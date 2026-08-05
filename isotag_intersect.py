#!/usr/bin/env python3
"""
isotag_intersect - Create Isotag Match Index (BAM-to-BAM)

Creates a lightweight SQLite index mapping isotag values to file offsets
in two BAM files. Enables efficient downstream queries and extraction.

Usage:
    python3 isotag_intersect.py -a experiment.bam -b reference.bam -o matches.db
    python3 isotag_intersect.py -a exp.bam -b ref.bam --tag XS -o matches.db

Index format:
    - SQLite database with pickled offset arrays
    - Minimal storage: isotag + [A offsets] + [B offsets]
    - ~10 MB for 1M reads with 45K unique isotags

Author: IsoTag Team
Version: 1.0.0
"""

import sqlite3
import pickle
import pysam
import os
import sys
import time
import json
import argparse
from collections import defaultdict, Counter
from datetime import datetime


class SimpleIsotagIndex:
    """
    Ultra-simple isotag index using SQLite + Pickle

    Stores only: isotag -> [A_offsets] + [B_offsets]
    No coordinates, no read names - maximum compactness!
    """

    def __init__(self, db_path):
        self.db_path = db_path
        self.conn = None
        self.cursor = None

    def create(self, bam_a_path, bam_b_path, tag_type):
        """Create new index database"""
        # Remove existing database if present
        if os.path.exists(self.db_path):
            os.remove(self.db_path)

        self.conn = sqlite3.connect(self.db_path)
        self.cursor = self.conn.cursor()

        # Create metadata table
        self.cursor.execute('''
            CREATE TABLE metadata (
                key TEXT PRIMARY KEY,
                value TEXT
            )
        ''')

        # Create matches table with tag_type
        # NOTE: Indexes are created AFTER data insertion for 10-100× speedup!
        self.cursor.execute('''
            CREATE TABLE matches (
                isotag TEXT NOT NULL,
                tag_type TEXT NOT NULL,
                a_offsets BLOB NOT NULL,
                b_offsets BLOB NOT NULL,
                PRIMARY KEY (isotag, tag_type)
            )
        ''')

        # Store metadata
        metadata = {
            'version': '2.0.0',
            'created': datetime.now().isoformat(),
            'bam_a': os.path.abspath(bam_a_path),
            'bam_b': os.path.abspath(bam_b_path),
            'multi_tag': 'true'
        }

        for key, value in metadata.items():
            self.cursor.execute(
                "INSERT INTO metadata (key, value) VALUES (?, ?)",
                (key, value)
            )

        self.conn.commit()

    def add_match(self, isotag, tag_type, a_offsets, b_offsets):
        """Add match record for this isotag and tag type"""
        # Pickle the offset arrays
        a_blob = pickle.dumps(a_offsets, protocol=pickle.HIGHEST_PROTOCOL)
        b_blob = pickle.dumps(b_offsets, protocol=pickle.HIGHEST_PROTOCOL)

        self.cursor.execute('''
            INSERT OR REPLACE INTO matches (isotag, tag_type, a_offsets, b_offsets)
            VALUES (?, ?, ?, ?)
        ''', (isotag, tag_type, a_blob, b_blob))

    def add_matches_batch(self, matches_dict, tag_type):
        """Add multiple matches efficiently for a specific tag type"""
        batch = []
        for isotag, (a_offsets, b_offsets) in matches_dict.items():
            a_blob = pickle.dumps(a_offsets, protocol=pickle.HIGHEST_PROTOCOL)
            b_blob = pickle.dumps(b_offsets, protocol=pickle.HIGHEST_PROTOCOL)
            batch.append((isotag, tag_type, a_blob, b_blob))

        self.cursor.executemany('''
            INSERT OR REPLACE INTO matches (isotag, tag_type, a_offsets, b_offsets)
            VALUES (?, ?, ?, ?)
        ''', batch)

    def begin_transaction(self):
        """Begin explicit transaction for bulk inserts"""
        self.cursor.execute("BEGIN TRANSACTION")

    def end_transaction(self):
        """End transaction and commit"""
        self.cursor.execute("END TRANSACTION")

    def create_indexes(self):
        """
        Create indexes AFTER data insertion for 10-100× speedup!

        Why? SQLite updates indexes on every INSERT, even within a transaction.
        Creating indexes after bulk insertion is much faster for large datasets.
        """
        print("Creating indexes (this may take a few minutes for large datasets)...")

        # Create indexes for fast lookups
        print("  Creating idx_isotag...", end='', flush=True)
        self.cursor.execute('CREATE INDEX idx_isotag ON matches(isotag)')
        print(" Done!")

        print("  Creating idx_tag_type...", end='', flush=True)
        self.cursor.execute('CREATE INDEX idx_tag_type ON matches(tag_type)')
        print(" Done!")

        print("  Creating idx_isotag_tag...", end='', flush=True)
        self.cursor.execute('CREATE INDEX idx_isotag_tag ON matches(isotag, tag_type)')
        print(" Done!")

        self.conn.commit()
        print("All indexes created successfully!")

    def update_metadata(self, key, value):
        """Update metadata value"""
        self.cursor.execute(
            "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
            (key, str(value))
        )

    def commit(self):
        """Commit changes"""
        if self.conn:
            self.conn.commit()

    def close(self):
        """Close database connection"""
        if self.conn:
            self.conn.commit()
            self.conn.close()
            self.conn = None

    def get_match(self, isotag, tag_type):
        """Get offsets for specific isotag and tag type (for query tool later)"""
        self.cursor.execute(
            "SELECT a_offsets, b_offsets FROM matches WHERE isotag = ? AND tag_type = ?",
            (isotag, tag_type)
        )
        result = self.cursor.fetchone()

        if result:
            a_offsets = pickle.loads(result[0])
            b_offsets = pickle.loads(result[1])
            return a_offsets, b_offsets

        return None, None

    def get_all_isotags(self, tag_type=None):
        """Get all isotag values (for query tool later)"""
        if tag_type:
            self.cursor.execute("SELECT DISTINCT isotag FROM matches WHERE tag_type = ?", (tag_type,))
        else:
            self.cursor.execute("SELECT DISTINCT isotag FROM matches")
        return [row[0] for row in self.cursor.fetchall()]

    def get_tag_types(self):
        """Get all tag types stored in the index"""
        self.cursor.execute("SELECT DISTINCT tag_type FROM matches")
        return [row[0] for row in self.cursor.fetchall()]

    def get_metadata(self):
        """Get all metadata (for query tool later)"""
        self.cursor.execute("SELECT key, value FROM metadata")
        return dict(self.cursor.fetchall())


class StatisticsCollector:
    """Collect statistics during index creation"""

    def __init__(self):
        self.start_time = time.time()

        # BAM A statistics
        self.total_reads_a = 0
        self.reads_with_tag_a = 0
        self.unique_isotags_a = set()

        # BAM B statistics
        self.total_reads_b = 0
        self.reads_with_tag_b = 0
        self.unique_isotags_b = set()

        # Matching statistics
        self.isotags_in_both = 0
        self.isotags_only_in_a = 0
        self.isotags_only_in_b = 0
        self.total_isotags = 0

    def update_bam_a(self, isotags_a):
        """Update statistics for BAM A"""
        self.unique_isotags_a = isotags_a

    def update_bam_b(self, isotags_b):
        """Update statistics for BAM B"""
        self.unique_isotags_b = isotags_b

    def finalize(self):
        """Calculate final statistics"""
        self.isotags_in_both = len(self.unique_isotags_a & self.unique_isotags_b)
        self.isotags_only_in_a = len(self.unique_isotags_a - self.unique_isotags_b)
        self.isotags_only_in_b = len(self.unique_isotags_b - self.unique_isotags_a)
        self.total_isotags = len(self.unique_isotags_a | self.unique_isotags_b)

    def get_elapsed_time(self):
        """Get elapsed time in seconds"""
        return time.time() - self.start_time

    def to_dict(self):
        """Return statistics as dictionary"""
        return {
            'bam_a': {
                'total_reads': self.total_reads_a,
                'reads_with_tag': self.reads_with_tag_a,
                'unique_isotags': len(self.unique_isotags_a)
            },
            'bam_b': {
                'total_reads': self.total_reads_b,
                'reads_with_tag': self.reads_with_tag_b,
                'unique_isotags': len(self.unique_isotags_b)
            },
            'matching': {
                'isotags_in_both': self.isotags_in_both,
                'isotags_only_in_a': self.isotags_only_in_a,
                'isotags_only_in_b': self.isotags_only_in_b,
                'total_unique_isotags': self.total_isotags
            },
            'performance': {
                'elapsed_time_seconds': round(self.get_elapsed_time(), 2),
                'reads_per_second': int((self.total_reads_a + self.total_reads_b) / self.get_elapsed_time())
            }
        }


def index_bam_by_all_isotags(bam_path, description="BAM"):
    """
    Index a BAM file by all isotag types (XI, XB, XS, XT, XV)

    Returns:
        dict: {tag_type: {isotag -> [list of file offsets]}}
        int: total reads processed
        dict: {tag_type: reads_with_tag}
    """
    print(f"Indexing {description}: {os.path.basename(bam_path)}")

    bam = pysam.AlignmentFile(bam_path, 'rb')

    # All isotag types we support
    TAG_TYPES = ['XI', 'XB', 'XS', 'XT', 'XV', 'XC', 'XR']

    # Create separate index for each tag type
    isotag_indexes = {tag: defaultdict(list) for tag in TAG_TYPES}
    reads_with_tag_counts = {tag: 0 for tag in TAG_TYPES}

    total_reads = 0

    while True:
        # CRITICAL: Get offset BEFORE reading
        offset = bam.tell()

        try:
            read = next(bam)
        except StopIteration:
            break

        total_reads += 1

        # Progress update (reduced frequency for better performance)
        if total_reads % 50000 == 0:
            print(f"  Processed {total_reads:,} reads...", end='\r')

        # Optimize: Get all tags at once to avoid multiple KeyError exceptions
        # This is 2-3× faster than 7 separate get_tag() calls with try/except
        tags = read.get_tags()
        tags_dict = {tag: value for tag, value in tags if tag in TAG_TYPES}

        # Extract all isotag types (no try/except needed - much faster!)
        for tag_type in TAG_TYPES:
            if tag_type in tags_dict:
                isotag_value = tags_dict[tag_type]
                isotag_indexes[tag_type][isotag_value].append(offset)
                reads_with_tag_counts[tag_type] += 1

    bam.close()

    print(f"  Processed {total_reads:,} reads - Done!     ")

    # Print summary for each tag type
    for tag_type in TAG_TYPES:
        count = reads_with_tag_counts[tag_type]
        unique = len(isotag_indexes[tag_type])
        if count > 0:
            print(f"  {tag_type} tag: {count:,} reads, {unique:,} unique values")

    return isotag_indexes, total_reads, reads_with_tag_counts


def create_intersect_index(bam_a_path, bam_b_path, output_db,
                           stats_path=None, quiet=False):
    """
    Create isotag intersection index for ALL isotag types

    Args:
        bam_a_path: Path to query BAM file
        bam_b_path: Path to reference BAM file
        output_db: Output SQLite database path
        stats_path: Optional path to write statistics JSON
        quiet: Suppress progress messages
    """

    # Redirect output if quiet
    if quiet:
        import io
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()

    stats = StatisticsCollector()

    print("=" * 70)
    print("isotag_intersect - Create Isotag Match Index (Multi-Tag)")
    print("=" * 70)
    print(f"Query BAM (A):     {bam_a_path}")
    print(f"Reference BAM (B): {bam_b_path}")
    print(f"Storing tags:      XI, XB, XS, XT, XV, XC, XR (all tags)")
    print(f"Output database:   {output_db}")
    print("=" * 70)
    print()

    # Step 1: Index BAM B (reference) - all tags
    print("STEP 1: Indexing reference BAM (B) - all tags")
    print("-" * 70)
    b_indexes, total_b, with_tag_b_dict = index_bam_by_all_isotags(bam_b_path, "reference BAM")
    stats.total_reads_b = total_b
    print()

    # Step 2: Index BAM A (query) - all tags
    print("STEP 2: Indexing query BAM (A) - all tags")
    print("-" * 70)
    a_indexes, total_a, with_tag_a_dict = index_bam_by_all_isotags(bam_a_path, "query BAM")
    stats.total_reads_a = total_a
    print()

    # Step 3: Create database and write matches
    print("STEP 3: Creating index database")
    print("-" * 70)

    db_index = SimpleIsotagIndex(output_db)
    db_index.create(bam_a_path, bam_b_path, None)  # No single tag_type anymore

    # Get all tag types present in at least one BAM file
    TAG_TYPES = ['XI', 'XB', 'XS', 'XT', 'XV', 'XC', 'XR']
    active_tags = []
    for tag in TAG_TYPES:
        if len(a_indexes[tag]) > 0 or len(b_indexes[tag]) > 0:
            active_tags.append(tag)

    print(f"Active tag types: {', '.join(active_tags)}")
    print()

    # Begin transaction BEFORE all inserts (10-100× speedup!)
    db_index.begin_transaction()

    # Process each tag type separately
    total_written = 0
    for tag_type in active_tags:
        a_index = a_indexes[tag_type]
        b_index = b_indexes[tag_type]

        # Combine all isotags for this tag type
        # Optimize: Convert to sorted list for better memory locality
        all_isotags = sorted(set(a_index.keys()) | set(b_index.keys()))

        if len(all_isotags) == 0:
            continue

        print(f"Writing {tag_type} tag: {len(all_isotags):,} unique values")

        # Batch write for efficiency (increased from 1000 to 5000 for better performance)
        batch_size = 5000
        batch = {}
        written = 0

        for isotag in all_isotags:
            a_offsets = a_index.get(isotag, [])
            b_offsets = b_index.get(isotag, [])

            batch[isotag] = (a_offsets, b_offsets)

            if len(batch) >= batch_size:
                db_index.add_matches_batch(batch, tag_type)
                written += len(batch)
                total_written += len(batch)
                # Progress update every batch (now every 5000 instead of 1000)
                if written % 5000 == 0 or written == len(all_isotags):
                    print(f"  {tag_type}: {written:,} / {len(all_isotags):,} records", end='\r')
                batch = {}

        # Write remaining batch
        if batch:
            db_index.add_matches_batch(batch, tag_type)
            written += len(batch)
            total_written += len(batch)

        print(f"  {tag_type}: {written:,} / {len(all_isotags):,} records - Done!     ")

    # End transaction AFTER all inserts (commit once!)
    db_index.end_transaction()

    print()
    print(f"Total records written: {total_written:,}")
    print()

    # STEP 4: Create indexes AFTER data insertion (10-100× speedup!)
    print("STEP 4: Creating indexes")
    print("-" * 70)
    db_index.create_indexes()
    print()

    # Update metadata with statistics
    stats_dict = stats.to_dict()
    for category, values in stats_dict.items():
        for key, value in values.items():
            db_index.update_metadata(f"{category}_{key}", value)

    db_index.commit()
    db_index.close()

    # Get database size
    db_size_mb = os.path.getsize(output_db) / (1024 * 1024)

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"BAM A (query):")
    print(f"  Total reads:        {stats.total_reads_a:,}")
    for tag in active_tags:
        count = with_tag_a_dict[tag]
        unique = len(a_indexes[tag])
        print(f"  {tag} tag:             {count:,} reads, {unique:,} unique")
    print()
    print(f"BAM B (reference):")
    print(f"  Total reads:        {stats.total_reads_b:,}")
    for tag in active_tags:
        count = with_tag_b_dict[tag]
        unique = len(b_indexes[tag])
        print(f"  {tag} tag:             {count:,} reads, {unique:,} unique")
    print()
    print(f"Performance:")
    print(f"  Processing time:    {stats.get_elapsed_time():.1f} seconds")
    print(f"  Processing speed:   {stats_dict['performance']['reads_per_second']:,} reads/sec")
    print(f"  Database size:      {db_size_mb:.2f} MB")
    print(f"  Total records:      {total_written:,}")
    print("=" * 70)
    print(f"✓ Index created successfully: {output_db}")
    print(f"✓ Stored {len(active_tags)} tag types: {', '.join(active_tags)}")
    print()

    # Write statistics JSON if requested
    if stats_path:
        stats_output = {
            'created': datetime.now().isoformat(),
            'version': '1.0.0',
            'input_files': {
                'bam_a': os.path.abspath(bam_a_path),
                'bam_b': os.path.abspath(bam_b_path)
            },
            'tag_type': tag_type,
            'statistics': stats_dict,
            'database': {
                'path': os.path.abspath(output_db),
                'size_mb': round(db_size_mb, 2)
            }
        }

        with open(stats_path, 'w') as f:
            json.dump(stats_output, f, indent=2)

        print(f"✓ Statistics written to: {stats_path}")

    # Restore stdout if quiet
    if quiet:
        sys.stdout = old_stdout


def main():
    """Main entry point"""

    parser = argparse.ArgumentParser(
        description='isotag_intersect - Create isotag match index (BAM-to-BAM)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - stores ALL isotag types automatically
  isotag_intersect.py -a experiment.bam -b reference.bam -o matches.db

  # With statistics output
  isotag_intersect.py -a exp.bam -b ref.bam -o matches.db --stats stats.json

Database format:
  - SQLite with pickled offset arrays
  - Stores ALL tag types: XI, XB, XS, XT, XV, XC, XR
  - Minimal storage: (isotag, tag_type) -> [A offsets] + [B offsets]
  - ~30-40 MB for 1M reads with 45K unique isotags per tag

Notes:
  - Input BAM files must be tagged with isotag.py first
  - Automatically stores all available isotag types
  - Output is a SQLite database for efficient queries
  - Use isotag_query_index.py to query by specific tag type
        """
    )

    parser.add_argument('-a', required=True, metavar='FILE',
                       help='Input BAM file A (query)')
    parser.add_argument('-b', required=True, metavar='FILE',
                       help='Input BAM file B (reference)')
    parser.add_argument('-o', required=True, metavar='FILE',
                       help='Output SQLite database (.db extension)')
    parser.add_argument('--stats', metavar='FILE',
                       help='Write statistics to JSON file')
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='Suppress progress messages')
    parser.add_argument('-v', '--version', action='version',
                       version='isotag_intersect 2.0.0')

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.a):
        print(f"Error: Input file not found: {args.a}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.b):
        print(f"Error: Reference file not found: {args.b}", file=sys.stderr)
        sys.exit(1)

    # Check BAM files are readable
    try:
        bam = pysam.AlignmentFile(args.a, 'rb')
        bam.close()
    except Exception as e:
        print(f"Error: Cannot open BAM file A: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        bam = pysam.AlignmentFile(args.b, 'rb')
        bam.close()
    except Exception as e:
        print(f"Error: Cannot open BAM file B: {e}", file=sys.stderr)
        sys.exit(1)

    # Check output directory is writable
    output_dir = os.path.dirname(os.path.abspath(args.o))
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except Exception as e:
            print(f"Error: Cannot create output directory: {e}", file=sys.stderr)
            sys.exit(1)

    # Create index
    try:
        create_intersect_index(
            args.a,
            args.b,
            args.o,
            stats_path=args.stats,
            quiet=args.quiet
        )
    except KeyboardInterrupt:
        print("\n\nInterrupted by user", file=sys.stderr)
        # Clean up partial database
        if os.path.exists(args.o):
            os.remove(args.o)
        sys.exit(1)
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        # Clean up partial database
        if os.path.exists(args.o):
            os.remove(args.o)
        sys.exit(1)


if __name__ == '__main__':
    main()
