#!/usr/bin/env python3
"""
isotag_query_index - Query Isotag Intersection Index

Query and display information from isotag intersection index created
by isotag_intersect.py. Shows read IDs, sequences, and other information
for specific isotags.

Usage:
    # Query specific isotag
    python3 isotag_query_index.py --index matches.db --isotag fuIF7PN...

    # Show sequences
    python3 isotag_query_index.py --index matches.db --isotag fuIF7PN... --format sequence

    # Expanded format (cartesian product)
    python3 isotag_query_index.py --index matches.db --isotag fuIF7PN... --expand

    # List all isotags
    python3 isotag_query_index.py --index matches.db --list

Author: IsoTag Team
Version: 1.0.0
"""

import sqlite3
import pickle
import pysam
import os
import sys
import argparse
from collections import defaultdict


class IsotagIndexReader:
    """Read and query isotag intersection index"""

    def __init__(self, db_path):
        self.db_path = db_path
        self.conn = None
        self.cursor = None
        self.metadata = {}

        # BAM file handles (opened on demand)
        self.bam_a = None
        self.bam_b = None
        self.bam_a_path = None
        self.bam_b_path = None

    def open(self):
        """Open database and read metadata"""
        if not os.path.exists(self.db_path):
            raise FileNotFoundError(f"Index database not found: {self.db_path}")

        self.conn = sqlite3.connect(self.db_path)
        self.cursor = self.conn.cursor()

        # Read metadata
        self.cursor.execute("SELECT key, value FROM metadata")
        self.metadata = dict(self.cursor.fetchall())

        self.bam_a_path = self.metadata.get('bam_a')
        self.bam_b_path = self.metadata.get('bam_b')

        # Validate BAM files exist
        if not os.path.exists(self.bam_a_path):
            print(f"Warning: BAM A not found: {self.bam_a_path}", file=sys.stderr)
            print(f"  (Index created with this path, but file may have moved)", file=sys.stderr)

        if not os.path.exists(self.bam_b_path):
            print(f"Warning: BAM B not found: {self.bam_b_path}", file=sys.stderr)
            print(f"  (Index created with this path, but file may have moved)", file=sys.stderr)

    def open_bam_files(self):
        """Open BAM files for reading"""
        if self.bam_a is None and os.path.exists(self.bam_a_path):
            self.bam_a = pysam.AlignmentFile(self.bam_a_path, 'rb')

        if self.bam_b is None and os.path.exists(self.bam_b_path):
            self.bam_b = pysam.AlignmentFile(self.bam_b_path, 'rb')

    def close(self):
        """Close database and BAM files"""
        if self.conn:
            self.conn.close()
        if self.bam_a:
            self.bam_a.close()
        if self.bam_b:
            self.bam_b.close()

    def get_match(self, isotag, tag_type='XI'):
        """Get offsets for specific isotag and tag type"""
        self.cursor.execute(
            "SELECT a_offsets, b_offsets FROM matches WHERE isotag = ? AND tag_type = ?",
            (isotag, tag_type)
        )
        result = self.cursor.fetchone()

        if result is None:
            return None, None

        # Unpickle offset arrays
        a_offsets = pickle.loads(result[0])
        b_offsets = pickle.loads(result[1])

        return a_offsets, b_offsets

    def get_all_isotags(self, tag_type='XI'):
        """Get all isotags with counts for specific tag type"""
        self.cursor.execute(
            "SELECT isotag, a_offsets, b_offsets FROM matches WHERE tag_type = ?",
            (tag_type,)
        )

        isotags = []
        for row in self.cursor.fetchall():
            isotag = row[0]
            a_count = len(pickle.loads(row[1]))
            b_count = len(pickle.loads(row[2]))
            isotags.append((isotag, a_count, b_count))

        return isotags

    def get_tag_types(self):
        """Get all tag types stored in the index"""
        self.cursor.execute("SELECT DISTINCT tag_type FROM matches")
        return [row[0] for row in self.cursor.fetchall()]

    def get_read_by_offset(self, offset, from_bam='a'):
        """Get read from BAM file using offset"""
        self.open_bam_files()

        if from_bam == 'a':
            if self.bam_a is None:
                raise FileNotFoundError(f"Cannot open BAM A: {self.bam_a_path}")
            bam = self.bam_a
        else:
            if self.bam_b is None:
                raise FileNotFoundError(f"Cannot open BAM B: {self.bam_b_path}")
            bam = self.bam_b

        # Seek to offset and read
        bam.seek(offset)
        read = next(bam)
        return read

    def get_reads_by_offsets(self, offsets, from_bam='a'):
        """Get multiple reads from BAM file"""
        reads = []
        for offset in offsets:
            try:
                read = self.get_read_by_offset(offset, from_bam)
                reads.append(read)
            except Exception as e:
                print(f"Warning: Could not read at offset {offset}: {e}", file=sys.stderr)
                reads.append(None)
        return reads


def format_read_id(read):
    """Get read ID"""
    if read is None:
        return "NULL"
    return read.query_name


def format_read_sequence(read):
    """Get read sequence"""
    if read is None:
        return "NULL"
    return read.query_sequence if read.query_sequence else "NULL"


def format_read_cigar(read):
    """Get CIGAR string"""
    if read is None:
        return "NULL"
    return read.cigarstring if read.cigarstring else "NULL"


def print_compact_format(isotag, a_reads, b_reads, format_fields):
    """Print compact format (comma-separated)"""

    # Build output based on format fields
    output = [isotag]

    if 'id' in format_fields:
        a_ids = ','.join([format_read_id(r) for r in a_reads])
        b_ids = ','.join([format_read_id(r) for r in b_reads])
        output.extend([str(len(a_reads)), str(len(b_reads)), a_ids, b_ids])

    if 'sequence' in format_fields:
        a_seqs = ','.join([format_read_sequence(r) for r in a_reads])
        b_seqs = ','.join([format_read_sequence(r) for r in b_reads])
        output.extend([a_seqs, b_seqs])

    if 'cigar' in format_fields:
        a_cigars = ','.join([format_read_cigar(r) for r in a_reads])
        b_cigars = ','.join([format_read_cigar(r) for r in b_reads])
        output.extend([a_cigars, b_cigars])

    print('\t'.join(output))


def print_expanded_format(isotag, a_reads, b_reads, format_fields):
    """Print expanded format (cartesian product)"""

    # Generate all pairs
    for a_read in a_reads:
        for b_read in b_reads:
            output = [isotag]

            if 'id' in format_fields:
                output.extend([format_read_id(a_read), format_read_id(b_read)])

            if 'sequence' in format_fields:
                output.extend([format_read_sequence(a_read), format_read_sequence(b_read)])

            if 'cigar' in format_fields:
                output.extend([format_read_cigar(a_read), format_read_cigar(b_read)])

            print('\t'.join(output))


def print_header_compact(format_fields):
    """Print header for compact format"""
    header = ['isotag']

    if 'id' in format_fields:
        header.extend(['a_count', 'b_count', 'a_ids', 'b_ids'])

    if 'sequence' in format_fields:
        header.extend(['a_sequences', 'b_sequences'])

    if 'cigar' in format_fields:
        header.extend(['a_cigars', 'b_cigars'])

    print('\t'.join(header))


def print_header_expanded(format_fields):
    """Print header for expanded format"""
    header = ['isotag']

    if 'id' in format_fields:
        header.extend(['a_id', 'b_id'])

    if 'sequence' in format_fields:
        header.extend(['a_sequence', 'b_sequence'])

    if 'cigar' in format_fields:
        header.extend(['a_cigar', 'b_cigar'])

    print('\t'.join(header))


def query_isotag(index_reader, isotag, tag_type, format_str, expand):
    """Query and print information for specific isotag"""

    # Parse format string
    format_fields = [f.strip() for f in format_str.split(',')]

    # Validate format fields
    valid_fields = {'id', 'sequence', 'cigar'}
    for field in format_fields:
        if field not in valid_fields:
            print(f"Error: Invalid format field '{field}'. Valid: {valid_fields}", file=sys.stderr)
            sys.exit(1)

    # Get offsets
    a_offsets, b_offsets = index_reader.get_match(isotag, tag_type)

    if a_offsets is None:
        print(f"Error: Isotag not found in index for tag type {tag_type}: {isotag}", file=sys.stderr)
        print(f"Available tag types: {', '.join(index_reader.get_tag_types())}", file=sys.stderr)
        sys.exit(1)

    # Fetch reads if needed (for sequence/cigar formats)
    if 'id' in format_fields and len(format_fields) == 1:
        # ID only - no need to fetch reads, can use offsets directly
        # But for consistency, we'll fetch read names
        # Actually, we need to fetch to get read names
        pass

    # Fetch reads from BAM files
    a_reads = index_reader.get_reads_by_offsets(a_offsets, from_bam='a')
    b_reads = index_reader.get_reads_by_offsets(b_offsets, from_bam='b')

    # Print header
    if expand:
        print_header_expanded(format_fields)
        print_expanded_format(isotag, a_reads, b_reads, format_fields)
    else:
        print_header_compact(format_fields)
        print_compact_format(isotag, a_reads, b_reads, format_fields)


def list_all_isotags(index_reader, tag_type, top_n=None):
    """List all isotags with counts for specific tag type"""

    print(f"Getting all isotags from index (tag type: {tag_type})...", file=sys.stderr)
    isotags = index_reader.get_all_isotags(tag_type)

    # Sort by total count (A + B) descending
    isotags.sort(key=lambda x: x[1] + x[2], reverse=True)

    # Limit to top N if specified
    if top_n:
        isotags = isotags[:top_n]

    # Print header
    print('\t'.join(['isotag', 'tag_type', 'a_count', 'b_count', 'total', 'status']))

    # Print isotags
    for isotag, a_count, b_count in isotags:
        total = a_count + b_count

        # Determine status
        if a_count > 0 and b_count > 0:
            status = "in_both"
        elif a_count > 0:
            status = "only_in_a"
        else:
            status = "only_in_b"

        print('\t'.join([isotag, tag_type, str(a_count), str(b_count), str(total), status]))

    print(f"\nTotal isotags for {tag_type}: {len(isotags)}", file=sys.stderr)


def compare_all_tags(index_reader):
    """Compare matching rates across all tag types"""

    tag_types = sorted(index_reader.get_tag_types())

    if not tag_types:
        print("Error: No tag types found in index", file=sys.stderr)
        return

    print("="*80)
    print("Multi-Tag Comparison Report")
    print("="*80)
    print()

    # Collect statistics for all tags
    stats = []
    for tag in tag_types:
        isotags = index_reader.get_all_isotags(tag)

        in_both = sum(1 for _, a, b in isotags if a > 0 and b > 0)
        only_a = sum(1 for _, a, b in isotags if a > 0 and b == 0)
        only_b = sum(1 for _, a, b in isotags if a == 0 and b > 0)
        total = len(isotags)
        match_rate = (in_both / total * 100) if total > 0 else 0

        # Count total reads (not just unique isotags)
        reads_in_both_a = sum(a for _, a, b in isotags if a > 0 and b > 0)
        reads_in_both_b = sum(b for _, a, b in isotags if a > 0 and b > 0)
        reads_only_a = sum(a for _, a, b in isotags if a > 0 and b == 0)
        reads_only_b = sum(b for _, a, b in isotags if a == 0 and b > 0)

        stats.append({
            'tag': tag,
            'in_both': in_both,
            'only_a': only_a,
            'only_b': only_b,
            'total': total,
            'match_rate': match_rate,
            'reads_in_both_a': reads_in_both_a,
            'reads_in_both_b': reads_in_both_b,
            'reads_only_a': reads_only_a,
            'reads_only_b': reads_only_b
        })

    # Find best tag by match rate
    best_tag = max(stats, key=lambda x: x['match_rate'])

    # Print header
    print(f"{'Tag':<6} {'Isotags_Both':>13} {'Only_A':>10} {'Only_B':>10} {'Total':>10} {'Match_Rate':>12}")
    print("-"*80)

    # Print stats
    for s in stats:
        marker = " ← BEST" if s['tag'] == best_tag['tag'] else ""
        print(f"{s['tag']:<6} {s['in_both']:>13,} {s['only_a']:>10,} {s['only_b']:>10,} "
              f"{s['total']:>10,} {s['match_rate']:>11.1f}%{marker}")

    print()
    print("Read counts (not unique isotags):")
    print("-"*80)
    print(f"{'Tag':<6} {'Reads_Both_A':>14} {'Reads_Both_B':>14} {'Reads_Only_A':>14} {'Reads_Only_B':>14}")
    print("-"*80)

    for s in stats:
        print(f"{s['tag']:<6} {s['reads_in_both_a']:>14,} {s['reads_in_both_b']:>14,} "
              f"{s['reads_only_a']:>14,} {s['reads_only_b']:>14,}")

    print()
    print("="*80)
    print("ANALYSIS")
    print("="*80)

    # Recommendation
    print(f"✓ Best match rate: {best_tag['tag']} ({best_tag['match_rate']:.1f}%)")

    # Compare XI vs XS (if both present)
    xi_stats = next((s for s in stats if s['tag'] == 'XI'), None)
    xs_stats = next((s for s in stats if s['tag'] == 'XS'), None)

    if xi_stats and xs_stats:
        rescued = xs_stats['in_both'] - xi_stats['in_both']
        if rescued > 0:
            print(f"✓ XS fuzzy matching rescued {rescued:,} additional isotags vs XI exact matching")
            print(f"  ({rescued / xi_stats['total'] * 100:.1f}% increase in matches)")
        elif rescued < 0:
            print(f"⚠ XI has {abs(rescued):,} more matches than XS (unexpected!)")

    # Compare XT vs XI (biological vs structural)
    xt_stats = next((s for s in stats if s['tag'] == 'XT'), None)
    if xi_stats and xt_stats:
        diff = xi_stats['in_both'] - xt_stats['in_both']
        if diff > 0:
            print(f"✓ XI found {diff:,} more structural matches than XT transcript groups")
            print(f"  (same gene may have multiple isoforms)")

    print()
    print("RECOMMENDATIONS:")
    print("-"*80)
    if best_tag['match_rate'] > 90:
        print(f"• High match rate ({best_tag['match_rate']:.1f}%) - excellent data quality!")
    elif best_tag['match_rate'] > 75:
        print(f"• Good match rate ({best_tag['match_rate']:.1f}%) - acceptable results")
    else:
        print(f"⚠ Low match rate ({best_tag['match_rate']:.1f}%) - check data quality")

    # Tag-specific recommendations
    if xs_stats and xs_stats['tag'] == best_tag['tag']:
        print(f"• Use XS for maximum sensitivity (fuzzy ±5bp splice junction matching)")
    elif xi_stats and xi_stats['tag'] == best_tag['tag']:
        print(f"• Use XI for maximum precision (exact exon structure matching)")
    elif xt_stats and xt_stats['tag'] == best_tag['tag']:
        print(f"• Use XT for biological transcript grouping (gene-level analysis)")

    print("="*80)
    print()


def show_index_info(index_reader):
    """Show information about the index"""

    print("="*70)
    print("Index Information")
    print("="*70)

    # Basic metadata
    print(f"Database:     {index_reader.db_path}")
    print(f"Version:      {index_reader.metadata.get('version', 'unknown')}")
    print(f"Created:      {index_reader.metadata.get('created', 'unknown')}")

    # Get tag types
    tag_types = index_reader.get_tag_types()
    print(f"Tag types:    {', '.join(tag_types)} ({len(tag_types)} types)")
    print()

    print("Input files:")
    print(f"  BAM A:      {index_reader.bam_a_path}")
    print(f"              {'✓ exists' if os.path.exists(index_reader.bam_a_path) else '✗ not found'}")
    print(f"  BAM B:      {index_reader.bam_b_path}")
    print(f"              {'✓ exists' if os.path.exists(index_reader.bam_b_path) else '✗ not found'}")
    print()

    # Statistics from metadata
    print("Statistics:")
    stats_keys = [
        ('bam_a_total_reads', 'BAM A total reads'),
        ('bam_a_reads_with_tag', 'BAM A with tag'),
        ('bam_a_unique_isotags', 'BAM A unique isotags'),
        ('bam_b_total_reads', 'BAM B total reads'),
        ('bam_b_reads_with_tag', 'BAM B with tag'),
        ('bam_b_unique_isotags', 'BAM B unique isotags'),
        ('matching_isotags_in_both', 'Isotags in both'),
        ('matching_isotags_only_in_a', 'Isotags only in A'),
        ('matching_isotags_only_in_b', 'Isotags only in B'),
    ]

    for key, label in stats_keys:
        value = index_reader.metadata.get(key, 'N/A')
        if value != 'N/A':
            value = f"{int(value):,}"
        print(f"  {label:25} {value}")

    print("="*70)


def main():
    """Main entry point"""

    parser = argparse.ArgumentParser(
        description='isotag_query_index - Query isotag intersection index',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # List all isotags (default: XI tag)
  isotag_query_index.py --index matches.db --list

  # List XS (splice) isotags
  isotag_query_index.py --index matches.db --list --tag XS

  # Show top 10 most abundant XI isotags
  isotag_query_index.py --index matches.db --list --tag XI --top 10

  # Query specific isotag by XI (structure)
  isotag_query_index.py --index matches.db --isotag fuIF7PN23g2gq9sFxqhUNGnfOCZhkQJS

  # Query same isotag by XS (splice pattern)
  isotag_query_index.py --index matches.db --isotag aKF498dAp... --tag XS

  # Show sequences (compact)
  isotag_query_index.py --index matches.db --isotag fuIF7PN... --format sequence

  # Expanded format (cartesian product)
  isotag_query_index.py --index matches.db --isotag fuIF7PN... --expand

  # Expanded with sequences and custom tag
  isotag_query_index.py --index matches.db --isotag fuIF7PN... --tag XT --format sequence --expand

  # Show index information (all tag types)
  isotag_query_index.py --index matches.db --info

Tag types:
  XI - Exact exon structure (default)
  XS - Splice junctions
  XB - 5'/3' boundaries
  XT - Transcript groups
  XV - Variants
  XC - Cluster IDs
  XR - Representative tags

Format options:
  id         - Read names only (default)
  sequence   - Read sequences only
  cigar      - CIGAR strings only
  id,sequence       - Both IDs and sequences
  id,sequence,cigar - All three

Display modes:
  Default (compact):  Comma-separated lists
  --expand:           Cartesian product (one line per A-B pair)

        """
    )

    parser.add_argument('--index', required=True, metavar='FILE',
                       help='Input index database (.db)')

    # Query modes
    query_group = parser.add_mutually_exclusive_group()
    query_group.add_argument('--isotag', metavar='ID',
                            help='Query specific isotag')
    query_group.add_argument('--list', action='store_true',
                            help='List all isotags with counts')
    query_group.add_argument('--info', action='store_true',
                            help='Show index information')
    query_group.add_argument('--compare', action='store_true',
                            help='Compare matching rates across all tag types')

    # Tag selection
    parser.add_argument('--tag', default='XI',
                       choices=['XI', 'XB', 'XS', 'XT', 'XV', 'XC', 'XR'],
                       help='Tag type to query (default: XI)')

    # Format options
    parser.add_argument('--format', default='id',
                       help='Output format: id, sequence, cigar, id,sequence, etc (default: id)')
    parser.add_argument('--expand', action='store_true',
                       help='Expanded format (cartesian product)')

    # List options
    parser.add_argument('--top', type=int, metavar='N',
                       help='Show top N isotags (with --list)')

    parser.add_argument('-v', '--version', action='version',
                       version='isotag_query_index 2.0.0')

    args = parser.parse_args()

    # Validate index file exists
    if not os.path.exists(args.index):
        print(f"Error: Index file not found: {args.index}", file=sys.stderr)
        sys.exit(1)

    # Open index
    try:
        reader = IsotagIndexReader(args.index)
        reader.open()
    except Exception as e:
        print(f"Error: Cannot open index: {e}", file=sys.stderr)
        sys.exit(1)

    # Execute query
    try:
        if args.info:
            # Show index information
            show_index_info(reader)

        elif args.compare:
            # Compare all tag types
            compare_all_tags(reader)

        elif args.list:
            # List all isotags for specific tag type
            list_all_isotags(reader, args.tag, top_n=args.top)

        elif args.isotag:
            # Query specific isotag with tag type
            query_isotag(reader, args.isotag, args.tag, args.format, args.expand)

        else:
            # No query specified
            print("Error: Must specify --isotag, --list, --info, or --compare", file=sys.stderr)
            parser.print_help()
            sys.exit(1)

    except KeyboardInterrupt:
        print("\n\nInterrupted by user", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
    finally:
        reader.close()


if __name__ == '__main__':
    main()
