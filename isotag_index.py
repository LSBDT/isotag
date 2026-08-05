#!/usr/bin/env python3
"""
ISO-Tools Index - Create and Query Isotag Databases

Create searchable databases of isotag data for fast lookup and retrieval.
Support for complex queries and batch operations.
"""

import subprocess
import click
import sys
import sqlite3
import json
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import re
import time


class IsotogDatabase:
    """SQLite-based isotag database for fast queries"""
    
    def __init__(self, db_path: str):
        self.db_path = db_path
        self.conn = None
        self.cursor = None
    
    def connect(self):
        """Connect to the database"""
        self.conn = sqlite3.connect(self.db_path)
        self.cursor = self.conn.cursor()
        self.conn.execute("PRAGMA journal_mode=WAL")  # Better concurrent access
        self.conn.execute("PRAGMA synchronous=NORMAL")  # Better performance
    
    def disconnect(self):
        """Disconnect from the database"""
        if self.conn:
            self.conn.close()
    
    def create_schema(self):
        """Create database schema"""
        # Main isotag table
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS isotags (
                isotag_id TEXT PRIMARY KEY,
                chrom TEXT NOT NULL,
                strand TEXT NOT NULL,
                start_pos INTEGER NOT NULL,
                end_pos INTEGER NOT NULL,
                exon_count INTEGER NOT NULL,
                total_length INTEGER NOT NULL,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        
        # Exons table
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS exons (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                isotag_id TEXT NOT NULL,
                exon_number INTEGER NOT NULL,
                start_pos INTEGER NOT NULL,
                end_pos INTEGER NOT NULL,
                FOREIGN KEY (isotag_id) REFERENCES isotags (isotag_id)
            )
        """)
        
        # Sample data table
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS sample_data (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                isotag_id TEXT NOT NULL,
                sample_name TEXT NOT NULL,
                read_count INTEGER NOT NULL,
                total_reads INTEGER NOT NULL,
                bam_file TEXT,
                FOREIGN KEY (isotag_id) REFERENCES isotags (isotag_id)
            )
        """)
        
        # Variant data table
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS variants (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                isotag_id TEXT NOT NULL,
                variant_id TEXT,
                variant_count INTEGER DEFAULT 1,
                FOREIGN KEY (isotag_id) REFERENCES isotags (isotag_id)
            )
        """)
        
        # Annotations table
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS annotations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                isotag_id TEXT NOT NULL,
                annotation_type TEXT NOT NULL,
                gene_id TEXT,
                transcript_id TEXT,
                classification TEXT,
                overlap_score REAL,
                FOREIGN KEY (isotag_id) REFERENCES isotags (isotag_id)
            )
        """)
        
        # Create indexes for faster queries
        indexes = [
            "CREATE INDEX IF NOT EXISTS idx_isotags_chrom_strand ON isotags (chrom, strand)",
            "CREATE INDEX IF NOT EXISTS idx_isotags_position ON isotags (chrom, start_pos, end_pos)",
            "CREATE INDEX IF NOT EXISTS idx_exons_isotag ON exons (isotag_id)",
            "CREATE INDEX IF NOT EXISTS idx_sample_data_isotag ON sample_data (isotag_id)",
            "CREATE INDEX IF NOT EXISTS idx_sample_data_sample ON sample_data (sample_name)",
            "CREATE INDEX IF NOT EXISTS idx_variants_isotag ON variants (isotag_id)",
            "CREATE INDEX IF NOT EXISTS idx_annotations_isotag ON annotations (isotag_id)",
            "CREATE INDEX IF NOT EXISTS idx_annotations_gene ON annotations (gene_id)",
        ]
        
        for index_sql in indexes:
            self.cursor.execute(index_sql)
        
        self.conn.commit()


class IsotogIndexer:
    """Index and query isotag databases"""
    
    def __init__(self):
        self.db = None
        self.stats = {
            'isotags_indexed': 0,
            'samples_processed': 0,
            'exons_indexed': 0,
            'variants_indexed': 0,
            'annotations_indexed': 0
        }
    
    def create_database(self, db_path: str):
        """Create new isotag database"""
        click.echo(f"🗄️  Creating isotag database: {db_path}")
        
        # Remove existing database
        db_file = Path(db_path)
        if db_file.exists():
            db_file.unlink()
            click.echo(f"   ⚠️  Removed existing database")
        
        self.db = IsotogDatabase(db_path)
        self.db.connect()
        self.db.create_schema()
        
        click.echo(f"   ✅ Database created with schema")
    
    def index_bam_files(self, bam_files: List[str], sample_names: Optional[List[str]] = None):
        """Index isotag data from BAM files"""
        if not sample_names:
            sample_names = [Path(f).stem for f in bam_files]
        
        if len(bam_files) != len(sample_names):
            click.echo("❌ Number of BAM files and sample names must match")
            sys.exit(1)
        
        click.echo(f"📋 Indexing {len(bam_files)} BAM files...")
        
        isotag_structures = {}  # isotag_id -> structure_info
        
        for bam_file, sample_name in zip(bam_files, sample_names):
            click.echo(f"   📁 Processing {sample_name}...")
            
            isotag_counts = Counter()
            isotag_reads = defaultdict(list)
            
            view_cmd = ['samtools', 'view', bam_file]
            
            try:
                process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)
                
                reads_processed = 0
                for line in process.stdout:
                    reads_processed += 1
                    
                    if reads_processed % 10000 == 0:
                        click.echo(f"      ⏳ Processed {reads_processed:,} reads...")
                    
                    fields = line.strip().split('\t')
                    
                    if len(fields) >= 11:
                        flag = int(fields[1])
                        rname = fields[2]
                        pos = int(fields[3])
                        cigar = fields[5]
                        
                        # Skip unmapped reads
                        if flag & 0x4 or cigar == '*':
                            continue
                        
                        # Extract isotag and variant info
                        isotag_id = None
                        variant_id = None
                        
                        for field in fields[11:]:
                            if field.startswith('XI:Z:'):
                                isotag_id = field[5:]
                            elif field.startswith('XV:Z:'):
                                variant_id = field[5:]
                        
                        if isotag_id:
                            isotag_counts[isotag_id] += 1
                            
                            # Store structure info (first time we see this isotag)
                            if isotag_id not in isotag_structures:
                                strand = '-' if flag & 0x10 else '+'
                                exons = self.parse_cigar_to_blocks(pos, cigar)
                                
                                if exons:
                                    isotag_structures[isotag_id] = {
                                        'chrom': rname,
                                        'strand': strand,
                                        'exons': exons,
                                        'start': min(e[0] for e in exons),
                                        'end': max(e[1] for e in exons)
                                    }
                            
                            # Store variant info
                            if variant_id:
                                isotag_reads[isotag_id].append(variant_id)
                
                process.wait()
                
                # Insert sample data into database
                self.insert_sample_data(isotag_counts, sample_name, bam_file, len(isotag_reads))
                
                click.echo(f"      ✅ Indexed {len(isotag_counts)} isotags, {reads_processed:,} total reads")
                self.stats['samples_processed'] += 1
                
            except subprocess.CalledProcessError as e:
                click.echo(f"❌ Error reading {bam_file}: {e}")
                sys.exit(1)
        
        # Insert isotag structures
        self.insert_isotag_structures(isotag_structures)
        
        # Insert variant data
        self.insert_variant_data(isotag_reads)
    
    def parse_cigar_to_blocks(self, pos: int, cigar: str) -> List[Tuple[int, int]]:
        """Parse CIGAR to genomic blocks (exons)"""
        blocks = []
        current_pos = pos
        block_start = pos
        in_block = True
        
        pattern = r'(\d+)([MIDNSHP=X])'
        matches = re.findall(pattern, cigar)
        
        for length_str, op_char in matches:
            length = int(length_str)
            
            if op_char in ['M', '=', 'X']:  # Match/alignment
                if not in_block:
                    block_start = current_pos
                    in_block = True
                current_pos += length
            
            elif op_char == 'D':  # Deletion from reference
                current_pos += length
            
            elif op_char == 'I':  # Insertion to reference
                pass  # Don't advance reference position
            
            elif op_char == 'N':  # Intron/skipped region
                if in_block:
                    blocks.append((block_start, current_pos - 1))
                    in_block = False
                current_pos += length
            
            elif op_char in ['S', 'H']:  # Soft/hard clipping
                pass  # Don't affect block structure
        
        # Add final block
        if in_block and block_start < current_pos:
            blocks.append((block_start, current_pos - 1))
        
        return blocks
    
    def insert_isotag_structures(self, structures: Dict):
        """Insert isotag structure data"""
        click.echo(f"   📊 Inserting {len(structures)} isotag structures...")
        
        for isotag_id, structure in structures.items():
            exons = structure['exons']
            
            # Insert main isotag record
            self.db.cursor.execute("""
                INSERT OR IGNORE INTO isotags 
                (isotag_id, chrom, strand, start_pos, end_pos, exon_count, total_length)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (
                isotag_id,
                structure['chrom'],
                structure['strand'],
                structure['start'],
                structure['end'],
                len(exons),
                sum(end - start + 1 for start, end in exons)
            ))
            
            # Insert exon records
            for i, (start, end) in enumerate(exons, 1):
                self.db.cursor.execute("""
                    INSERT INTO exons (isotag_id, exon_number, start_pos, end_pos)
                    VALUES (?, ?, ?, ?)
                """, (isotag_id, i, start, end))
        
        self.db.conn.commit()
        self.stats['isotags_indexed'] = len(structures)
        self.stats['exons_indexed'] = sum(len(s['exons']) for s in structures.values())
    
    def insert_sample_data(self, counts: Counter, sample_name: str, bam_file: str, total_isotags: int):
        """Insert sample count data"""
        for isotag_id, read_count in counts.items():
            self.db.cursor.execute("""
                INSERT INTO sample_data (isotag_id, sample_name, read_count, total_reads, bam_file)
                VALUES (?, ?, ?, ?, ?)
            """, (isotag_id, sample_name, read_count, total_isotags, bam_file))
        
        self.db.conn.commit()
    
    def insert_variant_data(self, variant_data: Dict):
        """Insert variant data"""
        click.echo(f"   🧬 Inserting variant data...")
        
        for isotag_id, variant_list in variant_data.items():
            variant_counts = Counter(variant_list)
            
            for variant_id, count in variant_counts.items():
                self.db.cursor.execute("""
                    INSERT INTO variants (isotag_id, variant_id, variant_count)
                    VALUES (?, ?, ?)
                """, (isotag_id, variant_id, count))
        
        self.db.conn.commit()
        self.stats['variants_indexed'] = sum(len(v) for v in variant_data.values())
    
    def query_isotags(self, **kwargs) -> List[Dict]:
        """Query isotags with various filters"""
        conditions = []
        params = []
        
        # Basic filters
        if 'chrom' in kwargs:
            conditions.append("chrom = ?")
            params.append(kwargs['chrom'])
        
        if 'strand' in kwargs:
            conditions.append("strand = ?")
            params.append(kwargs['strand'])
        
        if 'min_exons' in kwargs:
            conditions.append("exon_count >= ?")
            params.append(kwargs['min_exons'])
        
        if 'max_exons' in kwargs:
            conditions.append("exon_count <= ?")
            params.append(kwargs['max_exons'])
        
        if 'min_length' in kwargs:
            conditions.append("total_length >= ?")
            params.append(kwargs['min_length'])
        
        # Position filters
        if 'region' in kwargs:
            # Format: chr:start-end
            region = kwargs['region']
            if ':' in region and '-' in region:
                chrom, pos_range = region.split(':', 1)
                start, end = map(int, pos_range.split('-'))
                conditions.extend([
                    "chrom = ?",
                    "start_pos <= ?",
                    "end_pos >= ?"
                ])
                params.extend([chrom, end, start])
        
        # Build query
        base_query = """
            SELECT isotag_id, chrom, strand, start_pos, end_pos, 
                   exon_count, total_length, created_at
            FROM isotags
        """
        
        if conditions:
            query = base_query + " WHERE " + " AND ".join(conditions)
        else:
            query = base_query
        
        # Add ordering and limit
        query += " ORDER BY chrom, start_pos"
        
        if 'limit' in kwargs:
            query += " LIMIT ?"
            params.append(kwargs['limit'])
        
        self.db.cursor.execute(query, params)
        
        columns = [desc[0] for desc in self.db.cursor.description]
        results = []
        
        for row in self.db.cursor.fetchall():
            result_dict = dict(zip(columns, row))
            results.append(result_dict)
        
        return results
    
    def get_isotag_details(self, isotag_id: str) -> Optional[Dict]:
        """Get detailed information for a specific isotag"""
        # Get basic info
        self.db.cursor.execute("""
            SELECT * FROM isotags WHERE isotag_id = ?
        """, (isotag_id,))
        
        isotag_info = self.db.cursor.fetchone()
        if not isotag_info:
            return None
        
        # Get exons
        self.db.cursor.execute("""
            SELECT exon_number, start_pos, end_pos
            FROM exons WHERE isotag_id = ?
            ORDER BY exon_number
        """, (isotag_id,))
        
        exons = self.db.cursor.fetchall()
        
        # Get sample data
        self.db.cursor.execute("""
            SELECT sample_name, read_count, total_reads, bam_file
            FROM sample_data WHERE isotag_id = ?
        """, (isotag_id,))
        
        samples = self.db.cursor.fetchall()
        
        # Get variants
        self.db.cursor.execute("""
            SELECT variant_id, variant_count
            FROM variants WHERE isotag_id = ?
        """, (isotag_id,))
        
        variants = self.db.cursor.fetchall()
        
        # Get annotations
        self.db.cursor.execute("""
            SELECT annotation_type, gene_id, transcript_id, classification, overlap_score
            FROM annotations WHERE isotag_id = ?
        """, (isotag_id,))
        
        annotations = self.db.cursor.fetchall()
        
        # Compile results
        result = {
            'isotag_info': dict(zip([desc[0] for desc in self.db.cursor.description], isotag_info)) if isotag_info else None,
            'exons': [{'exon': e[0], 'start': e[1], 'end': e[2]} for e in exons],
            'samples': [{'sample': s[0], 'reads': s[1], 'total': s[2], 'bam': s[3]} for s in samples],
            'variants': [{'variant_id': v[0], 'count': v[1]} for v in variants],
            'annotations': [{'type': a[0], 'gene': a[1], 'transcript': a[2], 'class': a[3], 'score': a[4]} for a in annotations]
        }
        
        return result
    
    def get_database_stats(self) -> Dict:
        """Get database statistics"""
        stats = {}
        
        # Count tables
        tables = ['isotags', 'exons', 'sample_data', 'variants', 'annotations']
        for table in tables:
            self.db.cursor.execute(f"SELECT COUNT(*) FROM {table}")
            stats[f'{table}_count'] = self.db.cursor.fetchone()[0]
        
        # Unique samples
        self.db.cursor.execute("SELECT COUNT(DISTINCT sample_name) FROM sample_data")
        stats['unique_samples'] = self.db.cursor.fetchone()[0]
        
        # Chromosome distribution
        self.db.cursor.execute("""
            SELECT chrom, COUNT(*) 
            FROM isotags 
            GROUP BY chrom 
            ORDER BY COUNT(*) DESC 
            LIMIT 10
        """)
        stats['top_chromosomes'] = dict(self.db.cursor.fetchall())
        
        return stats


@click.group()
def index():
    """ISO-Tools Index - Create and Query Isotag Databases"""
    pass


@index.command()
@click.argument('database_path', type=click.Path())
@click.argument('bam_files', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('--samples', help='Comma-separated sample names')
def create(database_path, bam_files, samples):
    """
    Create isotag database from BAM files
    
    Examples:
        # Create database from BAM files
        iso-tools index create isotag.db sample1.bam sample2.bam
        
        # With custom sample names
        iso-tools index create isotag.db *.bam --samples "Control1,Control2,Treatment1,Treatment2"
    """
    
    sample_names = None
    if samples:
        sample_names = [s.strip() for s in samples.split(',')]
    
    indexer = IsotogIndexer()
    
    # Create database
    indexer.create_database(database_path)
    
    # Index BAM files
    indexer.index_bam_files(list(bam_files), sample_names)
    
    # Display final stats
    click.echo("\n" + "="*50)
    click.echo("✅ INDEXING COMPLETE")
    click.echo("="*50)
    click.echo(f"📊 Isotags indexed:     {indexer.stats['isotags_indexed']:,}")
    click.echo(f"📁 Samples processed:   {indexer.stats['samples_processed']:,}")
    click.echo(f"🧬 Exons indexed:       {indexer.stats['exons_indexed']:,}")
    click.echo(f"🔬 Variants indexed:    {indexer.stats['variants_indexed']:,}")
    click.echo(f"💾 Database: {database_path}")
    
    # Show database stats
    db_stats = indexer.get_database_stats()
    click.echo(f"\n📈 Database contains:")
    for key, value in db_stats.items():
        if not key.startswith('top_'):
            click.echo(f"   {key}: {value:,}")
    
    indexer.db.disconnect()


@index.command()
@click.argument('database_path', type=click.Path(exists=True))
@click.option('--chrom', help='Filter by chromosome')
@click.option('--strand', help='Filter by strand (+/-)')
@click.option('--region', help='Filter by region (chr:start-end)')
@click.option('--min-exons', type=int, help='Minimum number of exons')
@click.option('--max-exons', type=int, help='Maximum number of exons')
@click.option('--min-length', type=int, help='Minimum total length')
@click.option('--limit', type=int, default=100, help='Maximum results to return')
@click.option('--output', '-o', help='Output results to file')
def query(database_path, chrom, strand, region, min_exons, max_exons, min_length, limit, output):
    """
    Query isotag database
    
    Examples:
        # Query all isotags on chromosome 1
        iso-tools index query isotag.db --chrom 1 --limit 50
        
        # Query by region
        iso-tools index query isotag.db --region "1:1000000-2000000"
        
        # Filter by exon count
        iso-tools index query isotag.db --min-exons 3 --max-exons 10
        
        # Export results
        iso-tools index query isotag.db --chrom X --output results.json
    """
    
    indexer = IsotogIndexer()
    indexer.db = IsotogDatabase(database_path)
    indexer.db.connect()
    
    # Build query parameters
    query_params = {}
    if chrom:
        query_params['chrom'] = chrom
    if strand:
        query_params['strand'] = strand
    if region:
        query_params['region'] = region
    if min_exons:
        query_params['min_exons'] = min_exons
    if max_exons:
        query_params['max_exons'] = max_exons
    if min_length:
        query_params['min_length'] = min_length
    if limit:
        query_params['limit'] = limit
    
    # Execute query
    results = indexer.query_isotags(**query_params)
    
    if output:
        # Export to file
        with open(output, 'w') as f:
            json.dump(results, f, indent=2)
        click.echo(f"📄 Results exported to {output}")
    else:
        # Display results
        click.echo(f"Found {len(results)} isotags:")
        click.echo("-" * 80)
        
        for result in results:
            click.echo(f"{result['isotag_id'][:32]}... "
                      f"{result['chrom']}:{result['start_pos']}-{result['end_pos']} "
                      f"({result['strand']}) "
                      f"{result['exon_count']} exons, "
                      f"{result['total_length']} bp")
    
    indexer.db.disconnect()


@index.command()
@click.argument('database_path', type=click.Path(exists=True))
@click.argument('isotag_id')
def details(database_path, isotag_id):
    """
    Get detailed information for a specific isotag
    
    Examples:
        iso-tools index details isotag.db Uy3v73FzY4ZhB-w0ZLsDwQLJfMl34Hzx
    """
    
    indexer = IsotogIndexer()
    indexer.db = IsotogDatabase(database_path)
    indexer.db.connect()
    
    details = indexer.get_isotag_details(isotag_id)
    
    if not details:
        click.echo(f"❌ Isotag not found: {isotag_id}")
        sys.exit(1)
    
    # Display detailed information
    info = details['isotag_info']
    click.echo(f"🆔 Isotag ID: {info['isotag_id']}")
    click.echo(f"📍 Location: {info['chrom']}:{info['start_pos']}-{info['end_pos']} ({info['strand']})")
    click.echo(f"🧬 Structure: {info['exon_count']} exons, {info['total_length']} bp")
    
    if details['exons']:
        click.echo(f"\n📋 EXONS")
        for exon in details['exons']:
            click.echo(f"   Exon {exon['exon']}: {exon['start']}-{exon['end']}")
    
    if details['samples']:
        click.echo(f"\n📊 SAMPLE DATA")
        for sample in details['samples']:
            click.echo(f"   {sample['sample']}: {sample['reads']} reads")
    
    if details['variants']:
        click.echo(f"\n🔬 VARIANTS")
        for variant in details['variants'][:5]:  # Show first 5
            click.echo(f"   {variant['variant_id'][:32]}...: {variant['count']} reads")
    
    if details['annotations']:
        click.echo(f"\n🔗 ANNOTATIONS")
        for ann in details['annotations']:
            click.echo(f"   {ann['type']}: {ann['gene']} ({ann['class']})")
    
    indexer.db.disconnect()


if __name__ == '__main__':
    index()