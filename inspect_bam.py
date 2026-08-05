#!/usr/bin/env python3
"""
BAM file inspection utility using samtools
Memory-efficient approach without pysam
"""

import subprocess
import click
import json
import sys
from collections import defaultdict

VERSION = "2.0.0"


def parse_sam_line(sam_line: str):
    """Parse a single SAM format line"""
    if sam_line.startswith('@') or not sam_line.strip():
        return None
    
    fields = sam_line.strip().split('\t')
    if len(fields) < 11:
        return None
    
    # Parse basic fields
    qname = fields[0]
    flag = int(fields[1])
    rname = fields[2]
    pos = int(fields[3])
    cigar = fields[5]
    seq = fields[9]
    
    # Parse optional fields for MD tag
    md_tag = None
    for field in fields[11:]:
        if field.startswith('MD:Z:'):
            md_tag = field[5:]
            break
    
    return {
        'qname': qname,
        'flag': flag,
        'rname': rname,
        'pos': pos,
        'cigar': cigar,
        'seq': seq,
        'md_tag': md_tag,
        'is_unmapped': flag & 0x4,
        'is_reverse': flag & 0x10
    }


@click.command()
@click.option('--bam', help='BAM file to inspect')
@click.option('--limit', default=10, help='Number of reads to show in detail')
@click.option('--stats-limit', default=10000, help='Number of reads to analyze for statistics')
@click.option('--version', is_flag=True, help='Show version and exit')
def inspect_bam(bam, limit, stats_limit, version):
    """Inspect BAM file structure and MD tags using streaming samtools"""

    if version:
        click.echo(f"inspect_bam.py v{VERSION}")
        sys.exit(0)

    if not bam:
        click.echo("Error: --bam is required", err=True)
        sys.exit(1)

    click.echo(f"Inspecting BAM file: {bam}")
    click.echo(f"Will show detailed info for first {limit} reads")
    click.echo(f"Will analyze {stats_limit} reads for statistics")
    
    # First, get header information
    try:
        # Get BAM header
        header_cmd = ['samtools', 'view', '-H', bam]
        header_result = subprocess.run(header_cmd, capture_output=True, text=True, check=True)
        
        click.echo("\n=== BAM Header Info ===")
        reference_count = 0
        for line in header_result.stdout.split('\n'):
            if line.startswith('@SQ'):
                reference_count += 1
                if reference_count <= 3:  # Show first 3 references
                    parts = line.split('\t')
                    seq_name = next((p[3:] for p in parts if p.startswith('SN:')), 'Unknown')
                    seq_len = next((p[3:] for p in parts if p.startswith('LN:')), 'Unknown')
                    click.echo(f"  Reference {reference_count}: {seq_name} (length: {seq_len})")
        
        if reference_count > 3:
            click.echo(f"  ... and {reference_count - 3} more references")
        
        click.echo(f"Total references: {reference_count}")
        
    except subprocess.CalledProcessError as e:
        click.echo(f"Error reading BAM header: {e}")
        return
    
    # Stream through reads for detailed inspection and statistics
    try:
        # Use samtools view to stream reads
        cmd = ['samtools', 'view', bam]
        
        stats = defaultdict(int)
        chromosomes = defaultdict(int)
        cigar_patterns = defaultdict(int)
        
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
            click.echo(f"\n=== Sample Reads (first {limit}) ===")
            
            for line_num, line in enumerate(process.stdout):
                read_info = parse_sam_line(line)
                if not read_info:
                    continue
                
                # Collect statistics
                stats['total_reads'] += 1
                
                if not read_info['is_unmapped']:
                    stats['mapped_reads'] += 1
                    chromosomes[read_info['rname']] += 1
                
                if read_info['md_tag']:
                    stats['reads_with_md'] += 1
                
                # Check for splicing (N in CIGAR)
                if 'N' in read_info['cigar']:
                    stats['reads_with_splicing'] += 1
                
                # Count CIGAR patterns
                if read_info['cigar'] != '*':
                    # Simplify CIGAR for pattern counting
                    simplified_cigar = ''.join(c for c in read_info['cigar'] if c.isalpha())
                    cigar_patterns[simplified_cigar] += 1
                
                # Show detailed info for first few reads
                if line_num < limit:
                    click.echo(f"\nRead {line_num + 1}: {read_info['qname']}")
                    click.echo(f"  Reference: {read_info['rname']}")
                    click.echo(f"  Position: {read_info['pos']}")
                    click.echo(f"  Strand: {'-' if read_info['is_reverse'] else '+'}")
                    click.echo(f"  CIGAR: {read_info['cigar']}")
                    click.echo(f"  MD tag: {read_info['md_tag'] or 'None'}")
                    click.echo(f"  Mapped: {not read_info['is_unmapped']}")
                    click.echo(f"  Sequence length: {len(read_info['seq'])}")
                
                # Stop after stats_limit for statistics
                if stats['total_reads'] >= stats_limit:
                    break
        
        # Display statistics
        click.echo(f"\n=== Summary Statistics (from {stats['total_reads']} reads) ===")
        click.echo(f"Total reads analyzed: {stats['total_reads']}")
        click.echo(f"Mapped reads: {stats['mapped_reads']} ({100*stats['mapped_reads']/stats['total_reads']:.1f}%)")
        click.echo(f"Reads with MD tags: {stats['reads_with_md']} ({100*stats['reads_with_md']/stats['total_reads']:.1f}%)")
        click.echo(f"Reads with splicing: {stats['reads_with_splicing']} ({100*stats['reads_with_splicing']/stats['total_reads']:.1f}%)")
        
        # Top chromosomes
        if chromosomes:
            click.echo(f"\n=== Top Chromosomes ===")
            sorted_chroms = sorted(chromosomes.items(), key=lambda x: x[1], reverse=True)
            for chrom, count in sorted_chroms[:5]:
                click.echo(f"  {chrom}: {count} reads")
        
        # Common CIGAR patterns
        if cigar_patterns:
            click.echo(f"\n=== Common CIGAR Patterns ===")
            sorted_cigars = sorted(cigar_patterns.items(), key=lambda x: x[1], reverse=True)
            for cigar, count in sorted_cigars[:5]:
                click.echo(f"  {cigar}: {count} reads")
        
        # Warnings
        click.echo(f"\n=== Analysis ===")
        if stats['reads_with_md'] == 0:
            click.echo("⚠️  WARNING: No MD tags found!")
            click.echo("   You need to run 'samtools calmd' to add MD tags for variant detection.")
        else:
            click.echo("✅ MD tags found - variant detection will work")
        
        if stats['reads_with_splicing'] > 0:
            click.echo(f"✅ Found spliced reads ({stats['reads_with_splicing']}) - this appears to be RNA-seq data")
        else:
            click.echo("ℹ️  No spliced reads found - might be DNA-seq or unspliced RNA")
        
        if stats['mapped_reads'] / stats['total_reads'] < 0.5:
            click.echo("⚠️  WARNING: Low mapping rate - check reference genome compatibility")
    
    except subprocess.CalledProcessError as e:
        click.echo(f"Error reading BAM file: {e}")
    except Exception as e:
        click.echo(f"Unexpected error: {e}")


if __name__ == "__main__":
    inspect_bam()