#!/usr/bin/env python3
"""
IsoTag BED Index Creator v1.0

Creates compact reference index from tagged BAM files for fast bedtools-based intersection.
Groups reads by clustering tag (XC, XT, or XI) and outputs genomic regions with tag lists.

This solves the storage bloat problem in isotag_intersect.py by creating a compact index
(~1-5 MB) instead of storing all file offsets (~500-1000 MB).

Workflow:
    1. Create compact reference: isotag_create_bed_index.py -i ref.bam -o ref.bed
    2. Fast intersection: bedtools intersect -a exp.bam -b ref.bed -wa -wb > overlaps.txt
    3. Check tag matches: isotag_query_bed_index.py -i overlaps.txt -o matches.tsv

Output Formats:
    - BED: Standard BED format with metadata in name field
    - GTF: GTF format with attributes for tag lists
    - TSV: Tab-delimited with dedicated columns (most flexible)

Usage:
    # Group by XC (gene/locus level), include XI and XS tags
    python3 isotag_create_bed_index.py -i ref.bam -o ref.bed --group-by XC --include-tags XI,XS

    # GTF output with all tags
    python3 isotag_create_bed_index.py -i ref.bam -o ref.gtf --format gtf --include-tags XI,XS,XT

    # TSV with minimum 5 reads per cluster
    python3 isotag_create_bed_index.py -i ref.bam -o ref.tsv --format tsv --min-reads 5

Author: RIKEN IMS - Laboratory for Large-Scale Biomedical Data Technology
Version: 1.0
Date: 2026-01-27
"""

import sys
import click
import pysam
from collections import defaultdict
from typing import Dict, List, Set, Tuple
import os

__version__ = "1.0"


class ClusterGroup:
    """Represents a group of reads with the same clustering tag"""

    def __init__(self, chrom: str, strand: str, cluster_id: str):
        self.chrom = chrom
        self.strand = strand
        self.cluster_id = cluster_id
        self.start = float('inf')
        self.end = 0
        self.xi_tags = set()
        self.xb_tags = set()
        self.xs_tags = set()
        self.xt_tags = set()
        self.xv_tags = set()
        self.read_count = 0

    def add_read(self, read, include_tags: Set[str]):
        """Add a read to this cluster group"""
        self.read_count += 1

        # Update genomic span
        if read.reference_start < self.start:
            self.start = read.reference_start
        if read.reference_end > self.end:
            self.end = read.reference_end

        # Collect tags
        if 'XI' in include_tags and read.has_tag('XI'):
            self.xi_tags.add(read.get_tag('XI'))
        if 'XB' in include_tags and read.has_tag('XB'):
            self.xb_tags.add(read.get_tag('XB'))
        if 'XS' in include_tags and read.has_tag('XS'):
            self.xs_tags.add(read.get_tag('XS'))
        if 'XT' in include_tags and read.has_tag('XT'):
            self.xt_tags.add(read.get_tag('XT'))
        if 'XV' in include_tags and read.has_tag('XV'):
            self.xv_tags.add(read.get_tag('XV'))

    def to_bed(self, name_format: str = "compact") -> str:
        """Convert to BED format"""
        # BED format: chrom start end name score strand
        if name_format == "compact":
            # Compact: cluster_id;n=count;XI=count;XS=count
            name = f"{self.cluster_id[:8]};n={self.read_count}"
            if self.xi_tags:
                name += f";XI={len(self.xi_tags)}"
            if self.xs_tags:
                name += f";XS={len(self.xs_tags)}"
            if self.xt_tags:
                name += f";XT={len(self.xt_tags)}"
        else:
            # Full: include actual tag lists (may be long)
            parts = [self.cluster_id[:8]]
            if self.xi_tags:
                parts.append(f"XI={','.join(sorted(self.xi_tags)[:5])}")  # Limit to 5
            parts.append(f"n={self.read_count}")
            name = ";".join(parts)

        return f"{self.chrom}\t{self.start}\t{self.end}\t{name}\t.\t{self.strand}"

    def to_gtf(self) -> str:
        """Convert to GTF format"""
        # GTF: seqname source feature start end score strand frame attributes
        attributes = [
            f'cluster_id "{self.cluster_id}"',
            f'read_count "{self.read_count}"'
        ]

        if self.xi_tags:
            xi_list = ','.join(sorted(self.xi_tags))
            attributes.append(f'XI "{xi_list}"')
            attributes.append(f'XI_count "{len(self.xi_tags)}"')

        if self.xs_tags:
            xs_list = ','.join(sorted(self.xs_tags))
            attributes.append(f'XS "{xs_list}"')
            attributes.append(f'XS_count "{len(self.xs_tags)}"')

        if self.xt_tags:
            xt_list = ','.join(sorted(self.xt_tags))
            attributes.append(f'XT "{xt_list}"')
            attributes.append(f'XT_count "{len(self.xt_tags)}"')

        if self.xb_tags:
            attributes.append(f'XB_count "{len(self.xb_tags)}"')

        if self.xv_tags:
            attributes.append(f'XV_count "{len(self.xv_tags)}"')

        attr_str = "; ".join(attributes) + ";"

        # GTF is 1-based, but our coordinates are 0-based from BAM
        gtf_start = self.start + 1

        return f"{self.chrom}\tisotag\tcluster\t{gtf_start}\t{self.end}\t.\t{self.strand}\t.\t{attr_str}"

    def to_tsv(self) -> str:
        """Convert to TSV format"""
        xi_list = ','.join(sorted(self.xi_tags)) if self.xi_tags else '.'
        xb_list = ','.join(sorted(self.xb_tags)) if self.xb_tags else '.'
        xs_list = ','.join(sorted(self.xs_tags)) if self.xs_tags else '.'
        xt_list = ','.join(sorted(self.xt_tags)) if self.xt_tags else '.'
        xv_list = ','.join(sorted(self.xv_tags)) if self.xv_tags else '.'

        return (f"{self.chrom}\t{self.start}\t{self.end}\t{self.strand}\t"
                f"{self.cluster_id}\t{self.read_count}\t"
                f"{len(self.xi_tags)}\t{xi_list}\t"
                f"{len(self.xs_tags)}\t{xs_list}\t"
                f"{len(self.xt_tags)}\t{xt_list}\t"
                f"{len(self.xb_tags)}\t{xb_list}\t"
                f"{len(self.xv_tags)}\t{xv_list}")


def parse_include_tags(include_tags: str) -> Set[str]:
    """Parse comma-separated tag list"""
    if not include_tags:
        return {'XI', 'XS', 'XT'}  # Default tags

    valid_tags = {'XI', 'XB', 'XS', 'XT', 'XV'}
    tags = set(tag.strip().upper() for tag in include_tags.split(','))

    invalid = tags - valid_tags
    if invalid:
        click.echo(f"Warning: Invalid tags ignored: {invalid}", err=True)

    return tags & valid_tags


def get_strand(read) -> str:
    """Get strand from read"""
    return '-' if read.is_reverse else '+'


@click.command()
@click.option('-i', '--input', 'input_bam', required=True, type=click.Path(exists=True),
              help='Input BAM file with isotag tags (XI, XB, XS, XT, XC)')
@click.option('-o', '--output', required=True, type=click.Path(),
              help='Output file (.bed, .gtf, or .tsv)')
@click.option('--group-by', default='XC', type=click.Choice(['XC', 'XT', 'XI'], case_sensitive=False),
              help='Tag to group reads by (default: XC for gene/locus level)')
@click.option('--include-tags', default='XI,XS,XT',
              help='Comma-separated list of tags to include (default: XI,XS,XT)')
@click.option('--format', 'output_format', type=click.Choice(['bed', 'gtf', 'tsv'], case_sensitive=False),
              help='Output format (auto-detected from extension if not specified)')
@click.option('--min-reads', default=1, type=int,
              help='Minimum reads per cluster to include (default: 1)')
@click.option('--name-format', type=click.Choice(['compact', 'full'], case_sensitive=False), default='compact',
              help='BED name field format: compact (counts only) or full (tag lists)')
@click.version_option(version=__version__)
def main(input_bam, output, group_by, include_tags, output_format, min_reads, name_format):
    """
    Create compact BED/GTF/TSV index from tagged BAM file.

    This tool groups reads by a clustering tag (XC, XT, or XI) and creates a compact
    reference file with genomic coordinates and tag lists. This enables fast bedtools
    intersection followed by precise tag matching.

    Storage savings: ~100-1000× smaller than isotag_intersect.py database!

    Examples:

        # Create BED index grouped by XC (gene/locus level)
        isotag_create_bed_index.py -i reference.bam -o ref.bed --group-by XC

        # Create GTF with all tags, minimum 5 reads per cluster
        isotag_create_bed_index.py -i ref.bam -o ref.gtf --format gtf --min-reads 5

        # TSV format for maximum flexibility
        isotag_create_bed_index.py -i ref.bam -o ref.tsv --format tsv
    """

    # Auto-detect format from extension if not specified
    if not output_format:
        ext = os.path.splitext(output)[1].lower()
        if ext == '.bed':
            output_format = 'bed'
        elif ext == '.gtf':
            output_format = 'gtf'
        elif ext in ['.tsv', '.txt']:
            output_format = 'tsv'
        else:
            click.echo("Error: Cannot auto-detect format. Please specify --format", err=True)
            sys.exit(1)

    output_format = output_format.lower()
    group_by = group_by.upper()

    # Parse included tags
    include_tag_set = parse_include_tags(include_tags)

    click.echo(f"Reading BAM file: {input_bam}")
    click.echo(f"Grouping by: {group_by}")
    click.echo(f"Including tags: {', '.join(sorted(include_tag_set))}")
    click.echo(f"Output format: {output_format}")
    click.echo(f"Minimum reads: {min_reads}")

    # Open BAM file
    try:
        bam = pysam.AlignmentFile(input_bam, 'rb')
    except Exception as e:
        click.echo(f"Error opening BAM file: {e}", err=True)
        sys.exit(1)

    # Group reads by clustering tag
    clusters: Dict[Tuple[str, str, str], ClusterGroup] = {}
    skipped_no_tag = 0
    skipped_unmapped = 0
    reads_processed = 0

    click.echo("\nProcessing reads...")

    for read in bam:
        reads_processed += 1

        if reads_processed % 10000 == 0:
            click.echo(f"  Processed {reads_processed:,} reads, found {len(clusters):,} clusters...", nl=False)
            click.echo('\r', nl=False)

        # Skip unmapped
        if read.is_unmapped:
            skipped_unmapped += 1
            continue

        # Check for clustering tag
        if not read.has_tag(group_by):
            skipped_no_tag += 1
            continue

        cluster_id = read.get_tag(group_by)
        chrom = read.reference_name
        strand = get_strand(read)

        # Create cluster group key
        key = (chrom, strand, cluster_id)

        # Create or update cluster
        if key not in clusters:
            clusters[key] = ClusterGroup(chrom, strand, cluster_id)

        clusters[key].add_read(read, include_tag_set)

    bam.close()

    click.echo(f"\nProcessed {reads_processed:,} reads")
    click.echo(f"  Skipped {skipped_unmapped:,} unmapped reads")
    click.echo(f"  Skipped {skipped_no_tag:,} reads without {group_by} tag")
    click.echo(f"  Found {len(clusters):,} unique {group_by} clusters")

    # Filter by minimum reads
    filtered_clusters = [c for c in clusters.values() if c.read_count >= min_reads]

    if min_reads > 1:
        click.echo(f"  {len(filtered_clusters):,} clusters pass min-reads filter (>= {min_reads})")

    # Write output
    click.echo(f"\nWriting {output_format.upper()} output: {output}")

    with open(output, 'w') as out:
        # Write header for TSV
        if output_format == 'tsv':
            out.write("chrom\tstart\tend\tstrand\tcluster_id\tread_count\t")
            out.write("XI_count\tXI_tags\tXS_count\tXS_tags\tXT_count\tXT_tags\t")
            out.write("XB_count\tXB_tags\tXV_count\tXV_tags\n")

        # Sort clusters by chromosome and position
        sorted_clusters = sorted(filtered_clusters, key=lambda c: (c.chrom, c.start))

        for cluster in sorted_clusters:
            if output_format == 'bed':
                out.write(cluster.to_bed(name_format) + '\n')
            elif output_format == 'gtf':
                out.write(cluster.to_gtf() + '\n')
            elif output_format == 'tsv':
                out.write(cluster.to_tsv() + '\n')

    # Statistics
    total_reads = sum(c.read_count for c in filtered_clusters)
    avg_reads = total_reads / len(filtered_clusters) if filtered_clusters else 0

    click.echo(f"\n✓ Success!")
    click.echo(f"  Clusters written: {len(filtered_clusters):,}")
    click.echo(f"  Total reads represented: {total_reads:,}")
    click.echo(f"  Average reads per cluster: {avg_reads:.1f}")

    # Tag diversity statistics
    if include_tag_set:
        click.echo("\n  Tag diversity:")
        if 'XI' in include_tag_set:
            total_xi = sum(len(c.xi_tags) for c in filtered_clusters)
            avg_xi = total_xi / len(filtered_clusters) if filtered_clusters else 0
            click.echo(f"    XI tags: {total_xi:,} total, {avg_xi:.1f} avg per cluster")
        if 'XS' in include_tag_set:
            total_xs = sum(len(c.xs_tags) for c in filtered_clusters)
            avg_xs = total_xs / len(filtered_clusters) if filtered_clusters else 0
            click.echo(f"    XS tags: {total_xs:,} total, {avg_xs:.1f} avg per cluster")
        if 'XT' in include_tag_set:
            total_xt = sum(len(c.xt_tags) for c in filtered_clusters)
            avg_xt = total_xt / len(filtered_clusters) if filtered_clusters else 0
            click.echo(f"    XT tags: {total_xt:,} total, {avg_xt:.1f} avg per cluster")

    click.echo(f"\nNext step: Use bedtools intersect with this reference file")
    click.echo(f"  bedtools intersect -a experiment.bam -b {output} -wa -wb -s > overlaps.txt")


if __name__ == '__main__':
    main()
