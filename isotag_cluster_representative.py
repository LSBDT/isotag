#!/usr/bin/env python3
"""
isotag_cluster_representative - Select Representative Reads from XC Clusters

Selects one representative read per XC cluster based on user-specified criteria.
This is useful for:
- Reducing redundancy in datasets
- Creating non-redundant isoform catalogs
- Picking the "best" example from each cluster

Available selection criteria:
- mapping_quality: Highest MAPQ (default)
- read_length: Longest read
- first: First encountered read (fastest)

Usage:
    # Select by mapping quality (default)
    python3 isotag_cluster_representative.py -i tagged.bam -o representatives.bam

    # Select by read length
    python3 isotag_cluster_representative.py -i tagged.bam -o representatives.bam --by read_length

    # Output read IDs only (TSV)
    python3 isotag_cluster_representative.py -i tagged.bam -o representatives.tsv --format tsv

Author: IsoTag Team
Version: 1.0.0
"""

import pysam
import click
import sys
from collections import defaultdict
from typing import Dict, List, Tuple, Optional


class ClusterRepresentativeSelector:
    """Selects representative reads from XC clusters"""

    def __init__(self, selection_method: str = "mapping_quality"):
        """
        Initialize selector

        Args:
            selection_method: How to pick representative
                - "mapping_quality": Highest MAPQ
                - "read_length": Longest read
                - "first": First encountered
        """
        self.selection_method = selection_method

    def collect_clusters(self, bam_path: str) -> Dict[str, List[Tuple]]:
        """
        Collect reads grouped by XC cluster

        Returns:
            Dict of {xc_tag: [(read_name, mapq, read_length, file_offset), ...]}
        """
        clusters = defaultdict(list)

        try:
            bam = pysam.AlignmentFile(bam_path, 'rb')
        except Exception as e:
            click.echo(f"Error opening BAM file: {e}", err=True)
            sys.exit(1)

        read_count = 0
        xc_count = 0

        for read in bam.fetch(until_eof=True):
            read_count += 1

            if read_count % 100000 == 0:
                click.echo(f"  Processed {read_count:,} reads...", err=True)

            # Get XC tag
            try:
                xc_tag = read.get_tag('XC')
                xc_count += 1
            except KeyError:
                continue

            # Store read info for selection
            read_info = (
                read.query_name,
                read.mapping_quality,
                read.query_length or 0,
                bam.tell()  # File offset for extraction
            )
            clusters[xc_tag].append(read_info)

        bam.close()

        click.echo(f"  Total reads: {read_count:,}", err=True)
        click.echo(f"  Reads with XC: {xc_count:,}", err=True)
        click.echo(f"  Unique XC clusters: {len(clusters):,}", err=True)

        return clusters

    def select_representatives(self, clusters: Dict[str, List[Tuple]]) -> Dict[str, Tuple]:
        """
        Select one representative per cluster

        Returns:
            Dict of {xc_tag: (read_name, mapq, read_length, file_offset)}
        """
        representatives = {}

        for xc_tag, reads in clusters.items():
            if self.selection_method == "mapping_quality":
                # Sort by MAPQ (descending), then read length (descending), then name (for determinism)
                best = max(reads, key=lambda x: (x[1], x[2], x[0]))
            elif self.selection_method == "read_length":
                # Sort by read length (descending), then MAPQ (descending), then name
                best = max(reads, key=lambda x: (x[2], x[1], x[0]))
            else:  # first
                best = reads[0]

            representatives[xc_tag] = best

        return representatives


def write_bam_output(input_bam: str, output_bam: str, representatives: Dict[str, Tuple]):
    """Write representative reads to output BAM"""

    # Get file offsets of representatives
    target_offsets = {rep[3] for rep in representatives.values()}

    click.echo(f"  Extracting {len(target_offsets):,} representative reads...", err=True)

    bam_in = pysam.AlignmentFile(input_bam, 'rb')
    bam_out = pysam.AlignmentFile(output_bam, 'wb', template=bam_in)

    written = 0
    for read in bam_in.fetch(until_eof=True):
        offset = bam_in.tell()
        if offset in target_offsets:
            bam_out.write(read)
            written += 1

            if written % 10000 == 0:
                click.echo(f"  Written {written:,} reads...", err=True)

    bam_in.close()
    bam_out.close()

    # Index output BAM
    click.echo(f"  Indexing output BAM...", err=True)
    pysam.index(output_bam)

    return written


def write_tsv_output(output_tsv: str, representatives: Dict[str, Tuple]):
    """Write representative read IDs to TSV"""

    with open(output_tsv, 'w') as f:
        # Header
        f.write("xc_cluster\tread_name\tmapping_quality\tread_length\n")

        # Sort by XC tag for consistent output
        for xc_tag in sorted(representatives.keys()):
            read_name, mapq, read_length, _ = representatives[xc_tag]
            f.write(f"{xc_tag}\t{read_name}\t{mapq}\t{read_length}\n")

    return len(representatives)


@click.command()
@click.option('-i', '--input', 'input_bam', required=True,
              type=click.Path(exists=True),
              help='Input BAM file with XC tags')
@click.option('-o', '--output', 'output_file', required=True,
              type=click.Path(),
              help='Output file (BAM or TSV)')
@click.option('--by', 'selection_method',
              type=click.Choice(['mapping_quality', 'read_length', 'first']),
              default='mapping_quality',
              help='Selection criterion (default: mapping_quality)')
@click.option('--format', 'output_format',
              type=click.Choice(['bam', 'tsv', 'auto']),
              default='auto',
              help='Output format (default: auto-detect from extension)')
@click.option('--stats', 'stats_file',
              type=click.Path(),
              help='Write cluster statistics to this file')
def main(input_bam, output_file, selection_method, output_format, stats_file):
    """
    Select Representative Reads from XC Clusters

    Picks one "best" read per XC cluster based on the specified criterion.
    This reduces redundancy while preserving cluster diversity.

    Examples:

        # Select by mapping quality (default)
        isotag_cluster_representative.py -i tagged.bam -o reps.bam

        # Select by read length
        isotag_cluster_representative.py -i tagged.bam -o reps.bam --by read_length

        # Output as TSV (read IDs only)
        isotag_cluster_representative.py -i tagged.bam -o reps.tsv --format tsv

    The XC tag must be present in the input BAM (added by isotag.py v9.0+).
    """

    click.echo(f"XC Cluster Representative Selection")
    click.echo(f"=" * 60)
    click.echo(f"Input:  {input_bam}")
    click.echo(f"Output: {output_file}")
    click.echo(f"Method: {selection_method}")
    click.echo(f"=" * 60)

    # Determine output format
    if output_format == 'auto':
        if output_file.endswith('.bam'):
            output_format = 'bam'
        elif output_file.endswith('.tsv') or output_file.endswith('.txt'):
            output_format = 'tsv'
        else:
            output_format = 'bam'

    click.echo(f"Format: {output_format}")
    click.echo()

    # Collect clusters
    click.echo("Phase 1: Collecting XC clusters...")
    selector = ClusterRepresentativeSelector(selection_method)
    clusters = selector.collect_clusters(input_bam)

    if not clusters:
        click.echo("Error: No XC tags found in input BAM", err=True)
        click.echo("  Make sure input was tagged with isotag.py v9.0+", err=True)
        sys.exit(1)

    # Select representatives
    click.echo()
    click.echo("Phase 2: Selecting representatives...")
    representatives = selector.select_representatives(clusters)
    click.echo(f"  Selected {len(representatives):,} representatives")

    # Write output
    click.echo()
    click.echo("Phase 3: Writing output...")

    if output_format == 'bam':
        written = write_bam_output(input_bam, output_file, representatives)
    else:
        written = write_tsv_output(output_file, representatives)

    click.echo(f"  Written {written:,} records")

    # Optional: Write statistics
    if stats_file:
        click.echo()
        click.echo(f"Writing cluster statistics to {stats_file}...")

        with open(stats_file, 'w') as f:
            f.write("xc_cluster\tcluster_size\trep_read_name\trep_mapq\trep_length\n")

            for xc_tag in sorted(clusters.keys()):
                cluster_size = len(clusters[xc_tag])
                read_name, mapq, read_length, _ = representatives[xc_tag]
                f.write(f"{xc_tag}\t{cluster_size}\t{read_name}\t{mapq}\t{read_length}\n")

    # Summary
    click.echo()
    click.echo("=" * 60)
    click.echo("Summary:")
    click.echo(f"  Input clusters: {len(clusters):,}")
    click.echo(f"  Representatives: {len(representatives):,}")

    # Calculate compression
    total_reads = sum(len(reads) for reads in clusters.values())
    compression = (1 - len(representatives) / total_reads) * 100 if total_reads > 0 else 0
    click.echo(f"  Compression: {compression:.1f}% ({total_reads:,} -> {len(representatives):,})")

    click.echo(f"  Output: {output_file}")
    click.echo("=" * 60)


if __name__ == '__main__':
    main()
