#!/usr/bin/env python3
"""
IsoTag Stats - Comprehensive Multi-Tag Statistics and Analysis

Analyzes all IsoTag tags (XI, XB, XS, XT, XV, XC, XR) with:
- Tag presence statistics
- Distribution analysis
- Clustering metrics
- Tag co-occurrence analysis
- Export capabilities (JSON, TSV)

Version: v2.0 (Multi-tag support)
"""

import subprocess
import click
import json
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

# Import shared tag definitions
try:
    from tag_definitions import (
        TAG_TYPES, TagExtractor, TagValidator, TagStatistics,
        get_all_tag_names, format_tag_for_display
    )
except ImportError:
    click.echo("❌ Error: tag_definitions.py not found. Ensure it's in the same directory.")
    sys.exit(1)


class MultiTagStatsAnalyzer:
    """Comprehensive multi-tag statistics analyzer"""

    def __init__(self):
        self.stats = {
            'total_reads': 0,
            'version': None,
            'tag_counts': {},
            'unique_counts': {},
            'complexity_metrics': {},
            'clustering_metrics': {},
            'tag_co_occurrence': defaultdict(lambda: defaultdict(int))
        }

        # Counters for each tag type
        self.tag_counters = {tag: Counter() for tag in get_all_tag_names()}

        # Tag statistics helper
        self.tag_statistics = TagStatistics()

        # Read-level tracking
        self.structure_read_counts = defaultdict(list)
        self.variant_read_counts = defaultdict(list)

    def analyze_bam(self, bam_file, limit=None):
        """Analyze all tags in BAM file"""
        click.echo(f"🔍 Analyzing IsoTag tags in {bam_file}")

        view_cmd = ['samtools', 'view', bam_file]

        try:
            process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

            for line_num, line in enumerate(process.stdout):
                if limit and line_num >= limit:
                    break

                self.stats['total_reads'] += 1

                if self.stats['total_reads'] % 10000 == 0:
                    click.echo(f"⏳ Processed {self.stats['total_reads']:,} reads...")

                fields = line.strip().split('\t')

                if len(fields) >= 11:
                    qname = fields[0]

                    # Extract all tags
                    tags = TagExtractor.extract_all_tags(fields)

                    # Update statistics
                    self.tag_statistics.update(tags)

                    # Count each tag type
                    for tag_name, tag_value in tags.items():
                        self.tag_counters[tag_name][tag_value] += 1

                        # Track reads for structure and variants
                        if tag_name == 'XI':
                            self.structure_read_counts[tag_value].append(qname)
                        elif tag_name == 'XV':
                            self.variant_read_counts[tag_value].append(qname)

                    # Tag co-occurrence analysis
                    present_tags = list(tags.keys())
                    for i, tag1 in enumerate(present_tags):
                        for tag2 in present_tags[i+1:]:
                            key = f"{tag1}+{tag2}"
                            self.stats['tag_co_occurrence'][tag1][tag2] += 1

            process.wait()

        except subprocess.CalledProcessError as e:
            click.echo(f"❌ Error reading BAM file: {e}")
            sys.exit(1)

    def compute_statistics(self):
        """Compute comprehensive statistics"""

        # Get tag presence statistics
        tag_summary = self.tag_statistics.get_summary()
        self.stats['version'] = tag_summary.get('version', 'unknown')
        self.stats['tag_counts'] = {
            tag: tag_summary['tag_presence'][tag]['count']
            for tag in get_all_tag_names()
        }

        # Unique counts
        self.stats['unique_counts'] = {
            tag: len(self.tag_counters[tag])
            for tag in get_all_tag_names()
        }

        # Complexity metrics
        if self.stats['total_reads'] > 0:
            self.stats['complexity_metrics'] = {
                'structure_diversity': self.stats['unique_counts']['XI'] / self.stats['total_reads'],
                'splicetag_diversity': self.stats['unique_counts']['XS'] / self.stats['tag_counts'].get('XS', 1),
                'transcript_diversity': self.stats['unique_counts']['XT'] / self.stats['tag_counts'].get('XT', 1),
                'variant_diversity': self.stats['unique_counts']['XV'] / self.stats['tag_counts'].get('XV', 1) if self.stats['tag_counts'].get('XV', 0) > 0 else 0,
                'tagging_efficiency': self.stats['tag_counts']['XI'] / self.stats['total_reads'],
                'multi_exon_rate': self.stats['tag_counts'].get('XS', 0) / self.stats['total_reads'],
                'variant_rate': self.stats['tag_counts'].get('XV', 0) / self.stats['total_reads'],
                'clustering_rate': self.stats['tag_counts'].get('XC', 0) / self.stats['total_reads']
            }

        # Clustering metrics (if clustered)
        if self.stats['tag_counts'].get('XC', 0) > 0:
            # Calculate clustering efficiency
            clustered_reads = self.stats['tag_counts']['XC']
            unique_clusters = self.stats['unique_counts']['XC']

            # Determine what was clustered (XI or XS based)
            if self.stats['tag_counts'].get('XS', 0) > 0:
                # XS-based clustering
                unique_before = self.stats['unique_counts']['XS']
                basis = 'XS (splicetag)'
            else:
                # XI-based clustering
                unique_before = self.stats['unique_counts']['XI']
                basis = 'XI (structure)'

            compression = ((unique_before - unique_clusters) / unique_before * 100) if unique_before > 0 else 0

            self.stats['clustering_metrics'] = {
                'clustered_reads': clustered_reads,
                'unique_clusters': unique_clusters,
                'avg_cluster_size': clustered_reads / unique_clusters if unique_clusters > 0 else 0,
                'compression_ratio': compression,
                'clustering_basis': basis,
                'unique_before_clustering': unique_before
            }

    def display_summary(self):
        """Display comprehensive summary statistics"""
        click.echo("\n" + "="*70)
        click.echo("📊 ISOTAG MULTI-TAG STATISTICS")
        click.echo("="*70)

        click.echo(f"📋 Total reads:  {self.stats['total_reads']:,}")
        click.echo(f"🏷️  Tag version:  {self.stats['version']}")

        # Tag presence summary
        click.echo(f"\n📊 TAG PRESENCE SUMMARY")
        click.echo("-" * 70)

        for tag_name in ['XI', 'XB', 'XS', 'XT', 'XV', 'XC', 'XR']:
            count = self.stats['tag_counts'].get(tag_name, 0)
            unique = self.stats['unique_counts'].get(tag_name, 0)
            pct = 100 * count / self.stats['total_reads'] if self.stats['total_reads'] > 0 else 0
            tag_desc = TAG_TYPES[tag_name]['name']

            status = "✅" if count > 0 else "  "
            click.echo(f"{status} {tag_name} ({tag_desc:20s}): {count:8,} reads ({pct:5.1f}%), {unique:6,} unique")

        # Complexity metrics
        if self.stats['complexity_metrics']:
            metrics = self.stats['complexity_metrics']
            click.echo(f"\n📈 COMPLEXITY METRICS")
            click.echo("-" * 70)
            click.echo(f"Structure diversity:     {metrics['structure_diversity']:.4f}")
            click.echo(f"Splicetag diversity:     {metrics['splicetag_diversity']:.4f}")
            click.echo(f"Transcript diversity:    {metrics['transcript_diversity']:.4f}")
            click.echo(f"Variant diversity:       {metrics['variant_diversity']:.4f}")
            click.echo(f"Multi-exon rate:         {metrics['multi_exon_rate']:.4f} ({100*metrics['multi_exon_rate']:.1f}%)")
            click.echo(f"Variant rate:            {metrics['variant_rate']:.4f} ({100*metrics['variant_rate']:.1f}%)")

        # Clustering metrics
        if self.stats['clustering_metrics']:
            cm = self.stats['clustering_metrics']
            click.echo(f"\n🔗 CLUSTERING METRICS")
            click.echo("-" * 70)
            click.echo(f"Clustering basis:        {cm['clustering_basis']}")
            click.echo(f"Clustered reads:         {cm['clustered_reads']:,}")
            click.echo(f"Unique clusters:         {cm['unique_clusters']:,}")
            click.echo(f"Avg cluster size:        {cm['avg_cluster_size']:.2f}")
            click.echo(f"Compression ratio:       {cm['compression_ratio']:.1f}%")
            click.echo(f"Before clustering:       {cm['unique_before_clustering']:,} unique")

        # Top structures
        if self.tag_counters['XI']:
            click.echo(f"\n🏆 TOP 5 STRUCTURES (XI)")
            click.echo("-" * 70)
            for i, (structure_id, count) in enumerate(self.tag_counters['XI'].most_common(5), 1):
                pct = 100 * count / self.stats['total_reads']
                click.echo(f"{i:2}. {structure_id[:40]:40s} {count:6,} reads ({pct:5.2f}%)")

        # Top splicetags
        if self.tag_counters['XS']:
            click.echo(f"\n🔗 TOP 5 SPLICETAGS (XS)")
            click.echo("-" * 70)
            for i, (splicetag, count) in enumerate(self.tag_counters['XS'].most_common(5), 1):
                pct = 100 * count / self.stats['tag_counts']['XS']
                click.echo(f"{i:2}. {splicetag[:40]:40s} {count:6,} reads ({pct:5.2f}%)")

        # Top transcript groups
        if self.tag_counters['XT']:
            click.echo(f"\n🎯 TOP 5 TRANSCRIPT GROUPS (XT)")
            click.echo("-" * 70)
            for i, (transcript, count) in enumerate(self.tag_counters['XT'].most_common(5), 1):
                pct = 100 * count / self.stats['tag_counts']['XT']
                click.echo(f"{i:2}. {transcript[:40]:40s} {count:6,} reads ({pct:5.2f}%)")

        # Top clusters
        if self.tag_counters['XC']:
            click.echo(f"\n⭐ TOP 5 CLUSTERS (XC)")
            click.echo("-" * 70)
            for i, (cluster, count) in enumerate(self.tag_counters['XC'].most_common(5), 1):
                pct = 100 * count / self.stats['tag_counts']['XC']
                click.echo(f"{i:2}. {cluster[:40]:40s} {count:6,} reads ({pct:5.2f}%)")

    def export_json(self, output_file):
        """Export statistics to JSON"""
        # Convert Counter objects to regular dicts for JSON serialization
        export_data = dict(self.stats)

        # Add full tag counts (limited to top 1000 per tag to avoid huge files)
        export_data['tag_distributions'] = {
            tag: dict(counter.most_common(1000))
            for tag, counter in self.tag_counters.items()
            if len(counter) > 0
        }

        # Convert defaultdict to regular dict
        export_data['tag_co_occurrence'] = {
            k: dict(v) for k, v in self.stats['tag_co_occurrence'].items()
        }

        with open(output_file, 'w') as f:
            json.dump(export_data, f, indent=2)

        click.echo(f"💾 Statistics exported to JSON: {output_file}")

    def export_tsv(self, output_file):
        """Export statistics to TSV"""
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')

            # Write summary statistics
            writer.writerow(['Metric', 'Value'])
            writer.writerow(['Total_Reads', self.stats['total_reads']])
            writer.writerow(['Version', self.stats['version']])

            # Tag counts
            for tag in get_all_tag_names():
                writer.writerow([f'{tag}_Reads', self.stats['tag_counts'].get(tag, 0)])
                writer.writerow([f'{tag}_Unique', self.stats['unique_counts'].get(tag, 0)])

            # Complexity metrics
            if self.stats['complexity_metrics']:
                for metric, value in self.stats['complexity_metrics'].items():
                    writer.writerow([metric, value])

            # Clustering metrics
            if self.stats['clustering_metrics']:
                for metric, value in self.stats['clustering_metrics'].items():
                    writer.writerow([f'clustering_{metric}', value])

            writer.writerow([])  # Empty row

            # Tag distributions (top 100 per tag)
            for tag in get_all_tag_names():
                if self.tag_counters[tag]:
                    writer.writerow([f'{tag}_Tag', 'Count'])
                    for tag_value, count in self.tag_counters[tag].most_common(100):
                        writer.writerow([tag_value, count])
                    writer.writerow([])

        click.echo(f"📋 Statistics exported to TSV: {output_file}")


@click.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.option('--output', '-o', help='Output file (JSON or TSV based on extension)')
@click.option('--format', '-f', type=click.Choice(['json', 'tsv', 'auto']), default='auto', help='Output format')
@click.option('--limit', '-l', type=int, help='Limit number of reads to analyze (for testing)')
@click.option('--display/--no-display', default=True, help='Display summary to console')
def stats(bam_file, output, format, limit, display):
    """
    IsoTag Stats - Comprehensive Multi-Tag Analysis

    Analyze all IsoTag tags (XI, XB, XS, XT, XV, XC, XR) in BAM files.
    Provides statistics, distributions, and clustering metrics.

    Examples:
        # Basic analysis
        isotag_stats.py tagged.bam

        # Export to JSON
        isotag_stats.py tagged.bam -o stats.json

        # Export to TSV
        isotag_stats.py tagged.bam -o stats.tsv

        # Analyze first 1000 reads only
        isotag_stats.py tagged.bam --limit 1000
    """

    bam_path = Path(bam_file)
    if not bam_path.exists():
        click.echo(f"❌ BAM file not found: {bam_file}")
        sys.exit(1)

    analyzer = MultiTagStatsAnalyzer()

    # Analyze BAM file
    analyzer.analyze_bam(bam_file, limit)
    analyzer.compute_statistics()

    # Display results
    if display:
        analyzer.display_summary()

    # Export results
    if output:
        output_path = Path(output)

        # Determine format
        if format == 'auto':
            if output_path.suffix.lower() == '.json':
                format = 'json'
            elif output_path.suffix.lower() in ['.tsv', '.txt']:
                format = 'tsv'
            else:
                format = 'json'  # Default

        if format == 'json':
            analyzer.export_json(output_path)
        else:
            analyzer.export_tsv(output_path)

    click.echo(f"\n✅ Analysis complete!")


if __name__ == '__main__':
    stats()
