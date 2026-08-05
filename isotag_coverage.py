#!/usr/bin/env python3
"""
ISO-Tools Coverage - Coverage Analysis at Isotag Level (v2.1+)

Analyze coverage patterns at the isotag level.
Calculate per-isoform coverage, identify lowly covered regions,
and generate coverage statistics for each isotag.
"""

import subprocess
import click
import sys
import json
import csv
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import re


class IsotagCoverageAnalyzer:
    """Analyze coverage at isotag level"""

    def __init__(self, use_tag: str = 'XI'):
        self.use_tag = use_tag
        self.isotag_coverage = defaultdict(lambda: {'count': 0, 'total_bases': 0, 'reads': []})
        self.stats = {
            'total_reads': 0,
            'tagged_reads': 0,
            'unique_isotags': 0
        }

    def parse_cigar(self, cigar_string: str) -> List[Tuple[str, int]]:
        """Parse CIGAR string"""
        if not cigar_string or cigar_string == "*":
            return []

        operations = []
        pattern = r'(\d+)([MIDNSHP=X])'
        matches = re.findall(pattern, cigar_string)

        for length_str, op_char in matches:
            operations.append((op_char, int(length_str)))

        return operations

    def calculate_alignment_length(self, cigar_ops: List[Tuple[str, int]]) -> int:
        """Calculate aligned length from CIGAR"""
        length = 0
        for op, op_len in cigar_ops:
            if op in ['M', '=', 'X', 'D', 'N']:  # Operations that consume reference
                length += op_len
        return length

    def analyze_coverage(self, bam_file: str):
        """Analyze coverage from BAM file"""
        click.echo(f"🔍 Analyzing {self.use_tag} coverage from {bam_file}...")

        tag_prefix = f'{self.use_tag}:Z:'

        view_cmd = ['samtools', 'view', bam_file]

        try:
            process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

            for line in process.stdout:
                self.stats['total_reads'] += 1

                if self.stats['total_reads'] % 100000 == 0:
                    click.echo(f"   ⏳ Processed {self.stats['total_reads']:,} reads...")

                fields = line.strip().split('\t')

                if len(fields) < 11:
                    continue

                # Skip unmapped
                flag = int(fields[1])
                if flag & 0x4:
                    continue

                pos = int(fields[3])
                cigar = fields[5]

                # Extract tag
                tag_id = None
                for field in fields[11:]:
                    if field.startswith(tag_prefix):
                        tag_id = field[5:]
                        self.stats['tagged_reads'] += 1
                        break

                if not tag_id:
                    continue

                # Calculate coverage
                cigar_ops = self.parse_cigar(cigar)
                if cigar_ops:
                    aligned_length = self.calculate_alignment_length(cigar_ops)

                    self.isotag_coverage[tag_id]['count'] += 1
                    self.isotag_coverage[tag_id]['total_bases'] += aligned_length
                    self.isotag_coverage[tag_id]['reads'].append({
                        'qname': fields[0],
                        'pos': pos,
                        'length': aligned_length
                    })

            process.wait()

            self.stats['unique_isotags'] = len(self.isotag_coverage)

            click.echo(f"   ✅ Analyzed {self.stats['tagged_reads']:,} tagged reads")
            click.echo(f"   🆔 Found {self.stats['unique_isotags']:,} unique {self.use_tag} tags")

        except subprocess.CalledProcessError as e:
            click.echo(f"❌ Error reading BAM file: {e}")
            sys.exit(1)

    def compute_statistics(self):
        """Compute coverage statistics for each isotag"""
        click.echo(f"📊 Computing coverage statistics...")

        for tag_id, coverage_data in self.isotag_coverage.items():
            count = coverage_data['count']
            total_bases = coverage_data['total_bases']

            # Calculate mean coverage
            mean_coverage = total_bases / count if count > 0 else 0

            # Calculate length statistics
            lengths = [r['length'] for r in coverage_data['reads']]
            min_length = min(lengths) if lengths else 0
            max_length = max(lengths) if lengths else 0
            mean_length = sum(lengths) / len(lengths) if lengths else 0

            # Store computed stats
            coverage_data['mean_coverage'] = mean_coverage
            coverage_data['min_length'] = min_length
            coverage_data['max_length'] = max_length
            coverage_data['mean_length'] = mean_length

        click.echo(f"   ✅ Statistics computed for {len(self.isotag_coverage)} isotags")

    def export_coverage_tsv(self, output_file: str, min_count: int = 1):
        """Export coverage statistics to TSV"""
        click.echo(f"📄 Exporting coverage to {output_file}...")

        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')

            # Write header
            header = [
                f'{self.use_tag}_Tag',
                'Read_Count',
                'Total_Bases',
                'Mean_Coverage',
                'Min_Length',
                'Max_Length',
                'Mean_Length'
            ]
            writer.writerow(header)

            # Sort by read count (descending)
            sorted_isotags = sorted(
                self.isotag_coverage.items(),
                key=lambda x: x[1]['count'],
                reverse=True
            )

            exported = 0
            for tag_id, coverage_data in sorted_isotags:
                if coverage_data['count'] < min_count:
                    continue

                row = [
                    tag_id,
                    coverage_data['count'],
                    coverage_data['total_bases'],
                    f"{coverage_data['mean_coverage']:.2f}",
                    coverage_data['min_length'],
                    coverage_data['max_length'],
                    f"{coverage_data['mean_length']:.2f}"
                ]
                writer.writerow(row)
                exported += 1

        click.echo(f"   ✅ Exported coverage for {exported} isotags")

    def export_coverage_json(self, output_file: str, top_n: int = 100):
        """Export coverage statistics to JSON"""
        click.echo(f"📄 Exporting coverage summary to {output_file}...")

        # Sort by read count
        sorted_isotags = sorted(
            self.isotag_coverage.items(),
            key=lambda x: x[1]['count'],
            reverse=True
        )

        summary = {
            'statistics': self.stats,
            'tag_type': self.use_tag,
            'top_isotags': []
        }

        for tag_id, coverage_data in sorted_isotags[:top_n]:
            summary['top_isotags'].append({
                'tag_id': tag_id,
                'read_count': coverage_data['count'],
                'total_bases': coverage_data['total_bases'],
                'mean_coverage': round(coverage_data['mean_coverage'], 2),
                'min_length': coverage_data['min_length'],
                'max_length': coverage_data['max_length'],
                'mean_length': round(coverage_data['mean_length'], 2)
            })

        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2)

        click.echo(f"   ✅ Exported summary with top {min(top_n, len(sorted_isotags))} isotags")

    def display_summary(self, top_n: int = 10):
        """Display coverage summary"""
        click.echo(f"\n{'='*60}")
        click.echo(f"📊 ISO-TOOLS COVERAGE SUMMARY")
        click.echo(f"{'='*60}")

        click.echo(f"Tag type:          {self.use_tag}")
        click.echo(f"Total reads:       {self.stats['total_reads']:,}")
        click.echo(f"Tagged reads:      {self.stats['tagged_reads']:,}")
        click.echo(f"Unique isotags:    {self.stats['unique_isotags']:,}")

        # Calculate global statistics
        total_count = sum(data['count'] for data in self.isotag_coverage.values())
        total_bases = sum(data['total_bases'] for data in self.isotag_coverage.values())

        if total_count > 0:
            click.echo(f"Total bases:       {total_bases:,}")
            click.echo(f"Mean bases/read:   {total_bases / total_count:.2f}")

        # Top isotags
        sorted_isotags = sorted(
            self.isotag_coverage.items(),
            key=lambda x: x[1]['count'],
            reverse=True
        )

        if sorted_isotags:
            click.echo(f"\n🏆 TOP {top_n} {self.use_tag} TAGS BY COVERAGE")
            click.echo("-" * 60)

            for i, (tag_id, coverage_data) in enumerate(sorted_isotags[:top_n], 1):
                percentage = 100 * coverage_data['count'] / self.stats['tagged_reads']
                click.echo(f"{i:2d}. {tag_id[:35]:<35} "
                          f"{coverage_data['count']:>6,} reads ({percentage:5.2f}%) "
                          f"avg: {coverage_data['mean_length']:.0f}bp")


@click.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.option('--output', '-o', required=True, help='Output file (TSV or JSON)')
@click.option('--use-tag', '-t', type=click.Choice(['XI', 'XB', 'XS', 'XT', 'XC']),
              default='XI', help='Tag to use for coverage analysis (default: XI)')
@click.option('--format', '-f', type=click.Choice(['tsv', 'json', 'both']),
              default='tsv', help='Output format (default: tsv)')
@click.option('--min-count', default=1, help='Minimum read count per isotag (for TSV export)')
@click.option('--top', default=100, help='Number of top isotags in JSON summary')
def coverage(bam_file, output, use_tag, format, min_count, top):
    """
    ISO-Tools Coverage - Coverage Analysis at Isotag Level (v2.1+)

    Analyze coverage patterns at the isotag level.
    Calculate per-isoform coverage, identify lowly covered regions,
    and generate coverage statistics for each isotag.

    Analysis Modes:
        XI: Coverage per structure (full exon coordinates)
        XB: Coverage per boundary (5'/3' ends)
        XS: Coverage per splicetag (splice junctions)
        XT: Coverage per transcript (biological groups)
        XC: Coverage per cluster (after clustering)

    Examples:
        # Basic coverage analysis (XI, TSV output)
        isotag_coverage.py tagged.bam -o coverage.tsv

        # Coverage by splicetag (XS)
        isotag_coverage.py tagged.bam -o junction_coverage.tsv -t XS

        # Coverage by cluster (XC) with JSON output
        isotag_coverage.py clustered.bam -o cluster_coverage.json -t XC -f json

        # Filter low-count isotags
        isotag_coverage.py tagged.bam -o coverage.tsv --min-count 10

        # Both TSV and JSON output
        isotag_coverage.py tagged.bam -o coverage -f both

    Use Cases:
        - Identify most abundant isoforms
        - Filter lowly covered isoforms
        - Quality control for sequencing depth
        - Compare coverage across samples
    """

    analyzer = IsotagCoverageAnalyzer(use_tag=use_tag)

    # Analyze coverage
    analyzer.analyze_coverage(bam_file)

    # Compute statistics
    analyzer.compute_statistics()

    # Export results
    output_path = Path(output)

    if format == 'tsv' or format == 'both':
        if format == 'both':
            tsv_output = output_path.parent / f"{output_path.stem}.tsv"
        else:
            tsv_output = output

        analyzer.export_coverage_tsv(str(tsv_output), min_count)

    if format == 'json' or format == 'both':
        if format == 'both':
            json_output = output_path.parent / f"{output_path.stem}.json"
        else:
            json_output = output

        analyzer.export_coverage_json(str(json_output), top)

    # Display summary
    analyzer.display_summary()

    click.echo(f"\n✅ Coverage analysis complete!")


if __name__ == '__main__':
    coverage()
