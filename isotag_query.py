#!/usr/bin/env python3
"""
ISO-Tools Query - Query Interface for Isotag Data (v2.1+)

Interactive and command-line query interface for isotag-tagged BAM files.
Search, filter, and extract reads by tag values, genomic regions, or patterns.
"""

import subprocess
import click
import sys
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Optional


class IsotagQuery:
    """Query interface for isotag-tagged BAM files"""

    def __init__(self, bam_file: str):
        self.bam_file = bam_file
        self.results = []
        self.stats = {
            'reads_scanned': 0,
            'reads_matched': 0
        }

    def query_by_tag(self, tag_type: str, tag_values: List[str], output_bam: Optional[str] = None):
        """Query reads by specific tag values"""
        click.echo(f"🔍 Querying {tag_type} tags: {', '.join(tag_values[:5])}{'...' if len(tag_values) > 5 else ''}")

        tag_set = set(tag_values)
        tag_prefix = f'{tag_type}:Z:'

        # If output_bam specified, write to BAM, else collect results
        if output_bam:
            return self._query_to_bam(tag_prefix, tag_set, output_bam)
        else:
            return self._query_to_list(tag_prefix, tag_set)

    def _query_to_bam(self, tag_prefix: str, tag_set: Set[str], output_bam: str):
        """Query and write results to BAM file"""
        # Read header
        header_cmd = ['samtools', 'view', '-H', self.bam_file]
        header_result = subprocess.run(header_cmd, capture_output=True, text=True, check=True)
        header = header_result.stdout

        # Open output BAM
        output_cmd = ['samtools', 'view', '-b', '-o', output_bam]
        output_process = subprocess.Popen(output_cmd, stdin=subprocess.PIPE, text=True)
        output_process.stdin.write(header)

        # Query reads
        view_cmd = ['samtools', 'view', self.bam_file]
        view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

        for line in view_process.stdout:
            self.stats['reads_scanned'] += 1

            if self.stats['reads_scanned'] % 100000 == 0:
                click.echo(f"   ⏳ Scanned {self.stats['reads_scanned']:,} reads, matched {self.stats['reads_matched']:,}...")

            fields = line.split('\t')
            matched = False

            for field in fields[11:]:
                if field.startswith(tag_prefix):
                    tag_value = field[5:]
                    if tag_value in tag_set:
                        matched = True
                        break

            if matched:
                output_process.stdin.write(line)
                self.stats['reads_matched'] += 1

        view_process.wait()
        output_process.stdin.close()
        output_process.wait()

        click.echo(f"   ✅ Found {self.stats['reads_matched']:,} reads matching query")
        click.echo(f"   📄 Output: {output_bam}")

        return self.stats['reads_matched']

    def _query_to_list(self, tag_prefix: str, tag_set: Set[str]):
        """Query and collect results in memory"""
        view_cmd = ['samtools', 'view', self.bam_file]
        view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

        for line in view_process.stdout:
            self.stats['reads_scanned'] += 1

            if self.stats['reads_scanned'] % 100000 == 0:
                click.echo(f"   ⏳ Scanned {self.stats['reads_scanned']:,} reads...")

            fields = line.split('\t')

            for field in fields[11:]:
                if field.startswith(tag_prefix):
                    tag_value = field[5:]
                    if tag_value in tag_set:
                        self.results.append(line)
                        self.stats['reads_matched'] += 1
                        break

        view_process.wait()

        click.echo(f"   ✅ Found {self.stats['reads_matched']:,} reads matching query")

        return self.results

    def query_by_region(self, region: str, tag_filter: Optional[Dict[str, List[str]]] = None,
                       output_bam: Optional[str] = None):
        """Query reads by genomic region with optional tag filtering"""
        click.echo(f"🔍 Querying region: {region}")

        # Use samtools view with region
        view_cmd = ['samtools', 'view', self.bam_file, region]

        if output_bam:
            # Write to output BAM
            header_cmd = ['samtools', 'view', '-H', self.bam_file]
            header_result = subprocess.run(header_cmd, capture_output=True, text=True, check=True)
            header = header_result.stdout

            output_cmd = ['samtools', 'view', '-b', '-o', output_bam]
            output_process = subprocess.Popen(output_cmd, stdin=subprocess.PIPE, text=True)
            output_process.stdin.write(header)

            view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

            for line in view_process.stdout:
                self.stats['reads_scanned'] += 1

                # Apply tag filter if specified
                if tag_filter:
                    matched = self._check_tag_filter(line, tag_filter)
                    if not matched:
                        continue

                output_process.stdin.write(line)
                self.stats['reads_matched'] += 1

            view_process.wait()
            output_process.stdin.close()
            output_process.wait()

            click.echo(f"   ✅ Found {self.stats['reads_matched']:,} reads in region")
            click.echo(f"   📄 Output: {output_bam}")

        else:
            # Collect in memory
            view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

            for line in view_process.stdout:
                self.stats['reads_scanned'] += 1

                if tag_filter:
                    matched = self._check_tag_filter(line, tag_filter)
                    if not matched:
                        continue

                self.results.append(line)
                self.stats['reads_matched'] += 1

            view_process.wait()

            click.echo(f"   ✅ Found {self.stats['reads_matched']:,} reads in region")

        return self.stats['reads_matched']

    def _check_tag_filter(self, line: str, tag_filter: Dict[str, List[str]]) -> bool:
        """Check if read matches tag filter"""
        fields = line.split('\t')

        for tag_type, required_values in tag_filter.items():
            tag_prefix = f'{tag_type}:Z:'
            matched = False

            for field in fields[11:]:
                if field.startswith(tag_prefix):
                    tag_value = field[5:]
                    if tag_value in required_values:
                        matched = True
                        break

            if not matched:
                return False

        return True

    def count_by_tag(self, tag_type: str, top_n: int = 20):
        """Count and rank tags"""
        click.echo(f"🔢 Counting {tag_type} tags...")

        tag_prefix = f'{tag_type}:Z:'
        tag_counter = Counter()

        view_cmd = ['samtools', 'view', self.bam_file]
        view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

        for line in view_process.stdout:
            self.stats['reads_scanned'] += 1

            if self.stats['reads_scanned'] % 100000 == 0:
                click.echo(f"   ⏳ Processed {self.stats['reads_scanned']:,} reads...")

            fields = line.split('\t')

            for field in fields[11:]:
                if field.startswith(tag_prefix):
                    tag_value = field[5:]
                    tag_counter[tag_value] += 1
                    break

        view_process.wait()

        click.echo(f"   ✅ Found {len(tag_counter)} unique {tag_type} tags")
        click.echo(f"\n🏆 TOP {top_n} {tag_type} TAGS")
        click.echo("-" * 60)

        for i, (tag, count) in enumerate(tag_counter.most_common(top_n), 1):
            percentage = 100 * count / self.stats['reads_scanned']
            click.echo(f"{i:2d}. {tag[:40]:<40} {count:>8,} ({percentage:5.2f}%)")

        return tag_counter

    def find_reads_with_pattern(self, tag_type: str, pattern: str, output_bam: Optional[str] = None):
        """Find reads with tags matching a pattern (substring match)"""
        click.echo(f"🔍 Finding {tag_type} tags containing: '{pattern}'")

        tag_prefix = f'{tag_type}:Z:'

        if output_bam:
            return self._find_pattern_to_bam(tag_prefix, pattern, output_bam)
        else:
            return self._find_pattern_to_list(tag_prefix, pattern)

    def _find_pattern_to_bam(self, tag_prefix: str, pattern: str, output_bam: str):
        """Find pattern and write to BAM"""
        header_cmd = ['samtools', 'view', '-H', self.bam_file]
        header_result = subprocess.run(header_cmd, capture_output=True, text=True, check=True)
        header = header_result.stdout

        output_cmd = ['samtools', 'view', '-b', '-o', output_bam]
        output_process = subprocess.Popen(output_cmd, stdin=subprocess.PIPE, text=True)
        output_process.stdin.write(header)

        view_cmd = ['samtools', 'view', self.bam_file]
        view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

        for line in view_process.stdout:
            self.stats['reads_scanned'] += 1

            if self.stats['reads_scanned'] % 100000 == 0:
                click.echo(f"   ⏳ Scanned {self.stats['reads_scanned']:,} reads...")

            fields = line.split('\t')

            for field in fields[11:]:
                if field.startswith(tag_prefix):
                    tag_value = field[5:]
                    if pattern in tag_value:
                        output_process.stdin.write(line)
                        self.stats['reads_matched'] += 1
                        break

        view_process.wait()
        output_process.stdin.close()
        output_process.wait()

        click.echo(f"   ✅ Found {self.stats['reads_matched']:,} reads with pattern")
        click.echo(f"   📄 Output: {output_bam}")

        return self.stats['reads_matched']

    def _find_pattern_to_list(self, tag_prefix: str, pattern: str):
        """Find pattern and collect in memory"""
        view_cmd = ['samtools', 'view', self.bam_file]
        view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

        for line in view_process.stdout:
            self.stats['reads_scanned'] += 1

            fields = line.split('\t')

            for field in fields[11:]:
                if field.startswith(tag_prefix):
                    tag_value = field[5:]
                    if pattern in tag_value:
                        self.results.append(line)
                        self.stats['reads_matched'] += 1
                        break

        view_process.wait()

        click.echo(f"   ✅ Found {self.stats['reads_matched']:,} reads with pattern")

        return self.results


@click.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.option('--query-type', '-q', type=click.Choice(['tag', 'region', 'count', 'pattern']),
              required=True, help='Type of query to perform')
@click.option('--tag', '-t', type=click.Choice(['XI', 'XB', 'XS', 'XT', 'XV', 'XC', 'XR']),
              help='Tag type to query')
@click.option('--values', '-v', help='Comma-separated tag values to search for (for tag/pattern query)')
@click.option('--pattern', '-p', help='Pattern to search within tags (substring match)')
@click.option('--region', '-r', help='Genomic region (chr:start-end) for region query')
@click.option('--output', '-o', help='Output BAM file for results')
@click.option('--top', default=20, help='Number of top results to show (for count query)')
def query(bam_file, query_type, tag, values, pattern, region, output, top):
    """
    ISO-Tools Query - Query Interface for Isotag Data (v2.1+)

    Interactive and command-line query interface for isotag-tagged BAM files.
    Search, filter, and extract reads by tag values, genomic regions, or patterns.

    Query Types:
        tag:     Query by specific tag values
        region:  Query by genomic region
        count:   Count and rank tags
        pattern: Find tags containing pattern (substring match)

    Examples:
        # Query reads with specific XI tags
        isotag_query.py tagged.bam -q tag -t XI -v "tag1,tag2,tag3" -o subset.bam

        # Query reads in genomic region
        isotag_query.py tagged.bam -q region -r "chr1:1000-5000" -o region.bam

        # Count and rank XC (cluster) tags
        isotag_query.py tagged.bam -q count -t XC --top 50

        # Find XI tags containing pattern
        isotag_query.py tagged.bam -q pattern -t XI -p "aKF498" -o pattern_matches.bam

        # Count XS (splicetag) usage
        isotag_query.py tagged.bam -q count -t XS --top 100

        # Query reads with specific clusters
        isotag_query.py clustered.bam -q tag -t XC -v "cluster_id1,cluster_id2" -o clusters.bam

    Use Cases:
        - Extract reads for specific isoforms
        - Find all reads in a genomic region
        - Rank most abundant tags
        - Search for tag patterns
        - Filter by cluster assignment
    """

    querier = IsotagQuery(bam_file)

    if query_type == 'tag':
        if not tag or not values:
            click.echo("❌ --tag and --values required for tag query")
            sys.exit(1)

        value_list = [v.strip() for v in values.split(',')]
        querier.query_by_tag(tag, value_list, output)

    elif query_type == 'region':
        if not region:
            click.echo("❌ --region required for region query")
            sys.exit(1)

        querier.query_by_region(region, output_bam=output)

    elif query_type == 'count':
        if not tag:
            click.echo("❌ --tag required for count query")
            sys.exit(1)

        querier.count_by_tag(tag, top)

    elif query_type == 'pattern':
        if not tag or not pattern:
            click.echo("❌ --tag and --pattern required for pattern query")
            sys.exit(1)

        querier.find_reads_with_pattern(tag, pattern, output)

    # Display stats
    click.echo(f"\n📊 Query Statistics")
    click.echo(f"   Reads scanned: {querier.stats['reads_scanned']:,}")
    click.echo(f"   Reads matched: {querier.stats['reads_matched']:,}")

    if querier.stats['reads_scanned'] > 0:
        match_rate = 100 * querier.stats['reads_matched'] / querier.stats['reads_scanned']
        click.echo(f"   Match rate:    {match_rate:.2f}%")

    click.echo(f"\n✅ Query complete!")


if __name__ == '__main__':
    query()
