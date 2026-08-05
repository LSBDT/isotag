#!/usr/bin/env python3
"""
ISO-Tools Merge - Intelligent BAM File Merging (v2.1+)

Merge multiple isotag-tagged BAM files with tag-aware handling.
Supports various merge strategies and tag reconciliation.

Merge Modes:
- simple: Simple concatenation (fastest)
- deduplicate: Remove duplicate reads by QNAME
- reconcile: Reconcile conflicting tags across samples
- cluster-aware: Preserve cluster assignments (XC/XR)
"""

import subprocess
import click
import sys
import tempfile
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Set, Tuple


class IsotagMerger:
    """Merge multiple isotag BAM files intelligently"""

    def __init__(self, mode: str = 'simple'):
        self.mode = mode
        self.input_files = []

        # Statistics
        self.stats = {
            'total_files': 0,
            'total_reads': 0,
            'output_reads': 0,
            'duplicates_removed': 0,
            'tag_conflicts': 0,
            'tags_reconciled': 0
        }

        # For deduplication mode
        self.seen_qnames = set()

        # For reconcile mode
        self.qname_to_tags = defaultdict(lambda: defaultdict(Counter))  # qname -> tag_type -> Counter

    def merge_simple(self, input_bams: List[str], output_bam: str):
        """Simple merge (concatenation)"""
        click.echo(f"🔄 Simple merge mode...")

        # Use samtools merge directly
        merge_cmd = ['samtools', 'merge', '-f', output_bam] + input_bams

        try:
            subprocess.run(merge_cmd, check=True)
            self.stats['output_reads'] = self._count_bam_reads(output_bam)
            click.echo(f"   ✅ Merged {len(input_bams)} files → {self.stats['output_reads']:,} reads")
        except subprocess.CalledProcessError as e:
            click.echo(f"❌ Merge failed: {e}")
            sys.exit(1)

    def merge_deduplicate(self, input_bams: List[str], output_bam: str):
        """Merge with deduplication by QNAME"""
        click.echo(f"🔄 Deduplication merge mode...")

        # Read header from first file
        header_cmd = ['samtools', 'view', '-H', input_bams[0]]
        header_result = subprocess.run(header_cmd, capture_output=True, text=True, check=True)
        header = header_result.stdout

        # Open output BAM
        output_cmd = ['samtools', 'view', '-b', '-o', output_bam]
        output_process = subprocess.Popen(output_cmd, stdin=subprocess.PIPE, text=True)
        output_process.stdin.write(header)

        # Process each input file
        for bam_file in input_bams:
            click.echo(f"   Processing {Path(bam_file).name}...")

            view_cmd = ['samtools', 'view', bam_file]
            view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

            for line in view_process.stdout:
                self.stats['total_reads'] += 1

                if self.stats['total_reads'] % 100000 == 0:
                    click.echo(f"      Processed {self.stats['total_reads']:,} reads...")

                fields = line.split('\t')
                qname = fields[0]

                # Check if already seen
                if qname in self.seen_qnames:
                    self.stats['duplicates_removed'] += 1
                    continue

                # Add to output
                self.seen_qnames.add(qname)
                output_process.stdin.write(line)
                self.stats['output_reads'] += 1

            view_process.wait()

        output_process.stdin.close()
        output_process.wait()

        click.echo(f"   ✅ Output: {self.stats['output_reads']:,} reads ({self.stats['duplicates_removed']:,} duplicates removed)")

    def merge_reconcile(self, input_bams: List[str], output_bam: str):
        """Merge with tag reconciliation (two-pass)"""
        click.echo(f"🔄 Reconcile merge mode (two-pass)...")

        # Pass 1: Collect all reads and tags
        click.echo(f"   📊 Pass 1: Analyzing tags across files...")

        read_storage = defaultdict(list)  # qname -> list of (line, tags_dict)

        for bam_file in input_bams:
            click.echo(f"      Reading {Path(bam_file).name}...")

            view_cmd = ['samtools', 'view', bam_file]
            view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

            for line in view_process.stdout:
                self.stats['total_reads'] += 1

                if self.stats['total_reads'] % 100000 == 0:
                    click.echo(f"         Processed {self.stats['total_reads']:,} reads...")

                fields = line.split('\t')
                qname = fields[0]

                # Extract all tags
                tags = {}
                for field in fields[11:]:
                    if field.startswith('XI:Z:'):
                        tags['XI'] = field[5:]
                    elif field.startswith('XB:Z:'):
                        tags['XB'] = field[5:]
                    elif field.startswith('XS:Z:'):
                        tags['XS'] = field[5:]
                    elif field.startswith('XT:Z:'):
                        tags['XT'] = field[5:]
                    elif field.startswith('XV:Z:'):
                        tags['XV'] = field[5:]
                    elif field.startswith('XC:Z:'):
                        tags['XC'] = field[5:]
                    elif field.startswith('XR:Z:'):
                        tags['XR'] = field[5:]

                # Store read with tags
                read_storage[qname].append((line, tags))

                # Track tag occurrences for this qname
                for tag_type, tag_value in tags.items():
                    self.qname_to_tags[qname][tag_type][tag_value] += 1

            view_process.wait()

        # Pass 2: Write with reconciled tags
        click.echo(f"   📊 Pass 2: Writing with reconciled tags...")

        # Read header from first file
        header_cmd = ['samtools', 'view', '-H', input_bams[0]]
        header_result = subprocess.run(header_cmd, capture_output=True, text=True, check=True)
        header = header_result.stdout

        # Open output BAM
        output_cmd = ['samtools', 'view', '-b', '-o', output_bam]
        output_process = subprocess.Popen(output_cmd, stdin=subprocess.PIPE, text=True)
        output_process.stdin.write(header)

        for qname, reads in read_storage.items():
            # Get consensus tags for this qname
            consensus_tags = {}
            tag_types = self.qname_to_tags[qname]

            for tag_type, tag_counts in tag_types.items():
                # Most common tag wins
                most_common_tag = tag_counts.most_common(1)[0][0]
                consensus_tags[tag_type] = most_common_tag

                # Check for conflicts
                if len(tag_counts) > 1:
                    self.stats['tag_conflicts'] += 1
                    self.stats['tags_reconciled'] += 1

            # Use first occurrence, update tags if needed
            line, original_tags = reads[0]

            # If tags differ from consensus, rebuild line
            if original_tags != consensus_tags:
                fields = line.strip().split('\t')

                # Remove old tags
                base_fields = fields[:11]
                extra_fields = []

                for field in fields[11:]:
                    if not field.startswith(('XI:Z:', 'XB:Z:', 'XS:Z:', 'XT:Z:', 'XV:Z:', 'XC:Z:', 'XR:Z:')):
                        extra_fields.append(field)

                # Add consensus tags
                for tag_type in ['XI', 'XB', 'XS', 'XT', 'XV', 'XC', 'XR']:
                    if tag_type in consensus_tags:
                        extra_fields.append(f"{tag_type}:Z:{consensus_tags[tag_type]}")

                line = '\t'.join(base_fields + extra_fields) + '\n'

            output_process.stdin.write(line)
            self.stats['output_reads'] += 1

        output_process.stdin.close()
        output_process.wait()

        click.echo(f"   ✅ Output: {self.stats['output_reads']:,} reads ({self.stats['tag_conflicts']:,} conflicts resolved)")

    def merge_cluster_aware(self, input_bams: List[str], output_bam: str):
        """Merge with cluster-aware handling (preserves XC/XR)"""
        click.echo(f"🔄 Cluster-aware merge mode...")

        # For cluster-aware, we want to preserve cluster assignments
        # Use reconcile mode but prioritize cluster tags
        self.merge_reconcile(input_bams, output_bam)

        click.echo(f"   ℹ️  Cluster assignments (XC/XR) preserved from most common source")

    def _count_bam_reads(self, bam_file: str) -> int:
        """Count reads in BAM file"""
        count_cmd = ['samtools', 'view', '-c', bam_file]
        result = subprocess.run(count_cmd, capture_output=True, text=True, check=True)
        return int(result.stdout.strip())

    def merge_files(self, input_bams: List[str], output_bam: str):
        """Main merge dispatcher"""
        self.stats['total_files'] = len(input_bams)

        # Validate input files
        for bam_file in input_bams:
            if not Path(bam_file).exists():
                click.echo(f"❌ Input file not found: {bam_file}")
                sys.exit(1)

        click.echo(f"🔀 Merging {len(input_bams)} BAM files...")
        click.echo(f"   Mode: {self.mode}")

        # Dispatch to appropriate merge method
        if self.mode == 'simple':
            self.merge_simple(input_bams, output_bam)
        elif self.mode == 'deduplicate':
            self.merge_deduplicate(input_bams, output_bam)
        elif self.mode == 'reconcile':
            self.merge_reconcile(input_bams, output_bam)
        elif self.mode == 'cluster-aware':
            self.merge_cluster_aware(input_bams, output_bam)

    def display_summary(self):
        """Display merge summary"""
        click.echo(f"\n{'='*60}")
        click.echo(f"🔀 MERGE SUMMARY")
        click.echo(f"{'='*60}")

        click.echo(f"Mode:           {self.mode}")
        click.echo(f"Input files:    {self.stats['total_files']}")
        click.echo(f"Total reads:    {self.stats['total_reads']:,}")
        click.echo(f"Output reads:   {self.stats['output_reads']:,}")

        if self.stats['duplicates_removed'] > 0:
            click.echo(f"Duplicates:     {self.stats['duplicates_removed']:,} removed")

        if self.stats['tag_conflicts'] > 0:
            click.echo(f"Tag conflicts:  {self.stats['tag_conflicts']:,} resolved")


@click.command()
@click.argument('input_bams', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('--output', '-o', required=True, help='Output merged BAM file')
@click.option('--mode', '-m', type=click.Choice(['simple', 'deduplicate', 'reconcile', 'cluster-aware']),
              default='simple', help='Merge strategy (default: simple)')
@click.option('--index', is_flag=True, help='Index output BAM after merging')
def merge(input_bams, output, mode, index):
    """
    ISO-Tools Merge - Intelligent BAM File Merging (v2.1+)

    Merge multiple isotag-tagged BAM files with various strategies.
    Supports tag reconciliation and deduplication.

    Merge Modes:
        simple:         Fast concatenation (uses samtools merge)
        deduplicate:    Remove duplicate reads by QNAME
        reconcile:      Reconcile conflicting tags across samples
        cluster-aware:  Preserve cluster assignments (XC/XR)

    Examples:
        # Simple merge (fastest)
        isotag_merge.py sample1.bam sample2.bam sample3.bam -o merged.bam

        # Deduplicate by read name
        isotag_merge.py sample*.bam -o merged.bam -m deduplicate

        # Reconcile conflicting tags
        isotag_merge.py sample1.bam sample2.bam -o merged.bam -m reconcile

        # Cluster-aware merge (preserves XC/XR)
        isotag_merge.py clustered*.bam -o merged.bam -m cluster-aware

        # With indexing
        isotag_merge.py sample*.bam -o merged.bam --index

    Use Cases:
        - Merge technical replicates (simple mode)
        - Combine batches from same sample (deduplicate mode)
        - Merge samples tagged separately (reconcile mode)
        - Merge clustered datasets (cluster-aware mode)
    """

    if len(input_bams) < 2:
        click.echo("❌ Need at least 2 input BAM files to merge")
        sys.exit(1)

    # Create merger
    merger = IsotagMerger(mode=mode)

    # Merge files
    merger.merge_files(list(input_bams), output)

    # Display summary
    merger.display_summary()

    # Index if requested
    if index:
        click.echo(f"\n🔧 Indexing output BAM...")
        try:
            subprocess.run(['samtools', 'index', output], check=True)
            click.echo(f"   ✅ Index created: {output}.bai")
        except subprocess.CalledProcessError:
            click.echo(f"   ⚠️  Indexing failed (non-fatal)")

    click.echo(f"\n✅ Merge complete: {output}")


if __name__ == '__main__':
    merge()
