#!/usr/bin/env python3
"""
ISO-Tools Subset - Smart BAM File Subsetting (v2.1+)

Create representative subsets of large isotag-tagged BAM files.
Supports multiple subsetting strategies optimized for different analyses.

Subsetting Modes:
- random: Random sampling (fastest)
- diverse: Maximize isotag diversity (structure-aware)
- proportional: Maintain isotag frequency distribution
- cluster: Sample from each cluster (XC-aware)
- boundary: Sample diverse boundaries (XB-aware)
- splicetag: Sample diverse splice patterns (XS-aware)
"""

import subprocess
import click
import sys
import random
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple


class IsotagSubsetter:
    """Create representative subsets of isotag BAM files"""

    def __init__(self, mode: str = 'random', target_reads: int = 10000):
        self.mode = mode
        self.target_reads = target_reads

        # Statistics
        self.total_reads = 0
        self.selected_reads = 0

        # Tag tracking for smart sampling
        self.xi_counts = Counter()
        self.xb_counts = Counter()
        self.xs_counts = Counter()
        self.xc_counts = Counter()

        # Selected tags (for diversity tracking)
        self.selected_xi = set()
        self.selected_xb = set()
        self.selected_xs = set()
        self.selected_xc = set()

        # Read storage for two-pass modes
        self.read_storage = []  # [(read_line, tags_dict)]

    def _extract_tags(self, line: str) -> Dict[str, str]:
        """Extract all tags from a read line"""
        fields = line.split('\t')
        tags = {}

        if len(fields) < 11:
            return tags

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

        return tags

    def _should_select_random(self) -> bool:
        """Random sampling decision"""
        if self.total_reads == 0:
            return False

        # Reservoir sampling
        if self.selected_reads < self.target_reads:
            return True
        else:
            # Probability: target_reads / total_reads
            return random.random() < (self.target_reads / self.total_reads)

    def _should_select_diverse(self, tags: Dict[str, str]) -> bool:
        """Diversity-based sampling (prefer unseen isotags)"""
        if self.selected_reads >= self.target_reads:
            return False

        xi_tag = tags.get('XI')
        if not xi_tag:
            return False

        # Prioritize new XI tags
        if xi_tag not in self.selected_xi:
            return True

        # Once we've seen enough unique isotags, sample proportionally
        if len(self.selected_xi) > self.target_reads * 0.8:
            return self._should_select_random()

        # Otherwise, sample less frequently seen isotags
        current_freq = self.xi_counts[xi_tag]
        return random.random() < (1.0 / (1 + current_freq))

    def _should_select_cluster(self, tags: Dict[str, str]) -> bool:
        """Cluster-aware sampling (sample from each cluster)"""
        if self.selected_reads >= self.target_reads:
            return False

        xc_tag = tags.get('XC')
        if not xc_tag:
            # If no cluster tag, fallback to diverse mode
            return self._should_select_diverse(tags)

        # Try to sample from each cluster
        if xc_tag not in self.selected_xc:
            return True

        # Sample proportionally within clusters
        cluster_count = self.xc_counts[xc_tag]
        return random.random() < (1.0 / (1 + cluster_count))

    def _should_select_boundary(self, tags: Dict[str, str]) -> bool:
        """Boundary-aware sampling (diverse XB tags)"""
        if self.selected_reads >= self.target_reads:
            return False

        xb_tag = tags.get('XB')
        if not xb_tag:
            return self._should_select_diverse(tags)

        # Prioritize new boundaries
        if xb_tag not in self.selected_xb:
            return True

        boundary_count = self.xb_counts[xb_tag]
        return random.random() < (1.0 / (1 + boundary_count))

    def _should_select_splicetag(self, tags: Dict[str, str]) -> bool:
        """Splicetag-aware sampling (diverse XS tags)"""
        if self.selected_reads >= self.target_reads:
            return False

        xs_tag = tags.get('XS')
        if not xs_tag:
            return self._should_select_diverse(tags)

        # Prioritize new splice patterns
        if xs_tag not in self.selected_xs:
            return True

        splicetag_count = self.xs_counts[xs_tag]
        return random.random() < (1.0 / (1 + splicetag_count))

    def _compute_selection_probability(self, tags: Dict[str, str]) -> float:
        """Compute selection probability for proportional mode"""
        xi_tag = tags.get('XI')
        if not xi_tag:
            return 0.0

        # Use frequency from first pass
        xi_freq = self.xi_counts[xi_tag]
        total_xi = sum(self.xi_counts.values())

        if total_xi == 0:
            return 0.0

        # Proportional probability
        target_for_this_xi = (xi_freq / total_xi) * self.target_reads
        selected_for_this_xi = sum(1 for t in self.selected_xi if t == xi_tag)

        if selected_for_this_xi >= target_for_this_xi:
            return 0.0

        return min(1.0, (target_for_this_xi - selected_for_this_xi) / xi_freq)

    def subset_bam(self, input_bam: str, output_bam: str):
        """Create subset BAM file"""
        click.echo(f"🔄 Creating subset using '{self.mode}' mode...")
        click.echo(f"   Target: {self.target_reads:,} reads")

        # For proportional mode, need two passes
        if self.mode == 'proportional':
            self._subset_bam_proportional(input_bam, output_bam)
        else:
            self._subset_bam_streaming(input_bam, output_bam)

        # Display statistics
        self._display_statistics()

    def _subset_bam_streaming(self, input_bam: str, output_bam: str):
        """Single-pass streaming subset (all modes except proportional)"""
        # Read header
        header_cmd = ['samtools', 'view', '-H', input_bam]
        header_result = subprocess.run(header_cmd, capture_output=True, text=True, check=True)
        header = header_result.stdout

        # Open output BAM
        output_cmd = ['samtools', 'view', '-b', '-o', output_bam]
        output_process = subprocess.Popen(output_cmd, stdin=subprocess.PIPE, text=True)

        # Write header
        output_process.stdin.write(header)

        # Process reads
        view_cmd = ['samtools', 'view', input_bam]
        view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

        for line in view_process.stdout:
            self.total_reads += 1

            if self.total_reads % 100000 == 0:
                click.echo(f"   ⏳ Processed {self.total_reads:,} reads, selected {self.selected_reads:,}...")

            # Extract tags
            tags = self._extract_tags(line)

            # Update counters (for statistics)
            if 'XI' in tags:
                self.xi_counts[tags['XI']] += 1
            if 'XB' in tags:
                self.xb_counts[tags['XB']] += 1
            if 'XS' in tags:
                self.xs_counts[tags['XS']] += 1
            if 'XC' in tags:
                self.xc_counts[tags['XC']] += 1

            # Selection logic based on mode
            should_select = False
            if self.mode == 'random':
                should_select = self._should_select_random()
            elif self.mode == 'diverse':
                should_select = self._should_select_diverse(tags)
            elif self.mode == 'cluster':
                should_select = self._should_select_cluster(tags)
            elif self.mode == 'boundary':
                should_select = self._should_select_boundary(tags)
            elif self.mode == 'splicetag':
                should_select = self._should_select_splicetag(tags)

            if should_select:
                output_process.stdin.write(line)
                self.selected_reads += 1

                # Track selected tags
                if 'XI' in tags:
                    self.selected_xi.add(tags['XI'])
                if 'XB' in tags:
                    self.selected_xb.add(tags['XB'])
                if 'XS' in tags:
                    self.selected_xs.add(tags['XS'])
                if 'XC' in tags:
                    self.selected_xc.add(tags['XC'])

        view_process.wait()
        output_process.stdin.close()
        output_process.wait()

        click.echo(f"   ✅ Selected {self.selected_reads:,} reads from {self.total_reads:,} total")

    def _subset_bam_proportional(self, input_bam: str, output_bam: str):
        """Two-pass proportional subset (maintains frequency distribution)"""
        click.echo(f"   📊 Pass 1: Analyzing isotag frequencies...")

        # First pass: count frequencies
        view_cmd = ['samtools', 'view', input_bam]
        view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

        for line in view_process.stdout:
            self.total_reads += 1

            if self.total_reads % 100000 == 0:
                click.echo(f"      Processed {self.total_reads:,} reads...")

            tags = self._extract_tags(line)

            if 'XI' in tags:
                self.xi_counts[tags['XI']] += 1

            # Store reads for second pass
            self.read_storage.append((line, tags))

        view_process.wait()

        click.echo(f"   ✅ Found {len(self.xi_counts)} unique isotags")
        click.echo(f"   📊 Pass 2: Proportional sampling...")

        # Second pass: proportional sampling
        header_cmd = ['samtools', 'view', '-H', input_bam]
        header_result = subprocess.run(header_cmd, capture_output=True, text=True, check=True)
        header = header_result.stdout

        output_cmd = ['samtools', 'view', '-b', '-o', output_bam]
        output_process = subprocess.Popen(output_cmd, stdin=subprocess.PIPE, text=True)
        output_process.stdin.write(header)

        for line, tags in self.read_storage:
            probability = self._compute_selection_probability(tags)

            if random.random() < probability:
                output_process.stdin.write(line)
                self.selected_reads += 1

                if 'XI' in tags:
                    self.selected_xi.add(tags['XI'])

        output_process.stdin.close()
        output_process.wait()

        click.echo(f"   ✅ Selected {self.selected_reads:,} reads proportionally")

    def _display_statistics(self):
        """Display subsetting statistics"""
        click.echo(f"\n{'='*60}")
        click.echo(f"📊 SUBSETTING STATISTICS")
        click.echo(f"{'='*60}")

        click.echo(f"Mode:                  {self.mode}")
        click.echo(f"Total reads:           {self.total_reads:,}")
        click.echo(f"Selected reads:        {self.selected_reads:,} ({100*self.selected_reads/self.total_reads:.1f}%)")
        click.echo(f"Target reads:          {self.target_reads:,}")

        click.echo(f"\n🆔 DIVERSITY METRICS")
        click.echo(f"-" * 40)
        click.echo(f"Total unique XI:       {len(self.xi_counts):,}")
        click.echo(f"Selected unique XI:    {len(self.selected_xi):,} ({100*len(self.selected_xi)/len(self.xi_counts):.1f}%)")

        if len(self.xb_counts) > 0:
            click.echo(f"Total unique XB:       {len(self.xb_counts):,}")
            click.echo(f"Selected unique XB:    {len(self.selected_xb):,} ({100*len(self.selected_xb)/len(self.xb_counts):.1f}%)")

        if len(self.xs_counts) > 0:
            click.echo(f"Total unique XS:       {len(self.xs_counts):,}")
            click.echo(f"Selected unique XS:    {len(self.selected_xs):,} ({100*len(self.selected_xs)/len(self.xs_counts):.1f}%)")

        if len(self.xc_counts) > 0:
            click.echo(f"Total unique XC:       {len(self.xc_counts):,}")
            click.echo(f"Selected unique XC:    {len(self.selected_xc):,} ({100*len(self.selected_xc)/len(self.xc_counts):.1f}%)")


@click.command()
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('output_bam', type=click.Path())
@click.option('--mode', '-m', type=click.Choice(['random', 'diverse', 'proportional', 'cluster', 'boundary', 'splicetag']),
              default='diverse', help='Subsetting strategy (default: diverse)')
@click.option('--reads', '-n', default=10000, help='Target number of reads (default: 10000)')
@click.option('--seed', '-s', type=int, help='Random seed for reproducibility')
def subset(input_bam, output_bam, mode, reads, seed):
    """
    ISO-Tools Subset - Smart BAM File Subsetting (v2.1+)

    Create representative subsets of large isotag-tagged BAM files.
    Supports multiple subsetting strategies optimized for different analyses.

    Subsetting Modes:
        random:       Random sampling (fastest, no bias)
        diverse:      Maximize isotag diversity (default, structure-aware)
        proportional: Maintain isotag frequency distribution (two-pass)
        cluster:      Sample from each cluster (XC-aware)
        boundary:     Sample diverse boundaries (XB-aware)
        splicetag:    Sample diverse splice patterns (XS-aware)

    Examples:
        # Diverse subset (10K reads, maximize unique isotags)
        isotag_subset.py large.bam subset.bam

        # Random subset (50K reads)
        isotag_subset.py large.bam subset.bam -m random -n 50000

        # Proportional subset (maintains frequency distribution)
        isotag_subset.py large.bam subset.bam -m proportional -n 20000

        # Cluster-aware subset (samples from each cluster)
        isotag_subset.py clustered.bam subset.bam -m cluster -n 10000

        # Boundary-aware subset (diverse 5'/3' boundaries)
        isotag_subset.py tagged.bam subset.bam -m boundary -n 10000

        # Splicetag-aware subset (diverse splice junctions)
        isotag_subset.py tagged.bam subset.bam -m splicetag -n 10000

        # Reproducible subset (fixed seed)
        isotag_subset.py large.bam subset.bam --seed 42
    """

    # Set random seed if provided
    if seed is not None:
        random.seed(seed)
        click.echo(f"🎲 Random seed: {seed}")

    # Create subsetter
    subsetter = IsotagSubsetter(mode=mode, target_reads=reads)

    # Create subset
    subsetter.subset_bam(input_bam, output_bam)

    # Index output BAM
    click.echo(f"\n🔧 Indexing output BAM...")
    try:
        subprocess.run(['samtools', 'index', output_bam], check=True)
        click.echo(f"   ✅ Index created: {output_bam}.bai")
    except subprocess.CalledProcessError:
        click.echo(f"   ⚠️  Indexing failed (non-fatal)")

    click.echo(f"\n✅ Subset complete: {output_bam}")


if __name__ == '__main__':
    subset()
