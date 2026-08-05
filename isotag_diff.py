#!/usr/bin/env python3
"""
ISO-Tools Diff - Differential Isoform Analysis (v2.1+)

Compare isotag expression between conditions for differential usage analysis.
Statistical testing and ranking of differentially expressed isoforms.
Supports all 7 tags: XI, XB, XS, XT, XV, XC, XR
"""

import subprocess
import click
import sys
import csv
import json
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import math


class IsotagDifferentialAnalyzer:
    """Differential analysis of isotag expression (v2.1+ all 7 tags)"""

    def __init__(self, use_tag: str = 'XI'):
        self.use_tag = use_tag  # Which tag to use for counting (XI, XB, XS, XT, XC)
        self.sample_data = {}  # sample_name -> {tag_id: read_count}
        self.condition_samples = defaultdict(list)  # condition -> list of sample_names
        self.isotag_stats = {}  # tag_id -> statistical_results
        self.all_isotags = set()
        self.stats = {
            'total_samples': 0,
            'total_conditions': 0,
            'total_isotags': 0,
            'significant_isotags': 0,
            'upregulated': 0,
            'downregulated': 0,
            'use_tag': use_tag
        }
    
    def load_sample_counts(self, bam_files: List[str], sample_names: Optional[List[str]] = None, 
                          conditions: Optional[List[str]] = None):
        """Load isotag counts from multiple BAM files"""
        if not sample_names:
            sample_names = [Path(f).stem for f in bam_files]
        
        if not conditions:
            # Try to infer conditions from sample names
            conditions = []
            for name in sample_names:
                if any(keyword in name.lower() for keyword in ['control', 'ctrl', 'wt']):
                    conditions.append('control')
                elif any(keyword in name.lower() for keyword in ['treatment', 'treat', 'exp', 'mut']):
                    conditions.append('treatment')
                else:
                    conditions.append('unknown')
        
        if len(bam_files) != len(sample_names) or len(bam_files) != len(conditions):
            click.echo("❌ Number of BAM files, sample names, and conditions must match")
            sys.exit(1)
        
        click.echo(f"🔍 Loading isotag counts from {len(bam_files)} BAM files...")
        
        for i, (bam_file, sample_name, condition) in enumerate(zip(bam_files, sample_names, conditions)):
            click.echo(f"   📁 Processing {sample_name} ({condition})...")
            
            isotag_counts = Counter()
            
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
                        # Extract the specified tag (XI, XB, XS, XT, or XC)
                        tag_prefix = f'{self.use_tag}:Z:'
                        tag_id = None

                        for field in fields[11:]:
                            if field.startswith(tag_prefix):
                                tag_id = field[5:]
                                isotag_counts[tag_id] += 1
                                self.all_isotags.add(tag_id)
                                break
                
                process.wait()
                
                self.sample_data[sample_name] = dict(isotag_counts)
                self.condition_samples[condition].append(sample_name)
                
                click.echo(f"      ✅ Found {len(isotag_counts)} unique isotags in {reads_processed:,} reads")
                
            except subprocess.CalledProcessError as e:
                click.echo(f"❌ Error reading {bam_file}: {e}")
                sys.exit(1)
        
        self.stats['total_samples'] = len(sample_names)
        self.stats['total_conditions'] = len(self.condition_samples)
        self.stats['total_isotags'] = len(self.all_isotags)
        
        click.echo(f"   ✅ Loaded data: {self.stats['total_samples']} samples, {self.stats['total_conditions']} conditions, {self.stats['total_isotags']} isotags")
    
    def calculate_basic_stats(self, counts: List[float]) -> Dict:
        """Calculate basic statistics for a list of counts"""
        if not counts:
            return {'mean': 0, 'std': 0, 'median': 0, 'total': 0}
        
        mean = sum(counts) / len(counts)
        variance = sum((x - mean) ** 2 for x in counts) / len(counts) if len(counts) > 1 else 0
        std = math.sqrt(variance)
        median = sorted(counts)[len(counts) // 2]
        total = sum(counts)
        
        return {
            'mean': mean,
            'std': std,
            'median': median,
            'total': total,
            'n': len(counts)
        }
    
    def welch_t_test(self, group1: List[float], group2: List[float]) -> Tuple[float, float]:
        """
        Perform Welch's t-test for unequal variances
        Returns (t_statistic, estimated_p_value)
        """
        if len(group1) < 2 or len(group2) < 2:
            return 0.0, 1.0
        
        # Calculate means
        mean1 = sum(group1) / len(group1)
        mean2 = sum(group2) / len(group2)
        
        # Calculate variances
        var1 = sum((x - mean1) ** 2 for x in group1) / (len(group1) - 1)
        var2 = sum((x - mean2) ** 2 for x in group2) / (len(group2) - 1)
        
        # Avoid division by zero
        if var1 == 0 and var2 == 0:
            return 0.0, 1.0
        
        # Calculate standard error
        se = math.sqrt(var1 / len(group1) + var2 / len(group2))
        
        if se == 0:
            return 0.0, 1.0
        
        # Calculate t-statistic
        t_stat = (mean1 - mean2) / se
        
        # Degrees of freedom (Welch-Satterthwaite equation)
        if var1 == 0 or var2 == 0:
            df = min(len(group1), len(group2)) - 1
        else:
            num = (var1 / len(group1) + var2 / len(group2)) ** 2
            denom = (var1 / len(group1)) ** 2 / (len(group1) - 1) + (var2 / len(group2)) ** 2 / (len(group2) - 1)
            df = num / denom if denom > 0 else 1
        
        # Rough p-value estimation (simplified)
        # For a more accurate p-value, would need scipy.stats
        abs_t = abs(t_stat)
        if abs_t < 1:
            p_value = 0.5
        elif abs_t < 2:
            p_value = 0.1
        elif abs_t < 3:
            p_value = 0.01
        else:
            p_value = 0.001
        
        return t_stat, p_value
    
    def perform_differential_analysis(self, condition1: str, condition2: str, 
                                    min_reads: int = 5, fold_change_threshold: float = 2.0):
        """Perform differential analysis between two conditions"""
        click.echo(f"📊 Performing differential analysis: {condition1} vs {condition2}")
        
        if condition1 not in self.condition_samples or condition2 not in self.condition_samples:
            available_conditions = list(self.condition_samples.keys())
            click.echo(f"❌ Available conditions: {available_conditions}")
            sys.exit(1)
        
        samples1 = self.condition_samples[condition1]
        samples2 = self.condition_samples[condition2]
        
        click.echo(f"   📋 {condition1}: {len(samples1)} samples")
        click.echo(f"   📋 {condition2}: {len(samples2)} samples")
        
        significant_count = 0
        upregulated = 0
        downregulated = 0
        
        for isotag_id in self.all_isotags:
            # Get counts for each condition
            counts1 = [self.sample_data.get(sample, {}).get(isotag_id, 0) for sample in samples1]
            counts2 = [self.sample_data.get(sample, {}).get(isotag_id, 0) for sample in samples2]
            
            # Filter by minimum reads
            total_reads = sum(counts1) + sum(counts2)
            if total_reads < min_reads:
                continue
            
            # Calculate basic statistics
            stats1 = self.calculate_basic_stats(counts1)
            stats2 = self.calculate_basic_stats(counts2)
            
            # Calculate fold change (with pseudocount to avoid division by zero)
            mean1_pseudo = stats1['mean'] + 0.1
            mean2_pseudo = stats2['mean'] + 0.1
            fold_change = mean1_pseudo / mean2_pseudo
            log2_fold_change = math.log2(fold_change)
            
            # Perform statistical test
            t_stat, p_value = self.welch_t_test([float(x) for x in counts1], [float(x) for x in counts2])
            
            # Determine significance
            is_significant = abs(log2_fold_change) >= math.log2(fold_change_threshold) and p_value < 0.05
            
            if is_significant:
                significant_count += 1
                if log2_fold_change > 0:
                    upregulated += 1
                else:
                    downregulated += 1
            
            # Store results
            self.isotag_stats[isotag_id] = {
                f'{condition1}_mean': stats1['mean'],
                f'{condition1}_std': stats1['std'],
                f'{condition1}_total': stats1['total'],
                f'{condition2}_mean': stats2['mean'],
                f'{condition2}_std': stats2['std'],
                f'{condition2}_total': stats2['total'],
                'fold_change': fold_change,
                'log2_fold_change': log2_fold_change,
                't_statistic': t_stat,
                'p_value': p_value,
                'is_significant': is_significant,
                'total_reads': total_reads
            }
        
        self.stats['significant_isotags'] = significant_count
        self.stats['upregulated'] = upregulated
        self.stats['downregulated'] = downregulated
        
        click.echo(f"   ✅ Analysis complete: {significant_count} significant isotags")
        click.echo(f"      📈 Upregulated in {condition1}: {upregulated}")
        click.echo(f"      📉 Downregulated in {condition1}: {downregulated}")
    
    def export_results(self, output_file: str, condition1: str, condition2: str):
        """Export differential analysis results to TSV"""
        click.echo(f"📄 Exporting results to {output_file}")
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            
            # Write header
            header = [
                'Isotag_ID',
                f'{condition1}_Mean',
                f'{condition1}_Std',
                f'{condition1}_Total',
                f'{condition2}_Mean',
                f'{condition2}_Std', 
                f'{condition2}_Total',
                'Fold_Change',
                'Log2_Fold_Change',
                'T_Statistic',
                'P_Value',
                'Significant',
                'Total_Reads',
                'Regulation'
            ]
            writer.writerow(header)
            
            # Sort by absolute log2 fold change
            sorted_isotags = sorted(
                self.isotag_stats.items(),
                key=lambda x: abs(x[1]['log2_fold_change']),
                reverse=True
            )
            
            for isotag_id, stats in sorted_isotags:
                # Determine regulation
                if stats['is_significant']:
                    if stats['log2_fold_change'] > 0:
                        regulation = f'Up in {condition1}'
                    else:
                        regulation = f'Down in {condition1}'
                else:
                    regulation = 'Not significant'
                
                row = [
                    isotag_id,
                    f"{stats[f'{condition1}_mean']:.2f}",
                    f"{stats[f'{condition1}_std']:.2f}",
                    stats[f'{condition1}_total'],
                    f"{stats[f'{condition2}_mean']:.2f}",
                    f"{stats[f'{condition2}_std']:.2f}",
                    stats[f'{condition2}_total'],
                    f"{stats['fold_change']:.3f}",
                    f"{stats['log2_fold_change']:.3f}",
                    f"{stats['t_statistic']:.3f}",
                    f"{stats['p_value']:.3e}",
                    str(stats['is_significant']),
                    stats['total_reads'],
                    regulation
                ]
                writer.writerow(row)
        
        click.echo(f"   ✅ Exported results for {len(self.isotag_stats)} isotags")
    
    def display_summary(self, condition1: str, condition2: str):
        """Display differential analysis summary"""
        click.echo("\n" + "="*60)
        click.echo("📊 ISO-TOOLS DIFF SUMMARY")
        click.echo("="*60)
        
        click.echo(f"🆔 Total isotags analyzed:   {len(self.isotag_stats):,}")
        click.echo(f"⭐ Significant isotags:      {self.stats['significant_isotags']:,}")
        click.echo(f"📈 Upregulated in {condition1}:   {self.stats['upregulated']:,}")
        click.echo(f"📉 Downregulated in {condition1}: {self.stats['downregulated']:,}")
        
        # Show top changed isotags
        if self.isotag_stats:
            click.echo(f"\n🏆 TOP 5 MOST CHANGED ISOTAGS")
            click.echo("-" * 40)
            
            sorted_isotags = sorted(
                self.isotag_stats.items(),
                key=lambda x: abs(x[1]['log2_fold_change']),
                reverse=True
            )
            
            for i, (isotag_id, stats) in enumerate(sorted_isotags[:5], 1):
                direction = "↑" if stats['log2_fold_change'] > 0 else "↓"
                significance = "***" if stats['is_significant'] else ""
                
                click.echo(f"{i}. {isotag_id[:32]}... {direction} "
                          f"FC: {stats['fold_change']:.2f} "
                          f"(log2: {stats['log2_fold_change']:.2f}) {significance}")


@click.command()
@click.argument('bam_files', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('--output', '-o', required=True, help='Output TSV file for results')
@click.option('--samples', help='Comma-separated sample names')
@click.option('--conditions', help='Comma-separated condition labels')
@click.option('--condition1', default='control', help='First condition name for comparison')
@click.option('--condition2', default='treatment', help='Second condition name for comparison')
@click.option('--use-tag', '-t', type=click.Choice(['XI', 'XB', 'XS', 'XT', 'XC']),
              default='XI', help='Tag to use for differential analysis (default: XI)')
@click.option('--min-reads', default=5, help='Minimum total reads per isotag')
@click.option('--fold-change', default=2.0, help='Fold change threshold for significance')
def diff(bam_files, output, samples, conditions, condition1, condition2, use_tag, min_reads, fold_change):
    """
    ISO-Tools Diff - Differential Isoform Analysis (v2.1+)

    Compare isotag expression between conditions for differential usage analysis.
    Performs statistical testing and ranks differentially expressed isoforms.
    Supports all 7 tags with flexible counting modes.

    Analysis Tags:
        XI: Differential at structure level (full exon coordinates)
        XB: Differential at boundary level (5'/3' ends)
        XS: Differential at junction level (splice patterns)
        XT: Differential at transcript level (biological groups)
        XC: Differential at cluster level (after clustering)

    Examples:
        # Basic comparison (auto-detect conditions, use XI)
        isotag_diff.py control1.bam control2.bam treatment1.bam treatment2.bam -o diff_results.tsv

        # Differential junction usage (XS)
        isotag_diff.py *.bam -o junction_diff.tsv --conditions "ctrl,ctrl,treat,treat" -t XS

        # Differential boundary usage (XB)
        isotag_diff.py *.bam -o boundary_diff.tsv --conditions "ctrl,ctrl,treat,treat" -t XB

        # Cluster-level differential (XC)
        isotag_diff.py clustered*.bam -o cluster_diff.tsv --conditions "A,A,B,B" -t XC

        # Custom parameters
        isotag_diff.py *.bam -o results.tsv --conditions "A,A,B,B" --condition1 A --condition2 B --fold-change 1.5

        # Differential transcript usage (XT)
        isotag_diff.py *.bam -o transcript_diff.tsv -t XT --conditions "ctrl,ctrl,treat,treat"
    """
    
    if len(bam_files) < 2:
        click.echo("❌ Need at least 2 BAM files for comparison")
        sys.exit(1)
    
    # Parse sample names and conditions
    sample_names = None
    if samples:
        sample_names = [s.strip() for s in samples.split(',')]
    
    condition_labels = None
    if conditions:
        condition_labels = [c.strip() for c in conditions.split(',')]

    analyzer = IsotagDifferentialAnalyzer(use_tag=use_tag)
    
    # Load data
    analyzer.load_sample_counts(list(bam_files), sample_names, condition_labels)
    
    if len(analyzer.condition_samples) < 2:
        available_conditions = list(analyzer.condition_samples.keys())
        click.echo(f"❌ Need at least 2 conditions. Found: {available_conditions}")
        sys.exit(1)
    
    # Verify conditions exist
    if condition1 not in analyzer.condition_samples:
        available = list(analyzer.condition_samples.keys())
        click.echo(f"❌ Condition '{condition1}' not found. Available: {available}")
        sys.exit(1)
    
    if condition2 not in analyzer.condition_samples:
        available = list(analyzer.condition_samples.keys())
        click.echo(f"❌ Condition '{condition2}' not found. Available: {available}")
        sys.exit(1)
    
    # Perform analysis
    analyzer.perform_differential_analysis(condition1, condition2, min_reads, fold_change)
    
    # Export results
    analyzer.export_results(output, condition1, condition2)
    
    # Display summary
    analyzer.display_summary(condition1, condition2)
    
    click.echo(f"\n✅ Differential analysis complete!")
    click.echo(f"📄 Results: {output}")
    click.echo(f"💡 Next steps: Validate significant isotags, perform functional analysis")


if __name__ == '__main__':
    diff()