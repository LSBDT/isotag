#!/usr/bin/env python3
"""
ISO-Tools Compare - Compare Isotags Between Multiple BAM Files

Compare isotag overlap and differences between BAM files.
Generate Venn diagram data, intersection statistics, and differential analysis.
"""

import subprocess
import click
import json
import csv
from collections import Counter, defaultdict
from pathlib import Path
import sys


class IsotogComparer:
    """Compare isotags between multiple BAM files (v2.1+ - all 7 tags)"""

    def __init__(self):
        self.file_data = {}  # filename -> {all tag sets and counts}
        self.comparison_results = {}
    
    def extract_isotags_from_bam(self, bam_file, sample_name=None):
        """Extract all isotags from a BAM file"""
        if sample_name is None:
            sample_name = Path(bam_file).stem
        
        click.echo(f"🔍 Extracting isotags from {sample_name}...")

        # v1.0 tags
        xi_tags = set()
        xv_tags = set()
        xi_counts = Counter()
        xv_counts = Counter()
        # v2.0+ tags
        xb_tags = set()
        xs_tags = set()
        xt_tags = set()
        xc_tags = set()
        xr_tags = set()
        xb_counts = Counter()
        xs_counts = Counter()
        xt_counts = Counter()
        xc_counts = Counter()
        xr_counts = Counter()

        total_reads = 0
        reads_with_tags = 0
        
        view_cmd = ['samtools', 'view', bam_file]
        
        try:
            process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)
            
            for line in process.stdout:
                total_reads += 1
                
                if total_reads % 10000 == 0:
                    click.echo(f"   ⏳ Processed {total_reads:,} reads...")
                
                fields = line.strip().split('\t')
                
                if len(fields) >= 11:
                    has_isotag = False

                    for field in fields[11:]:
                        if field.startswith('XI:Z:'):
                            xi_tag = field[5:]
                            xi_tags.add(xi_tag)
                            xi_counts[xi_tag] += 1
                            has_isotag = True
                        elif field.startswith('XB:Z:'):
                            xb_tag = field[5:]
                            xb_tags.add(xb_tag)
                            xb_counts[xb_tag] += 1
                            has_isotag = True
                        elif field.startswith('XS:Z:'):
                            xs_tag = field[5:]
                            xs_tags.add(xs_tag)
                            xs_counts[xs_tag] += 1
                            has_isotag = True
                        elif field.startswith('XT:Z:'):
                            xt_tag = field[5:]
                            xt_tags.add(xt_tag)
                            xt_counts[xt_tag] += 1
                            has_isotag = True
                        elif field.startswith('XV:Z:'):
                            xv_tag = field[5:]
                            xv_tags.add(xv_tag)
                            xv_counts[xv_tag] += 1
                            has_isotag = True
                        elif field.startswith('XC:Z:'):
                            xc_tag = field[5:]
                            xc_tags.add(xc_tag)
                            xc_counts[xc_tag] += 1
                            has_isotag = True
                        elif field.startswith('XR:Z:'):
                            xr_tag = field[5:]
                            xr_tags.add(xr_tag)
                            xr_counts[xr_tag] += 1
                            has_isotag = True

                    if has_isotag:
                        reads_with_tags += 1
            
            process.wait()
            
            self.file_data[sample_name] = {
                'xi_tags': xi_tags,
                'xb_tags': xb_tags,
                'xs_tags': xs_tags,
                'xt_tags': xt_tags,
                'xv_tags': xv_tags,
                'xc_tags': xc_tags,
                'xr_tags': xr_tags,
                'xi_counts': xi_counts,
                'xb_counts': xb_counts,
                'xs_counts': xs_counts,
                'xt_counts': xt_counts,
                'xv_counts': xv_counts,
                'xc_counts': xc_counts,
                'xr_counts': xr_counts,
                'total_reads': total_reads,
                'reads_with_tags': reads_with_tags
            }

            click.echo(f"   ✅ {sample_name}: XI={len(xi_tags)}, XB={len(xb_tags)}, XS={len(xs_tags)}, XT={len(xt_tags)}, XV={len(xv_tags)}, XC={len(xc_tags)}, XR={len(xr_tags)} from {reads_with_tags:,}/{total_reads:,} reads")
            
        except subprocess.CalledProcessError as e:
            click.echo(f"❌ Error reading {bam_file}: {e}")
            sys.exit(1)
    
    def compute_pairwise_comparisons(self):
        """Compute pairwise comparisons between all samples"""
        sample_names = list(self.file_data.keys())
        
        self.comparison_results = {
            'samples': sample_names,
            'pairwise': {},
            'overall': {
                'all_structures': set(),
                'all_boundaries': set(),
                'all_splicetags': set(),
                'all_transcripts': set(),
                'all_variants': set(),
                'all_clusters': set(),
                'all_representatives': set(),
                'intersection_matrix': {}
            }
        }

        # Collect all unique isotags (v2.1+ all tags)
        for sample_data in self.file_data.values():
            self.comparison_results['overall']['all_structures'].update(sample_data['xi_tags'])
            self.comparison_results['overall']['all_boundaries'].update(sample_data['xb_tags'])
            self.comparison_results['overall']['all_splicetags'].update(sample_data['xs_tags'])
            self.comparison_results['overall']['all_transcripts'].update(sample_data['xt_tags'])
            self.comparison_results['overall']['all_variants'].update(sample_data['xv_tags'])
            self.comparison_results['overall']['all_clusters'].update(sample_data['xc_tags'])
            self.comparison_results['overall']['all_representatives'].update(sample_data['xr_tags'])
        
        # Pairwise comparisons
        for i, sample1 in enumerate(sample_names):
            for j, sample2 in enumerate(sample_names[i+1:], i+1):
                self._compare_pair(sample1, sample2)
        
        # Convert sets to lists for JSON serialization
        self.comparison_results['overall']['all_structures'] = list(self.comparison_results['overall']['all_structures'])
        self.comparison_results['overall']['all_boundaries'] = list(self.comparison_results['overall']['all_boundaries'])
        self.comparison_results['overall']['all_splicetags'] = list(self.comparison_results['overall']['all_splicetags'])
        self.comparison_results['overall']['all_transcripts'] = list(self.comparison_results['overall']['all_transcripts'])
        self.comparison_results['overall']['all_variants'] = list(self.comparison_results['overall']['all_variants'])
        self.comparison_results['overall']['all_clusters'] = list(self.comparison_results['overall']['all_clusters'])
        self.comparison_results['overall']['all_representatives'] = list(self.comparison_results['overall']['all_representatives'])
    
    def _compare_pair(self, sample1, sample2):
        """Compare two samples"""
        data1 = self.file_data[sample1]
        data2 = self.file_data[sample2]
        
        # Structure comparison
        xi_intersection = data1['xi_tags'] & data2['xi_tags']
        xi_union = data1['xi_tags'] | data2['xi_tags']
        xi_only1 = data1['xi_tags'] - data2['xi_tags']
        xi_only2 = data2['xi_tags'] - data1['xi_tags']
        
        # Variant comparison
        xv_intersection = data1['xv_tags'] & data2['xv_tags']
        xv_union = data1['xv_tags'] | data2['xv_tags']
        xv_only1 = data1['xv_tags'] - data2['xv_tags']
        xv_only2 = data2['xv_tags'] - data1['xv_tags']
        
        # Jaccard similarity
        xi_jaccard = len(xi_intersection) / len(xi_union) if xi_union else 0
        xv_jaccard = len(xv_intersection) / len(xv_union) if xv_union else 0
        
        pair_key = f"{sample1}_vs_{sample2}"
        
        self.comparison_results['pairwise'][pair_key] = {
            'sample1': sample1,
            'sample2': sample2,
            'structures': {
                'intersection': list(xi_intersection),
                'union': list(xi_union),
                'only_sample1': list(xi_only1),
                'only_sample2': list(xi_only2),
                'jaccard_similarity': xi_jaccard,
                'counts': {
                    'intersection': len(xi_intersection),
                    'sample1_total': len(data1['xi_tags']),
                    'sample2_total': len(data2['xi_tags']),
                    'union': len(xi_union)
                }
            },
            'variants': {
                'intersection': list(xv_intersection),
                'union': list(xv_union),
                'only_sample1': list(xv_only1),
                'only_sample2': list(xv_only2),
                'jaccard_similarity': xv_jaccard,
                'counts': {
                    'intersection': len(xv_intersection),
                    'sample1_total': len(data1['xv_tags']),
                    'sample2_total': len(data2['xv_tags']),
                    'union': len(xv_union)
                }
            }
        }
    
    def display_summary(self):
        """Display comparison summary"""
        click.echo("\n" + "="*60)
        click.echo("🔍 ISO-TOOLS COMPARE SUMMARY")
        click.echo("="*60)
        
        sample_names = self.comparison_results['samples']
        
        # Sample overview (v2.1+ all tags)
        click.echo(f"📂 Samples analyzed: {len(sample_names)}")
        for name in sample_names:
            data = self.file_data[name]
            click.echo(f"   • {name}:")
            click.echo(f"     XI={len(data['xi_tags'])}, XB={len(data['xb_tags'])}, XS={len(data['xs_tags'])}, XT={len(data['xt_tags'])}")
            click.echo(f"     XV={len(data['xv_tags'])}, XC={len(data['xc_tags'])}, XR={len(data['xr_tags'])}")

        # Overall statistics
        overall = self.comparison_results['overall']

        click.echo(f"\n🆔 Total Unique Tags (v2.1+):")
        click.echo(f"   XI (Structures):      {len(overall['all_structures'])}")
        click.echo(f"   XB (Boundaries):      {len(overall['all_boundaries'])}")
        click.echo(f"   XS (Splicetags):      {len(overall['all_splicetags'])}")
        click.echo(f"   XT (Transcripts):     {len(overall['all_transcripts'])}")
        click.echo(f"   XV (Variants):        {len(overall['all_variants'])}")
        click.echo(f"   XC (Clusters):        {len(overall['all_clusters'])}")
        click.echo(f"   XR (Representatives): {len(overall['all_representatives'])}")
        
        # Pairwise comparisons (show first few)
        click.echo(f"\n🔗 PAIRWISE COMPARISONS")
        click.echo("-" * 40)
        
        pairwise_items = list(self.comparison_results['pairwise'].items())
        for pair_key, comparison in pairwise_items[:5]:  # Show first 5 comparisons
            sample1 = comparison['sample1'] 
            sample2 = comparison['sample2']
            
            xi_jaccard = comparison['structures']['jaccard_similarity']
            xv_jaccard = comparison['variants']['jaccard_similarity']
            xi_shared = comparison['structures']['counts']['intersection']
            xv_shared = comparison['variants']['counts']['intersection']
            
            click.echo(f"📊 {sample1} vs {sample2}:")
            click.echo(f"   Structures: {xi_shared} shared (Jaccard: {xi_jaccard:.3f})")
            click.echo(f"   Variants:   {xv_shared} shared (Jaccard: {xv_jaccard:.3f})")
            click.echo()
        
        if len(pairwise_items) > 5:
            click.echo(f"   ... and {len(pairwise_items) - 5} more comparisons")
    
    def export_results(self, output_file, format='json'):
        """Export comparison results"""
        output_path = Path(output_file)
        
        if format == 'json':
            # Convert Counter objects to dicts for JSON serialization
            export_data = dict(self.comparison_results)
            
            # Add sample details
            export_data['sample_details'] = {}
            for sample_name, data in self.file_data.items():
                export_data['sample_details'][sample_name] = {
                    'total_reads': data['total_reads'],
                    'reads_with_tags': data['reads_with_tags'],
                    'unique_tags': {
                        'XI': len(data['xi_tags']),
                        'XB': len(data['xb_tags']),
                        'XS': len(data['xs_tags']),
                        'XT': len(data['xt_tags']),
                        'XV': len(data['xv_tags']),
                        'XC': len(data['xc_tags']),
                        'XR': len(data['xr_tags'])
                    },
                    'top_tags': {
                        'XI': dict(data['xi_counts'].most_common(10)),
                        'XB': dict(data['xb_counts'].most_common(10)),
                        'XS': dict(data['xs_counts'].most_common(10)),
                        'XT': dict(data['xt_counts'].most_common(10)),
                        'XV': dict(data['xv_counts'].most_common(10)),
                        'XC': dict(data['xc_counts'].most_common(10)),
                        'XR': dict(data['xr_counts'].most_common(10))
                    }
                }
            
            with open(output_path, 'w') as f:
                json.dump(export_data, f, indent=2)
            
            click.echo(f"💾 Comparison results exported to JSON: {output_file}")
        
        elif format == 'tsv':
            with open(output_path, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                
                # Write header
                writer.writerow([
                    'Sample1', 'Sample2', 'Structures_Shared', 'Structures_Sample1_Only',
                    'Structures_Sample2_Only', 'Structures_Jaccard',
                    'Variants_Shared', 'Variants_Sample1_Only', 'Variants_Sample2_Only',
                    'Variants_Jaccard'
                ])
                
                # Write pairwise comparisons
                for comparison in self.comparison_results['pairwise'].values():
                    writer.writerow([
                        comparison['sample1'],
                        comparison['sample2'],
                        comparison['structures']['counts']['intersection'],
                        len(comparison['structures']['only_sample1']),
                        len(comparison['structures']['only_sample2']),
                        f"{comparison['structures']['jaccard_similarity']:.4f}",
                        comparison['variants']['counts']['intersection'],
                        len(comparison['variants']['only_sample1']),
                        len(comparison['variants']['only_sample2']),
                        f"{comparison['variants']['jaccard_similarity']:.4f}"
                    ])
            
            click.echo(f"📋 Comparison results exported to TSV: {output_file}")


@click.command()
@click.argument('bam_files', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('--output', '-o', help='Output file for comparison results')
@click.option('--format', '-f', type=click.Choice(['json', 'tsv', 'auto']), default='auto', help='Output format')
@click.option('--names', help='Comma-separated sample names (default: use filenames)')
def compare(bam_files, output, format, names):
    """
    ISO-Tools Compare - Compare Isotags Between BAM Files
    
    Compare isotag overlap and differences between multiple BAM files.
    Generates intersection statistics and similarity metrics.
    
    Examples:
        # Compare two BAM files
        iso-tools compare sample1.bam sample2.bam
        
        # Compare multiple files with custom names
        iso-tools compare *.bam --names "Control,Treatment1,Treatment2" -o comparison.json
        
        # Export to TSV format
        iso-tools compare sample1.bam sample2.bam -o results.tsv -f tsv
    """
    
    if len(bam_files) < 2:
        click.echo("❌ At least 2 BAM files are required for comparison")
        sys.exit(1)
    
    # Parse sample names
    sample_names = None
    if names:
        sample_names = [name.strip() for name in names.split(',')]
        if len(sample_names) != len(bam_files):
            click.echo(f"❌ Number of sample names ({len(sample_names)}) doesn't match number of files ({len(bam_files)})")
            sys.exit(1)
    
    comparer = IsotogComparer()
    
    # Extract isotags from each file
    for i, bam_file in enumerate(bam_files):
        sample_name = sample_names[i] if sample_names else None
        comparer.extract_isotags_from_bam(bam_file, sample_name)
    
    # Compute comparisons
    comparer.compute_pairwise_comparisons()
    
    # Display results
    comparer.display_summary()
    
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
        
        comparer.export_results(output, format)
    
    click.echo(f"\n✅ Comparison complete!")


if __name__ == '__main__':
    compare()