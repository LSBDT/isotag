#!/usr/bin/env python3
"""
ISO-Tools Count - Generate Expression Matrices with Isotag IDs

Generate expression count matrices using isotag IDs as features.
Supports both bulk and single-cell RNA-seq data with barcode support.
"""

import subprocess
import click
import csv
import json
import gzip
from collections import Counter, defaultdict
from pathlib import Path
import sys
import re


class IsotogCounter:
    """Generate expression matrices from isotag-annotated BAM files"""
    
    def __init__(self):
        self.isotag_counts = defaultdict(Counter)  # isotag_id -> {sample/barcode: count}
        self.sample_stats = {}  # sample -> {total_reads, tagged_reads, unique_isotags}
        self.all_isotags = set()
        self.all_samples = set()
    
    def extract_barcode_from_read(self, fields, barcode_tag='CB'):
        """Extract cell barcode from BAM fields"""
        for field in fields[11:]:
            if field.startswith(f'{barcode_tag}:Z:'):
                return field[5:]
        return None
    
    def process_bam_file(self, bam_file, sample_name=None, mode='bulk', barcode_tag='CB', umi_tag='UB', count_by='XI'):
        """Process BAM file and count isotags (v2.1+ supports all tags)"""
        if sample_name is None:
            sample_name = Path(bam_file).stem

        click.echo(f"🔢 Counting isotags in {sample_name} ({mode} mode, counting by {count_by})...")
        
        total_reads = 0
        tagged_reads = 0
        isotag_sample_counts = Counter()
        
        view_cmd = ['samtools', 'view', bam_file]
        
        try:
            process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)
            
            for line in process.stdout:
                total_reads += 1
                
                if total_reads % 10000 == 0:
                    click.echo(f"   ⏳ Processed {total_reads:,} reads...")
                
                fields = line.strip().split('\t')
                
                if len(fields) >= 11:
                    # Extract all tags (v2.1+)
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

                    # Select tag to count by
                    isotag_id = tags.get(count_by)

                    if isotag_id:
                        tagged_reads += 1
                        self.all_isotags.add(isotag_id)
                        
                        if mode == 'bulk':
                            # Bulk RNA-seq: count by sample
                            count_key = sample_name
                            self.isotag_counts[isotag_id][count_key] += 1
                            isotag_sample_counts[isotag_id] += 1
                            
                        elif mode == 'single-cell':
                            # Single-cell: count by cell barcode
                            barcode = self.extract_barcode_from_read(fields, barcode_tag)
                            if barcode:
                                count_key = f"{sample_name}_{barcode}"
                                self.isotag_counts[isotag_id][count_key] += 1
                                isotag_sample_counts[isotag_id] += 1
                                self.all_samples.add(count_key)
                            
                        elif mode == 'umi':
                            # UMI-based counting (collapse by UMI)
                            barcode = self.extract_barcode_from_read(fields, barcode_tag)
                            umi = self.extract_barcode_from_read(fields, umi_tag)
                            
                            if barcode and umi:
                                # For UMI mode, we'd need to collapse identical UMIs
                                # This is a simplified version - would need UMI deduplication
                                count_key = f"{sample_name}_{barcode}"
                                umi_key = f"{isotag_id}_{barcode}_{umi}"
                                # In practice, you'd store UMI info and deduplicate later
                                self.isotag_counts[isotag_id][count_key] += 1
                                self.all_samples.add(count_key)
            
            process.wait()
            
            # Store sample statistics
            if mode == 'bulk':
                self.all_samples.add(sample_name)
            
            self.sample_stats[sample_name] = {
                'total_reads': total_reads,
                'tagged_reads': tagged_reads,
                'unique_isotags': len(isotag_sample_counts),
                'tagging_rate': tagged_reads / total_reads if total_reads > 0 else 0
            }
            
            click.echo(f"   ✅ {sample_name}: {tagged_reads:,}/{total_reads:,} tagged reads, {len(isotag_sample_counts)} unique isotags")
            
        except subprocess.CalledProcessError as e:
            click.echo(f"❌ Error reading {bam_file}: {e}")
            sys.exit(1)
    
    def export_matrix_tsv(self, output_file):
        """Export count matrix as TSV"""
        click.echo(f"📊 Exporting count matrix to {output_file}")
        
        # Sort samples and isotags for consistent output
        sorted_samples = sorted(self.all_samples)
        sorted_isotags = sorted(self.all_isotags)
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            
            # Write header
            header = ['Isotag_ID'] + sorted_samples
            writer.writerow(header)
            
            # Write count matrix
            for isotag_id in sorted_isotags:
                row = [isotag_id]
                for sample in sorted_samples:
                    count = self.isotag_counts[isotag_id].get(sample, 0)
                    row.append(count)
                writer.writerow(row)
    
    def export_matrix_csv(self, output_file):
        """Export count matrix as CSV"""
        click.echo(f"📊 Exporting count matrix to {output_file}")
        
        sorted_samples = sorted(self.all_samples)
        sorted_isotags = sorted(self.all_isotags)
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Write header
            header = ['Isotag_ID'] + sorted_samples
            writer.writerow(header)
            
            # Write count matrix
            for isotag_id in sorted_isotags:
                row = [isotag_id]
                for sample in sorted_samples:
                    count = self.isotag_counts[isotag_id].get(sample, 0)
                    row.append(count)
                writer.writerow(row)
    
    def export_matrix_mtx(self, output_file):
        """Export count matrix in Matrix Market format (for single-cell analysis)"""
        click.echo(f"📊 Exporting count matrix to Matrix Market format: {output_file}")
        
        sorted_samples = sorted(self.all_samples)
        sorted_isotags = sorted(self.all_isotags)
        
        # Count non-zero entries
        nnz = sum(1 for isotag_counts in self.isotag_counts.values() 
                 for count in isotag_counts.values() if count > 0)
        
        output_path = Path(output_file)
        
        # Write matrix file
        with open(output_path, 'w') as f:
            # Header
            f.write("%%MatrixMarket matrix coordinate integer general\n")
            f.write(f"{len(sorted_isotags)} {len(sorted_samples)} {nnz}\n")
            
            # Data (1-indexed)
            for isotag_idx, isotag_id in enumerate(sorted_isotags, 1):
                for sample_idx, sample in enumerate(sorted_samples, 1):
                    count = self.isotag_counts[isotag_id].get(sample, 0)
                    if count > 0:
                        f.write(f"{isotag_idx} {sample_idx} {count}\n")
        
        # Write features file (isotag IDs)
        features_file = output_path.parent / f"{output_path.stem}_features.tsv"
        with open(features_file, 'w') as f:
            for isotag_id in sorted_isotags:
                f.write(f"{isotag_id}\t{isotag_id}\tIsoform\n")
        
        # Write barcodes file (sample/cell IDs)
        barcodes_file = output_path.parent / f"{output_path.stem}_barcodes.tsv"
        with open(barcodes_file, 'w') as f:
            for sample in sorted_samples:
                f.write(f"{sample}\n")
        
        click.echo(f"   📄 Features: {features_file}")
        click.echo(f"   📄 Barcodes: {barcodes_file}")
    
    def export_summary_json(self, output_file):
        """Export summary statistics as JSON"""
        summary = {
            'total_isotags': len(self.all_isotags),
            'total_samples': len(self.all_samples),
            'sample_statistics': self.sample_stats,
            'top_isotags': {},
            'matrix_dimensions': {
                'features': len(self.all_isotags),
                'samples': len(self.all_samples)
            }
        }
        
        # Calculate top isotags across all samples
        total_counts = Counter()
        for isotag_id, sample_counts in self.isotag_counts.items():
            total_counts[isotag_id] = sum(sample_counts.values())
        
        summary['top_isotags'] = dict(total_counts.most_common(20))
        
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        click.echo(f"📋 Summary exported to {output_file}")
    
    def display_summary(self):
        """Display counting summary"""
        click.echo("\n" + "="*60)
        click.echo("🔢 ISO-TOOLS COUNT SUMMARY")
        click.echo("="*60)
        
        click.echo(f"🆔 Total unique isotags:  {len(self.all_isotags)}")
        click.echo(f"📂 Total samples/cells:   {len(self.all_samples)}")
        
        # Sample statistics
        click.echo(f"\n📊 SAMPLE STATISTICS")
        click.echo("-" * 40)
        for sample, stats in self.sample_stats.items():
            click.echo(f"• {sample}:")
            click.echo(f"   Reads: {stats['tagged_reads']:,}/{stats['total_reads']:,} ({stats['tagging_rate']:.1%})")
            click.echo(f"   Unique isotags: {stats['unique_isotags']:,}")
        
        # Top isotags
        total_counts = Counter()
        for isotag_id, sample_counts in self.isotag_counts.items():
            total_counts[isotag_id] = sum(sample_counts.values())
        
        click.echo(f"\n🏆 TOP 5 ISOTAGS")
        click.echo("-" * 40)
        for i, (isotag_id, count) in enumerate(total_counts.most_common(5), 1):
            click.echo(f"{i}. {isotag_id[:32]}... {count:,} reads")


@click.command()
@click.argument('bam_files', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('--output', '-o', required=True, help='Output file prefix')
@click.option('--format', '-f', type=click.Choice(['tsv', 'csv', 'mtx', 'all']), default='tsv', help='Output format')
@click.option('--mode', '-m', type=click.Choice(['bulk', 'single-cell', 'umi']), default='bulk', help='Analysis mode')
@click.option('--count-by', type=click.Choice(['XI', 'XB', 'XS', 'XT', 'XC']), default='XI',
              help='Tag to count by: XI=structure, XB=boundary, XS=splicetag, XT=transcript, XC=cluster (default: XI)')
@click.option('--barcode-tag', default='CB', help='Cell barcode tag (default: CB)')
@click.option('--umi-tag', default='UB', help='UMI tag (default: UB)')
@click.option('--names', help='Comma-separated sample names (default: use filenames)')
def count(bam_files, output, format, mode, count_by, barcode_tag, umi_tag, names):
    """
    ISO-Tools Count - Generate Expression Matrices with Isotag IDs (v2.1+)

    Create count matrices using isotag IDs as features. Supports counting by:
    - XI (structure): Full exon structure (default)
    - XB (boundary): 5'/3' boundaries only
    - XS (splicetag): Splice junctions only
    - XT (transcript): Biological transcript groups
    - XC (cluster): Cluster assignments

    Supports bulk RNA-seq, single-cell, and UMI-based counting.

    Examples:
        # Count by structure (default)
        iso-tools count sample1.bam sample2.bam -o expression -f tsv

        # Count by splicetag (junction level)
        iso-tools count sample1.bam -o junction_matrix -f tsv --count-by XS

        # Count by transcript group
        iso-tools count sample1.bam -o transcript_matrix -f tsv --count-by XT

        # Count by cluster (after clustering)
        iso-tools count clustered.bam -o cluster_matrix -f tsv --count-by XC

        # Single-cell analysis (10x format)
        iso-tools count cells.bam -o sc_matrix -f mtx -m single-cell

        # Multiple formats
        iso-tools count *.bam -o results -f all --names "Ctrl,Treat1,Treat2"
    """
    
    # Parse sample names
    sample_names = None
    if names:
        sample_names = [name.strip() for name in names.split(',')]
        if len(sample_names) != len(bam_files):
            click.echo(f"❌ Number of sample names ({len(sample_names)}) doesn't match number of files ({len(bam_files)})")
            sys.exit(1)
    
    counter = IsotogCounter()
    
    # Process each BAM file
    for i, bam_file in enumerate(bam_files):
        sample_name = sample_names[i] if sample_names else None
        counter.process_bam_file(bam_file, sample_name, mode, barcode_tag, umi_tag, count_by)
    
    # Display summary
    counter.display_summary()
    
    # Export results
    output_path = Path(output)
    
    if format == 'tsv' or format == 'all':
        counter.export_matrix_tsv(f"{output_path}.tsv")
    
    if format == 'csv' or format == 'all':
        counter.export_matrix_csv(f"{output_path}.csv")
    
    if format == 'mtx' or format == 'all':
        counter.export_matrix_mtx(f"{output_path}.mtx")
    
    # Always export summary
    counter.export_summary_json(f"{output_path}_summary.json")
    
    click.echo(f"\n✅ Count matrix generation complete!")


if __name__ == '__main__':
    count()