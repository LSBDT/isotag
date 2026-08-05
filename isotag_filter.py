#!/usr/bin/env python3
"""
ISO-Tools Filter - Filter BAM Files by Isotag IDs

Filter BAM files to include/exclude specific isotag IDs.
Useful for focused analysis of specific isoforms or variants.
"""

import subprocess
import click
import sys
import tempfile
import os
from pathlib import Path
from typing import Set, Optional


class IsotogFilter:
    """Filter BAM files by isotag IDs (v2.1+ - all 7 tags)"""

    def __init__(self):
        # v1.0 tags
        self.include_xi = set()  # Structure IDs to include
        self.exclude_xi = set()  # Structure IDs to exclude
        self.include_xv = set()  # Variant IDs to include
        self.exclude_xv = set()  # Variant IDs to exclude
        # v2.0+ tags
        self.include_xb = set()  # Boundary tags to include
        self.exclude_xb = set()  # Boundary tags to exclude
        self.include_xs = set()  # Splicetags to include
        self.exclude_xs = set()  # Splicetags to exclude
        self.include_xt = set()  # Transcript groups to include
        self.exclude_xt = set()  # Transcript groups to exclude
        self.include_xc = set()  # Cluster IDs to include
        self.exclude_xc = set()  # Cluster IDs to exclude
        self.include_xr = set()  # Representative tags to include
        self.exclude_xr = set()  # Representative tags to exclude

        self.stats = {
            'total_reads': 0,
            'passed_reads': 0,
            'filtered_reads': 0,
            'reads_with_xi': 0,
            'reads_with_xb': 0,
            'reads_with_xs': 0,
            'reads_with_xt': 0,
            'reads_with_xv': 0,
            'reads_with_xc': 0,
            'reads_with_xr': 0
        }
    
    def load_isotag_list(self, file_path: str) -> Set[str]:
        """Load isotag IDs from file (one per line)"""
        isotag_set = set()
        
        with open(file_path, 'r') as f:
            for line in f:
                isotag_id = line.strip()
                if isotag_id and not isotag_id.startswith('#'):
                    isotag_set.add(isotag_id)
        
        return isotag_set
    
    def should_keep_read(self, tags: dict) -> bool:
        """Determine if a read should be kept based on filter criteria (v2.1+ all tags)"""

        # Extract all tags
        xi_tag = tags.get('XI')
        xb_tag = tags.get('XB')
        xs_tag = tags.get('XS')
        xt_tag = tags.get('XT')
        xv_tag = tags.get('XV')
        xc_tag = tags.get('XC')
        xr_tag = tags.get('XR')

        # Check each tag type
        passes = []

        # XI (structure) filters
        if self.include_xi or self.exclude_xi:
            xi_pass = True
            if self.include_xi:
                xi_pass = xi_tag in self.include_xi
            if self.exclude_xi and xi_tag:
                xi_pass = xi_pass and (xi_tag not in self.exclude_xi)
            passes.append(xi_pass)

        # XB (boundary) filters
        if self.include_xb or self.exclude_xb:
            xb_pass = True
            if self.include_xb:
                xb_pass = xb_tag in self.include_xb
            if self.exclude_xb and xb_tag:
                xb_pass = xb_pass and (xb_tag not in self.exclude_xb)
            passes.append(xb_pass)

        # XS (splicetag) filters
        if self.include_xs or self.exclude_xs:
            xs_pass = True
            if self.include_xs:
                xs_pass = xs_tag in self.include_xs
            if self.exclude_xs and xs_tag:
                xs_pass = xs_pass and (xs_tag not in self.exclude_xs)
            passes.append(xs_pass)

        # XT (transcript group) filters
        if self.include_xt or self.exclude_xt:
            xt_pass = True
            if self.include_xt:
                xt_pass = xt_tag in self.include_xt
            if self.exclude_xt and xt_tag:
                xt_pass = xt_pass and (xt_tag not in self.exclude_xt)
            passes.append(xt_pass)

        # XV (variant) filters
        if self.include_xv or self.exclude_xv:
            xv_pass = True
            if self.include_xv:
                xv_pass = xv_tag in self.include_xv
            if self.exclude_xv and xv_tag:
                xv_pass = xv_pass and (xv_tag not in self.exclude_xv)
            passes.append(xv_pass)

        # XC (cluster) filters
        if self.include_xc or self.exclude_xc:
            xc_pass = True
            if self.include_xc:
                xc_pass = xc_tag in self.include_xc
            if self.exclude_xc and xc_tag:
                xc_pass = xc_pass and (xc_tag not in self.exclude_xc)
            passes.append(xc_pass)

        # XR (representative) filters
        if self.include_xr or self.exclude_xr:
            xr_pass = True
            if self.include_xr:
                xr_pass = xr_tag in self.include_xr
            if self.exclude_xr and xr_tag:
                xr_pass = xr_pass and (xr_tag not in self.exclude_xr)
            passes.append(xr_pass)

        # All filters must pass
        return all(passes) if passes else True
    
    def filter_bam(self, input_bam: str, output_bam: str, keep_untagged: bool = False):
        """Filter BAM file based on isotag criteria"""
        
        click.echo(f"🔍 Filtering {input_bam} -> {output_bam}")
        click.echo(f"   Include XI: {len(self.include_xi)} IDs")
        click.echo(f"   Exclude XI: {len(self.exclude_xi)} IDs") 
        click.echo(f"   Include XV: {len(self.include_xv)} IDs")
        click.echo(f"   Exclude XV: {len(self.exclude_xv)} IDs")
        click.echo(f"   Keep untagged: {keep_untagged}")
        
        input_path = Path(input_bam)
        output_path = Path(output_bam)
        
        # Determine input/output formats
        input_is_bam = input_path.suffix.lower() == '.bam'
        output_is_bam = output_path.suffix.lower() == '.bam'
        
        try:
            # Create temporary SAM file for processing
            with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as temp_output:
                temp_output_path = temp_output.name
                
                # Read input BAM/SAM
                if input_is_bam:
                    samtools_cmd = ['samtools', 'view', '-h', input_bam]
                else:
                    samtools_cmd = ['cat', input_bam]
                
                with subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE, text=True) as proc:
                    for line in proc.stdout:
                        # Pass through header lines
                        if line.startswith('@'):
                            temp_output.write(line)
                            continue
                        
                        self.stats['total_reads'] += 1
                        
                        if self.stats['total_reads'] % 10000 == 0:
                            click.echo(f"   ⏳ Processed {self.stats['total_reads']:,} reads...")
                        
                        # Parse read line
                        fields = line.strip().split('\t')

                        if len(fields) >= 11:
                            # Extract all isotag tags (v2.1+)
                            tags = {}

                            for field in fields[11:]:
                                if field.startswith('XI:Z:'):
                                    tags['XI'] = field[5:]
                                    self.stats['reads_with_xi'] += 1
                                elif field.startswith('XB:Z:'):
                                    tags['XB'] = field[5:]
                                    self.stats['reads_with_xb'] += 1
                                elif field.startswith('XS:Z:'):
                                    tags['XS'] = field[5:]
                                    self.stats['reads_with_xs'] += 1
                                elif field.startswith('XT:Z:'):
                                    tags['XT'] = field[5:]
                                    self.stats['reads_with_xt'] += 1
                                elif field.startswith('XV:Z:'):
                                    tags['XV'] = field[5:]
                                    self.stats['reads_with_xv'] += 1
                                elif field.startswith('XC:Z:'):
                                    tags['XC'] = field[5:]
                                    self.stats['reads_with_xc'] += 1
                                elif field.startswith('XR:Z:'):
                                    tags['XR'] = field[5:]
                                    self.stats['reads_with_xr'] += 1

                            # Apply filtering logic
                            should_keep = False

                            if tags:
                                # Read has isotags - apply filters
                                should_keep = self.should_keep_read(tags)
                            else:
                                # Untagged read - keep if specified
                                should_keep = keep_untagged
                            
                            if should_keep:
                                temp_output.write(line)
                                self.stats['passed_reads'] += 1
                            else:
                                self.stats['filtered_reads'] += 1
                        else:
                            # Malformed line - keep if untagged reads are kept
                            if keep_untagged:
                                temp_output.write(line)
                                self.stats['passed_reads'] += 1
                            else:
                                self.stats['filtered_reads'] += 1
            
            # Convert to final output format
            if output_is_bam:
                click.echo("🔄 Converting to BAM format...")
                subprocess.run(['samtools', 'view', '-b', temp_output_path, '-o', str(output_path)], 
                             check=True)
                os.unlink(temp_output_path)
            else:
                # Output is SAM, just move the file
                import shutil
                shutil.move(temp_output_path, str(output_path))
            
            # Display results
            self.display_results(output_path)
            
        except subprocess.CalledProcessError as e:
            click.echo(f"❌ Error running samtools: {e}")
            sys.exit(1)
        except Exception as e:
            click.echo(f"❌ Error: {e}")
            sys.exit(1)
    
    def display_results(self, output_file):
        """Display filtering results"""
        click.echo("\n" + "="*50)
        click.echo("✅ FILTERING COMPLETE")
        click.echo("="*50)
        
        click.echo(f"📊 Total reads processed:    {self.stats['total_reads']:,}")
        click.echo(f"✅ Reads passed filter:      {self.stats['passed_reads']:,} ({100*self.stats['passed_reads']/self.stats['total_reads']:.1f}%)")
        click.echo(f"❌ Reads filtered out:       {self.stats['filtered_reads']:,} ({100*self.stats['filtered_reads']/self.stats['total_reads']:.1f}%)")
        click.echo(f"\n🏷️  Tag Coverage (v2.1+):")
        click.echo(f"   XI (Structure):      {self.stats['reads_with_xi']:,}")
        click.echo(f"   XB (Boundary):       {self.stats['reads_with_xb']:,}")
        click.echo(f"   XS (Splicetag):      {self.stats['reads_with_xs']:,}")
        click.echo(f"   XT (Transcript):     {self.stats['reads_with_xt']:,}")
        click.echo(f"   XV (Variants):       {self.stats['reads_with_xv']:,}")
        click.echo(f"   XC (Cluster):        {self.stats['reads_with_xc']:,}")
        click.echo(f"   XR (Representative): {self.stats['reads_with_xr']:,}")
        
        click.echo(f"\n💾 Output file: {output_file}")
        
        if output_file.suffix.lower() == '.bam':
            click.echo(f"🚀 View results: samtools view {output_file} | head")
        else:
            click.echo(f"🚀 View results: head {output_file}")


@click.command()
@click.option('--input', '-i', 'input_file', required=True, type=click.Path(exists=True), help='Input BAM/SAM file')
@click.option('--output', '-o', required=True, help='Output BAM/SAM file')
# v1.0 tags
@click.option('--include-xi', type=click.Path(exists=True), help='File with XI (structure) IDs to include')
@click.option('--exclude-xi', type=click.Path(exists=True), help='File with XI (structure) IDs to exclude')
@click.option('--include-xv', type=click.Path(exists=True), help='File with XV (variant) IDs to include')
@click.option('--exclude-xv', type=click.Path(exists=True), help='File with XV (variant) IDs to exclude')
# v2.0+ tags
@click.option('--include-xb', type=click.Path(exists=True), help='File with XB (boundary) tags to include')
@click.option('--exclude-xb', type=click.Path(exists=True), help='File with XB (boundary) tags to exclude')
@click.option('--include-xs', type=click.Path(exists=True), help='File with XS (splicetag) tags to include')
@click.option('--exclude-xs', type=click.Path(exists=True), help='File with XS (splicetag) tags to exclude')
@click.option('--include-xt', type=click.Path(exists=True), help='File with XT (transcript group) IDs to include')
@click.option('--exclude-xt', type=click.Path(exists=True), help='File with XT (transcript group) IDs to exclude')
@click.option('--include-xc', type=click.Path(exists=True), help='File with XC (cluster) IDs to include')
@click.option('--exclude-xc', type=click.Path(exists=True), help='File with XC (cluster) IDs to exclude')
@click.option('--include-xr', type=click.Path(exists=True), help='File with XR (representative) tags to include')
@click.option('--exclude-xr', type=click.Path(exists=True), help='File with XR (representative) tags to exclude')
# Individual IDs
@click.option('--include-xi-id', multiple=True, help='Individual XI ID to include (can be used multiple times)')
@click.option('--exclude-xi-id', multiple=True, help='Individual XI ID to exclude (can be used multiple times)')
@click.option('--include-xv-id', multiple=True, help='Individual XV ID to include (can be used multiple times)')
@click.option('--exclude-xv-id', multiple=True, help='Individual XV ID to exclude (can be used multiple times)')
@click.option('--keep-untagged/--no-untagged', default=False, help='Keep reads without isotags')
def filter(input_file, output, include_xi, exclude_xi, include_xv, exclude_xv,
          include_xb, exclude_xb, include_xs, exclude_xs, include_xt, exclude_xt,
          include_xc, exclude_xc, include_xr, exclude_xr,
          include_xi_id, exclude_xi_id, include_xv_id, exclude_xv_id, keep_untagged):
    """
    ISO-Tools Filter - Filter BAM Files by Isotag IDs
    
    Filter reads based on isotag IDs. Can include/exclude specific structures
    or variants. Supports both file-based and command-line ID specification.
    
    Examples:
        # Include specific structure IDs
        iso-tools filter -i input.bam -o filtered.bam --include-xi structures.txt
        
        # Exclude specific variants
        iso-tools filter -i input.bam -o filtered.bam --exclude-xv variants.txt
        
        # Include specific IDs from command line
        iso-tools filter -i input.bam -o filtered.bam --include-xi-id ABC123 --include-xi-id DEF456
        
        # Complex filtering with multiple criteria
        iso-tools filter -i input.bam -o filtered.bam --include-xi structures.txt --exclude-xv variants.txt --keep-untagged
        
        # Keep only reads with specific structure, exclude untagged
        iso-tools filter -i input.bam -o filtered.bam --include-xi-id Uy3v73FzY4ZhB-w0ZLsDwQLJfMl34Hzx
    """
    
    filterer = IsotogFilter()

    # Load include/exclude lists from files (v1.0 tags)
    if include_xi:
        filterer.include_xi = filterer.load_isotag_list(include_xi)
        click.echo(f"📥 Loaded {len(filterer.include_xi)} XI IDs to include from {include_xi}")

    if exclude_xi:
        filterer.exclude_xi = filterer.load_isotag_list(exclude_xi)
        click.echo(f"📥 Loaded {len(filterer.exclude_xi)} XI IDs to exclude from {exclude_xi}")

    if include_xv:
        filterer.include_xv = filterer.load_isotag_list(include_xv)
        click.echo(f"📥 Loaded {len(filterer.include_xv)} XV IDs to include from {include_xv}")

    if exclude_xv:
        filterer.exclude_xv = filterer.load_isotag_list(exclude_xv)
        click.echo(f"📥 Loaded {len(filterer.exclude_xv)} XV IDs to exclude from {exclude_xv}")

    # Load include/exclude lists from files (v2.0+ tags)
    if include_xb:
        filterer.include_xb = filterer.load_isotag_list(include_xb)
        click.echo(f"📥 Loaded {len(filterer.include_xb)} XB tags to include from {include_xb}")

    if exclude_xb:
        filterer.exclude_xb = filterer.load_isotag_list(exclude_xb)
        click.echo(f"📥 Loaded {len(filterer.exclude_xb)} XB tags to exclude from {exclude_xb}")

    if include_xs:
        filterer.include_xs = filterer.load_isotag_list(include_xs)
        click.echo(f"📥 Loaded {len(filterer.include_xs)} XS tags to include from {include_xs}")

    if exclude_xs:
        filterer.exclude_xs = filterer.load_isotag_list(exclude_xs)
        click.echo(f"📥 Loaded {len(filterer.exclude_xs)} XS tags to exclude from {exclude_xs}")

    if include_xt:
        filterer.include_xt = filterer.load_isotag_list(include_xt)
        click.echo(f"📥 Loaded {len(filterer.include_xt)} XT IDs to include from {include_xt}")

    if exclude_xt:
        filterer.exclude_xt = filterer.load_isotag_list(exclude_xt)
        click.echo(f"📥 Loaded {len(filterer.exclude_xt)} XT IDs to exclude from {exclude_xt}")

    if include_xc:
        filterer.include_xc = filterer.load_isotag_list(include_xc)
        click.echo(f"📥 Loaded {len(filterer.include_xc)} XC IDs to include from {include_xc}")

    if exclude_xc:
        filterer.exclude_xc = filterer.load_isotag_list(exclude_xc)
        click.echo(f"📥 Loaded {len(filterer.exclude_xc)} XC IDs to exclude from {exclude_xc}")

    if include_xr:
        filterer.include_xr = filterer.load_isotag_list(include_xr)
        click.echo(f"📥 Loaded {len(filterer.include_xr)} XR tags to include from {include_xr}")

    if exclude_xr:
        filterer.exclude_xr = filterer.load_isotag_list(exclude_xr)
        click.echo(f"📥 Loaded {len(filterer.exclude_xr)} XR tags to exclude from {exclude_xr}")
    
    # Add individual IDs from command line
    if include_xi_id:
        filterer.include_xi.update(include_xi_id)
        click.echo(f"➕ Added {len(include_xi_id)} XI IDs to include from command line")
    
    if exclude_xi_id:
        filterer.exclude_xi.update(exclude_xi_id)
        click.echo(f"➖ Added {len(exclude_xi_id)} XI IDs to exclude from command line")
    
    if include_xv_id:
        filterer.include_xv.update(include_xv_id)
        click.echo(f"➕ Added {len(include_xv_id)} XV IDs to include from command line")
    
    if exclude_xv_id:
        filterer.exclude_xv.update(exclude_xv_id)
        click.echo(f"➖ Added {len(exclude_xv_id)} XV IDs to exclude from command line")
    
    # Validate that we have some filtering criteria
    total_filters = (len(filterer.include_xi) + len(filterer.exclude_xi) +
                    len(filterer.include_xb) + len(filterer.exclude_xb) +
                    len(filterer.include_xs) + len(filterer.exclude_xs) +
                    len(filterer.include_xt) + len(filterer.exclude_xt) +
                    len(filterer.include_xv) + len(filterer.exclude_xv) +
                    len(filterer.include_xc) + len(filterer.exclude_xc) +
                    len(filterer.include_xr) + len(filterer.exclude_xr))

    if total_filters == 0 and not keep_untagged:
        click.echo("❌ No filtering criteria specified. Use --include-*, --exclude-*, or --keep-untagged options.")
        sys.exit(1)
    
    # Run filtering
    filterer.filter_bam(input_file, output, keep_untagged)


if __name__ == '__main__':
    filter()