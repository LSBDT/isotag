#!/usr/bin/env python3
"""
ISO-Tools Convert - Convert Isotag Data to Standard Formats (v2.1+)

Convert isotag-annotated reads to standard genomic formats like BED, GTF, GFF3.
Supports all 7 tags: XI, XB, XS, XT, XV, XC, XR
Enables integration with genome browsers and annotation tools.
"""

import subprocess
import click
import sys
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import re
import json


class IsotagConverter:
    """Convert isotag data to standard genomic formats (v2.1+ all 7 tags)"""

    def __init__(self, use_tag: str = 'XI'):
        self.use_tag = use_tag  # Which tag to use for grouping (XI, XB, XS, XT, XC)
        self.isoform_structures = defaultdict(list)  # tag_id -> list of structures
        self.isotag_to_reads = defaultdict(list)     # tag_id -> list of read_info

        # Tag-specific storage
        self.xi_structures = defaultdict(list)
        self.xb_structures = defaultdict(list)
        self.xs_structures = defaultdict(list)
        self.xt_structures = defaultdict(list)
        self.xc_structures = defaultdict(list)

        # Junction tracking (for junction matrix)
        self.junction_counts = Counter()  # (chrom, start, end, strand) -> count
        self.junction_to_tags = defaultdict(set)  # junction -> set of tags

        self.stats = {
            'total_reads': 0,
            'tagged_reads': 0,
            'unique_structures': 0,
            'reads_with_xi': 0,
            'reads_with_xb': 0,
            'reads_with_xs': 0,
            'reads_with_xt': 0,
            'reads_with_xv': 0,
            'reads_with_xc': 0,
            'reads_with_xr': 0
        }
    
    def parse_cigar_operations(self, cigar_string: str) -> List[Tuple[str, int]]:
        """Parse CIGAR string into operations"""
        if not cigar_string or cigar_string == "*":
            return []
        
        operations = []
        pattern = r'(\d+)([MIDNSHP=X])'
        matches = re.findall(pattern, cigar_string)
        
        for length_str, op_char in matches:
            operations.append((op_char, int(length_str)))
        
        return operations
    
    def cigar_to_blocks(self, pos: int, cigar_ops: List[Tuple[str, int]]) -> List[Tuple[int, int]]:
        """Convert CIGAR operations to genomic blocks (exons)"""
        blocks = []
        current_pos = pos
        block_start = pos
        in_block = True
        
        for op, length in cigar_ops:
            if op in ['M', '=', 'X']:  # Match/alignment
                if not in_block:
                    block_start = current_pos
                    in_block = True
                current_pos += length
            
            elif op == 'D':  # Deletion from reference
                current_pos += length
            
            elif op == 'I':  # Insertion to reference
                pass  # Don't advance reference position
            
            elif op == 'N':  # Intron/skipped region
                if in_block:
                    blocks.append((block_start, current_pos - 1))
                    in_block = False
                current_pos += length
            
            elif op in ['S', 'H']:  # Soft/hard clipping
                pass  # Don't affect block structure
        
        # Add final block
        if in_block and block_start < current_pos:
            blocks.append((block_start, current_pos - 1))
        
        return blocks
    
    def extract_isoform_structures(self, bam_file: str):
        """Extract isoform structures from BAM file"""
        click.echo(f"🔍 Extracting isoform structures from {bam_file}...")
        
        view_cmd = ['samtools', 'view', bam_file]
        
        try:
            process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)
            
            for line in process.stdout:
                self.stats['total_reads'] += 1
                
                if self.stats['total_reads'] % 10000 == 0:
                    click.echo(f"   ⏳ Processed {self.stats['total_reads']:,} reads...")
                
                fields = line.strip().split('\t')
                
                if len(fields) >= 11:
                    qname = fields[0]
                    flag = int(fields[1])
                    rname = fields[2]
                    pos = int(fields[3])
                    cigar = fields[5]
                    
                    # Skip unmapped reads
                    if flag & 0x4 or cigar == '*':
                        continue
                    
                    # Determine strand
                    strand = '-' if flag & 0x10 else '+'

                    # Extract all 7 tags
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

                    # Get the tag we're using for grouping
                    group_tag = tags.get(self.use_tag)

                    if group_tag or tags.get('XI'):  # At minimum need XI or specified tag
                        self.stats['tagged_reads'] += 1

                        # Parse CIGAR to get exon blocks
                        cigar_ops = self.parse_cigar_operations(cigar)
                        if cigar_ops:
                            blocks = self.cigar_to_blocks(pos, cigar_ops)

                            if blocks:
                                # Store read information with all tags
                                read_info = {
                                    'qname': qname,
                                    'chrom': rname,
                                    'strand': strand,
                                    'blocks': blocks,
                                    'start': min(b[0] for b in blocks),
                                    'end': max(b[1] for b in blocks),
                                    'tags': tags  # Store all tags
                                }

                                # Add to appropriate storage
                                if group_tag:
                                    self.isotag_to_reads[group_tag].append(read_info)

                                # Also add to tag-specific storage
                                if 'XI' in tags:
                                    self.xi_structures[tags['XI']].append(read_info)
                                if 'XB' in tags:
                                    self.xb_structures[tags['XB']].append(read_info)
                                if 'XS' in tags:
                                    self.xs_structures[tags['XS']].append(read_info)
                                if 'XT' in tags:
                                    self.xt_structures[tags['XT']].append(read_info)
                                if 'XC' in tags:
                                    self.xc_structures[tags['XC']].append(read_info)

                                # Track junctions for junction matrix
                                if len(blocks) > 1:
                                    for i in range(len(blocks) - 1):
                                        junction = (rname, blocks[i][1], blocks[i+1][0], strand)
                                        self.junction_counts[junction] += 1
                                        if group_tag:
                                            self.junction_to_tags[junction].add(group_tag)
            
            process.wait()
            
            # Compute consensus structures for each isotag
            self._compute_consensus_structures()
            
            self.stats['unique_structures'] = len(self.isoform_structures)
            
            click.echo(f"   ✅ Extracted {self.stats['tagged_reads']:,} tagged reads")
            click.echo(f"   🆔 Found {self.stats['unique_structures']:,} unique isoform structures")
            
        except subprocess.CalledProcessError as e:
            click.echo(f"❌ Error reading BAM file: {e}")
            sys.exit(1)
    
    def _compute_consensus_structures(self):
        """Compute consensus structure for each isotag"""
        for isotag_id, reads in self.isotag_to_reads.items():
            if not reads:
                continue
            
            # Group reads by chromosome and strand
            chrom_strand_groups = defaultdict(list)
            for read in reads:
                key = (read['chrom'], read['strand'])
                chrom_strand_groups[key].append(read)
            
            # For each chromosome/strand combination, find consensus
            for (chrom, strand), group_reads in chrom_strand_groups.items():
                # Simple consensus: use the most common block structure
                block_patterns = defaultdict(int)
                
                for read in group_reads:
                    pattern = tuple(read['blocks'])
                    block_patterns[pattern] += 1
                
                # Get most common pattern
                if block_patterns:
                    consensus_blocks = max(block_patterns.keys(), key=lambda x: block_patterns[x])
                    read_count = block_patterns[consensus_blocks]
                    
                    start = min(b[0] for b in consensus_blocks)
                    end = max(b[1] for b in consensus_blocks)
                    
                    self.isoform_structures[isotag_id].append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'blocks': list(consensus_blocks),
                        'read_count': read_count,
                        'total_reads': len(group_reads)
                    })
    
    def export_bed12(self, output_file: str, include_variants: bool = False):
        """Export to BED12 format"""
        click.echo(f"📄 Exporting to BED12: {output_file}")
        
        with open(output_file, 'w') as f:
            for isotag_id, structures in self.isoform_structures.items():
                for i, struct in enumerate(structures):
                    chrom = struct['chrom']
                    start = struct['start'] - 1  # BED is 0-based
                    end = struct['end']
                    strand = struct['strand']
                    blocks = struct['blocks']
                    read_count = struct['read_count']
                    
                    # Create name
                    if len(structures) > 1:
                        name = f"{isotag_id}_{i+1}"
                    else:
                        name = isotag_id
                    
                    # Calculate block information
                    block_count = len(blocks)
                    block_sizes = []
                    block_starts = []
                    
                    for block_start, block_end in blocks:
                        block_sizes.append(block_end - block_start + 1)
                        block_starts.append(block_start - struct['start'])
                    
                    block_sizes_str = ','.join(map(str, block_sizes)) + ','
                    block_starts_str = ','.join(map(str, block_starts)) + ','
                    
                    # BED12 format
                    f.write(f"{chrom}\t{start}\t{end}\t{name}\t{read_count}\t{strand}\t"
                           f"{start}\t{end}\t0\t{block_count}\t{block_sizes_str}\t{block_starts_str}\n")
        
        click.echo(f"   ✅ Exported {sum(len(s) for s in self.isoform_structures.values())} isoform entries")
    
    def export_gtf(self, output_file: str, source: str = "isotag"):
        """Export to GTF format"""
        click.echo(f"📄 Exporting to GTF: {output_file}")
        
        with open(output_file, 'w') as f:
            # Write GTF header
            f.write(f'#gtf-version 2.2\n')
            f.write(f'#source: iso-tools convert\n')
            
            for isotag_id, structures in self.isoform_structures.items():
                for i, struct in enumerate(structures):
                    chrom = struct['chrom']
                    strand = struct['strand']
                    blocks = struct['blocks']
                    read_count = struct['read_count']
                    
                    # Create gene and transcript IDs
                    if len(structures) > 1:
                        transcript_id = f"{isotag_id}_{i+1}"
                    else:
                        transcript_id = isotag_id
                    
                    gene_id = f"GENE_{isotag_id.split('_')[0] if '_' in isotag_id else isotag_id[:8]}"
                    
                    # Write transcript feature
                    start = struct['start']
                    end = struct['end']
                    attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; read_count "{read_count}";'
                    
                    f.write(f"{chrom}\t{source}\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n")
                    
                    # Write exon features
                    for exon_num, (exon_start, exon_end) in enumerate(blocks, 1):
                        exon_attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; exon_number "{exon_num}";'
                        f.write(f"{chrom}\t{source}\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t{exon_attributes}\n")
        
        total_transcripts = sum(len(s) for s in self.isoform_structures.values())
        total_exons = sum(len(struct['blocks']) for structures in self.isoform_structures.values() 
                         for struct in structures)
        
        click.echo(f"   ✅ Exported {total_transcripts} transcripts with {total_exons} exons")
    
    def export_gff3(self, output_file: str, source: str = "isotag"):
        """Export to GFF3 format"""
        click.echo(f"📄 Exporting to GFF3: {output_file}")
        
        with open(output_file, 'w') as f:
            # Write GFF3 header
            f.write('##gff-version 3\n')
            f.write(f'##source: iso-tools convert\n')
            
            gene_counter = 1
            
            for isotag_id, structures in self.isoform_structures.items():
                gene_id = f"gene_{gene_counter:06d}"
                gene_counter += 1
                
                # Calculate gene boundaries
                all_starts = [struct['start'] for struct in structures]
                all_ends = [struct['end'] for struct in structures]
                gene_start = min(all_starts)
                gene_end = max(all_ends)
                gene_chrom = structures[0]['chrom']  # Assume same chromosome
                gene_strand = structures[0]['strand']  # Assume same strand
                
                # Write gene feature
                gene_attributes = f"ID={gene_id};Name={isotag_id}_gene"
                f.write(f"{gene_chrom}\t{source}\tgene\t{gene_start}\t{gene_end}\t.\t{gene_strand}\t.\t{gene_attributes}\n")
                
                for i, struct in enumerate(structures):
                    chrom = struct['chrom']
                    strand = struct['strand'] 
                    blocks = struct['blocks']
                    read_count = struct['read_count']
                    
                    # Create mRNA ID
                    if len(structures) > 1:
                        mrna_id = f"{gene_id}.{i+1}"
                    else:
                        mrna_id = f"{gene_id}.1"
                    
                    # Write mRNA feature
                    start = struct['start']
                    end = struct['end']
                    mrna_attributes = f"ID={mrna_id};Parent={gene_id};Name={isotag_id};read_count={read_count}"
                    
                    f.write(f"{chrom}\t{source}\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t{mrna_attributes}\n")
                    
                    # Write exon features
                    for exon_num, (exon_start, exon_end) in enumerate(blocks, 1):
                        exon_id = f"{mrna_id}.exon{exon_num}"
                        exon_attributes = f"ID={exon_id};Parent={mrna_id}"
                        f.write(f"{chrom}\t{source}\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t{exon_attributes}\n")
        
        total_genes = len(self.isoform_structures)
        total_transcripts = sum(len(s) for s in self.isoform_structures.values())
        total_exons = sum(len(struct['blocks']) for structures in self.isoform_structures.values() 
                         for struct in structures)
        
        click.echo(f"   ✅ Exported {total_genes} genes, {total_transcripts} transcripts, {total_exons} exons")
    
    def export_junction_matrix(self, output_file: str):
        """Export junction usage matrix"""
        click.echo(f"📄 Exporting junction matrix: {output_file}")

        with open(output_file, 'w') as f:
            # Write header
            f.write("Chromosome\tStart\tEnd\tStrand\tCount\tAssociated_Tags\n")

            # Sort junctions by count (descending)
            sorted_junctions = sorted(self.junction_counts.items(),
                                    key=lambda x: x[1], reverse=True)

            for junction, count in sorted_junctions:
                chrom, start, end, strand = junction
                associated_tags = self.junction_to_tags[junction]
                tags_str = ','.join(list(associated_tags)[:10])  # Limit to 10 tags

                if len(associated_tags) > 10:
                    tags_str += f" (+{len(associated_tags) - 10} more)"

                f.write(f"{chrom}\t{start}\t{end}\t{strand}\t{count}\t{tags_str}\n")

        click.echo(f"   ✅ Exported {len(self.junction_counts)} junctions")

    def export_tag_summary_json(self, output_file: str):
        """Export tag summary as JSON (all 7 tags)"""
        click.echo(f"📄 Exporting tag summary: {output_file}")

        summary = {
            'statistics': self.stats,
            'unique_counts': {
                'XI': len(self.xi_structures),
                'XB': len(self.xb_structures),
                'XS': len(self.xs_structures),
                'XT': len(self.xt_structures),
                'XC': len(self.xc_structures)
            },
            'junctions': {
                'total_junctions': len(self.junction_counts),
                'top_junctions': []
            }
        }

        # Add top 20 junctions
        sorted_junctions = sorted(self.junction_counts.items(),
                                key=lambda x: x[1], reverse=True)[:20]
        for junction, count in sorted_junctions:
            chrom, start, end, strand = junction
            summary['junctions']['top_junctions'].append({
                'chromosome': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'count': count
            })

        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2)

        click.echo(f"   ✅ Tag summary exported")

    def display_summary(self):
        """Display conversion summary (v2.1+ all 7 tags)"""
        click.echo("\n" + "="*60)
        click.echo(f"🔄 ISO-TOOLS CONVERT SUMMARY (v2.1+)")
        click.echo("="*60)

        click.echo(f"📊 Total reads processed:     {self.stats['total_reads']:,}")
        click.echo(f"🏷️  Tagged reads:              {self.stats['tagged_reads']:,}")
        click.echo(f"🆔 Unique isoform structures: {self.stats['unique_structures']:,}")

        click.echo(f"\n📊 TAG PRESENCE")
        click.echo(f"-" * 40)
        click.echo(f"Reads with XI: {self.stats['reads_with_xi']:,} ({len(self.xi_structures)} unique)")
        click.echo(f"Reads with XB: {self.stats['reads_with_xb']:,} ({len(self.xb_structures)} unique)")
        click.echo(f"Reads with XS: {self.stats['reads_with_xs']:,} ({len(self.xs_structures)} unique)")
        click.echo(f"Reads with XT: {self.stats['reads_with_xt']:,} ({len(self.xt_structures)} unique)")
        click.echo(f"Reads with XV: {self.stats['reads_with_xv']:,}")
        click.echo(f"Reads with XC: {self.stats['reads_with_xc']:,} ({len(self.xc_structures)} unique)")
        click.echo(f"Reads with XR: {self.stats['reads_with_xr']:,}")

        # Junction statistics
        if len(self.junction_counts) > 0:
            click.echo(f"\n🔗 JUNCTION STATISTICS")
            click.echo(f"-" * 40)
            click.echo(f"Total junctions: {len(self.junction_counts):,}")
            top_junction_count = max(self.junction_counts.values()) if self.junction_counts else 0
            click.echo(f"Top junction count: {top_junction_count:,}")

        # Structure complexity
        total_isoforms = sum(len(structures) for structures in self.isoform_structures.values())
        if total_isoforms > 0:
            avg_isoforms_per_structure = total_isoforms / len(self.isoform_structures)
            click.echo(f"\n📈 STRUCTURE COMPLEXITY")
            click.echo(f"-" * 40)
            click.echo(f"Total isoform entries:     {total_isoforms:,}")
            click.echo(f"Avg entries per structure: {avg_isoforms_per_structure:.1f}")


@click.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.option('--output', '-o', required=True, help='Output file')
@click.option('--format', '-f', type=click.Choice(['bed', 'bed12', 'gtf', 'gff3', 'junction', 'summary', 'auto']),
              default='auto', help='Output format (auto-detected from extension)')
@click.option('--source', '-s', default='isotag', help='Source field for GTF/GFF3')
@click.option('--use-tag', '-t', type=click.Choice(['XI', 'XB', 'XS', 'XT', 'XC']),
              default='XI', help='Tag to use for grouping (default: XI)')
@click.option('--include-variants', is_flag=True, help='Include variant information in output')
def convert(bam_file, output, format, source, use_tag, include_variants):
    """
    ISO-Tools Convert - Convert to Standard Genomic Formats (v2.1+)

    Convert isotag-annotated reads to BED, GTF, GFF3, or custom formats.
    Supports all 7 tags: XI, XB, XS, XT, XV, XC, XR

    Output Formats:
        bed/bed12:  BED12 format for genome browsers
        gtf:        GTF format for annotation
        gff3:       GFF3 format for annotation
        junction:   Junction usage matrix (TSV)
        summary:    Tag summary statistics (JSON)
        auto:       Auto-detect from file extension

    Grouping Tags:
        XI: Group by structure (full exon coordinates)
        XB: Group by boundary (5'/3' ends)
        XS: Group by splicetag (splice junctions)
        XT: Group by transcript (biological groups)
        XC: Group by cluster (clustered isoforms)

    Examples:
        # Convert to BED12 format (default: XI grouping)
        isotag_convert.py tagged.bam -o isoforms.bed

        # Convert using XC (cluster) grouping
        isotag_convert.py clustered.bam -o clusters.bed -t XC

        # Convert using XS (splicetag) grouping
        isotag_convert.py tagged.bam -o junctions.bed -t XS

        # Export junction matrix
        isotag_convert.py tagged.bam -o junctions.tsv -f junction

        # Export tag summary
        isotag_convert.py tagged.bam -o summary.json -f summary

        # Convert to GTF (auto-detected from extension)
        isotag_convert.py tagged.bam -o isoforms.gtf

        # Convert to GFF3 with custom source
        isotag_convert.py tagged.bam -o isoforms.gff3 -s "MyStudy"

        # Include variant information
        isotag_convert.py tagged.bam -o isoforms.bed --include-variants
    """

    converter = IsotagConverter(use_tag=use_tag)

    # Extract isoform structures
    converter.extract_isoform_structures(bam_file)

    if converter.stats['unique_structures'] == 0 and format not in ['junction', 'summary']:
        click.echo("❌ No isotag-annotated reads found in BAM file")
        sys.exit(1)

    # Auto-detect format from extension
    output_path = Path(output)
    if format == 'auto':
        ext = output_path.suffix.lower()
        if ext == '.bed':
            format = 'bed12'
        elif ext == '.gtf':
            format = 'gtf'
        elif ext in ['.gff', '.gff3']:
            format = 'gff3'
        elif ext == '.tsv':
            format = 'junction'
        elif ext == '.json':
            format = 'summary'
        else:
            click.echo("❌ Cannot auto-detect format. Please specify --format")
            sys.exit(1)

    # Convert to specified format
    if format in ['bed', 'bed12']:
        converter.export_bed12(output, include_variants)
    elif format == 'gtf':
        converter.export_gtf(output, source)
    elif format == 'gff3':
        converter.export_gff3(output, source)
    elif format == 'junction':
        converter.export_junction_matrix(output)
    elif format == 'summary':
        converter.export_tag_summary_json(output)

    # Display summary
    converter.display_summary()

    click.echo(f"\n✅ Conversion complete!")
    click.echo(f"📄 Output file: {output}")

    # Suggest viewing commands
    if format in ['bed', 'bed12']:
        click.echo(f"🚀 View in UCSC: Upload to genome browser as BED track")
        click.echo(f"🚀 Command line: head -5 {output}")
    elif format in ['gtf', 'gff3']:
        click.echo(f"🚀 View with IGV: Load as annotation track")
        click.echo(f"🚀 Command line: head -10 {output}")
    elif format == 'junction':
        click.echo(f"🚀 Command line: head -20 {output}")
    elif format == 'summary':
        click.echo(f"🚀 Command line: cat {output} | jq .")


if __name__ == '__main__':
    convert()