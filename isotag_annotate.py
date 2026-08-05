#!/usr/bin/env python3
"""
ISO-Tools Annotate - Match Isotags to Reference Annotations

Match isotag-annotated reads to reference transcriptome annotations.
Classify isoforms as exact matches, partial matches, or novel.
"""

import subprocess
import click
import sys
import re
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import gzip


class ReferenceAnnotation:
    """Represents a reference transcript annotation"""
    
    def __init__(self, transcript_id: str, gene_id: str, chrom: str, strand: str, exons: List[Tuple[int, int]]):
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.chrom = chrom
        self.strand = strand
        self.exons = sorted(exons)  # List of (start, end) tuples
        self.start = min(e[0] for e in exons)
        self.end = max(e[1] for e in exons)


class IsotogAnnotator:
    """Annotate isotags against reference transcriptome"""
    
    def __init__(self):
        self.reference_transcripts = {}  # transcript_id -> ReferenceAnnotation
        self.gene_transcripts = defaultdict(list)  # gene_id -> list of transcript_ids
        self.isotag_structures = {}  # isotag_id -> (chrom, strand, exons, read_count)
        self.annotation_results = {}  # isotag_id -> annotation_info
        self.stats = {
            'total_isotags': 0,
            'exact_matches': 0,
            'partial_matches': 0,
            'novel_isoforms': 0,
            'intergenic': 0,
            'antisense': 0,
            'reference_transcripts': 0,
            'reference_genes': 0
        }
    
    def load_gtf_annotation(self, gtf_file: str):
        """Load reference annotation from GTF file"""
        click.echo(f"📋 Loading reference annotation from {gtf_file}")
        
        transcripts = defaultdict(lambda: {'gene_id': None, 'chrom': None, 'strand': None, 'exons': []})
        
        # Handle gzipped files
        if gtf_file.endswith('.gz'):
            file_handle = gzip.open(gtf_file, 'rt')
        else:
            file_handle = open(gtf_file, 'r')
        
        try:
            for line in file_handle:
                # Skip comments and headers
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                chrom, source, feature, start, end, score, strand, frame, attributes = fields
                
                if feature not in ['transcript', 'exon']:
                    continue
                
                # Parse attributes
                attr_dict = {}
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr:
                        # Handle both quoted and unquoted attributes
                        match = re.match(r'(\w+)\s*["\s]*([^"]*?)["\s]*$', attr)
                        if match:
                            key, value = match.groups()
                            attr_dict[key] = value.strip('"').strip()
                
                transcript_id = attr_dict.get('transcript_id')
                gene_id = attr_dict.get('gene_id')
                
                if not transcript_id:
                    continue
                
                if feature == 'transcript':
                    transcripts[transcript_id]['gene_id'] = gene_id
                    transcripts[transcript_id]['chrom'] = chrom
                    transcripts[transcript_id]['strand'] = strand
                
                elif feature == 'exon':
                    exon_start = int(start)
                    exon_end = int(end)
                    transcripts[transcript_id]['exons'].append((exon_start, exon_end))
                    
                    # Set transcript info if not set by transcript feature
                    if not transcripts[transcript_id]['chrom']:
                        transcripts[transcript_id]['gene_id'] = gene_id
                        transcripts[transcript_id]['chrom'] = chrom
                        transcripts[transcript_id]['strand'] = strand
        
        finally:
            file_handle.close()
        
        # Convert to ReferenceAnnotation objects
        for transcript_id, info in transcripts.items():
            if info['exons'] and info['chrom'] and info['strand']:
                ref_annotation = ReferenceAnnotation(
                    transcript_id=transcript_id,
                    gene_id=info['gene_id'] or 'unknown',
                    chrom=info['chrom'],
                    strand=info['strand'],
                    exons=info['exons']
                )
                
                self.reference_transcripts[transcript_id] = ref_annotation
                self.gene_transcripts[ref_annotation.gene_id].append(transcript_id)
        
        self.stats['reference_transcripts'] = len(self.reference_transcripts)
        self.stats['reference_genes'] = len(self.gene_transcripts)
        
        click.echo(f"   ✅ Loaded {self.stats['reference_transcripts']} transcripts from {self.stats['reference_genes']} genes")
    
    def extract_isotag_structures(self, bam_file: str):
        """Extract isoform structures from isotag-annotated BAM file"""
        click.echo(f"🔍 Extracting isotag structures from {bam_file}")
        
        isotag_reads = defaultdict(list)
        
        view_cmd = ['samtools', 'view', bam_file]
        
        try:
            process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)
            
            reads_processed = 0
            for line in process.stdout:
                reads_processed += 1
                
                if reads_processed % 10000 == 0:
                    click.echo(f"   ⏳ Processed {reads_processed:,} reads...")
                
                fields = line.strip().split('\t')
                
                if len(fields) >= 11:
                    flag = int(fields[1])
                    rname = fields[2]
                    pos = int(fields[3])
                    cigar = fields[5]

                    # Skip unmapped reads
                    if flag & 0x4 or cigar == '*':
                        continue

                    # Extract all tags (v2.1+)
                    isotag_id = None
                    tags = {}
                    for field in fields[11:]:
                        if field.startswith('XI:Z:'):
                            isotag_id = field[5:]
                            tags['XI'] = isotag_id
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
                    
                    if isotag_id:
                        strand = '-' if flag & 0x10 else '+'
                        exons = self.parse_cigar_to_blocks(pos, cigar)

                        if exons:
                            isotag_reads[isotag_id].append({
                                'chrom': rname,
                                'strand': strand,
                                'exons': exons,
                                'tags': tags  # Store v2.1+ tags
                            })
            
            process.wait()
            
            # Generate consensus structures for each isotag
            for isotag_id, reads in isotag_reads.items():
                if reads:
                    # Use most common structure
                    structure_counter = Counter()
                    for read in reads:
                        structure_key = (read['chrom'], read['strand'], tuple(read['exons']))
                        structure_counter[structure_key] += 1

                    most_common = structure_counter.most_common(1)[0]
                    chrom, strand, exons = most_common[0]
                    read_count = most_common[1]

                    # Get additional tags from first read (v2.1+)
                    extra_tags = reads[0].get('tags', {})

                    self.isotag_structures[isotag_id] = {
                        'chrom': chrom,
                        'strand': strand,
                        'exons': list(exons),
                        'read_count': read_count,
                        'total_reads': len(reads),
                        'tags': extra_tags  # v2.1+ tags
                    }
            
            self.stats['total_isotags'] = len(self.isotag_structures)
            click.echo(f"   ✅ Extracted {self.stats['total_isotags']} unique isotag structures")
            
        except subprocess.CalledProcessError as e:
            click.echo(f"❌ Error reading BAM file: {e}")
            sys.exit(1)
    
    def parse_cigar_to_blocks(self, pos: int, cigar: str) -> List[Tuple[int, int]]:
        """Parse CIGAR to genomic blocks (exons)"""
        blocks = []
        current_pos = pos
        block_start = pos
        in_block = True
        
        pattern = r'(\d+)([MIDNSHP=X])'
        matches = re.findall(pattern, cigar)
        
        for length_str, op_char in matches:
            length = int(length_str)
            
            if op_char in ['M', '=', 'X']:  # Match/alignment
                if not in_block:
                    block_start = current_pos
                    in_block = True
                current_pos += length
            
            elif op_char == 'D':  # Deletion from reference
                current_pos += length
            
            elif op_char == 'I':  # Insertion to reference
                pass  # Don't advance reference position
            
            elif op_char == 'N':  # Intron/skipped region
                if in_block:
                    blocks.append((block_start, current_pos - 1))
                    in_block = False
                current_pos += length
            
            elif op_char in ['S', 'H']:  # Soft/hard clipping
                pass  # Don't affect block structure
        
        # Add final block
        if in_block and block_start < current_pos:
            blocks.append((block_start, current_pos - 1))
        
        return blocks
    
    def calculate_overlap(self, exons1: List[Tuple[int, int]], exons2: List[Tuple[int, int]]) -> float:
        """Calculate overlap between two exon sets"""
        if not exons1 or not exons2:
            return 0.0
        
        # Calculate total exonic length for each
        length1 = sum(end - start + 1 for start, end in exons1)
        length2 = sum(end - start + 1 for start, end in exons2)
        
        # Calculate overlapping bases
        overlap = 0
        for start1, end1 in exons1:
            for start2, end2 in exons2:
                overlap_start = max(start1, start2)
                overlap_end = min(end1, end2)
                if overlap_start <= overlap_end:
                    overlap += overlap_end - overlap_start + 1
        
        # Return Jaccard similarity
        union = length1 + length2 - overlap
        return overlap / union if union > 0 else 0.0
    
    def classify_isotag(self, isotag_id: str, structure: Dict) -> Dict:
        """Classify isotag against reference annotation"""
        chrom = structure['chrom']
        strand = structure['strand']
        exons = structure['exons']
        
        best_match = None
        best_overlap = 0.0
        best_classification = 'novel'
        
        # Find overlapping reference transcripts on same chromosome
        candidates = []
        for transcript_id, ref in self.reference_transcripts.items():
            if ref.chrom == chrom:
                # Calculate overlap
                overlap = self.calculate_overlap(exons, ref.exons)
                if overlap > 0.1:  # At least 10% overlap
                    candidates.append((transcript_id, ref, overlap))
        
        # Sort by overlap
        candidates.sort(key=lambda x: x[2], reverse=True)
        
        if candidates:
            best_transcript_id, best_ref, best_overlap = candidates[0]
            best_match = {
                'transcript_id': best_transcript_id,
                'gene_id': best_ref.gene_id,
                'overlap': best_overlap,
                'same_strand': strand == best_ref.strand
            }
            
            # Classify based on overlap and strand
            if best_overlap >= 0.95 and strand == best_ref.strand:
                if exons == best_ref.exons:
                    best_classification = 'exact_match'
                else:
                    best_classification = 'near_exact_match'
            elif best_overlap >= 0.5 and strand == best_ref.strand:
                best_classification = 'partial_match'
            elif best_overlap >= 0.3 and strand != best_ref.strand:
                best_classification = 'antisense'
            else:
                best_classification = 'novel'
        else:
            best_classification = 'intergenic'
        
        return {
            'classification': best_classification,
            'best_match': best_match,
            'all_candidates': [(t_id, overlap) for t_id, ref, overlap in candidates[:5]]
        }
    
    def annotate_isotags(self):
        """Annotate all isotags against reference"""
        click.echo(f"🔗 Annotating {len(self.isotag_structures)} isotags against reference...")
        
        for isotag_id, structure in self.isotag_structures.items():
            annotation = self.classify_isotag(isotag_id, structure)
            self.annotation_results[isotag_id] = annotation
            
            # Update stats
            classification = annotation['classification']
            if classification == 'exact_match':
                self.stats['exact_matches'] += 1
            elif classification in ['partial_match', 'near_exact_match']:
                self.stats['partial_matches'] += 1
            elif classification == 'antisense':
                self.stats['antisense'] += 1
            elif classification == 'intergenic':
                self.stats['intergenic'] += 1
            else:
                self.stats['novel_isoforms'] += 1
    
    def export_annotations(self, output_file: str):
        """Export annotation results to TSV"""
        click.echo(f"📄 Exporting annotations to {output_file}")
        
        with open(output_file, 'w') as f:
            # Write header (v2.1+ with additional tags)
            f.write("Isotag_ID\tClassification\tBest_Match_Transcript\tBest_Match_Gene\t"
                   "Overlap\tSame_Strand\tRead_Count\tTotal_Reads\tChrom\tStrand\tExon_Count\t"
                   "XB_Tag\tXS_Tag\tXT_Tag\tXC_Tag\n")
            
            for isotag_id, annotation in self.annotation_results.items():
                structure = self.isotag_structures[isotag_id]
                
                classification = annotation['classification']
                best_match = annotation['best_match']
                
                if best_match:
                    transcript_id = best_match['transcript_id']
                    gene_id = best_match['gene_id']
                    overlap = f"{best_match['overlap']:.3f}"
                    same_strand = str(best_match['same_strand'])
                else:
                    transcript_id = gene_id = overlap = same_strand = 'N/A'

                # Get v2.1+ tags
                tags = structure.get('tags', {})
                xb_tag = tags.get('XB', 'N/A')
                xs_tag = tags.get('XS', 'N/A')
                xt_tag = tags.get('XT', 'N/A')
                xc_tag = tags.get('XC', 'N/A')

                f.write(f"{isotag_id}\t{classification}\t{transcript_id}\t{gene_id}\t"
                       f"{overlap}\t{same_strand}\t{structure['read_count']}\t"
                       f"{structure['total_reads']}\t{structure['chrom']}\t{structure['strand']}\t"
                       f"{len(structure['exons'])}\t{xb_tag}\t{xs_tag}\t{xt_tag}\t{xc_tag}\n")
        
        click.echo(f"   ✅ Exported annotations for {len(self.annotation_results)} isotags")
    
    def display_summary(self):
        """Display annotation summary"""
        click.echo("\n" + "="*60)
        click.echo("🔗 ISO-TOOLS ANNOTATE SUMMARY")
        click.echo("="*60)
        
        click.echo(f"📚 Reference: {self.stats['reference_transcripts']:,} transcripts, {self.stats['reference_genes']:,} genes")
        click.echo(f"🆔 Isotags analyzed: {self.stats['total_isotags']:,}")
        
        click.echo(f"\n📊 CLASSIFICATION RESULTS")
        click.echo("-" * 40)
        
        total = self.stats['total_isotags']
        if total > 0:
            click.echo(f"✅ Exact matches:     {self.stats['exact_matches']:,} ({100*self.stats['exact_matches']/total:.1f}%)")
            click.echo(f"🔗 Partial matches:   {self.stats['partial_matches']:,} ({100*self.stats['partial_matches']/total:.1f}%)")
            click.echo(f"🆕 Novel isoforms:    {self.stats['novel_isoforms']:,} ({100*self.stats['novel_isoforms']/total:.1f}%)")
            click.echo(f"🔄 Antisense:         {self.stats['antisense']:,} ({100*self.stats['antisense']/total:.1f}%)")
            click.echo(f"🌌 Intergenic:        {self.stats['intergenic']:,} ({100*self.stats['intergenic']/total:.1f}%)")
            
            # Summary categories
            known = self.stats['exact_matches'] + self.stats['partial_matches']
            novel = total - known
            click.echo(f"\n📈 SUMMARY")
            click.echo(f"Known isoforms:  {known:,} ({100*known/total:.1f}%)")
            click.echo(f"Novel elements:  {novel:,} ({100*novel/total:.1f}%)")


@click.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.option('--gtf', '-g', required=True, type=click.Path(exists=True), help='Reference GTF annotation file')
@click.option('--output', '-o', required=True, help='Output TSV file with annotations')
def annotate(bam_file, gtf, output):
    """
    ISO-Tools Annotate - Match Isotags to Reference Annotations
    
    Match isotag structures against reference transcriptome annotations.
    Classify isoforms as exact matches, partial matches, or novel elements.
    
    Examples:
        # Annotate against GENCODE
        iso-tools annotate tagged.bam -g gencode.v44.annotation.gtf -o annotations.tsv
        
        # Use gzipped GTF file
        iso-tools annotate tagged.bam -g ensembl.gtf.gz -o annotations.tsv
    """
    
    annotator = IsotogAnnotator()
    
    # Load reference annotation
    annotator.load_gtf_annotation(gtf)
    
    if annotator.stats['reference_transcripts'] == 0:
        click.echo("❌ No transcripts loaded from GTF file")
        sys.exit(1)
    
    # Extract isotag structures
    annotator.extract_isotag_structures(bam_file)
    
    if annotator.stats['total_isotags'] == 0:
        click.echo("❌ No isotag-annotated reads found in BAM file")
        sys.exit(1)
    
    # Perform annotation
    annotator.annotate_isotags()
    
    # Export results
    annotator.export_annotations(output)
    
    # Display summary
    annotator.display_summary()
    
    click.echo(f"\n✅ Annotation complete!")
    click.echo(f"📄 Results: {output}")
    click.echo(f"💡 Next steps: Analyze novel isoforms, validate exact matches")


if __name__ == '__main__':
    annotate()