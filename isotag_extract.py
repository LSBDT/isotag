#!/usr/bin/env python3
"""
ISO-Tools Extract - Extract Sequences for Specific Isotags

Extract FASTA sequences for specific isotag IDs from BAM files.
Generate consensus sequences and individual read sequences.
"""

import subprocess
import click
import sys
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import re


class IsotogExtractor:
    """Extract sequences for specific isotags"""
    
    def __init__(self):
        self.isotag_reads = defaultdict(list)  # isotag_id -> list of read_info
        self.isotag_sequences = defaultdict(list)  # isotag_id -> list of sequences
        self.stats = {
            'total_reads': 0,
            'tagged_reads': 0,
            'target_reads': 0,
            'extracted_sequences': 0
        }
        self.target_isotags = set()
        
    def load_target_isotags(self, isotag_file: Optional[str], isotag_ids: List[str]):
        """Load target isotag IDs from file or command line"""
        if isotag_file:
            with open(isotag_file, 'r') as f:
                for line in f:
                    isotag_id = line.strip()
                    if isotag_id and not isotag_id.startswith('#'):
                        self.target_isotags.add(isotag_id)
        
        if isotag_ids:
            self.target_isotags.update(isotag_ids)
        
        if not self.target_isotags:
            click.echo("❌ No target isotag IDs specified")
            sys.exit(1)
        
        click.echo(f"🎯 Targeting {len(self.target_isotags)} isotag IDs")
    
    def reverse_complement(self, seq: str) -> str:
        """Generate reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(seq.upper()))
    
    def parse_cigar_for_sequence(self, seq: str, cigar: str) -> str:
        """Extract the aligned portion of sequence based on CIGAR"""
        if not cigar or cigar == '*':
            return seq
        
        # Parse CIGAR operations
        operations = []
        pattern = r'(\d+)([MIDNSHP=X])'
        matches = re.findall(pattern, cigar)
        
        result_seq = ""
        seq_pos = 0
        
        for length_str, op_char in matches:
            length = int(length_str)
            
            if op_char in ['M', '=', 'X']:  # Match/mismatch - keep sequence
                result_seq += seq[seq_pos:seq_pos + length]
                seq_pos += length
            elif op_char == 'I':  # Insertion - keep sequence
                result_seq += seq[seq_pos:seq_pos + length] 
                seq_pos += length
            elif op_char == 'D':  # Deletion - skip (nothing from read sequence)
                pass
            elif op_char in ['S', 'H']:  # Clipping - skip clipped bases
                seq_pos += length
            elif op_char == 'N':  # Intron - skip (nothing from read sequence)
                pass
        
        return result_seq
    
    def extract_isotag_reads(self, bam_file: str):
        """Extract reads matching target isotags"""
        click.echo(f"🔍 Extracting reads for target isotags from {bam_file}...")
        
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
                    seq = fields[9]
                    qual = fields[10]
                    
                    # Skip unmapped reads or reads without sequence
                    if flag & 0x4 or seq == '*':
                        continue
                    
                    # Extract isotag
                    isotag_id = None
                    variant_id = None
                    
                    for field in fields[11:]:
                        if field.startswith('XI:Z:'):
                            isotag_id = field[5:]
                            self.stats['tagged_reads'] += 1
                        elif field.startswith('XV:Z:'):
                            variant_id = field[5:]
                    
                    # Check if this is a target isotag
                    if isotag_id and isotag_id in self.target_isotags:
                        self.stats['target_reads'] += 1
                        
                        # Determine strand and sequence orientation
                        strand = '-' if flag & 0x10 else '+'
                        
                        # Get the aligned sequence portion
                        aligned_seq = self.parse_cigar_for_sequence(seq, cigar)
                        
                        # Reverse complement if on negative strand
                        if strand == '-':
                            final_seq = self.reverse_complement(aligned_seq)
                        else:
                            final_seq = aligned_seq
                        
                        read_info = {
                            'qname': qname,
                            'chrom': rname,
                            'pos': pos,
                            'strand': strand,
                            'cigar': cigar,
                            'sequence': final_seq,
                            'quality': qual,
                            'variant_id': variant_id,
                            'raw_sequence': seq,
                            'flag': flag
                        }
                        
                        self.isotag_reads[isotag_id].append(read_info)
            
            process.wait()
            
            click.echo(f"   ✅ Found {self.stats['target_reads']:,} reads matching target isotags")
            
            # Report per-isotag counts
            for isotag_id in sorted(self.target_isotags):
                count = len(self.isotag_reads[isotag_id])
                if count > 0:
                    click.echo(f"      {isotag_id}: {count} reads")
                else:
                    click.echo(f"      {isotag_id}: 0 reads ⚠️")
            
        except subprocess.CalledProcessError as e:
            click.echo(f"❌ Error reading BAM file: {e}")
            sys.exit(1)
    
    def generate_consensus_sequence(self, sequences: List[str]) -> str:
        """Generate consensus sequence from multiple reads"""
        if not sequences:
            return ""
        
        if len(sequences) == 1:
            return sequences[0]
        
        # Find the longest sequence as template
        max_len = max(len(seq) for seq in sequences)
        
        consensus = []
        for pos in range(max_len):
            base_counts = Counter()
            
            for seq in sequences:
                if pos < len(seq):
                    base_counts[seq[pos].upper()] += 1
            
            # Choose most common base, or N if tie
            if base_counts:
                most_common = base_counts.most_common(1)[0]
                consensus.append(most_common[0])
            else:
                consensus.append('N')
        
        return ''.join(consensus)
    
    def export_individual_sequences(self, output_file: str):
        """Export individual read sequences"""
        click.echo(f"📄 Exporting individual sequences to {output_file}")
        
        sequence_count = 0
        
        with open(output_file, 'w') as f:
            for isotag_id in sorted(self.target_isotags):
                reads = self.isotag_reads[isotag_id]
                
                for i, read in enumerate(reads):
                    header = f">{isotag_id}_{read['qname']}"
                    if read['variant_id']:
                        header += f"_var:{read['variant_id'][:16]}..."
                    
                    header += f" {read['chrom']}:{read['pos']} {read['strand']}"
                    
                    f.write(f"{header}\n")
                    f.write(f"{read['sequence']}\n")
                    sequence_count += 1
        
        self.stats['extracted_sequences'] = sequence_count
        click.echo(f"   ✅ Exported {sequence_count} individual sequences")
    
    def export_consensus_sequences(self, output_file: str):
        """Export consensus sequences for each isotag"""
        click.echo(f"📄 Exporting consensus sequences to {output_file}")
        
        consensus_count = 0
        
        with open(output_file, 'w') as f:
            for isotag_id in sorted(self.target_isotags):
                reads = self.isotag_reads[isotag_id]
                
                if reads:
                    # Generate consensus from all sequences
                    sequences = [read['sequence'] for read in reads]
                    consensus_seq = self.generate_consensus_sequence(sequences)
                    
                    # Create header with metadata
                    header = f">{isotag_id}_consensus"
                    header += f" reads:{len(reads)}"
                    
                    # Add location info from first read
                    first_read = reads[0]
                    header += f" {first_read['chrom']}:{first_read['pos']}"
                    header += f" {first_read['strand']}"
                    
                    # Add variant info if consistent
                    variant_ids = set(r['variant_id'] for r in reads if r['variant_id'])
                    if len(variant_ids) == 1:
                        header += f" var:{list(variant_ids)[0][:16]}..."
                    elif len(variant_ids) > 1:
                        header += f" variants:{len(variant_ids)}"
                    
                    f.write(f"{header}\n")
                    f.write(f"{consensus_seq}\n")
                    consensus_count += 1
        
        click.echo(f"   ✅ Exported {consensus_count} consensus sequences")
    
    def export_read_info(self, output_file: str):
        """Export detailed read information as TSV"""
        click.echo(f"📄 Exporting read details to {output_file}")
        
        with open(output_file, 'w') as f:
            # Write header
            f.write("Isotag_ID\tRead_Name\tChrom\tPos\tStrand\tCIGAR\t"
                   "Sequence_Length\tVariant_ID\tFlag\n")
            
            total_rows = 0
            for isotag_id in sorted(self.target_isotags):
                reads = self.isotag_reads[isotag_id]
                
                for read in reads:
                    f.write(f"{isotag_id}\t{read['qname']}\t{read['chrom']}\t"
                           f"{read['pos']}\t{read['strand']}\t{read['cigar']}\t"
                           f"{len(read['sequence'])}\t{read['variant_id'] or 'None'}\t"
                           f"{read['flag']}\n")
                    total_rows += 1
        
        click.echo(f"   ✅ Exported details for {total_rows} reads")
    
    def display_summary(self):
        """Display extraction summary"""
        click.echo("\n" + "="*60)
        click.echo("🧬 ISO-TOOLS EXTRACT SUMMARY")
        click.echo("="*60)
        
        click.echo(f"📊 Total reads processed:     {self.stats['total_reads']:,}")
        click.echo(f"🏷️  Tagged reads:              {self.stats['tagged_reads']:,}")
        click.echo(f"🎯 Target isotag reads:       {self.stats['target_reads']:,}")
        click.echo(f"🧬 Sequences extracted:       {self.stats['extracted_sequences']:,}")
        
        click.echo(f"\n🎯 TARGET ISOTAGS")
        click.echo("-" * 40)
        
        found_isotags = 0
        for isotag_id in sorted(self.target_isotags):
            read_count = len(self.isotag_reads[isotag_id])
            status = "✅" if read_count > 0 else "❌"
            click.echo(f"{status} {isotag_id}: {read_count} reads")
            
            if read_count > 0:
                found_isotags += 1
        
        click.echo(f"\n📈 Found sequences for {found_isotags}/{len(self.target_isotags)} target isotags")


@click.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.option('--output', '-o', required=True, help='Output FASTA file prefix')
@click.option('--isotag-file', type=click.Path(exists=True), help='File with isotag IDs to extract (one per line)')
@click.option('--isotag-id', multiple=True, help='Individual isotag ID to extract (can be used multiple times)')
@click.option('--mode', '-m', type=click.Choice(['individual', 'consensus', 'both']), 
              default='both', help='Extraction mode')
@click.option('--details', is_flag=True, help='Export detailed read information as TSV')
def extract(bam_file, output, isotag_file, isotag_id, mode, details):
    """
    ISO-Tools Extract - Extract Sequences for Specific Isotags
    
    Extract FASTA sequences for specific isotag IDs from BAM files.
    Can generate individual read sequences or consensus sequences.
    
    Examples:
        # Extract sequences for specific isotags
        iso-tools extract tagged.bam -o sequences --isotag-id ABC123 --isotag-id DEF456
        
        # Extract from file list
        iso-tools extract tagged.bam -o sequences --isotag-file isotags.txt
        
        # Generate only consensus sequences
        iso-tools extract tagged.bam -o consensus --isotag-file isotags.txt -m consensus
        
        # Include detailed read information
        iso-tools extract tagged.bam -o results --isotag-file isotags.txt --details
    """
    
    extractor = IsotogExtractor()
    
    # Load target isotags
    extractor.load_target_isotags(isotag_file, list(isotag_id))
    
    # Extract reads
    extractor.extract_isotag_reads(bam_file)
    
    if extractor.stats['target_reads'] == 0:
        click.echo("❌ No reads found matching target isotags")
        sys.exit(1)
    
    # Export sequences
    output_path = Path(output)
    
    if mode in ['individual', 'both']:
        individual_file = f"{output_path}_individual.fasta"
        extractor.export_individual_sequences(individual_file)
    
    if mode in ['consensus', 'both']:
        consensus_file = f"{output_path}_consensus.fasta"
        extractor.export_consensus_sequences(consensus_file)
    
    # Export details if requested
    if details:
        details_file = f"{output_path}_details.tsv"
        extractor.export_read_info(details_file)
    
    # Display summary
    extractor.display_summary()
    
    click.echo(f"\n✅ Sequence extraction complete!")
    
    # Suggest next steps
    click.echo(f"💡 Next steps:")
    if mode in ['individual', 'both']:
        click.echo(f"   • View individual: head {output_path}_individual.fasta")
    if mode in ['consensus', 'both']:
        click.echo(f"   • View consensus: head {output_path}_consensus.fasta")
        click.echo(f"   • Align consensus sequences with BLAST or similar")
    if details:
        click.echo(f"   • Analyze details: head {output_path}_details.tsv")


if __name__ == '__main__':
    extract()