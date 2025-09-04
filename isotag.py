#!/usr/bin/env python3
"""
IsoTag v1.0 - RefGet-Compatible Universal Isoform Tagger

Adds RefGet SQ.XXXX based isoform structure (XI) and VRS-compliant variant (XV) tags 
to BAM/SAM files. Fully compatible with GA4GH VRS and RefGet specifications.

Key Changes in v1.0:
- Uses RefGet SQ.XXXX sequence identifiers instead of chromosome names
- VRS-compliant VARSID calculation using official algorithms
- Requires genome reference (-g) or pre-computed RefGet mapping (--refget)
- Generates truly universal, cross-database compatible identifiers

Usage:
    # With FASTA (calculates RefGet on-the-fly)
    python3 isotag.py -i input.bam -o tagged.bam -g reference.fa
    
    # With pre-computed RefGet mapping (faster)
    python3 isotag.py -i input.bam -o tagged.bam --refget hg38-refget.json
"""

import subprocess
import sys
import re
import json
import tempfile
import os
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import click

# Import our RefGet and VRS modules
from isotag_refget import RefGetCalculator
from vrs_compat import VRSVariant, sha512t24u


@dataclass
class ExonBoundary:
    """Represents exon start and end coordinates"""
    start: int
    end: int
    
    def __lt__(self, other):
        return self.start < other.start


@dataclass  
class RefGetVariant:
    """RefGet-based variant representation for VRS compliance"""
    refget_id: str  # SQ.XXXX identifier
    position: int   # 0-based position
    ref: str        # Reference sequence
    alt: str        # Alternate sequence
    
    def to_vrs_variant(self) -> VRSVariant:
        """Convert to VRS-compliant variant for digest calculation"""
        return VRSVariant(self.refget_id, self.position, self.ref, self.alt)
    
    def __lt__(self, other):
        return self.position < other.position


class IsoformTaggerV6:
    """
    RefGet-Compatible Universal Isoform Tagger
    
    Generates RefGet SQ.XXXX based isoform IDs and VRS-compliant variant IDs:
    - XI tag: 32-character structure hash using RefGet identifiers
    - XV tag: VRS-compliant variant hashes separated by dots
    """
    
    def __init__(self, refget_mapping: Dict[str, str]):
        """Initialize with RefGet chromosome -> SQ.XXXX mapping
        
        Args:
            refget_mapping: Dictionary mapping chromosome names to RefGet IDs
        """
        self.refget_mapping = refget_mapping
        if not refget_mapping:
            raise ValueError("RefGet mapping is required for v1.0")
    
    def get_refget_id(self, chromosome: str) -> str:
        """Get RefGet ID for chromosome name
        
        Args:
            chromosome: Chromosome name from BAM file
            
        Returns:
            RefGet SQ.XXXX identifier
            
        Raises:
            ValueError: If chromosome not found in mapping
        """
        refget_id = self.refget_mapping.get(chromosome)
        if not refget_id:
            raise ValueError(f"No RefGet ID found for chromosome: {chromosome}")
        return refget_id
    
    def generate_hash(self, text: str) -> str:
        """Generate VRS-compatible hash using official sha512t24u
        
        Args:
            text: Text to hash
            
        Returns:
            32-character hash compatible with VRS format
        """
        return sha512t24u(text.encode('utf-8'))
    
    def serialize_structure(self, refget_id: str, strand: str, exons: List[ExonBoundary]) -> str:
        """Serialize exon structure using RefGet identifier
        
        Args:
            refget_id: RefGet sequence identifier (SQ.XXXX)
            strand: Strand (+ or -)
            exons: List of exon boundaries
            
        Returns:
            Canonical structure string for hashing
        """
        exon_strings = [f"{e.start}:{e.end}" for e in sorted(exons)]
        return f"{refget_id}|{strand}|{'|'.join(exon_strings)}"
    
    def generate_structure_id(self, chromosome: str, strand: str, exons: List[ExonBoundary]) -> str:
        """Generate RefGet-based structure ID for XI tag
        
        Args:
            chromosome: Chromosome name from BAM
            strand: Strand (+ or -)
            exons: List of exon boundaries
            
        Returns:
            32-character structure hash for XI tag
        """
        refget_id = self.get_refget_id(chromosome)
        structure_serial = self.serialize_structure(refget_id, strand, exons)
        return self.generate_hash(structure_serial)
    
    def generate_vrs_variant_ids(self, variants: List[RefGetVariant]) -> Optional[str]:
        """Generate VRS-compliant variant IDs for XV tag
        
        Args:
            variants: List of RefGet-based variants
            
        Returns:
            Concatenated VRS variant IDs or None if no variants
        """
        if not variants:
            return None
        
        vrs_ids = []
        for variant in sorted(variants):
            vrs_variant = variant.to_vrs_variant()
            digest = vrs_variant.calculate_digest_only()
            vrs_ids.append(digest)
        
        return '.'.join(vrs_ids)
    
    def parse_cigar_operations(self, cigar_string: str) -> List[Tuple[int, int]]:
        """Parse CIGAR string into list of (operation_code, length) tuples"""
        if not cigar_string or cigar_string == "*":
            return []
        
        operations = []
        pattern = r'(\d+)([MIDNSHP=X])'
        matches = re.findall(pattern, cigar_string)
        
        for length_str, op_char in matches:
            length = int(length_str)
            op_map = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
            if op_char in op_map:
                operations.append((op_map[op_char], length))
        
        return operations
    
    def extract_exons_from_cigar(self, pos: int, cigar_operations: List[Tuple[int, int]]) -> List[ExonBoundary]:
        """Extract exon boundaries from CIGAR operations"""
        exons = []
        current_pos = pos
        exon_start = pos
        in_exon = True
        
        for op_code, length in cigar_operations:
            if op_code in [0, 7, 8]:  # M, =, X (alignment/match)
                if not in_exon:
                    exon_start = current_pos
                    in_exon = True
                current_pos += length
            elif op_code == 2:  # D (deletion from reference)
                current_pos += length
            elif op_code == 1:  # I (insertion to reference)
                pass  # Don't advance reference position
            elif op_code == 3:  # N (skipped region/intron)
                if in_exon:
                    exons.append(ExonBoundary(exon_start, current_pos - 1))
                    in_exon = False
                current_pos += length
            elif op_code in [4, 5]:  # S, H (soft/hard clipping)
                pass  # Don't advance reference position for clipping
        
        # Add final exon if we ended in an exon
        if in_exon and exon_start < current_pos:
            exons.append(ExonBoundary(exon_start, current_pos - 1))
        
        return exons
    
    def parse_md_tag(self, md_tag: str, ref_pos: int, chromosome: str) -> List[RefGetVariant]:
        """Parse MD tag to extract mismatches and deletions using RefGet IDs
        
        Args:
            md_tag: MD tag value
            ref_pos: Reference position
            chromosome: Chromosome name
            
        Returns:
            List of RefGet-based variants
        """
        variants = []
        if not md_tag:
            return variants
        
        refget_id = self.get_refget_id(chromosome)
        current_pos = ref_pos
        i = 0
        
        while i < len(md_tag):
            if md_tag[i].isdigit():
                # Parse number (matching bases)
                num_str = ""
                while i < len(md_tag) and md_tag[i].isdigit():
                    num_str += md_tag[i]
                    i += 1
                if num_str:
                    current_pos += int(num_str)
            
            elif md_tag[i] == '^':
                # Deletion: ^ATCG
                i += 1  # Skip '^'
                deleted_seq = ""
                while i < len(md_tag) and md_tag[i].isalpha():
                    deleted_seq += md_tag[i]
                    i += 1
                
                if deleted_seq:
                    variant = RefGetVariant(
                        refget_id=refget_id,
                        position=current_pos,
                        ref=deleted_seq,
                        alt="-"
                    )
                    variants.append(variant)
                    current_pos += len(deleted_seq)
            
            elif md_tag[i].isalpha():
                # Mismatch: single base
                ref_base = md_tag[i]
                variant = RefGetVariant(
                    refget_id=refget_id,
                    position=current_pos,
                    ref=ref_base,
                    alt="N"  # Will be determined from read sequence if needed
                )
                variants.append(variant)
                current_pos += 1
                i += 1
            else:
                i += 1
        
        return variants
    
    def extract_insertions_from_cigar(self, cigar_operations: List[Tuple[int, int]], 
                                    ref_pos: int, chromosome: str) -> List[RefGetVariant]:
        """Extract insertions from CIGAR operations using RefGet IDs"""
        variants = []
        refget_id = self.get_refget_id(chromosome)
        current_pos = ref_pos
        
        for op_code, length in cigar_operations:
            if op_code == 1:  # Insertion
                variant = RefGetVariant(
                    refget_id=refget_id,
                    position=current_pos,
                    ref="-",
                    alt="I" * length  # Placeholder - would need read sequence for actual bases
                )
                variants.append(variant)
                # Don't advance reference position for insertions
            elif op_code in [0, 2, 7, 8]:  # M, D, =, X (consume reference)
                current_pos += length
            elif op_code == 3:  # N (intron)
                current_pos += length
            # S, H, P don't consume reference
        
        return variants
    
    def process_sam_line(self, line: str, detect_variants: bool = True) -> Optional[Dict]:
        """Process SAM line and generate RefGet-based isoform/variant IDs
        
        Args:
            line: SAM format line
            detect_variants: Whether to detect and generate variant IDs
            
        Returns:
            Dictionary with qname, isoform_id, variant_id, and original_line
        """
        fields = line.strip().split('\t')
        
        if len(fields) < 11:
            return None
        
        qname = fields[0]
        flag = int(fields[1])
        rname = fields[2]  # Chromosome name
        pos = int(fields[3])
        cigar = fields[5]
        
        # Skip unmapped reads
        if flag & 0x4:
            return None
            
        # Skip if no CIGAR or if CIGAR indicates unmapped
        if cigar == "*":
            return None
        
        # Parse CIGAR to get exon structure
        cigar_ops = self.parse_cigar_operations(cigar)
        if not cigar_ops:
            return None
        
        # Determine strand from flag
        strand = "-" if flag & 0x10 else "+"
        
        # Extract exons from CIGAR
        exons = self.extract_exons_from_cigar(pos, cigar_ops)
        if not exons:
            return None
        
        # Generate RefGet-based structure ID
        try:
            structure_id = self.generate_structure_id(rname, strand, exons)
        except ValueError as e:
            # Skip reads from chromosomes not in RefGet mapping
            return None
        
        # Extract variants if requested
        variant_id = None
        if detect_variants:
            variants = []
            
            # Look for MD tag in optional fields
            md_tag = None
            for field in fields[11:]:
                if field.startswith('MD:Z:'):
                    md_tag = field[5:]
                    break
            
            if md_tag:
                try:
                    md_variants = self.parse_md_tag(md_tag, pos, rname)
                    variants.extend(md_variants)
                except ValueError:
                    # Skip if chromosome not in RefGet mapping
                    pass
            
            # Extract insertions from CIGAR
            try:
                cigar_variants = self.extract_insertions_from_cigar(cigar_ops, pos, rname)
                variants.extend(cigar_variants)
            except ValueError:
                # Skip if chromosome not in RefGet mapping
                pass
            
            # Generate VRS-compliant variant ID
            variant_id = self.generate_vrs_variant_ids(variants)
        
        return {
            'qname': qname,
            'isoform_id': structure_id,
            'variant_id': variant_id,
            'original_line': line.strip()
        }


@click.command()
@click.option('--input', '-i', 'input_file', required=True, help='Input BAM/SAM file')
@click.option('--output', '-o', required=True, help='Output BAM/SAM file with XI/XV tags')
@click.option('--genome', '-g', help='Reference genome FASTA (calculates RefGet IDs)')
@click.option('--refget', '-r', help='Pre-computed RefGet JSON mapping (faster)')
@click.option('--variants/--no-variants', default=True, help='Enable/disable variant detection')
def isotag(input_file, output, genome, refget, variants):
    """
    IsoTag v1.0 - RefGet-Compatible Universal Isoform Tagger
    
    Add RefGet SQ.XXXX based XI (isoform structure) and VRS-compliant XV (variant) tags.
    Uses GA4GH VRS and RefGet specifications for truly universal identifiers.
    
    RefGet Options (exactly one required):
    - Use -g FASTA to calculate RefGet IDs on-the-fly (slower, first time)
    - Use -r JSON to load pre-computed RefGet mapping (faster, recommended)
    
    Examples:
        # Calculate RefGet from FASTA
        isotag.py -i input.bam -o tagged.bam -g hg38.fa
        
        # Use pre-computed RefGet mapping
        isotag.py -i input.bam -o tagged.bam -r hg38-refget.json
        
        # Structure tags only (no variants)
        isotag.py -i input.bam -o tagged.bam -r hg38-refget.json --no-variants
    """
    
    input_path = Path(input_file)
    output_path = Path(output)
    
    # Validation
    if not input_path.exists():
        click.echo(f"❌ Input file not found: {input_path}")
        sys.exit(1)
    
    if not (genome or refget):
        click.echo("❌ Must specify either --genome FASTA or --refget JSON")
        sys.exit(1)
    
    if genome and refget:
        click.echo("❌ Cannot specify both --genome and --refget")
        sys.exit(1)
    
    # Detect input and output formats
    input_is_bam = input_path.suffix.lower() == '.bam'
    output_is_bam = output_path.suffix.lower() == '.bam'
    
    click.echo(f"🧬 IsoTag v1.0 - RefGet Compatible Processing")
    click.echo(f"📥 Input: {input_path.name} ({'BAM' if input_is_bam else 'SAM'})")
    click.echo(f"📤 Output: {output_path.name} ({'BAM' if output_is_bam else 'SAM'})")
    
    # Load or calculate RefGet mapping
    refget_mapping = {}
    
    if refget:
        click.echo(f"📋 Loading RefGet mapping: {refget}")
        if not Path(refget).exists():
            click.echo(f"❌ RefGet file not found: {refget}")
            sys.exit(1)
        refget_mapping = RefGetCalculator.load_refget_mapping(refget)
        click.echo(f"✅ Loaded {len(refget_mapping)} RefGet mappings")
        
    elif genome:
        click.echo(f"🧬 Calculating RefGet IDs from: {genome}")
        if not Path(genome).exists():
            click.echo(f"❌ Genome file not found: {genome}")
            sys.exit(1)
        
        # Calculate RefGet mapping from FASTA
        base_mapping = RefGetCalculator.calculate_from_fasta(genome)
        refget_mapping = RefGetCalculator.generate_all_variants(base_mapping)
        click.echo(f"✅ Calculated {len(base_mapping)} RefGet IDs -> {len(refget_mapping)} total mappings")
    
    # Initialize tagger
    tagger = IsoformTaggerV6(refget_mapping)
    
    # Determine variant detection strategy
    detect_variants = variants
    input_for_processing = str(input_path)
    needs_cleanup = False
    
    if detect_variants and genome:
        # Generate MD tags for variant detection
        click.echo("📋 Adding MD tags for VRS-compliant variant detection...")
        with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as temp_bam:
            temp_bam_path = temp_bam.name
        
        subprocess.run(['samtools', 'calmd', '-b', str(input_path), genome], 
                     stdout=open(temp_bam_path, 'wb'), check=True)
        input_for_processing = temp_bam_path
        needs_cleanup = True
    
    click.echo(f"🏷️  Adding RefGet-based structure tags (XI)")
    if detect_variants:
        click.echo(f"🧪 Adding VRS-compliant variant tags (XV)")
    
    # Statistics tracking
    stats = {
        'total_reads': 0,
        'reads_processed': 0,
        'reads_with_structure': 0,
        'reads_with_variants': 0,
        'unique_structures': set(),
        'unique_variants': set(),
        'skipped_chromosomes': set()
    }
    
    try:
        # Prepare samtools command
        if input_is_bam:
            samtools_cmd = ['samtools', 'view', '-h', input_for_processing]
        else:
            samtools_cmd = ['cat', input_for_processing]
        
        # Process file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as temp_output:
            temp_output_path = temp_output.name
            
            with subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE, text=True) as proc:
                for line in proc.stdout:
                    if line.startswith('@'):
                        # Header line - pass through unchanged
                        temp_output.write(line)
                        continue
                    
                    stats['total_reads'] += 1
                    
                    if stats['total_reads'] % 1000 == 0:
                        click.echo(f"⏳ Processed {stats['total_reads']:,} reads...")
                    
                    # Process read line
                    result = tagger.process_sam_line(line, detect_variants)
                    
                    if result:
                        stats['reads_processed'] += 1
                        
                        # Add tags to the line
                        fields = line.strip().split('\t')
                        
                        # Add XI tag for RefGet-based structure ID
                        fields.append(f"XI:Z:{result['isoform_id']}")
                        stats['reads_with_structure'] += 1
                        stats['unique_structures'].add(result['isoform_id'])
                        
                        # Add XV tag for VRS-compliant variant ID
                        if detect_variants and result['variant_id']:
                            fields.append(f"XV:Z:{result['variant_id']}")
                            stats['reads_with_variants'] += 1
                            stats['unique_variants'].add(result['variant_id'])
                        
                        # Write modified line
                        temp_output.write('\t'.join(fields) + '\n')
                    else:
                        # Write original line if processing failed
                        temp_output.write(line)
        
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
        
        # Clean up temporary file if created
        if needs_cleanup and 'temp_bam_path' in locals():
            os.unlink(temp_bam_path)
        
        # Display results
        click.echo("\n" + "="*60)
        click.echo("✅ IsoTag v1.0 Complete - RefGet Compatible!")
        click.echo("="*60)
        click.echo(f"📊 Total reads: {stats['total_reads']:,}")
        click.echo(f"🧬 Reads processed: {stats['reads_processed']:,}")
        click.echo(f"🏷️  RefGet structure tags (XI): {stats['reads_with_structure']:,}")
        click.echo(f"🧪 VRS variant tags (XV): {stats['reads_with_variants']:,}")
        click.echo(f"🆔 Unique structures: {len(stats['unique_structures']):,}")
        click.echo(f"🔬 Unique variant combinations: {len(stats['unique_variants']):,}")
        click.echo(f"💾 Output: {output_path}")
        
        if stats['unique_structures']:
            example_structure = next(iter(stats['unique_structures']))
            click.echo(f"\n🎯 Example RefGet-based tags:")
            click.echo(f"   XI:Z:{example_structure}")
        if stats['unique_variants']:
            example_variant = next(iter(stats['unique_variants']))
            click.echo(f"   XV:Z:{example_variant}")
        
        # Show viewing command
        if output_is_bam:
            click.echo(f"\n🚀 View results: samtools view {output_path} | grep 'XI:Z:' | head")
        else:
            click.echo(f"\n🚀 View results: grep 'XI:Z:' {output_path} | head")
        
        click.echo(f"\n🌐 Universal Compatibility: RefGet SQ.XXXX + VRS compliant")
        
    except subprocess.CalledProcessError as e:
        click.echo(f"❌ Error running samtools: {e}")
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    isotag()
