#!/usr/bin/env python3
"""
IsoTag - Universal Isoform Tagger for BAM/SAM Files

Adds compact isoform structure (XI), reversible splicetag (XS), transcript group (XT),
and variant (XV) tags to BAM/SAM files. Uses RefGet-based chromosome hashing for
universal compatibility across chr1/Chr1/CHR1/1 naming conventions.

Tags added:
    XI:Z: Isoform structure ID (32-char hash, full exon coordinates)
    XB:Z: Reversible boundary tag (8-char chr hash + hex 5'/3' ends)
    XS:Z: Reversible splicetag (8-char chr hash + hex coordinates)
    XT:Z: Transcript group ID (32-char hash, position-based clustering)
    XV:Z: Variant ID (32-char hashes, if variant detection enabled)

Usage:
    python3 isotag.py -i input.bam -o tagged.bam
    python3 isotag.py -i input.bam -o tagged.bam --clustermode 5prime
    python3 isotag.py -i input.bam -o tagged.bam -r genome-refget.json
"""

import subprocess
import sys
import re
import hashlib
import base64
import click
import json
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import tempfile
import os

# Import VRS-compatible sha512t24u function
try:
    from vrs_compat import sha512t24u
except ImportError:
    # Fallback implementation if vrs_compat not available
    def sha512t24u(blob: bytes) -> str:
        """Generate base64url-encoded, truncated SHA-512 digest"""
        digest_size = 24
        digest = hashlib.sha512(blob).digest()
        tdigest_b64us = base64.urlsafe_b64encode(digest[:digest_size])
        return tdigest_b64us.decode("ascii")


@dataclass
class ExonBoundary:
    """Represents exon start and end coordinates"""
    start: int
    end: int

    def __lt__(self, other):
        return self.start < other.start


@dataclass
class SimpleVariant:
    """Simple variant representation using chr_hash:pos:ref>alt format"""
    chromosome: str
    position: int
    ref: str
    alt: str

    def to_string(self, chr_hash: str) -> str:
        """Convert to canonical string format with chromosome hash"""
        return f"{chr_hash}:{self.position}:{self.ref}>{self.alt}"

    def __lt__(self, other):
        return self.position < other.position


def mask_ambiguous_bases(sequence: str, keep_ambiguous: bool = False) -> str:
    """
    Mask ambiguous IUPAC codes with 'N' for consistent RefGet hashing

    This ensures UCSC hg38 and NCBI GRCh38 produce identical chromosome hashes
    despite differences in ambiguous base representation (R, Y, S, W, K, M, etc.)

    Args:
        sequence: DNA sequence string
        keep_ambiguous: If True, keep R/Y/etc as-is (may cause cross-genome incompatibility)

    Returns:
        Normalized sequence for hashing
    """
    if keep_ambiguous:
        return sequence.upper()

    # Convert ambiguous IUPAC bases to N
    canonical_bases = set('ACGTN')
    normalized = ''.join(base if base in canonical_bases else 'N'
                        for base in sequence.upper())

    return normalized


class RefGetCache:
    """Manages RefGet chromosome hash caching"""

    DEFAULT_CACHE_DIR = Path.home() / ".isotag_cache"

    @staticmethod
    def get_cache_path(genome_file: Optional[str]) -> Path:
        """Get cache file path for given genome"""
        if genome_file:
            genome_name = Path(genome_file).stem
        else:
            genome_name = "default"

        cache_dir = RefGetCache.DEFAULT_CACHE_DIR
        cache_dir.mkdir(exist_ok=True)
        return cache_dir / f"{genome_name}_refget.json"

    @staticmethod
    def load_or_generate(genome_file: Optional[str], refget_file: Optional[str],
                        keep_ambiguous: bool = False) -> Dict[str, str]:
        """Load RefGet mapping from file or generate from genome FASTA"""

        # If user provided RefGet JSON, use it
        if refget_file:
            click.echo(f"📂 Loading RefGet mapping from: {refget_file}")
            with open(refget_file, 'r') as f:
                data = json.load(f)
                return data.get("refget_mapping", {})

        # If no genome file, return empty mapping (will use legacy normalization)
        if not genome_file:
            click.echo("⚠️  No RefGet mapping or genome FASTA provided")
            click.echo("⚠️  Chromosome hashes will use legacy normalization (not recommended)")
            return {}

        # Check cache
        cache_path = RefGetCache.get_cache_path(genome_file)
        if cache_path.exists():
            click.echo(f"📂 Loading cached RefGet mapping from: {cache_path}")
            with open(cache_path, 'r') as f:
                data = json.load(f)
                return data.get("refget_mapping", {})

        # Generate RefGet mapping from genome FASTA
        click.echo(f"🧬 Generating RefGet mapping from genome: {genome_file}")
        click.echo(f"💾 This will be cached at: {cache_path}")

        mapping = RefGetCache.generate_from_fasta(genome_file, keep_ambiguous)

        # Save to cache
        RefGetCache.save_cache(mapping, cache_path, Path(genome_file).stem, keep_ambiguous)

        return mapping

    @staticmethod
    def generate_from_fasta(fasta_file: str, keep_ambiguous: bool = False) -> Dict[str, str]:
        """Generate RefGet mapping from FASTA file with optional ambiguous base masking"""
        mapping = {}
        current_chrom = None
        current_seq = []
        ambiguous_count_total = 0

        if not keep_ambiguous:
            click.echo("⏳ Processing chromosomes (masking ambiguous bases with 'N')...")
        else:
            click.echo("⏳ Processing chromosomes (keeping ambiguous IUPAC codes)...")

        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()

                if line.startswith('>'):
                    # Process previous chromosome
                    if current_chrom and current_seq:
                        seq = ''.join(current_seq)

                        # Normalize sequence (mask ambiguous bases unless --keep-ambiguous-bases)
                        normalized_seq = mask_ambiguous_bases(seq, keep_ambiguous)

                        # Count ambiguous bases if masking
                        if not keep_ambiguous:
                            ambiguous_in_chrom = sum(1 for orig, norm in zip(seq.upper(), normalized_seq)
                                                    if orig != norm)
                            ambiguous_count_total += ambiguous_in_chrom
                            if ambiguous_in_chrom > 0:
                                click.echo(f"   ⚠️  {current_chrom}: {ambiguous_in_chrom} ambiguous bases masked")

                        refget_id = sha512t24u(normalized_seq.encode('ascii'))
                        mapping[current_chrom] = f"SQ.{refget_id}"
                        click.echo(f"   ✅ {current_chrom} -> SQ.{refget_id[:8]}...")

                    # Start new chromosome
                    header = line[1:]
                    current_chrom = header.split()[0]
                    current_seq = []

                elif current_chrom:
                    current_seq.append(line)

            # Process last chromosome
            if current_chrom and current_seq:
                seq = ''.join(current_seq)

                # Normalize sequence
                normalized_seq = mask_ambiguous_bases(seq, keep_ambiguous)

                # Count ambiguous bases if masking
                if not keep_ambiguous:
                    ambiguous_in_chrom = sum(1 for orig, norm in zip(seq.upper(), normalized_seq)
                                            if orig != norm)
                    ambiguous_count_total += ambiguous_in_chrom
                    if ambiguous_in_chrom > 0:
                        click.echo(f"   ⚠️  {current_chrom}: {ambiguous_in_chrom} ambiguous bases masked")

                refget_id = sha512t24u(normalized_seq.encode('ascii'))
                mapping[current_chrom] = f"SQ.{refget_id}"
                click.echo(f"   ✅ {current_chrom} -> SQ.{refget_id[:8]}...")

        # Summary
        if not keep_ambiguous and ambiguous_count_total > 0:
            click.echo(f"   📊 Total ambiguous bases masked: {ambiguous_count_total:,}")
            click.echo(f"   ℹ️  This ensures UCSC hg38 and NCBI GRCh38 produce identical hashes")

        # Generate chromosome name variants
        extended = RefGetCache.generate_variants(mapping)
        click.echo(f"📊 Generated {len(mapping)} chromosomes -> {len(extended)} total mappings")

        return extended

    @staticmethod
    def generate_variants(base_mapping: Dict[str, str]) -> Dict[str, str]:
        """Generate all chromosome name variants (chr1, Chr1, CHR1, 1)"""
        extended = {}

        for chrom_name, refget_id in base_mapping.items():
            extended[chrom_name] = refget_id

            # Handle chr1 -> 1, Chr1, CHR1
            if chrom_name.lower().startswith('chr'):
                base = chrom_name[3:]
                extended[base] = refget_id
                extended[f"chr{base}"] = refget_id
                extended[f"Chr{base}"] = refget_id
                extended[f"CHR{base}"] = refget_id

            # Handle 1 -> chr1, Chr1, CHR1
            elif chrom_name.isalnum():
                extended[f"chr{chrom_name}"] = refget_id
                extended[f"Chr{chrom_name}"] = refget_id
                extended[f"CHR{chrom_name}"] = refget_id

        return extended

    @staticmethod
    def save_cache(mapping: Dict[str, str], cache_path: Path, genome_name: str,
                  keep_ambiguous: bool = False):
        """Save RefGet mapping to cache file"""
        from datetime import datetime

        data = {
            "metadata": {
                "genome": genome_name,
                "generated": datetime.now().isoformat(),
                "total_mappings": len(mapping),
                "ambiguous_bases_masked": not keep_ambiguous,
                "description": "IsoTag RefGet chromosome hash cache"
            },
            "refget_mapping": mapping
        }

        with open(cache_path, 'w') as f:
            json.dump(data, f, indent=2, sort_keys=True)

        click.echo(f"💾 Cached RefGet mapping: {cache_path}")


class IsoformTagger:
    """
    Universal Isoform Tagger with RefGet-based chromosome hashing

    Generates universal isoform IDs and adds them as BAM/SAM tags:
    - XI tag: 32-char structure hash (full RefGet hash + exon coordinates)
    - XB tag: Reversible boundary tag (8-char chr hash + hex 5'/3' ends)
    - XS tag: Reversible splicetag (8-char chr hash + hex coordinates)
    - XT tag: 32-char transcript group hash (mode-based clustering)
    - XV tag: 32-char variant hashes (full RefGet hash + variant info)
    """

    def __init__(self,
                 refget_mapping: Optional[Dict[str, str]] = None,
                 clustermode: str = "middle",
                 position_quantum: int = 10000,
                 span_quantum: int = 1000,
                 exon_quantum: int = 1000):
        """
        Initialize IsoformTagger

        Args:
            refget_mapping: RefGet chromosome hash mapping
            clustermode: Position for XT clustering (5prime, middle, 3prime)
            position_quantum: Bin size for position quantization (bp)
            span_quantum: Bin size for genomic span quantization (bp)
            exon_quantum: Bin size for exon length quantization (bp)
        """
        self.refget_mapping = refget_mapping or {}
        self.clustermode = clustermode
        self.position_quantum = position_quantum
        self.span_quantum = span_quantum
        self.exon_quantum = exon_quantum

    def get_chromosome_hash(self, chrom_name: str, hash_length: int = 32) -> str:
        """
        Get chromosome hash from RefGet ID

        Args:
            chrom_name: Chromosome name (chr1, 1, Chr1, etc.)
            hash_length: Hash length (8 for XS tag, 32 for others)

        Returns:
            Chromosome hash of specified length
        """
        # Try RefGet mapping first
        refget_id = self.refget_mapping.get(chrom_name)

        if refget_id and refget_id.startswith("SQ."):
            # Extract hash: "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2" -> "aKF498dA" or full
            full_hash = refget_id.split('.')[1]
            return full_hash[:hash_length]

        # Fallback: legacy normalization (not recommended)
        normalized = chrom_name.upper().replace('CHR', '')
        fallback_hash = self.generate_hash(normalized)
        return fallback_hash[:hash_length]

    def generate_hash(self, text: str) -> str:
        """Generate URL-safe base64 hash, 24 bytes for 32-char IDs (VRS format)"""
        return sha512t24u(text.encode('ascii'))

    def serialize_structure(self, chromosome: str, strand: str, exons: List[ExonBoundary]) -> str:
        """
        Serialize exon structure with 32-char chromosome hash
        Format: chr_hash_32|strand|exon1_start:exon1_end|exon2_start:exon2_end|...
        """
        chr_hash = self.get_chromosome_hash(chromosome, hash_length=32)
        exon_strings = [f"{e.start}:{e.end}" for e in sorted(exons)]
        return f"{chr_hash}|{strand}|{'|'.join(exon_strings)}"

    def encode_reversible_splicetag(self, chromosome: str, strand: str, exons: List[ExonBoundary]) -> Optional[str]:
        """
        Encode reversible splicetag with 8-char chromosome hash + hex coordinates
        Format: [8-CHAR-CHR-HASH][STRAND].hex1.hex2.hex3...

        Example: aKF498dAp.4b0.7d0.866.bb8

        Returns None for single-exon transcripts (no splice junctions)
        """
        sorted_exons = sorted(exons)

        # Need at least 2 exons for splice junctions
        if len(sorted_exons) < 2:
            return None

        # Get 8-character chromosome hash
        chr_hash_8 = self.get_chromosome_hash(chromosome, hash_length=8)

        # Encode strand: p(+) or m(-)
        strand_char = 'p' if strand == '+' else 'm'

        # Extract splice junction coordinates
        splice_coords = []
        for i in range(len(sorted_exons) - 1):
            exon_end = sorted_exons[i].end
            next_exon_start = sorted_exons[i + 1].start
            splice_coords.extend([exon_end, next_exon_start])

        # Convert to hex (lowercase, no padding)
        hex_coords = '.'.join(f'{coord:x}' for coord in splice_coords)

        # Build reversible splicetag
        return f"{chr_hash_8}{strand_char}.{hex_coords}"

    def decode_reversible_splicetag(self, encoded: str) -> Tuple[str, str, List[int]]:
        """
        Decode reversible splicetag back to components

        Args:
            encoded: Reversible splicetag string

        Returns:
            (chr_hash_8, strand, splice_coordinates)
        """
        chr_hash_8 = encoded[:8]
        strand = '+' if encoded[8] == 'p' else '-'
        hex_coords = encoded[10:]  # Skip chr_hash + strand + first dot

        coordinates = [int(coord, 16) for coord in hex_coords.split('.')]
        return chr_hash_8, strand, coordinates

    def encode_reversible_boundarytag(self, chromosome: str, strand: str, exons: List[ExonBoundary]) -> str:
        """
        Encode reversible boundary tag with 8-char chromosome hash + hex 5'/3' ends
        Format: [8-CHAR-CHR-HASH][STRAND].start_hex.end_hex

        Example: aKF498dAp.3e8.1004 (5' end = 1000, 3' end = 4100)

        This tag stores transcript termini and enables:
        - Full coordinate reconstruction when combined with XS tag
        - Cross-validation of XI structure hash
        - 5'/3' boundary clustering and degradation analysis
        """
        sorted_exons = sorted(exons)

        # Get 8-character chromosome hash
        chr_hash_8 = self.get_chromosome_hash(chromosome, hash_length=8)

        # Encode strand: p(+) or m(-)
        strand_char = 'p' if strand == '+' else 'm'

        # Get 5' and 3' transcript ends
        five_prime_end = sorted_exons[0].start
        three_prime_end = sorted_exons[-1].end

        # Convert to hex (lowercase, no padding)
        start_hex = f'{five_prime_end:x}'
        end_hex = f'{three_prime_end:x}'

        # Build reversible boundary tag
        return f"{chr_hash_8}{strand_char}.{start_hex}.{end_hex}"

    def decode_reversible_boundarytag(self, encoded: str) -> Tuple[str, str, int, int]:
        """
        Decode reversible boundary tag back to components

        Args:
            encoded: Reversible boundary tag string

        Returns:
            (chr_hash_8, strand, five_prime_end, three_prime_end)
        """
        chr_hash_8 = encoded[:8]
        strand = '+' if encoded[8] == 'p' else '-'
        hex_coords = encoded[10:]  # Skip chr_hash + strand + first dot

        coords = hex_coords.split('.')
        five_prime_end = int(coords[0], 16)
        three_prime_end = int(coords[1], 16)

        return chr_hash_8, strand, five_prime_end, three_prime_end

    def generate_structure_id(self, chromosome: str, strand: str, exons: List[ExonBoundary]) -> str:
        """Generate structure ID for XI tag (32-char hash)"""
        structure_serial = self.serialize_structure(chromosome, strand, exons)
        return self.generate_hash(structure_serial)

    def generate_splicetag_id(self, chromosome: str, strand: str, exons: List[ExonBoundary]) -> Optional[str]:
        """
        Generate reversible splicetag for XS tag
        Returns None for single-exon transcripts
        """
        return self.encode_reversible_splicetag(chromosome, strand, exons)

    def generate_boundarytag_id(self, chromosome: str, strand: str, exons: List[ExonBoundary]) -> str:
        """
        Generate reversible boundary tag for XB tag
        Always returns a value (works for single-exon and multi-exon)
        """
        return self.encode_reversible_boundarytag(chromosome, strand, exons)

    # Quantization functions
    def round_down_to_bin(self, value: int, bin_size: int) -> int:
        """Round value down to nearest bin boundary"""
        return (value // bin_size) * bin_size

    def round_up_to_bin(self, value: int, bin_size: int) -> int:
        """Round value up to nearest bin boundary"""
        return ((value + bin_size - 1) // bin_size) * bin_size

    def round_to_nearest_bin(self, value: int, bin_size: int) -> int:
        """Round value to nearest bin boundary"""
        return round(value / bin_size) * bin_size

    def calculate_transcript_metrics(self, exons: List[ExonBoundary]) -> Tuple[int, int, int, int]:
        """
        Calculate transcript metrics
        Returns: (genomic_start, genomic_end, exon_total_length, genomic_span)
        """
        sorted_exons = sorted(exons)

        genomic_start = sorted_exons[0].start
        genomic_end = sorted_exons[-1].end
        genomic_span = genomic_end - genomic_start
        exon_total_length = sum(exon.end - exon.start + 1 for exon in exons)

        return genomic_start, genomic_end, exon_total_length, genomic_span

    def calculate_cluster_position(self, genomic_start: int, genomic_end: int, strand: str) -> int:
        """
        Calculate clustering position based on mode

        Args:
            genomic_start: Transcript start position
            genomic_end: Transcript end position
            strand: Strand (+ or -)

        Returns:
            Position to use for clustering
        """
        if self.clustermode == "5prime":
            # 5prime (TSS): + strand uses start, - strand uses end
            if strand == "+":
                return self.round_down_to_bin(genomic_start, self.position_quantum)
            else:
                return self.round_up_to_bin(genomic_end, self.position_quantum)

        elif self.clustermode == "3prime":
            # 3prime (TES): + strand uses end, - strand uses start
            if strand == "+":
                return self.round_up_to_bin(genomic_end, self.position_quantum)
            else:
                return self.round_down_to_bin(genomic_start, self.position_quantum)

        else:  # middle (default)
            # Middle position: use average for both strands
            middle = (genomic_start + genomic_end) / 2
            return self.round_to_nearest_bin(int(middle), self.position_quantum)

    def serialize_transcript_group_xt(self, chromosome: str, strand: str, exons: List[ExonBoundary]) -> str:
        """
        Serialize transcript for XT tag with mode-based clustering
        Format: chr_hash_32|strand|quantized_position|quantized_exon_length|quantized_span|splice_junctions
        """
        chr_hash = self.get_chromosome_hash(chromosome, hash_length=32)

        # Calculate transcript metrics
        genomic_start, genomic_end, exon_total_length, genomic_span = self.calculate_transcript_metrics(exons)

        # Calculate clustering position based on mode
        cluster_position = self.calculate_cluster_position(genomic_start, genomic_end, strand)

        # Quantize metrics
        quantized_exon_length = self.round_to_nearest_bin(exon_total_length, self.exon_quantum)
        quantized_span = self.round_to_nearest_bin(genomic_span, self.span_quantum)

        # Build serialization
        base_parts = [
            chr_hash,
            strand,
            str(cluster_position),
            str(quantized_exon_length),
            str(quantized_span)
        ]

        # Add splice junctions if multi-exon
        sorted_exons = sorted(exons)
        if len(sorted_exons) > 1:
            splice_coords = []
            for i in range(len(sorted_exons) - 1):
                exon_end = sorted_exons[i].end
                next_exon_start = sorted_exons[i + 1].start
                splice_coords.extend([str(exon_end), str(next_exon_start)])
            base_parts.extend(splice_coords)

        return '|'.join(base_parts)

    def generate_transcript_group_xt_id(self, chromosome: str, strand: str, exons: List[ExonBoundary]) -> str:
        """Generate XT tag ID with mode-based clustering"""
        xt_serial = self.serialize_transcript_group_xt(chromosome, strand, exons)
        return self.generate_hash(xt_serial)

    def generate_individual_variant_ids(self, chromosome: str, variants: List[SimpleVariant]) -> Optional[str]:
        """
        Generate individual IDs for each variant with 32-char chromosome hash
        Returns: concatenated variant IDs like '32chars.32chars.32chars' or None if no variants
        """
        if not variants:
            return None

        chr_hash = self.get_chromosome_hash(chromosome, hash_length=32)

        individual_ids = []
        for variant in sorted(variants):
            # Create ID with chromosome hash: chr_hash:pos:ref>alt
            variant_string = variant.to_string(chr_hash)
            variant_hash = self.generate_hash(variant_string)
            individual_ids.append(variant_hash)

        return '.'.join(individual_ids)

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
                pass
            elif op_code == 3:  # N (skipped region/intron)
                if in_exon:
                    exons.append(ExonBoundary(exon_start, current_pos - 1))
                    in_exon = False
                current_pos += length
            elif op_code in [4, 5]:  # S, H (soft/hard clipping)
                pass

        # Add final exon
        if in_exon and exon_start < current_pos:
            exons.append(ExonBoundary(exon_start, current_pos - 1))

        return exons

    def parse_md_tag(self, md_tag: str, ref_pos: int) -> List[SimpleVariant]:
        """Parse MD tag to extract mismatches and deletions"""
        variants = []
        if not md_tag:
            return variants

        current_pos = ref_pos
        i = 0

        while i < len(md_tag):
            if md_tag[i].isdigit():
                num_str = ""
                while i < len(md_tag) and md_tag[i].isdigit():
                    num_str += md_tag[i]
                    i += 1
                if num_str:
                    current_pos += int(num_str)

            elif md_tag[i] == '^':
                i += 1
                deleted_seq = ""
                while i < len(md_tag) and md_tag[i].isalpha():
                    deleted_seq += md_tag[i]
                    i += 1

                if deleted_seq:
                    variant = SimpleVariant(
                        chromosome="unknown",
                        position=current_pos,
                        ref=deleted_seq,
                        alt="-"
                    )
                    variants.append(variant)
                    current_pos += len(deleted_seq)

            elif md_tag[i].isalpha():
                ref_base = md_tag[i]
                variant = SimpleVariant(
                    chromosome="unknown",
                    position=current_pos,
                    ref=ref_base,
                    alt="N"
                )
                variants.append(variant)
                current_pos += 1
                i += 1
            else:
                i += 1

        return variants

    def extract_insertions_from_cigar(self, cigar_operations: List[Tuple[int, int]], ref_pos: int) -> List[SimpleVariant]:
        """Extract insertions from CIGAR operations"""
        variants = []
        current_pos = ref_pos

        for op_code, length in cigar_operations:
            if op_code == 1:  # Insertion
                variant = SimpleVariant(
                    chromosome="unknown",
                    position=current_pos,
                    ref="-",
                    alt="I" * length
                )
                variants.append(variant)
            elif op_code in [0, 2, 7, 8]:  # M, D, =, X
                current_pos += length
            elif op_code == 3:  # N
                current_pos += length

        return variants

    def process_sam_line(self, line: str, detect_variants: bool = True) -> Optional[Dict]:
        """
        Process a single SAM line and generate universal isoform/variant IDs

        Args:
            line: SAM format line
            detect_variants: Whether to detect and generate variant IDs

        Returns:
            Dictionary with qname, isoform_id, splicetag_id, xt_group_id, variant_id, original_line
        """
        fields = line.strip().split('\t')

        if len(fields) < 11:
            return None

        qname = fields[0]
        flag = int(fields[1])
        rname = fields[2]
        pos = int(fields[3])
        cigar = fields[5]

        # Skip unmapped reads
        if flag & 0x4 or cigar == "*":
            return None

        # Parse CIGAR to get exon structure
        cigar_ops = self.parse_cigar_operations(cigar)
        if not cigar_ops:
            return None

        # Determine strand
        strand = "-" if flag & 0x10 else "+"

        # Extract exons
        exons = self.extract_exons_from_cigar(pos, cigar_ops)
        if not exons:
            return None

        # Generate structure ID (XI tag) - 32-char hash
        structure_id = self.generate_structure_id(rname, strand, exons)

        # Generate reversible boundary tag (XB tag) - 8-char chr hash + hex 5'/3' ends
        boundarytag_id = self.generate_boundarytag_id(rname, strand, exons)

        # Generate reversible splicetag (XS tag) - 8-char chr hash + hex coords
        splicetag_id = self.generate_splicetag_id(rname, strand, exons)

        # Generate XT tag ID (mode-based clustering) - 32-char hash
        xt_group_id = self.generate_transcript_group_xt_id(rname, strand, exons)

        # Extract variants if enabled
        variant_id = None
        if detect_variants:
            variants = []

            # Look for MD tag
            md_tag = None
            for field in fields[11:]:
                if field.startswith('MD:Z:'):
                    md_tag = field[5:]
                    break

            if md_tag:
                md_variants = self.parse_md_tag(md_tag, pos)
                for variant in md_variants:
                    variant.chromosome = rname
                variants.extend(md_variants)

            # Extract insertions from CIGAR
            cigar_variants = self.extract_insertions_from_cigar(cigar_ops, pos)
            for variant in cigar_variants:
                variant.chromosome = rname
            variants.extend(cigar_variants)

            # Generate variant ID with 32-char chromosome hash
            variant_id = self.generate_individual_variant_ids(rname, variants)

        return {
            'qname': qname,
            'isoform_id': structure_id,
            'boundarytag_id': boundarytag_id,
            'splicetag_id': splicetag_id,
            'xt_group_id': xt_group_id,
            'variant_id': variant_id,
            'original_line': line.strip()
        }

    def check_for_md_tags(self, bam_file: str, sample_size: int = 100) -> bool:
        """Check if BAM file contains MD tags by sampling reads"""
        try:
            is_bam = Path(bam_file).suffix.lower() == '.bam'

            if is_bam:
                cmd = ['samtools', 'view', bam_file]
            else:
                cmd = ['cat', bam_file]

            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

            reads_checked = 0
            md_tags_found = 0

            for line in process.stdout:
                if line.startswith('@'):
                    continue

                reads_checked += 1
                if reads_checked > sample_size:
                    break

                fields = line.strip().split('\t')
                if len(fields) >= 11:
                    for field in fields[11:]:
                        if field.startswith('MD:Z:'):
                            md_tags_found += 1
                            break

            process.terminate()

            return md_tags_found > 0 and (md_tags_found / reads_checked) >= 0.1 if reads_checked > 0 else False

        except (subprocess.SubprocessError, FileNotFoundError):
            return False


@click.command()
@click.option('--input', '-i', 'input_file', required=True, help='Input BAM/SAM file')
@click.option('--output', '-o', required=True, help='Output BAM/SAM file with XI/XB/XS/XT/XV tags')
@click.option('--genome', '-g', help='Reference genome FASTA (for RefGet cache generation and/or variant detection)')
@click.option('--refget', '-r', help='RefGet JSON mapping file (optional, will auto-generate from genome if not provided)')
@click.option('--keep-ambiguous-bases', is_flag=True,
              help='Keep ambiguous IUPAC bases (R,Y,S,W,K,M,etc) in RefGet hashing (may cause incompatibility across genome versions)')
@click.option('--clustermode', type=click.Choice(['5prime', 'middle', '3prime']), default='middle',
              help='Position for XT clustering: 5prime=CAGE/TSS, middle=RNA-seq, 3prime=polyA/TES (default: middle)')
@click.option('--position-quantum', type=int, default=10000,
              help='Bin size for position quantization in bp (default: 10000)')
@click.option('--span-quantum', type=int, default=1000,
              help='Bin size for genomic span quantization in bp (default: 1000)')
@click.option('--exon-quantum', type=int, default=1000,
              help='Bin size for exon length quantization in bp (default: 1000)')
def isotag(input_file, output, genome, refget, keep_ambiguous_bases, clustermode, position_quantum, span_quantum, exon_quantum):
    """
    IsoTag - Universal Isoform Tagger

    Add RefGet-based universal tags to BAM/SAM files:
    - XI: Isoform structure (32-char hash with full RefGet chr hash)
    - XB: Reversible boundary tag (8-char chr hash + hex 5'/3' ends)
    - XS: Reversible splicetag (8-char chr hash + hex coordinates)
    - XT: Transcript group (32-char hash, mode-based clustering)
    - XV: Variant IDs (32-char hashes with full RefGet chr hash)

    RefGet Behavior:
    - If --refget provided: Use specified RefGet JSON mapping
    - If --genome provided: Auto-generate and cache RefGet mapping from FASTA
    - If neither: Use legacy chromosome normalization (not recommended)

    Clustering Modes (--clustermode):
    - 5prime: CAGE/TSS clustering (+ uses start, - uses end)
    - middle: RNA-seq clustering (uses transcript midpoint)
    - 3prime: polyA/TES clustering (+ uses end, - uses start)

    Variant Detection:
    - No MD tags + No genome: Only structure tags (XI/XS/XT)
    - Has MD tags: Add variant tags (XV)
    - No MD tags + Has genome: Generate MD tags, add XV tags

    Examples:
        # Basic tagging with auto RefGet cache
        isotag.py -i input.bam -o tagged.bam -g genome.fa

        # CAGE data (5' clustering)
        isotag.py -i cage.bam -o tagged.bam -g genome.fa --clustermode 5prime

        # Use pre-computed RefGet mapping
        isotag.py -i input.bam -o tagged.bam -r genome-refget.json

        # Custom quantization
        isotag.py -i input.bam -o tagged.bam -g genome.fa --position-quantum 5000
    """

    input_path = Path(input_file)
    output_path = Path(output)

    if not input_path.exists():
        click.echo(f"❌ Input file not found: {input_path}")
        sys.exit(1)

    # Detect file formats
    input_is_bam = input_path.suffix.lower() == '.bam'
    output_is_bam = output_path.suffix.lower() == '.bam'

    click.echo(f"🧬 IsoTag v8.2 - Reversible Boundary Tags Edition")
    click.echo(f"📥 Input: {input_path.name} ({'BAM' if input_is_bam else 'SAM'})")
    click.echo(f"📤 Output: {output_path.name} ({'BAM' if output_is_bam else 'SAM'})")
    click.echo(f"🎯 Cluster mode: {clustermode}")
    click.echo(f"📏 Quantization: position={position_quantum}bp, span={span_quantum}bp, exon={exon_quantum}bp")

    # Load or generate RefGet mapping
    refget_mapping = RefGetCache.load_or_generate(genome, refget, keep_ambiguous_bases)

    if refget_mapping:
        click.echo(f"✅ RefGet mapping loaded: {len(refget_mapping)} chromosome mappings")
    else:
        click.echo("⚠️  Using legacy chromosome normalization (consider providing genome FASTA)")

    # Initialize tagger
    tagger = IsoformTagger(
        refget_mapping=refget_mapping,
        clustermode=clustermode,
        position_quantum=position_quantum,
        span_quantum=span_quantum,
        exon_quantum=exon_quantum
    )

    # Statistics
    stats = {
        'total_reads': 0,
        'reads_processed': 0,
        'reads_with_structure': 0,
        'reads_with_boundarytags': 0,
        'reads_with_splicetags': 0,
        'reads_with_xt_groups': 0,
        'reads_with_variants': 0,
        'unique_structures': set(),
        'unique_boundarytags': set(),
        'unique_splicetags': set(),
        'unique_xt_groups': set(),
        'unique_variants': set()
    }

    try:
        # Determine variant detection strategy
        detect_variants = False
        input_for_processing = str(input_path)
        needs_cleanup = False

        if genome:
            # Generate MD tags for variant detection
            click.echo("📋 Preparing BAM with MD tags for variant detection...")
            with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as temp_bam:
                temp_bam_path = temp_bam.name

            subprocess.run(['samtools', 'calmd', '-b', str(input_path), genome],
                         stdout=open(temp_bam_path, 'wb'), check=True)
            input_for_processing = temp_bam_path
            detect_variants = True
            needs_cleanup = True
            click.echo("🔬 Variant detection: Enabled")
        else:
            # Check for existing MD tags
            has_md_tags = tagger.check_for_md_tags(str(input_path))

            if has_md_tags:
                detect_variants = True
                click.echo("🔬 Variant detection: Enabled (found MD tags)")
            else:
                detect_variants = False
                click.echo("🔬 Variant detection: Disabled")

        click.echo("🏷️  Adding isoform tags (XI, XB, XS, XT)...")
        if detect_variants:
            click.echo("🏷️  Adding variant tags (XV)...")

        # Process input
        if input_is_bam:
            samtools_cmd = ['samtools', 'view', '-h', input_for_processing]
        else:
            samtools_cmd = ['cat', input_for_processing]

        # Process and add tags
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as temp_output:
            temp_output_path = temp_output.name

            with subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE, text=True) as proc:
                for line in proc.stdout:
                    if line.startswith('@'):
                        temp_output.write(line)
                        continue

                    stats['total_reads'] += 1

                    if stats['total_reads'] % 10000 == 0:
                        click.echo(f"⏳ Processed {stats['total_reads']:,} reads...")

                    # Process read
                    result = tagger.process_sam_line(line, detect_variants)

                    if result:
                        stats['reads_processed'] += 1
                        fields = line.strip().split('\t')

                        # Add XI tag (structure, 32-char hash)
                        fields.append(f"XI:Z:{result['isoform_id']}")
                        stats['reads_with_structure'] += 1
                        stats['unique_structures'].add(result['isoform_id'])

                        # Add XB tag (reversible boundary tag, 8-char chr hash + hex ends)
                        fields.append(f"XB:Z:{result['boundarytag_id']}")
                        stats['reads_with_boundarytags'] += 1
                        stats['unique_boundarytags'].add(result['boundarytag_id'])

                        # Add XS tag (reversible splicetag, 8-char chr hash + hex)
                        if result['splicetag_id']:
                            fields.append(f"XS:Z:{result['splicetag_id']}")
                            stats['reads_with_splicetags'] += 1
                            stats['unique_splicetags'].add(result['splicetag_id'])

                        # Add XT tag (transcript group, 32-char hash)
                        fields.append(f"XT:Z:{result['xt_group_id']}")
                        stats['reads_with_xt_groups'] += 1
                        stats['unique_xt_groups'].add(result['xt_group_id'])

                        # Add XV tag (variants, 32-char hashes)
                        if detect_variants and result['variant_id']:
                            fields.append(f"XV:Z:{result['variant_id']}")
                            stats['reads_with_variants'] += 1
                            stats['unique_variants'].add(result['variant_id'])

                        temp_output.write('\t'.join(fields) + '\n')
                    else:
                        temp_output.write(line)

        # Convert to final format
        if output_is_bam:
            click.echo("🔄 Converting to BAM format...")
            subprocess.run(['samtools', 'view', '-b', temp_output_path, '-o', str(output_path)],
                         check=True)
            os.unlink(temp_output_path)
        else:
            import shutil
            shutil.move(temp_output_path, str(output_path))

        # Cleanup
        if needs_cleanup and 'temp_bam_path' in locals():
            os.unlink(temp_bam_path)

        # Display results
        click.echo("\n" + "="*60)
        click.echo("✅ IsoTag Complete!")
        click.echo("="*60)
        click.echo(f"📊 Total reads: {stats['total_reads']:,}")
        click.echo(f"🧬 Reads processed: {stats['reads_processed']:,}")
        click.echo(f"🏷️  XI tags (structure): {stats['reads_with_structure']:,}")
        click.echo(f"🔚 XB tags (boundaries): {stats['reads_with_boundarytags']:,}")
        click.echo(f"🔗 XS tags (splicetag): {stats['reads_with_splicetags']:,}")
        click.echo(f"🎯 XT tags (transcript group): {stats['reads_with_xt_groups']:,}")
        click.echo(f"🔬 XV tags (variants): {stats['reads_with_variants']:,}")
        click.echo(f"🆔 Unique structures: {len(stats['unique_structures']):,}")
        click.echo(f"🔚 Unique boundarytags: {len(stats['unique_boundarytags']):,}")
        click.echo(f"🧬 Unique splicetags: {len(stats['unique_splicetags']):,}")
        click.echo(f"🎯 Unique XT groups: {len(stats['unique_xt_groups']):,}")
        click.echo(f"🧪 Unique variant combos: {len(stats['unique_variants']):,}")
        click.echo(f"💾 Output: {output_path}")

        # Show example tags
        click.echo(f"\n🎯 Example tags added:")
        if stats['unique_structures']:
            example = next(iter(stats['unique_structures']))
            click.echo(f"   XI:Z:{example} (32-char structure)")
        if stats['unique_boundarytags']:
            example = next(iter(stats['unique_boundarytags']))
            click.echo(f"   XB:Z:{example} (reversible boundary tag)")
        if stats['unique_splicetags']:
            example = next(iter(stats['unique_splicetags']))
            click.echo(f"   XS:Z:{example[:50]}... (reversible splicetag)")
        if stats['unique_xt_groups']:
            example = next(iter(stats['unique_xt_groups']))
            click.echo(f"   XT:Z:{example} (32-char {clustermode} group)")
        if stats['unique_variants']:
            example = next(iter(stats['unique_variants']))
            click.echo(f"   XV:Z:{example[:50]}... (variant IDs)")

        # Viewing command
        if output_is_bam:
            click.echo(f"\n🚀 View results: samtools view {output_path} | grep 'XI:Z:' | head")
        else:
            click.echo(f"\n🚀 View results: grep 'XI:Z:' {output_path} | head")

    except subprocess.CalledProcessError as e:
        click.echo(f"❌ Error running samtools: {e}")
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    isotag()
