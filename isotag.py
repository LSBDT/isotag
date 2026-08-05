#!/usr/bin/env python3
"""
IsoTag - Universal Isoform Tagger for BAM/SAM Files

Adds compact isoform structure (XI), reversible splicetag (XS), transcript group (XT),
cluster (XC), and variant (XV) tags to BAM/SAM files. Uses RefGet-based chromosome hashing
for universal compatibility across chr1/Chr1/CHR1/1 naming conventions.

Tags added:
    XI:Z: Isoform structure ID (32-char hash, full exon coordinates)
    XB:Z: Reversible boundary tag (8-char chr hash + hex 5'/3' ends)
    XS:Z: Reversible splicetag (8-char chr hash + hex coordinates)
    XT:Z: Transcript group ID (32-char hash, position-based clustering)
    XC:Z: Locus ID (spatial bin tag - chr+strand+binned positions only)
    X5:Z: TSS cluster ID (32-char hash, binned 5' transcript end)
    X3:Z: PolyA cluster ID (32-char hash, binned 3' transcript end)
    XV:Z: GA4GH VRS v2 Allele IDs (SNVs/substitutions, comma-separated, if variant detection enabled)

Usage:
    python3 isotag.py -i input.bam -o tagged.bam
    python3 isotag.py -i input.bam -o tagged.bam --clustermode 5prime
    python3 isotag.py -i input.bam -o tagged.bam -r genome-refget.json
    python3 isotag.py -i input.bam -o tagged.bam --xc-bin-size 50
"""

import subprocess
import sys
import re
import click
import json
import pysam
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import tempfile
import os

# Import VRS-compatible sha512t24u function + VRS v2 Allele identifier support
try:
    from vrs_compat import sha512t24u, VRSAlleleV2
except ImportError:
    from isotag_utils import sha512t24u_fallback as sha512t24u  # type: ignore[assignment]
    VRSAlleleV2 = None  # VRS v2 XV generation unavailable without vrs_compat.py

from isotag_utils import mask_ambiguous_bases  # shared with isotag_refget.py


@dataclass
class ExonBoundary:
    """Represents exon start and end coordinates"""
    start: int
    end: int

    def __lt__(self, other):
        return self.start < other.start


@dataclass
class SimpleVariant:
    """Simple variant representation: 0-based inter-residue position + literal ref/alt"""
    chromosome: str
    position: int
    ref: str
    alt: str

    def __lt__(self, other):
        return self.position < other.position


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
            if not data.get("metadata", {}).get("ambiguous_bases_masked", True):
                click.echo("⚠️  WARNING: This RefGet JSON lacks ambiguous base masking (generated with v2.2.0 or earlier)")
                click.echo("⚠️  Tags may be silently wrong for GRCh38 and other genomes with IUPAC codes")
                click.echo("⚠️  Regenerate: python3 isotag_refget.py -f genome.fa -o genome-refget.json -g genome_name")
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
            if not data.get("metadata", {}).get("ambiguous_bases_masked", True):
                click.echo(f"⚠️  Stale cache lacks ambiguous base masking — deleting and regenerating: {cache_path}")
                cache_path.unlink()
            else:
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
            click.echo("   ℹ️  This ensures UCSC hg38 and NCBI GRCh38 produce identical hashes")

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
    - XV tag: Full GA4GH VRS v2 Allele IDs (SNVs/substitutions, comma-separated)
    """

    def __init__(self,
                 refget_mapping: Optional[Dict[str, str]] = None,
                 clustermode: str = "middle",
                 position_quantum: int = 10000,
                 span_quantum: int = 1000,
                 exon_quantum: int = 1000,
                 xc_bin_size: int = 10,
                 xc_position_quantum: int = 1000,
                 x5x3_bin_size: int = 200,
                 xc_midpoint_only: bool = False):
        """
        Initialize IsoformTagger

        Args:
            refget_mapping: RefGet chromosome hash mapping
            clustermode: Position for XT clustering (5prime, middle, 3prime)
            position_quantum: Bin size for position quantization (bp)
            span_quantum: Bin size for genomic span quantization (bp)
            exon_quantum: Bin size for exon length quantization (bp)
            xc_bin_size: Bin size for XC exon length bins (bp, default: 10).
                         Absorbs ±1bp Nanopore splice site wobble.
            xc_position_quantum: Bin size for XC midpoint position (bp, default: 1000).
                                 Coarser than xc_bin_size for locus anchoring.
            x5x3_bin_size: Bin size for X5/X3 end position bins (bp, default: 200).
                           Groups reads by TSS (X5) or polyA site (X3) within this window.
        """
        self.refget_mapping = refget_mapping or {}
        self.clustermode = clustermode
        self.position_quantum = position_quantum
        self.span_quantum = span_quantum
        self.exon_quantum = exon_quantum
        self.xc_bin_size = xc_bin_size
        self.xc_position_quantum = xc_position_quantum
        self.x5x3_bin_size = x5x3_bin_size
        self.xc_midpoint_only = xc_midpoint_only

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
        Serialize transcript for XT tag with mode-based clustering.
        NOTE: splice_junctions are exact integer coordinates (not binned). Any ±1bp
        Nanopore wobble at any junction produces a different XT hash — this tag is
        NOT wobble-tolerant despite prior documentation claiming "fuzzy matching".
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

    def serialize_cluster_xc(self, chromosome: str, strand: str, exons: List[ExonBoundary]) -> str:
        """
        Serialize transcript for XC cluster tag - midpoint + binned exon lengths (v11.0).

        XC groups transcripts with similar exon structure at the same genomic locus,
        tolerating ±1bp Nanopore splice site wobble without requiring junction coordinates.

        Format: chr_hash_32|strand|middle_bin|exon_len_bin1|exon_len_bin2|...

        Two bin sizes serve different purposes:
          - xc_position_quantum (default 1kb): coarse locus anchor via transcript midpoint.
            Midpoint is more truncation-tolerant than either terminal end alone.
          - xc_bin_size (default 10bp): fine exon length bins absorb ±1bp wobble.
            Exon count is implicit in the number of length bins (no explicit count needed).

        Args:
            chromosome: Chromosome name
            strand: Strand (+ or -)
            exons: List of exon boundaries

        Returns:
            Serialized string for XC hashing (NO splice junction coordinates)
        """
        chr_hash = self.get_chromosome_hash(chromosome, hash_length=32)
        sorted_exons = sorted(exons)

        # Midpoint position — more truncation-tolerant than start or end alone
        genomic_start = sorted_exons[0].start
        genomic_end = sorted_exons[-1].end
        midpoint = (genomic_start + genomic_end) // 2
        middle_bin = midpoint // self.xc_position_quantum

        # Bin exon lengths — absorbs ±1bp Nanopore wobble at splice sites
        # Exon count is implicit in the number of bins (no explicit count needed)
        exon_length_bins = []
        for exon in sorted_exons:
            exon_length = exon.end - exon.start + 1
            binned_length = exon_length // self.xc_bin_size
            exon_length_bins.append(str(binned_length))

        # Build serialization — NO splice junction coordinates, NO explicit exon count
        # --xc-midpoint-only: skip exon lengths for midpoint-only baseline comparison
        parts = [
            chr_hash,
            strand,
            str(middle_bin)
        ]
        if not self.xc_midpoint_only:
            parts += exon_length_bins

        return '|'.join(parts)

    def generate_cluster_xc_id(self, chromosome: str, strand: str, exons: List[ExonBoundary]) -> str:
        """
        Generate XC cluster tag ID - midpoint + binned exon lengths (v11.0).

        XC captures:
        - Chromosome (via RefGet hash)
        - Strand
        - Transcript midpoint (binned by --xc-position-quantum, default 1kb)
        - Exon lengths (binned by --xc-bin-size, default 10bp)
          Exon count is implicit in the number of length bins.

        XC ignores (unlike XI/XS/XT):
        - Precise splice junction coordinates
        - Exact exon boundary positions
        - Explicit exon count (implicit via length bin count)

        This tolerates ±1bp Nanopore splice site wobble while distinguishing
        transcripts with different exon sizes or substantially different loci.
        Midpoint is more robust to terminal truncation than start or end alone.
        """
        xc_serial = self.serialize_cluster_xc(chromosome, strand, exons)
        return self.generate_hash(xc_serial)

    def generate_x5_id(self, chromosome: str, strand: str, exons: List[ExonBoundary]) -> str:
        """
        Generate X5 tag — binned 5' transcript end for TSS/CAGE clustering.

        X5 groups reads sharing the same transcription start site within
        self.x5x3_bin_size bp windows (default 200 bp).

        Convention (same as XB tag):
            + strand: 5' end = leftmost exon start
            - strand: 5' end = rightmost exon end

        Format hashed: chr_hash_32|strand|bin(5prime, x5x3_bin_size)
        """
        chr_hash = self.get_chromosome_hash(chromosome, hash_length=32)
        sorted_exons = sorted(exons)
        if strand == '+':
            five_prime = sorted_exons[0].start
        else:
            five_prime = sorted_exons[-1].end
        pos_bin = five_prime // self.x5x3_bin_size
        serial = f"{chr_hash}|{strand}|{pos_bin}"
        return self.generate_hash(serial)

    def generate_x3_id(self, chromosome: str, strand: str, exons: List[ExonBoundary]) -> str:
        """
        Generate X3 tag — binned 3' transcript end for polyA site clustering.

        X3 groups reads sharing the same cleavage/polyadenylation site within
        self.x5x3_bin_size bp windows (default 200 bp).

        Convention (same as XB tag):
            + strand: 3' end = rightmost exon end
            - strand: 3' end = leftmost exon start

        Format hashed: chr_hash_32|strand|bin(3prime, x5x3_bin_size)
        """
        chr_hash = self.get_chromosome_hash(chromosome, hash_length=32)
        sorted_exons = sorted(exons)
        if strand == '+':
            three_prime = sorted_exons[-1].end
        else:
            three_prime = sorted_exons[0].start
        pos_bin = three_prime // self.x5x3_bin_size
        serial = f"{chr_hash}|{strand}|{pos_bin}"
        return self.generate_hash(serial)

    def generate_individual_variant_ids(self, chromosome: str, variants: List[SimpleVariant]) -> Optional[str]:
        """
        Generate canonical GA4GH VRS v2 Allele IDs for reference-independent substitutions.

        Returns comma-separated `ga4gh:VA.<digest>` identifiers (sorted, deduplicated),
        or None if there are no emittable variants. Pure insertions/deletions are
        deliberately omitted: their canonical VRS representation requires reference-aware
        normalization not yet implemented here (see vrs_compat.normalize_substitution).

        Raises:
            RuntimeError: vrs_compat.VRSAlleleV2 is unavailable (sha512t24u-only fallback mode)
            ValueError: no RefGet mapping is available for `chromosome`
        """
        if not variants:
            return None
        if VRSAlleleV2 is None:
            raise RuntimeError("vrs_compat.py (with VRSAlleleV2) is required for VRS v2 XV generation")

        refget_accession = self.refget_mapping.get(chromosome)
        if not refget_accession:
            raise ValueError(f"VRS XV requires a RefGet mapping for reference {chromosome!r}")

        individual_ids = []
        for variant in sorted(variants):
            try:
                allele_id = VRSAlleleV2(
                    refget_accession=refget_accession,
                    start=variant.position,
                    ref=variant.ref,
                    alt=variant.alt,
                ).identifier()
            except ValueError:
                # Indel pending reference-aware normalization — omit, don't fail the read.
                continue
            individual_ids.append(allele_id)

        # A comma is unambiguous because GA4GH identifiers contain a period.
        return ",".join(sorted(set(individual_ids))) or None

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

    def parse_md_events(self, md_tag: str) -> Dict[int, Tuple[str, str]]:
        """Map the MD reference axis to mismatch/deletion events.

        The MD axis includes aligned reference bases and CIGAR deletions, but
        excludes CIGAR N skipped regions. Keeping it separate from the genomic
        cursor prevents coordinate drift after introns.
        """
        events: Dict[int, Tuple[str, str]] = {}
        md_offset = 0
        i = 0
        while i < len(md_tag):
            if md_tag[i].isdigit():
                j = i
                while j < len(md_tag) and md_tag[j].isdigit():
                    j += 1
                md_offset += int(md_tag[i:j])
                i = j
            elif md_tag[i] == "^":
                i += 1
                while i < len(md_tag) and md_tag[i].isalpha():
                    events[md_offset] = ("deletion", md_tag[i].upper())
                    md_offset += 1
                    i += 1
            elif md_tag[i].isalpha():
                events[md_offset] = ("mismatch", md_tag[i].upper())
                md_offset += 1
                i += 1
            else:
                i += 1
        return events

    def extract_variants_from_alignment(
        self,
        cigar_operations: List[Tuple[int, int]],
        sam_position: int,
        query_sequence: str,
        md_tag: str,
    ) -> List[SimpleVariant]:
        """Extract read-supported alleles while synchronizing CIGAR, MD, and SEQ.

        Walks the genomic cursor (advanced by every CIGAR op, including N introns)
        separately from the MD-axis cursor (advanced only by M/D, per parse_md_events),
        so mismatches downstream of a splice junction land at the correct genomic
        position. Positions are 0-based inter-residue, matching VRS convention.
        """
        if not md_tag or not query_sequence or query_sequence == "*":
            return []

        events = self.parse_md_events(md_tag)
        variants: List[SimpleVariant] = []
        reference_cursor = sam_position - 1
        query_cursor = 0
        md_cursor = 0

        for op_code, length in cigar_operations:
            if op_code in (0, 7, 8):  # M, =, X
                for _ in range(length):
                    event = events.get(md_cursor)
                    if event and event[0] == "mismatch":
                        ref = event[1]
                        alt = query_sequence[query_cursor].upper()
                        if ref in "ACGT" and alt in "ACGT" and ref != alt:
                            variants.append(
                                SimpleVariant("unknown", reference_cursor, ref, alt)
                            )
                    reference_cursor += 1
                    query_cursor += 1
                    md_cursor += 1
            elif op_code == 1:  # I
                inserted = query_sequence[query_cursor : query_cursor + length].upper()
                variants.append(
                    SimpleVariant("unknown", reference_cursor, "", inserted)
                )
                query_cursor += length
            elif op_code == 2:  # D
                deleted = "".join(
                    events.get(md_cursor + offset, ("deletion", "N"))[1]
                    for offset in range(length)
                )
                variants.append(
                    SimpleVariant("unknown", reference_cursor, deleted, "")
                )
                reference_cursor += length
                md_cursor += length
            elif op_code == 3:  # N: excluded from the MD coordinate axis
                reference_cursor += length
            elif op_code == 4:  # S
                query_cursor += length
            elif op_code in (5, 6):  # H, P
                pass

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

        # Generate XC tag ID (location-based cluster, no splice precision) - 32-char hash
        xc_cluster_id = self.generate_cluster_xc_id(rname, strand, exons)

        # Generate X5/X3 end-site tags (TSS and polyA site clustering)
        x5_id = self.generate_x5_id(rname, strand, exons)
        x3_id = self.generate_x3_id(rname, strand, exons)

        # Extract variants if enabled
        variant_id = None
        variant_candidate_count = 0
        variant_emitted_count = 0
        if detect_variants:
            variants = []

            # Look for MD tag
            md_tag = None
            for field in fields[11:]:
                if field.startswith('MD:Z:'):
                    md_tag = field[5:]
                    break

            if md_tag:
                seq = fields[9]
                extracted_variants = self.extract_variants_from_alignment(cigar_ops, pos, seq, md_tag)
                for variant in extracted_variants:
                    variant.chromosome = rname
                variants.extend(extracted_variants)

            # Generate full GA4GH VRS v2 identifiers.
            variant_candidate_count = len(variants)
            variant_id = self.generate_individual_variant_ids(rname, variants)
            variant_emitted_count = len(variant_id.split(",")) if variant_id else 0

        return {
            'qname': qname,
            'isoform_id': structure_id,
            'boundarytag_id': boundarytag_id,
            'splicetag_id': splicetag_id,
            'xt_group_id': xt_group_id,
            'xc_cluster_id': xc_cluster_id,
            'x5_id': x5_id,
            'x3_id': x3_id,
            'variant_id': variant_id,
            'variant_candidate_count': variant_candidate_count,
            'variant_emitted_count': variant_emitted_count,
            'original_line': line.strip()
        }

    def process_read(self, read: pysam.AlignedSegment, detect_variants: bool = True) -> Optional[Dict]:
        """Process a pysam AlignedSegment and generate universal isoform/variant IDs.

        Pysam-native equivalent of process_sam_line — no SAM text parsing.
        """
        if read.is_unmapped or not read.cigartuples:
            return None

        rname = read.reference_name
        pos = read.reference_start + 1  # pysam is 0-based; SAM convention is 1-based
        cigar_ops = list(read.cigartuples)  # same (op_code, length) format as parse_cigar_operations

        strand = '-' if read.is_reverse else '+'

        exons = self.extract_exons_from_cigar(pos, cigar_ops)
        if not exons:
            return None

        structure_id = self.generate_structure_id(rname, strand, exons)
        boundarytag_id = self.generate_boundarytag_id(rname, strand, exons)
        splicetag_id = self.generate_splicetag_id(rname, strand, exons)
        xt_group_id = self.generate_transcript_group_xt_id(rname, strand, exons)
        xc_cluster_id = self.generate_cluster_xc_id(rname, strand, exons)
        x5_id = self.generate_x5_id(rname, strand, exons)
        x3_id = self.generate_x3_id(rname, strand, exons)

        variant_id = None
        variant_candidate_count = 0
        variant_emitted_count = 0
        if detect_variants:
            variants = []

            if read.has_tag('MD'):
                md_tag = read.get_tag('MD')
                seq = read.query_sequence or ''
                extracted_variants = self.extract_variants_from_alignment(cigar_ops, pos, seq, md_tag)
                for variant in extracted_variants:
                    variant.chromosome = rname
                variants.extend(extracted_variants)

            variant_candidate_count = len(variants)
            variant_id = self.generate_individual_variant_ids(rname, variants)
            variant_emitted_count = len(variant_id.split(",")) if variant_id else 0

        return {
            'isoform_id': structure_id,
            'boundarytag_id': boundarytag_id,
            'splicetag_id': splicetag_id,
            'xt_group_id': xt_group_id,
            'xc_cluster_id': xc_cluster_id,
            'x5_id': x5_id,
            'x3_id': x3_id,
            'variant_id': variant_id,
            'variant_candidate_count': variant_candidate_count,
            'variant_emitted_count': variant_emitted_count,
        }

    def check_for_md_tags(self, bam_file: str, sample_size: int = 100) -> bool:
        """Check if BAM file contains MD tags by sampling reads"""
        try:
            mode = 'rb' if Path(bam_file).suffix.lower() == '.bam' else 'r'
            with pysam.AlignmentFile(bam_file, mode) as bam:
                reads_checked = 0
                md_tags_found = 0
                for read in bam:
                    if read.is_unmapped:
                        continue
                    reads_checked += 1
                    if reads_checked > sample_size:
                        break
                    if read.has_tag('MD'):
                        md_tags_found += 1

            return md_tags_found > 0 and (md_tags_found / reads_checked) >= 0.1 if reads_checked > 0 else False

        except Exception:
            return False


@click.command()
@click.option('--input', '-i', 'input_file', required=True, help='Input BAM/SAM file')
@click.option('--output', '-o', required=True, help='Output BAM/SAM file with XI/XB/XS/XT/XV/XC/X5/X3 tags')
@click.option('--genome', '-g', help='Reference genome FASTA (for RefGet cache generation and/or variant detection)')
@click.option('--refget', '-r', help='RefGet JSON mapping file (optional, will auto-generate from genome if not provided)')
@click.option('--no-variants', is_flag=True,
              help='Disable XV variant tag detection even when genome is provided or MD tags are present')
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
@click.option('--xc-bin-size', type=int, default=10,
              help='Bin size for XC exon length bins in bp (default: 10). Absorbs ±1bp Nanopore wobble.')
@click.option('--xc-position-quantum', type=int, default=1000,
              help='Bin size for XC midpoint position in bp (default: 1000). Coarser locus anchor.')
@click.option('--x5x3-bin-size', type=int, default=200,
              help='Bin size for X5/X3 end-site tags in bp (default: 200). '
                   'X5 clusters reads by TSS; X3 clusters reads by polyA site.')
@click.option('--xc-midpoint-only', is_flag=True,
              help='XC uses midpoint position only (no exon lengths) — for baseline comparison vs v11.0')
def isotag(input_file, output, genome, refget, no_variants, keep_ambiguous_bases, clustermode, position_quantum, span_quantum, exon_quantum, xc_bin_size, xc_position_quantum, x5x3_bin_size, xc_midpoint_only):
    """
    IsoTag - Universal Isoform Tagger

    Add RefGet-based universal tags to BAM/SAM files:
    - XI: Isoform structure (32-char hash with full RefGet chr hash)
    - XB: Reversible boundary tag (8-char chr hash + hex 5'/3' ends)
    - XS: Reversible splicetag (8-char chr hash + hex coordinates)
    - XT: Transcript group (32-char hash, mode-based clustering)
    - XC: Cluster ID (32-char hash, midpoint + binned exon lengths, wobble-tolerant)
    - XV: Full GA4GH VRS v2 Allele IDs (SNVs/substitutions, comma-separated)

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

    VRS v2 XV Safety Boundary:
    - Requires --refget or --genome (chromosome name hashing has no RefGet accession for VRS)
    - Rejects a RefGet mapping generated with ambiguous bases masked (default masking behavior
      of this tool and of any refget-vX.json bundled with IsoTag) — VRS Allele IDs must be
      computed from the exact canonical reference sequence to match ClinVar/ClinGen. Regenerate
      with --keep-ambiguous-bases from the exact canonical FASTA if you need real XV values.
    - Pure insertions/deletions are omitted from XV (reference-aware VRS normalization not
      yet implemented); XI/XB/XS/XT/XC/X5/X3 are unaffected.

    Examples:
        # Basic tagging with auto RefGet cache
        isotag.py -i input.bam -o tagged.bam -g genome.fa

        # CAGE data (5' clustering)
        isotag.py -i cage.bam -o tagged.bam -g genome.fa --clustermode 5prime

        # Use pre-computed RefGet mapping
        isotag.py -i input.bam -o tagged.bam -r genome-refget.json

        # Custom quantization
        isotag.py -i input.bam -o tagged.bam -g genome.fa --position-quantum 5000

        # Gene-level XC clustering with larger bins (50kb instead of 10kb default)
        isotag.py -i input.bam -o tagged.bam -g genome.fa --xc-bin-size 50000
    """

    input_path = Path(input_file)
    output_path = Path(output)

    if not input_path.exists():
        click.echo(f"❌ Input file not found: {input_path}")
        sys.exit(1)

    # Detect file formats
    input_is_bam = input_path.suffix.lower() == '.bam'
    output_is_bam = output_path.suffix.lower() == '.bam'

    click.echo("🧬 IsoTag v11.1 - XC Cluster Tag (midpoint + binned exon lengths) + GA4GH VRS v2 XV")
    click.echo(f"📥 Input: {input_path.name} ({'BAM' if input_is_bam else 'SAM'})")
    click.echo(f"📤 Output: {output_path.name} ({'BAM' if output_is_bam else 'SAM'})")
    click.echo(f"🎯 Cluster mode: {clustermode}")
    click.echo(f"📏 XT quantization: position={position_quantum}bp, span={span_quantum}bp, exon={exon_quantum}bp")
    click.echo(f"🔗 XC: midpoint/{xc_position_quantum}bp bins, exon lengths/{xc_bin_size}bp bins")

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
        exon_quantum=exon_quantum,
        xc_bin_size=xc_bin_size,
        xc_position_quantum=xc_position_quantum,
        x5x3_bin_size=x5x3_bin_size,
        xc_midpoint_only=xc_midpoint_only,
    )

    # Statistics
    stats = {
        'total_reads': 0,
        'reads_processed': 0,
        'reads_with_structure': 0,
        'reads_with_boundarytags': 0,
        'reads_with_splicetags': 0,
        'reads_single_exon_xs': 0,
        'reads_with_xt_groups': 0,
        'reads_with_xc_clusters': 0,
        'reads_with_x5_tags': 0,
        'reads_with_x3_tags': 0,
        'reads_with_variants': 0,
        'variant_candidates': 0,
        'variant_ids_emitted': 0,
        'unique_structures': {'count': 0, 'example': None},
        'unique_boundarytags': {'count': 0, 'example': None},
        'unique_splicetags': {'count': 0, 'example': None},
        'unique_xt_groups': {'count': 0, 'example': None},
        'unique_xc_clusters': {'count': 0, 'example': None},
        'unique_x5_tags': {'count': 0, 'example': None},
        'unique_x3_tags': {'count': 0, 'example': None},
        'unique_variants': {'count': 0, 'example': None}
    }

    temp_bam_path = None

    try:
        # Determine variant detection strategy
        detect_variants = False
        input_for_processing = str(input_path)
        needs_cleanup = False

        if no_variants:
            click.echo("🔬 Variant detection: Disabled (--no-variants)")
        elif genome:
            # Generate MD tags for variant detection
            click.echo("📋 Preparing BAM with MD tags for variant detection...")
            with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as temp_bam:
                temp_bam_path = temp_bam.name

            with open(temp_bam_path, 'wb') as calmd_out:
                subprocess.run(['samtools', 'calmd', '-b', str(input_path), genome],
                             stdout=calmd_out, check=True)
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

        if detect_variants and not refget_mapping:
            raise click.ClickException(
                "VRS v2 XV generation requires --refget or --genome "
                "(legacy chromosome-name hashing has no RefGet accession for VRS)"
            )
        if detect_variants:
            # Ambiguous-base masking (this tool's default, and the default for any bundled
            # refget/*.json) changes the chromosome sequence hash and therefore produces
            # XV identifiers that will not match ClinVar/ClinGen VRS-indexed data.
            mapping_is_masked = True
            if refget:
                with open(refget, "r") as refget_handle:
                    refget_document = json.load(refget_handle)
                mapping_is_masked = bool(
                    refget_document.get("metadata", {}).get("ambiguous_bases_masked", True)
                )
            elif genome:
                mapping_is_masked = not keep_ambiguous_bases
            if mapping_is_masked:
                raise click.ClickException(
                    "This RefGet mapping was generated after masking ambiguous bases and cannot "
                    "be used for external VRS matching. Regenerate it from the exact canonical "
                    "FASTA with --keep-ambiguous-bases (XI/XB/XS/XT/XC/X5/X3 are unaffected by "
                    "this check and may still use a masked mapping)."
                )

        click.echo("🏷️  Adding isoform tags (XI, XB, XS, XT, XC, X5, X3)...")
        if detect_variants:
            click.echo("🏷️  Adding variant tags (XV)...")

        # Process BAM/SAM directly with pysam — no temp SAM file, no samtools subprocess
        mode_r = 'rb' if input_is_bam else 'r'
        mode_w = 'wb' if output_is_bam else 'w'

        with pysam.AlignmentFile(input_for_processing, mode_r) as in_bam:
            header_dict = in_bam.header.to_dict()
            header_dict.setdefault('PG', []).append({
                'ID': 'isotag', 'PN': 'isotag', 'VN': '2.3.0',
                'CL': ' '.join(sys.argv),
            })
            out_header = pysam.AlignmentHeader.from_dict(header_dict)
        with pysam.AlignmentFile(input_for_processing, mode_r) as in_bam, \
             pysam.AlignmentFile(str(output_path), mode_w, header=out_header) as out_bam:

            for read in in_bam:
                stats['total_reads'] += 1

                if stats['total_reads'] % 10000 == 0:
                    click.echo(f"⏳ Processed {stats['total_reads']:,} reads...")

                result = tagger.process_read(read, detect_variants)

                if result:
                    stats['reads_processed'] += 1

                    # Set tags directly on the pysam read object
                    read.set_tag('XI', result['isoform_id'])
                    stats['reads_with_structure'] += 1
                    u = stats['unique_structures']
                    u['count'] += 1
                    if u['example'] is None:
                        u['example'] = result['isoform_id']

                    read.set_tag('XB', result['boundarytag_id'])
                    stats['reads_with_boundarytags'] += 1
                    u = stats['unique_boundarytags']
                    u['count'] += 1
                    if u['example'] is None:
                        u['example'] = result['boundarytag_id']

                    if result['splicetag_id']:
                        read.set_tag('XS', result['splicetag_id'])
                        stats['reads_with_splicetags'] += 1
                        u = stats['unique_splicetags']
                        u['count'] += 1
                        if u['example'] is None:
                            u['example'] = result['splicetag_id']
                    else:
                        # Sentinel for single-exon reads (no splice junctions).
                        # Prevents silent GROUP BY XS grouping of all single-exon reads as NULL.
                        read.set_tag('XS', 'single')
                        stats['reads_single_exon_xs'] += 1

                    read.set_tag('XT', result['xt_group_id'])
                    stats['reads_with_xt_groups'] += 1
                    u = stats['unique_xt_groups']
                    u['count'] += 1
                    if u['example'] is None:
                        u['example'] = result['xt_group_id']

                    read.set_tag('XC', result['xc_cluster_id'])
                    stats['reads_with_xc_clusters'] += 1
                    u = stats['unique_xc_clusters']
                    u['count'] += 1
                    if u['example'] is None:
                        u['example'] = result['xc_cluster_id']

                    read.set_tag('X5', result['x5_id'])
                    stats['reads_with_x5_tags'] += 1
                    u = stats['unique_x5_tags']
                    u['count'] += 1
                    if u['example'] is None:
                        u['example'] = result['x5_id']

                    read.set_tag('X3', result['x3_id'])
                    stats['reads_with_x3_tags'] += 1
                    u = stats['unique_x3_tags']
                    u['count'] += 1
                    if u['example'] is None:
                        u['example'] = result['x3_id']

                    if detect_variants:
                        stats['variant_candidates'] += result['variant_candidate_count']
                        stats['variant_ids_emitted'] += result['variant_emitted_count']
                        if result['variant_id']:
                            read.set_tag('XV', result['variant_id'])
                            stats['reads_with_variants'] += 1
                            u = stats['unique_variants']
                            u['count'] += 1
                            if u['example'] is None:
                                u['example'] = result['variant_id']

                out_bam.write(read)

        # Cleanup
        if needs_cleanup and temp_bam_path is not None:
            os.unlink(temp_bam_path)
            temp_bam_path = None

        # Display results
        click.echo("\n" + "="*60)
        click.echo("✅ IsoTag Complete!")
        click.echo("="*60)
        click.echo(f"📊 Total reads: {stats['total_reads']:,}")
        click.echo(f"🧬 Reads processed: {stats['reads_processed']:,}")
        click.echo(f"🏷️  XI tags (structure): {stats['reads_with_structure']:,}")
        click.echo(f"🔚 XB tags (boundaries): {stats['reads_with_boundarytags']:,}")
        click.echo(f"🔗 XS tags (splicetag): {stats['reads_with_splicetags']:,} multi-exon, {stats['reads_single_exon_xs']:,} single-exon (sentinel 'single')")
        click.echo(f"🎯 XT tags (transcript group): {stats['reads_with_xt_groups']:,}")
        click.echo(f"📍 XC tags (cluster): {stats['reads_with_xc_clusters']:,}")
        click.echo(f"⬆️  X5 tags (TSS/5' end): {stats['reads_with_x5_tags']:,}")
        click.echo(f"⬇️  X3 tags (polyA/3' end): {stats['reads_with_x3_tags']:,}")
        click.echo(f"🔬 XV tags (variants): {stats['reads_with_variants']:,}")
        if detect_variants:
            skipped = stats['variant_candidates'] - stats['variant_ids_emitted']
            click.echo(f"🧬 Variant candidates (MD mismatches/indels): {stats['variant_candidates']:,}")
            click.echo(f"✅ VRS v2 IDs emitted: {stats['variant_ids_emitted']:,}")
            click.echo(f"⏸️  Indel/unsupported candidates skipped: {skipped:,}")
        click.echo(f"🆔 Unique structures: {stats['unique_structures']['count']:,}")
        click.echo(f"🔚 Unique boundarytags: {stats['unique_boundarytags']['count']:,}")
        click.echo(f"🧬 Unique splicetags: {stats['unique_splicetags']['count']:,}")
        click.echo(f"🎯 Unique XT groups: {stats['unique_xt_groups']['count']:,}")
        click.echo(f"📍 Unique XC clusters: {stats['unique_xc_clusters']['count']:,}")
        click.echo(f"⬆️  Unique X5 TSS sites: {stats['unique_x5_tags']['count']:,}")
        click.echo(f"⬇️  Unique X3 polyA sites: {stats['unique_x3_tags']['count']:,}")
        click.echo(f"🧪 Unique variant combos: {stats['unique_variants']['count']:,}")
        click.echo(f"💾 Output: {output_path}")

        # Show example tags
        click.echo("\n🎯 Example tags added:")
        if stats['unique_structures']['example']:
            click.echo(f"   XI:Z:{stats['unique_structures']['example']} (32-char structure)")
        if stats['unique_boundarytags']['example']:
            click.echo(f"   XB:Z:{stats['unique_boundarytags']['example']} (reversible boundary tag)")
        if stats['unique_splicetags']['example']:
            click.echo(f"   XS:Z:{stats['unique_splicetags']['example'][:50]}... (reversible splicetag)")
        if stats['unique_xt_groups']['example']:
            click.echo(f"   XT:Z:{stats['unique_xt_groups']['example']} (32-char {clustermode} group)")
        if stats['unique_xc_clusters']['example']:
            click.echo(f"   XC:Z:{stats['unique_xc_clusters']['example']} (32-char cluster, midpoint/{xc_position_quantum}bp + exon_lengths/{xc_bin_size}bp)")
        if stats['unique_variants']['example']:
            click.echo(f"   XV:Z:{stats['unique_variants']['example'][:50]}... (variant IDs)")

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
    finally:
        if temp_bam_path and os.path.exists(temp_bam_path):
            os.unlink(temp_bam_path)


if __name__ == '__main__':
    isotag()
