#!/usr/bin/env python3
"""
Splicetag Decoder - Decode reversible splicetag IDs to coordinates

Decodes reversible splicetag format back to splice junction coordinates.
Format: [8-CHAR-CHR-HASH][STRAND].hex1.hex2.hex3...

Usage:
    python3 decode_splicetag.py -s "aKF498dAp.4b0.7d0.866.bb8"
    python3 decode_splicetag.py -s "aKF498dAp.4b0.7d0.866.bb8" -r genome-refget.json
"""

import click
import json
from typing import List, Tuple, Optional, Dict


def decode_reversible_splicetag(encoded: str) -> Tuple[str, str, List[int]]:
    """
    Decode reversible splicetag back to components

    Args:
        encoded: Reversible splicetag string

    Returns:
        (chr_hash_8, strand, splice_coordinates)
    """
    if len(encoded) < 10:
        raise ValueError(f"Invalid splicetag format: {encoded}")

    chr_hash_8 = encoded[:8]
    strand = '+' if encoded[8] == 'p' else '-'

    if encoded[9] != '.':
        raise ValueError(f"Invalid splicetag format: expected '.' at position 9")

    hex_coords = encoded[10:]  # Skip chr_hash + strand + first dot

    if not hex_coords:
        raise ValueError(f"No coordinates found in splicetag: {encoded}")

    try:
        coordinates = [int(coord, 16) for coord in hex_coords.split('.')]
    except ValueError as e:
        raise ValueError(f"Invalid hex coordinate in splicetag: {e}")

    return chr_hash_8, strand, coordinates


def lookup_chromosome(chr_hash_8: str, refget_file: Optional[str] = None) -> List[str]:
    """
    Look up chromosome names from 8-character hash

    Args:
        chr_hash_8: 8-character chromosome hash
        refget_file: RefGet JSON mapping file

    Returns:
        List of matching chromosome names
    """
    if not refget_file:
        return [f"chr? (hash: {chr_hash_8})"]

    try:
        with open(refget_file, 'r') as f:
            data = json.load(f)
            refget_mapping = data.get("refget_mapping", {})

        matches = []
        for chrom, refget_id in refget_mapping.items():
            if refget_id.startswith(f"SQ.{chr_hash_8}"):
                matches.append(chrom)

        return matches if matches else [f"Unknown (hash: {chr_hash_8})"]

    except (FileNotFoundError, json.JSONDecodeError) as e:
        return [f"Error reading RefGet file: {e}"]


def format_splice_junctions(coordinates: List[int]) -> List[Tuple[int, int]]:
    """
    Format splice junction coordinates as (exon_end, next_exon_start) pairs

    Args:
        coordinates: List of splice junction coordinates

    Returns:
        List of (exon_end, next_exon_start) tuples
    """
    if len(coordinates) % 2 != 0:
        raise ValueError(f"Invalid number of coordinates: {len(coordinates)} (must be even)")

    junctions = []
    for i in range(0, len(coordinates), 2):
        junctions.append((coordinates[i], coordinates[i + 1]))

    return junctions


@click.command()
@click.option('--splicetag', '-s', required=True, help='Reversible splicetag ID to decode')
@click.option('--refget', '-r', help='RefGet JSON mapping file (optional, for chromosome name lookup)')
def decode(splicetag, refget):
    """
    Splicetag Decoder - Decode reversible splicetag IDs

    Decodes reversible splicetag format back to chromosome hash, strand,
    and splice junction coordinates.

    Examples:
        # Basic decoding
        python3 decode_splicetag.py -s "aKF498dAp.4b0.7d0.866.bb8"

        # With chromosome name lookup
        python3 decode_splicetag.py -s "aKF498dAp.4b0.7d0.866.bb8" -r genome-refget.json
    """
    click.echo(f"🔍 Decoding splicetag: {splicetag}")
    click.echo("="*60)

    try:
        # Decode splicetag
        chr_hash_8, strand, coordinates = decode_reversible_splicetag(splicetag)

        click.echo(f"\n📊 Decoded Components:")
        click.echo(f"   Chromosome hash (8-char): {chr_hash_8}")
        click.echo(f"   Strand: {strand}")
        click.echo(f"   Number of coordinates: {len(coordinates)}")
        click.echo(f"   Number of junctions: {len(coordinates) // 2}")

        # Look up chromosome names if RefGet provided
        if refget:
            chromosomes = lookup_chromosome(chr_hash_8, refget)
            click.echo(f"\n🧬 Matching Chromosomes:")
            for chrom in chromosomes:
                click.echo(f"   {chrom}")
        else:
            click.echo(f"\n⚠️  Chromosome name: Unknown (provide --refget for lookup)")
            click.echo(f"   Chromosome hash: {chr_hash_8}")

        # Format and display splice junctions
        junctions = format_splice_junctions(coordinates)

        click.echo(f"\n🔗 Splice Junctions:")
        for i, (exon_end, next_exon_start) in enumerate(junctions, 1):
            intron_length = next_exon_start - exon_end - 1
            click.echo(f"   Junction {i}: {exon_end} -> {next_exon_start} (intron: {intron_length:,} bp)")

        # Display raw coordinates
        click.echo(f"\n📍 Raw Coordinates (for IGV/BED):")
        click.echo(f"   {', '.join(map(str, coordinates))}")

        # Generate summary
        genomic_start = coordinates[0]
        genomic_end = coordinates[-1]
        genomic_span = genomic_end - genomic_start

        click.echo(f"\n📐 Transcript Summary:")
        click.echo(f"   Genomic start: {genomic_start:,}")
        click.echo(f"   Genomic end: {genomic_end:,}")
        click.echo(f"   Genomic span: {genomic_span:,} bp")
        click.echo(f"   Exon count: {len(junctions) + 1}")

    except ValueError as e:
        click.echo(f"❌ Error: {e}")
        return 1

    click.echo("\n✅ Decoding complete!")


if __name__ == '__main__':
    decode()
