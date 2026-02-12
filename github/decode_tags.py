#!/usr/bin/env python3
"""
IsoTag Tag Decoder - Decode XB and XS reversible tags

Decodes:
- XB:Z: Reversible boundary tag (8-char chr hash + hex 5'/3' ends)
- XS:Z: Reversible splicetag (8-char chr hash + hex coordinates)

Usage:
    python3 decode_tags.py -b "aKF498dAp.3e8.1004"
    python3 decode_tags.py -s "aKF498dAp.4b0.7d0.866.bb8"
    python3 decode_tags.py -b "aKF498dAp.3e8.1004" -s "aKF498dAp.4b0.7d0.866.bb8"
    python3 decode_tags.py -b "aKF498dAp.3e8.1004" -s "aKF498dAp.4b0.7d0" -r hg38-refget.json
"""

import click
import json
from typing import Dict, List, Tuple, Optional


def decode_boundarytag(encoded: str) -> Tuple[str, str, int, int]:
    """
    Decode XB reversible boundary tag

    Args:
        encoded: XB tag value (e.g., "aKF498dAp.3e8.1004")

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


def decode_splicetag(encoded: str) -> Tuple[str, str, List[int]]:
    """
    Decode XS reversible splicetag

    Args:
        encoded: XS tag value (e.g., "aKF498dAp.4b0.7d0.866.bb8")

    Returns:
        (chr_hash_8, strand, splice_coordinates)
    """
    chr_hash_8 = encoded[:8]
    strand = '+' if encoded[8] == 'p' else '-'
    hex_coords = encoded[10:]  # Skip chr_hash + strand + first dot

    coordinates = [int(coord, 16) for coord in hex_coords.split('.')]
    return chr_hash_8, strand, coordinates


def reconstruct_exons(boundarytag: str, splicetag: Optional[str] = None) -> List[Tuple[int, int]]:
    """
    Reconstruct full exon coordinates from XB and XS tags

    Args:
        boundarytag: XB tag value
        splicetag: XS tag value (optional for single-exon)

    Returns:
        List of (start, end) exon coordinates
    """
    chr_hash_8, strand, five_prime, three_prime = decode_boundarytag(boundarytag)

    # Single-exon transcript
    if not splicetag or splicetag == "None":
        return [(five_prime, three_prime)]

    # Multi-exon transcript
    _, _, splice_coords = decode_splicetag(splicetag)

    # Reconstruct exons from boundaries and splice junctions
    exons = []

    # First exon: 5' end to first junction
    exons.append((five_prime, splice_coords[0]))

    # Internal exons: junction pairs
    for i in range(1, len(splice_coords) - 1, 2):
        exons.append((splice_coords[i], splice_coords[i + 1]))

    # Last exon: last junction to 3' end
    exons.append((splice_coords[-1], three_prime))

    return exons


def lookup_chromosome_names(chr_hash_8: str, refget_file: Optional[str] = None) -> List[str]:
    """
    Look up chromosome names from 8-char hash using RefGet mapping

    Args:
        chr_hash_8: 8-character chromosome hash
        refget_file: Path to RefGet JSON mapping file

    Returns:
        List of matching chromosome names
    """
    if not refget_file:
        return [f"<chr_hash:{chr_hash_8}>"]

    with open(refget_file, 'r') as f:
        data = json.load(f)
        mapping = data.get("refget_mapping", {})

    matches = []
    for chrom, refget_id in mapping.items():
        if refget_id.startswith(f"SQ.{chr_hash_8}"):
            matches.append(chrom)

    return matches if matches else [f"<chr_hash:{chr_hash_8}>"]


@click.command()
@click.option('--boundarytag', '-b', help='XB tag value to decode')
@click.option('--splicetag', '-s', help='XS tag value to decode')
@click.option('--refget', '-r', help='RefGet JSON mapping file for chromosome name lookup')
@click.option('--reconstruct', is_flag=True, help='Reconstruct full exon coordinates from XB + XS')
def decode_tags(boundarytag, splicetag, refget, reconstruct):
    """
    Decode IsoTag reversible tags (XB and XS)

    Examples:
        # Decode XB boundary tag
        decode_tags.py -b "aKF498dAp.3e8.1004"

        # Decode XS splicetag
        decode_tags.py -s "aKF498dAp.4b0.7d0.866.bb8"

        # Reconstruct full exons from both tags
        decode_tags.py -b "aKF498dAp.3e8.1004" -s "aKF498dAp.4b0.7d0.866.bb8" --reconstruct

        # With chromosome name lookup
        decode_tags.py -b "aKF498dAp.3e8.1004" -r hg38-refget.json
    """

    if not boundarytag and not splicetag:
        click.echo("❌ Error: Provide at least one tag to decode (-b or -s)")
        return

    click.echo("🔍 IsoTag Tag Decoder v8.2\n")

    # Decode XB boundary tag
    if boundarytag:
        chr_hash_8, strand, five_prime, three_prime = decode_boundarytag(boundarytag)

        click.echo("📍 XB Boundary Tag Decoded:")
        click.echo(f"   Chromosome hash: {chr_hash_8}")
        click.echo(f"   Strand: {strand}")
        click.echo(f"   5' end: {five_prime:,} (0x{five_prime:x})")
        click.echo(f"   3' end: {three_prime:,} (0x{three_prime:x})")
        click.echo(f"   Genomic span: {three_prime - five_prime + 1:,} bp")

        # Chromosome name lookup
        if refget:
            chrom_names = lookup_chromosome_names(chr_hash_8, refget)
            click.echo(f"   Chromosomes: {', '.join(chrom_names)}")

        click.echo()

    # Decode XS splicetag
    if splicetag and splicetag != "None":
        chr_hash_8, strand, splice_coords = decode_splicetag(splicetag)

        click.echo("🔗 XS Splicetag Decoded:")
        click.echo(f"   Chromosome hash: {chr_hash_8}")
        click.echo(f"   Strand: {strand}")
        click.echo(f"   Splice junctions: {len(splice_coords) // 2}")
        click.echo(f"   Coordinates:")

        for i, coord in enumerate(splice_coords):
            junction_type = "exon end" if i % 2 == 0 else "exon start"
            click.echo(f"      [{i}] {coord:,} (0x{coord:x}) - {junction_type}")

        # Chromosome name lookup
        if refget:
            chrom_names = lookup_chromosome_names(chr_hash_8, refget)
            click.echo(f"   Chromosomes: {', '.join(chrom_names)}")

        click.echo()

    # Reconstruct full exons
    if reconstruct and boundarytag:
        click.echo("🧬 Reconstructed Exon Structure:")

        exons = reconstruct_exons(boundarytag, splicetag)

        chr_hash_8, strand, _, _ = decode_boundarytag(boundarytag)
        chrom_display = chr_hash_8

        if refget:
            chrom_names = lookup_chromosome_names(chr_hash_8, refget)
            if chrom_names and not chrom_names[0].startswith('<'):
                chrom_display = chrom_names[0]

        click.echo(f"   Chromosome: {chrom_display}")
        click.echo(f"   Strand: {strand}")
        click.echo(f"   Exons: {len(exons)}")
        click.echo(f"   Structure:")

        for i, (start, end) in enumerate(exons, 1):
            exon_length = end - start + 1
            click.echo(f"      Exon {i}: {start:,}-{end:,} ({exon_length:,} bp)")

        total_exon_length = sum(end - start + 1 for start, end in exons)
        genomic_span = exons[-1][1] - exons[0][0] + 1

        click.echo(f"\n   Total exon length: {total_exon_length:,} bp")
        click.echo(f"   Genomic span: {genomic_span:,} bp")

        if len(exons) > 1:
            total_intron_length = genomic_span - total_exon_length
            click.echo(f"   Total intron length: {total_intron_length:,} bp")


if __name__ == '__main__':
    decode_tags()
