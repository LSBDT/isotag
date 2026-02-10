# IsoTag Tag Format Specification v2.2

This document provides the complete technical specification for all IsoTag BAM/SAM tags.

## Table of Contents
- [Overview](#overview)
- [Tag Summary](#tag-summary)
- [XI Tag - Structure ID](#xi-tag---structure-id)
- [XB Tag - Boundary Tag](#xb-tag---boundary-tag)
- [XS Tag - Splicetag](#xs-tag---splicetag)
- [XT Tag - Transcript Group](#xt-tag---transcript-group)
- [XC Tag - Gene/Locus ID](#xc-tag---genelocus-id)
- [XV Tag - Variants](#xv-tag---variants)
- [RefGet Chromosome Hashing](#refget-chromosome-hashing)
- [Encoding Algorithms](#encoding-algorithms)
- [Decoding Algorithms](#decoding-algorithms)

---

## Overview

IsoTag v2.2 uses **six custom BAM tags** to encode transcript isoform structure, splice junctions, biological clustering, gene/locus identity, and variants:

- **XI**: Unique isoform structure identifier (32-char hash)
- **XB**: Reversible boundary tag (8-char chr hash + hex-encoded ends)
- **XS**: Reversible splicetag (8-char chr hash + hex-encoded splice junctions)
- **XT**: Biological transcript group (32-char hash with mode-based clustering)
- **XC**: Gene/locus cluster identifier (32-char hash, pure location-based)
- **XV**: Individual variant identifiers (32-char hashes, optional)

### Design Principles

1. **RefGet-compatible**: Uses GA4GH RefGet chromosome identifiers
2. **VRS-compliant**: Official `sha512t24u` algorithm throughout
3. **Reversible encoding**: XB and XS tags can be decoded to exact coordinates
4. **Universal compatibility**: Same chromosome sequence = same hash regardless of naming
5. **Compact storage**: Hex encoding minimizes BAM file size overhead

---

## Tag Summary

| Tag | Type | Length | Reversible | Purpose |
|-----|------|--------|------------|---------|
| **XI** | Z (string) | 32 chars | No | Unique isoform structure ID |
| **XB** | Z (string) | Variable | **Yes** | 5'/3' transcript boundary coordinates |
| **XS** | Z (string) | Variable | **Yes** | Internal splice junction coordinates |
| **XT** | Z (string) | 32 chars | No | Transcript group (fuzzy clustering) |
| **XC** | Z (string) | 32 chars | No | Gene/locus cluster (pure location) |
| **XV** | Z (string) | Variable | No | Individual variant IDs |

### Tag Presence Rules

- **XI**: Always present (required for all reads with CIGAR)
- **XB**: Always present (required for all reads with CIGAR)
- **XS**: Present only for multi-exon transcripts (2+ exons)
- **XT**: Always present (required for all reads with CIGAR)
- **XC**: Always present (required for all reads with CIGAR)
- **XV**: Present only when variants detected AND variant detection enabled

---

## XI Tag - Structure ID

### Purpose
Unique identifier for complete isoform structure (all exon coordinates).

### Format
```
XI:Z:[32-character-hash]
```

### Serialization (before hashing)
```
[32-chr-refget-hash]|[strand]|[exon1_start]:[exon1_end]|[exon2_start]:[exon2_end]|...
```

### Example
```bash
# Serialization
"aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2|+|1000:1200|2000:2150|3000:3500"

# Final XI tag
XI:Z:fuIF7PN23g2gq9sFxqhUNGnfOCZhkQJS
```

### Encoding Algorithm
```python
def generate_structure_id(chr_refget_32: str, strand: str, exons: List[Tuple[int, int]]) -> str:
    """Generate XI tag structure ID"""
    exon_strs = [f"{start}:{end}" for start, end in sorted(exons)]
    serialization = f"{chr_refget_32}|{strand}|{'|'.join(exon_strs)}"
    hash_bytes = serialization.encode('utf-8')
    return sha512t24u(hash_bytes)  # Returns 32-char base64url string
```

### Properties
- **Length**: Always 32 characters
- **Character set**: Base64URL (A-Z, a-z, 0-9, -, _)
- **Collision probability**: < 10^-50 for millions of structures
- **Reversibility**: No (one-way hash)

---

## XB Tag - Boundary Tag

### Purpose
Reversible encoding of transcript 5' and 3' ends for full coordinate reconstruction.

### Format
```
XB:Z:[8-chr-hash][strand].[5'-end-hex].[3'-end-hex]
```

### Components
- **8-chr-hash**: First 8 characters of RefGet chromosome hash
- **strand**: `p` (plus) or `m` (minus)
- **5'-end-hex**: 5' transcript end in hexadecimal (lowercase)
- **3'-end-hex**: 3' transcript end in hexadecimal (lowercase)

### Example
```bash
# Input
Chromosome: chr1 (RefGet hash: aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2)
Strand: +
5' end: 1000 (0x3e8)
3' end: 4100 (0x1004)

# XB tag
XB:Z:aKF498dAp.3e8.1004
```

### Visual Breakdown
```
aKF498dAp.3e8.1004
│       │ │   │
│       │ │   └─── 3' end: 4100 (hex: 1004)
│       │ └─────── 5' end: 1000 (hex: 3e8)
│       └───────── Strand: + (p=plus, m=minus)
└─────────────── Chromosome: aKF498dA (8-char RefGet hash)
```

### Encoding Algorithm
```python
def generate_boundarytag(chr_refget_32: str, strand: str,
                         five_prime_end: int, three_prime_end: int) -> str:
    """Generate XB tag boundary tag"""
    chr_hash_8 = chr_refget_32[:8]
    strand_char = 'p' if strand == '+' else 'm'
    five_hex = format(five_prime_end, 'x')  # Lowercase hex
    three_hex = format(three_prime_end, 'x')
    return f"{chr_hash_8}{strand_char}.{five_hex}.{three_hex}"
```

### Decoding Algorithm
```python
def decode_boundarytag(xb_tag: str) -> Tuple[str, str, int, int]:
    """Decode XB tag to coordinates"""
    chr_hash_8 = xb_tag[:8]
    strand = '+' if xb_tag[8] == 'p' else '-'
    coords = xb_tag[10:].split('.')  # Skip 8-char hash + strand + dot
    five_prime_end = int(coords[0], 16)
    three_prime_end = int(coords[1], 16)
    return chr_hash_8, strand, five_prime_end, three_prime_end
```

### Properties
- **Length**: Variable (typical: 20-25 characters)
- **Reversibility**: **Yes** (exact coordinate reconstruction)
- **Encoding**: Hexadecimal (lowercase)
- **Cross-validation**: Can verify XI tag correctness

---

## XS Tag - Splicetag

### Purpose
Reversible encoding of internal splice junction coordinates.

### Format
```
XS:Z:[8-chr-hash][strand].[coord1-hex].[coord2-hex].[coord3-hex]...
```

### Components
- **8-chr-hash**: First 8 characters of RefGet chromosome hash
- **strand**: `p` (plus) or `m` (minus)
- **coord-hex**: Splice junction coordinates in hexadecimal (lowercase)
  - For + strand: exon ends, then exon starts (alternating)
  - For - strand: exon starts, then exon ends (alternating)

### Example
```bash
# Input (3 exons on + strand)
Exon 1: 1000-1200
Exon 2: 2000-2150
Exon 3: 3000-3500

# Internal coordinates (splice junctions)
1200 (exon1 end), 2000 (exon2 start), 2150 (exon2 end), 3000 (exon3 start)
= 0x4b0, 0x7d0, 0x866, 0xbb8

# XS tag
XS:Z:aKF498dAp.4b0.7d0.866.bb8
```

### Visual Breakdown
```
XS:Z:aKF498dAp.4b0.7d0.866.bb8
     │       │ │   │   │   │
     │       │ └───┴───┴───┴─── Splice junction coordinates (hex)
     │       └───────────────── Strand: +
     └─────────────────────── Chromosome: aKF498dA (8-char RefGet hash)
```

### Encoding Algorithm
```python
def generate_splicetag(chr_refget_32: str, strand: str,
                       exons: List[Tuple[int, int]]) -> str:
    """Generate XS tag splicetag (multi-exon only)"""
    if len(exons) < 2:
        return "None"  # Single-exon transcripts have no splice junctions

    chr_hash_8 = chr_refget_32[:8]
    strand_char = 'p' if strand == '+' else 'm'

    # Extract internal splice junction coordinates
    coords = []
    for i in range(len(exons) - 1):
        coords.append(exons[i][1])      # Current exon end
        coords.append(exons[i+1][0])    # Next exon start

    # Convert to hex
    hex_coords = [format(coord, 'x') for coord in coords]
    return f"{chr_hash_8}{strand_char}.{'.'.join(hex_coords)}"
```

### Decoding Algorithm
```python
def decode_splicetag(xs_tag: str) -> Tuple[str, str, List[int]]:
    """Decode XS tag to splice junction coordinates"""
    chr_hash_8 = xs_tag[:8]
    strand = '+' if xs_tag[8] == 'p' else '-'
    hex_coords = xs_tag[10:].split('.')  # Skip 8-char hash + strand + dot
    coords = [int(coord, 16) for coord in hex_coords]
    return chr_hash_8, strand, coords
```

### Full Exon Reconstruction (XB + XS)
```python
def reconstruct_exons(xb_tag: str, xs_tag: str) -> List[Tuple[int, int]]:
    """Reconstruct full exon structure from XB and XS tags"""
    # Decode boundary tag
    chr_hash, strand, five_prime, three_prime = decode_boundarytag(xb_tag)

    # Single-exon transcript
    if xs_tag == "None":
        return [(five_prime, three_prime)]

    # Multi-exon transcript
    _, _, splice_coords = decode_splicetag(xs_tag)

    # Reconstruct exons
    exons = []
    exons.append((five_prime, splice_coords[0]))  # First exon

    # Internal exons
    for i in range(1, len(splice_coords) - 1, 2):
        exons.append((splice_coords[i], splice_coords[i+1]))

    # Last exon
    exons.append((splice_coords[-1], three_prime))

    return exons
```

### Properties
- **Length**: Variable (depends on exon count)
- **Reversibility**: **Yes** (exact splice junction reconstruction)
- **Encoding**: Hexadecimal (lowercase)
- **Presence**: Multi-exon transcripts only

---

## XT Tag - Transcript Group

### Purpose
Biological clustering of transcripts with fuzzy boundaries (TSS variation, degradation, etc).

### Format
```
XT:Z:[32-character-hash]
```

### Clustering Modes

IsoTag v2.0 supports three biological clustering modes:

| Mode | Position Used | Use Case | Data Type |
|------|---------------|----------|-----------|
| **5prime** | 5' transcript end | TSS clustering | CAGE, Cap-seq |
| **middle** | Middle position | General clustering | RNA-seq, Iso-Seq |
| **3prime** | 3' transcript end | TES clustering | PolyA-seq, 3'-seq |

### Serialization (before hashing)

**Mode: middle (default)**
```
[32-chr-refget]|[strand]|[rounded_middle]|[rounded_exon_total]|[rounded_span]|[splice_junctions]
```

**Mode: 5prime**
```
[32-chr-refget]|[strand]|[rounded_5prime]|[rounded_exon_total]|[rounded_span]|[splice_junctions]
```

**Mode: 3prime**
```
[32-chr-refget]|[strand]|[rounded_3prime]|[rounded_exon_total]|[rounded_span]|[splice_junctions]
```

### Example
```bash
# Input (middle mode)
Chromosome: chr1 (RefGet: aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2)
Strand: +
Exons: 1000-1200, 2000-2150, 3000-3500
Middle position: 2250 → rounded to 0 (quantum 10000)
Exon total: 850 → rounded to 1000 (quantum 1000)
Genomic span: 2501 → rounded to 0 (quantum 10000)
Splice junctions: 1200,2000,2150,3000

# Serialization
"aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2|+|0|1000|0|1200|2000|2150|3000"

# Final XT tag
XT:Z:266CbPqmZz8eS-EzT4xtnYtmm-SoIhnL
```

### Quantization Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| **position-quantum** | 10000 bp | Round middle/5'/3' positions |
| **span-quantum** | 10000 bp | Round genomic span |
| **exon-quantum** | 1000 bp | Round total exon length |

### Encoding Algorithm
```python
def generate_xt_group_id(chr_refget_32: str, strand: str, exons: List[Tuple[int, int]],
                         mode: str = 'middle', pos_quantum: int = 10000,
                         span_quantum: int = 10000, exon_quantum: int = 1000) -> str:
    """Generate XT tag transcript group ID"""
    # Calculate position based on mode
    if mode == '5prime':
        position = exons[0][0] if strand == '+' else exons[-1][1]
    elif mode == '3prime':
        position = exons[-1][1] if strand == '+' else exons[0][0]
    else:  # middle
        start = min(e[0] for e in exons)
        end = max(e[1] for e in exons)
        position = (start + end) // 2

    # Round position
    rounded_position = round(position / pos_quantum) * pos_quantum

    # Calculate exon total length
    exon_total = sum(end - start + 1 for start, end in exons)
    rounded_exon_total = round(exon_total / exon_quantum) * exon_quantum

    # Calculate genomic span
    genomic_span = max(e[1] for e in exons) - min(e[0] for e in exons) + 1
    rounded_span = round(genomic_span / span_quantum) * span_quantum

    # Extract splice junctions
    splice_junctions = []
    for i in range(len(exons) - 1):
        splice_junctions.append(str(exons[i][1]))
        splice_junctions.append(str(exons[i+1][0]))

    # Build serialization
    parts = [chr_refget_32, strand, str(rounded_position),
             str(rounded_exon_total), str(rounded_span)]
    parts.extend(splice_junctions)

    serialization = '|'.join(parts)
    return sha512t24u(serialization.encode('utf-8'))
```

### Properties
- **Length**: Always 32 characters
- **Reversibility**: No (one-way hash)
- **Purpose**: Fuzzy boundary clustering for biological analysis
- **Sample-independent**: Same quantization = same group across experiments

---

## XC Tag - Gene/Locus ID

### Purpose
Pure location-based gene/locus identifier that groups ALL transcripts at the same genomic location, regardless of isoform structure. Suitable for use as a gene ID or locus ID without requiring annotation databases.

### Format
```
XC:Z:[32-character-hash]
```

### Serialization (before hashing)
```
[32-chr-refget-hash]|[strand]|[start_bin]|[end_bin]
```

Where:
- **start_bin** = `transcript_start // xc_bin_size`
- **end_bin** = `transcript_end // xc_bin_size`
- **xc_bin_size** = configurable (default: 10000bp)

### Example
```bash
# Input
Chromosome: chr1 (RefGet: aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2)
Strand: +
Exons: 1000000-1001200, 1005000-1005150, 1010000-1020000
# start_bin = 1000000 // 10000 = 100
# end_bin = 1020000 // 10000 = 102

# Serialization
"aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2|+|100|102"

# Final XC tag
XC:Z:a7Bf9xK2mP3qR5tN8wY1zC4dF6hJ0lO
```

### XC vs XT Comparison

| Feature | XT (Transcript Group) | XC (Gene/Locus) |
|---------|----------------------|-----------------|
| Position | Binned (mode-based) | Binned (start/end) |
| Strand | ✓ | ✓ |
| Exon count/lengths | Quantized | **Ignored** |
| Splice junctions | ✓ | **Ignored** |
| Genomic span | Quantized | **Ignored** |
| Same gene, different isoforms | May differ | **Always same** ✅ |
| Clustering level | Isoform group | Gene/locus |
| Use case | Isoform clustering | Gene-level ID |

### Bin Size Parameter

| Bin Size | Clustering Level | Use Case |
|----------|-----------------|----------|
| **10bp** | Very fine | Splice wobble detection |
| **100bp** | Fine | Isoform discovery with tolerance |
| **1kb** | Moderate | Locus-level grouping |
| **10kb** (default) | Gene-level | Gene-level analysis ⭐ |
| **100kb** | Very coarse | May merge nearby genes |

### Encoding Algorithm
```python
def generate_xc_cluster_id(chr_refget_32: str, strand: str,
                           exons: List[Tuple[int, int]],
                           xc_bin_size: int = 10000) -> str:
    """Generate XC tag gene/locus cluster ID"""
    sorted_exons = sorted(exons)
    start_bin = sorted_exons[0][0] // xc_bin_size
    end_bin = sorted_exons[-1][1] // xc_bin_size

    serialization = f"{chr_refget_32}|{strand}|{start_bin}|{end_bin}"
    return sha512t24u(serialization.encode('utf-8'))
```

### Properties
- **Length**: Always 32 characters
- **Character set**: Base64URL (A-Z, a-z, 0-9, -, _)
- **Reversibility**: No (one-way hash)
- **Key property**: All isoforms at the same genomic locus produce the same XC tag

---

## XV Tag - Variants

### Purpose
Individual variant identifiers for linking isoforms to genomic variations.

### Format
```
XV:Z:[variant1-32chars].[variant2-32chars].[variant3-32chars]...
```

### Variant Serialization (before hashing)
```
[32-chr-refget]:[position]:[ref]>[alt]
```

### Example
```bash
# Variant 1: chr1:1050:A>G
Serialization: "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2:1050:A>G"
Hash: Q4fUfjJXgQwpSxFgeGVowhJaLTVg3Fqk

# Variant 2: chr1:2075:C>T
Serialization: "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2:2075:C>T"
Hash: 85IT2XfuQ7nofo_xMFc4PzCZNtgAcZ96

# Final XV tag
XV:Z:Q4fUfjJXgQwpSxFgeGVowhJaLTVg3Fqk.85IT2XfuQ7nofo_xMFc4PzCZNtgAcZ96
```

### Variant Format Rules
- **SNPs**: `pos:ref>alt` (e.g., `1050:A>G`)
- **Deletions**: `pos:seq>-` (e.g., `2075:ATG>-`)
- **Insertions**: `pos:->seq` (e.g., `3100:->TACG`)
- **MNVs**: `pos:ref>alt` (e.g., `4000:AC>GT`)

### Encoding Algorithm
```python
def generate_variant_id(chr_refget_32: str, position: int,
                        ref: str, alt: str) -> str:
    """Generate individual variant ID"""
    variant_str = f"{chr_refget_32}:{position}:{ref}>{alt}"
    return sha512t24u(variant_str.encode('utf-8'))
```

### Properties
- **Length**: Variable (32 chars × number of variants + dots)
- **Reversibility**: No (one-way hash)
- **Presence**: Optional (only when variants detected)
- **Separator**: Dot (`.`) between variant IDs

---

## RefGet Chromosome Hashing

### Overview
IsoTag v2.0 uses **RefGet-based chromosome hashing** to solve naming inconsistencies (chr1 vs Chr1 vs CHR1 vs 1).

### Algorithm
```python
def calculate_refget_id(chromosome_sequence: str) -> str:
    """Calculate RefGet ID for chromosome sequence"""
    seq_bytes = chromosome_sequence.upper().encode('ascii')
    return sha512t24u(seq_bytes)  # Returns 32-char base64url string
```

### Hash Lengths Used

| Tag Type | Hash Length | Purpose |
|----------|-------------|---------|
| **XB, XS** | 8 chars | Compact reversible encoding |
| **XI, XT, XC, XV** | 32 chars | Full chromosome hash in serialization |

### RefGet Cache Format
```json
{
  "metadata": {
    "genome": "hg38",
    "generated": "2025-10-02T12:00:00",
    "total_mappings": 150
  },
  "refget_mapping": {
    "chr1": "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2",
    "Chr1": "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2",
    "CHR1": "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2",
    "1": "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"
  }
}
```

### Universal Compatibility
```bash
# Different chromosome names → Same RefGet hash
chr1  → SQ.aKF498dA... → Same tags!
Chr1  → SQ.aKF498dA... → Same tags!
CHR1  → SQ.aKF498dA... → Same tags!
1     → SQ.aKF498dA... → Same tags!
```

---

## Encoding Algorithms

### VRS-Compatible sha512t24u
```python
def sha512t24u(blob: bytes) -> str:
    """
    Generate base64url-encoded, truncated SHA-512 digest

    This is the official GA4GH VRS algorithm for generating
    compact, collision-resistant identifiers.

    Args:
        blob: Input bytes to hash

    Returns:
        32-character base64url-encoded string
    """
    import hashlib
    import base64

    digest_size = 24  # 24 bytes = 192 bits
    digest = hashlib.sha512(blob).digest()
    truncated = digest[:digest_size]
    encoded = base64.urlsafe_b64encode(truncated)
    return encoded.decode("ascii").rstrip('=')  # Remove padding
```

### Coordinate Encoding (Hex)
```python
def encode_coordinate(position: int) -> str:
    """Encode genomic coordinate as lowercase hexadecimal"""
    return format(position, 'x')

# Examples
encode_coordinate(1000)  # "3e8"
encode_coordinate(4100)  # "1004"
encode_coordinate(65535) # "ffff"
```

---

## Decoding Algorithms

### XB Tag Decoder
```python
def decode_xb_tag(xb_value: str) -> dict:
    """
    Decode XB boundary tag

    Args:
        xb_value: XB tag value (e.g., "aKF498dAp.3e8.1004")

    Returns:
        {
            'chr_hash_8': '8-character chromosome hash',
            'strand': '+' or '-',
            'five_prime_end': int,
            'three_prime_end': int,
            'genomic_span': int
        }
    """
    chr_hash_8 = xb_value[:8]
    strand = '+' if xb_value[8] == 'p' else '-'
    coords = xb_value[10:].split('.')
    five_prime = int(coords[0], 16)
    three_prime = int(coords[1], 16)

    return {
        'chr_hash_8': chr_hash_8,
        'strand': strand,
        'five_prime_end': five_prime,
        'three_prime_end': three_prime,
        'genomic_span': three_prime - five_prime + 1
    }
```

### XS Tag Decoder
```python
def decode_xs_tag(xs_value: str) -> dict:
    """
    Decode XS splicetag

    Args:
        xs_value: XS tag value (e.g., "aKF498dAp.4b0.7d0.866.bb8")

    Returns:
        {
            'chr_hash_8': '8-character chromosome hash',
            'strand': '+' or '-',
            'splice_coordinates': [int, int, ...]
        }
    """
    if xs_value == "None":
        return None

    chr_hash_8 = xs_value[:8]
    strand = '+' if xs_value[8] == 'p' else '-'
    hex_coords = xs_value[10:].split('.')
    coordinates = [int(coord, 16) for coord in hex_coords]

    return {
        'chr_hash_8': chr_hash_8,
        'strand': strand,
        'splice_coordinates': coordinates
    }
```

### Full Exon Reconstruction
```python
def reconstruct_exons_from_tags(xb_value: str, xs_value: str) -> List[Tuple[int, int]]:
    """
    Reconstruct complete exon structure from XB and XS tags

    Args:
        xb_value: XB tag value
        xs_value: XS tag value

    Returns:
        List of (start, end) tuples for each exon
    """
    xb_decoded = decode_xb_tag(xb_value)
    five_prime = xb_decoded['five_prime_end']
    three_prime = xb_decoded['three_prime_end']

    # Single-exon transcript
    if xs_value == "None":
        return [(five_prime, three_prime)]

    # Multi-exon transcript
    xs_decoded = decode_xs_tag(xs_value)
    splice_coords = xs_decoded['splice_coordinates']

    exons = []
    # First exon
    exons.append((five_prime, splice_coords[0]))

    # Internal exons (coordinates come in pairs: start, end)
    for i in range(1, len(splice_coords) - 1, 2):
        exons.append((splice_coords[i], splice_coords[i+1]))

    # Last exon
    exons.append((splice_coords[-1], three_prime))

    return exons
```

---

## Reference Implementation

Full implementation available in:
- `isotag.py`: Tag generation
- `decode_tags.py`: Tag decoding and reconstruction
- `vrs_compat.py`: VRS-compliant sha512t24u algorithm

---

## Version History

- **v2.2.0** (2026-02-09): Added XC gene/locus tag (pure location-based clustering)
- **v2.0.0** (2025-10-02): Added XB/XS reversible tags, universal chromosome hashing, biological clustering
- **v1.0.0** (2025-08-20): Initial release with XI/XV tags

---

## References

- **GA4GH RefGet**: https://samtools.github.io/hts-specs/refget.html
- **VRS Specification**: https://vrs.ga4gh.org
- **SAM/BAM Format**: https://samtools.github.io/hts-specs/SAMv1.pdf

---

**IsoTag Tag Format Specification v2.2**
**Last Updated**: February 9, 2026
**Maintained by**: LSBDT Team
