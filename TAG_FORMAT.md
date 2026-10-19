# IsoTag Tag Format Specification v2.5

This document provides the complete technical specification for all IsoTag BAM/SAM tags.

## Table of Contents
- [Overview](#overview)
- [Tag Summary](#tag-summary)
- [XI Tag - Structure ID](#xi-tag---structure-id)
- [XB Tag - Boundary Tag](#xb-tag---boundary-tag)
- [XS Tag - Splicetag](#xs-tag---splicetag)
- [XT Tag - Transcript Group](#xt-tag---transcript-group)
- [XC Tag - Cluster ID (v12.0, end-free)](#xc-tag---cluster-id-v120-end-free)
- [XV Tag - Variants](#xv-tag---variants)
- [RefGet Chromosome Hashing](#refget-chromosome-hashing)
- [Encoding Algorithms](#encoding-algorithms)
- [Decoding Algorithms](#decoding-algorithms)

---

## Overview

IsoTag v2.5 uses **eight custom BAM tags** to encode transcript isoform structure, splice junctions, biological clustering, locus identity, TSS/polyA sites, and variants:

- **XI**: Unique isoform structure identifier (32-char hash)
- **XB**: Reversible boundary tag (8-char chr hash + hex-encoded ends)
- **XS**: Reversible splicetag (8-char chr hash + hex-encoded splice junctions)
- **XT**: Biological transcript group (32-char hash, exact splice junction matching)
- **XC**: Splice-structure cluster identifier (32-char hash, end-free: internal span + binned intron/internal-exon lengths; v12.0)
- **X5**: TSS cluster (32-char hash, binned 5′ end; default 200bp bins)
- **X3**: PolyA cluster (32-char hash, binned 3′ end; default 200bp bins)
- **XV**: GA4GH VRS v2 Allele identifiers (SNVs/substitutions, comma-separated, optional; v2.4+)

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
| **XB** | Z (string) | Variable | **Yes** | 5′/3′ transcript boundary coordinates |
| **XS** | Z (string) | Variable | **Yes** | Internal splice junction coordinates |
| **XT** | Z (string) | 32 chars | No | Transcript group (exact splice junctions) |
| **XC** | Z (string) | 32 chars | No | Splice-structure cluster, end-free (v12.0) |
| **X5** | Z (string) | 32 chars | No | TSS cluster (binned 5′ end) |
| **X3** | Z (string) | 32 chars | No | PolyA cluster (binned 3′ end) |
| **XV** | Z (string) | Variable | No | GA4GH VRS v2 Allele IDs (comma-separated; v2.4+) |

### Tag Presence Rules

- **XI**: Always present (required for all reads with CIGAR)
- **XB**: Always present (required for all reads with CIGAR)
- **XS**: Present only for multi-exon transcripts (2+ exons)
- **XT**: Always present (required for all reads with CIGAR)
- **XC**: Always present (required for all reads with CIGAR)
- **X5**: Always present (required for all reads with CIGAR)
- **X3**: Always present (required for all reads with CIGAR)
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

> **⚠️ Tag name collision warning:** `XS` is also used by BWA-MEM (`XS:i`, insert size for paired-end reads) and STAR (`XS:A`, strand of the gene). IsoTag uses `XS:Z` (string type), which is a distinct SAM type. These tags coexist safely in a BAM if written by separate tools, but downstream parsers that query by tag name `XS` without checking the type may return unexpected results. Check your aligner's tag usage before running IsoTag on pre-aligned BAMs.

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
Biological transcript group using exact splice junction coordinates. Transcripts with identical splice junctions (same exon boundaries) share the same XT tag. **Not wobble-tolerant**: any ±1bp Nanopore junction wobble produces a different XT hash. Use XC (v12.0) for end-independent structure grouping with 10 bp bin tolerance.

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

## XC Tag - Cluster ID (v12.0, end-free)

> **Status: Experimental** — XC is available and functional but has known limitations documented below. It is a splice-structure cluster id, not a gene id.
>
> **Breaking change in v12.0 (October 2026, IsoTag v2.5.0):** the XC serialization was redesigned again. v11.0 hashed the binned length of *every* exon, including the first and last one. Those two lengths are set by where the read starts and ends, so v11.0 silently required 5′/3′ agreement within a 10 bp bin: on HPRC Kinnex HG00097 (10.9 M multi-exon reads) only 12.4% of reads shared an XC v11.0 value with a GENCODE v50 transcript, versus 46.9% for XT and 76.1% for the exact junction chain XS. v12.0 drops the terminal exon lengths and anchors on the internal span instead of the read midpoint; the same measurement gives 76.9%. v11.0 and v12.0 hashes are not comparable — re-tag with `isotag.py` (v11.0 remains available as `--xc-version 11`).

### Purpose
Groups reads that share the same **internal splice structure** at the same locus, independent of their 5′/3′ ends. Only junction coordinates enter the serialization: terminal exon lengths, TSS and polyA positions do not (they are carried by X5/X3/XB). Compared with XS (exact junction coordinates), XC tolerates small junction shifts through 10 bp length bins; compared with XT, it is not split by read ends.

### Format
```
XC:Z:[32-character-hash]
```

### Serialization (before hashing) — v12.0

Multi-exon reads (K junctions, K ≥ 1). Exons are sorted by coordinate (1-based, inclusive, as in XS); `d_i` is the donor (end) of exon *i*, `a_i` the acceptor (start) of exon *i+1*:
```
[8-chr-refget][s]|[ispan_mid_bin]|[ispan_bin]|[intron_1_bin]|[exon_2_bin]|[intron_2_bin]|...|[exon_K_bin]|[intron_K_bin]
```
Where:
- `[8-chr-refget]` = first 8 characters of the chromosome RefGet digest (same prefix XS/XB use); `[s]` = `p` (+ strand) or `m` (− strand)
- `ispan_mid_bin = ((d_1 + a_K) // 2) // xc_position_quantum` — midpoint of the **internal span** (first donor to last acceptor), default 1000 bp bins
- `ispan_bin = (a_K − d_1) // xc_bin_size` — internal span length, default 10 bp bins
- `intron_i_bin = (a_i − d_i − 1) // xc_bin_size` — every intron length
- `exon_i_bin = (d_i − a_{i−1} + 1) // xc_bin_size` — every **internal** exon length (exons 2..K); the first and last exon are not encoded
- The intron terms keep "exon-shift" structures apart (same internal exon lengths, different intron partition); the span term separates alternative first/last exons that leave the internal exons unchanged.

Mono-exon reads have no internal structure; XC falls back to a coarse locus id:
```
[8-chr-refget][s]|S|[midpoint // xc_position_quantum]|[length // xc_position_quantum]
```

### Example
```bash
# Input (v12.0 defaults: xc_position_quantum=1000, xc_bin_size=10)
Chromosome: chr1 (test-suite fallback hash; with a RefGet mapping the digest prefix is used), strand +
Exons: 1000-1200, 5000-5150, 10000-20000
Donors d = (1200, 5150); acceptors a = (5000, 10000)
ispan_mid_bin = ((1200 + 10000) // 2) // 1000 = 5600 // 1000 = 5
ispan_bin     = (10000 - 1200) // 10 = 880
intron_1_bin  = (5000 - 1200 - 1) // 10 = 379
exon_2_bin    = (5150 - 5000 + 1) // 10 = 15
intron_2_bin  = (10000 - 5150 - 1) // 10 = 484
# First exon (201 bp) and last exon (10001 bp) are NOT encoded.

# Serialization
"Tf9Oo0Dwp|5|880|379|15|484"

# Final XC tag
XC:Z:H2CGaj_IMJtHrWMK1vcy-5cHZQqUmoVg
```
A read of the same transcript that starts 300 bp downstream or ends 2 kb earlier has the same XC; a read missing the last junction (truncated before exon 3) does not — XC is chain-level, not gene-level.

### What changes XC and what does not (v12.0)
| Change | XC v12.0 | XC v11.0 | XS | XT |
|--------|----------|----------|----|----|
| 5′ or 3′ end moves (any amount, same junctions) | same | differs if a terminal exon crosses a 10 bp bin | same | differs if a 10 kb / 1 kb bin is crossed |
| One junction shifts by 1 bp (wobble) | same with probability ≈ 0.81 (two length terms, each crosses a 10 bp bin with p = 0.1) | similar | differs | differs |
| Junction added / removed (truncation past a junction, skipped exon) | differs | differs | differs | differs |
| Internal exon length differs by ≥ 10 bp | differs | differs | differs | differs |

### Measured behaviour (HPRC Kinnex HG00097, GENCODE v50; 2026-10-12)
| Metric | XC v11.0 | XT | XS | XC v12.0 |
|--------|----------|----|----|----------|
| Multi-exon reads whose tag value equals an annotated transcript's (10.9 M reads) | 12.4% | 46.9% | 76.1% | 76.9% |
| GENCODE v50 multi-exon transcripts hit by ≥ 1 read | 3.1% | 10.5% | 18.4% | 20.1% |
| Same, all 206 HPRC samples (2.9 G reads) | — | — | 54.6% | 55.8% |
| Distinct XC values per XS junction chain (chr22, 400,680 multi-exon reads) | 2.70 (max 267) | — | 1 | **1.00 (max 1)** |
| Same, full HG00097 (10.9 M multi-exon reads, 706,712 chains) | 2.84 (max 412) | — | 1 | **1.00 (max 1)** |
| GENCODE v50 values shared by ≥ 2 genes (of ~600 k) | 15 | 74 | 154 | 162 |
| GENCODE v50 values shared by ≥ 2 transcripts | 15,099 | 10,615 | 16,546 | 29,417 |
| Gene-level purity of multi-read XC clusters (full HG00097 vs GENCODE v50 gene bodies) | 75.9% | — | — | 77.2% |
| Gene completeness (genes held by a single XC value) | 19.0% | — | — | 20.7% |
| ARI / V-measure vs genes (full HG00097) | 0.080 / 0.558 | — | — | 0.087 / 0.598 |

The 29,417 shared v12.0 values on the annotation side split into 12,129 that hold transcripts with the **same junction chain** (they differ only in their ends — the intended merge), 17,126 that hold different chains of the **same gene** (internal differences smaller than one 10 bp bin, e.g. NAGNAG acceptor shifts, binned together) and 162 that span two genes (same order as XS's 154: overlapping or read-through loci). Use XS when isoforms differing by < 10 bp at a junction must be kept apart.

### Legacy v11.0 (`--xc-version 11`)
```
[32-chr-refget]|[strand]|[middle_bin]|[exon_len_bin_1]|...|[exon_len_bin_N]
middle_bin = (transcript_start + transcript_end) // 2 // xc_position_quantum ; exon_len_bin_i = exon_i_length // xc_bin_size
```
Every exon length is encoded, so the first and last bins move with the read's ends. Kept only for re-tagging legacy datasets; the position-only baseline `--xc-midpoint-only` also uses this formula without the length bins.

### Known Limitations (v12.0)
1. **Not a gene id.** Long genes carry many XC values (one per internal splice structure); isoforms of one gene that differ in any junction get different XC values. Group by locus with X5/X3 or by annotation if gene-level ids are needed.
2. **Bin boundaries merge near-identical isoforms.** Internal differences smaller than `xc_bin_size` (10 bp) can fall into the same bin — 17,126 GENCODE v50 values hold two or more distinct chains of the same gene. Use XS for exact structure.
3. **Mono-exon reads are locus-level only** (1 kb midpoint and length bins). They have no internal structure to encode.
4. **Cross-gene collisions** are inherited from the junction chain itself (overlapping genes, read-through): 162 of ~600 k GENCODE values, the same order as XS.
5. **Hash versions are not comparable.** BAMs tagged with v9.0, v10.0 or v11.0 must be re-tagged.

### Properties
- **Length**: Always 32 characters
- **Character set**: Base64URL (A-Z, a-z, 0-9, -, _)
- **Reversibility**: No (one-way hash); the junction chain itself is recoverable from XS
- **CLI flags**: `--xc-version 12` (default) / `11` (legacy); `--xc-position-quantum 1000` (internal-span midpoint bin); `--xc-bin-size 10` (length bins); `--xc-midpoint-only` (position-only baseline, v11 formula)

---

## X5 Tag - TSS Cluster

### Purpose
Transcription start site (TSS) cluster — groups reads sharing the same binned 5′ end position.
Enables TSS discovery and CAGE data comparison without external annotation.

### Format
```
X5:Z:[32-character-hash]
```

### Serialization (before hashing)
```
[32-chr-refget-hash]|[strand]|[tss_bin]
```

Where:
- **tss_bin** = `transcript_5prime_end // x5x3_bin_size`
- **x5x3_bin_size** = configurable (default: 200bp; use 35bp for finer resolution)
- 5′ end = leftmost exon start on + strand; rightmost exon end on − strand

### Validation
- 43.3% of X5 clusters (200bp bins) overlap FANTOM5 CAGE peaks within ±100bp
- 54.9% at 35bp bins with ±100bp window; 45.3% at 35bp bins with matched ±35bp window
- Genic-space permutation null: 43.1× enrichment over random genic positions (Z=1903)

---

## X3 Tag - PolyA Cluster

### Purpose
Poly-adenylation site (polyA) cluster — groups reads sharing the same binned 3′ end position.
Enables polyA site discovery and PolyA_DB / PolyASite comparison without external annotation.

### Format
```
X3:Z:[32-character-hash]
```

### Serialization (before hashing)
```
[32-chr-refget-hash]|[strand]|[polya_bin]
```

Where:
- **polya_bin** = `transcript_3prime_end // x5x3_bin_size`
- 3′ end = rightmost exon end on + strand; leftmost exon start on − strand

### Validation
- 29.16% of X3 clusters overlap PolyA_DB v4 sites within ±100bp

---

## XV Tag - Variants

**Changed in v2.4 (VRS v2 migration):** XV now emits full GA4GH VRS v2 Allele
identifiers instead of 32-char IsoTag-internal hashes. This is a breaking format change —
v2.3-and-earlier XV values are not comparable to v2.4+ XV values.

### Purpose
Per-read variant identifiers, computed as canonical GA4GH VRS v2 Allele IDs so they can be
matched directly against ClinVar/ClinGen and other VRS-indexed external databases.

### Format
```
XV:Z:ga4gh:VA.<32-char-digest>[,ga4gh:VA.<32-char-digest>...]
```
Multiple variants on one read are comma-separated, sorted, and deduplicated. A comma is
unambiguous as a separator because GA4GH identifiers already contain a period (`ga4gh:VA.`).

### Example
```bash
# NC_000017.11:g.43106487A>G (BRCA1, ClinVar golden fixture)
XV:Z:ga4gh:VA.gbvJw0s4OeAvloCeAM6BNNvrjFC_Dhc8
```

### Computation
1. Extract read-supported alleles from CIGAR + MD + SEQ, synchronizing the genomic cursor
   (advanced by every CIGAR op, including `N` introns) against the MD-axis cursor (advanced
   only by `M`/`D`) so mismatches downstream of a splice junction land at the correct position.
2. Convert to 0-based inter-residue coordinates (VRS convention): `start = SAM_POS - 1`.
3. Trim shared literal ref/alt flanks. If trimming would produce a pure insertion or deletion,
   the variant is **omitted** — canonical VRS normalization for indels requires reference-aware
   left-alignment not yet implemented here.
4. Digest per GA4GH's `sha512t24u` algorithm over canonical JSON: a `SequenceLocation` digest
   (chr RefGet accession + start/end), then an `Allele` digest over `{location: <location
   digest>, state: {type: LiteralSequenceExpression, sequence: alt}}`.
5. Prefix with `ga4gh:VA.`.

### Safety Boundary
- **Requires an unmasked RefGet mapping.** IsoTag's default ambiguous-base masking (see below)
  changes the chromosome sequence hash and therefore produces XV identifiers that will not
  match ClinVar/ClinGen. `isotag.py` refuses to emit XV against a masked mapping — regenerate
  with `--keep-ambiguous-bases` from the exact canonical FASTA. XI/XB/XS/XT/XC/X5/X3 are
  unaffected by this restriction and may still use a masked mapping.
- **Requires `--refget` or `--genome`.** Legacy chromosome-name hashing (no genome/refget
  provided) has no RefGet accession to build a VRS `SequenceReference` from.
- **Indels are silently omitted from XV**, not emitted with a wrong or placeholder identifier.

### Properties
- **Length**: Variable (one or more `ga4gh:VA.<32-char>` IDs, comma-separated)
- **Reversibility**: No (one-way digest)
- **Presence**: Optional — only for SNV/substitution candidates, only when variant detection
  is enabled and an unmasked RefGet mapping is available

---

## RefGet Chromosome Hashing

### Overview
IsoTag uses a **RefGet-inspired chromosome hashing scheme** to solve naming inconsistencies (chr1 vs Chr1 vs CHR1 vs 1).

**Important difference from canonical GA4GH RefGet:** IsoTag masks IUPAC ambiguous bases (R, Y, S, W, K, M, B, D, H, V) to 'N' before hashing. The GA4GH RefGet specification hashes sequences as-is. As a result, IsoTag chromosome hashes **will not match** canonical GA4GH RefGet API IDs for any chromosome that contains ambiguous bases.

**Why the difference exists:** UCSC hg38 and NCBI GRCh38 represent the same genome but use different ambiguous base encodings at some positions. Masking to 'N' before hashing ensures both assemblies produce identical XC/XI/XT tags, enabling cross-assembly comparison. Canonical GA4GH compatibility was sacrificed for cross-assembly reproducibility.

**Exception: XV requires the unmasked mapping.** Since the v2.4 VRS v2 migration, XV Allele
IDs must be computed from the exact canonical reference sequence to match ClinVar/ClinGen —
masking breaks that match. `isotag.py` refuses to emit XV against a masked mapping; regenerate
with `--keep-ambiguous-bases` from the exact canonical FASTA if you need real XV values. This
does not change the default masking behavior for XI/XB/XS/XT/XC/X5/X3.

### Algorithm
```python
def calculate_refget_id(chromosome_sequence: str) -> str:
    """
    Calculate IsoTag chromosome hash (RefGet-inspired, with ambiguous base masking).

    NOTE: This does NOT produce canonical GA4GH RefGet IDs.
    Ambiguous IUPAC bases are masked to 'N' before hashing to ensure
    UCSC hg38 and NCBI GRCh38 produce identical hashes.
    """
    import re
    AMBIGUOUS_BASES = re.compile(r'[RYSWKMBDHV]', re.IGNORECASE)
    normalized = AMBIGUOUS_BASES.sub('N', chromosome_sequence.upper())
    return sha512t24u(normalized.encode('ascii'))  # Returns 32-char base64url string
```

> **v2.2.0 → v2.2.1 migration warning:** RefGet JSONs generated with IsoTag v2.0–v2.2.0 did NOT mask ambiguous bases. Hashes from those versions are incompatible with v2.2.1+ hashes. Any `.json` file generated before v2.2.1 must be regenerated with `isotag_refget.py` to produce correct cross-assembly tags.

### Hash Lengths Used

| Tag Type | Hash Length | Purpose |
|----------|-------------|---------|
| **XB, XS, XC (v12.0)** | 8 chars | Compact reversible encoding (XB/XS); XC v12.0 serialization prefix |
| **XI, XT, XC (v11.0 legacy)** | 32 chars | Full chromosome hash in serialization |
| **XV** | 32-char base64url per digest (`ga4gh:VA.<digest>`, comma-separated) | GA4GH VRS v2 Allele ID (see XV Tag section) |

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

- **v2.5.0** (2026-10-12): XC redesigned to v12.0 (end-free: internal-span anchor + binned intron/internal-exon lengths; terminal exon lengths dropped); `--xc-version 11` keeps the legacy formula
- **v2.4.0** (2026-08-05): X5/X3 tags; XC v11.0 (midpoint + binned exon lengths); XV migrated to GA4GH VRS v2
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
