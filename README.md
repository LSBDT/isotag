# IsoTag v2.3 - Universal Isoform Identification System 🚀

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![GA4GH](https://img.shields.io/badge/GA4GH-RefGet%20Compatible-green.svg)](https://samtools.github.io/hts-specs/refget.html)

A production-ready system for generating standardized, RefGet-compatible isoform identifiers for long-read transcript sequencing data. **v2.3** merges UCSC+NCBI RefGet JSONs and fixes a chromosome alias bug.

## 🎯 What's New in v2.3

### **Merged RefGet JSONs** 🔀 NEW!
- **Unified genome files**: UCSC and NCBI versions combined into single JSON (e.g., `mm39.GRCm39-refget.json`)
- **All naming conventions**: Look up by `chr1`, `1`, `CHR1`, `Chr1`, or `NC_000067.7` — all in one file
- **Bug fix**: Removed 280 bogus chromosome aliases incorrectly derived from NCBI accession numbers (e.g., `NC_000067.7` was generating `chr67` instead of being mouse chr1)
- **Merge tool**: `python3 isotag_refget.py -m file1.json -m file2.json -o merged.json`

#### Pre-built Merged RefGet JSONs (5 files)

| File | Genomes | Sequences | Mappings |
|------|---------|-----------|----------|
| `hg38.GRCh38-refget.json` | Human (UCSC+NCBI) | 705 | 2,525 |
| `hg19.GRCh37-refget.json` | Human (UCSC+NCBI) | 297 | 669 |
| `hs1.T2T-CHM13v2-refget.json` | Human T2T (UCSC+NCBI) | 25 | 124 |
| `mm39.GRCm39-refget.json` | Mouse (UCSC+NCBI) | 61 | 305 |
| `mm10.GRCm38-refget.json` | Mouse (UCSC+NCBI) | 239 | 503 |

## 🎯 What's New in v2.2

### **Gene/Locus Tag (XC Tag)** 🧬 NEW!
- **Pure location-based clustering**: Groups all transcripts at the same genomic locus into one XC ID
- **No annotation required**: Automatically assigns gene-level IDs to novel and known transcripts alike
- **Configurable resolution**: `--xc-bin-size` parameter (default: 10kb for gene-level, 10bp for splice wobble)
- **Same gene, different isoforms → Same XC**: Unlike XI/XT, XC ignores exon structure entirely
- **32-char hash**: Same format as XI/XT for consistency

```bash
# XC groups all isoforms at the same locus
Isoform A (5 exons): XI=aaa... XT=bbb... XC=zzz...  ← Same XC
Isoform B (3 exons): XI=ccc... XT=ddd... XC=zzz...  ← Same XC
Isoform C (4 exons): XI=eee... XT=fff... XC=zzz...  ← Same XC
```

## 🎯 What's New in v2.1

### **Unified Genome Build Support** 🔧 NEW!
- **Ambiguous base masking**: Automatically converts ambiguous IUPAC codes (R, Y, S, W, K, M, etc.) to 'N'
- **Cross-build compatibility**: UCSC hg38 and NCBI GRCh38 now produce identical chromosome hashes
- **Informative output**: Reports ambiguous bases masked per chromosome during RefGet cache generation
- **Optional override**: `--keep-ambiguous-bases` flag available if needed (not recommended)

## 🎯 What's New in v2.0

### **Reversible Splicetags (XS Tag)** ✨
- **Compact splice junction encoding**: hex-encoded coordinates for minimal storage
- **Full coordinate reconstruction**: decode exact splice sites from tag alone
- **Universal compatibility**: chr1/Chr1/CHR1/1 all map to same 8-char RefGet hash

### **Boundary Tags (XB Tag)** ✨
- **5'/3' transcript end encoding**: captures full transcript span
- **Cross-validation**: XB + XS = complete exon structure verification
- **No reference needed**: reconstruct coordinates without genome FASTA

### **Biological Clustering (XT Tag)** ✨
- **Three clustering modes**: 5prime (CAGE/TSS), middle (RNA-seq), 3prime (polyA/TES)
- **Fuzzy boundary handling**: groups transcripts with variable ends
- **Sample-independent**: consistent grouping across experiments

### **Universal Chromosome Hashing** 🌍
- **RefGet-based**: Uses GA4GH sequence identifiers
- **Auto-caching**: Generates RefGet cache from genome FASTA automatically
- **Cross-database**: Same chromosome sequence = same hash regardless of naming

## 📋 Tag Format (v2.2)

```bash
XI:Z:fuIF7PN23g2gq9sFxqhUNGnfOCZhkQJS              # Structure ID (32-char)
XB:Z:aKF498dAp.3e8.1004                           # Boundary tag (8-char chr + hex ends)
XS:Z:aKF498dAp.4b0.7d0.866.bb8                    # Splicetag (8-char chr + hex coords)
XT:Z:266CbPqmZz8eS-EzT4xtnYtmm-SoIhnL              # Transcript group (32-char)
XC:Z:a7Bf9xK2mP3qR5tN8wY1zC4dF6hJ0lO              # Gene/locus cluster (32-char)
XV:Z:Q4fUfjJXgQwpSxFgeGVowhJaLTVg3Fqk              # Variants (32-char, optional)
```

### Tag Breakdown

| Tag | Name | Format | Purpose |
|-----|------|--------|---------|
| **XI** | Structure ID | 32-char hash | Unique isoform structure identifier |
| **XB** | Boundary Tag | `[8-chr][s].[hex1].[hex2]` | Reversible 5'/3' transcript ends |
| **XS** | Splicetag | `[8-chr][s].[hex1].[hex2]...` | Reversible splice junction coordinates |
| **XT** | Transcript Group | 32-char hash | Biological clustering with fuzzy boundaries |
| **XC** | Gene/Locus ID | 32-char hash | Pure location-based gene/locus cluster |
| **XV** | Variants | 32-char hashes | Individual variant IDs (optional) |

**Legend**: `[8-chr]` = 8-char RefGet chromosome hash, `[s]` = strand (p/m), `[hex]` = hexadecimal coordinates

## 🚀 Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/LSBDT/isotag.git
cd isotag

# Install dependencies
pip install click pysam
```

**Requirements**: Python 3.8+, samtools 1.18+, click

### Basic Usage

```bash
# Standard tagging with auto RefGet cache generation
python3 isotag.py -i input.bam -o tagged.bam -g reference.fa

# Fast processing with pre-computed RefGet mapping (recommended)
python3 isotag.py -i input.bam -o tagged.bam -r hg38.GRCh38-refget.json

# Structure tags only (no variants)
python3 isotag.py -i input.bam -o tagged.bam -r hg38.GRCh38-refget.json --no-variants

# CAGE data (5' TSS clustering)
python3 isotag.py -i cage.bam -o tagged.bam -g reference.fa --clustermode 5prime

# PolyA data (3' TES clustering)
python3 isotag.py -i polya.bam -o tagged.bam -g reference.fa --clustermode 3prime

# Fine-grained gene/locus clustering (100bp XC bins)
python3 isotag.py -i input.bam -o tagged.bam -g reference.fa --xc-bin-size 100
```

### Decode Tags

```bash
# Decode boundary tag (XB)
python3 decode_tags.py -b "aKF498dAp.3e8.1004"

# Decode splicetag (XS)
python3 decode_tags.py -s "aKF498dAp.4b0.7d0.866.bb8"

# Reconstruct full exon structure
python3 decode_tags.py -b "aKF498dAp.3e8.1004" -s "aKF498dAp.4b0.7d0.866.bb8" --reconstruct

# With chromosome name lookup
python3 decode_tags.py -b "aKF498dAp.3e8.1004" -r hg38.GRCh38-refget.json
```

## 📖 Detailed Examples

### Example 1: Standard RNA-seq Processing

```bash
# Process RNA-seq BAM with automatic RefGet cache
python3 isotag.py -i rnaseq.bam -o tagged.bam -g hg38.fa

# Cache created at: ~/.isotag_cache/hg38_refget.json
# Subsequent runs use cache automatically
```

### Example 2: CAGE TSS Clustering

```bash
# Use 5' position for clustering (TSS-focused)
python3 isotag.py -i cage.bam -o tagged.bam -g mm39.fa --clustermode 5prime
```

### Example 3: Cross-Database Integration

```bash
# Lab A (UCSC naming: chr1, chr2)
python3 isotag.py -i labA.bam -o taggedA.bam -g ucsc_hg38.fa

# Lab B (Ensembl naming: 1, 2)
python3 isotag.py -i labB.bam -o taggedB.bam -g ensembl_hg38.fa

# Same chromosome sequences → Same RefGet hashes → Perfect compatibility! ✅
```

### Example 4: Decode and Reconstruct

```bash
# Extract tags from BAM
samtools view tagged.bam | grep "XB:Z:" | head -1 > tags.txt

# Decode to get coordinates
python3 decode_tags.py -b "aKF498dAp.3e8.1004" -s "aKF498dAp.4b0.7d0"

# Output:
# Chromosome: chr1 (from RefGet mapping)
# Strand: +
# Exon 1: 1,000-1,200
# Exon 2: 2,000-2,150
```

## 🧬 Use Cases

### Genomics Research
- **Cross-database integration**: Same isoform = same ID across Ensembl/RefSeq/GENCODE
- **Long-read sequencing**: Standardize novel isoform identification
- **Differential expression**: Compare isoform usage between samples
- **Splice variant analysis**: Track alternative splicing patterns

### Clinical Applications
- **Disease isoforms**: Identify disease-specific transcript variants
- **Biomarker discovery**: Find diagnostic isoform signatures
- **Drug target identification**: Map isoform-specific therapeutic targets

### Data Sharing
- **Universal compatibility**: Share data without chromosome naming conflicts
- **Reproducibility**: Same input = same output across labs
- **Database integration**: Link to external resources via RefGet IDs

## 🔬 Technical Details

### RefGet Chromosome Hashing

IsoTag solves chromosome naming inconsistencies (chr1 vs Chr1 vs CHR1 vs 1) by hashing **chromosome sequences** instead of names:

```python
# Step 1: Extract chromosome sequence
chr_sequence = "ACGTACGTACGT..."  # Full chromosome from FASTA

# Step 2: Calculate RefGet ID (GA4GH sha512t24u algorithm)
refget_id = sha512t24u(chr_sequence.encode('ascii'))
# Result: "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2" (32-char)

# Step 3: Use appropriate hash length
chr_hash_8 = refget_id[:8]   # "aKF498dA" - for XB/XS tags
chr_hash_32 = refget_id       # Full 32 chars - for XI/XT/XC/XV tags
```

**Result**: chr1, Chr1, CHR1, and 1 all map to same hash if they have identical sequences!

### Auto RefGet Cache System

First run with genome FASTA:
```bash
python3 isotag.py -i input.bam -o tagged.bam -g hg38.fa
# → Generates ~/.isotag_cache/hg38_refget.json automatically
```

Subsequent runs:
```bash
python3 isotag.py -i input2.bam -o tagged2.bam -g hg38.fa
# → Uses cached RefGet mapping (instant lookup, no FASTA parsing)
```

### Reversible Tag Encoding

**XB Tag (Boundary)**: `[8-chr-hash][strand].[5'-hex].[3'-hex]`
```
aKF498dAp.3e8.1004
│       │ │   │
│       │ │   └─ 3' end: 4100 (0x1004)
│       │ └───── 5' end: 1000 (0x3e8)
│       └─────── Strand: + (p=plus, m=minus)
└─────────────── Chromosome: aKF498dA (8-char RefGet hash)
```

**XS Tag (Splicetag)**: `[8-chr-hash][strand].[coord1].[coord2].[coord3]...`
```
aKF498dAp.4b0.7d0.866.bb8
│       │ │   │   │   │
│       │ └───┴───┴───┴─ Splice coordinates in hex
│       └─────────────── Strand: +
└───────────────────── Chromosome: aKF498dA
```

## 📊 Performance

- **Processing Speed**: ~1000 reads/second (tested on real data)
- **Memory Usage**: Streaming processing, constant memory
- **Storage Overhead**: Minimal (compact hex encoding)
- **Cache Generation**: One-time per genome (~30 seconds for hg38)
- **Tag Lookup**: Instant with cached RefGet mapping

## 🛠️ Advanced Options

### Clustering Parameters

```bash
# Custom position quantization (default: 10000bp)
python3 isotag.py -i input.bam -o tagged.bam -g genome.fa --position-quantum 5000

# Custom genomic span quantization (default: 10000bp)
python3 isotag.py -i input.bam -o tagged.bam -g genome.fa --span-quantum 5000

# Custom exon length quantization (default: 1000bp)
python3 isotag.py -i input.bam -o tagged.bam -g genome.fa --exon-quantum 500
```

### Output Control

```bash
# Quiet mode (minimal output)
python3 isotag.py -i input.bam -o tagged.bam -q

# Structure tags only (no variants, faster)
python3 isotag.py -i input.bam -o tagged.bam --no-variants

# Show progress updates
python3 isotag.py -i input.bam -o tagged.bam -g genome.fa
```

## 📁 File Formats

### RefGet Cache File

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
    "1": "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"
  }
}
```

### Tagged BAM Output

```
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422
read_001	0	chr1	1000	60	200M2000N150M	*	0	0	ATCG...	IIII...	XI:Z:fuIF7PN23g2gq9sFxqhUNGnfOCZhkQJS	XB:Z:aKF498dAp.3e8.1004	XS:Z:aKF498dAp.4b0.7d0	XT:Z:266CbPqmZz8eS-EzT4xtnYtmm-SoIhnL
```

## 🔍 Validation & Testing

```bash
# Test with sample data
samtools view -h input.bam | head -1000 | samtools view -b > sample.bam
python3 isotag.py -i sample.bam -o tagged_sample.bam -g genome.fa

# Verify tags added
samtools view tagged_sample.bam | head -5 | grep "XI:Z:"

# Decode random tag
samtools view tagged_sample.bam | grep "XS:Z:" | head -1 | \
  sed 's/.*XS:Z:\([^ ]*\).*/\1/' | \
  xargs -I {} python3 decode_tags.py -s {}

# Count unique isoforms
samtools view tagged_sample.bam | grep -o "XI:Z:[^ ]*" | sort -u | wc -l
```

## 🤝 Contributing

We welcome contributions! Areas for improvement:
- Additional clustering modes
- Performance optimizations
- Integration with annotation databases
- Support for fusion transcripts

See [CHANGELOG.md](CHANGELOG.md) for version history and [TAG_FORMAT.md](TAG_FORMAT.md) for detailed tag specifications.

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details.

## 🆘 Support & Citation

For issues and questions:
- GitHub Issues: https://github.com/LSBDT/isotag/issues
- Documentation: See [TAG_FORMAT.md](TAG_FORMAT.md) for technical details

If you use IsoTag in your research, please cite:
```
IsoTag: Universal Isoform Identification System using RefGet-compatible Identifiers
GitHub: https://github.com/LSBDT/isotag
Version: 2.3.0 (2026)
```

## 🔗 Related Resources

- **GA4GH RefGet Specification**: https://samtools.github.io/hts-specs/refget.html
- **VRS (Variation Representation Specification)**: https://vrs.ga4gh.org
- **SAM/BAM Format**: https://samtools.github.io/hts-specs/SAMv1.pdf

---

**Status**: ✅ Production Ready | **Version**: 2.3.0 | **Last Updated**: February 2026
