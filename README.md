# IsoTag v2.4 - Universal Isoform Identification System 🚀

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![GA4GH](https://img.shields.io/badge/GA4GH-RefGet--inspired-blue.svg)](https://samtools.github.io/hts-specs/refget.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18676569.svg)](https://doi.org/10.5281/zenodo.18676569)

A production-ready system for generating standardized, RefGet-compatible isoform identifiers for long-read transcript sequencing data. **v2.4** adds TSS (X5) and PolyA (X3) cluster tags, rewrites the XC locus tag with validated isoform-length encoding, and ships the full 38-tool analysis toolkit.

## 🎯 What's New in v2.4

### **TSS Cluster Tag (X5)** 🧬 NEW!
- **Binned 5' end**: Groups reads by transcription start site (default 200bp bins; use `--x5x3-bin-size 35` for higher resolution)
- **Validated against FANTOM5 CAGE**: 54.9% of X5 clusters (35bp bins) overlap a FANTOM5 TSS within ±100bp; 45.3% with matched ±35bp window — **43× above genic-space null** (Z=1903, 1000 permutations within GENCODE v47 gene bodies)
- **Use case**: Group reads by promoter identity; TSS-polyA coupling analysis with X3

### **PolyA Cluster Tag (X3)** 🧬 NEW!
- **Binned 3' end**: Groups reads by polyadenylation site (default 200bp bins)
- **Validated against PolyA_DB v4**: 29.16% of X3 clusters overlap a known polyA site within ±100bp — **14× above genic-space null** (Z=811, 1000 permutations)
- **Use case**: Characterize alternative polyadenylation; coupled TSS/polyA analysis across the same read

### **X5/X3 Validation Summary**

| Tag | Bin size | Reference | Overlap | Null (genic-space) | Enrichment |
|-----|----------|-----------|---------|-------------------|-----------|
| X5 | 200bp | FANTOM5 CAGE (209,911 peaks) | 43.3% | 1.27% | 34× |
| X5 | 35bp | FANTOM5 CAGE | **54.9%** (±100bp); 45.3% (±35bp) | 1.27% | **43×** |
| X3 | 200bp | PolyA_DB v4 (303,312 sites) | **29.16%** | 2.06% | **14×** |
| X3 | 200bp | PolyASite v3.0 high-conf (88,284 sites) | 18.5% | — | — |

### **XC Tag v11.0 Rewrite** 🔄 (replaces v10.0)
- **Algorithm**: `chr_hash | strand | middle_bin(1kb) | exon_len_bin1(10bp) | ...` — midpoint + binned exon lengths
- **Isoform-sensitive**: Reads from different isoforms (different exon structure) receive different XC tags; terminal truncations usually share the same XC
- **Validated**: 76.0% multi-read cluster purity; ARI=0.4238; V-measure H=0.913 vs GENCODE v47 gene annotation (MPC human dataset, hg38)
- ⚠️ **Wobble caveat**: 10bp length bins tolerate ±1bp for single-junction reads (~80% same tag), but 7-junction reads fail ~79% of the time at any junction. Use XC for locus-level grouping, not exact isoform identity — use **XI** for exact matches.
- **Breaking change from v10.0**: Re-tag with `isotag.py` — v10.0 and v11.0 hashes are incompatible

### **Full Toolkit (38 tools)** 🛠️ NEW!
All 38 analysis scripts now included in this release. See [Full Toolkit Reference](#-full-toolkit-reference) below.

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

### **Locus Tag (XC Tag v10.0)** *(superseded by v11.0 in v2.4)*
- ~~**Spatial bin tag**: Groups transcripts at the same genomic location using 10kb position bins~~
- **v10.0 is replaced by XC v11.0** in v2.4: the location-only algorithm was insufficient; see v2.4 section above for the validated midpoint + exon-length approach
- Tags generated with v10.0 are incompatible with v11.0; re-tag using `isotag.py` from v2.4

## 🎯 What's New in v2.1

### **Unified Genome Build Support** 🔧 NEW!
- **Ambiguous base masking**: Automatically converts ambiguous IUPAC codes (R, Y, S, W, K, M, etc.) to 'N'
- **Cross-build compatibility**: UCSC hg38 and NCBI GRCh38 now produce identical chromosome hashes
- **Informative output**: Reports ambiguous bases masked per chromosome during RefGet cache generation
- **Optional override**: `--keep-ambiguous-bases` flag available if needed (not recommended)

> **⚠️ Migration warning for v2.2.0 users:** RefGet JSON files generated with IsoTag v2.0–v2.2.0 are **incompatible** with v2.2.1+. The ambiguous base masking introduced in v2.2.1 changed how chromosome hashes are computed. Any `.json` file from v2.0–v2.2.0 must be regenerated:
> ```bash
> python3 isotag_refget.py -f genome.fa -o genome-refget.json -g genome_name
> ```
> Using an old JSON with a new `isotag.py` will silently produce wrong XI/XB/XS/XT/XC tags (hashes will not match between datasets).

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
- **Exact splice junction matching**: groups reads sharing identical junction coordinates into the same transcript group
- **Not wobble-tolerant**: ±1bp Nanopore junction wobble changes the XT hash; use **XC** for approximate locus-level clustering
- **Sample-independent**: consistent grouping across experiments

### **Universal Chromosome Hashing** 🌍
- **RefGet-based**: Uses GA4GH sequence identifiers
- **Auto-caching**: Generates RefGet cache from genome FASTA automatically
- **Cross-database**: Same chromosome sequence = same hash regardless of naming

## 📋 Tag Format (v2.4.0)

```bash
XI:Z:fuIF7PN23g2gq9sFxqhUNGnfOCZhkQJS              # Structure ID (32-char)
XB:Z:aKF498dAp.3e8.1004                           # Boundary tag (8-char chr + hex ends)
XS:Z:aKF498dAp.4b0.7d0.866.bb8                    # Splicetag (8-char chr + hex coords)
XT:Z:266CbPqmZz8eS-EzT4xtnYtmm-SoIhnL              # Transcript group (32-char, exact junctions)
XC:Z:yQ61R9bjPTau4Rm13SDg3D0MVr_EIos3              # Cluster ID (midpoint + exon lengths, v11.0)
X5:Z:H6hLexXgCSsl8SDqLZs8zbJzUaaqTNif              # TSS cluster (binned 5' end, 200bp default)
X3:Z:9ho2wDAU3-1r3o8SPNK_5FhTlbWGKIQC              # PolyA cluster (binned 3' end, 200bp default)
XV:Z:ga4gh:VA.gbvJw0s4OeAvloCeAM6BNNvrjFC_Dhc8     # GA4GH VRS v2 Allele ID (comma-separated, optional; v2.4+)
```

### Tag Breakdown

| Tag | Name | Format | Purpose |
|-----|------|--------|---------|
| **XI** | Structure ID | 32-char hash | Unique isoform structure identifier |
| **XB** | Boundary Tag | `[8-chr][s].[hex1].[hex2]` | Reversible 5'/3' transcript ends |
| **XS** | Splicetag | `[8-chr][s].[hex1].[hex2]...` | Reversible splice junction coordinates |
| **XT** | Transcript Group | 32-char hash | Biological clustering (exact splice junctions; not wobble-tolerant) |
| **XC** | Cluster ID | 32-char hash | Locus cluster (midpoint 1kb + exon lengths 10bp bins; v11.0) |
| **X5** | TSS Cluster | 32-char hash | TSS cluster by binned 5' end (default 200bp) |
| **X3** | PolyA Cluster | 32-char hash | PolyA cluster by binned 3' end (default 200bp) |
| **XV** | Variants | `ga4gh:VA.<32-char>`, comma-separated | GA4GH VRS v2 Allele IDs, SNVs/substitutions only (optional; v2.4+) |

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

**Requirements**: Python 3.8+, samtools 1.18+, click, pysam

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

## 📦 Example Datasets

Pre-tagged BAM files are available on Zenodo for both hg38 and T2Tv2 genomes, annotated against GENCODE v49 and RefSeq:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18676569.svg)](https://doi.org/10.5281/zenodo.18676569)

**Download**: https://doi.org/10.5281/zenodo.18676569

| File | Genome | Annotation | Size |
|------|--------|------------|------|
| `hg38_GENCODEv49_sorted_isotagged.bam` | hg38 | GENCODE v49 | 280 MB |
| `hg38_RefSeq_sorted_isotagged.bam` | hg38 | RefSeq | 249 MB |
| `T2Tv2_GENCODEv49_sorted_isotagged.bam` | T2Tv2 | GENCODE v49 | 304 MB |
| `T2Tv2_RefSeq_sorted_isotagged.bam` | T2Tv2 | RefSeq | 241 MB |

All files contain the full IsoTag set (XI, XB, XS, XT, XC, XV tags) and can be used as benchmarks for validating your own tagged BAMs.

> **Note:** Zenodo BAMs use GENCODE v49 annotations. The validation metrics reported in this README (purity, ARI, V-measure) were computed on an internal dataset against GENCODE v47. These metrics are not directly comparable to results derived from the Zenodo BAMs.

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

## Background

Long-read RNA sequencing has enabled the detection of tens of thousands of novel transcripts, yet standardizing this molecular variability across tools, datasets, and databases remains an open challenge. As noted by Monzó, Frankish & Conesa (*Genome Research*, 2025), the community must "agree on the strategy to capture molecular variability while still defining reference annotations that are useful for the genomics community" ([doi:10.1101/gr.279865.124](https://doi.org/10.1101/gr.279865.124)).

IsoTag addresses this by providing read-level deterministic tags written directly into BAM files. Unlike transcript assembly tools (IsoQuant, StringTie2) or quantification tools (Bambu, oarfish, TranSigner) — which output GTF files or expression matrices — IsoTag encodes transcript identity as permanent, portable BAM tags. The result is an analysis-ready BAM that any downstream tool can consume without re-running alignment or reassembly.

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

IsoTag solves chromosome naming inconsistencies (chr1 vs Chr1 vs CHR1 vs 1) by hashing **chromosome sequences** instead of names, using the GA4GH `sha512t24u` algorithm.

**One important difference from canonical GA4GH RefGet:** IsoTag masks IUPAC ambiguous bases (R, Y, S, W, K, M, etc.) to `N` before hashing. This means IsoTag hashes will **not match** canonical GA4GH RefGet server IDs. The reason: UCSC hg38 and NCBI GRCh38 encode some positions differently; masking to `N` first ensures both assemblies produce identical tags.

```python
# Step 1: Extract and normalize chromosome sequence
chr_sequence = "ACGTACGTACGT..."  # Full chromosome from FASTA
normalized = re.sub(r'[RYSWKMBDHV]', 'N', chr_sequence.upper())  # mask ambiguous bases

# Step 2: Hash with GA4GH sha512t24u algorithm
refget_id = sha512t24u(normalized.encode('ascii'))
# Result: "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2" (32-char)

# Step 3: Use appropriate hash length
chr_hash_8 = refget_id[:8]   # "aKF498dA" - for XB/XS tags
chr_hash_32 = refget_id       # Full 32 chars - for XI/XT/XC/XV tags
```

**Result**: chr1, Chr1, CHR1, and 1 all map to the same hash if they have the same underlying sequence.

> **⚠️ Hash collision probability**: IsoTag uses 24-byte (192-bit) truncated SHA-512. The probability of any two different sequences producing the same hash is <10⁻⁴⁶ — negligible for any realistic transcript dataset.

### XC Tag Limitations (v11.0)

The XC locus tag groups transcripts by genomic midpoint (1kb bins) + binned exon lengths (10bp bins). Key limitations:

- **Hard bin boundaries**: A gene spanning a bin boundary produces two different XC tags for reads at each end. Two reads from the same gene can receive different XC tags if their midpoints fall in adjacent 1kb bins.
- **Wobble tolerance degrades with junction count**: Single-junction reads tolerate ±1bp (80% same tag); 7-junction reads have 79% probability of bin boundary cross. Median human transcript has 6–8 junctions.
- **XC is a locus/isoform tag, not a gene ID.** Empirical validation (hg38, MPC dataset, GENCODE v47): multi-read cluster purity = **76.0%**; false merge rate = **24.0%**; ARI vs gene annotation = 0.4238. Both XC v11.0 and midpoint-only XC (73.5%) remain below 1kb genomic background (80.4%) — enrichment 0.946×.
- **Alternative 3' ends**: Transcripts with different polyadenylation sites in different bins get different XC tags even if they share the same 5' end.

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

- **Processing Speed**: ~1,000 reads/second (single-threaded Python; tested on a 2.8 GHz Intel Xeon with a 9,968-read ONT Nanopore BAM; throughput will vary with read length, BAM compression, and I/O speed)
- **Memory Usage**: Streaming processing — constant memory regardless of BAM size; no full-file loading
- **Storage Overhead**: Minimal (compact hex encoding adds ~100–200 bytes per read)
- **Cache Generation**: One-time per genome (~30 seconds for hg38 on a standard server)
- **Tag Lookup**: Instant with pre-built RefGet JSON (no FASTA parsing on subsequent runs)

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

## 🧰 Full Toolkit Reference

IsoTag ships 38 tools for end-to-end transcript analysis:

| Tool | Purpose |
|------|---------|
| `isotag.py` | **Core** — tag BAM with all 8 IsoTag tags (XI/XB/XS/XT/XV/XC/X5/X3) |
| `isotag_refget.py` | Build and merge RefGet JSON chromosome maps |
| `isotag_create_bed_index.py` | BED index — 100–5000× smaller than full BAM, 4–5× faster queries |
| `isotag_query_bed_index.py` | Query BED index for novel/matching isoforms |
| `isotag_intersect.py` | Find reads matching tags in a reference BAM |
| `decode_tags.py` | Decode XB/XS reversible tags back to coordinates |
| `decode_splicetag.py` | Decode XS splicetag to junction coordinates |
| `isotag_stats.py` | Comprehensive per-tag statistics and distribution analysis |
| `isotag_filter.py` | Filter BAM by specific isotag IDs |
| `isotag_bam_filter.py` | Filter BAM by tag presence/absence/value criteria |
| `isotag_count.py` | Count reads per tag value or cluster |
| `isotag_subset.py` | Subset BAM to reads matching a tag criterion |
| `isotag_extract.py` | Extract reads by tag criteria to BED/TSV |
| `isotag_annotate.py` | Annotate reads with gene/transcript information |
| `isotag_validate.py` | Validate tag format and internal consistency |
| `isotag_merge.py` | Merge two tagged BAMs (union or intersection of tag sets) |
| `isotag_diff.py` | Differential tag analysis between two BAMs |
| `isotag_compare.py` | Pairwise experiment comparison (shared/unique isoforms) |
| `isotag_clustering.py` | Cluster analysis: purity, entropy, cluster size distribution |
| `isotag_cluster_representative.py` | Select representative read per XC/XI cluster |
| `isotag_coverage.py` | Per-tag coverage profiles across genomic regions |
| `isotag_index.py` | Build tag-to-read index for fast lookup |
| `isotag_query.py` | Query reads by tag value |
| `isotag_query_index.py` | Query reads via pre-built index |
| `isotag_novel_ranking.py` | Rank novel isoforms by evidence and uniqueness |
| `isotag_isoform_variants.py` | Characterize variant patterns per isoform |
| `isotag_boundary_analysis.py` | 5'/3' end heterogeneity across isoform clusters |
| `isotag_fuzzy_rescue.py` | Rescue reads with ±1bp wobble across bin boundaries |
| `isotag_convert.py` | Convert IsoTag output to GTF/GFF/BED formats |
| `isotag_stats_compare.py` | Compare tag statistics across two experiments |
| `isotag_xc_validate.py` | Validate XC clusters vs GENCODE (ARI/V-measure/purity) |
| `isotag_xc_gencode_compare.py` | XC version stability across GENCODE releases |
| `inspect_bam.py` | Inspect BAM header, tag presence, and alignment stats |
| `inspect_tags.py` | Inspect tag values and decode summaries |
| `tag_definitions.py` | Shared tag schema and hash utilities (library) |
| `tag_summary.py` | Summarize all tag types in a BAM to TSV |
| `isotag_utils.py` | Shared utilities: masking, hashing, encoding (library) |
| `vrs_compat.py` | GA4GH VRS-compatible sequence identifier utilities |

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
Version: 2.4.0 (2026)
```

## 🔗 Downstream Tools

IsoTag-tagged BAM files are compatible with existing long-read RNA analysis tools:

- **[TranSigner](https://github.com/Loren1994/TranSigner)** — isoform quantification using BAM-native tags
- **[oarfish](https://github.com/COMBINE-lab/oarfish)** — Bayesian isoform quantification for long-read RNA-seq
- **[scywalker](https://github.com/jksr/scywalker)** — single-cell long-read isoform analysis; XC tags can substitute for UMI-based cell grouping

## 🔗 Related Resources

- **IsoTag Example Datasets (Zenodo)**: https://doi.org/10.5281/zenodo.18676569
- **GA4GH RefGet Specification**: https://samtools.github.io/hts-specs/refget.html
- **GA4GH SeqCol Specification v1.0**: https://ga4gh.github.io/refget/seqcols/ (bioRxiv: https://doi.org/10.1101/2025.10.06.680641)
- **FANTOM5 CAGE Atlas**: https://fantom.gsc.riken.jp/5/
- **PolyA_DB v4**: https://exon.njms.rutgers.edu/polya_db/v4/
- **VRS (Variation Representation Specification)**: https://vrs.ga4gh.org
- **SAM/BAM Format**: https://samtools.github.io/hts-specs/SAMv1.pdf
- **Poster AGBT General Meeting 2026**: https://zenodo.org/records/18653366
- **Monzó et al. 2025** (long-read RNA standardization rationale): https://doi.org/10.1101/gr.279865.124

---

**Status**: ✅ Production Ready | **Version**: 2.4.0 | **Last Updated**: June 2026
