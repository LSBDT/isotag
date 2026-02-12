# Changelog

All notable changes to IsoTag will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.3.0] - 2026-02-12

### 🔀 Minor Release: Merged RefGet JSONs + NC_ Alias Bug Fix

This release merges UCSC and NCBI RefGet JSON files into unified files and fixes a bug where NCBI accession numbers were incorrectly converted to chromosome aliases.

### Added

#### Merged RefGet JSONs
- **Unified genome files**: UCSC and NCBI versions combined into single JSON per genome
- **5 merged files** replace previous 10 separate files:
  - `hg38.GRCh38-refget.json` (2,525 mappings)
  - `hg19.GRCh37-refget.json` (669 mappings)
  - `hs1.T2T-CHM13v2-refget.json` (124 mappings)
  - `mm39.GRCm39-refget.json` (305 mappings)
  - `mm10.GRCm38-refget.json` (503 mappings)
- **All naming conventions in one file**: `chr1`, `1`, `CHR1`, `Chr1`, `NC_000067.7`

#### Merge Tool
- **New `--merge` option** for `isotag_refget.py`: combine multiple RefGet JSONs
- Usage: `python3 isotag_refget.py -m file1.json -m file2.json -o merged.json -g name`
- Automatic cleanup of bogus aliases during merge
- Conflict detection for mismatched RefGet IDs

### Fixed

#### NC_ Accession Alias Bug 🐛
- **Bug**: NCBI RefSeq accession numbers (e.g., `NC_000067.7`) were incorrectly converted to chromosome aliases by extracting the numeric portion
- **Impact**: Mouse chr1 (`NC_000067.7`) generated bogus aliases `67`, `chr67`, `Chr67`, `CHR67`
- **Affected genomes**: GRCm39, GRCm38 (88 bogus each), GRCh38, GRCh37 (4 bogus each), T2T-CHM13v2 (96 bogus)
- **Total**: 280 bogus aliases removed across all genomes
- **Fix**: NC_ accessions are now kept as-is without generating chromosome number aliases

### Migration

- Replace individual RefGet JSONs with merged versions:
  - `hg38-refget.json` or `GRCh38-refget.json` → `hg38.GRCh38-refget.json`
  - `mm39-refget.json` or `GRCm39-refget.json` → `mm39.GRCm39-refget.json`
- Update `--refget` paths in your scripts accordingly

---

## [2.2.0] - 2026-02-09

### 🧬 Minor Release: Gene/Locus Tag (XC)

This release adds the **XC tag** - a pure location-based gene/locus identifier that groups all transcripts at the same genomic location regardless of isoform structure. Ideal for gene-level analysis without requiring annotation databases.

### Added

#### XC Tag - Gene/Locus Cluster ID
- **New tag**: `XC:Z:` - 32-character hash identifying genomic locus
- **Pure location-based**: Uses only chromosome, strand, and binned start/end positions
- **Gene-level grouping**: All isoforms at the same locus get the same XC, regardless of alternative splicing
- **No annotation required**: Automatically assigns gene-level IDs to novel and known transcripts
- **Configurable resolution**: `--xc-bin-size` parameter (default: 10000bp = 10kb)

#### XC vs Existing Tags

| Feature | XI (Structure) | XT (Transcript Group) | **XC (Gene/Locus)** |
|---------|---------------|----------------------|---------------------|
| Exon coordinates | ✓ | Quantized | **✗** |
| Splice junctions | ✓ | ✓ | **✗** |
| Exon count/lengths | ✓ | Quantized | **✗** |
| Position | Exact | Binned | **Binned** |
| Strand | ✓ | ✓ | **✓** |
| Same gene, different isoforms | Different | May differ | **Same** ✅ |
| Use case | Exact isoform | Isoform clustering | **Gene-level ID** |

#### CLI Options
- **`--xc-bin-size`**: Bin size for XC clustering (default: 10000bp)
  - 10kb (default): Gene-level analysis
  - 100bp: Isoform discovery with tolerance
  - 10bp: Splice wobble detection

### Unchanged
- All 5 existing tags (XI, XB, XS, XT, XV) are fully backward compatible
- `decode_tags.py`, `isotag_refget.py`, `vrs_compat.py` unchanged (XC is a non-reversible hash)
- RefGet cache format unchanged

### Technical Details

```bash
# XC tag format
XC:Z:a7Bf9xK2mP3qR5tN8wY1zC4dF6hJ0lO    # 32-char hash (gene/locus ID)

# XC serialization (before hashing)
# chr_hash_32|strand|start_bin|end_bin
"aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2|+|100|102"

# Same-location isoforms → Same XC
Isoform A (5 exons, 1M-1.02M): XC = abc123...  ← SAME
Isoform B (3 exons, 1M-1.02M): XC = abc123...  ← SAME
Isoform C (4 exons, 1M-1.02M): XC = abc123...  ← SAME
```

### Usage

```bash
# Standard tagging (XC included by default)
python3 isotag.py -i input.bam -o tagged.bam -g genome.fa

# Fine-grained XC clustering (100bp bins)
python3 isotag.py -i input.bam -o tagged.bam -g genome.fa --xc-bin-size 100

# Splice wobble detection (10bp bins)
python3 isotag.py -i input.bam -o tagged.bam -g genome.fa --xc-bin-size 10
```

---

## [2.1.0] - 2025-10-23

### 🔧 Minor Release: Unified Genome Build Support

This release adds ambiguous base masking to ensure that UCSC hg38 and NCBI GRCh38 (and other genome build variants) produce identical RefGet chromosome hashes, improving cross-database compatibility.

### Added

#### Ambiguous Base Masking
- **Automatic masking**: Ambiguous IUPAC nucleotide codes (R, Y, S, W, K, M, etc.) are converted to 'N' before RefGet hashing
- **Cross-build compatibility**: UCSC hg38 and NCBI GRCh38 now produce identical chromosome hashes
- **Informative output**: Reports count of ambiguous bases masked per chromosome
- **New CLI option**: `--keep-ambiguous-bases` flag to preserve old behavior if needed (not recommended)
- **Cache metadata**: RefGet cache now includes `ambiguous_bases_masked` flag

#### Benefits
- Universal isotags work across genome build variants (hg38, GRCh38, etc.)
- Eliminates hash mismatches caused by different ambiguous base representations
- Maintains backward compatibility (existing caches will be regenerated with new behavior)

### Changed

#### RefGet Cache Generation
- Sequences are normalized (ambiguous bases → 'N') before hashing by default
- Enhanced progress output shows ambiguous base masking statistics
- Cache files include metadata about masking status

### Technical Details

```bash
# Example: UCSC hg38 chr1 has ambiguous bases (R, Y, etc.)
# NCBI GRCh38 chr1 may have different ambiguous base representation

# v2.0.0: Different hashes for UCSC vs NCBI
UCSC hg38  chr1 → SQ.aKF498dA... (includes R, Y, etc.)
NCBI GRCh38 chr1 → SQ.bXG723fB... (different ambiguous bases)

# v2.1.0: Identical hashes (ambiguous bases → N)
UCSC hg38  chr1 → SQ.aKF498dA... (R,Y,etc → N)
NCBI GRCh38 chr1 → SQ.aKF498dA... (R,Y,etc → N)  ✅ SAME!
```

### Migration

No action required! The change is backward compatible:
- Existing RefGet caches will be regenerated automatically on next use
- Isotags will remain consistent within your analysis
- For cross-database comparisons, regenerate RefGet cache with v2.1.0

### Performance

No performance impact - masking is done during cache generation only (one-time operation).

---

## [2.0.0] - 2025-10-02

### 🚀 Major Release: Reversible Tags & Universal Chromosome Hashing

This is a **major version** release introducing reversible tag encoding, universal chromosome hashing, and biological clustering modes. These features enable full coordinate reconstruction without reference genomes and solve cross-database compatibility issues.

### Added

#### Reversible Splicetag System (XS Tag)
- **XS Tag**: Reversible splicetag with 8-char chromosome hash + hex-encoded coordinates
- Format: `[8-chr-hash][strand].[hex1].[hex2].[hex3]...`
- Full splice junction reconstruction from tag alone (no reference needed)
- Compact hex encoding for minimal storage overhead
- Example: `XS:Z:aKF498dAp.4b0.7d0.866.bb8`

#### Boundary Tag System (XB Tag)
- **XB Tag**: Reversible boundary tag encoding transcript 5'/3' ends
- Format: `[8-chr-hash][strand].[5'-hex].[3'-hex]`
- Captures full transcript span in compact format
- Cross-validation: XB + XS = complete exon structure verification
- Example: `XB:Z:aKF498dAp.3e8.1004`

#### Biological Clustering (XT Tag Enhancement)
- **Three clustering modes**: 5prime, middle, 3prime
  - `5prime`: CAGE/TSS data (clusters by transcription start site)
  - `middle`: RNA-seq data (clusters by middle position)
  - `3prime`: PolyA/TES data (clusters by transcription end site)
- Fuzzy boundary handling for variable transcript ends
- Sample-independent clustering (same groups across experiments)
- Configurable quantization parameters

#### Universal Chromosome Hashing
- **RefGet-based chromosome hashing**: Hash sequences instead of names
- **Auto RefGet cache system**: Generates `~/.isotag_cache/[genome]_refget.json` automatically
- **Cross-database compatibility**: chr1/Chr1/CHR1/1 all map to same hash
- **8-char chromosome hashes**: Used in XB and XS tags for compact encoding
- **32-char chromosome hashes**: Used in XI, XT, XV tags for full hashing

#### Tag Decoder Utility
- **decode_tags.py**: Comprehensive tag decoding utility
- Decode XB boundary tags to get 5'/3' ends
- Decode XS splicetags to get splice junction coordinates
- Reconstruct full exon structure from XB + XS
- Chromosome name lookup via RefGet mapping
- Command-line interface with multiple options

### Changed

#### Core Tag System
- **XI Tag**: Now uses 32-char RefGet chromosome hash in serialization
- **XT Tag**: Enhanced with mode-based clustering (5prime/middle/3prime)
- **XV Tag**: Now uses 32-char RefGet chromosome hash for variants

#### Performance & Caching
- First run with genome FASTA: Generates RefGet cache automatically
- Subsequent runs: Uses cached RefGet mapping (instant lookup)
- Cache location: `~/.isotag_cache/[genome_name]_refget.json`
- No FASTA parsing needed after initial cache generation

#### Command-Line Interface
- Added `--clustermode` option: Choose clustering mode (5prime/middle/3prime)
- Added `--position-quantum`: Custom position quantization (default: 10000bp)
- Added `--span-quantum`: Custom genomic span quantization (default: 10000bp)
- Added `--exon-quantum`: Custom exon length quantization (default: 1000bp)

### Technical Details

#### Tag Format Changes
```bash
# v1.0 Format
XI:Z:Uy3v73FzY4ZhB-w0ZLsDwQLJfMl34Hzx    # Structure ID only
XV:Z:variant1.variant2                   # Variant IDs

# v2.0 Format (NEW)
XI:Z:fuIF7PN23g2gq9sFxqhUNGnfOCZhkQJS    # Structure ID (RefGet chr hash)
XB:Z:aKF498dAp.3e8.1004                  # Boundary tag ✨ NEW
XS:Z:aKF498dAp.4b0.7d0.866.bb8           # Splicetag ✨ NEW
XT:Z:266CbPqmZz8eS-EzT4xtnYtmm-SoIhnL    # Transcript group (mode-based) ✨ ENHANCED
XV:Z:Q4fUfjJXgQwpSxFgeGVowhJaLTVg3Fqk    # Variants (RefGet chr hash)
```

#### Serialization Format Changes
- **v1.0**: `chr|strand|exon1_start:exon1_end|...`
- **v2.0**: `[32-chr-refget-hash]|strand|exon1_start:exon1_end|...`

#### Backward Compatibility
- ⚠️ **Breaking Change**: Tag serialization format changed (includes RefGet hashes)
- v1.0 tags are **not compatible** with v2.0 (different hashing)
- Re-tagging required for existing BAM files
- Migration: Simply re-run isotag.py v2.0 on original BAM files

### Migration Guide (v1.0 → v2.0)

#### Step 1: Update Dependencies
```bash
# Ensure you have the latest version
git pull origin main
pip install --upgrade click
```

#### Step 2: Re-tag Existing BAM Files
```bash
# v1.0 tagged BAM
old_tagged.bam (with v1.0 XI/XV tags)

# Re-tag with v2.0
python3 isotag.py -i original.bam -o new_tagged.bam -g genome.fa

# Result: v2.0 tags (XI/XB/XS/XT/XV)
```

#### Step 3: Update Analysis Pipelines
```bash
# Old workflow (v1.0)
samtools view tagged.bam | grep "XI:Z:"

# New workflow (v2.0) - same commands still work!
samtools view tagged.bam | grep "XI:Z:"

# New capabilities (v2.0)
python3 decode_tags.py -b "XB_tag_value" -s "XS_tag_value" --reconstruct
```

### Performance Improvements
- **Cache-based processing**: RefGet cache eliminates repeated FASTA parsing
- **Compact hex encoding**: XB/XS tags use minimal storage
- **Streaming processing**: Constant memory usage regardless of BAM size
- **Parallel-ready**: Chromosome-based processing enables parallelization

### Documentation
- Comprehensive README.md update with v2.0 features
- TAG_FORMAT.md specification document
- CHANGELOG.md (this file)
- Detailed examples for all clustering modes
- Migration guide from v1.0 to v2.0

---

## [1.0.0] - 2025-08-20

### 🎉 Initial Public Release

First public release of IsoTag on GitHub (https://github.com/LSBDT/isotag).

### Added

#### Core Functionality
- **XI Tag**: Structure ID using 32-char VRS-compatible hash
- **XV Tag**: Individual variant IDs with dot-separated concatenation
- RefGet-compatible identifier generation
- VRS-compliant sha512t24u algorithm
- BAM/SAM format support (input and output)

#### Features
- RefGet SQ.XXXX integration for universal compatibility
- VRS-compliant VARSID for variant identification
- Cross-database universal identifiers (Ensembl, GENCODE, RefSeq, UCSC)
- Pre-computed RefGet tables support (JSON format)
- Memory-efficient streaming processing
- Individual variant tracking with simple chr:pos:ref>alt format

#### Tools & Utilities
- `isotag.py`: Main tagging script
- `isotag_refget.py`: RefGet SQ.XXXX calculator
- `vrs_compat.py`: VRS-compliant algorithms

#### Documentation
- README.md with quick start guide
- Installation instructions
- Basic usage examples
- RefGet specification validation

### Technical Details

#### Tag Format (v1.0)
```bash
XI:Z:Uy3v73FzY4ZhB-w0ZLsDwQLJfMl34Hzx              # Structure ID
XV:Z:85IT2XfuQ7nofo_xMFc4PzCZNtgAcZ96              # Variant IDs
```

#### Serialization Format (v1.0)
```
Structure: chr|strand|exon1_start:exon1_end|exon2_start:exon2_end|...
Variants:  chr:pos:ref>alt
```

#### Performance (v1.0)
- Processing speed: ~1000 reads/second
- Memory usage: Constant (streaming)
- File size overhead: Minimal

### Validation
- GA4GH RefGet specification test vectors: ✅ PASSED
- Real-world testing: mm39 mouse genome (22 chromosomes, 2.8GB)
- Production testing: SRR10540253.bam (1.4GB mouse RNA-seq)

---

## Version Numbering Scheme

- **Major version** (X.0.0): Breaking changes, new core features, incompatible with previous versions
- **Minor version** (0.X.0): New features, backward compatible
- **Patch version** (0.0.X): Bug fixes, documentation updates

## Links

- **Repository**: https://github.com/LSBDT/isotag
- **Issues**: https://github.com/LSBDT/isotag/issues
- **Releases**: https://github.com/LSBDT/isotag/releases

## Contributors

IsoTag is developed and maintained by the LSBDT team.

---

**Legend**:
- ✨ New feature
- 🔧 Enhancement
- 🐛 Bug fix
- ⚠️ Breaking change
- 📝 Documentation
