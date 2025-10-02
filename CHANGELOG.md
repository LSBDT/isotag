# Changelog

All notable changes to IsoTag will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
