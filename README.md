# IsoTag v1.0 - RefGet-Compatible Universal Isoform Tagger

A RefGet-compatible, VRS-compliant universal isoform identification system for long-read transcript sequencing data. IsoTag v1.0 generates deterministic, globally unique identifiers using GA4GH RefGet SQ.XXXX sequence identifiers and official VRS algorithms, ensuring true universal compatibility across all genomic databases.

## Features

- **🌐 RefGet Compatible**: Uses GA4GH RefGet SQ.XXXX sequence identifiers for true universality
- **🧬 VRS Compliant**: Official GA4GH VRS algorithms ensure ecosystem integration
- **🔄 Cross-Database Universal**: Same sequence content = same RefGet ID across Ensembl, GENCODE, RefSeq, UCSC
- **✅ Specification Validated**: Passes official RefGet test vectors
- **⚡ Memory Efficient**: Stream-based processing using samtools for large BAM files
- **📁 Flexible Input**: Supports both BAM and SAM formats with automatic format detection
- **🏷️ Comprehensive Tagging**: Adds both RefGet-based structure (XI) and VRS-compliant variant (XV) tags

## Quick Start

### Installation

Requirements:
- Python 3.6+
- samtools (for BAM processing)
- click library

```bash
pip install click
```

### Basic Usage

**Important**: v1.0 requires either a reference genome FASTA or pre-computed RefGet mapping.

```bash
# With FASTA (calculates RefGet on-the-fly)
python3 isotag.py -i input.bam -o tagged.bam -g reference.fa

# With pre-computed RefGet mapping (recommended, faster)
python3 isotag.py -i input.bam -o tagged.bam -r hg38-refget.json

# Structure tags only (no variants)
python3 isotag.py -i input.bam -o tagged.bam -r hg38-refget.json --no-variants

# SAM format support
python3 isotag.py -i input.sam -o tagged.sam -g reference.fa
```

## Output Tags

IsoTag v1.0 adds two RefGet/VRS-compliant tags to BAM/SAM files:

- **XI:Z:** - 32-character RefGet-based structure hash (isoform structure ID)
- **XV:Z:** - VRS-compliant variant IDs separated by dots (when variants are present)

### Example Tags
```
XI:Z:fuIF7PN23g2gq9sFxqhUNGnfOCZhkQJS    # RefGet-based Structure ID
XV:Z:Q4fUfjJXgQwpSxFgeGVowhJaLTVg3Fqk.P8rTn3mKwQxYs5LcHgVjRkNpBzAd    # VRS-compliant Variant IDs
```

## Algorithm Overview

### RefGet-Based Structure ID Generation (v1.0)
1. Extract exon boundaries from CIGAR operations
2. Map chromosome name to RefGet SQ.XXXX identifier
3. Serialize as: `SQ.XXXX|strand|exon1_start:exon1_end|exon2_start:exon2_end|...`
4. Generate VRS-compatible hash using official `sha512t24u` algorithm (32 chars)

### VRS-Compliant Variant Detection (v1.0)
1. Parse MD tags for mismatches and deletions
2. Extract insertions from CIGAR operations
3. Create RefGet-based variants with SQ.XXXX identifiers
4. Generate individual VRS-compliant hashes using official GA4GH algorithms
5. Concatenate with dots for the XV tag

## Use Cases

- **Long-read single-cell analysis**: Consistent isoform identification across samples
- **Cross-database integration**: Map isoforms between Ensembl, RefSeq, GENCODE
- **Comparative transcriptomics**: Match isoforms across studies and conditions
- **Disease variant tracking**: Link specific variants to transcript structures

## Technical Details

### Variant Detection Modes

1. **Structure Only**: No MD tags, no reference genome
   - Only XI (structure) tags added
   - Fastest processing mode

2. **With Existing MD Tags**: Input contains MD tags
   - Both XI and XV tags added
   - Uses existing variant information

3. **Generate MD Tags**: Reference genome provided
   - Generates MD tags using `samtools calmd`
   - Both XI and XV tags added
   - Most comprehensive variant detection

### Memory Efficiency

- Stream-based processing avoids loading entire BAM into memory
- Uses samtools for efficient BAM/SAM handling
- Temporary files for format conversion only

## Command Line Options

```
Options:
  -i, --input TEXT      Input BAM/SAM file [required]
  -o, --output TEXT     Output BAM/SAM file with XI/XV tags [required]  
  -g, --genome TEXT     Reference genome FASTA (calculates RefGet IDs)
  -r, --refget TEXT     Pre-computed RefGet JSON mapping (recommended)
  --variants            Enable variant detection (default: True)
  --no-variants         Disable variant detection
  --help               Show this message and exit
```

**Note**: Exactly one of `--genome` or `--refget` must be specified in v1.0.

## Example Workflow

```bash
# 1. Generate RefGet mapping (one-time per genome)
python3 -c "from our_refget import RefGetCalculator; RefGetCalculator.calculate_from_fasta('genome.fa')"

# 2. RefGet-compatible isoform tagging
python3 isotag.py -i sample.bam -o sample_tagged.bam -g genome.fa

# 3. View RefGet-based results
samtools view sample_tagged.bam | grep 'XI:Z:' | head

# 4. With pre-computed RefGet mapping (faster)
python3 isotag.py -i sample.bam -o sample_tagged.bam -r hg38-refget.json

# 5. Extract RefGet-based isoform statistics
samtools view sample_tagged.bam | cut -f12- | grep 'XI:Z:' | sort | uniq -c
```

## Output Statistics

IsoTag v1.0 provides comprehensive processing statistics:

```
============================================================
✅ IsoTag v1.0 Complete - RefGet Compatible!
============================================================
📊 Total reads: 10,000
🧬 Reads processed: 9,847
🏷️ RefGet structure tags (XI): 9,847
🧪 VRS variant tags (XV): 8,123
🆔 Unique structures: 1,234
🔬 Unique variant combinations: 856
💾 Output: sample_tagged.bam

🎯 Example RefGet-based tags:
   XI:Z:fuIF7PN23g2gq9sFxqhUNGnfOCZhkQJS
   XV:Z:Q4fUfjJXgQwpSxFgeGVowhJaLTVg3Fqk

🌐 Universal Compatibility: RefGet SQ.XXXX + VRS compliant
```

## File Structure

The repository contains the following core files:

- **`isotag.py`** - Main v1.0 script (RefGet-compatible, VRS-compliant)
- **`isotag_refget.py`** - RefGet SQ.XXXX calculator module
- **`vrs_compat.py`** - VRS-compliant variant handling module
- **`README.md`** - This documentation
- **`LICENSE`** - MIT license

## GA4GH Standards Compliance

IsoTag v1.0 is fully compliant with:

- **RefGet**: Uses official SQ.XXXX sequence identifiers
- **VRS**: Official `sha512t24u` algorithm from GA4GH VRS-Python
- **Cross-database universal**: Works with Ensembl, GENCODE, RefSeq, UCSC

## Contributing

This is Phase 1 of the IsoTag system focusing on core RefGet/VRS compatibility. Future enhancements may include:

- Integration with annotation databases
- REST API for ID lookup services  
- Support for fusion transcripts
- Advanced variant normalization

## License

MIT License - see LICENSE file for details.

## Citation

If you use IsoTag v1.0 in your research, please cite:

```
IsoTag v1.0: RefGet-Compatible Universal Isoform Tagger
LSBDT Team (2025)
GitHub: https://github.com/LSBDT/isotag
```

## Contact

For questions and support, please open an issue on GitHub or contact the LSBDT team.