#!/usr/bin/env python3
"""
isotag_fuzzy_rescue - Find Reads Recovered by Fuzzy Matching

Identifies reads that match between A and B using XS (fuzzy splice junctions ±5bp)
but do NOT match using XI (exact structure). These represent sequencing errors
or minor variants that fuzzy matching successfully rescued.

This tool uses OFFSET-BASED comparison to correctly identify rescued reads.

Usage:
    python3 isotag_fuzzy_rescue.py --index matches.db --output rescued.tsv

    # With BAM extraction
    python3 isotag_fuzzy_rescue.py --index matches.db --output rescued.tsv \
        --bam-a experiment.bam --bam-out rescued.bam

Related Tools:
    - XC tag (isotag.py v9.0+): Provides location-based clustering without splice
      precision. Use --xc-bin-size to control clustering granularity.
      XC groups transcripts with similar location+exon structure, ignoring
      precise splice junction coordinates.

    - This tool (fuzzy_rescue): Compares XI vs XS matches in intersection database.
      Use to evaluate the benefit of fuzzy matching on your data.

Author: IsoTag Team
Version: 2.1.0 (Updated with XC reference)
"""

import sqlite3
import pickle
import pysam
import os
import sys
import argparse
from collections import defaultdict


def load_offsets_in_both(db_path, tag_type):
    """
    Load file offsets for reads that match between A and B for specific tag type

    Returns:
        tuple: (a_offsets_set, b_offsets_set, isotag_to_offsets)
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute(
        "SELECT isotag, a_offsets, b_offsets FROM matches WHERE tag_type = ?",
        (tag_type,)
    )

    a_offsets_in_both = set()
    b_offsets_in_both = set()
    isotag_to_offsets = {}

    for row in cursor.fetchall():
        isotag = row[0]
        a_offsets = pickle.loads(row[1])
        b_offsets = pickle.loads(row[2])

        a_count = len(a_offsets)
        b_count = len(b_offsets)

        # Only consider reads that match (in_both)
        if a_count > 0 and b_count > 0:
            a_offsets_in_both.update(a_offsets)
            b_offsets_in_both.update(b_offsets)
            isotag_to_offsets[isotag] = (a_offsets, b_offsets)

    conn.close()
    return a_offsets_in_both, b_offsets_in_both, isotag_to_offsets


def find_fuzzy_rescued(db_path):
    """
    Find reads rescued by fuzzy matching (matched via XS but not XI)

    Returns:
        dict: {
            'rescued_a_offsets': set of offsets in A,
            'rescued_b_offsets': set of offsets in B,
            'xi_offsets_a': set of XI offsets in A,
            'xi_offsets_b': set of XI offsets in B,
            'xs_offsets_a': set of XS offsets in A,
            'xs_offsets_b': set of XS offsets in B
        }
    """

    print("Loading XI matches (exact structure)...", file=sys.stderr)
    xi_a, xi_b, xi_isotags = load_offsets_in_both(db_path, 'XI')
    print(f"  XI: {len(xi_a):,} reads in A, {len(xi_b):,} reads in B", file=sys.stderr)

    print("Loading XS matches (fuzzy splice ±5bp)...", file=sys.stderr)
    xs_a, xs_b, xs_isotags = load_offsets_in_both(db_path, 'XS')
    print(f"  XS: {len(xs_a):,} reads in A, {len(xs_b):,} reads in B", file=sys.stderr)

    # Rescued reads = matched via XS but NOT via XI
    rescued_a = xs_a - xi_a
    rescued_b = xs_b - xi_b

    print(f"Fuzzy rescued: {len(rescued_a):,} reads in A, {len(rescued_b):,} reads in B", file=sys.stderr)
    print()

    return {
        'rescued_a_offsets': rescued_a,
        'rescued_b_offsets': rescued_b,
        'xi_offsets_a': xi_a,
        'xi_offsets_b': xi_b,
        'xs_offsets_a': xs_a,
        'xs_offsets_b': xs_b,
        'xs_isotags': xs_isotags
    }


def analyze_rescued_reads(rescued_data):
    """Analyze characteristics of rescued reads"""

    rescued_a = rescued_data['rescued_a_offsets']
    rescued_b = rescued_data['rescued_b_offsets']
    xi_a = rescued_data['xi_offsets_a']
    xi_b = rescued_data['xi_offsets_b']
    xs_a = rescued_data['xs_offsets_a']
    xs_b = rescued_data['xs_offsets_b']

    # Calculate rescue rates
    rescue_rate_a = (len(rescued_a) / len(xi_a) * 100) if len(xi_a) > 0 else 0
    rescue_rate_b = (len(rescued_b) / len(xi_b) * 100) if len(xi_b) > 0 else 0

    # Fuzzy fraction
    fuzzy_frac_a = (len(rescued_a) / len(xs_a) * 100) if len(xs_a) > 0 else 0
    fuzzy_frac_b = (len(rescued_b) / len(xs_b) * 100) if len(xs_b) > 0 else 0

    return {
        'rescued_reads_a': len(rescued_a),
        'rescued_reads_b': len(rescued_b),
        'total_rescued': len(rescued_a) + len(rescued_b),
        'rescue_rate_a': rescue_rate_a,
        'rescue_rate_b': rescue_rate_b,
        'fuzzy_fraction_a': fuzzy_frac_a,
        'fuzzy_fraction_b': fuzzy_frac_b,
        'xi_reads_a': len(xi_a),
        'xi_reads_b': len(xi_b),
        'xs_reads_a': len(xs_a),
        'xs_reads_b': len(xs_b)
    }


def write_rescued_report(rescued_data, output_path):
    """Write rescued reads statistics to TSV file"""

    rescued_a = rescued_data['rescued_a_offsets']
    rescued_b = rescued_data['rescued_b_offsets']
    xs_isotags = rescued_data['xs_isotags']

    # Find which XS isotags contain rescued reads
    rescued_isotags = {}

    for isotag, (a_offsets, b_offsets) in xs_isotags.items():
        rescued_in_a = [off for off in a_offsets if off in rescued_a]
        rescued_in_b = [off for off in b_offsets if off in rescued_b]

        if rescued_in_a or rescued_in_b:
            rescued_isotags[isotag] = {
                'rescued_a': len(rescued_in_a),
                'rescued_b': len(rescued_in_b),
                'total_a': len(a_offsets),
                'total_b': len(b_offsets)
            }

    with open(output_path, 'w') as f:
        # Header
        f.write('\t'.join([
            'isotag_xs',
            'rescued_reads_a',
            'rescued_reads_b',
            'total_rescued',
            'total_reads_a',
            'total_reads_b',
            'rescue_fraction_a',
            'rescue_fraction_b'
        ]) + '\n')

        # Sort by total rescued (most important first)
        sorted_isotags = sorted(
            rescued_isotags.items(),
            key=lambda x: x[1]['rescued_a'] + x[1]['rescued_b'],
            reverse=True
        )

        for isotag, stats in sorted_isotags:
            rescued_total = stats['rescued_a'] + stats['rescued_b']
            total_reads = stats['total_a'] + stats['total_b']

            frac_a = (stats['rescued_a'] / stats['total_a'] * 100) if stats['total_a'] > 0 else 0
            frac_b = (stats['rescued_b'] / stats['total_b'] * 100) if stats['total_b'] > 0 else 0

            f.write('\t'.join([
                isotag,
                str(stats['rescued_a']),
                str(stats['rescued_b']),
                str(rescued_total),
                str(stats['total_a']),
                str(stats['total_b']),
                f"{frac_a:.1f}",
                f"{frac_b:.1f}"
            ]) + '\n')


def extract_rescued_reads(rescued_data, bam_a_path, bam_out_path):
    """Extract rescued reads from BAM A to output BAM"""

    rescued_offsets = rescued_data['rescued_a_offsets']

    if not os.path.exists(bam_a_path):
        print(f"Error: BAM file not found: {bam_a_path}", file=sys.stderr)
        return False

    print(f"Extracting {len(rescued_offsets):,} rescued reads to {bam_out_path}...", file=sys.stderr)

    bam_in = pysam.AlignmentFile(bam_a_path, 'rb')
    bam_out = pysam.AlignmentFile(bam_out_path, 'wb', template=bam_in)

    extracted = 0
    for offset in sorted(rescued_offsets):
        try:
            bam_in.seek(offset)
            read = next(bam_in)
            bam_out.write(read)
            extracted += 1

            if extracted % 10000 == 0:
                print(f"  Extracted {extracted:,} reads...", file=sys.stderr)
        except Exception as e:
            print(f"Warning: Could not extract read at offset {offset}: {e}", file=sys.stderr)

    bam_in.close()
    bam_out.close()

    print(f"✓ Extracted {extracted:,} reads to {bam_out_path}", file=sys.stderr)
    return True


def print_summary(stats):
    """Print human-readable summary"""

    print("="*80)
    print("Fuzzy Match Rescue Analysis (v2.0 - Offset-Based)")
    print("="*80)
    print()

    print("Exact Matching (XI):")
    print("-"*80)
    print(f"Reads matched in A:      {stats['xi_reads_a']:>10,} reads")
    print(f"Reads matched in B:      {stats['xi_reads_b']:>10,} reads")
    print(f"Total XI matches:        {stats['xi_reads_a'] + stats['xi_reads_b']:>10,} reads")
    print()

    print("Fuzzy Matching (XS ±5bp):")
    print("-"*80)
    print(f"Reads matched in A:      {stats['xs_reads_a']:>10,} reads")
    print(f"Reads matched in B:      {stats['xs_reads_b']:>10,} reads")
    print(f"Total XS matches:        {stats['xs_reads_a'] + stats['xs_reads_b']:>10,} reads")
    print()

    print("Fuzzy Rescue (XS - XI):")
    print("-"*80)
    print(f"Rescued reads in A:      {stats['rescued_reads_a']:>10,} reads ({stats['rescue_rate_a']:.1f}% of XI)")
    print(f"Rescued reads in B:      {stats['rescued_reads_b']:>10,} reads ({stats['rescue_rate_b']:.1f}% of XI)")
    print(f"Total rescued:           {stats['total_rescued']:>10,} reads")
    print()

    print("Fuzzy-Only Fraction:")
    print("-"*80)
    print(f"A: {stats['fuzzy_fraction_a']:.1f}% of XS matches are fuzzy-only (not in XI)")
    print(f"B: {stats['fuzzy_fraction_b']:.1f}% of XS matches are fuzzy-only (not in XI)")
    print()

    print("="*80)
    print("INTERPRETATION")
    print("="*80)
    print()

    avg_rescue = (stats['rescue_rate_a'] + stats['rescue_rate_b']) / 2

    if avg_rescue < 5:
        quality = "EXCELLENT"
        interpretation = "Very few sequencing errors. Fuzzy matching provides minimal benefit."
    elif avg_rescue < 15:
        quality = "GOOD"
        interpretation = "Typical sequencing quality. Fuzzy matching rescues ~10% additional reads."
    elif avg_rescue < 30:
        quality = "MODERATE"
        interpretation = "Some sequencing errors. Fuzzy matching provides significant benefit."
    else:
        quality = "POOR"
        interpretation = "Many sequencing errors or highly divergent samples. Check data quality!"

    print(f"Data Quality: {quality}")
    print(f"Rescue Rate:  {avg_rescue:.1f}%")
    print()
    print(f"• {interpretation}")
    print()

    print("Rescued reads likely represent:")
    print("• Sequencing errors (1-5bp shifts in splice junctions)")
    print("• PCR or RT artifacts")
    print("• Minor sequence variants")
    print("• Alignment boundary differences")
    print()

    print("RECOMMENDATIONS:")
    print("-"*80)
    if avg_rescue < 5:
        print("✓ Your data is high quality - exact matching (XI) is sufficient")
    elif avg_rescue < 15:
        print("✓ Fuzzy matching (XS) recommended - rescues ~10% additional reads")
    elif avg_rescue < 30:
        print("⚠ Consider using XS instead of XI for better sensitivity")
    else:
        print("⚠ High rescue rate suggests data quality issues")
        print("  - Check sequencing quality metrics")
        print("  - Verify alignment parameters")
        print("  - Consider resequencing if critical")

    print()
    print("="*80)
    print()


def main():
    """Main entry point"""

    parser = argparse.ArgumentParser(
        description='isotag_fuzzy_rescue - Find reads recovered by fuzzy matching',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  isotag_fuzzy_rescue.py --index matches.db --output rescued.tsv

  # With BAM extraction
  isotag_fuzzy_rescue.py --index matches.db --output rescued.tsv \\
      --bam-a experiment.bam --bam-out rescued_reads.bam

Algorithm:
  1. Load offsets for reads matched via XI (exact structure)
  2. Load offsets for reads matched via XS (fuzzy ±5bp splice)
  3. Rescued = (XS matched) - (XI matched)

Interpretation:
  • Low rescue rate (<5%) = High quality data
  • Medium rescue rate (5-15%) = Typical sequencing quality
  • High rescue rate (>30%) = Data quality issues

Use cases:
  • Evaluate benefit of fuzzy matching
  • Assess sequencing quality
  • Decide between XI (exact) vs XS (fuzzy) matching
  • Identify reads to validate manually

Version 2.0 Changes:
  • Fixed: Now uses OFFSET-BASED comparison (correct!)
  • Old v1.0 compared isotag VALUES across tag types (wrong!)
  • See archived/isotag_fuzzy_rescue_broken_v1.py for old version
        """
    )

    parser.add_argument('--index', required=True, metavar='FILE',
                       help='Input index database (.db)')
    parser.add_argument('--output', metavar='FILE',
                       help='Output TSV file (optional)')
    parser.add_argument('--bam-a', metavar='FILE',
                       help='BAM file A (experiment) for extraction')
    parser.add_argument('--bam-out', metavar='FILE',
                       help='Output BAM with rescued reads (requires --bam-a)')
    parser.add_argument('-v', '--version', action='version',
                       version='isotag_fuzzy_rescue 2.1.0 (updated with XC reference)')

    args = parser.parse_args()

    # Validate input
    if not os.path.exists(args.index):
        print(f"Error: Index file not found: {args.index}", file=sys.stderr)
        sys.exit(1)

    if args.bam_out and not args.bam_a:
        print("Error: --bam-out requires --bam-a", file=sys.stderr)
        sys.exit(1)

    print("="*80, file=sys.stderr)
    print("Fuzzy Match Rescue Analysis (v2.0)", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(f"Index: {args.index}", file=sys.stderr)
    if args.output:
        print(f"Output: {args.output}", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(file=sys.stderr)

    # Find rescued reads
    rescued_data = find_fuzzy_rescued(args.index)

    # Analyze
    stats = analyze_rescued_reads(rescued_data)

    # Write TSV if requested
    if args.output:
        try:
            write_rescued_report(rescued_data, args.output)
            print(f"✓ Report written to: {args.output}", file=sys.stderr)
        except Exception as e:
            print(f"Error writing report: {e}", file=sys.stderr)
            sys.exit(1)

    # Extract BAM if requested
    if args.bam_out:
        success = extract_rescued_reads(rescued_data, args.bam_a, args.bam_out)
        if not success:
            sys.exit(1)

    print(file=sys.stderr)

    # Print summary
    print_summary(stats)


if __name__ == '__main__':
    main()
