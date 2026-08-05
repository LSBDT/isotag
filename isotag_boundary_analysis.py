#!/usr/bin/env python3
"""
isotag_boundary_analysis.py - Analyze TSS/PolyA Site Variation

Identifies boundary variation by comparing XI (structure) and XB (boundaries).
Finds alternative TSS and polyadenylation sites.

Usage:
    python3 isotag_boundary_analysis.py --index matches.db --output boundaries.tsv

Author: IsoTag Team
Version: 1.0.0
"""

import sqlite3
import pickle
import os
import sys
import argparse
from collections import defaultdict


def load_tag_data(db_path, tag_type):
    """Load isotag data for specific tag type"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute(
        "SELECT isotag, a_offsets, b_offsets FROM matches WHERE tag_type = ?",
        (tag_type,)
    )

    data = {}
    for row in cursor.fetchall():
        isotag = row[0]
        a_offsets = pickle.loads(row[1])
        b_offsets = pickle.loads(row[2])
        data[isotag] = (a_offsets, b_offsets)

    conn.close()
    return data


def parse_xb_boundary(xb_tag):
    """
    Parse XB tag to extract boundary coordinates
    Format: chr_hash_8 + strand + 5'_hex + 3'_hex
    Example: aKF498dAp.3e8.1004
    """
    try:
        parts = xb_tag.split('.')
        if len(parts) < 3:
            return None

        strand = parts[0][-1]  # Last char is strand
        start_hex = parts[1]
        end_hex = parts[2]

        start = int(start_hex, 16)
        end = int(end_hex, 16)

        return (start, end, strand)
    except:
        return None


def build_offset_to_isotag_map(tag_data):
    """Build mapping from file offset to isotag value"""
    offset_to_isotag = {}

    for isotag, (a_offsets, b_offsets) in tag_data.items():
        for offset in a_offsets:
            offset_to_isotag[offset] = isotag

    return offset_to_isotag


def find_boundary_variants(xi_data, xb_data):
    """
    Find structures with multiple boundary variants

    Returns:
        dict: {xi_isotag: {xb_isotag: [offsets]}}
    """
    print("Building offset mappings...", file=sys.stderr)

    offset_to_xi = build_offset_to_isotag_map(xi_data)
    offset_to_xb = build_offset_to_isotag_map(xb_data)

    print(f"  XI mapping: {len(offset_to_xi):,} offsets", file=sys.stderr)
    print(f"  XB mapping: {len(offset_to_xb):,} offsets", file=sys.stderr)

    # Group by XI, collect XB variants
    xi_to_xb_variants = defaultdict(lambda: defaultdict(list))

    for offset in offset_to_xi:
        if offset in offset_to_xb:
            xi = offset_to_xi[offset]
            xb = offset_to_xb[offset]
            xi_to_xb_variants[xi][xb].append(offset)

    # Filter to only XI with multiple XB
    variants = {}
    for xi, xb_dict in xi_to_xb_variants.items():
        if len(xb_dict) > 1:  # Multiple boundaries
            variants[xi] = dict(xb_dict)

    print(f"Found {len(variants):,} structures with boundary variation", file=sys.stderr)

    return variants


def analyze_boundary_shifts(variants):
    """Analyze characteristics of boundary shifts"""

    tss_variants = 0
    polya_variants = 0
    both_variants = 0

    tss_shifts = []
    polya_shifts = []

    for xi, xb_dict in variants.items():
        xb_coords = []
        for xb in xb_dict.keys():
            parsed = parse_xb_boundary(xb)
            if parsed:
                xb_coords.append(parsed)

        if len(xb_coords) < 2:
            continue

        # Analyze shifts
        starts = [c[0] for c in xb_coords]
        ends = [c[1] for c in xb_coords]

        start_variation = max(starts) - min(starts)
        end_variation = max(ends) - min(ends)

        if start_variation > 0:
            tss_shifts.append(start_variation)
        if end_variation > 0:
            polya_shifts.append(end_variation)

        # Classify
        if start_variation > 0 and end_variation > 0:
            both_variants += 1
        elif start_variation > 0:
            tss_variants += 1
        elif end_variation > 0:
            polya_variants += 1

    return {
        'tss_variants': tss_variants,
        'polya_variants': polya_variants,
        'both_variants': both_variants,
        'tss_shifts': tss_shifts,
        'polya_shifts': polya_shifts
    }


def write_boundary_report(variants, output_path):
    """Write boundary variants to TSV"""

    with open(output_path, 'w') as f:
        # Header
        f.write('\t'.join([
            'structure_xi',
            'num_boundaries',
            'boundary_xb_list',
            'reads_per_boundary',
            'total_reads',
            'dominant_boundary',
            'dominant_fraction',
            'has_tss_variation',
            'has_polya_variation'
        ]) + '\n')

        # Sort by number of boundaries
        sorted_variants = sorted(
            variants.items(),
            key=lambda x: len(x[1]),
            reverse=True
        )

        for xi, xb_dict in sorted_variants:
            num_boundaries = len(xb_dict)

            # Get read counts
            xb_counts = {xb: len(offsets) for xb, offsets in xb_dict.items()}
            total_reads = sum(xb_counts.values())

            # Find dominant
            dominant_xb = max(xb_counts, key=xb_counts.get)
            dominant_count = xb_counts[dominant_xb]
            dominant_fraction = dominant_count / total_reads if total_reads > 0 else 0

            # Check for TSS/polyA variation
            xb_coords = []
            for xb in xb_dict.keys():
                parsed = parse_xb_boundary(xb)
                if parsed:
                    xb_coords.append(parsed)

            has_tss = False
            has_polya = False
            if len(xb_coords) >= 2:
                starts = [c[0] for c in xb_coords]
                ends = [c[1] for c in xb_coords]
                if max(starts) - min(starts) > 0:
                    has_tss = True
                if max(ends) - min(ends) > 0:
                    has_polya = True

            # Format lists
            xb_list = list(xb_dict.keys())
            if len(xb_list) <= 5:
                xb_str = ', '.join(xb_list)
            else:
                xb_str = ', '.join(xb_list[:5]) + f', ... ({len(xb_list)-5} more)'

            reads_per = ', '.join([str(xb_counts[xb]) for xb in xb_list[:5]])
            if len(xb_list) > 5:
                reads_per += ', ...'

            f.write('\t'.join([
                xi,
                str(num_boundaries),
                xb_str,
                reads_per,
                str(total_reads),
                dominant_xb,
                f"{dominant_fraction:.3f}",
                'YES' if has_tss else 'NO',
                'YES' if has_polya else 'NO'
            ]) + '\n')


def print_summary(variants, stats):
    """Print summary statistics"""

    if not variants:
        print("No boundary variants found")
        return

    print()
    print("="*80)
    print("Boundary Variation Analysis")
    print("="*80)
    print()

    print(f"Structures with boundary variation: {len(variants):,}")
    print()

    # Variation types
    print("Variation types:")
    print(f"  TSS variation only:      {stats['tss_variants']:>10,} ({stats['tss_variants']/len(variants)*100:>5.1f}%)")
    print(f"  PolyA variation only:    {stats['polya_variants']:>10,} ({stats['polya_variants']/len(variants)*100:>5.1f}%)")
    print(f"  Both TSS and PolyA:      {stats['both_variants']:>10,} ({stats['both_variants']/len(variants)*100:>5.1f}%)")
    print()

    # Shift statistics
    if stats['tss_shifts']:
        import statistics
        print("TSS shifts:")
        print(f"  Median:  {statistics.median(stats['tss_shifts']):>6.0f} bp")
        print(f"  Mean:    {statistics.mean(stats['tss_shifts']):>6.1f} bp")
        print(f"  Range:   {min(stats['tss_shifts']):>6.0f} - {max(stats['tss_shifts']):>6.0f} bp")
        print()

    if stats['polya_shifts']:
        import statistics
        print("PolyA shifts:")
        print(f"  Median:  {statistics.median(stats['polya_shifts']):>6.0f} bp")
        print(f"  Mean:    {statistics.mean(stats['polya_shifts']):>6.1f} bp")
        print(f"  Range:   {min(stats['polya_shifts']):>6.0f} - {max(stats['polya_shifts']):>6.0f} bp")
        print()

    # Top variants
    sorted_by_boundaries = sorted(
        variants.items(),
        key=lambda x: len(x[1]),
        reverse=True
    )[:10]

    print("Top 10 structures by boundary variation:")
    print(f"{'Rank':<6} {'Structure_XI':<35} {'Boundaries':>12} {'Total_Reads':>12}")
    print("-"*80)

    for rank, (xi, xb_dict) in enumerate(sorted_by_boundaries, 1):
        num_boundaries = len(xb_dict)
        total_reads = sum(len(offsets) for offsets in xb_dict.values())
        xi_short = xi[:32] + "..." if len(xi) > 35 else xi
        print(f"{rank:<6} {xi_short:<35} {num_boundaries:>12} {total_reads:>12,}")

    print()
    print("="*80)
    print("BIOLOGICAL INTERPRETATION")
    print("="*80)

    if stats['tss_variants'] > stats['polya_variants']:
        print("• TSS variation dominant - suggests alternative promoter usage")
        print("  → Different transcription start sites")
        print("  → May indicate tissue-specific or condition-specific regulation")
    elif stats['polya_variants'] > stats['tss_variants']:
        print("• PolyA variation dominant - suggests alternative polyadenylation")
        print("  → Different 3' UTR lengths")
        print("  → May affect mRNA stability and localization")

    if stats['both_variants'] > 0:
        pct = stats['both_variants'] / len(variants) * 100
        if pct > 10:
            print(f"• {pct:.1f}% have both TSS and polyA variation")
            print("  → May represent different transcript isoforms rather than just boundary shifts")

    print()
    print("RECOMMENDATIONS:")
    print("-"*80)
    print("• Validate top candidates with 5' and 3' RACE")
    print("• Check if boundary shifts correspond to known regulatory elements")
    print("• Compare across tissues/conditions for functional insights")
    print("• Use XB tags to precisely map alternative start/end sites")
    print("="*80)
    print()


def main():
    """Main entry point"""

    parser = argparse.ArgumentParser(
        description='isotag_boundary_analysis - Analyze TSS/PolyA variation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  isotag_boundary_analysis.py --index matches.db --output boundaries.tsv

  # Just print summary
  isotag_boundary_analysis.py --index matches.db

Output format (TSV):
  structure_xi         - XI isotag (exon structure)
  num_boundaries       - Number of different boundaries
  boundary_xb_list     - List of XB tags
  reads_per_boundary   - Read counts
  total_reads          - Total reads
  dominant_boundary    - Most common boundary
  dominant_fraction    - Fraction with dominant boundary
  has_tss_variation    - YES/NO
  has_polya_variation  - YES/NO

Interpretation:
  • Same XI (structure) with different XB (boundaries) = boundary variation
  • TSS variation = alternative transcription start sites
  • PolyA variation = alternative polyadenylation signals
  • Both = complex regulation or different transcripts

Use cases:
  • Find alternative promoters
  • Discover alternative polyadenylation
  • Map precise TSS and polyA sites
  • Compare boundary usage across conditions
        """
    )

    parser.add_argument('--index', required=True, metavar='FILE',
                       help='Input index database (.db)')
    parser.add_argument('--output', metavar='FILE',
                       help='Output TSV file (optional)')
    parser.add_argument('-v', '--version', action='version',
                       version='isotag_boundary_analysis 1.0.0')

    args = parser.parse_args()

    # Validate input
    if not os.path.exists(args.index):
        print(f"Error: Index file not found: {args.index}", file=sys.stderr)
        sys.exit(1)

    # Check tags
    conn = sqlite3.connect(args.index)
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT tag_type FROM matches")
    tag_types = set(row[0] for row in cursor.fetchall())
    conn.close()

    if 'XI' not in tag_types or 'XB' not in tag_types:
        print("Error: Index must contain both XI and XB tags", file=sys.stderr)
        print(f"Available tags: {', '.join(sorted(tag_types))}", file=sys.stderr)
        sys.exit(1)

    print("="*80, file=sys.stderr)
    print("Boundary Variation Analysis", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(f"Index: {args.index}", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(file=sys.stderr)

    # Load data
    print("Loading XI (structures)...", file=sys.stderr)
    xi_data = load_tag_data(args.index, 'XI')
    print(f"  Loaded {len(xi_data):,} XI isotags", file=sys.stderr)

    print("Loading XB (boundaries)...", file=sys.stderr)
    xb_data = load_tag_data(args.index, 'XB')
    print(f"  Loaded {len(xb_data):,} XB isotags", file=sys.stderr)
    print(file=sys.stderr)

    # Find variants
    print("Finding boundary variants...", file=sys.stderr)
    variants = find_boundary_variants(xi_data, xb_data)
    print(file=sys.stderr)

    if not variants:
        print("No boundary variants found", file=sys.stderr)
        sys.exit(0)

    # Analyze
    print("Analyzing boundary shifts...", file=sys.stderr)
    stats = analyze_boundary_shifts(variants)
    print(file=sys.stderr)

    # Print summary
    print_summary(variants, stats)

    # Write TSV
    if args.output:
        try:
            write_boundary_report(variants, args.output)
            print(f"✓ Boundary variants written to: {args.output}")
        except Exception as e:
            print(f"Error: Cannot write output: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    main()
