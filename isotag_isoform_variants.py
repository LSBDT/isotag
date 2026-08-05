#!/usr/bin/env python3
"""
isotag_isoform_variants.py - Discover Alternative Isoforms

Finds genes with multiple isoforms by grouping reads by XT (transcript group)
and identifying different XI (structures). This reveals alternative splicing.

Usage:
    python3 isotag_isoform_variants.py --index matches.db --output variants.tsv

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


def build_offset_to_isotag_map(tag_data):
    """Build mapping from file offset to isotag value"""

    offset_to_isotag = {}

    for isotag, (a_offsets, b_offsets) in tag_data.items():
        for offset in a_offsets:
            offset_to_isotag[offset] = isotag

    return offset_to_isotag


def find_isoform_variants(xt_data, xi_data):
    """
    Find genes with multiple isoforms

    Returns:
        dict: {xt_isotag: {xi_isotag: [offsets]}}
    """

    print("Building offset mappings...", file=sys.stderr)

    # Build offset → isotag mappings
    offset_to_xt = build_offset_to_isotag_map(xt_data)
    offset_to_xi = build_offset_to_isotag_map(xi_data)

    print(f"  XT mapping: {len(offset_to_xt):,} offsets", file=sys.stderr)
    print(f"  XI mapping: {len(offset_to_xi):,} offsets", file=sys.stderr)

    # Group by XT, collect XI variants
    xt_to_xi_variants = defaultdict(lambda: defaultdict(list))

    for offset in offset_to_xt:
        if offset in offset_to_xi:
            xt = offset_to_xt[offset]
            xi = offset_to_xi[offset]
            xt_to_xi_variants[xt][xi].append(offset)

    # Filter to only XT with multiple XI
    variants = {}
    for xt, xi_dict in xt_to_xi_variants.items():
        if len(xi_dict) > 1:  # Multiple structures
            variants[xt] = dict(xi_dict)

    print(f"Found {len(variants):,} genes with multiple isoforms", file=sys.stderr)

    return variants


def classify_variant_type(num_structures):
    """Classify variant based on number of structures"""

    if num_structures >= 5:
        return "major_alternative_splicing"
    elif num_structures >= 3:
        return "minor_alternative_splicing"
    elif num_structures == 2:
        return "simple_alternative"
    else:
        return "single_isoform"


def write_variants_report(variants, output_path):
    """Write variants to TSV file"""

    with open(output_path, 'w') as f:
        # Header
        f.write('\t'.join([
            'transcript_group_xt',
            'num_structures',
            'structures_xi',
            'reads_per_structure',
            'total_reads',
            'dominant_structure',
            'dominant_fraction',
            'type'
        ]) + '\n')

        # Sort by number of structures (most complex first)
        sorted_variants = sorted(
            variants.items(),
            key=lambda x: len(x[1]),
            reverse=True
        )

        for xt, xi_dict in sorted_variants:
            num_structures = len(xi_dict)

            # Get read counts
            xi_counts = {xi: len(offsets) for xi, offsets in xi_dict.items()}
            total_reads = sum(xi_counts.values())

            # Find dominant structure
            dominant_xi = max(xi_counts, key=xi_counts.get)
            dominant_count = xi_counts[dominant_xi]
            dominant_fraction = dominant_count / total_reads if total_reads > 0 else 0

            # Format structures list (truncate if too many)
            structures_list = list(xi_dict.keys())
            if len(structures_list) <= 5:
                structures_str = ', '.join(structures_list)
            else:
                structures_str = ', '.join(structures_list[:5]) + f', ... ({len(structures_list)-5} more)'

            # Format read counts
            reads_per_str = ', '.join([f"{xi_counts[xi]}" for xi in structures_list[:5]])
            if len(structures_list) > 5:
                reads_per_str += ', ...'

            # Classify
            variant_type = classify_variant_type(num_structures)

            f.write('\t'.join([
                xt,
                str(num_structures),
                structures_str,
                reads_per_str,
                str(total_reads),
                dominant_xi,
                f"{dominant_fraction:.3f}",
                variant_type
            ]) + '\n')


def print_summary(variants):
    """Print summary statistics"""

    if not variants:
        print("No genes with multiple isoforms found")
        return

    print()
    print("="*80)
    print("Alternative Isoform Discovery Summary")
    print("="*80)
    print()

    # Count structures
    structure_counts = [len(xi_dict) for xi_dict in variants.values()]
    total_variants = len(variants)
    max_structures = max(structure_counts)
    avg_structures = sum(structure_counts) / len(structure_counts)

    print(f"Genes with multiple isoforms: {total_variants:,}")
    print(f"Average isoforms per gene:    {avg_structures:.2f}")
    print(f"Maximum isoforms (one gene):  {max_structures}")
    print()

    # Distribution
    print("Isoform distribution:")
    print(f"{'Num Isoforms':<15} {'Genes':>10} {'Percentage':>12}")
    print("-"*80)

    dist = defaultdict(int)
    for count in structure_counts:
        dist[count] += 1

    for num_iso in sorted(dist.keys(), reverse=True):
        gene_count = dist[num_iso]
        pct = gene_count / total_variants * 100
        print(f"{num_iso:<15} {gene_count:>10,} {pct:>11.1f}%")

    print()

    # Top 10 most complex
    sorted_by_complexity = sorted(
        variants.items(),
        key=lambda x: len(x[1]),
        reverse=True
    )[:10]

    print("Top 10 most complex genes (by isoform count):")
    print(f"{'Rank':<6} {'Transcript_Group':<35} {'Isoforms':>10} {'Total_Reads':>12}")
    print("-"*80)

    for rank, (xt, xi_dict) in enumerate(sorted_by_complexity, 1):
        num_iso = len(xi_dict)
        total_reads = sum(len(offsets) for offsets in xi_dict.values())
        xt_short = xt[:32] + "..." if len(xt) > 35 else xt
        print(f"{rank:<6} {xt_short:<35} {num_iso:>10} {total_reads:>12,}")

    print()
    print("="*80)
    print("BIOLOGICAL INSIGHTS")
    print("="*80)

    if max_structures >= 5:
        print(f"• Found gene with {max_structures} isoforms - highly complex regulation!")
        print("  → Likely tissue-specific or condition-specific expression")
        print("  → Consider as candidate for detailed validation")

    major_alt = sum(1 for c in structure_counts if c >= 5)
    if major_alt > 0:
        print(f"• {major_alt} genes with ≥5 isoforms (major alternative splicing)")
        print("  → These represent the most complex transcriptional programs")

    minor_alt = sum(1 for c in structure_counts if 3 <= c < 5)
    if minor_alt > 0:
        print(f"• {minor_alt} genes with 3-4 isoforms (minor alternative splicing)")

    simple_alt = sum(1 for c in structure_counts if c == 2)
    if simple_alt > 0:
        print(f"• {simple_alt} genes with 2 isoforms (simple alternative)")
        print("  → May represent cassette exons or alternative TSS/polyA sites")

    print()
    print("RECOMMENDATIONS:")
    print("-"*80)
    print("• Focus on genes with ≥3 isoforms for detailed analysis")
    print("• Validate top candidates with qPCR or targeted sequencing")
    print("• Check if isoforms correspond to known tissue-specific variants")
    print("• Use XI isotag values to retrieve exact exon structures")
    print("="*80)
    print()


def main():
    """Main entry point"""

    parser = argparse.ArgumentParser(
        description='isotag_isoform_variants - Discover alternative isoforms',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  isotag_isoform_variants.py --index matches.db --output variants.tsv

  # Just print summary (no file output)
  isotag_isoform_variants.py --index matches.db

Output format (TSV):
  transcript_group_xt  - XT isotag (transcript group)
  num_structures       - Number of different XI structures
  structures_xi        - List of XI isotag values
  reads_per_structure  - Read counts for each structure
  total_reads          - Total reads for this gene
  dominant_structure   - Most abundant XI isotag
  dominant_fraction    - Fraction of reads with dominant structure
  type                 - Classification (major/minor/simple alternative)

Interpretation:
  • XT groups reads from same gene/transcript
  • Multiple XI within one XT = alternative splicing
  • More isoforms = more complex regulation
  • Dominant fraction < 0.5 = balanced isoform usage

Use cases:
  • Discover alternative splicing events
  • Find most complex genes
  • Prioritize candidates for validation
  • Compare isoform complexity across samples
        """
    )

    parser.add_argument('--index', required=True, metavar='FILE',
                       help='Input index database (.db)')
    parser.add_argument('--output', metavar='FILE',
                       help='Output TSV file (optional)')
    parser.add_argument('--min-isoforms', type=int, default=2, metavar='N',
                       help='Minimum isoforms to report (default: 2)')
    parser.add_argument('-v', '--version', action='version',
                       version='isotag_isoform_variants 1.0.0')

    args = parser.parse_args()

    # Validate input
    if not os.path.exists(args.index):
        print(f"Error: Index file not found: {args.index}", file=sys.stderr)
        sys.exit(1)

    # Check if index has XT and XI tags
    conn = sqlite3.connect(args.index)
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT tag_type FROM matches")
    tag_types = set(row[0] for row in cursor.fetchall())
    conn.close()

    if 'XT' not in tag_types:
        print("Error: Index must contain XT (transcript group) tag", file=sys.stderr)
        print(f"Available tags: {', '.join(sorted(tag_types))}", file=sys.stderr)
        sys.exit(1)

    if 'XI' not in tag_types:
        print("Error: Index must contain XI (structure) tag", file=sys.stderr)
        print(f"Available tags: {', '.join(sorted(tag_types))}", file=sys.stderr)
        sys.exit(1)

    print("="*80, file=sys.stderr)
    print("Alternative Isoform Discovery", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(f"Index: {args.index}", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(file=sys.stderr)

    # Load data
    print("Loading XT (transcript groups)...", file=sys.stderr)
    try:
        xt_data = load_tag_data(args.index, 'XT')
        print(f"  Loaded {len(xt_data):,} XT isotags", file=sys.stderr)
    except Exception as e:
        print(f"Error loading XT data: {e}", file=sys.stderr)
        sys.exit(1)

    print("Loading XI (structures)...", file=sys.stderr)
    try:
        xi_data = load_tag_data(args.index, 'XI')
        print(f"  Loaded {len(xi_data):,} XI isotags", file=sys.stderr)
    except Exception as e:
        print(f"Error loading XI data: {e}", file=sys.stderr)
        sys.exit(1)

    print(file=sys.stderr)

    # Find variants
    print("Finding genes with multiple isoforms...", file=sys.stderr)
    variants = find_isoform_variants(xt_data, xi_data)

    # Filter by minimum isoforms
    if args.min_isoforms > 2:
        variants = {xt: xi_dict for xt, xi_dict in variants.items()
                   if len(xi_dict) >= args.min_isoforms}
        print(f"Filtered to {len(variants):,} genes with ≥{args.min_isoforms} isoforms", file=sys.stderr)

    print(file=sys.stderr)

    # Print summary
    print_summary(variants)

    # Write TSV if requested
    if args.output:
        try:
            write_variants_report(variants, args.output)
            print(f"✓ Variants written to: {args.output}")
        except Exception as e:
            print(f"Error: Cannot write output file: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    main()
