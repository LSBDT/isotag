#!/usr/bin/env python3
"""
isotag_stats_compare.py - Unified Statistical Comparison Tool

Replaces and combines three separate tools:
  - isotag_tag_compare.py (matching rates across tags)
  - isotag_venn_ab.py (dataset A vs B Venn diagrams)
  - isotag_concordance.py (tag agreement matrix)

Three modes of operation:
  1. tags      - Compare matching rates across tag types (TSV output)
  2. venn      - Generate Venn diagram data for A vs B (JSON output)
  3. concordance - Calculate tag agreement matrix (TSV output)

Usage:
    # Tag comparison mode
    python3 isotag_stats_compare.py --index matches.db --mode tags -o compare.tsv

    # Venn diagram mode
    python3 isotag_stats_compare.py --index matches.db --mode venn --tag XI -o venn.json

    # Concordance matrix mode
    python3 isotag_stats_compare.py --index matches.db --mode concordance -o concordance.tsv

Author: IsoTag Team
Version: 2.0.0 (Unified tool)
"""

import sqlite3
import pickle
import json
import os
import sys
import argparse
from collections import defaultdict


# ==============================================================================
# MODE 1: TAG COMPARISON (replaces isotag_tag_compare.py)
# ==============================================================================

def load_tag_comparison_stats(db_path):
    """Load statistics for all tag types from index"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("SELECT DISTINCT tag_type FROM matches")
    tag_types = sorted([row[0] for row in cursor.fetchall()])

    stats = {}

    for tag in tag_types:
        cursor.execute(
            "SELECT isotag, a_offsets, b_offsets FROM matches WHERE tag_type = ?",
            (tag,)
        )

        isotags_in_both = 0
        isotags_only_a = 0
        isotags_only_b = 0
        reads_in_both_a = 0
        reads_in_both_b = 0
        reads_only_a = 0
        reads_only_b = 0

        for row in cursor.fetchall():
            isotag = row[0]
            a_offsets = pickle.loads(row[1])
            b_offsets = pickle.loads(row[2])

            a_count = len(a_offsets)
            b_count = len(b_offsets)

            if a_count > 0 and b_count > 0:
                isotags_in_both += 1
                reads_in_both_a += a_count
                reads_in_both_b += b_count
            elif a_count > 0:
                isotags_only_a += 1
                reads_only_a += a_count
            else:
                isotags_only_b += 1
                reads_only_b += b_count

        total_isotags = isotags_in_both + isotags_only_a + isotags_only_b
        match_rate = (isotags_in_both / total_isotags) if total_isotags > 0 else 0

        sensitivity = (isotags_in_both / (isotags_in_both + isotags_only_a)) if (isotags_in_both + isotags_only_a) > 0 else 0
        specificity = (isotags_in_both / (isotags_in_both + isotags_only_b)) if (isotags_in_both + isotags_only_b) > 0 else 0

        stats[tag] = {
            'isotags_in_both': isotags_in_both,
            'isotags_only_a': isotags_only_a,
            'isotags_only_b': isotags_only_b,
            'total_isotags': total_isotags,
            'match_rate': match_rate,
            'reads_in_both_a': reads_in_both_a,
            'reads_in_both_b': reads_in_both_b,
            'reads_only_a': reads_only_a,
            'reads_only_b': reads_only_b,
            'sensitivity': sensitivity,
            'specificity': specificity
        }

    conn.close()
    return stats


def write_tag_comparison_tsv(stats, output_path):
    """Write comparison report to TSV file"""
    with open(output_path, 'w') as f:
        f.write('\t'.join([
            'tag',
            'isotags_in_both',
            'isotags_only_a',
            'isotags_only_b',
            'total_isotags',
            'match_rate',
            'reads_in_both_a',
            'reads_in_both_b',
            'reads_only_a',
            'reads_only_b',
            'sensitivity',
            'specificity'
        ]) + '\n')

        for tag in sorted(stats.keys()):
            s = stats[tag]
            f.write('\t'.join([
                tag,
                str(s['isotags_in_both']),
                str(s['isotags_only_a']),
                str(s['isotags_only_b']),
                str(s['total_isotags']),
                f"{s['match_rate']:.6f}",
                str(s['reads_in_both_a']),
                str(s['reads_in_both_b']),
                str(s['reads_only_a']),
                str(s['reads_only_b']),
                f"{s['sensitivity']:.6f}",
                f"{s['specificity']:.6f}"
            ]) + '\n')


def print_tag_comparison_summary(stats):
    """Print human-readable summary"""
    print("="*80)
    print("Tag Comparison Summary")
    print("="*80)
    print()

    best_tag = max(stats.keys(), key=lambda t: stats[t]['match_rate'])

    print(f"{'Tag':<6} {'Match_Rate':>12} {'Sensitivity':>12} {'Specificity':>12} {'Isotags':>10}")
    print("-"*80)

    for tag in sorted(stats.keys()):
        s = stats[tag]
        marker = " ← BEST" if tag == best_tag else ""
        print(f"{tag:<6} {s['match_rate']:>11.1%} {s['sensitivity']:>11.1%} {s['specificity']:>11.1%} "
              f"{s['total_isotags']:>10,}{marker}")

    print()
    print(f"✓ Best match rate: {best_tag} ({stats[best_tag]['match_rate']:.1%})")

    if 'XI' in stats and 'XS' in stats:
        xi_matches = stats['XI']['isotags_in_both']
        xs_matches = stats['XS']['isotags_in_both']
        rescued = xs_matches - xi_matches
        if rescued > 0:
            print(f"✓ XS fuzzy matching rescued {rescued:,} isotags vs XI ({rescued/xi_matches:.1%} increase)")

    print("="*80)
    print()


# ==============================================================================
# MODE 2: VENN DIAGRAMS (replaces isotag_venn_ab.py)
# ==============================================================================

def load_venn_ab_sets(db_path, tag_type):
    """Load sets of isotags present in A and B for specific tag type"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute(
        "SELECT isotag, a_offsets, b_offsets FROM matches WHERE tag_type = ?",
        (tag_type,)
    )

    set_a = set()
    set_b = set()
    in_both = set()
    only_a = set()
    only_b = set()

    reads_in_a = 0
    reads_in_b = 0
    reads_in_both_a = 0
    reads_in_both_b = 0

    for row in cursor.fetchall():
        isotag = row[0]
        a_offsets = pickle.loads(row[1])
        b_offsets = pickle.loads(row[2])

        a_count = len(a_offsets)
        b_count = len(b_offsets)

        if a_count > 0:
            set_a.add(isotag)
            reads_in_a += a_count

        if b_count > 0:
            set_b.add(isotag)
            reads_in_b += b_count

        if a_count > 0 and b_count > 0:
            in_both.add(isotag)
            reads_in_both_a += a_count
            reads_in_both_b += b_count
        elif a_count > 0:
            only_a.add(isotag)
        elif b_count > 0:
            only_b.add(isotag)

    conn.close()

    return {
        'set_a': set_a,
        'set_b': set_b,
        'in_both': in_both,
        'only_a': only_a,
        'only_b': only_b,
        'reads_in_a': reads_in_a,
        'reads_in_b': reads_in_b,
        'reads_in_both_a': reads_in_both_a,
        'reads_in_both_b': reads_in_both_b
    }


def calculate_venn_stats(stats):
    """Calculate Venn diagram statistics"""
    n_a = len(stats['set_a'])
    n_b = len(stats['set_b'])
    n_both = len(stats['in_both'])
    n_only_a = len(stats['only_a'])
    n_only_b = len(stats['only_b'])
    n_union = len(stats['set_a'] | stats['set_b'])

    pct_a_in_both = (n_both / n_a * 100) if n_a > 0 else 0
    pct_b_in_both = (n_both / n_b * 100) if n_b > 0 else 0
    pct_a_only = (n_only_a / n_a * 100) if n_a > 0 else 0
    pct_b_only = (n_only_b / n_b * 100) if n_b > 0 else 0

    jaccard = (n_both / n_union) if n_union > 0 else 0

    return {
        'n_a': n_a,
        'n_b': n_b,
        'n_both': n_both,
        'n_only_a': n_only_a,
        'n_only_b': n_only_b,
        'n_union': n_union,
        'pct_a_in_both': pct_a_in_both,
        'pct_b_in_both': pct_b_in_both,
        'pct_a_only': pct_a_only,
        'pct_b_only': pct_b_only,
        'jaccard': jaccard,
        'reads_in_a': stats['reads_in_a'],
        'reads_in_b': stats['reads_in_b'],
        'reads_in_both_a': stats['reads_in_both_a'],
        'reads_in_both_b': stats['reads_in_both_b']
    }


def write_venn_json(tag_stats, output_path, metadata):
    """Write Venn diagram data to JSON"""
    data = {
        'format': 'venn_diagram_ab',
        'version': '2.0',
        'description': 'Venn diagram comparing dataset A vs B for each isotag type',
        'metadata': metadata,
        'tags': {}
    }

    for tag, stats in tag_stats.items():
        data['tags'][tag] = stats

    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)


def print_venn_summary(tag_stats, metadata):
    """Print human-readable Venn summary"""
    print("="*80)
    print("Venn Diagram: Dataset A vs B")
    print("="*80)
    print()
    print(f"Dataset A: {metadata['bam_a']}")
    print(f"Dataset B: {metadata['bam_b']}")
    print()

    for tag in sorted(tag_stats.keys()):
        stats = tag_stats[tag]

        print("="*80)
        print(f"Tag Type: {tag}")
        print("="*80)
        print()

        print("ISOTAG COUNTS:")
        print(f"  Dataset A total:     {stats['n_a']:>10,} isotags")
        print(f"  Dataset B total:     {stats['n_b']:>10,} isotags")
        print(f"  In both (overlap):   {stats['n_both']:>10,} isotags")
        print(f"  Only in A (novel):   {stats['n_only_a']:>10,} isotags")
        print(f"  Only in B:           {stats['n_only_b']:>10,} isotags")
        print(f"  Union (A ∪ B):       {stats['n_union']:>10,} isotags")
        print()

        print("PERCENTAGES:")
        print(f"  A matched in B:      {stats['pct_a_in_both']:>10.1f}% ({stats['n_both']:,} / {stats['n_a']:,})")
        print(f"  B matched in A:      {stats['pct_b_in_both']:>10.1f}% ({stats['n_both']:,} / {stats['n_b']:,})")
        print(f"  A-specific (novel):  {stats['pct_a_only']:>10.1f}% ({stats['n_only_a']:,} / {stats['n_a']:,})")
        print()

        print("SIMILARITY:")
        print(f"  Jaccard index:       {stats['jaccard']:>10.3f} (overlap/union)")
        print()

    print("="*80)
    print()


def generate_r_code(output_path, tag_stats, metadata):
    """Generate R code for plotting Venn diagrams"""
    r_code = f"""# R code to plot Venn diagrams (Dataset A vs B)
# Install package if needed: install.packages("VennDiagram")

library(VennDiagram)
library(jsonlite)
library(grid)

# Load data
data <- fromJSON("{os.path.basename(output_path)}")

# Metadata
bam_a <- "{os.path.basename(metadata['bam_a'])}"
bam_b <- "{os.path.basename(metadata['bam_b'])}"

"""

    for tag in sorted(tag_stats.keys()):
        stats = tag_stats[tag]

        r_code += f"""
# Venn diagram for {tag}
venn_{tag.lower()} <- draw.pairwise.venn(
  area1 = {stats['n_a']},
  area2 = {stats['n_b']},
  cross.area = {stats['n_both']},
  category = c(paste0("A: ", bam_a), paste0("B: ", bam_b)),
  fill = c("lightblue", "pink"),
  alpha = 0.5,
  cat.pos = c(0, 0),
  cat.dist = 0.05,
  main = "{tag} Tag - Dataset A vs B"
)

# Save to file
png("venn_{tag.lower()}.png", width=800, height=800)
grid.draw(venn_{tag.lower()})
dev.off()

cat("Saved venn_{tag.lower()}.png\\n")
"""

    r_code += """
# Summary statistics
cat("\\n=== Venn Diagram Summary ===\\n")
for (tag in names(data$tags)) {
  stats <- data$tags[[tag]]
  cat(sprintf("\\n%s:\\n", tag))
  cat(sprintf("  A: %d isotags\\n", stats$n_a))
  cat(sprintf("  B: %d isotags\\n", stats$n_b))
  cat(sprintf("  Overlap: %d isotags (%.1f%% of A, %.1f%% of B)\\n",
              stats$n_both, stats$pct_a_in_both, stats$pct_b_in_both))
  cat(sprintf("  Jaccard: %.3f\\n", stats$jaccard))
}
"""

    return r_code


# ==============================================================================
# MODE 3: CONCORDANCE (replaces isotag_concordance.py)
# ==============================================================================

def load_tag_data_concordance(db_path, tag_type):
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


def calculate_concordance(tag1_map, tag2_map):
    """Calculate concordance between two tags"""
    all_offsets = set(tag1_map.keys()) | set(tag2_map.keys())

    both_present = 0
    for offset in all_offsets:
        if offset in tag1_map and offset in tag2_map:
            both_present += 1

    concordance = both_present / len(all_offsets) if all_offsets else 0

    return concordance, both_present, len(all_offsets)


def calculate_all_concordances(tag_maps):
    """Calculate pairwise concordance for all tag combinations"""
    tags = sorted(tag_maps.keys())

    concordance_matrix = {}
    for tag1 in tags:
        concordance_matrix[tag1] = {}
        for tag2 in tags:
            if tag1 == tag2:
                concordance_matrix[tag1][tag2] = 1.0
            else:
                concord, both, total = calculate_concordance(tag_maps[tag1], tag_maps[tag2])
                concordance_matrix[tag1][tag2] = concord

    return concordance_matrix


def write_concordance_tsv(concordance_matrix, output_path):
    """Write concordance matrix to TSV"""
    tags = sorted(concordance_matrix.keys())

    with open(output_path, 'w') as f:
        f.write('tag\t' + '\t'.join(tags) + '\n')

        for tag1 in tags:
            row = [tag1]
            for tag2 in tags:
                value = concordance_matrix[tag1][tag2]
                row.append(f"{value:.4f}")
            f.write('\t'.join(row) + '\n')


def print_concordance_summary(concordance_matrix):
    """Print formatted concordance matrix"""
    tags = sorted(concordance_matrix.keys())

    print("="*80)
    print("Tag Concordance Matrix")
    print("="*80)
    print()

    print("Pairwise agreement (fraction of reads present in both tags):")
    print()

    header = "       " + "".join([f"{tag:>8}" for tag in tags])
    print(header)
    print("-"*80)

    for tag1 in tags:
        row = f"{tag1:<7}"
        for tag2 in tags:
            value = concordance_matrix[tag1][tag2]
            if tag1 == tag2:
                row += f"{'100%':>8}"
            else:
                row += f"{value*100:>7.1f}%"
        print(row)

    print()
    print("="*80)
    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    """Main entry point"""

    parser = argparse.ArgumentParser(
        description='isotag_stats_compare - Unified statistical comparison tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This tool replaces three separate tools:
  • isotag_tag_compare.py   → --mode tags
  • isotag_venn_ab.py        → --mode venn
  • isotag_concordance.py    → --mode concordance

Examples:
  # Tag comparison (matching rates across tag types)
  isotag_stats_compare.py --index matches.db --mode tags -o compare.tsv

  # Venn diagram (dataset A vs B overlap)
  isotag_stats_compare.py --index matches.db --mode venn --tag XI -o venn.json
  isotag_stats_compare.py --index matches.db --mode venn --tags XI,XS,XB -o venn_multi.json

  # Venn with R plotting code
  isotag_stats_compare.py --index matches.db --mode venn --tag XI -o venn.json --r-code plot.R

  # Concordance matrix (tag agreement)
  isotag_stats_compare.py --index matches.db --mode concordance -o concordance.tsv

Mode descriptions:
  tags        - Compare matching rates across tag types (TSV output)
                Shows isotags_in_both, only_a, only_b, match_rate, sensitivity, specificity

  venn        - Generate Venn diagram data for A vs B comparison (JSON output)
                Shows overlap, novel isoforms, Jaccard index
                Optionally generates R plotting code

  concordance - Calculate tag agreement matrix (TSV output)
                Shows pairwise concordance between tag types (QC metric)

Output formats:
  • tags mode: TSV with columns per tag type
  • venn mode: JSON with nested tag statistics (+ optional R code)
  • concordance mode: TSV matrix (tags × tags)
        """
    )

    parser.add_argument('--index', required=True, metavar='FILE',
                       help='Input index database (.db)')
    parser.add_argument('--mode', required=True,
                       choices=['tags', 'venn', 'concordance'],
                       help='Comparison mode: tags, venn, or concordance')
    parser.add_argument('-o', '--output', metavar='FILE',
                       help='Output file (TSV or JSON depending on mode)')

    # Venn-specific options
    parser.add_argument('--tag', metavar='TAG',
                       help='Single tag type for venn mode (e.g., XI)')
    parser.add_argument('--tags', metavar='TAGS',
                       help='Comma-separated tag types for venn mode (e.g., XI,XS,XB)')
    parser.add_argument('--r-code', metavar='FILE',
                       help='Generate R plotting code to file (venn mode only)')

    parser.add_argument('-v', '--version', action='version',
                       version='isotag_stats_compare 2.0.0')

    args = parser.parse_args()

    # Validate input
    if not os.path.exists(args.index):
        print(f"Error: Index file not found: {args.index}", file=sys.stderr)
        sys.exit(1)

    # Mode-specific validation
    if args.mode == 'venn':
        if not args.tag and not args.tags:
            print("Error: venn mode requires --tag or --tags", file=sys.stderr)
            sys.exit(1)
        if args.tag and args.tags:
            print("Error: Cannot specify both --tag and --tags", file=sys.stderr)
            sys.exit(1)

    print("="*80, file=sys.stderr)
    print(f"Statistical Comparison - Mode: {args.mode.upper()}", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(f"Index: {args.index}", file=sys.stderr)
    if args.output:
        print(f"Output: {args.output}", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(file=sys.stderr)

    # ==============================================================================
    # MODE: TAGS
    # ==============================================================================
    if args.mode == 'tags':
        print("Loading tag comparison statistics...", file=sys.stderr)
        try:
            stats = load_tag_comparison_stats(args.index)
        except Exception as e:
            print(f"Error: Cannot load index: {e}", file=sys.stderr)
            sys.exit(1)

        if not stats:
            print("Error: No tag types found in index", file=sys.stderr)
            sys.exit(1)

        print(f"Loaded statistics for {len(stats)} tag types", file=sys.stderr)
        print(file=sys.stderr)

        print_tag_comparison_summary(stats)

        if args.output:
            try:
                write_tag_comparison_tsv(stats, args.output)
                print(f"✓ TSV report written to: {args.output}")
            except Exception as e:
                print(f"Error: Cannot write output file: {e}", file=sys.stderr)
                sys.exit(1)

    # ==============================================================================
    # MODE: VENN
    # ==============================================================================
    elif args.mode == 'venn':
        if not args.output:
            print("Error: venn mode requires --output", file=sys.stderr)
            sys.exit(1)

        # Parse tags
        if args.tag:
            tags = [args.tag]
        else:
            tags = [t.strip() for t in args.tags.split(',')]

        valid_tags = {'XI', 'XB', 'XS', 'XT', 'XV', 'XC', 'XR'}
        for tag in tags:
            if tag not in valid_tags:
                print(f"Error: Invalid tag '{tag}'. Valid: {valid_tags}", file=sys.stderr)
                sys.exit(1)

        print(f"Tags: {', '.join(tags)}", file=sys.stderr)
        print(file=sys.stderr)

        # Load metadata
        conn = sqlite3.connect(args.index)
        cursor = conn.cursor()
        cursor.execute("SELECT key, value FROM metadata")
        metadata = dict(cursor.fetchall())
        conn.close()

        # Load A/B sets for each tag
        print("Loading isotag sets...", file=sys.stderr)
        tag_stats = {}

        for tag in tags:
            try:
                print(f"  Processing {tag}...", file=sys.stderr)
                stats = load_venn_ab_sets(args.index, tag)
                venn_stats = calculate_venn_stats(stats)
                tag_stats[tag] = venn_stats
                print(f"    A: {venn_stats['n_a']:,}, B: {venn_stats['n_b']:,}, Both: {venn_stats['n_both']:,}", file=sys.stderr)
            except Exception as e:
                print(f"Error loading {tag}: {e}", file=sys.stderr)
                sys.exit(1)

        print(file=sys.stderr)

        # Write JSON
        try:
            write_venn_json(tag_stats, args.output, metadata)
            print(f"✓ Venn data written to: {args.output}", file=sys.stderr)
        except Exception as e:
            print(f"Error writing output: {e}", file=sys.stderr)
            sys.exit(1)

        # Generate R code if requested
        if args.r_code:
            try:
                r_code = generate_r_code(args.output, tag_stats, metadata)
                with open(args.r_code, 'w') as f:
                    f.write(r_code)
                print(f"✓ R code written to: {args.r_code}", file=sys.stderr)
            except Exception as e:
                print(f"Error writing R code: {e}", file=sys.stderr)

        print(file=sys.stderr)

        # Print summary
        print_venn_summary(tag_stats, metadata)

    # ==============================================================================
    # MODE: CONCORDANCE
    # ==============================================================================
    elif args.mode == 'concordance':
        # Get available tags
        conn = sqlite3.connect(args.index)
        cursor = conn.cursor()
        cursor.execute("SELECT DISTINCT tag_type FROM matches")
        tags = sorted([row[0] for row in cursor.fetchall()])
        conn.close()

        if len(tags) < 2:
            print("Error: Need at least 2 tag types for concordance analysis", file=sys.stderr)
            sys.exit(1)

        print(f"Tags: {', '.join(tags)}", file=sys.stderr)
        print(file=sys.stderr)

        # Load data for all tags
        print("Loading tag data...", file=sys.stderr)
        tag_maps = {}

        for tag in tags:
            print(f"  Loading {tag}...", file=sys.stderr)
            data = load_tag_data_concordance(args.index, tag)
            tag_maps[tag] = build_offset_to_isotag_map(data)
            print(f"    {len(data):,} isotags, {len(tag_maps[tag]):,} reads", file=sys.stderr)

        print(file=sys.stderr)

        # Calculate concordances
        print("Calculating pairwise concordances...", file=sys.stderr)
        concordance_matrix = calculate_all_concordances(tag_maps)
        print(f"  Computed {len(tags)*(len(tags)-1)//2} pairwise concordances", file=sys.stderr)
        print(file=sys.stderr)

        # Print matrix
        print_concordance_summary(concordance_matrix)

        # Write TSV
        if args.output:
            try:
                write_concordance_tsv(concordance_matrix, args.output)
                print(f"✓ Concordance matrix written to: {args.output}")
            except Exception as e:
                print(f"Error writing output: {e}", file=sys.stderr)
                sys.exit(1)


if __name__ == '__main__':
    main()
