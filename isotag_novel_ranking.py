#!/usr/bin/env python3
"""
isotag_novel_ranking.py - Rank Novel Isoforms by Confidence

Ranks novel isoforms based on how many tag types they're novel in.
Higher novelty score = more confident discovery.

Usage:
    python3 isotag_novel_ranking.py --index matches.db --output ranked_novel.tsv

Author: IsoTag Team
Version: 1.0.0
"""

import sqlite3
import pickle
import os
import sys
import argparse
from collections import defaultdict


def load_novel_isotags(db_path, tag_type):
    """Load isotags that are novel (only_in_a)"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute(
        "SELECT isotag, a_offsets FROM matches WHERE tag_type = ?",
        (tag_type,)
    )

    novel = {}
    for row in cursor.fetchall():
        isotag = row[0]
        a_offsets = pickle.loads(row[1])

        if len(a_offsets) > 0:  # Has reads in A
            # Check if it's truly novel (not in B)
            cursor.execute(
                "SELECT b_offsets FROM matches WHERE tag_type = ? AND isotag = ?",
                (tag_type, isotag)
            )
            result = cursor.fetchone()
            if result:
                b_offsets = pickle.loads(result[0])
                if len(b_offsets) == 0:  # Novel!
                    novel[isotag] = len(a_offsets)

    conn.close()
    return novel


def build_offset_to_isotag_map(db_path, tag_type):
    """Build mapping from offset to isotag for a tag type"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute(
        "SELECT isotag, a_offsets FROM matches WHERE tag_type = ?",
        (tag_type,)
    )

    offset_map = {}
    for row in cursor.fetchall():
        isotag = row[0]
        a_offsets = pickle.loads(row[1])
        for offset in a_offsets:
            offset_map[offset] = isotag

    conn.close()
    return offset_map


def calculate_novelty_scores(db_path, tag_types=['XI', 'XS', 'XB', 'XT']):
    """
    Calculate novelty score for each unique read/isotag combination

    Returns:
        dict: {xi_isotag: {'score': N, 'tags': {...}, 'reads': N}}
    """

    print("Loading novel isotags for each tag type...", file=sys.stderr)

    # Load novel isotags for each tag
    novel_by_tag = {}
    for tag in tag_types:
        novel = load_novel_isotags(db_path, tag)
        novel_by_tag[tag] = set(novel.keys())
        print(f"  {tag}: {len(novel):,} novel isotags", file=sys.stderr)

    # Build offset mappings
    print("Building offset mappings...", file=sys.stderr)
    offset_maps = {}
    for tag in tag_types:
        offset_maps[tag] = build_offset_to_isotag_map(db_path, tag)

    # For each XI isotag, calculate novelty score
    print("Calculating novelty scores...", file=sys.stderr)
    xi_novelty = {}

    if 'XI' not in offset_maps:
        print("Error: XI tag required for ranking", file=sys.stderr)
        return {}

    xi_isotags = load_novel_isotags(db_path, 'XI')

    for xi, read_count in xi_isotags.items():
        # Get offsets for this XI
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute(
            "SELECT a_offsets FROM matches WHERE tag_type = 'XI' AND isotag = ?",
            (xi,)
        )
        result = cursor.fetchone()
        if not result:
            conn.close()
            continue

        offsets = pickle.loads(result[0])
        conn.close()

        # Check novelty in each tag
        novelty_tags = {}
        for tag in tag_types:
            if tag not in offset_maps:
                continue

            # Check if these reads are novel in this tag
            novel_in_tag = False
            for offset in offsets:
                if offset in offset_maps[tag]:
                    tag_isotag = offset_maps[tag][offset]
                    if tag_isotag in novel_by_tag[tag]:
                        novel_in_tag = True
                        break

            novelty_tags[tag] = 'NEW' if novel_in_tag else 'KNOWN'

        # Calculate score (number of NEW tags)
        score = sum(1 for status in novelty_tags.values() if status == 'NEW')

        xi_novelty[xi] = {
            'score': score,
            'tags': novelty_tags,
            'reads': read_count
        }

    return xi_novelty


def classify_priority(score):
    """Classify priority based on novelty score"""
    if score >= 4:
        return "CRITICAL"
    elif score == 3:
        return "HIGH"
    elif score == 2:
        return "MEDIUM"
    elif score == 1:
        return "LOW"
    else:
        return "KNOWN"


def write_ranking_report(xi_novelty, output_path, tag_types):
    """Write ranked novel isoforms to TSV"""

    with open(output_path, 'w') as f:
        # Header
        header = ['isotag_xi', 'novelty_score', 'priority', 'reads']
        header.extend([f"{tag}_status" for tag in tag_types])
        f.write('\t'.join(header) + '\n')

        # Sort by score (highest first), then by read count
        sorted_isoforms = sorted(
            xi_novelty.items(),
            key=lambda x: (x[1]['score'], x[1]['reads']),
            reverse=True
        )

        for xi, data in sorted_isoforms:
            priority = classify_priority(data['score'])

            row = [
                xi,
                str(data['score']),
                priority,
                str(data['reads'])
            ]

            for tag in tag_types:
                status = data['tags'].get(tag, 'N/A')
                row.append(status)

            f.write('\t'.join(row) + '\n')


def print_summary(xi_novelty, tag_types):
    """Print ranking summary"""

    if not xi_novelty:
        print("No novel isoforms found")
        return

    print()
    print("="*80)
    print("Novel Isoform Ranking Summary")
    print("="*80)
    print()

    # Score distribution
    score_dist = defaultdict(int)
    score_reads = defaultdict(int)

    for data in xi_novelty.values():
        score = data['score']
        score_dist[score] += 1
        score_reads[score] += data['reads']

    print("Novelty score distribution:")
    print(f"{'Score':<10} {'Priority':<12} {'Isoforms':>12} {'Reads':>15}")
    print("-"*80)

    for score in sorted(score_dist.keys(), reverse=True):
        count = score_dist[score]
        reads = score_reads[score]
        priority = classify_priority(score)
        print(f"{score:<10} {priority:<12} {count:>12,} {reads:>15,}")

    print()

    # Top candidates
    sorted_isoforms = sorted(
        xi_novelty.items(),
        key=lambda x: (x[1]['score'], x[1]['reads']),
        reverse=True
    )[:20]

    print("Top 20 novel isoforms (by novelty score and read count):")
    print(f"{'Rank':<6} {'XI_Isotag':<35} {'Score':>7} {'Reads':>10} {'Priority':<10}")
    print("-"*80)

    for rank, (xi, data) in enumerate(sorted_isoforms, 1):
        xi_short = xi[:32] + "..." if len(xi) > 35 else xi
        priority = classify_priority(data['score'])
        print(f"{rank:<6} {xi_short:<35} {data['score']:>7} {data['reads']:>10,} {priority:<10}")

    print()
    print("="*80)
    print("INTERPRETATION")
    print("="*80)

    # Completely novel
    completely_novel = sum(1 for d in xi_novelty.values() if d['score'] == len(tag_types))
    if completely_novel > 0:
        pct = completely_novel / len(xi_novelty) * 100
        print(f"• Completely novel (all {len(tag_types)} tags): {completely_novel:,} isoforms ({pct:.1f}%)")
        print("  → These are the highest confidence discoveries")
        print("  → May represent novel genes or transcripts")

    high_conf = sum(1 for d in xi_novelty.values() if d['score'] >= 3)
    if high_conf > 0:
        pct = high_conf / len(xi_novelty) * 100
        print(f"• High confidence (≥3 tags): {high_conf:,} isoforms ({pct:.1f}%)")
        print("  → Strong evidence for novelty")
        print("  → Prioritize these for validation")

    low_conf = sum(1 for d in xi_novelty.values() if d['score'] == 1)
    if low_conf > 0:
        pct = low_conf / len(xi_novelty) * 100
        print(f"• Low confidence (1 tag only): {low_conf:,} isoforms ({pct:.1f}%)")
        print("  → May be technical artifacts or rare variants")
        print("  → Lower priority for validation")

    print()
    print("RECOMMENDATIONS:")
    print("-"*80)
    print("• Focus on novelty score ≥3 (high confidence)")
    print("• Validate top 10-20 by read count")
    print("• Check if completely novel isoforms match unannotated loci")
    print("• Use XI isotag to retrieve exact exon structures")
    print("• Compare across samples to confirm reproducibility")
    print("="*80)
    print()


def main():
    """Main entry point"""

    parser = argparse.ArgumentParser(
        description='isotag_novel_ranking - Rank novel isoforms by confidence',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  isotag_novel_ranking.py --index matches.db --output ranked_novel.tsv

  # Custom tag types
  isotag_novel_ranking.py --index matches.db --tags XI,XS,XB --output ranked.tsv

Output format (TSV):
  • isotag_xi: XI isotag value
  • novelty_score: Number of tags where this is novel (0-4+)
  • priority: CRITICAL/HIGH/MEDIUM/LOW/KNOWN
  • reads: Read count
  • *_status: NEW or KNOWN for each tag

Novelty scoring:
  • Score 4: Novel in all tags (completely novel)
  • Score 3: Novel in 3 tags (high confidence)
  • Score 2: Novel in 2 tags (medium confidence)
  • Score 1: Novel in 1 tag only (low confidence)
  • Score 0: Known in all tags (not novel)

Use cases:
  • Prioritize novel isoforms for validation
  • Filter technical artifacts (score 1)
  • Find highest confidence discoveries (score ≥3)
  • Guide experimental validation efforts
        """
    )

    parser.add_argument('--index', required=True, metavar='FILE',
                       help='Input index database (.db)')
    parser.add_argument('--output', metavar='FILE',
                       help='Output TSV file (optional)')
    parser.add_argument('--tags', default='XI,XS,XB,XT',
                       help='Comma-separated tag types (default: XI,XS,XB,XT)')
    parser.add_argument('--min-score', type=int, default=1, metavar='N',
                       help='Minimum novelty score to report (default: 1)')
    parser.add_argument('-v', '--version', action='version',
                       version='isotag_novel_ranking 1.0.0')

    args = parser.parse_args()

    # Validate
    if not os.path.exists(args.index):
        print(f"Error: Index file not found: {args.index}", file=sys.stderr)
        sys.exit(1)

    # Parse tags
    tag_types = [t.strip() for t in args.tags.split(',')]

    print("="*80, file=sys.stderr)
    print("Novel Isoform Ranking", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(f"Index: {args.index}", file=sys.stderr)
    print(f"Tags:  {', '.join(tag_types)}", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(file=sys.stderr)

    # Calculate novelty scores
    try:
        xi_novelty = calculate_novelty_scores(args.index, tag_types)
    except Exception as e:
        print(f"Error calculating novelty scores: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    if not xi_novelty:
        print("No novel isoforms found", file=sys.stderr)
        sys.exit(0)

    # Filter by minimum score
    if args.min_score > 1:
        xi_novelty = {k: v for k, v in xi_novelty.items() if v['score'] >= args.min_score}
        print(f"Filtered to {len(xi_novelty):,} isoforms with score ≥{args.min_score}", file=sys.stderr)
        print(file=sys.stderr)

    # Print summary
    print_summary(xi_novelty, tag_types)

    # Write TSV
    if args.output:
        try:
            write_ranking_report(xi_novelty, args.output, tag_types)
            print(f"✓ Ranked novel isoforms written to: {args.output}")
        except Exception as e:
            print(f"Error writing output: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    main()
