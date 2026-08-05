#!/usr/bin/env python3
"""
isotag_bam_filter.py - Unified BAM Filtering Tool

Replaces and combines two separate tools:
  - isotag_multi_extract.py (Boolean logic: require/exclude tags)
  - isotag_confidence_filter.py (Confidence scoring: count matching tags)

Supports both filtering strategies in one unified interface:
  1. Boolean filtering: --require XI,XS --exclude XT
  2. Confidence scoring: --min-confidence 3
  3. Combined: --min-confidence 2 --require XI --exclude XT

Usage:
    # Boolean filtering (require multiple tags)
    python3 isotag_bam_filter.py --index matches.db --require XI,XS \
        --bam-a exp.bam -o high_conf.bam

    # Confidence scoring (minimum N tags must match)
    python3 isotag_bam_filter.py --index matches.db --min-confidence 3 \
        --bam-a exp.bam -o confident.bam

    # Combined filtering (confidence + Boolean logic)
    python3 isotag_bam_filter.py --index matches.db --min-confidence 2 --require XI \
        --exclude XT --bam-a exp.bam -o filtered.bam

    # Extract novel isoforms only (status=only_a)
    python3 isotag_bam_filter.py --index matches.db --min-confidence 3 \
        --status only_a --bam-a exp.bam -o novel_high_conf.bam

Author: IsoTag Team
Version: 2.0.0 (Unified tool)
"""

import sqlite3
import pickle
import pysam
import os
import sys
import argparse
from collections import defaultdict


# ==============================================================================
# CORE: Load tag matches
# ==============================================================================

def load_tag_matches(db_path, tag_type, status='in_both'):
    """Load isotag matches for specific tag and status"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute(
        "SELECT isotag, a_offsets, b_offsets FROM matches WHERE tag_type = ?",
        (tag_type,)
    )

    matches = {}
    for row in cursor.fetchall():
        isotag = row[0]
        a_offsets = pickle.loads(row[1])
        b_offsets = pickle.loads(row[2])

        a_count = len(a_offsets)
        b_count = len(b_offsets)

        # Filter by status
        if status == 'in_both' and a_count > 0 and b_count > 0:
            matches[isotag] = a_offsets
        elif status == 'only_a' and a_count > 0 and b_count == 0:
            matches[isotag] = a_offsets
        elif status == 'only_b' and a_count == 0 and b_count > 0:
            matches[isotag] = []  # No A offsets
        elif status == 'any':
            if a_count > 0:
                matches[isotag] = a_offsets

    conn.close()
    return matches


def load_all_tag_matches(db_path, status='in_both'):
    """Load matches for all tag types"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("SELECT DISTINCT tag_type FROM matches")
    tag_types = [row[0] for row in cursor.fetchall()]

    all_matches = {}

    for tag in tag_types:
        cursor.execute(
            "SELECT isotag, a_offsets, b_offsets FROM matches WHERE tag_type = ?",
            (tag,)
        )

        tag_matches = {}
        for row in cursor.fetchall():
            isotag = row[0]
            a_offsets = pickle.loads(row[1])
            b_offsets = pickle.loads(row[2])

            a_count = len(a_offsets)
            b_count = len(b_offsets)

            # Filter by status
            if status == 'in_both' and a_count > 0 and b_count > 0:
                tag_matches[isotag] = (a_offsets, b_offsets)
            elif status == 'only_a' and a_count > 0 and b_count == 0:
                tag_matches[isotag] = (a_offsets, b_offsets)
            elif status == 'only_b' and a_count == 0 and b_count > 0:
                tag_matches[isotag] = (a_offsets, b_offsets)
            elif status == 'any':
                tag_matches[isotag] = (a_offsets, b_offsets)

        all_matches[tag] = tag_matches

    conn.close()
    return all_matches


# ==============================================================================
# STRATEGY 1: Boolean Filtering (replaces isotag_multi_extract.py)
# ==============================================================================

def build_offset_sets(tag_matches_dict):
    """Build set of offsets for each tag"""
    offset_sets = {}

    for tag, matches in tag_matches_dict.items():
        offsets = set()
        for isotag, offset_list in matches.items():
            offsets.update(offset_list)
        offset_sets[tag] = offsets

    return offset_sets


def apply_boolean_logic(offset_sets, require_tags, exclude_tags):
    """
    Apply Boolean logic to filter offsets

    Logic:
        (require_tag1 AND require_tag2 AND ...) AND NOT (exclude_tag1 OR exclude_tag2 OR ...)
    """

    # Start with all offsets if no requirements
    if require_tags:
        # Intersection of all required tags
        passing = offset_sets[require_tags[0]].copy()
        for tag in require_tags[1:]:
            passing &= offset_sets[tag]
    else:
        # Union of all tags
        passing = set()
        for offsets in offset_sets.values():
            passing |= offsets

    # Subtract excluded tags
    if exclude_tags:
        for tag in exclude_tags:
            if tag in offset_sets:
                passing -= offset_sets[tag]

    return passing


# ==============================================================================
# STRATEGY 2: Confidence Scoring (replaces isotag_confidence_filter.py)
# ==============================================================================

def calculate_confidence_scores(all_matches, tag_priority=['XI', 'XS', 'XB', 'XT']):
    """
    Calculate confidence score for each unique isotag combination

    Returns:
        dict: {offset: {'confidence': N, 'tags': {tag: isotag}}}
    """
    # Build mapping of read offsets to isotag values for each tag
    offset_to_tags = defaultdict(lambda: {})

    for tag in tag_priority:
        if tag not in all_matches:
            continue

        for isotag, (a_offsets, b_offsets) in all_matches[tag].items():
            for offset in a_offsets:
                if tag not in offset_to_tags[offset]:
                    offset_to_tags[offset][tag] = isotag

    # Calculate confidence scores
    offset_scores = {}

    for offset, tags in offset_to_tags.items():
        confidence = len(tags)
        offset_scores[offset] = {
            'confidence': confidence,
            'tags': tags
        }

    return offset_scores


def filter_by_confidence(offset_scores, min_confidence=2):
    """Filter offsets by minimum confidence score"""
    passing_offsets = set()

    for offset, data in offset_scores.items():
        if data['confidence'] >= min_confidence:
            passing_offsets.add(offset)

    return passing_offsets


# ==============================================================================
# COMBINED: Confidence + Boolean
# ==============================================================================

def filter_combined(offset_scores, min_confidence, require_tags, exclude_tags):
    """
    Apply both confidence and Boolean filtering

    Args:
        offset_scores: Output from calculate_confidence_scores()
        min_confidence: Minimum number of matching tags
        require_tags: List of tags that must be present
        exclude_tags: List of tags that must NOT be present

    Returns:
        set: Offsets that pass all filters
    """
    passing_offsets = set()

    for offset, data in offset_scores.items():
        # Check minimum confidence
        if data['confidence'] < min_confidence:
            continue

        tags = data['tags']

        # Check required tags
        if require_tags:
            if not all(tag in tags for tag in require_tags):
                continue

        # Check excluded tags
        if exclude_tags:
            if any(tag in tags for tag in exclude_tags):
                continue

        # Passed all filters
        passing_offsets.add(offset)

    return passing_offsets


# ==============================================================================
# BAM EXTRACTION
# ==============================================================================

def extract_reads_to_bam(input_bam_path, output_bam_path, passing_offsets):
    """Extract reads at specific offsets to output BAM"""

    input_bam = pysam.AlignmentFile(input_bam_path, 'rb')
    output_bam = pysam.AlignmentFile(output_bam_path, 'wb', template=input_bam)

    total_reads = 0
    extracted_reads = 0

    while True:
        offset = input_bam.tell()

        try:
            read = next(input_bam)
        except StopIteration:
            break

        total_reads += 1

        if total_reads % 10000 == 0:
            print(f"  Processed {total_reads:,} reads, extracted {extracted_reads:,}...", end='\r', file=sys.stderr)

        if offset in passing_offsets:
            output_bam.write(read)
            extracted_reads += 1

    input_bam.close()
    output_bam.close()

    print(f"  Processed {total_reads:,} reads, extracted {extracted_reads:,} - Done!     ", file=sys.stderr)

    return total_reads, extracted_reads


# ==============================================================================
# SUMMARY & REPORTING
# ==============================================================================

def print_filter_summary(mode, all_matches, passing_offsets, min_confidence,
                        require_tags, exclude_tags, status):
    """Print summary of filtering results"""

    print()
    print("="*80)
    print(f"BAM Filtering Summary - Mode: {mode}")
    print("="*80)
    print()

    print("Available tag types:")
    for tag in sorted(all_matches.keys()):
        print(f"  {tag}: {len(all_matches[tag]):,} isotags")
    print()

    print("Filter criteria:")
    print(f"  Status:         {status}")
    if min_confidence:
        print(f"  Min confidence: {min_confidence} tags")
    if require_tags:
        print(f"  Required tags:  {', '.join(require_tags)}")
    if exclude_tags:
        print(f"  Excluded tags:  {', '.join(exclude_tags)}")
    print()

    # Boolean logic explanation
    if mode == 'combined':
        print("Combined logic:")
        parts = []
        if min_confidence:
            parts.append(f"confidence ≥ {min_confidence}")
        if require_tags:
            parts.append(f"({' AND '.join(require_tags)})")
        if exclude_tags:
            parts.append(f"NOT ({' OR '.join(exclude_tags)})")
        print(f"  {' AND '.join(parts)}")
    elif mode == 'boolean' and require_tags and exclude_tags:
        req_str = ' AND '.join(require_tags)
        exc_str = ' OR '.join(exclude_tags)
        print(f"Boolean logic: ({req_str}) AND NOT ({exc_str})")
    elif mode == 'boolean' and require_tags:
        print(f"Boolean logic: {' AND '.join(require_tags)}")
    elif mode == 'boolean' and exclude_tags:
        print(f"Boolean logic: NOT ({' OR '.join(exclude_tags)})")
    print()

    print(f"Reads passing filter: {len(passing_offsets):,}")
    print()

    print("="*80)
    print()


def print_confidence_distribution(offset_scores):
    """Print confidence score distribution"""
    confidence_dist = defaultdict(int)
    for data in offset_scores.values():
        conf = data['confidence']
        confidence_dist[conf] += 1

    print("Confidence distribution:")
    print(f"{'Confidence':<15} {'Reads':>15} {'Percentage':>15}")
    print("-"*80)

    total_reads = sum(confidence_dist.values())
    for conf in sorted(confidence_dist.keys(), reverse=True):
        count = confidence_dist[conf]
        pct = count / total_reads * 100 if total_reads > 0 else 0
        print(f"{conf} tags{' '*(15-7)} {count:>15,} {pct:>14.1f}%")

    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    """Main entry point"""

    parser = argparse.ArgumentParser(
        description='isotag_bam_filter - Unified BAM filtering with Boolean and confidence strategies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This tool replaces two separate tools:
  • isotag_multi_extract.py     → Boolean logic (--require/--exclude)
  • isotag_confidence_filter.py → Confidence scoring (--min-confidence)

Examples:
  # Boolean filtering: require both XI and XS
  isotag_bam_filter.py --index m.db --require XI,XS --bam-a exp.bam -o high_conf.bam

  # Boolean filtering: require XI but exclude XS (exact only)
  isotag_bam_filter.py --index m.db --require XI --exclude XS --bam-a exp.bam -o exact_only.bam

  # Confidence scoring: at least 3 tags must match
  isotag_bam_filter.py --index m.db --min-confidence 3 --bam-a exp.bam -o confident.bam

  # Combined: confidence ≥2 AND require XI AND exclude XT
  isotag_bam_filter.py --index m.db --min-confidence 2 --require XI --exclude XT \\
      --bam-a exp.bam -o filtered.bam

  # Novel isoforms with high confidence
  isotag_bam_filter.py --index m.db --min-confidence 3 --status only_a \\
      --bam-a exp.bam -o novel_high_conf.bam

  # rRNA removal (require XI,XS for high confidence, exclude matches)
  isotag_bam_filter.py --index rRNA.db --require XI,XS --status only_a \\
      --bam-a exp.bam -o clean.bam

Boolean logic:
  • --require: AND logic (all must be present)
  • --exclude: OR logic (none can be present)
  • Combined: (req1 AND req2) AND NOT (exc1 OR exc2)

Confidence scoring:
  • Each tag match adds 1 to confidence score
  • Score 4 = matches in XI, XS, XB, XT (highest confidence)
  • Score 3 = matches in 3 tags (high confidence)
  • Score 2 = matches in 2 tags (medium confidence)
  • Score 1 = matches in 1 tag only (low confidence)

Status options:
  • in_both: Reads matching in both A and B (known)
  • only_a: Reads only in A (novel)
  • only_b: Reads only in B (not in experiment)
  • any: Any reads with the tag

Use cases:
  • High confidence filtering: --min-confidence 3
  • Remove false positives: --min-confidence 2
  • Exact matching only: --require XI --exclude XS
  • Fuzzy rescued reads: --require XS --exclude XI
  • Novel discovery: --min-confidence 3 --status only_a
  • rRNA removal: --require XI,XS --status only_a
        """
    )

    parser.add_argument('--index', required=True, metavar='FILE',
                       help='Input index database (.db)')
    parser.add_argument('--bam-a', required=True, metavar='FILE',
                       help='Input BAM file A (to filter from)')
    parser.add_argument('-o', '--output', required=True, metavar='FILE',
                       help='Output BAM file')

    # Filtering options
    parser.add_argument('--min-confidence', type=int, metavar='N',
                       help='Minimum confidence score (1-4+)')
    parser.add_argument('--require', metavar='TAGS',
                       help='Comma-separated tags that must be present (AND logic)')
    parser.add_argument('--exclude', metavar='TAGS',
                       help='Comma-separated tags to exclude (OR logic)')
    parser.add_argument('--status', default='in_both',
                       choices=['in_both', 'only_a', 'only_b', 'any'],
                       help='Match status to consider (default: in_both)')

    parser.add_argument('-v', '--version', action='version',
                       version='isotag_bam_filter 2.0.0')

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.index):
        print(f"Error: Index file not found: {args.index}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.bam_a):
        print(f"Error: BAM file not found: {args.bam_a}", file=sys.stderr)
        sys.exit(1)

    # Must have at least one criterion
    if not args.min_confidence and not args.require and not args.exclude:
        print("Error: Must specify at least one filter:", file=sys.stderr)
        print("  --min-confidence N", file=sys.stderr)
        print("  --require TAGS", file=sys.stderr)
        print("  --exclude TAGS", file=sys.stderr)
        sys.exit(1)

    # Parse tag lists
    require_tags = args.require.split(',') if args.require else []
    exclude_tags = args.exclude.split(',') if args.exclude else []

    # Validate tags
    valid_tags = {'XI', 'XB', 'XS', 'XT', 'XV', 'XC', 'XR'}
    for tag in require_tags + exclude_tags:
        if tag not in valid_tags:
            print(f"Error: Invalid tag '{tag}'. Valid: {valid_tags}", file=sys.stderr)
            sys.exit(1)

    # Can't require and exclude the same tag
    overlap = set(require_tags) & set(exclude_tags)
    if overlap:
        print(f"Error: Cannot both require and exclude: {', '.join(overlap)}", file=sys.stderr)
        sys.exit(1)

    # Determine mode
    if args.min_confidence and (require_tags or exclude_tags):
        mode = 'combined'
    elif args.min_confidence:
        mode = 'confidence'
    else:
        mode = 'boolean'

    print("="*80, file=sys.stderr)
    print(f"BAM Filtering - Mode: {mode.upper()}", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(f"Index:      {args.index}", file=sys.stderr)
    print(f"Input BAM:  {args.bam_a}", file=sys.stderr)
    print(f"Output BAM: {args.output}", file=sys.stderr)
    if args.min_confidence:
        print(f"Min confidence: {args.min_confidence}", file=sys.stderr)
    if require_tags:
        print(f"Require:    {', '.join(require_tags)}", file=sys.stderr)
    if exclude_tags:
        print(f"Exclude:    {', '.join(exclude_tags)}", file=sys.stderr)
    print(f"Status:     {args.status}", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(file=sys.stderr)

    # ==============================================================================
    # LOAD DATA
    # ==============================================================================

    print("Loading tag matches from index...", file=sys.stderr)

    # For confidence mode, we need all tags
    if mode == 'confidence' or mode == 'combined':
        try:
            all_matches = load_all_tag_matches(args.index, status=args.status)
            for tag in sorted(all_matches.keys()):
                print(f"  Loaded {len(all_matches[tag]):,} isotags for {tag}", file=sys.stderr)
        except Exception as e:
            print(f"Error: Cannot load index: {e}", file=sys.stderr)
            sys.exit(1)

    # For boolean-only mode, load specific tags
    else:
        all_tags = list(set(require_tags + exclude_tags))
        all_matches = {}
        for tag in all_tags:
            matches = load_tag_matches(args.index, tag, status=args.status)
            all_matches[tag] = {isotag: offsets for isotag, offsets in matches.items()}
            print(f"  Loaded {len(matches):,} isotags for {tag}", file=sys.stderr)

    if not all_matches:
        print("Error: No tag matches found in index", file=sys.stderr)
        sys.exit(1)

    print(file=sys.stderr)

    # ==============================================================================
    # APPLY FILTERS
    # ==============================================================================

    print("Applying filters...", file=sys.stderr)

    if mode == 'confidence':
        # Confidence-only filtering
        offset_scores = calculate_confidence_scores(all_matches)
        print(f"  Calculated confidence scores for {len(offset_scores):,} reads", file=sys.stderr)

        passing_offsets = filter_by_confidence(offset_scores, args.min_confidence)
        print(f"  {len(passing_offsets):,} reads pass confidence filter", file=sys.stderr)

    elif mode == 'combined':
        # Combined confidence + Boolean filtering
        offset_scores = calculate_confidence_scores(all_matches)
        print(f"  Calculated confidence scores for {len(offset_scores):,} reads", file=sys.stderr)

        passing_offsets = filter_combined(
            offset_scores,
            args.min_confidence,
            require_tags,
            exclude_tags
        )
        print(f"  {len(passing_offsets):,} reads pass combined filter", file=sys.stderr)

    else:  # mode == 'boolean'
        # Boolean-only filtering
        offset_sets = build_offset_sets(all_matches)

        passing_offsets = apply_boolean_logic(offset_sets, require_tags, exclude_tags)
        print(f"  {len(passing_offsets):,} reads pass Boolean filter", file=sys.stderr)

    print(file=sys.stderr)

    # ==============================================================================
    # PRINT SUMMARY
    # ==============================================================================

    if mode == 'confidence' or mode == 'combined':
        print_confidence_distribution(offset_scores)

    print_filter_summary(
        mode, all_matches, passing_offsets,
        args.min_confidence, require_tags, exclude_tags, args.status
    )

    # ==============================================================================
    # EXTRACT TO BAM
    # ==============================================================================

    print("Extracting reads to output BAM...", file=sys.stderr)
    try:
        total, extracted = extract_reads_to_bam(args.bam_a, args.output, passing_offsets)
    except Exception as e:
        print(f"Error: Cannot extract reads: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    print(file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(f"✓ Complete! Extracted {extracted:,} / {total:,} reads ({extracted/total*100:.1f}%)", file=sys.stderr)
    print(f"✓ Output written to: {args.output}", file=sys.stderr)
    print("="*80, file=sys.stderr)


if __name__ == '__main__':
    main()
