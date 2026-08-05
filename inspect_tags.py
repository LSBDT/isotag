#!/usr/bin/env python3
"""
Inspect IsoTag Tags - Display All Tags in BAM File

Version: v2.0 (Multi-tag support)
Shows: XI, XB, XS, XT, XV, XC, XR tags
"""

import subprocess
import click
from collections import Counter

# Import shared tag definitions
try:
    from tag_definitions import TagExtractor, TAG_TYPES, get_all_tag_names
except ImportError:
    # Fallback if tag_definitions not available
    def extract_all_tags_fallback(fields):
        tags = {}
        for field in fields[11:]:
            if field.startswith('XI:Z:'): tags['XI'] = field[5:]
            elif field.startswith('XB:Z:'): tags['XB'] = field[5:]
            elif field.startswith('XS:Z:'): tags['XS'] = field[5:]
            elif field.startswith('XT:Z:'): tags['XT'] = field[5:]
            elif field.startswith('XV:Z:'): tags['XV'] = field[5:]
            elif field.startswith('XC:Z:'): tags['XC'] = field[5:]
            elif field.startswith('XR:Z:'): tags['XR'] = field[5:]
        return tags

    class TagExtractor:
        extract_all_tags = staticmethod(extract_all_tags_fallback)

    TAG_TYPES = {
        'XI': {'name': 'Structure ID'},
        'XB': {'name': 'Boundary'},
        'XS': {'name': 'Splicetag'},
        'XT': {'name': 'Transcript Group'},
        'XV': {'name': 'Variants'},
        'XC': {'name': 'Cluster ID'},
        'XR': {'name': 'Representative'}
    }

    def get_all_tag_names():
        return ['XI', 'XB', 'XS', 'XT', 'XV', 'XC', 'XR']


@click.command()
@click.option('--bam', required=True, help='BAM file with IsoTag tags')
@click.option('--limit', default=10, help='Number of reads to display')
@click.option('--show-only', help='Comma-separated tag names to show (e.g., XI,XS,XC)')
@click.option('--decode', is_flag=True, help='Decode XB and XS tags to coordinates')
def inspect_tags(bam, limit, show_only, decode):
    """
    Inspect IsoTag Tags - Display All Tags in BAM File

    Examples:
        # Show all tags for first 10 reads
        inspect_tags.py --bam tagged.bam

        # Show only specific tags
        inspect_tags.py --bam tagged.bam --show-only XI,XS,XC

        # Show first 20 reads
        inspect_tags.py --bam tagged.bam --limit 20

        # Decode reversible tags
        inspect_tags.py --bam tagged.bam --decode
    """

    # Parse show_only filter
    show_tags = None
    if show_only:
        show_tags = set(tag.strip().upper() for tag in show_only.split(','))

    # Count occurrences of each tag
    tag_counters = {tag: Counter() for tag in get_all_tag_names()}
    read_count = 0
    reads_with_tags = {tag: 0 for tag in get_all_tag_names()}

    # Stream through BAM
    view_cmd = ['samtools', 'view', bam]
    process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

    print(f"=== First {limit} reads with IsoTag tags ===\n")
    displayed = 0

    for line in process.stdout:
        read_count += 1
        fields = line.strip().split('\t')

        if len(fields) >= 11:
            qname = fields[0]

            # Extract all tags
            tags = TagExtractor.extract_all_tags(fields)

            # Update counters
            for tag_name, tag_value in tags.items():
                tag_counters[tag_name][tag_value] += 1
                reads_with_tags[tag_name] += 1

            # Display first few reads with tags
            if tags and displayed < limit:
                print(f"Read: {qname}")

                # Display tags in order
                for tag_name in ['XI', 'XB', 'XS', 'XT', 'XV', 'XC', 'XR']:
                    # Skip if not in show_tags filter
                    if show_tags and tag_name not in show_tags:
                        continue

                    tag_value = tags.get(tag_name)
                    tag_desc = TAG_TYPES.get(tag_name, {}).get('name', tag_name)

                    if tag_value:
                        # Truncate long values
                        display_value = tag_value if len(tag_value) <= 50 else tag_value[:47] + '...'

                        # Decode if requested
                        if decode and tag_name in ['XB', 'XS']:
                            try:
                                from decode_tags import decode_boundarytag, decode_splicetag
                                if tag_name == 'XB':
                                    chr_hash, strand, five_p, three_p = decode_boundarytag(tag_value)
                                    print(f"  {tag_name} ({tag_desc:18s}): {display_value}")
                                    print(f"      → {chr_hash} {strand} {five_p:,}-{three_p:,}")
                                elif tag_name == 'XS':
                                    chr_hash, strand, coords = decode_splicetag(tag_value)
                                    print(f"  {tag_name} ({tag_desc:18s}): {display_value}")
                                    print(f"      → {chr_hash} {strand} junctions: {coords}")
                            except:
                                print(f"  {tag_name} ({tag_desc:18s}): {display_value}")
                        else:
                            print(f"  {tag_name} ({tag_desc:18s}): {display_value}")
                    else:
                        if not show_tags:  # Only show "None" if not filtering
                            print(f"  {tag_name} ({tag_desc:18s}): None")

                print()
                displayed += 1

    process.wait()

    # Display statistics
    print(f"\n=== Tag Statistics ===")
    print(f"Total reads: {read_count:,}")
    print()

    for tag_name in ['XI', 'XB', 'XS', 'XT', 'XV', 'XC', 'XR']:
        count = reads_with_tags[tag_name]
        unique = len(tag_counters[tag_name])
        pct = 100 * count / read_count if read_count > 0 else 0
        tag_desc = TAG_TYPES.get(tag_name, {}).get('name', tag_name)
        status = "✅" if count > 0 else "  "
        print(f"{status} {tag_name} ({tag_desc:18s}): {count:8,} reads ({pct:5.1f}%), {unique:6,} unique")

    # Show most common values for each tag
    print(f"\n=== Top 3 Most Common Values Per Tag ===")
    for tag_name in ['XI', 'XB', 'XS', 'XT', 'XV', 'XC', 'XR']:
        if tag_counters[tag_name]:
            tag_desc = TAG_TYPES.get(tag_name, {}).get('name', tag_name)
            print(f"\n{tag_name} ({tag_desc}):")
            for value, count in tag_counters[tag_name].most_common(3):
                display_val = value if len(value) <= 40 else value[:37] + '...'
                print(f"  {display_val:40s}: {count:6,} reads")


if __name__ == '__main__':
    inspect_tags()
