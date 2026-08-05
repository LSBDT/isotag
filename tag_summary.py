#!/usr/bin/env python3
"""
IsoTag Tag Summary - Quick Tag Presence Checker

Quick one-liner to check which tags are present in a BAM file.
Useful for verifying tag version and completeness.
"""

import subprocess
import click
import sys
from collections import defaultdict

# Import shared tag definitions
try:
    from tag_definitions import TagExtractor, TAG_TYPES, TagValidator
except ImportError:
    # Fallback
    class TagExtractor:
        @staticmethod
        def extract_all_tags(fields):
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

    class TagValidator:
        @staticmethod
        def detect_tag_version(tags):
            if 'XI' in tags and 'XB' in tags and 'XT' in tags:
                return 'v2.0+'
            elif 'XI' in tags:
                return 'v1.0'
            return 'unknown'

    TAG_TYPES = {
        'XI': {'name': 'Structure ID', 'version': 'v1.0+'},
        'XB': {'name': 'Boundary Tag', 'version': 'v2.0+'},
        'XS': {'name': 'Splicetag', 'version': 'v2.0+'},
        'XT': {'name': 'Transcript Group', 'version': 'v2.0+'},
        'XV': {'name': 'Variants', 'version': 'v1.0+'},
        'XC': {'name': 'Cluster ID', 'version': 'v7.0+'},
        'XR': {'name': 'Representative', 'version': 'v7.0+'}
    }


@click.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.option('--sample', '-s', type=int, default=1000, help='Number of reads to sample (default: 1000)')
def tag_summary(bam_file, sample):
    """
    IsoTag Tag Summary - Quick Tag Presence Checker

    Quickly check which IsoTag tags are present in a BAM file.

    Examples:
        # Quick summary
        tag_summary.py tagged.bam

        # Sample 10,000 reads
        tag_summary.py tagged.bam --sample 10000
    """

    click.echo(f"🔍 Checking tags in: {bam_file}")
    click.echo(f"📊 Sampling: {sample:,} reads\n")

    tag_counts = defaultdict(int)
    total_reads = 0
    version_detected = None

    view_cmd = ['samtools', 'view', bam_file]

    try:
        process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)

        for line_num, line in enumerate(process.stdout):
            if line_num >= sample:
                break

            total_reads += 1
            fields = line.strip().split('\t')

            if len(fields) >= 11:
                # Extract all tags
                tags = TagExtractor.extract_all_tags(fields)

                # Detect version on first read
                if version_detected is None and tags:
                    version_detected = TagValidator.detect_tag_version(tags)

                # Count tag presence
                for tag_name in tags.keys():
                    tag_counts[tag_name] += 1

        process.kill()  # Stop after sampling

    except Exception as e:
        click.echo(f"❌ Error reading BAM: {e}")
        sys.exit(1)

    # Display results
    click.echo("="*60)
    click.echo("TAG PRESENCE SUMMARY")
    click.echo("="*60)

    # Tag presence
    for tag_name in ['XI', 'XB', 'XS', 'XT', 'XV', 'XC', 'XR']:
        count = tag_counts.get(tag_name, 0)
        pct = 100 * count / total_reads if total_reads > 0 else 0
        tag_info = TAG_TYPES.get(tag_name, {})
        tag_desc = tag_info.get('name', tag_name)
        tag_version = tag_info.get('version', '')

        if count > 0:
            status = "✅"
            detail = f"{count:6,}/{total_reads:6,} reads ({pct:5.1f}%)"
        else:
            status = "  "
            detail = "Not present"

        print(f"{status} {tag_name}:Z: ({tag_desc:18s}) - {detail:30s} [{tag_version}]")

    # Version detection
    click.echo("\n" + "="*60)
    if version_detected:
        click.echo(f"📋 Tag version detected: {version_detected}")
    else:
        click.echo(f"⚠️  Tag version: Unknown (no tags found)")

    # Quick diagnostics
    click.echo("\n" + "="*60)
    click.echo("QUICK DIAGNOSTICS")
    click.echo("="*60)

    has_xi = tag_counts.get('XI', 0) > 0
    has_xb = tag_counts.get('XB', 0) > 0
    has_xs = tag_counts.get('XS', 0) > 0
    has_xt = tag_counts.get('XT', 0) > 0
    has_xc = tag_counts.get('XC', 0) > 0

    if has_xi and has_xb and has_xt:
        click.echo("✅ Full v2.0+ tags present (XI, XB, XT)")
    elif has_xi:
        click.echo("⚠️  Only XI tag present (v1.0 format)")
    else:
        click.echo("❌ Missing required tags")

    if has_xs:
        multi_exon_pct = 100 * tag_counts['XS'] / total_reads
        click.echo(f"✅ Multi-exon transcripts: {multi_exon_pct:.1f}%")
    else:
        click.echo("⚠️  No multi-exon transcripts (no XS tags)")

    if has_xc:
        clustered_pct = 100 * tag_counts['XC'] / total_reads
        click.echo(f"✅ Clustered reads: {clustered_pct:.1f}%")
    else:
        click.echo("ℹ️  Not clustered (no XC tags)")

    click.echo()


if __name__ == '__main__':
    tag_summary()
