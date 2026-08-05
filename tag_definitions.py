#!/usr/bin/env python3
"""
IsoTag Tag Definitions - Shared Tag Infrastructure

Centralized tag definitions and utilities for all IsoTag tools.
Provides consistent tag handling across the entire toolkit.

Tags supported:
- XI:Z: Structure ID (full exon coordinates hash)
- XB:Z: Boundary tag (5'/3' ends) - REVERSIBLE
- XS:Z: Splicetag (splice junctions) - REVERSIBLE
- XT:Z: Transcript group (biological clustering)
- XV:Z: Variants (optional)
- XC:Z: Cluster ID (from clustering tools)
- XR:Z: Representative (from clustering tools)
"""

from typing import Dict, Set, Optional, List
from collections import defaultdict


# Tag type definitions
TAG_TYPES = {
    'XI': {
        'name': 'Structure ID',
        'description': 'Full exon structure hash (32 chars)',
        'type': 'hash',
        'reversible': False,
        'required': True,
        'version': 'v1.0+'
    },
    'XB': {
        'name': 'Boundary Tag',
        'description': 'Reversible 5\'/3\' boundary coordinates',
        'type': 'encoded',
        'reversible': True,
        'required': True,
        'version': 'v2.0+'
    },
    'XS': {
        'name': 'Splicetag',
        'description': 'Reversible splice junction coordinates',
        'type': 'encoded',
        'reversible': True,
        'required': False,
        'version': 'v2.0+'
    },
    'XT': {
        'name': 'Transcript Group',
        'description': 'Biological clustering by TSS/body/polyA',
        'type': 'hash',
        'reversible': False,
        'required': True,
        'version': 'v2.0+'
    },
    'XV': {
        'name': 'Variants',
        'description': 'Variant IDs (dot-separated)',
        'type': 'hash_list',
        'reversible': False,
        'required': False,
        'version': 'v1.0+'
    },
    'XC': {
        'name': 'Gene/Locus ID',
        'description': 'Location-based clustering (chr+strand+binned positions)',
        'type': 'hash',
        'reversible': False,
        'required': False,
        'version': 'v10.0+'
    },
    'XR': {
        'name': 'Representative',
        'description': 'Representative tag (XI or XS)',
        'type': 'reference',
        'reversible': False,
        'required': False,
        'version': 'v7.0+'
    }
}


# Tag categories
TAG_CATEGORIES = {
    'structure': ['XI', 'XB', 'XS'],  # Structural tags
    'clustering': ['XC', 'XR'],        # Clustering tags
    'biological': ['XT'],              # Biological grouping
    'variants': ['XV'],                # Variant tags
    'reversible': ['XB', 'XS'],        # Reversible tags
    'required_v1': ['XI'],             # Required in v1.0
    'required_v2': ['XI', 'XB', 'XT']  # Required in v2.0
}


class TagExtractor:
    """Extract tags from BAM fields"""

    @staticmethod
    def extract_all_tags(bam_fields: List[str]) -> Dict[str, str]:
        """
        Extract all IsoTag tags from BAM fields

        Args:
            bam_fields: List of BAM fields (from line.split('\\t'))

        Returns:
            Dictionary of tag_name -> tag_value
        """
        tags = {}

        # Start from field 11 (tags start after first 11 fields)
        for field in bam_fields[11:]:
            if field.startswith('XI:Z:'):
                tags['XI'] = field[5:]
            elif field.startswith('XB:Z:'):
                tags['XB'] = field[5:]
            elif field.startswith('XS:Z:'):
                tags['XS'] = field[5:]
            elif field.startswith('XT:Z:'):
                tags['XT'] = field[5:]
            elif field.startswith('XV:Z:'):
                tags['XV'] = field[5:]
            elif field.startswith('XC:Z:'):
                tags['XC'] = field[5:]
            elif field.startswith('XR:Z:'):
                tags['XR'] = field[5:]

        return tags

    @staticmethod
    def extract_tag(bam_fields: List[str], tag_name: str) -> Optional[str]:
        """
        Extract a specific tag from BAM fields

        Args:
            bam_fields: List of BAM fields
            tag_name: Tag name (e.g., 'XI', 'XS')

        Returns:
            Tag value or None if not found
        """
        tag_prefix = f'{tag_name}:Z:'

        for field in bam_fields[11:]:
            if field.startswith(tag_prefix):
                return field[5:]

        return None

    @staticmethod
    def has_tag(bam_fields: List[str], tag_name: str) -> bool:
        """Check if a tag is present"""
        return TagExtractor.extract_tag(bam_fields, tag_name) is not None


class TagValidator:
    """Validate tag presence and format"""

    @staticmethod
    def detect_tag_version(tags: Dict[str, str]) -> str:
        """
        Detect IsoTag version from tag presence

        Args:
            tags: Dictionary of extracted tags

        Returns:
            Version string ('v1.0', 'v2.0', 'unknown')
        """
        has_xi = 'XI' in tags
        has_xb = 'XB' in tags
        has_xt = 'XT' in tags

        if has_xi and has_xb and has_xt:
            return 'v2.0+'
        elif has_xi:
            return 'v1.0'
        else:
            return 'unknown'

    @staticmethod
    def validate_v2_tags(tags: Dict[str, str]) -> Dict[str, List[str]]:
        """
        Validate v2.0+ tag presence

        Returns:
            Dictionary with 'missing' and 'warnings' lists
        """
        issues = {
            'missing': [],
            'warnings': []
        }

        # Check required tags
        for tag in TAG_CATEGORIES['required_v2']:
            if tag not in tags:
                issues['missing'].append(f'{tag} tag missing (required in v2.0+)')

        # Check consistency
        if 'XC' in tags and 'XR' not in tags:
            issues['warnings'].append('XC (cluster) present but XR (representative) missing')

        if 'XR' in tags and 'XC' not in tags:
            issues['warnings'].append('XR (representative) present but XC (cluster) missing')

        return issues

    @staticmethod
    def is_multi_exon(tags: Dict[str, str]) -> bool:
        """Check if read is multi-exon (has XS tag)"""
        return 'XS' in tags and tags['XS'] != 'None'

    @staticmethod
    def is_clustered(tags: Dict[str, str]) -> bool:
        """Check if read is clustered (has XC tag)"""
        return 'XC' in tags


class TagStatistics:
    """Collect statistics about tag presence"""

    def __init__(self):
        self.tag_counts = defaultdict(int)
        self.total_reads = 0
        self.version_detected = None

    def update(self, tags: Dict[str, str]):
        """Update statistics with tags from one read"""
        self.total_reads += 1

        for tag_name in TAG_TYPES.keys():
            if tag_name in tags:
                self.tag_counts[tag_name] += 1

        # Detect version on first read
        if self.version_detected is None and tags:
            self.version_detected = TagValidator.detect_tag_version(tags)

    def get_summary(self) -> Dict[str, any]:
        """Get summary statistics"""
        if self.total_reads == 0:
            return {}

        return {
            'total_reads': self.total_reads,
            'version': self.version_detected or 'unknown',
            'tag_presence': {
                tag: {
                    'count': self.tag_counts[tag],
                    'percentage': 100 * self.tag_counts[tag] / self.total_reads
                }
                for tag in TAG_TYPES.keys()
            },
            'multi_exon_reads': self.tag_counts.get('XS', 0),
            'clustered_reads': self.tag_counts.get('XC', 0)
        }


def format_tag_for_display(tag_name: str, tag_value: Optional[str], max_length: int = 50) -> str:
    """
    Format tag value for display

    Args:
        tag_name: Tag name (e.g., 'XI')
        tag_value: Tag value or None
        max_length: Maximum display length

    Returns:
        Formatted string
    """
    if tag_value is None:
        return 'None'

    if len(tag_value) > max_length:
        return tag_value[:max_length] + '...'

    return tag_value


def get_tag_description(tag_name: str) -> str:
    """Get human-readable tag description"""
    if tag_name in TAG_TYPES:
        return TAG_TYPES[tag_name]['description']
    return 'Unknown tag'


def get_all_tag_names() -> List[str]:
    """Get list of all supported tag names"""
    return list(TAG_TYPES.keys())


def get_reversible_tags() -> List[str]:
    """Get list of reversible tag names"""
    return TAG_CATEGORIES['reversible']


def get_clustering_tags() -> List[str]:
    """Get list of clustering-related tag names"""
    return TAG_CATEGORIES['clustering']


def format_tag_table(tags_to_show: Optional[List[str]] = None, show_examples: bool = False) -> str:
    """
    Format tag definitions as a table

    Args:
        tags_to_show: List of tag names to show (None = all)
        show_examples: Include example values

    Returns:
        Formatted table string
    """
    if tags_to_show is None:
        tags_to_show = list(TAG_TYPES.keys())

    # Header
    lines = []
    lines.append("IsoTag Tag Definitions")
    lines.append("=" * 80)
    lines.append("")

    if show_examples:
        header = f"{'Tag':<4} {'Name':<20} {'Format':<25} {'Reversible':<12} {'Version':<10}"
    else:
        header = f"{'Tag':<4} {'Name':<20} {'Description':<40} {'Version':<10}"

    lines.append(header)
    lines.append("-" * 80)

    for tag in tags_to_show:
        if tag in TAG_TYPES:
            info = TAG_TYPES[tag]
            rev = "Yes" if info['reversible'] else "No"

            if show_examples:
                # Add format/example
                if tag == 'XI':
                    fmt = "32-char hash"
                elif tag == 'XB':
                    fmt = "8chr[s].hex.hex"
                elif tag == 'XS':
                    fmt = "8chr[s].hex.hex..."
                elif tag == 'XT':
                    fmt = "32-char hash"
                elif tag == 'XV':
                    fmt = "32-char.32-char..."
                elif tag == 'XC':
                    fmt = "32-char hash"
                elif tag == 'XR':
                    fmt = "XI or XS value"
                else:
                    fmt = info['type']

                line = f"{tag:<4} {info['name']:<20} {fmt:<25} {rev:<12} {info['version']:<10}"
            else:
                line = f"{tag:<4} {info['name']:<20} {info['description']:<40} {info['version']:<10}"

            lines.append(line)

    lines.append("")

    # Add examples if requested
    if show_examples:
        lines.append("Examples:")
        lines.append("-" * 80)
        examples = {
            'XI': 'fuIF7PN23g2gq9sFxqhUNGnfOCZhkQJS',
            'XB': 'aKF498dAp.3e8.1004',
            'XS': 'aKF498dAp.4b0.7d0.866.bb8',
            'XT': '266CbPqmZz8eS-EzT4xtnYtmm-SoIhnL',
            'XV': 'abc123def456ghi789jkl012.mno345pqr678',
            'XC': 'bBhGHqwiJQoVAnl-pmYtxA56uvsur7MH',
            'XR': 'dJP4XrGnm.1caf3.1cc94.1cd2d.1d00d'
        }

        for tag in tags_to_show:
            if tag in examples:
                lines.append(f"{tag}:Z:{examples[tag]}")
        lines.append("")

    # Add categories
    lines.append("Tag Categories:")
    lines.append("-" * 80)
    for category, tag_list in TAG_CATEGORIES.items():
        lines.append(f"{category:<20} {', '.join(tag_list)}")
    lines.append("")

    return '\n'.join(lines)


def format_tag_json(tags_to_show: Optional[List[str]] = None) -> str:
    """
    Format tag definitions as JSON

    Args:
        tags_to_show: List of tag names to show (None = all)

    Returns:
        JSON string
    """
    import json

    if tags_to_show is None:
        tags_to_show = list(TAG_TYPES.keys())

    output = {
        'tags': {tag: TAG_TYPES[tag] for tag in tags_to_show if tag in TAG_TYPES},
        'categories': TAG_CATEGORIES
    }

    return json.dumps(output, indent=2)


def format_tag_markdown(tags_to_show: Optional[List[str]] = None) -> str:
    """
    Format tag definitions as Markdown

    Args:
        tags_to_show: List of tag names to show (None = all)

    Returns:
        Markdown string
    """
    if tags_to_show is None:
        tags_to_show = list(TAG_TYPES.keys())

    lines = []
    lines.append("# IsoTag Tag Definitions")
    lines.append("")

    lines.append("## Tags")
    lines.append("")
    lines.append("| Tag | Name | Description | Type | Reversible | Version |")
    lines.append("|-----|------|-------------|------|-----------|---------|")

    for tag in tags_to_show:
        if tag in TAG_TYPES:
            info = TAG_TYPES[tag]
            rev = "✅" if info['reversible'] else "❌"
            lines.append(f"| {tag} | {info['name']} | {info['description']} | {info['type']} | {rev} | {info['version']} |")

    lines.append("")
    lines.append("## Categories")
    lines.append("")

    for category, tag_list in TAG_CATEGORIES.items():
        lines.append(f"- **{category}**: {', '.join(tag_list)}")

    lines.append("")

    return '\n'.join(lines)


# Export commonly used functions and classes
__all__ = [
    'TAG_TYPES',
    'TAG_CATEGORIES',
    'TagExtractor',
    'TagValidator',
    'TagStatistics',
    'format_tag_for_display',
    'get_tag_description',
    'get_all_tag_names',
    'get_reversible_tags',
    'get_clustering_tags',
    'format_tag_table',
    'format_tag_json',
    'format_tag_markdown'
]


# CLI Interface
if __name__ == '__main__':
    import click

    @click.command()
    @click.option('--tags', help='Comma-separated list of tags to show (e.g., XI,XB,XS)')
    @click.option('--format', 'output_format', type=click.Choice(['table', 'json', 'markdown']),
                  default='table', help='Output format (default: table)')
    @click.option('--show-examples', is_flag=True, help='Show example values (table format only)')
    @click.option('--output', '-o', help='Output file (default: stdout)')
    def main(tags, output_format, show_examples, output):
        """
        IsoTag Tag Definitions - Display tag information

        Shows definitions, formats, and examples for all IsoTag tags.

        Examples:
            # Show all tags
            python3 tag_definitions.py

            # Show specific tags
            python3 tag_definitions.py --tags XI,XB,XS

            # Export to JSON
            python3 tag_definitions.py --format json -o tags.json

            # Show with examples
            python3 tag_definitions.py --show-examples
        """

        # Parse tags if provided
        tags_to_show = None
        if tags:
            tags_to_show = [t.strip() for t in tags.split(',')]
            # Validate tags
            invalid = [t for t in tags_to_show if t not in TAG_TYPES]
            if invalid:
                click.echo(f"Error: Invalid tags: {', '.join(invalid)}", err=True)
                click.echo(f"Valid tags: {', '.join(TAG_TYPES.keys())}", err=True)
                return

        # Format output
        if output_format == 'json':
            content = format_tag_json(tags_to_show)
        elif output_format == 'markdown':
            content = format_tag_markdown(tags_to_show)
        else:  # table
            content = format_tag_table(tags_to_show, show_examples)

        # Write output
        if output:
            with open(output, 'w') as f:
                f.write(content)
            click.echo(f"✅ Tag definitions written to: {output}")
        else:
            click.echo(content)

    main()
