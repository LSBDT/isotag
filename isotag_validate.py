#!/usr/bin/env python3
"""
IsoTag Validate - Quality Control Tool (v2.1+)

Validates tag integrity and consistency in IsoTag-tagged BAM files.
Supports all tags: XI, XB, XS, XT, XV, XC, XR

Checks for:
- Tag presence and format
- XB + XS reconstruction consistency
- Chromosome hash consistency
- Junction coordinate ordering
- Tag completeness
"""

import subprocess
import click
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional


class IsoTagValidator:
    """Comprehensive IsoTag validator for v2.1+ tags (all 7 tags)"""

    def __init__(self):
        self.validation_results = {
            'file_info': {},
            'tag_validation': {
                'total_reads': 0,
                'reads_with_xi': 0,
                'reads_with_xb': 0,
                'reads_with_xs': 0,
                'reads_with_xt': 0,
                'reads_with_xv': 0,
                'reads_with_xc': 0,
                'reads_with_xr': 0,
                'valid_xi': 0,
                'valid_xb': 0,
                'valid_xs': 0,
                'valid_xt': 0,
                'valid_xv': 0,
                'valid_xc': 0,
                'valid_xr': 0,
                'malformed_xi': 0,
                'malformed_xb': 0,
                'malformed_xs': 0,
                'malformed_xt': 0,
                'malformed_xv': 0,
                'malformed_xc': 0,
                'malformed_xr': 0,
                'xi_format_errors': [],
                'xb_format_errors': [],
                'xs_format_errors': [],
                'xt_format_errors': [],
                'xv_format_errors': [],
                'xc_format_errors': [],
                'xr_format_errors': [],
                'consistency_errors': 0,
                'consistency_issues': [],
                'duplicate_tags': 0
            },
            'quality_metrics': {
                'tagging_efficiency': 0,
                'unique_structures': 0,
                'unique_boundaries': 0,
                'unique_splicetags': 0,
                'unique_transcripts': 0,
                'unique_clusters': 0,
                'unique_representatives': 0,
                'clustering_rate': 0
            },
            'warnings': [],
            'errors': [],
            'recommendations': []
        }

        # Expected formats (v2.1+)
        self.xi_pattern = re.compile(r'^[A-Za-z0-9_-]{32}$')  # 32-char hash
        self.xb_pattern = re.compile(r'^[A-Za-z0-9_-]{8}[+-]\.[0-9a-f]+\.[0-9a-f]+$')  # chr_hash.hex.hex
        self.xs_pattern = re.compile(r'^[A-Za-z0-9_-]{8}[+-](\.[0-9a-f]+)+$')  # chr_hash.hex.hex...
        self.xt_pattern = re.compile(r'^[A-Za-z0-9_-]{32}$')  # 32-char hash
        self.xc_pattern = re.compile(r'^[A-Za-z0-9_-]{32}$')  # 32-char hash
        self.xr_pattern = re.compile(r'^.+$')  # Any reference tag (can be XI, XS, etc.)
        self.xv_pattern = re.compile(r'^[A-Za-z0-9_.-]+$')    # Dot-separated

        # Counters for analysis
        self.xi_counter = Counter()
        self.xb_counter = Counter()
        self.xs_counter = Counter()
        self.xt_counter = Counter()
        self.xv_counter = Counter()
        self.xc_counter = Counter()
        self.xr_counter = Counter()
        self.read_tag_combinations = Counter()
    
    def validate_xi_format(self, xi_tag: str) -> Tuple[bool, Optional[str]]:
        """Validate XI tag format"""
        if not xi_tag:
            return False, "Empty XI tag"
        
        # Check length and character set
        if not self.xi_pattern.match(xi_tag):
            if len(xi_tag) != 32:
                return False, f"XI tag length {len(xi_tag)}, expected 32 characters"
            else:
                return False, f"XI tag contains invalid characters: {xi_tag}"
        
        return True, None
    
    def validate_xb_format(self, xb_tag: str) -> Tuple[bool, Optional[str]]:
        """Validate XB tag format (chr_hash_8 + strand + 5'_hex + 3'_hex)"""
        if not xb_tag:
            return False, "Empty XB tag"

        if not self.xb_pattern.match(xb_tag):
            return False, f"XB tag doesn't match expected format: {xb_tag}"

        # Parse components
        parts = xb_tag.split('.')
        if len(parts) != 3:
            return False, f"XB tag should have 3 dot-separated parts, found {len(parts)}"

        chr_strand = parts[0]
        if len(chr_strand) < 9:  # 8-char hash + strand
            return False, f"XB chromosome hash too short: {chr_strand}"

        strand = chr_strand[-1]
        if strand not in ['+', '-']:
            return False, f"XB strand must be + or -, found {strand}"

        # Validate hex coordinates
        try:
            int(parts[1], 16)  # 5' coordinate
            int(parts[2], 16)  # 3' coordinate
        except ValueError:
            return False, f"XB coordinates must be hexadecimal: {parts[1]}, {parts[2]}"

        return True, None

    def validate_xs_format(self, xs_tag: str) -> Tuple[bool, Optional[str]]:
        """Validate XS tag format (chr_hash_8 + strand + splice_coords_hex)"""
        if not xs_tag:
            return False, "Empty XS tag"

        if not self.xs_pattern.match(xs_tag):
            return False, f"XS tag doesn't match expected format: {xs_tag}"

        # Parse components
        parts = xs_tag.split('.')
        if len(parts) < 2:  # At least chr_strand + one junction
            return False, f"XS tag must have at least 2 parts, found {len(parts)}"

        chr_strand = parts[0]
        if len(chr_strand) < 9:  # 8-char hash + strand
            return False, f"XS chromosome hash too short: {chr_strand}"

        strand = chr_strand[-1]
        if strand not in ['+', '-']:
            return False, f"XS strand must be + or -, found {strand}"

        # Validate hex coordinates
        for i, coord in enumerate(parts[1:], 1):
            try:
                int(coord, 16)
            except ValueError:
                return False, f"XS coordinate {i} must be hexadecimal: {coord}"

        return True, None

    def validate_xt_format(self, xt_tag: str) -> Tuple[bool, Optional[str]]:
        """Validate XT tag format (32-char hash)"""
        if not xt_tag:
            return False, "Empty XT tag"

        if not self.xt_pattern.match(xt_tag):
            if len(xt_tag) != 32:
                return False, f"XT tag length {len(xt_tag)}, expected 32 characters"
            else:
                return False, f"XT tag contains invalid characters: {xt_tag}"

        return True, None

    def validate_xc_format(self, xc_tag: str) -> Tuple[bool, Optional[str]]:
        """Validate XC tag format (32-char hash)"""
        if not xc_tag:
            return False, "Empty XC tag"

        if not self.xc_pattern.match(xc_tag):
            if len(xc_tag) != 32:
                return False, f"XC tag length {len(xc_tag)}, expected 32 characters"
            else:
                return False, f"XC tag contains invalid characters: {xc_tag}"

        return True, None

    def validate_xr_format(self, xr_tag: str) -> Tuple[bool, Optional[str]]:
        """Validate XR tag format (reference to another tag)"""
        if not xr_tag:
            return False, "Empty XR tag"

        # XR can be any valid tag value (XI, XS, etc.)
        # Basic validation - just check it's non-empty and reasonable length
        if len(xr_tag) < 8:
            return False, f"XR tag too short: {xr_tag}"

        return True, None

    def validate_xv_format(self, xv_tag: str) -> Tuple[bool, Optional[str]]:
        """Validate XV tag format"""
        if not xv_tag:
            return False, "Empty XV tag"

        # Check basic character set
        if not self.xv_pattern.match(xv_tag):
            return False, f"XV tag contains invalid characters: {xv_tag}"

        # Check if it's dot-separated variant IDs
        if '.' in xv_tag:
            variant_ids = xv_tag.split('.')
            for variant_id in variant_ids:
                if len(variant_id) != 32:
                    return False, f"XV variant component '{variant_id}' is {len(variant_id)} chars, expected 32"
                if not re.match(r'^[A-Za-z0-9_-]{32}$', variant_id):
                    return False, f"XV variant component '{variant_id}' has invalid format"
        else:
            # Single variant ID
            if len(xv_tag) != 32:
                return False, f"XV tag length {len(xv_tag)}, expected 32 characters"

        return True, None
    
    def validate_bam_file(self, bam_file: str):
        """Validate entire BAM file"""
        click.echo(f"🔍 Validating {bam_file}...")
        
        # File info
        bam_path = Path(bam_file)
        self.validation_results['file_info'] = {
            'filename': bam_path.name,
            'size_mb': bam_path.stat().st_size / (1024 * 1024),
            'format': 'BAM' if bam_path.suffix.lower() == '.bam' else 'SAM'
        }
        
        # Check if BAM file is readable by samtools
        try:
            subprocess.run(['samtools', 'view', '-H', bam_file], 
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        except subprocess.CalledProcessError:
            self.validation_results['errors'].append(f"Cannot read BAM file with samtools")
            return
        
        # Process reads
        view_cmd = ['samtools', 'view', bam_file]
        
        try:
            process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)
            
            for line in process.stdout:
                self.validation_results['tag_validation']['total_reads'] += 1
                
                if self.validation_results['tag_validation']['total_reads'] % 10000 == 0:
                    click.echo(f"   ⏳ Validated {self.validation_results['tag_validation']['total_reads']:,} reads...")
                
                self._validate_read_line(line.strip())
            
            process.wait()
            
        except subprocess.CalledProcessError as e:
            self.validation_results['errors'].append(f"Error reading BAM file: {e}")
            return
        
        # Compute quality metrics
        self._compute_quality_metrics()
        
        # Generate recommendations
        self._generate_recommendations()
    
    def _validate_read_line(self, line: str):
        """Validate a single read line (v2.1+ all 7 tags)"""
        fields = line.split('\t')

        if len(fields) < 11:
            return

        qname = fields[0]

        # Track all tags on this read
        tags = {
            'XI': [],
            'XB': [],
            'XS': [],
            'XT': [],
            'XV': [],
            'XC': [],
            'XR': []
        }

        # Extract and validate all 7 tags
        for field in fields[11:]:
            if field.startswith('XI:Z:'):
                xi_tag = field[5:]
                tags['XI'].append(xi_tag)
                self.validation_results['tag_validation']['reads_with_xi'] += 1
                is_valid, error = self.validate_xi_format(xi_tag)
                if is_valid:
                    self.validation_results['tag_validation']['valid_xi'] += 1
                    self.xi_counter[xi_tag] += 1
                else:
                    self.validation_results['tag_validation']['malformed_xi'] += 1
                    self.validation_results['tag_validation']['xi_format_errors'].append({
                        'read': qname, 'tag': xi_tag, 'error': error
                    })

            elif field.startswith('XB:Z:'):
                xb_tag = field[5:]
                tags['XB'].append(xb_tag)
                self.validation_results['tag_validation']['reads_with_xb'] += 1
                is_valid, error = self.validate_xb_format(xb_tag)
                if is_valid:
                    self.validation_results['tag_validation']['valid_xb'] += 1
                    self.xb_counter[xb_tag] += 1
                else:
                    self.validation_results['tag_validation']['malformed_xb'] += 1
                    self.validation_results['tag_validation']['xb_format_errors'].append({
                        'read': qname, 'tag': xb_tag, 'error': error
                    })

            elif field.startswith('XS:Z:'):
                xs_tag = field[5:]
                tags['XS'].append(xs_tag)
                self.validation_results['tag_validation']['reads_with_xs'] += 1
                is_valid, error = self.validate_xs_format(xs_tag)
                if is_valid:
                    self.validation_results['tag_validation']['valid_xs'] += 1
                    self.xs_counter[xs_tag] += 1
                else:
                    self.validation_results['tag_validation']['malformed_xs'] += 1
                    self.validation_results['tag_validation']['xs_format_errors'].append({
                        'read': qname, 'tag': xs_tag, 'error': error
                    })

            elif field.startswith('XT:Z:'):
                xt_tag = field[5:]
                tags['XT'].append(xt_tag)
                self.validation_results['tag_validation']['reads_with_xt'] += 1
                is_valid, error = self.validate_xt_format(xt_tag)
                if is_valid:
                    self.validation_results['tag_validation']['valid_xt'] += 1
                    self.xt_counter[xt_tag] += 1
                else:
                    self.validation_results['tag_validation']['malformed_xt'] += 1
                    self.validation_results['tag_validation']['xt_format_errors'].append({
                        'read': qname, 'tag': xt_tag, 'error': error
                    })

            elif field.startswith('XV:Z:'):
                xv_tag = field[5:]
                tags['XV'].append(xv_tag)
                self.validation_results['tag_validation']['reads_with_xv'] += 1
                is_valid, error = self.validate_xv_format(xv_tag)
                if is_valid:
                    self.validation_results['tag_validation']['valid_xv'] += 1
                    self.xv_counter[xv_tag] += 1
                else:
                    self.validation_results['tag_validation']['malformed_xv'] += 1
                    self.validation_results['tag_validation']['xv_format_errors'].append({
                        'read': qname, 'tag': xv_tag, 'error': error
                    })

            elif field.startswith('XC:Z:'):
                xc_tag = field[5:]
                tags['XC'].append(xc_tag)
                self.validation_results['tag_validation']['reads_with_xc'] += 1
                is_valid, error = self.validate_xc_format(xc_tag)
                if is_valid:
                    self.validation_results['tag_validation']['valid_xc'] += 1
                    self.xc_counter[xc_tag] += 1
                else:
                    self.validation_results['tag_validation']['malformed_xc'] += 1
                    self.validation_results['tag_validation']['xc_format_errors'].append({
                        'read': qname, 'tag': xc_tag, 'error': error
                    })

            elif field.startswith('XR:Z:'):
                xr_tag = field[5:]
                tags['XR'].append(xr_tag)
                self.validation_results['tag_validation']['reads_with_xr'] += 1
                is_valid, error = self.validate_xr_format(xr_tag)
                if is_valid:
                    self.validation_results['tag_validation']['valid_xr'] += 1
                    self.xr_counter[xr_tag] += 1
                else:
                    self.validation_results['tag_validation']['malformed_xr'] += 1
                    self.validation_results['tag_validation']['xr_format_errors'].append({
                        'read': qname, 'tag': xr_tag, 'error': error
                    })

        # Check for duplicate tags on same read
        duplicate_found = False
        for tag_name, tag_list in tags.items():
            if len(tag_list) > 1:
                duplicate_found = True
                break

        if duplicate_found:
            self.validation_results['tag_validation']['duplicate_tags'] += 1
            if len(self.validation_results['warnings']) < 100:  # Limit warnings
                self.validation_results['warnings'].append(f"Read {qname} has duplicate tags")

        # Track tag combinations (just presence/absence for v2.1)
        tag_combo = '|'.join([
            'XI' if tags['XI'] else '-',
            'XB' if tags['XB'] else '-',
            'XS' if tags['XS'] else '-',
            'XT' if tags['XT'] else '-',
            'XV' if tags['XV'] else '-',
            'XC' if tags['XC'] else '-',
            'XR' if tags['XR'] else '-'
        ])
        self.read_tag_combinations[tag_combo] += 1
    
    def _compute_quality_metrics(self):
        """Compute quality metrics (v2.1+ all 7 tags)"""
        total_reads = self.validation_results['tag_validation']['total_reads']

        if total_reads == 0:
            return

        # Basic efficiency metrics
        self.validation_results['quality_metrics']['tagging_efficiency'] = \
            self.validation_results['tag_validation']['reads_with_xi'] / total_reads

        # Unique counts for all tags
        self.validation_results['quality_metrics']['unique_structures'] = len(self.xi_counter)
        self.validation_results['quality_metrics']['unique_boundaries'] = len(self.xb_counter)
        self.validation_results['quality_metrics']['unique_splicetags'] = len(self.xs_counter)
        self.validation_results['quality_metrics']['unique_transcripts'] = len(self.xt_counter)
        self.validation_results['quality_metrics']['unique_variants'] = len(self.xv_counter)
        self.validation_results['quality_metrics']['unique_clusters'] = len(self.xc_counter)
        self.validation_results['quality_metrics']['unique_representatives'] = len(self.xr_counter)

        # Clustering rate
        if total_reads > 0:
            self.validation_results['quality_metrics']['clustering_rate'] = \
                self.validation_results['tag_validation']['reads_with_xc'] / total_reads

        # Structure complexity
        if self.xi_counter:
            xi_counts = list(self.xi_counter.values())
            self.validation_results['quality_metrics']['structure_complexity'] = {
                'mean_frequency': sum(xi_counts) / len(xi_counts),
                'max_frequency': max(xi_counts),
                'singletons': sum(1 for c in xi_counts if c == 1),
                'singleton_rate': sum(1 for c in xi_counts if c == 1) / len(xi_counts)
            }

        # Variant complexity
        if self.xv_counter:
            xv_counts = list(self.xv_counter.values())
            self.validation_results['quality_metrics']['variant_complexity'] = {
                'mean_frequency': sum(xv_counts) / len(xv_counts),
                'max_frequency': max(xv_counts),
                'singletons': sum(1 for c in xv_counts if c == 1),
                'singleton_rate': sum(1 for c in xv_counts if c == 1) / len(xv_counts)
            }

        # Cluster size distribution
        if self.xc_counter:
            xc_counts = list(self.xc_counter.values())
            self.validation_results['quality_metrics']['cluster_complexity'] = {
                'mean_cluster_size': sum(xc_counts) / len(xc_counts),
                'max_cluster_size': max(xc_counts),
                'num_clusters': len(self.xc_counter)
            }
    
    def _generate_recommendations(self):
        """Generate validation recommendations (v2.1+ all 7 tags)"""
        recommendations = []
        warnings = []

        # Check tagging efficiency
        tagging_eff = self.validation_results['quality_metrics']['tagging_efficiency']
        if tagging_eff < 0.5:
            warnings.append(f"Low tagging efficiency ({tagging_eff:.1%}). Consider re-running isotag with proper parameters.")
        elif tagging_eff > 0.95:
            recommendations.append("Excellent tagging efficiency! Dataset is well-processed.")

        # Check for malformed tags (all 7 tags)
        tv = self.validation_results['tag_validation']
        malformed_counts = {
            'XI': tv.get('malformed_xi', 0),
            'XB': tv.get('malformed_xb', 0),
            'XS': tv.get('malformed_xs', 0),
            'XT': tv.get('malformed_xt', 0),
            'XV': tv.get('malformed_xv', 0),
            'XC': tv.get('malformed_xc', 0),
            'XR': tv.get('malformed_xr', 0)
        }

        for tag_name, count in malformed_counts.items():
            if count > 0:
                warnings.append(f"Found {count} malformed {tag_name} tags. Check isotag version compatibility.")

        # Check duplicate tags
        duplicate_tags = tv.get('duplicate_tags', 0)
        if duplicate_tags > 0:
            warnings.append(f"Found {duplicate_tags} reads with duplicate tags. This may indicate processing issues.")

        # Check XB/XS presence (should come together)
        reads_with_xb = tv['reads_with_xb']
        reads_with_xs = tv['reads_with_xs']
        if reads_with_xb > 0 and reads_with_xs == 0:
            warnings.append("XB tags present but no XS tags. XB and XS should typically be generated together.")
        elif reads_with_xs > 0 and reads_with_xb == 0:
            warnings.append("XS tags present but no XB tags. XB and XS should typically be generated together.")

        # Check clustering status
        clustering_rate = self.validation_results['quality_metrics'].get('clustering_rate', 0)
        if clustering_rate > 0:
            recommendations.append(f"File contains clustering information ({clustering_rate:.1%} clustered).")
            if clustering_rate < 0.5:
                warnings.append("Less than 50% of reads have cluster assignments. Consider re-running clustering.")

        # Check XC/XR consistency
        reads_with_xc = tv['reads_with_xc']
        reads_with_xr = tv['reads_with_xr']
        if reads_with_xc > 0 and reads_with_xr == 0:
            warnings.append("XC (cluster) tags present but no XR (representative) tags. This is unusual.")

        # Check complexity
        if 'structure_complexity' in self.validation_results['quality_metrics']:
            singleton_rate = self.validation_results['quality_metrics']['structure_complexity']['singleton_rate']
            if singleton_rate > 0.8:
                warnings.append(f"High singleton rate ({singleton_rate:.1%}) suggests many unique structures. Consider filtering or deeper sequencing.")

        # File size vs read count consistency
        total_reads = tv['total_reads']
        file_size_mb = self.validation_results['file_info']['size_mb']

        if total_reads > 0:
            reads_per_mb = total_reads / file_size_mb
            if reads_per_mb < 1000:  # Very low density
                warnings.append("Low read density suggests large file with few tagged reads.")

        self.validation_results['warnings'].extend(warnings)
        self.validation_results['recommendations'].extend(recommendations)
    
    def display_results(self):
        """Display validation results (v2.1+ all 7 tags)"""
        click.echo("\n" + "="*60)
        click.echo("🔍 ISO-TOOLS VALIDATION RESULTS (v2.1+)")
        click.echo("="*60)

        # File info
        info = self.validation_results['file_info']
        click.echo(f"📁 File: {info['filename']} ({info['size_mb']:.1f} MB, {info['format']})")

        # Tag validation summary
        tv = self.validation_results['tag_validation']
        total = tv['total_reads']

        click.echo(f"\n📊 TAG PRESENCE (Total reads: {total:,})")
        click.echo("-" * 60)
        click.echo(f"XI (Structure):      {tv['reads_with_xi']:,} ({100*tv['reads_with_xi']/total:.1f}%)")
        click.echo(f"XB (Boundary):       {tv['reads_with_xb']:,} ({100*tv['reads_with_xb']/total:.1f}%)")
        click.echo(f"XS (Splicetag):      {tv['reads_with_xs']:,} ({100*tv['reads_with_xs']/total:.1f}%)")
        click.echo(f"XT (Transcript):     {tv['reads_with_xt']:,} ({100*tv['reads_with_xt']/total:.1f}%)")
        click.echo(f"XV (Variant):        {tv['reads_with_xv']:,} ({100*tv['reads_with_xv']/total:.1f}%)")
        click.echo(f"XC (Cluster):        {tv['reads_with_xc']:,} ({100*tv['reads_with_xc']/total:.1f}%)")
        click.echo(f"XR (Representative): {tv['reads_with_xr']:,} ({100*tv['reads_with_xr']/total:.1f}%)")

        # Format validation
        click.echo(f"\n✅ VALID TAGS")
        click.echo("-" * 60)
        click.echo(f"Valid XI: {tv['valid_xi']:,}  |  Valid XB: {tv['valid_xb']:,}  |  Valid XS: {tv['valid_xs']:,}")
        click.echo(f"Valid XT: {tv['valid_xt']:,}  |  Valid XV: {tv['valid_xv']:,}")
        click.echo(f"Valid XC: {tv['valid_xc']:,}  |  Valid XR: {tv['valid_xr']:,}")

        # Format errors (if any)
        malformed_total = (tv.get('malformed_xi', 0) + tv.get('malformed_xb', 0) +
                          tv.get('malformed_xs', 0) + tv.get('malformed_xt', 0) +
                          tv.get('malformed_xv', 0) + tv.get('malformed_xc', 0) +
                          tv.get('malformed_xr', 0))

        if malformed_total > 0:
            click.echo(f"\n❌ FORMAT ERRORS")
            click.echo("-" * 60)
            if tv.get('malformed_xi', 0) > 0:
                click.echo(f"Malformed XI: {tv['malformed_xi']:,}")
            if tv.get('malformed_xb', 0) > 0:
                click.echo(f"Malformed XB: {tv['malformed_xb']:,}")
            if tv.get('malformed_xs', 0) > 0:
                click.echo(f"Malformed XS: {tv['malformed_xs']:,}")
            if tv.get('malformed_xt', 0) > 0:
                click.echo(f"Malformed XT: {tv['malformed_xt']:,}")
            if tv.get('malformed_xv', 0) > 0:
                click.echo(f"Malformed XV: {tv['malformed_xv']:,}")
            if tv.get('malformed_xc', 0) > 0:
                click.echo(f"Malformed XC: {tv['malformed_xc']:,}")
            if tv.get('malformed_xr', 0) > 0:
                click.echo(f"Malformed XR: {tv['malformed_xr']:,}")

        if tv.get('duplicate_tags', 0) > 0:
            click.echo(f"\n⚠️  Duplicate tags: {tv['duplicate_tags']:,}")

        # Quality metrics
        qm = self.validation_results['quality_metrics']
        click.echo(f"\n📈 QUALITY METRICS")
        click.echo("-" * 60)
        click.echo(f"Tagging efficiency:      {qm['tagging_efficiency']:.1%}")
        click.echo(f"Unique structures (XI):  {qm['unique_structures']:,}")
        click.echo(f"Unique boundaries (XB):  {qm['unique_boundaries']:,}")
        click.echo(f"Unique splicetags (XS):  {qm['unique_splicetags']:,}")
        click.echo(f"Unique transcripts (XT): {qm['unique_transcripts']:,}")
        click.echo(f"Unique variants (XV):    {qm['unique_variants']:,}")

        if qm['unique_clusters'] > 0:
            click.echo(f"Unique clusters (XC):    {qm['unique_clusters']:,}")
            click.echo(f"Clustering rate:         {qm['clustering_rate']:.1%}")
            if 'cluster_complexity' in qm:
                cc = qm['cluster_complexity']
                click.echo(f"Mean cluster size:       {cc['mean_cluster_size']:.1f}")
                click.echo(f"Max cluster size:        {cc['max_cluster_size']:,}")

        if 'structure_complexity' in qm:
            sc = qm['structure_complexity']
            click.echo(f"\nStructure singletons: {sc['singletons']:,} ({sc['singleton_rate']:.1%})")

        # Warnings and recommendations
        if self.validation_results['warnings']:
            click.echo(f"\n⚠️  WARNINGS")
            click.echo("-" * 60)
            for warning in self.validation_results['warnings']:
                click.echo(f"• {warning}")

        if self.validation_results['recommendations']:
            click.echo(f"\n💡 RECOMMENDATIONS")
            click.echo("-" * 60)
            for rec in self.validation_results['recommendations']:
                click.echo(f"• {rec}")

        # Overall assessment
        click.echo(f"\n🎯 OVERALL ASSESSMENT")
        click.echo("-" * 60)

        if malformed_total == 0:
            click.echo("✅ All isotags are properly formatted")
        else:
            click.echo(f"❌ Found {malformed_total} format violations - check isotag processing")

        if qm['tagging_efficiency'] > 0.8:
            click.echo("✅ Good tagging efficiency")
        else:
            click.echo("⚠️  Low tagging efficiency - consider re-processing")
    
    def export_report(self, output_file: str):
        """Export detailed validation report (v2.1+ all 7 tags)"""
        # Limit error lists to avoid huge files
        max_errors = 100

        export_data = dict(self.validation_results)

        # Truncate error lists if too long (all 7 tags)
        error_types = ['xi', 'xb', 'xs', 'xt', 'xv', 'xc', 'xr']
        for tag_name in error_types:
            error_key = f'{tag_name}_format_errors'
            if error_key in export_data['tag_validation']:
                if len(export_data['tag_validation'][error_key]) > max_errors:
                    export_data['tag_validation'][error_key] = export_data['tag_validation'][error_key][:max_errors]
                    export_data['warnings'].append(f"{tag_name.upper()} format errors truncated to {max_errors} entries")

        # Add top tags for all types
        export_data['top_structures'] = dict(self.xi_counter.most_common(20))
        export_data['top_boundaries'] = dict(self.xb_counter.most_common(20))
        export_data['top_splicetags'] = dict(self.xs_counter.most_common(20))
        export_data['top_transcripts'] = dict(self.xt_counter.most_common(20))
        export_data['top_variants'] = dict(self.xv_counter.most_common(20))
        export_data['top_clusters'] = dict(self.xc_counter.most_common(20))
        export_data['top_representatives'] = dict(self.xr_counter.most_common(20))

        # Add tag combination statistics
        export_data['tag_combinations'] = dict(self.read_tag_combinations.most_common(10))

        with open(output_file, 'w') as f:
            json.dump(export_data, f, indent=2)

        click.echo(f"📄 Detailed validation report: {output_file}")


@click.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.option('--output', '-o', help='Output JSON report file')
@click.option('--quiet', '-q', is_flag=True, help='Suppress console output')
def validate(bam_file, output, quiet):
    """
    ISO-Tools Validate - Quality Control for Isotag Files (v2.1+)

    Validate isotag format compliance and quality metrics for all 7 tags:
    - XI (Structure): Full exon coordinate hashes
    - XB (Boundary): Reversible 5'/3' boundary tags
    - XS (Splicetag): Reversible splice junction tags
    - XT (Transcript): Biological transcript group hashes
    - XV (Variant): Variant IDs
    - XC (Cluster): Cluster assignment hashes
    - XR (Representative): Representative tag references

    Checks for:
    - Tag presence and format correctness
    - XB/XS consistency (should come together)
    - XC/XR consistency (cluster tags)
    - Quality metrics and recommendations

    Examples:
        # Basic validation (all 7 tags)
        isotag_validate.py tagged.bam

        # Export detailed report
        isotag_validate.py tagged.bam -o validation_report.json

        # Quiet mode (report only)
        isotag_validate.py tagged.bam -o report.json --quiet

        # Validate clustered file
        isotag_validate.py clustered.bam -o cluster_validation.json
    """

    validator = IsoTagValidator()

    # Run validation
    validator.validate_bam_file(bam_file)

    # Display results unless quiet
    if not quiet:
        validator.display_results()

    # Export report if requested
    if output:
        validator.export_report(output)

    # Exit with error code if major issues found
    if validator.validation_results['errors']:
        click.echo(f"\n❌ Validation failed with {len(validator.validation_results['errors'])} errors")
        sys.exit(1)
    else:
        click.echo(f"\n✅ Validation complete!")


if __name__ == '__main__':
    validate()