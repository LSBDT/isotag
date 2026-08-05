#!/usr/bin/env python3
"""
IsoTag Clustering - Multi-threaded BAM-to-BAM Processing

Clusters isotags by exon count with stable-end requirement.
Adds XC cluster tags and XR representative tags to output BAM for downstream analysis.

Usage:
    python3 isotag_clustering.py input.bam clustered.bam
    python3 isotag_clustering.py input.bam clustered.bam -j 8 -t 30
"""

import subprocess
import click
import json
import time
import os
import sys
import hashlib
import base64
import statistics
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Tuple, Set, Optional
import re
import multiprocessing


class IsotagRead:
    """Lightweight structure for isotag data"""
    def __init__(self, isotag_id: str, chrom: str, strand: str, start: int, end: int, exons: List[Tuple[int, int]]):
        self.isotag_id = isotag_id
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.end = end
        self.exons = exons
        self.exon_count = len(exons)


class ExonCountClusterer:
    """
    Exon-count-based clustering with deterioration awareness:
    - End exons: At least 1 end toward center must match (allows deterioration)
    - Middle exons: Both ends must match (preserves isoform structure)
    """
    
    def __init__(self, stable_end_tolerance=20):
        self.stable_end_tolerance = stable_end_tolerance
        self.stats = defaultdict(int)
    
    def has_stable_end(self, exon1: Tuple[int, int], exon2: Tuple[int, int], 
                      exon_position: int, total_exons: int) -> bool:
        """
        Position-aware exon matching with deterioration awareness:
        - End exons: At least 1 end toward center must match
        - Middle exons: Both ends must match
        """
        s1, e1 = exon1
        s2, e2 = exon2
        
        start_stable = abs(s1 - s2) <= self.stable_end_tolerance
        end_stable = abs(e1 - e2) <= self.stable_end_tolerance
        
        if total_exons == 1:
            # Single exon: at least one end must match (current behavior)
            return start_stable or end_stable
        elif exon_position == 0:
            # First exon: 3' end (toward center) must match - allows 5' deterioration
            return end_stable
        elif exon_position == total_exons - 1:
            # Last exon: 5' end (toward center) must match - allows 3' deterioration
            return start_stable
        else:
            # Middle exon: both ends must match - no deterioration allowed
            return start_stable and end_stable
    
    def can_cluster_together(self, isotag1: IsotagRead, isotag2: IsotagRead) -> bool:
        """Check if two isotags can be clustered together with deterioration awareness"""
        if isotag1.exon_count != isotag2.exon_count:
            return False
        
        total_exons = isotag1.exon_count
        
        for i in range(total_exons):
            exon1 = isotag1.exons[i]
            exon2 = isotag2.exons[i]
            if not self.has_stable_end(exon1, exon2, i, total_exons):
                return False
        
        return True
    
    def cluster_with_stable_ends(self, isotags: List[IsotagRead]) -> List[List[IsotagRead]]:
        """Cluster isotags with same exon count, requiring at least one stable end"""
        if not isotags:
            return []
        
        clusters = []
        unassigned = isotags.copy()
        
        while unassigned:
            seed = unassigned.pop(0)
            cluster = [seed]
            
            remaining = []
            for isotag in unassigned:
                if self.can_cluster_together(seed, isotag):
                    cluster.append(isotag)
                else:
                    remaining.append(isotag)
            
            unassigned = remaining
            clusters.append(cluster)
        
        return clusters
    
    def cluster_by_exon_count(self, isotags: List[IsotagRead]) -> List[List[IsotagRead]]:
        """Group by exon count, then cluster within each group"""
        exon_groups = defaultdict(list)
        for isotag in isotags:
            exon_groups[isotag.exon_count].append(isotag)
        
        all_clusters = []
        for exon_count, group_isotags in exon_groups.items():
            group_clusters = self.cluster_with_stable_ends(group_isotags)
            all_clusters.extend(group_clusters)
        
        return all_clusters


class ClusterIDGenerator:
    """Generate deterministic cluster IDs from isotag IDs"""
    
    @staticmethod
    def generate_cluster_id(isotag_ids: List[str]) -> str:
        """Generate 32-character cluster ID from sorted isotag IDs"""
        sorted_ids = sorted(isotag_ids)
        concatenated = '.'.join(sorted_ids)
        cluster_hash = hashlib.sha256(concatenated.encode()).digest()
        cluster_id = base64.urlsafe_b64encode(cluster_hash[:24]).decode().rstrip('=')
        return cluster_id
    
    @staticmethod
    def find_longest_isotag_representative(isotag_reads: List[IsotagRead]) -> str:
        """
        Find longest isotag as cluster representative with deterministic tie-breaking.
        
        Algorithm:
        1. Filter to isotags with maximum genomic span (end - start + 1)
        2. Sort alphabetically by isotag ID
        3. Return first one for deterministic results
        """
        if not isotag_reads:
            return None
        
        # Calculate genomic span for each isotag (end - start + 1)
        max_span = max(read.end - read.start + 1 for read in isotag_reads)
        
        # Filter to longest isotags
        longest_isotags = [
            read for read in isotag_reads 
            if (read.end - read.start + 1) == max_span
        ]
        
        # Deterministic tie-breaking: sort by isotag ID alphabetically
        longest_isotags.sort(key=lambda x: x.isotag_id)
        
        return longest_isotags[0].isotag_id


class ChromosomeProcessor:
    """Process individual chromosomes for clustering"""
    
    def __init__(self, stable_end_tolerance: int):
        self.stable_end_tolerance = stable_end_tolerance
    
    def parse_cigar_exons(self, pos: int, cigar: str) -> List[Tuple[int, int]]:
        """Parse CIGAR to extract exon positions"""
        exons = []
        current_pos = pos
        exon_start = pos
        in_exon = True
        
        pattern = r'(\d+)([MIDNSHP=X])'
        matches = re.findall(pattern, cigar)
        
        for length_str, op in matches:
            length = int(length_str)
            
            if op in ['M', '=', 'X']:
                if not in_exon:
                    exon_start = current_pos
                    in_exon = True
                current_pos += length
            elif op == 'D':
                current_pos += length
            elif op == 'I':
                pass
            elif op == 'N':
                if in_exon:
                    exons.append((exon_start, current_pos - 1))
                    in_exon = False
                current_pos += length
            elif op in ['S', 'H']:
                pass
        
        if in_exon and exon_start < current_pos:
            exons.append((exon_start, current_pos - 1))
        
        return exons
    
    def process_chromosome(self, bam_file: str, chromosome: str) -> Dict:
        """Process single chromosome and return cluster mapping"""
        cmd = f"samtools view {bam_file} {chromosome}"
        
        try:
            process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, text=True)
        except Exception as e:
            return {'chromosome': chromosome, 'error': str(e)}
        
        isotag_groups = defaultdict(list)
        valid_isotags = 0
        
        for line in process.stdout:
            if line.startswith('@'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            
            chrom = fields[2]
            if chrom == '*':
                continue
            
            pos = int(fields[3])
            cigar = fields[5]
            flag = int(fields[1])
            strand = '-' if flag & 16 else '+'
            
            # Extract isotag ID from XI tag
            isotag_id = None
            for field in fields[11:]:
                if field.startswith('XI:Z:'):
                    isotag_id = field[5:]
                    break
            
            if not isotag_id:
                continue
            
            exons = self.parse_cigar_exons(pos, cigar)
            if not exons:
                continue
            
            valid_isotags += 1
            start = exons[0][0]
            end = exons[-1][1]
            read = IsotagRead(isotag_id, chrom, strand, start, end, exons)
            
            key = (chrom, strand, read.exon_count)
            isotag_groups[key].append(read)
        
        process.wait()
        
        # Cluster each group
        clusterer = ExonCountClusterer(self.stable_end_tolerance)
        isotag_to_cluster = {}
        isotag_to_representative = {}
        cluster_stats = {'total_clusters': 0, 'total_isotags': valid_isotags}
        
        for (chrom, strand, exon_count), isotags in isotag_groups.items():
            if len(isotags) < 2:
                # Single isotag = single cluster
                if isotags:
                    cluster_id = ClusterIDGenerator.generate_cluster_id([isotags[0].isotag_id])
                    representative_id = ClusterIDGenerator.find_longest_isotag_representative(isotags)
                    isotag_to_cluster[isotags[0].isotag_id] = cluster_id
                    isotag_to_representative[isotags[0].isotag_id] = representative_id
                    cluster_stats['total_clusters'] += 1
                continue
            
            clusters = clusterer.cluster_by_exon_count(isotags)
            
            for cluster in clusters:
                isotag_ids = [read.isotag_id for read in cluster]
                cluster_id = ClusterIDGenerator.generate_cluster_id(isotag_ids)
                representative_id = ClusterIDGenerator.find_longest_isotag_representative(cluster)
                
                for isotag_id in isotag_ids:
                    isotag_to_cluster[isotag_id] = cluster_id
                    isotag_to_representative[isotag_id] = representative_id
                
                cluster_stats['total_clusters'] += 1
        
        return {
            'chromosome': chromosome,
            'isotag_to_cluster': isotag_to_cluster,
            'isotag_to_representative': isotag_to_representative,
            'stats': cluster_stats
        }


def extract_chromosomes_from_bam(bam_file: str) -> List[str]:
    """Extract chromosome list from BAM index"""
    cmd = f"samtools idxstats {bam_file}"
    
    try:
        result = subprocess.run(cmd.split(), capture_output=True, text=True, check=True)
        chromosomes = []
        
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) >= 3:
                chrom = fields[0]
                mapped_reads = int(fields[2])
                if chrom != '*' and mapped_reads > 0:
                    chromosomes.append(chrom)
        
        return chromosomes
    
    except subprocess.CalledProcessError as e:
        click.echo(f"❌ Error reading BAM index: {e}")
        sys.exit(1)


def parallel_chromosome_clustering(bam_file: str, threads: int, stable_end_tolerance: int) -> Dict:
    """Process chromosomes in parallel for clustering"""
    chromosomes = extract_chromosomes_from_bam(bam_file)
    click.echo(f"📊 Found {len(chromosomes)} chromosomes with reads")
    
    if threads == 0:
        threads = multiprocessing.cpu_count()
    
    max_threads = min(threads, len(chromosomes), multiprocessing.cpu_count())
    click.echo(f"⚡ Using {max_threads} threads for parallel processing")
    
    processor = ChromosomeProcessor(stable_end_tolerance)
    results = {}
    
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        future_to_chr = {
            executor.submit(processor.process_chromosome, bam_file, chr_name): chr_name 
            for chr_name in chromosomes
        }
        
        for future in as_completed(future_to_chr):
            chr_name = future_to_chr[future]
            try:
                chr_result = future.result()
                if 'error' not in chr_result:
                    results[chr_name] = chr_result
                    isotag_count = len(chr_result['isotag_to_cluster'])
                    cluster_count = chr_result['stats']['total_clusters']
                    click.echo(f"✅ {chr_name}: {isotag_count:,} isotags → {cluster_count:,} clusters")
                else:
                    click.echo(f"❌ {chr_name}: {chr_result['error']}")
            except Exception as exc:
                click.echo(f"❌ {chr_name}: {exc}")
    
    return results


def create_tagged_bam(input_bam: str, output_bam: str, cluster_mapping: Dict[str, str], representative_mapping: Dict[str, str]):
    """Create output BAM with XC cluster tags and XR representative tags"""
    click.echo(f"📝 Creating tagged BAM: {output_bam}")
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_bam)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    cmd_in = f"samtools view -h {input_bam}"
    cmd_out = f"samtools view -b -o {output_bam}"
    
    tagged_reads = 0
    total_reads = 0
    
    try:
        with subprocess.Popen(cmd_in.split(), stdout=subprocess.PIPE, text=True) as proc_in:
            with subprocess.Popen(cmd_out.split(), stdin=subprocess.PIPE, text=True) as proc_out:
                
                for line in proc_in.stdout:
                    if line.startswith('@'):  # Header
                        proc_out.stdin.write(line)
                        continue
                    
                    total_reads += 1
                    fields = line.strip().split('\t')
                    if len(fields) < 11:
                        proc_out.stdin.write(line)
                        continue
                    
                    # Look for XI tag to get isotag ID
                    isotag_id = None
                    for field in fields[11:]:
                        if field.startswith('XI:Z:'):
                            isotag_id = field[5:]
                            break
                    
                    # Add XC and XR tags if we have cluster mapping
                    if isotag_id and isotag_id in cluster_mapping:
                        cluster_id = cluster_mapping[isotag_id]
                        representative_id = representative_mapping[isotag_id]
                        line = line.rstrip() + f"\tXC:Z:{cluster_id}\tXR:Z:{representative_id}\n"
                        tagged_reads += 1
                    
                    proc_out.stdin.write(line)
        
        # Index the output BAM
        click.echo(f"🔗 Indexing output BAM...")
        subprocess.run(f"samtools index {output_bam}".split(), check=True)
        
        click.echo(f"✅ Tagged {tagged_reads:,}/{total_reads:,} reads with XC cluster tags and XR representative tags")
        
    except Exception as e:
        click.echo(f"❌ Error creating tagged BAM: {e}")
        sys.exit(1)


def create_json_stats(chromosome_results: Dict, json_file: str, total_time: float):
    """Create optional JSON statistics file"""
    total_isotags = sum(len(chr_result['isotag_to_cluster']) for chr_result in chromosome_results.values())
    total_clusters = sum(chr_result['stats']['total_clusters'] for chr_result in chromosome_results.values())
    
    stats = {
        'parameters': {
            'algorithm': 'exon_count_based_with_stable_ends',
            'processing_time': total_time
        },
        'summary': {
            'total_isotags': total_isotags,
            'total_clusters': total_clusters,
            'compression_ratio': (1 - total_clusters/total_isotags)*100 if total_isotags > 0 else 0,
            'chromosomes_processed': len(chromosome_results)
        },
        'by_chromosome': {
            chr_name: {
                'isotags': len(chr_result['isotag_to_cluster']),
                'clusters': chr_result['stats']['total_clusters']
            }
            for chr_name, chr_result in chromosome_results.items()
        }
    }
    
    with open(json_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    click.echo(f"📊 Statistics written to: {json_file}")


def create_tsv_output(input_bam: str, cluster_mapping: Dict[str, str], representative_mapping: Dict[str, str], tsv_file: str):
    """Create TSV output with isotag_id, cluster_id, representative_id, cluster_size, chromosome, strand, start, end"""
    click.echo(f"📋 Creating TSV output: {tsv_file}")
    
    # Count cluster sizes
    cluster_sizes = defaultdict(int)
    for cluster_id in cluster_mapping.values():
        cluster_sizes[cluster_id] += 1
    
    # Collect isotag details
    isotag_details = {}
    
    cmd = f"samtools view {input_bam}"
    try:
        process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, text=True)
        
        for line in process.stdout:
            if line.startswith('@'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            
            chrom = fields[2]
            if chrom == '*':
                continue
            
            pos = int(fields[3])
            flag = int(fields[1])
            strand = '-' if flag & 16 else '+'
            
            # Extract isotag ID from XI tag
            isotag_id = None
            for field in fields[11:]:
                if field.startswith('XI:Z:'):
                    isotag_id = field[5:]
                    break
            
            if not isotag_id or isotag_id not in cluster_mapping:
                continue
            
            # Parse CIGAR to get exon positions for start/end
            cigar = fields[5]
            processor = ChromosomeProcessor(20)  # Use same tolerance as main clustering
            exons = processor.parse_cigar_exons(pos, cigar)
            
            if exons:
                start = exons[0][0]
                end = exons[-1][1]
                isotag_details[isotag_id] = {
                    'chrom': chrom,
                    'strand': strand,
                    'start': start,
                    'end': end
                }
        
        process.wait()
    
    except Exception as e:
        click.echo(f"❌ Error reading BAM for TSV: {e}")
        return
    
    # Write TSV file
    with open(tsv_file, 'w') as f:
        # Header
        f.write("Isotag_ID\tCluster_ID\tRepresentative_ID\tCluster_Size\tChrom\tStrand\tStart\tEnd\n")
        
        # Sort by cluster ID for consistent output
        sorted_isotags = sorted(cluster_mapping.keys(), key=lambda x: cluster_mapping[x])
        
        for isotag_id in sorted_isotags:
            cluster_id = cluster_mapping[isotag_id]
            representative_id = representative_mapping[isotag_id]
            cluster_size = cluster_sizes[cluster_id]
            
            if isotag_id in isotag_details:
                details = isotag_details[isotag_id]
                f.write(f"{isotag_id}\t{cluster_id}\t{representative_id}\t{cluster_size}\t{details['chrom']}\t{details['strand']}\t{details['start']}\t{details['end']}\n")
            else:
                # Fallback if details not found
                f.write(f"{isotag_id}\t{cluster_id}\t{representative_id}\t{cluster_size}\t\t\t\t\n")
    
    total_entries = len(cluster_mapping)
    click.echo(f"📋 TSV written: {total_entries:,} entries → {tsv_file}")


def validate_inputs(input_bam: str, output_bam: str):
    """Validate input files and create output directory"""
    if not os.path.exists(input_bam):
        click.echo(f"❌ Input BAM not found: {input_bam}")
        sys.exit(1)
    
    # Check for BAM index
    index_file = f"{input_bam}.bai"
    if not os.path.exists(index_file):
        click.echo(f"⚠️  Creating BAM index: {index_file}")
        try:
            subprocess.run(f"samtools index {input_bam}".split(), check=True)
        except subprocess.CalledProcessError:
            click.echo(f"❌ Failed to create BAM index")
            sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_bam)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)


@click.command()
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('output_bam', type=click.Path())
@click.option('-t', '--stable-end-tolerance', type=int, default=20, help='BP tolerance for stable ends (default: 20)')
@click.option('-j', '--threads', type=int, default=4, help='Number of threads for chromosome processing (default: 4)')
@click.option('--json', help='Also output JSON statistics to this file')
@click.option('--tsv', help='Also output TSV file with isotag details to this file')
def main(input_bam, output_bam, stable_end_tolerance, threads, json, tsv):
    """
    IsoTag Clustering - Multi-threaded BAM-to-BAM Processing
    
    Clusters isotags by exon count with stable-end requirement.
    Adds XC cluster tags and XR representative tags to output BAM for downstream analysis.
    
    XC: Deterministic cluster ID (SHA-256 hash of sorted isotag IDs)
    XR: Longest isotag representative (biological meaning, deterministic tie-breaking)
    
    Examples:
        python3 isotag_clustering.py input.bam clustered.bam
        python3 isotag_clustering.py input.bam clustered.bam -j 8 -t 30
        python3 isotag_clustering.py input.bam clustered.bam --json stats.json
        python3 isotag_clustering.py input.bam clustered.bam --tsv results.tsv
        python3 isotag_clustering.py input.bam clustered.bam --json stats.json --tsv results.tsv
    """
    click.echo(f"🧬 IsoTag Clustering")
    click.echo(f"Input:  {input_bam}")
    click.echo(f"Output: {output_bam}")
    click.echo(f"Config: {threads} threads, {stable_end_tolerance}bp tolerance")
    click.echo(f"=" * 60)
    
    # Validate inputs
    validate_inputs(input_bam, output_bam)
    
    start_time = time.time()
    
    # Phase 1: Multi-threaded clustering by chromosome
    click.echo(f"🔄 Phase 1: Multi-threaded clustering...")
    chromosome_results = parallel_chromosome_clustering(input_bam, threads, stable_end_tolerance)
    clustering_time = time.time() - start_time
    
    # Combine results from all chromosomes
    all_isotag_to_cluster = {}
    all_isotag_to_representative = {}
    for chr_result in chromosome_results.values():
        all_isotag_to_cluster.update(chr_result['isotag_to_cluster'])
        all_isotag_to_representative.update(chr_result['isotag_to_representative'])
    
    click.echo(f"⚡ Clustering completed: {clustering_time:.1f}s")
    
    # Phase 2: Create tagged BAM with XC cluster IDs and XR representatives
    click.echo(f"🔄 Phase 2: Creating tagged BAM...")
    bam_start = time.time()
    create_tagged_bam(input_bam, output_bam, all_isotag_to_cluster, all_isotag_to_representative)
    bam_time = time.time() - bam_start
    
    total_time = time.time() - start_time
    
    # Phase 3: Optional JSON statistics
    if json:
        create_json_stats(chromosome_results, json, total_time)
    
    # Phase 4: Optional TSV output
    if tsv:
        create_tsv_output(input_bam, all_isotag_to_cluster, all_isotag_to_representative, tsv)
    
    # Summary
    total_isotags = len(all_isotag_to_cluster)
    total_clusters = len(set(all_isotag_to_cluster.values()))
    compression = (1 - total_clusters/total_isotags)*100 if total_isotags > 0 else 0
    
    click.echo(f"=" * 60)
    click.echo(f"✅ Complete: {total_time:.1f}s (clustering: {clustering_time:.1f}s, BAM: {bam_time:.1f}s)")
    click.echo(f"📊 {total_isotags:,} isotags → {total_clusters:,} clusters ({compression:.1f}% compression)")
    click.echo(f"🎯 Results: {output_bam} (XC cluster tags and XR representative tags added)")


if __name__ == '__main__':
    main()