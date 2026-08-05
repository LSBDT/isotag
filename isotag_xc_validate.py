#!/usr/bin/env python3
"""
IsoTag XC Validator - Empirical validation of XC clustering against GENCODE gene IDs

Validates that XC tags (pure location-based clustering) align with GENCODE gene assignments.

Metrics:
- Purity: % of XC clusters mapping to single gene
- Completeness: % of genes mapping to single XC cluster
- False merge rate: XC clusters incorrectly merging multiple genes
- Cluster size distribution

Input:
- BAM file with XC tags
- GENCODE GTF annotation

Output:
- TSV report with validation metrics
- Detailed problematic clusters

Memory: uses SQLite for all read/gene storage — constant RAM regardless of BAM size.
"""

import argparse
import hashlib
import math
import os
import random
import sys
import logging
import sqlite3
import statistics
import tempfile
import pysam
from pathlib import Path
from typing import Dict, List, Optional
import subprocess

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def create_validation_db(db_path: str) -> sqlite3.Connection:
    """
    Create a SQLite database for streaming validation data.

    Args:
        db_path: Path to SQLite file (use ':memory:' only for tiny datasets)

    Returns:
        Open SQLite connection
    """
    conn = sqlite3.connect(db_path)
    # Speed-up pragmas: we don't need durability for a temp analysis DB
    conn.execute('PRAGMA journal_mode = OFF')
    conn.execute('PRAGMA synchronous = OFF')
    conn.execute('PRAGMA cache_size = -32000')  # 32 MB page cache
    conn.execute('''
        CREATE TABLE read_xc (
            read_id TEXT PRIMARY KEY,
            xc      TEXT NOT NULL
        )
    ''')
    conn.execute('''
        CREATE TABLE read_gene (
            read_id TEXT NOT NULL,
            gene_id TEXT NOT NULL
        )
    ''')
    return conn


def stream_xc_to_db(bam_path: str, conn: sqlite3.Connection, bin_size: Optional[int] = None):
    """
    Stream XC tags from BAM into SQLite read_xc table.
    Skips secondary and supplementary alignments.

    Args:
        bam_path: Path to BAM with XC tags
        conn: Open SQLite connection
        bin_size: If set, ignore stored XC tags and recompute cluster IDs from read
                  coordinates using this bin size in bp (e.g. 5000, 10000, 50000).
                  XC format: chr:strand:b_start:b_end where b_start/b_end are coordinate // bin_size.
    """
    if bin_size is not None:
        logger.info(f"Extracting coordinates from {bam_path}, recomputing XC with bin_size={bin_size:,}")
    else:
        logger.info(f"Extracting XC tags from {bam_path}")

    total = 0
    skipped = 0
    n_xc = 0

    conn.execute('BEGIN')
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam:
            total += 1

            # Skip secondary/supplementary — multi-mappers inflate false merge rates
            if read.is_secondary or read.is_supplementary:
                skipped += 1
                continue

            if bin_size is not None:
                # Recompute XC from coordinates — ignore stored tag
                if read.is_unmapped:
                    continue
                chrom = read.reference_name
                strand = '-' if read.is_reverse else '+'
                b_start = read.reference_start // bin_size
                b_end = read.reference_end // bin_size
                xc = f"{chrom}:{strand}:{b_start}:{b_end}"
                conn.execute(
                    'INSERT OR REPLACE INTO read_xc (read_id, xc) VALUES (?, ?)',
                    (read.query_name, xc)
                )
                n_xc += 1
            elif read.has_tag('XC'):
                conn.execute(
                    'INSERT OR REPLACE INTO read_xc (read_id, xc) VALUES (?, ?)',
                    (read.query_name, read.get_tag('XC'))
                )
                n_xc += 1

            if total % 50000 == 0:
                conn.execute('COMMIT')
                conn.execute('BEGIN')
                logger.info(f"Processed {total:,} reads, {n_xc:,} with XC tags")

    conn.execute('COMMIT')

    primary = total - skipped
    if primary > 0:
        logger.info(
            f"Total: {total:,} reads ({skipped:,} secondary/supplementary skipped, "
            f"{primary:,} primary), {n_xc:,} with XC tags "
            f"({n_xc / primary * 100:.1f}% of primary)"
        )
    else:
        logger.info("Total: 0 primary reads")


def extract_genes_as_bed(gtf_path: str, bed_path: str) -> int:
    """
    Extract gene-level intervals from GTF as BED6 (chr, start, end, gene_id, 0, strand).

    Writing a gene-only BED6 (~60K lines) instead of passing the full GTF (~3M features)
    to bedtools reduces interval-tree RAM by ~50×, preventing OOM on shared servers.
    """
    import re
    count = 0
    with open(gtf_path) as gtf, open(bed_path, 'w') as bed:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            chrom = fields[0]
            start = int(fields[3]) - 1  # GTF is 1-based, BED is 0-based
            end = int(fields[4])
            strand = fields[6]
            m = re.search(r'gene_id "([^"]+)"', fields[8])
            if not m:
                continue
            gene_id = m.group(1).split('.')[0]  # strip version suffix
            bed.write(f'{chrom}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n')
            count += 1
    logger.info(f"Extracted {count:,} gene intervals from GTF as BED6")
    return count


def intersect_bam_with_gtf(bam_path: str, gtf_path: str, output_bed: str):
    """
    Intersect BAM with GTF using bedtools to assign genes to reads.

    Args:
        bam_path: Path to BAM file
        gtf_path: Path to GENCODE GTF
        output_bed: Path to output intersection file
    """
    logger.info("Intersecting BAM with GTF using bedtools")

    # Pre-extract gene-only intervals as BED6 to reduce bedtools RAM by ~50×.
    # Full GTF has ~3M features (genes + transcripts + exons + CDS); gene-only BED ~60K lines.
    # This prevents OOM on shared servers (prior failure: bedtools Killed at 20GB with full GTF).
    with tempfile.NamedTemporaryFile(mode='w', suffix='_genes.bed', delete=False) as tmp:
        gene_bed = tmp.name
    try:
        extract_genes_as_bed(gtf_path, gene_bed)

        # Pipe samtools view -F 0x900 | bedtools intersect to exclude secondary+supplementary
        # alignments before gene assignment (MAJOR-3 fix: bedtools -a BAM includes all alignment
        # types by default, inflating gene-assignment counts by ~71% from secondary reads)
        logger.info(
            f"Running: samtools view -F 0x900 -b {bam_path} | "
            f"bedtools intersect -a - -b {gene_bed} -wa -wb -bed -s"
        )

        with open(output_bed, 'w') as out:
            p_view = subprocess.Popen(['samtools', 'view', '-F', '0x900', '-b', bam_path],
                                       stdout=subprocess.PIPE)
            p_isect = subprocess.Popen(
                ['bedtools', 'intersect', '-a', '-', '-b', gene_bed, '-wa', '-wb', '-bed', '-s'],
                stdin=p_view.stdout, stdout=out, stderr=subprocess.PIPE, text=True,
            )
            p_view.stdout.close()
            _, isect_stderr = p_isect.communicate()
    finally:
        if os.path.exists(gene_bed):
            os.unlink(gene_bed)

    if p_isect.returncode != 0:
        logger.error(f"bedtools failed: {isect_stderr}")
        sys.exit(1)

    with open(output_bed) as f:
        n_overlaps = sum(1 for _ in f)

    logger.info(f"Found {n_overlaps:,} read-gene overlaps")


def stream_intersect_to_db(bed_path: str, conn: sqlite3.Connection):
    """
    Stream bedtools intersect output into SQLite read_gene table.
    Filters to gene-level GTF features only (avoids exon/transcript redundancy).

    Args:
        bed_path: Path to bedtools intersect output
        conn: Open SQLite connection
    """
    logger.info(f"Loading intersection results from {bed_path}")

    count = 0
    conn.execute('BEGIN')

    with open(bed_path) as f:
        for i, line in enumerate(f):
            fields = line.strip().split('\t')

            # bedtools -bed outputs BAM as BED12 (12 fields):
            #   0:chr 1:start 2:end 3:read_id 4:score 5:strand
            #   6:thickStart 7:thickEnd 8:itemRgb 9:blockCount 10:blockSizes 11:blockStarts
            # Then gene BED6 fields (6 fields):
            #   12:chr 13:start 14:end 15:gene_id 16:score 17:strand
            # (extract_genes_as_bed pre-filtered to gene-level, gene_id is already stripped)
            if len(fields) < 16:
                continue

            read_id = fields[3]
            gene_id = fields[15]  # gene_id column of BED6 (no GTF attribute parsing needed)

            if gene_id:
                conn.execute(
                    'INSERT INTO read_gene (read_id, gene_id) VALUES (?, ?)',
                    (read_id, gene_id)
                )
                count += 1

            if (i + 1) % 50000 == 0:
                conn.execute('COMMIT')
                conn.execute('BEGIN')
                logger.info(f"Parsed {i + 1:,} overlap lines, {count:,} gene assignments loaded")

    conn.execute('COMMIT')
    logger.info(f"Total: {count:,} read-gene pairs loaded into database")

    # Build indexes now that loading is complete (faster than indexing during insert)
    conn.execute('CREATE INDEX idx_rg_read ON read_gene(read_id)')
    conn.execute('CREATE INDEX idx_rx_xc ON read_xc(xc)')
    logger.info("Database indexes built")


def calculate_metrics_from_db(conn: sqlite3.Connection, xc_table: str = 'read_xc') -> Dict:
    """
    Calculate XC validation metrics from SQLite database.
    All computation done in SQL — no large structures loaded into RAM.

    Args:
        conn: Open SQLite connection with read_xc (or xc_table) and read_gene tables
        xc_table: Name of the table containing (read_id, xc) pairs.
                  Default 'read_xc'; use 'read_xc_shuffled' for negative-control permutations.

    Returns:
        Dict with validation metrics (summary numbers only, not per-cluster data)
    """
    logger.info(f"Calculating validation metrics from database (xc_table={xc_table})")

    # Basic read counts
    total_xc = conn.execute(f'SELECT COUNT(*) FROM {xc_table}').fetchone()[0]
    total_gene = conn.execute('SELECT COUNT(DISTINCT read_id) FROM read_gene').fetchone()[0]
    both = conn.execute(
        f'SELECT COUNT(*) FROM {xc_table} rx '
        'WHERE EXISTS (SELECT 1 FROM read_gene rg WHERE rg.read_id = rx.read_id)'
    ).fetchone()[0]

    if both == 0:
        logger.warning("No reads found with both XC tags and gene assignments")

    # XC cluster stats — only clusters that have at least one gene-assigned read
    conn.execute('DROP TABLE IF EXISTS xc_cluster_stats')
    conn.execute(f'''
        CREATE TEMP TABLE xc_cluster_stats AS
        SELECT rx.xc,
               COUNT(DISTINCT rx.read_id)  AS n_reads,
               COUNT(DISTINCT rg.gene_id)  AS n_genes
        FROM {xc_table} rx
        JOIN read_gene rg ON rx.read_id = rg.read_id
        GROUP BY rx.xc
    ''')
    logger.info(f"XC clusters with gene assignments: "
                f"{conn.execute('SELECT COUNT(*) FROM xc_cluster_stats').fetchone()[0]:,}")

    row = conn.execute('''
        SELECT
            COUNT(*)                                                     AS total,
            SUM(CASE WHEN n_reads = 1 THEN 1 ELSE 0 END)                AS singleton,
            SUM(CASE WHEN n_genes = 1 THEN 1 ELSE 0 END)                AS pure,
            SUM(CASE WHEN n_genes > 1 THEN 1 ELSE 0 END)                AS mixed,
            SUM(CASE WHEN n_genes = 1 AND n_reads > 1 THEN 1 ELSE 0 END) AS pure_multi,
            SUM(CASE WHEN n_genes > 1 AND n_reads > 1 THEN 1 ELSE 0 END) AS mixed_multi,
            MIN(n_reads),
            MAX(n_reads)
        FROM xc_cluster_stats
    ''').fetchone()
    xc_total, xc_singleton, xc_pure, xc_mixed, xc_pure_multi, xc_mixed_multi, size_min, size_max = row
    xc_total = xc_total or 0
    xc_singleton = xc_singleton or 0
    xc_multi = xc_total - xc_singleton

    # Cluster size median — fetch sorted sizes for statistics.median
    sizes = [r[0] for r in conn.execute('SELECT n_reads FROM xc_cluster_stats ORDER BY n_reads')]
    size_median = statistics.median(sizes) if sizes else 0

    # Cluster size distribution (top 10 most common sizes)
    size_dist = dict(conn.execute(
        'SELECT n_reads, COUNT(*) AS cnt FROM xc_cluster_stats '
        'GROUP BY n_reads ORDER BY cnt DESC LIMIT 10'
    ).fetchall())

    # Gene completeness — how many XC clusters per gene?
    conn.execute('DROP TABLE IF EXISTS gene_xc_stats')
    conn.execute(f'''
        CREATE TEMP TABLE gene_xc_stats AS
        SELECT rg.gene_id,
               COUNT(DISTINCT rx.xc) AS n_xc
        FROM read_gene rg
        JOIN {xc_table} rx ON rg.read_id = rx.read_id
        GROUP BY rg.gene_id
    ''')
    logger.info(f"Unique genes: "
                f"{conn.execute('SELECT COUNT(*) FROM gene_xc_stats').fetchone()[0]:,}")

    g_row = conn.execute('''
        SELECT COUNT(*),
               SUM(CASE WHEN n_xc = 1 THEN 1 ELSE 0 END),
               SUM(CASE WHEN n_xc > 1 THEN 1 ELSE 0 END)
        FROM gene_xc_stats
    ''').fetchone()
    gene_total, gene_pure, gene_split = g_row
    gene_total = gene_total or 0

    return {
        'total_reads_with_xc':     total_xc,
        'reads_with_gene_assignment': total_gene,
        'reads_with_both':         both,
        'xc_clusters_total':       xc_total,
        'xc_clusters_singleton':   xc_singleton,
        'xc_clusters_multi_read':  xc_multi,
        'xc_clusters_pure':        xc_pure or 0,
        'xc_clusters_mixed':       xc_mixed or 0,
        'xc_clusters_pure_multi':  xc_pure_multi or 0,
        'xc_clusters_mixed_multi': xc_mixed_multi or 0,
        'xc_purity_pct':           (xc_pure or 0) / xc_total * 100 if xc_total else 0,
        'xc_purity_multi_pct':     (xc_pure_multi or 0) / xc_multi * 100 if xc_multi else 0,
        'genes_total':             gene_total,
        'genes_pure':              gene_pure or 0,
        'genes_split':             gene_split or 0,
        'gene_completeness_pct':   (gene_pure or 0) / gene_total * 100 if gene_total else 0,
        'false_merge_rate':        (xc_mixed or 0) / xc_total * 100 if xc_total else 0,
        'cluster_size_min':        size_min or 0,
        'cluster_size_max':        size_max or 0,
        'cluster_size_median':     size_median,
        'cluster_size_distribution': size_dist,
    }


def calculate_ari_from_db(conn: sqlite3.Connection) -> float:
    """
    Compute Adjusted Rand Index (ARI) between XC clusters and GENCODE gene labels.

    Uses the contingency table formula — entirely in SQL, no large Python structures.
    ARI = 1.0 means XC clusters perfectly match gene assignments.
    ARI = 0.0 means XC clustering is no better than random.
    ARI < 0 means XC clustering is WORSE than random.

    This is a null-free metric: ARI mathematically adjusts for chance without
    needing permutation tests or spatial null models.

    Formula (from Hubert & Arabie 1985):
        ARI = (sum_ij C(n_ij,2) - expected) / (max_possible - expected)
        where expected = [sum_i C(a_i,2) * sum_j C(b_j,2)] / C(n,2)
    """
    logger.info("Computing Adjusted Rand Index (ARI) from contingency table")

    # Contingency table: n_ij = reads in both XC cluster i AND gene j
    conn.execute('DROP TABLE IF EXISTS ari_contingency')
    conn.execute('''
        CREATE TEMP TABLE ari_contingency AS
        SELECT rx.xc, rg.gene_id, COUNT(*) AS n_ij
        FROM read_xc rx
        JOIN read_gene rg ON rx.read_id = rg.read_id
        GROUP BY rx.xc, rg.gene_id
    ''')

    # sum_ij C(n_ij, 2) — use REAL to avoid integer overflow for large counts
    sum_comb_nij = conn.execute(
        'SELECT COALESCE(SUM(CAST(n_ij AS REAL) * (n_ij - 1) / 2), 0) FROM ari_contingency'
    ).fetchone()[0]

    # Row sums a_i (XC cluster sizes): sum_i C(a_i, 2)
    sum_comb_ai = conn.execute('''
        SELECT COALESCE(SUM(CAST(a AS REAL) * (a - 1) / 2), 0)
        FROM (SELECT SUM(n_ij) AS a FROM ari_contingency GROUP BY xc)
    ''').fetchone()[0]

    # Column sums b_j (gene sizes): sum_j C(b_j, 2)
    sum_comb_bj = conn.execute('''
        SELECT COALESCE(SUM(CAST(b AS REAL) * (b - 1) / 2), 0)
        FROM (SELECT SUM(n_ij) AS b FROM ari_contingency GROUP BY gene_id)
    ''').fetchone()[0]

    # n = total reads with both XC and gene assignment
    n = conn.execute(
        'SELECT COALESCE(SUM(n_ij), 0) FROM ari_contingency'
    ).fetchone()[0]

    conn.execute('DROP TABLE IF EXISTS ari_contingency')

    if n < 2:
        logger.warning("ARI: fewer than 2 reads with both XC and gene — returning 0.0")
        return 0.0

    cn2 = n * (n - 1) / 2
    expected = sum_comb_ai * sum_comb_bj / cn2
    numerator = sum_comb_nij - expected
    denominator = 0.5 * (sum_comb_ai + sum_comb_bj) - expected

    if denominator == 0:
        ari = 1.0 if numerator == 0 else 0.0
    else:
        ari = numerator / denominator

    logger.info(f"ARI = {ari:.4f} (n={n:,} reads, sum_comb_nij={sum_comb_nij:.0f})")
    return ari


def calculate_vmeasure_from_db(conn: sqlite3.Connection) -> Dict[str, float]:
    """
    Compute V-measure (homogeneity + completeness) between XC clusters and gene labels.

    V-measure = harmonic mean of:
    - Homogeneity: do all reads in an XC cluster map to a single gene?
    - Completeness: do all reads from a gene map to a single XC cluster?

    Note: high homogeneity + low completeness is expected for XC (it is a locus ID,
    not a gene ID — long genes span multiple XC bins by design).

    Reference: Rosenberg & Hirschberg (2007); adopted as standard in isONclust3
    (Bioinformatics 2025).
    """
    logger.info("Computing V-measure (homogeneity + completeness) from contingency table")

    rows = conn.execute('''
        SELECT rx.xc, rg.gene_id, COUNT(*) AS n_ij
        FROM read_xc rx
        JOIN read_gene rg ON rx.read_id = rg.read_id
        GROUP BY rx.xc, rg.gene_id
    ''').fetchall()

    if not rows:
        return {'homogeneity': 0.0, 'completeness': 0.0, 'v_measure': 0.0}

    cluster_totals: Dict[str, int] = {}
    gene_totals: Dict[str, int] = {}
    n = 0
    for xc, gene_id, n_ij in rows:
        cluster_totals[xc] = cluster_totals.get(xc, 0) + n_ij
        gene_totals[gene_id] = gene_totals.get(gene_id, 0) + n_ij
        n += n_ij

    if n == 0:
        return {'homogeneity': 0.0, 'completeness': 0.0, 'v_measure': 0.0}

    h_c = -sum((cnt / n) * math.log(cnt / n) for cnt in gene_totals.values())
    h_k = -sum((cnt / n) * math.log(cnt / n) for cnt in cluster_totals.values())
    h_c_given_k = -sum(
        (n_ij / n) * math.log(n_ij / cluster_totals[xc])
        for xc, gene_id, n_ij in rows
    )
    h_k_given_c = -sum(
        (n_ij / n) * math.log(n_ij / gene_totals[gene_id])
        for xc, gene_id, n_ij in rows
    )

    homogeneity = 1.0 - h_c_given_k / h_c if h_c > 0 else 1.0
    completeness = 1.0 - h_k_given_c / h_k if h_k > 0 else 1.0
    v_measure = (
        2 * homogeneity * completeness / (homogeneity + completeness)
        if (homogeneity + completeness) > 0 else 0.0
    )

    logger.info(f"V-measure = {v_measure:.4f} (homogeneity={homogeneity:.4f}, completeness={completeness:.4f})")
    return {'homogeneity': homogeneity, 'completeness': completeness, 'v_measure': v_measure}


def run_negative_control(conn: sqlite3.Connection, n_permutations: int) -> Dict[str, List[float]]:
    """
    Shuffle XC tag assignments randomly across reads N times and compute validation
    metrics for each permutation.  Reports observed vs. random to quantify how much
    better than chance the real XC clustering is.

    Memory: fetches all (read_id, xc) pairs into Python lists once (~400 MB for 8M
    reads), then processes one permutation at a time.  Does NOT accumulate all
    permutation results in memory.

    Args:
        conn: Open SQLite connection (read_xc and read_gene tables must exist)
        n_permutations: Number of shuffle permutations to run

    Returns:
        Dict mapping metric name → list of values (one per permutation):
            'xc_purity_multi_pct', 'gene_completeness_pct', 'false_merge_rate'
    """
    logger.info(f"Running {n_permutations} negative-control permutations")

    # Fetch all (read_id, xc) pairs — one-time cost
    rows = conn.execute('SELECT read_id, xc FROM read_xc').fetchall()
    if not rows:
        logger.warning("No rows in read_xc — skipping negative control")
        return {}

    read_ids = [r[0] for r in rows]
    xc_values = [r[1] for r in rows]
    logger.info(f"Loaded {len(read_ids):,} reads for shuffling")

    # Temp table for shuffled assignments — reused across permutations
    conn.execute('DROP TABLE IF EXISTS read_xc_shuffled')
    conn.execute('''
        CREATE TEMP TABLE read_xc_shuffled (
            read_id TEXT PRIMARY KEY,
            xc      TEXT NOT NULL
        )
    ''')

    results: Dict[str, List[float]] = {
        'xc_purity_multi_pct': [],
        'gene_completeness_pct': [],
        'false_merge_rate': [],
    }

    for perm_idx in range(n_permutations):
        logger.info(f"Negative control permutation {perm_idx + 1}/{n_permutations}")

        # Shuffle XC values in-place (fresh copy each time)
        shuffled = xc_values[:]
        random.shuffle(shuffled)

        # Populate temp table — commit the implicit DELETE transaction before explicit BEGIN
        conn.execute('DELETE FROM read_xc_shuffled')
        conn.commit()
        conn.execute('BEGIN')
        batch = []
        for i, (rid, xc) in enumerate(zip(read_ids, shuffled)):
            batch.append((rid, xc))
            if len(batch) == 50000:
                conn.executemany('INSERT INTO read_xc_shuffled VALUES (?, ?)', batch)
                batch = []
                conn.execute('COMMIT')
                conn.execute('BEGIN')
        if batch:
            conn.executemany('INSERT INTO read_xc_shuffled VALUES (?, ?)', batch)
        conn.execute('COMMIT')

        # Calculate metrics using the shuffled table
        m = calculate_metrics_from_db(conn, xc_table='read_xc_shuffled')
        results['xc_purity_multi_pct'].append(m['xc_purity_multi_pct'])
        results['gene_completeness_pct'].append(m['gene_completeness_pct'])
        results['false_merge_rate'].append(m['false_merge_rate'])

    conn.execute('DROP TABLE IF EXISTS read_xc_shuffled')
    logger.info("Negative control permutations complete")
    return results


def run_negative_control_chrom_local(
    conn: sqlite3.Connection,
    n_permutations: int,
) -> Dict[str, List[float]]:
    """
    Chromosome-local negative control: shuffle XC labels only within each chromosome.

    This destroys locus-level signal while preserving chromosomal distribution.
    Better than global shuffle (which moves reads to wrong chromosomes, trivially
    collapsing purity to ~0%) and better than spatial offset (which preserves
    neighbourhood gene density, confounding the null with local genomic structure).

    XC tags encode the chromosome hash as the first 32 chars of the hash, so we
    cannot recover chromosome from XC alone.  Instead we read chromosome directly
    from the XC tag prefix stored in read_xc: XC values computed by isotag.py
    are 32-char MD5 hashes.  We therefore JOIN with the BAM-derived chromosome
    stored in a temp table populated from read_xc via the read's chromosome field.

    Because we cannot recover chromosome from a hash, we require the caller to
    have previously stored (read_id → chromosome) in a 'read_chrom' table.
    If that table is absent this function falls back to global shuffle with a warning.
    """
    logger.info(f"Running {n_permutations} chromosome-local negative-control permutations")

    # Check if read_chrom table exists (populated by cmd_validate when chrom null requested)
    tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_temp_master WHERE type='table'")}
    if 'read_chrom' not in tables:
        logger.warning("read_chrom table missing — falling back to global shuffle for null model")
        return run_negative_control(conn, n_permutations)

    # Fetch (read_id, xc, chrom) — grouped by chromosome for within-chr shuffling
    rows_by_chrom: Dict[str, list] = {}
    for read_id, xc, chrom in conn.execute(
        'SELECT rx.read_id, rx.xc, rc.chrom '
        'FROM read_xc rx JOIN read_chrom rc ON rx.read_id = rc.read_id'
    ):
        rows_by_chrom.setdefault(chrom, []).append((read_id, xc))

    n_total = sum(len(v) for v in rows_by_chrom.values())
    logger.info(f"Loaded {n_total:,} reads across {len(rows_by_chrom)} chromosomes for chrom-local shuffle")

    conn.execute('DROP TABLE IF EXISTS read_xc_shuffled')
    conn.execute('''
        CREATE TEMP TABLE read_xc_shuffled (
            read_id TEXT PRIMARY KEY,
            xc      TEXT NOT NULL
        )
    ''')

    results: Dict[str, List[float]] = {
        'xc_purity_multi_pct': [],
        'gene_completeness_pct': [],
        'false_merge_rate': [],
        'type': 'chrom_local',
    }

    for perm_idx in range(n_permutations):
        logger.info(f"Chrom-local neg-control permutation {perm_idx + 1}/{n_permutations}")

        conn.execute('DELETE FROM read_xc_shuffled')
        conn.commit()
        conn.execute('BEGIN')
        batch = []

        for chrom_rows in rows_by_chrom.values():
            read_ids_chr = [r[0] for r in chrom_rows]
            xc_vals_chr  = [r[1] for r in chrom_rows]
            random.shuffle(xc_vals_chr)
            for rid, xc in zip(read_ids_chr, xc_vals_chr):
                batch.append((rid, xc))
                if len(batch) == 50000:
                    conn.executemany('INSERT INTO read_xc_shuffled VALUES (?, ?)', batch)
                    batch = []
                    conn.execute('COMMIT')
                    conn.execute('BEGIN')

        if batch:
            conn.executemany('INSERT INTO read_xc_shuffled VALUES (?, ?)', batch)
        conn.execute('COMMIT')

        m = calculate_metrics_from_db(conn, xc_table='read_xc_shuffled')
        results['xc_purity_multi_pct'].append(m['xc_purity_multi_pct'])
        results['gene_completeness_pct'].append(m['gene_completeness_pct'])
        results['false_merge_rate'].append(m['false_merge_rate'])

    conn.execute('DROP TABLE IF EXISTS read_xc_shuffled')
    logger.info("Chromosome-local negative control permutations complete")
    return results


def run_negative_control_spatial(
    bam_path: str,
    conn: sqlite3.Connection,
    n_permutations: int,
    offset_bins: int = 50,
    xc_position_quantum: int = 1000,
    xc_bin_size: int = 10,
) -> Dict[str, List[float]]:
    """
    Spatially-aware negative control: for each read, shift its midpoint bin by an
    independent random offset in [-offset_bins, +offset_bins] while preserving
    chromosome, strand, and exon length structure.

    This is a more realistic null than global label shuffling because XC is
    location-based; globally shuffling XC tags trivially yields near-0% purity
    (reads from chr1 land in chrX clusters).  Spatial shifts test whether XC
    outperforms random assignment *within the same genomic neighbourhood*.

    Args:
        bam_path: BAM file with XC tags and CIGAR strings
        conn: Open SQLite connection (read_xc and read_gene tables must exist)
        n_permutations: Number of permutations
        offset_bins: Max random offset in units of xc_position_quantum (default 50 → ±50 kb)
        xc_position_quantum: Midpoint bin size in bp (default 1000 = 1 kb)
        xc_bin_size: Exon-length bin size in bp (default 10)

    Returns:
        Dict mapping metric name → list of values (one per permutation), plus
        metadata keys 'type', 'offset_bins', 'xc_position_quantum'.
    """
    kb_range = offset_bins * xc_position_quantum // 1000
    logger.info(
        f"Running {n_permutations} spatially-aware negative-control permutations "
        f"(independent per-read offset ±{offset_bins} bins = ±{kb_range} kb; "
        f"preserves chr, strand, exon lengths)"
    )

    # ── Step 1: read BAM once, cache per-read spatial data ──────────────────────
    logger.info("Reading BAM coordinates + CIGAR for spatial null model")
    read_spatial = []  # list of (read_id, chrom, strand, midpoint_bin, exon_bins_str)
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            if not read.has_tag('XC') or not read.cigartuples:
                continue

            chrom = read.reference_name
            strand = '-' if read.is_reverse else '+'
            midpoint = (read.reference_start + read.reference_end) // 2
            midpoint_bin = midpoint // xc_position_quantum

            # Exon length bins from CIGAR — same logic as isotag.py extract_exons_from_cigar
            # M/=/X/D advance reference (exon); N ends exon (intron); I/S/H do not advance
            exon_bins = []
            exon_len = 0
            for op, length in read.cigartuples:
                if op in (0, 2, 7, 8):   # M, D, =, X — within exon, consumes reference
                    exon_len += length
                elif op == 3:             # N — intron: save exon and reset
                    if exon_len > 0:
                        exon_bins.append(exon_len // xc_bin_size)
                        exon_len = 0
                # op 1=I, 4=S, 5=H — do not consume reference; skip
            if exon_len > 0:
                exon_bins.append(exon_len // xc_bin_size)

            exon_bins_str = '|'.join(map(str, exon_bins))
            read_spatial.append((read.query_name, chrom, strand, midpoint_bin, exon_bins_str))

    logger.info(f"Loaded spatial data for {len(read_spatial):,} reads")
    if not read_spatial:
        logger.warning("No spatial data loaded — skipping spatial negative control")
        return {}

    # ── Step 2: temp table reused across permutations ───────────────────────────
    conn.execute('DROP TABLE IF EXISTS read_xc_shuffled')
    conn.execute('''
        CREATE TEMP TABLE read_xc_shuffled (
            read_id TEXT PRIMARY KEY,
            xc      TEXT NOT NULL
        )
    ''')

    results: Dict[str, List[float]] = {
        'xc_purity_multi_pct': [],
        'gene_completeness_pct': [],
        'false_merge_rate': [],
        'type': 'spatial',                        # metadata (not a list)
        'offset_bins': offset_bins,               # metadata
        'xc_position_quantum': xc_position_quantum,  # metadata
    }

    # ── Step 3: permutations ────────────────────────────────────────────────────
    for perm_idx in range(n_permutations):
        logger.info(f"Spatial neg-control permutation {perm_idx + 1}/{n_permutations}")

        conn.execute('DELETE FROM read_xc_shuffled')
        conn.commit()
        conn.execute('BEGIN')
        batch = []

        for read_id, chrom, strand, midpoint_bin, exon_bins_str in read_spatial:
            offset = random.randint(-offset_bins, offset_bins)
            shifted_bin = midpoint_bin + offset
            null_serial = f"{chrom}|{strand}|{shifted_bin}|{exon_bins_str}"
            null_xc = hashlib.sha256(null_serial.encode()).hexdigest()[:32]
            batch.append((read_id, null_xc))

            if len(batch) == 50000:
                conn.executemany('INSERT INTO read_xc_shuffled VALUES (?, ?)', batch)
                batch = []
                conn.execute('COMMIT')
                conn.execute('BEGIN')

        if batch:
            conn.executemany('INSERT INTO read_xc_shuffled VALUES (?, ?)', batch)
        conn.execute('COMMIT')

        m = calculate_metrics_from_db(conn, xc_table='read_xc_shuffled')
        results['xc_purity_multi_pct'].append(m['xc_purity_multi_pct'])
        results['gene_completeness_pct'].append(m['gene_completeness_pct'])
        results['false_merge_rate'].append(m['false_merge_rate'])

    conn.execute('DROP TABLE IF EXISTS read_xc_shuffled')
    logger.info("Spatial negative control permutations complete")
    return results


def write_report(
    metrics: Dict,
    conn: sqlite3.Connection,
    output_path: str,
    neg_control: Optional[Dict[str, List[float]]] = None,
):
    """
    Write validation report to TSV file.
    Queries SQLite for per-cluster and per-gene detail sections.

    Args:
        metrics: Summary validation metrics
        conn: Open SQLite connection (needs xc_cluster_stats, gene_xc_stats temp tables)
        output_path: Path to output TSV
        neg_control: Optional dict from run_negative_control(); if provided, appends
                     a negative-control section with observed vs. random comparison.
    """
    logger.info(f"Writing report to {output_path}")

    with open(output_path, 'w') as out:
        out.write("# IsoTag XC Validation Report\n")
        out.write("# ========================================\n\n")

        out.write("## Summary Metrics\n")
        out.write("Metric\tValue\n")
        out.write(f"Total reads with XC tags\t{metrics['total_reads_with_xc']}\n")
        out.write(f"Reads with gene assignment\t{metrics['reads_with_gene_assignment']}\n")
        out.write(f"Reads with both XC and gene\t{metrics['reads_with_both']}\n")
        out.write(f"XC clusters (total)\t{metrics['xc_clusters_total']}\n")
        out.write(f"XC clusters (singleton, 1 read)\t{metrics['xc_clusters_singleton']}\n")
        out.write(f"XC clusters (multi-read, >1 read)\t{metrics['xc_clusters_multi_read']}\n")
        out.write(f"XC clusters (pure, 1 gene)\t{metrics['xc_clusters_pure']}\n")
        out.write(f"XC clusters (mixed, >1 gene)\t{metrics['xc_clusters_mixed']}\n")
        out.write(f"XC purity - all clusters (%)\t{metrics['xc_purity_pct']:.1f}\n")
        out.write(f"XC clusters (pure multi-read)\t{metrics['xc_clusters_pure_multi']}\n")
        out.write(f"XC clusters (mixed multi-read)\t{metrics['xc_clusters_mixed_multi']}\n")
        out.write(f"XC purity - multi-read only (%)\t{metrics['xc_purity_multi_pct']:.1f}\n")
        out.write("NOTE: Singleton clusters are trivially pure (1 read cannot be mixed)\n")
        out.write(f"Genes (total)\t{metrics['genes_total']}\n")
        out.write(f"Genes (pure, 1 XC cluster)\t{metrics['genes_pure']}\n")
        out.write(f"Genes (split, >1 XC cluster)\t{metrics['genes_split']}\n")
        out.write(f"Gene completeness (%)\t{metrics['gene_completeness_pct']:.1f}\n")
        out.write(f"False merge rate (%)\t{metrics['false_merge_rate']:.1f}\n")
        if 'adjusted_rand_index' in metrics:
            out.write(f"Adjusted Rand Index (ARI)\t{metrics['adjusted_rand_index']:.4f}\n")
            out.write("NOTE: ARI=1.0=perfect clustering, ARI~0.0=random, ARI<0=worse than random\n")
        if 'v_measure' in metrics:
            out.write(f"V-measure\t{metrics['v_measure']:.4f}\n")
            out.write(f"V-measure homogeneity\t{metrics['homogeneity']:.4f}\n")
            out.write(f"V-measure completeness\t{metrics['completeness_entropy']:.4f}\n")
            out.write("NOTE: Low V-measure completeness expected — XC is a locus ID not gene ID\n")
        out.write(f"Cluster size (min)\t{metrics['cluster_size_min']}\n")
        out.write(f"Cluster size (max)\t{metrics['cluster_size_max']}\n")
        out.write(f"Cluster size (median)\t{metrics['cluster_size_median']}\n")
        out.write("\n")

        out.write("## Cluster Size Distribution\n")
        out.write("Size\tCount\n")
        for size, count in sorted(metrics['cluster_size_distribution'].items()):
            out.write(f"{size}\t{count}\n")
        out.write("\n")

        # Negative control section (optional)
        if neg_control:
            purity_list = neg_control['xc_purity_multi_pct']
            n_perms = len(purity_list)
            null_type = neg_control.get('type', 'global')
            if null_type == 'spatial':
                offset_bins = neg_control.get('offset_bins', 50)
                xc_pq = neg_control.get('xc_position_quantum', 1000)
                kb = offset_bins * xc_pq // 1000
                header = (f"## Negative Control ({n_perms} spatial permutations, "
                          f"±{offset_bins} bins = ±{kb} kb per-read midpoint shift)\n")
                out.write(header)
                out.write("# Each read's midpoint is shifted by an independent random offset;\n")
                out.write("# chromosome, strand, and exon lengths are preserved.\n")
            elif null_type == 'chrom_local':
                out.write(f"## Negative Control ({n_perms} chromosome-local shuffle permutations)\n")
                out.write("# XC labels shuffled within each chromosome — destroys locus signal,\n")
                out.write("# preserves chromosomal read distribution. Better null than spatial offset.\n")
            else:
                out.write(f"## Negative Control ({n_perms} global-shuffle permutations)\n")
                out.write("# XC labels shuffled globally across all reads.\n")
            out.write("# Z-score = (observed - random_mean) / random_std\n")
            out.write("# High Z-score means XC clustering is much better than random\n")
            out.write("Metric\tObserved\tRandom_Mean\tRandom_Std\tZ_Score\n")

            metric_map = [
                ('xc_purity_multi_pct',  'XC_Purity_Multi_Pct'),
                ('gene_completeness_pct', 'Gene_Completeness_Pct'),
                ('false_merge_rate',      'False_Merge_Rate'),
            ]
            # Observed values are in metrics dict
            observed_map = {
                'xc_purity_multi_pct':  metrics['xc_purity_multi_pct'],
                'gene_completeness_pct': metrics['gene_completeness_pct'],
                'false_merge_rate':      metrics['false_merge_rate'],
            }

            for key, label in metric_map:
                vals = neg_control[key]
                observed = observed_map[key]
                rand_mean = statistics.mean(vals) if vals else 0.0
                rand_std = statistics.stdev(vals) if len(vals) > 1 else 0.0
                z_score = (observed - rand_mean) / rand_std if rand_std > 0 else float('inf')
                out.write(f"{label}\t{observed:.2f}\t{rand_mean:.2f}\t{rand_std:.2f}\t{z_score:.1f}\n")
            out.write("\n")

        # Problematic XC clusters — query DB directly (avoids loading full dict)
        out.write("## Problematic XC Clusters (Multiple Genes)\n")
        out.write("XC_Tag\tNum_Genes\tGene_IDs\n")
        for (xc, n_genes) in conn.execute(
            'SELECT xc, n_genes FROM xc_cluster_stats WHERE n_genes > 1 ORDER BY n_genes DESC'
        ):
            genes = [r[0] for r in conn.execute(
                'SELECT DISTINCT rg.gene_id FROM read_gene rg '
                'JOIN read_xc rx ON rg.read_id = rx.read_id '
                'WHERE rx.xc = ? ORDER BY rg.gene_id',
                (xc,)
            )]
            out.write(f"{xc}\t{n_genes}\t{','.join(genes)}\n")
        out.write("\n")

        # Split genes — query DB directly
        out.write("## Split Genes (Multiple XC Clusters)\n")
        out.write("Gene_ID\tNum_XC_Clusters\tXC_Tags\n")
        for (gene_id, n_xc) in conn.execute(
            'SELECT gene_id, n_xc FROM gene_xc_stats WHERE n_xc > 1 ORDER BY n_xc DESC'
        ):
            xcs = [r[0] for r in conn.execute(
                'SELECT DISTINCT rx.xc FROM read_xc rx '
                'JOIN read_gene rg ON rx.read_id = rg.read_id '
                'WHERE rg.gene_id = ? ORDER BY rx.xc',
                (gene_id,)
            )]
            out.write(f"{gene_id}\t{n_xc}\t{','.join(xcs)}\n")

    logger.info("Report written successfully")


def cmd_validate(args):
    """Run XC validation"""

    output_dir = Path(args.output).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    intersect_bed = output_dir / "intersect_temp.bed"

    # SQLite temp file — use output_dir (not /tmp which has limited space)
    db_fd, db_path = tempfile.mkstemp(suffix='.db', prefix='isotag_xc_', dir=str(output_dir))
    os.close(db_fd)

    conn = None
    try:
        conn = create_validation_db(db_path)

        # Step 1: Stream XC tags from BAM into SQLite
        stream_xc_to_db(args.bam, conn, bin_size=args.recompute_xc)
        n_xc = conn.execute('SELECT COUNT(*) FROM read_xc').fetchone()[0]
        if n_xc == 0:
            logger.error("No XC tags found in BAM file")
            sys.exit(1)

        # Step 2: Intersect BAM with GTF
        intersect_bam_with_gtf(args.bam, args.gtf, str(intersect_bed))

        # Step 3: Stream intersection results into SQLite
        stream_intersect_to_db(str(intersect_bed), conn)
        n_pairs = conn.execute('SELECT COUNT(*) FROM read_gene').fetchone()[0]
        if n_pairs == 0:
            logger.error("No read-gene overlaps found")
            sys.exit(1)

        # Step 4: Calculate metrics from DB
        metrics = calculate_metrics_from_db(conn)

        # Step 4b: Compute ARI and V-measure (null-free clustering quality metrics)
        metrics['adjusted_rand_index'] = calculate_ari_from_db(conn)
        vmeasure = calculate_vmeasure_from_db(conn)
        metrics['homogeneity'] = vmeasure['homogeneity']
        metrics['completeness_entropy'] = vmeasure['completeness']
        metrics['v_measure'] = vmeasure['v_measure']

        # Step 5 (optional): Run negative control permutations
        neg_control = None
        if args.negative_control > 0:
            if args.negative_control_type == 'spatial':
                neg_control = run_negative_control_spatial(
                    args.bam, conn, args.negative_control,
                    offset_bins=args.neg_control_offset_bins,
                    xc_position_quantum=args.xc_position_quantum,
                    xc_bin_size=args.xc_bin_size,
                )
            elif args.negative_control_type == 'chrom':
                # Populate read_chrom temp table from BAM (needed by chrom-local null)
                logger.info("Building read_chrom table for chromosome-local null model")
                conn.execute('DROP TABLE IF EXISTS read_chrom')
                conn.execute('''
                    CREATE TEMP TABLE read_chrom (
                        read_id TEXT PRIMARY KEY,
                        chrom   TEXT NOT NULL
                    )
                ''')
                conn.execute('BEGIN')
                batch = []
                with pysam.AlignmentFile(args.bam, 'rb') as bam:
                    for read in bam:
                        if read.is_secondary or read.is_supplementary or read.is_unmapped:
                            continue
                        batch.append((read.query_name, read.reference_name))
                        if len(batch) == 50000:
                            conn.executemany('INSERT OR IGNORE INTO read_chrom VALUES (?, ?)', batch)
                            batch = []
                            conn.execute('COMMIT')
                            conn.execute('BEGIN')
                if batch:
                    conn.executemany('INSERT OR IGNORE INTO read_chrom VALUES (?, ?)', batch)
                conn.execute('COMMIT')
                n_chrom = conn.execute('SELECT COUNT(*) FROM read_chrom').fetchone()[0]
                logger.info(f"read_chrom populated: {n_chrom:,} reads")
                neg_control = run_negative_control_chrom_local(conn, args.negative_control)
            else:
                neg_control = run_negative_control(conn, args.negative_control)

        # Step 6: Write report (queries DB for detail sections)
        write_report(metrics, conn, args.output, neg_control=neg_control)

        # Print summary to console
        n_multi = metrics['xc_clusters_multi_read']
        print("\n" + "="*50)
        print("XC VALIDATION RESULTS")
        print("="*50)
        if args.recompute_xc:
            print(f"Mode: recomputed XC (bin_size={args.recompute_xc:,} bp)")
        print(f"XC Purity (all clusters):       {metrics['xc_purity_pct']:.1f}%"
              f" ({metrics['xc_clusters_pure']}/{metrics['xc_clusters_total']} pure)")
        print(f"  Singleton clusters:           {metrics['xc_clusters_singleton']}"
              " (trivially pure — 1 read cannot be mixed)")
        print(f"XC Purity (multi-read only):    {metrics['xc_purity_multi_pct']:.1f}%"
              f" ({metrics['xc_clusters_pure_multi']}/{n_multi} non-trivial clusters pure)")
        print(f"Gene Completeness:              {metrics['gene_completeness_pct']:.1f}%"
              f" ({metrics['genes_pure']}/{metrics['genes_total']} genes in single XC cluster)")
        print(f"False Merge Rate:               {metrics['false_merge_rate']:.1f}%"
              f" ({metrics['xc_clusters_mixed']} clusters merge multiple genes)")
        print(f"Adjusted Rand Index (ARI):      {metrics['adjusted_rand_index']:.4f}"
              " (1.0=perfect, 0.0=random, <0=worse than random)")
        print(f"V-measure:                      {metrics['v_measure']:.4f}"
              f" (homogeneity={metrics['homogeneity']:.4f},"
              f" completeness={metrics['completeness_entropy']:.4f})")
        print("  Note: low completeness expected — XC is a locus ID, not gene ID")

        if neg_control:
            purity_list = neg_control['xc_purity_multi_pct']
            n_perms = len(purity_list)
            null_type = neg_control.get('type', 'global')
            if null_type == 'spatial':
                ob = neg_control.get('offset_bins', 50)
                pq = neg_control.get('xc_position_quantum', 1000)
                print(f"\nNegative Control ({n_perms} spatial permutations, ±{ob} bins = ±{ob * pq // 1000} kb):")
            else:
                print(f"\nNegative Control ({n_perms} global-shuffle permutations):")
            obs_p = metrics['xc_purity_multi_pct']
            rand_p = statistics.mean(purity_list)
            print(f"  XC Purity (multi-read):  Observed {obs_p:.1f}% vs. Random {rand_p:.1f}%")
            obs_c = metrics['gene_completeness_pct']
            rand_c = statistics.mean(neg_control['gene_completeness_pct'])
            print(f"  Gene Completeness:       Observed {obs_c:.1f}% vs. Random {rand_c:.1f}%")

        print(f"\nDetailed report: {args.output}")
        print("="*50 + "\n")

    finally:
        if conn:
            conn.close()
        if os.path.exists(db_path):
            os.unlink(db_path)
        if intersect_bed.exists():
            intersect_bed.unlink()
            logger.info("Cleaned up temporary files")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Validate XC clustering against GENCODE gene IDs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate XC tags against GENCODE v47
  python3 isotag_xc_validate.py validate \\
      --bam test_xc_default.bam \\
      --gtf gencode.v47.annotation.gtf \\
      --output xc_validation_report.tsv

  # Negative control: 10 spatially-aware permutations (default, ±50 kb per-read shift)
  python3 isotag_xc_validate.py validate \\
      --bam test_xc_default.bam \\
      --gtf gencode.v47.annotation.gtf \\
      --negative-control 10 \\
      --output xc_validation_with_control.tsv

  # Negative control: global label shuffle (old behaviour, for comparison)
  python3 isotag_xc_validate.py validate \\
      --bam test_xc_default.bam \\
      --gtf gencode.v47.annotation.gtf \\
      --negative-control 10 --negative-control-type global \\
      --output xc_validation_global_control.tsv

  # Bin-size sensitivity: recompute XC at 20kb instead of stored 10kb
  python3 isotag_xc_validate.py validate \\
      --bam test_xc_default.bam \\
      --gtf gencode.v47.annotation.gtf \\
      --recompute-xc 20000 \\
      --output xc_validation_20kb.tsv

  # Bin-size sweep (run 5 times with different bin sizes)
  for BIN_SIZE in 5000 10000 20000 50000 100000; do
    python3 isotag_xc_validate.py validate \\
        --bam test_xc_default.bam \\
        --gtf gencode.v47.annotation.gtf \\
        --recompute-xc $BIN_SIZE \\
        -o results/xc_validation_binsize_${BIN_SIZE}.tsv
  done
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Subcommands')

    validate_parser = subparsers.add_parser('validate', help='Run XC validation')
    validate_parser.add_argument('--bam', required=True, help='BAM file with XC tags')
    validate_parser.add_argument('--gtf', required=True, help='GENCODE GTF annotation')
    validate_parser.add_argument('-o', '--output', required=True, help='Output report (TSV)')
    validate_parser.add_argument(
        '--negative-control', type=int, default=0, metavar='N',
        help='Run N random-shuffle permutations as negative control (default: 0 = skip). '
             'Reports observed vs. random purity/completeness with Z-scores.'
    )
    validate_parser.add_argument(
        '--recompute-xc', type=int, default=None, metavar='BIN_SIZE',
        help='Ignore stored XC tags; recompute cluster IDs from read coordinates using '
             'this bin size in bp (e.g. 5000, 10000, 20000, 50000, 100000). '
             'Use for bin-size sensitivity analysis.'
    )
    validate_parser.add_argument(
        '--negative-control-type', choices=['spatial', 'global', 'chrom'], default='spatial',
        help='Null model for --negative-control: "spatial" (default) shifts each read\'s '
             'midpoint by ±--neg-control-offset-bins bins while preserving chr/strand/exon '
             'lengths; "chrom" shuffles XC labels within each chromosome (recommended — '
             'destroys locus signal without gene-density confounding); '
             '"global" shuffles XC labels globally across all reads.'
    )
    validate_parser.add_argument(
        '--neg-control-offset-bins', type=int, default=50, metavar='N',
        help='Max per-read midpoint shift for spatial null, in units of --xc-position-quantum '
             '(default: 50 → ±50 kb at 1 kb quantum). Only used with --negative-control-type spatial.'
    )
    validate_parser.add_argument(
        '--xc-position-quantum', type=int, default=1000, metavar='BP',
        help='Midpoint bin size in bp used when computing spatial null XC values (default: 1000). '
             'Should match the --xc-position-quantum used during tagging.'
    )
    validate_parser.add_argument(
        '--xc-bin-size', type=int, default=10, metavar='BP',
        help='Exon-length bin size in bp used when computing spatial null XC values (default: 10). '
             'Should match the --xc-bin-size used during tagging.'
    )
    validate_parser.set_defaults(func=cmd_validate)

    args = parser.parse_args()

    if not hasattr(args, 'func'):
        parser.print_help()
        sys.exit(1)

    args.func(args)
