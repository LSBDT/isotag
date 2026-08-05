#!/usr/bin/env python3
"""
IsoTag XC GENCODE Version Comparison

Tests XC cluster stability across GENCODE versions (v42, v44, v45, v47).
Core question: do XC clusters map to the same gene regardless of which
annotation version is used? This validates XC as annotation-independent.

Metrics:
- Stability: % of XC clusters assigned to same gene in all versions
- Gain/Loss: clusters that gain or lose gene assignment across versions
- Remap: clusters reassigned to different gene between versions
- Never: clusters with no gene assignment in any version

Input:
- BAM file with XC tags
- 2-4 GENCODE GTF files specified as VERSION=PATH pairs

Output:
- xc_gencode_compare_summary.tsv — per-cluster classification
- xc_gencode_compare_matrix.tsv — pairwise version agreement matrix
- Console summary: stability %, gain/loss counts per version step

Memory: SQLite-backed, streaming bedtools — constant RAM regardless of BAM size.
"""

import argparse
import os
import re
import sys
import logging
import sqlite3
import tempfile
from pathlib import Path
import subprocess

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def normalize_gene_id(gene_id: str) -> str:
    """
    Strip GENCODE version suffix (.N) from gene IDs.

    ENSG00000290385.1 and ENSG00000290385.2 are the same gene model
    at different GENCODE update counts — not different genes.
    Without stripping, cross-version comparison falsely reports 'remapped'.
    """
    if '.' in gene_id:
        return gene_id.rsplit('.', 1)[0]
    return gene_id


def extract_genes_as_bed(gtf_path: str, bed_path: str) -> int:
    """
    Pre-extract gene-only BED6 from GTF to avoid loading full GTF interval tree.
    Full GTF (~3M features) causes bedtools OOM at 20GB; gene-only BED6 (~60K) is ~50x smaller.
    """
    count = 0
    with open(gtf_path) as gtf, open(bed_path, 'w') as bed:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            chrom = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]
            m = re.search(r'gene_id "([^"]+)"', fields[8])
            if not m:
                continue
            gene_id = m.group(1).split('.')[0]
            bed.write(f'{chrom}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n')
            count += 1
    return count


def create_compare_db(db_path: str) -> sqlite3.Connection:
    """
    Create SQLite database for multi-version comparison.

    Tables:
        read_xc      — read_id → xc tag (one row per primary read with XC)
        read_gene_v  — read_id × version → gene_id (one row per overlap)
    """
    conn = sqlite3.connect(db_path)
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
        CREATE TABLE read_gene_v (
            read_id TEXT NOT NULL,
            version TEXT NOT NULL,
            gene_id TEXT NOT NULL
        )
    ''')
    return conn


def stream_xc_to_db(bam_path: str, conn: sqlite3.Connection):
    """
    Stream XC tags from BAM into SQLite read_xc table via samtools view subprocess.
    Uses -F 0x900 to skip secondary/supplementary, -T SAM output parsed line-by-line.
    Avoids pysam page-cache OOM for large BAMs (16.9GB MPC_hg38_v11.bam).
    """
    logger.info(f"Extracting XC tags from {bam_path}")
    n_xc = 0
    total = 0

    cmd = f"samtools view -F 0x900 {bam_path}"
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, text=True)

    conn.execute('BEGIN')
    for line in proc.stdout:
        if line.startswith('@'):
            continue
        total += 1
        # SAM optional fields after column 11: search for XC:Z:HASH
        xc_tag = None
        for field in line.rstrip('\n').split('\t')[11:]:
            if field.startswith('XC:Z:'):
                xc_tag = field[5:]
                break
        if xc_tag:
            read_name = line.split('\t', 1)[0]
            conn.execute(
                'INSERT OR REPLACE INTO read_xc (read_id, xc) VALUES (?, ?)',
                (read_name, xc_tag)
            )
            n_xc += 1
        if total % 50000 == 0:
            conn.execute('COMMIT')
            conn.execute('BEGIN')
            logger.info(f"Processed {total:,} reads, {n_xc:,} with XC tags")
    conn.execute('COMMIT')
    proc.wait()
    if proc.returncode != 0:
        logger.error(f"samtools view failed with returncode {proc.returncode}")
        sys.exit(1)

    logger.info(
        f"Total: {total:,} primary reads, {n_xc:,} with XC tags "
        f"({n_xc / total * 100:.1f}% of primary)" if total else "Total: 0 primary reads"
    )


def intersect_bam_with_gtf(bam_path: str, gtf_path: str, output_bed: str):
    """
    Intersect BAM with gene-only BED6 (pre-extracted from GTF) using bedtools.
    Pre-extraction avoids loading full GTF (~3M features) interval tree → OOM fix.
    Same strand (-s) to avoid antisense overlaps.
    """
    logger.info(f"Intersecting BAM with GTF: {gtf_path}")
    with tempfile.NamedTemporaryFile(mode='w', suffix='_genes.bed', delete=False) as tmp:
        gene_bed = tmp.name
    try:
        n_genes = extract_genes_as_bed(gtf_path, gene_bed)
        logger.info(f"Extracted {n_genes:,} gene intervals from GTF as BED6")

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
        if p_isect.returncode != 0:
            logger.error(f"bedtools failed: {isect_stderr}")
            sys.exit(1)
    finally:
        if os.path.exists(gene_bed):
            os.unlink(gene_bed)

    with open(output_bed) as f:
        n_overlaps = sum(1 for _ in f)
    logger.info(f"Found {n_overlaps:,} read-gene overlaps")


def stream_intersect_to_db_versioned(
    bed_path: str, version: str, conn: sqlite3.Connection
):
    """
    Stream bedtools intersect output into SQLite read_gene_v table.

    Since intersect_bam_with_gtf now uses a gene-only BED6 (not full GTF),
    the output has 18 columns: BED12 (read) + BED6 (gene).
        BED12 (12): chr start end read_id score strand thickS thickE rgb blockN sizes starts
        BED6  ( 6): chr start end gene_id score strand
    read_id = fields[3]; gene_id = fields[15]  (index 12+3=15 in BED6)
    """
    logger.info(f"Loading intersection results for version {version} from {bed_path}")
    count = 0
    conn.execute('BEGIN')

    with open(bed_path) as f:
        for i, line in enumerate(f):
            fields = line.strip().split('\t')
            if len(fields) < 16:
                continue

            read_id = fields[3]
            gene_id = fields[15]  # BED6 name column = gene_id (pre-stripped of version)

            if gene_id:
                conn.execute(
                    'INSERT INTO read_gene_v (read_id, version, gene_id) VALUES (?, ?, ?)',
                    (read_id, version, normalize_gene_id(gene_id))
                )
                count += 1

            if (i + 1) % 50000 == 0:
                conn.execute('COMMIT')
                conn.execute('BEGIN')
                logger.info(
                    f"[{version}] Parsed {i + 1:,} overlap lines, {count:,} assignments"
                )

    conn.execute('COMMIT')
    logger.info(f"[{version}] Total: {count:,} read-gene pairs loaded")


def build_indexes(conn: sqlite3.Connection):
    """Build indexes after all data is loaded — faster than incremental indexing."""
    logger.info("Building database indexes")
    conn.execute('CREATE INDEX idx_rgv_read ON read_gene_v(read_id)')
    conn.execute('CREATE INDEX idx_rgv_ver ON read_gene_v(version)')
    conn.execute('CREATE INDEX idx_rx_xc ON read_xc(xc)')
    logger.info("Indexes built")


def build_xc_gene_view(conn: sqlite3.Connection, versions: list):
    """
    Build derived table: for each (xc, version), what gene(s) do reads map to?

    xc_gene_v: xc × version → gene_id × n_reads
    xc_versions: xc → set of versions where gene was assigned (as pipe-separated string)
    xc_stability: xc → stability classification
    """
    logger.info("Building XC-gene-version summary table")
    conn.execute('DROP TABLE IF EXISTS xc_gene_v')
    conn.execute('''
        CREATE TEMP TABLE xc_gene_v AS
        SELECT rx.xc,
               rgv.version,
               rgv.gene_id,
               COUNT(DISTINCT rx.read_id) AS n_reads
        FROM read_xc rx
        JOIN read_gene_v rgv ON rx.read_id = rgv.read_id
        GROUP BY rx.xc, rgv.version, rgv.gene_id
    ''')

    # Per-cluster version summary: how many distinct genes per version?
    conn.execute('DROP TABLE IF EXISTS xc_version_summary')
    conn.execute('''
        CREATE TEMP TABLE xc_version_summary AS
        SELECT xc,
               version,
               COUNT(DISTINCT gene_id) AS n_genes,
               GROUP_CONCAT(DISTINCT gene_id) AS genes
        FROM xc_gene_v
        GROUP BY xc, version
    ''')

    # Per-cluster overall: how many versions assigned, how many distinct genes?
    conn.execute('DROP TABLE IF EXISTS xc_overall')
    conn.execute('''
        CREATE TEMP TABLE xc_overall AS
        SELECT xc,
               COUNT(DISTINCT version)  AS n_versions_assigned,
               COUNT(DISTINCT gene_id)  AS n_genes_total
        FROM xc_gene_v
        GROUP BY xc
    ''')

    logger.info(
        f"XC clusters with any gene assignment: "
        f"{conn.execute('SELECT COUNT(*) FROM xc_overall').fetchone()[0]:,}"
    )


def classify_xc_clusters(conn: sqlite3.Connection, versions: list) -> dict:
    """
    Classify each XC cluster by stability across versions.

    Classifications:
        stable      — assigned to same single gene in all versions where assigned,
                      AND assigned in ≥2 versions
        stable_all  — subset of stable: assigned in ALL n versions
        remapped    — assigned in ≥2 versions but to different genes
        gained      — assigned only in later version(s), not in earliest
        lost        — assigned in earliest version(s), not in latest
        singleton   — assigned in exactly 1 version only
        never       — no gene assignment in any version

    Returns dict of summary counts + per-cluster details list.
    """
    logger.info("Classifying XC clusters by stability")

    n_versions = len(versions)
    earliest = versions[0]
    latest = versions[-1]

    # All XC clusters (including those with no gene assignment)
    all_xc = {r[0] for r in conn.execute('SELECT DISTINCT xc FROM read_xc')}
    assigned_xc = {r[0] for r in conn.execute('SELECT DISTINCT xc FROM xc_overall')}
    never_xc = all_xc - assigned_xc

    clusters = []

    for (xc, n_ver, n_genes_total) in conn.execute(
        'SELECT xc, n_versions_assigned, n_genes_total FROM xc_overall ORDER BY xc'
    ):
        # Get versions where this XC was assigned
        ver_genes = dict(conn.execute(
            'SELECT version, genes FROM xc_version_summary WHERE xc = ? ORDER BY version',
            (xc,)
        ).fetchall())

        in_earliest = earliest in ver_genes
        in_latest = latest in ver_genes

        # Collect all unique gene sets per version
        gene_sets = [frozenset(v.split(',')) for v in ver_genes.values()]
        all_genes_same = (len(set(frozenset(gs) for gs in gene_sets)) == 1)

        # Each version has only 1 gene?
        all_single = all(len(gs) == 1 for gs in gene_sets)

        if n_ver == 1:
            label = 'singleton'
        elif all_genes_same and all_single:
            # Same single gene in all assigned versions
            if n_ver == n_versions:
                label = 'stable_all'
            else:
                label = 'stable'
        elif not in_earliest and in_latest:
            label = 'gained'
        elif in_earliest and not in_latest:
            label = 'lost'
        else:
            label = 'remapped'

        # Consensus gene: the gene assigned in the most versions
        gene_version_count: dict = {}
        for gs in gene_sets:
            for g in gs:
                gene_version_count[g] = gene_version_count.get(g, 0) + 1
        consensus_gene = max(gene_version_count, key=lambda g: gene_version_count[g])

        # Build per-version assignment string: v42=GENE1,v44=GENE2,...
        ver_str = '|'.join(
            f"{v}={ver_genes[v]}" if v in ver_genes else f"{v}=none"
            for v in versions
        )

        clusters.append({
            'xc':             xc,
            'label':          label,
            'n_versions':     n_ver,
            'n_genes_total':  n_genes_total,
            'consensus_gene': consensus_gene,
            'per_version':    ver_str,
        })

    # Add 'never' clusters
    for xc in sorted(never_xc):
        clusters.append({
            'xc':             xc,
            'label':          'never',
            'n_versions':     0,
            'n_genes_total':  0,
            'consensus_gene': 'none',
            'per_version':    '|'.join(f"{v}=none" for v in versions),
        })

    # Count per label
    counts: dict = {}
    for c in clusters:
        counts[c['label']] = counts.get(c['label'], 0) + 1

    total = len(clusters)
    assigned = total - counts.get('never', 0)
    # Stable = stable_all + stable
    n_stable = counts.get('stable_all', 0) + counts.get('stable', 0)

    return {
        'total_xc_clusters':   total,
        'assigned_clusters':   assigned,
        'never_clusters':      counts.get('never', 0),
        'stable_all':          counts.get('stable_all', 0),
        'stable':              counts.get('stable', 0),
        'stable_total':        n_stable,
        'singleton':           counts.get('singleton', 0),
        'remapped':            counts.get('remapped', 0),
        'gained':              counts.get('gained', 0),
        'lost':                counts.get('lost', 0),
        'stability_pct':       n_stable / assigned * 100 if assigned else 0,
        'stability_all_pct':   counts.get('stable_all', 0) / assigned * 100 if assigned else 0,
        'clusters':            clusters,
    }


def calculate_version_matrix(conn: sqlite3.Connection, versions: list) -> list:
    """
    Compute pairwise agreement between consecutive version pairs.

    For each version pair (vA, vB):
        - n_both: XC clusters assigned in both versions
        - n_agree: of those, same gene assignment
        - pct_agree: n_agree / n_both * 100
        - n_only_a: assigned in vA but not vB (loss)
        - n_only_b: assigned in vB but not vA (gain)
    """
    matrix = []
    version_pairs = [(versions[i], versions[i + 1]) for i in range(len(versions) - 1)]

    for v_a, v_b in version_pairs:
        # Clusters assigned in both versions
        both = conn.execute('''
            SELECT a.xc
            FROM xc_version_summary a
            JOIN xc_version_summary b ON a.xc = b.xc
            WHERE a.version = ? AND b.version = ?
        ''', (v_a, v_b)).fetchall()
        n_both = len(both)

        if n_both == 0:
            matrix.append({
                'v_a': v_a, 'v_b': v_b,
                'n_both': 0, 'n_agree': 0, 'pct_agree': 0.0,
                'n_only_a': 0, 'n_only_b': 0,
            })
            continue

        # Agreement: both assigned to same single gene
        n_agree = conn.execute('''
            SELECT COUNT(*)
            FROM xc_version_summary a
            JOIN xc_version_summary b ON a.xc = b.xc
            WHERE a.version = ? AND b.version = ?
              AND a.genes = b.genes
              AND a.n_genes = 1 AND b.n_genes = 1
        ''', (v_a, v_b)).fetchone()[0]

        # Only in vA (gene lost in vB)
        n_only_a = conn.execute('''
            SELECT COUNT(*) FROM xc_version_summary
            WHERE version = ?
              AND xc NOT IN (SELECT xc FROM xc_version_summary WHERE version = ?)
        ''', (v_a, v_b)).fetchone()[0]

        # Only in vB (gene gained in vB)
        n_only_b = conn.execute('''
            SELECT COUNT(*) FROM xc_version_summary
            WHERE version = ?
              AND xc NOT IN (SELECT xc FROM xc_version_summary WHERE version = ?)
        ''', (v_b, v_a)).fetchone()[0]

        matrix.append({
            'v_a':      v_a,
            'v_b':      v_b,
            'n_both':   n_both,
            'n_agree':  n_agree,
            'pct_agree': n_agree / n_both * 100 if n_both else 0.0,
            'n_only_a': n_only_a,
            'n_only_b': n_only_b,
        })

    return matrix


def write_comparison_reports(
    metrics: dict,
    matrix: list,
    versions: list,
    output_dir: Path,
):
    """
    Write two output files:
        xc_gencode_compare_summary.tsv — per-cluster classification
        xc_gencode_compare_matrix.tsv  — pairwise version agreement
    """
    # --- Summary file ---
    summary_path = output_dir / 'xc_gencode_compare_summary.tsv'
    logger.info(f"Writing per-cluster summary to {summary_path}")

    with open(summary_path, 'w') as out:
        out.write("# IsoTag XC GENCODE Version Comparison — Per-Cluster Summary\n")
        out.write(f"# Versions: {', '.join(versions)}\n")
        out.write("#\n")
        out.write("# Labels:\n")
        out.write("#   stable_all  — same single gene in ALL versions\n")
        out.write("#   stable      — same single gene in all ASSIGNED versions (not all)\n")
        out.write("#   singleton   — gene assigned in exactly 1 version\n")
        out.write("#   remapped    — assigned in ≥2 versions but to different genes\n")
        out.write("#   gained      — not assigned in earliest version, assigned in latest\n")
        out.write("#   lost        — assigned in earliest version, not in latest\n")
        out.write("#   never       — no gene assignment in any version\n")
        out.write("#\n")
        # Header
        ver_cols = '\t'.join(versions)
        out.write(f"XC_Tag\tLabel\tN_Versions_Assigned\tN_Genes_Total\tConsensus_Gene\t{ver_cols}\n")

        for c in metrics['clusters']:
            # Parse per_version into per-column values
            ver_map = {}
            for token in c['per_version'].split('|'):
                v, g = token.split('=', 1)
                ver_map[v] = g
            ver_vals = '\t'.join(ver_map.get(v, 'none') for v in versions)
            out.write(
                f"{c['xc']}\t{c['label']}\t{c['n_versions']}\t"
                f"{c['n_genes_total']}\t{c['consensus_gene']}\t{ver_vals}\n"
            )

    # --- Matrix file ---
    matrix_path = output_dir / 'xc_gencode_compare_matrix.tsv'
    logger.info(f"Writing pairwise version matrix to {matrix_path}")

    with open(matrix_path, 'w') as out:
        out.write("# IsoTag XC GENCODE Version Comparison — Pairwise Agreement Matrix\n")
        out.write(f"# Versions compared: {', '.join(versions)}\n")
        out.write("#\n")
        out.write("Version_A\tVersion_B\t"
                  "N_Clusters_Both\tN_Agree\tPct_Agree\tN_Only_A\tN_Only_B\n")
        for row in matrix:
            out.write(
                f"{row['v_a']}\t{row['v_b']}\t"
                f"{row['n_both']}\t{row['n_agree']}\t{row['pct_agree']:.1f}\t"
                f"{row['n_only_a']}\t{row['n_only_b']}\n"
            )

    return summary_path, matrix_path


def print_console_summary(metrics: dict, matrix: list, versions: list):
    """Print human-readable summary to stdout."""
    print("\n" + "=" * 60)
    print("XC GENCODE VERSION STABILITY RESULTS")
    print("=" * 60)
    print(f"Versions analyzed:     {', '.join(versions)}")
    print(f"Total XC clusters:     {metrics['total_xc_clusters']:,}")
    print(f"  With gene assignment: {metrics['assigned_clusters']:,}")
    print(f"  No gene (any ver):    {metrics['never_clusters']:,}")
    print()
    print("Cluster classification (of assigned clusters):")
    print(f"  stable_all  (same gene, ALL versions):     "
          f"{metrics['stable_all']:>6,}  ({metrics['stability_all_pct']:.1f}%)")
    print(f"  stable      (same gene, assigned versions): "
          f"{metrics['stable']:>5,}")
    print(f"  --- Total stable: {metrics['stable_total']:,}  "
          f"({metrics['stability_pct']:.1f}% of assigned) ---")
    print(f"  singleton   (assigned in 1 version only):  "
          f"{metrics['singleton']:>6,}")
    print(f"  remapped    (different genes between vers): "
          f"{metrics['remapped']:>5,}")
    print(f"  gained      (not in earliest, in latest):  "
          f"{metrics['gained']:>6,}")
    print(f"  lost        (in earliest, not in latest):  "
          f"{metrics['lost']:>6,}")
    print()
    print("Version-to-version agreement:")
    for row in matrix:
        print(f"  {row['v_a']} → {row['v_b']}: "
              f"{row['n_agree']}/{row['n_both']} agree ({row['pct_agree']:.1f}%), "
              f"+{row['n_only_b']} gained, -{row['n_only_a']} lost")
    print("=" * 60 + "\n")


def cmd_compare(args):
    """Main entry point: compare XC stability across GENCODE versions."""

    # Parse --gtf VERSION=PATH arguments
    versions = []
    gtf_paths = {}
    for spec in args.gtf:
        if '=' not in spec:
            logger.error(f"--gtf argument must be VERSION=PATH, got: {spec!r}")
            sys.exit(1)
        version, path = spec.split('=', 1)
        if not os.path.exists(path):
            logger.error(f"GTF file not found: {path}")
            sys.exit(1)
        versions.append(version)
        gtf_paths[version] = path

    if len(versions) < 2:
        logger.error("At least 2 --gtf VERSION=PATH arguments required")
        sys.exit(1)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # SQLite temp file — use output_dir (not /tmp which has limited space)
    db_fd, db_path = tempfile.mkstemp(suffix='.db', prefix='isotag_xccompare_', dir=str(output_dir))
    os.close(db_fd)

    conn = None
    temp_beds = []
    try:
        conn = create_compare_db(db_path)

        # Step 1: Stream XC tags from BAM
        stream_xc_to_db(args.bam, conn)
        n_xc = conn.execute('SELECT COUNT(*) FROM read_xc').fetchone()[0]
        if n_xc == 0:
            logger.error("No XC tags found in BAM file")
            sys.exit(1)
        logger.info(f"Loaded {n_xc:,} reads with XC tags")

        # Step 2: For each version, intersect BAM with GTF and load into DB
        for version in versions:
            gtf_path = gtf_paths[version]
            bed_fd, bed_path = tempfile.mkstemp(suffix=f'_{version}.bed', prefix='isotag_xc_', dir=str(output_dir))
            os.close(bed_fd)
            temp_beds.append(bed_path)

            intersect_bam_with_gtf(args.bam, gtf_path, bed_path)
            stream_intersect_to_db_versioned(bed_path, version, conn)

        # Build indexes after all loading done
        build_indexes(conn)

        # Step 3: Build XC-gene-version summary
        build_xc_gene_view(conn, versions)

        # Step 4: Classify clusters
        metrics = classify_xc_clusters(conn, versions)

        # Step 5: Compute pairwise matrix
        matrix = calculate_version_matrix(conn, versions)

        # Step 6: Write output files
        summary_path, matrix_path = write_comparison_reports(
            metrics, matrix, versions, output_dir
        )

        # Step 7: Print console summary
        print_console_summary(metrics, matrix, versions)
        print("Detailed results:")
        print(f"  Per-cluster: {summary_path}")
        print(f"  Matrix:      {matrix_path}")

    finally:
        if conn:
            conn.close()
        if os.path.exists(db_path):
            os.unlink(db_path)
        for bed_path in temp_beds:
            if os.path.exists(bed_path):
                os.unlink(bed_path)
        logger.info("Cleaned up temporary files")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compare XC cluster stability across GENCODE annotation versions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare XC tags across GENCODE v42, v44, v45, v47
  python3 scripts/isotag_xc_gencode_compare.py compare \\
      --bam test_xc_default.bam \\
      --gtf v42=/path/to/gencode.v42.annotation.gtf \\
      --gtf v44=/path/to/gencode.v44.annotation.gtf \\
      --gtf v45=/path/to/gencode.v45.annotation.gtf \\
      --gtf v47=/path/to/gencode.v47.annotation.gtf \\
      --output-dir results/

  # Quick test with 2 versions
  python3 scripts/isotag_xc_gencode_compare.py compare \\
      --bam test_xc_tiny.bam \\
      --gtf v45=/path/to/gencode.v45.annotation.gtf \\
      --gtf v47=/path/to/gencode.v47.annotation.gtf \\
      --output-dir /tmp/
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Subcommands')

    compare_parser = subparsers.add_parser('compare', help='Run cross-version comparison')
    compare_parser.add_argument('--bam', required=True, help='BAM file with XC tags')
    compare_parser.add_argument(
        '--gtf', required=True, action='append', metavar='VERSION=PATH',
        help='GENCODE GTF as VERSION=PATH (repeat for each version, e.g. v42=/path/to/v42.gtf)'
    )
    compare_parser.add_argument(
        '--output-dir', default='results/', help='Output directory (default: results/)'
    )
    compare_parser.set_defaults(func=cmd_compare)

    args = parser.parse_args()

    if not hasattr(args, 'func'):
        parser.print_help()
        sys.exit(1)

    args.func(args)
