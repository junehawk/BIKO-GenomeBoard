#!/usr/bin/env python3
"""Build local ClinVar SQLite database from NCBI variant_summary.txt.gz"""

import gzip
import sqlite3
import logging
import sys
import time
from pathlib import Path
from datetime import datetime

logger = logging.getLogger(__name__)

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
DEFAULT_DB_PATH = "data/db/clinvar.sqlite3"

def download_clinvar(output_path: str) -> str:
    """Download ClinVar variant_summary.txt.gz"""
    import requests
    logger.info(f"Downloading ClinVar from {CLINVAR_URL}...")
    resp = requests.get(CLINVAR_URL, stream=True, timeout=300)
    resp.raise_for_status()
    with open(output_path, 'wb') as f:
        for chunk in resp.iter_content(chunk_size=8192):
            f.write(chunk)
    logger.info(f"Downloaded to {output_path}")
    return output_path

def build_db(tsv_gz_path: str, db_path: str = DEFAULT_DB_PATH) -> str:
    """Build SQLite DB from variant_summary.txt.gz

    Key columns from variant_summary.txt:
    #AlleleID, Type, Name, GeneID, GeneSymbol, HGNC_ID, ClinicalSignificance,
    ClinSigSimple, LastEvaluated, RS# (dbSNP), nsv/esv (dbVar), RCVaccession,
    PhenotypeIDS, PhenotypeList, Origin, OriginSimple, Assembly, ChromosomeAccession,
    Chromosome, Start, Stop, ReferenceAllele, AlternateAllele, Cytogenetic,
    ReviewStatus, NumberSubmitters, Guidelines, TestedInGTR, OtherIDs,
    SubmitterCategories, VariationID, PositionVCF, ReferenceAlleleVCF, AlternateAlleleVCF
    """
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")

    # Create table
    conn.execute("""
        CREATE TABLE IF NOT EXISTS variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            rsid TEXT,
            gene TEXT,
            clinical_significance TEXT,
            review_status TEXT,
            phenotype_list TEXT,
            variation_id TEXT,
            allele_id TEXT,
            origin TEXT,
            assembly TEXT DEFAULT 'GRCh38',
            last_evaluated TEXT,
            number_submitters INTEGER
        )
    """)

    # Create metadata table
    conn.execute("""
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT
        )
    """)

    conn.execute("DELETE FROM variants")  # Fresh build

    logger.info(f"Building ClinVar DB from {tsv_gz_path}...")
    count = 0
    skipped = 0

    open_func = gzip.open if tsv_gz_path.endswith('.gz') else open
    mode = 'rt' if tsv_gz_path.endswith('.gz') else 'r'

    with open_func(tsv_gz_path, mode) as f:
        header = None
        batch = []
        for line in f:
            if line.startswith('#'):
                # Parse header
                header = line.strip('#').strip().split('\t')
                continue
            if header is None:
                continue

            fields = line.strip().split('\t')
            if len(fields) < 33:
                skipped += 1
                continue

            row = dict(zip(header, fields))

            # Only keep GRCh38 assembly
            assembly = row.get('Assembly', '')
            if assembly != 'GRCh38':
                continue

            chrom = row.get('Chromosome', '')
            pos_str = row.get('PositionVCF', '')
            ref = row.get('ReferenceAlleleVCF', '')
            alt = row.get('AlternateAlleleVCF', '')

            if not chrom or not pos_str or pos_str == '-1' or not ref or not alt:
                skipped += 1
                continue

            try:
                pos = int(pos_str)
            except ValueError:
                skipped += 1
                continue

            rsid_raw = row.get('RS# (dbSNP)', '')
            rsid = f"rs{rsid_raw}" if rsid_raw and rsid_raw != '-1' else None

            batch.append((
                f"chr{chrom}" if not chrom.startswith('chr') else chrom,
                pos,
                ref,
                alt,
                rsid,
                row.get('GeneSymbol', ''),
                row.get('ClinicalSignificance', ''),
                row.get('ReviewStatus', ''),
                row.get('PhenotypeList', ''),
                row.get('VariationID', ''),
                row.get('#AlleleID', row.get('AlleleID', '')),
                row.get('OriginSimple', ''),
                assembly,
                row.get('LastEvaluated', ''),
                int(row.get('NumberSubmitters', '0') or '0'),
            ))

            count += 1
            if len(batch) >= 10000:
                conn.executemany("""
                    INSERT INTO variants (chrom, pos, ref, alt, rsid, gene,
                        clinical_significance, review_status, phenotype_list,
                        variation_id, allele_id, origin, assembly, last_evaluated,
                        number_submitters)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, batch)
                batch = []
                if count % 100000 == 0:
                    logger.info(f"  Processed {count:,} variants...")

    if batch:
        conn.executemany("""
            INSERT INTO variants (chrom, pos, ref, alt, rsid, gene,
                clinical_significance, review_status, phenotype_list,
                variation_id, allele_id, origin, assembly, last_evaluated,
                number_submitters)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, batch)

    # Create indexes
    logger.info("Creating indexes...")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_chrom_pos ON variants(chrom, pos)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_rsid ON variants(rsid)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_gene ON variants(gene)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_chrom_pos_ref_alt ON variants(chrom, pos, ref, alt)")

    # Store metadata
    now = datetime.utcnow().isoformat()
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (now,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('source', ?)", (CLINVAR_URL,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('variant_count', ?)", (str(count),))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('assembly', 'GRCh38')")

    conn.commit()
    conn.close()

    logger.info(f"ClinVar DB built: {count:,} variants, {skipped:,} skipped → {db_path}")
    return db_path

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    if len(sys.argv) > 1:
        # Build from existing file
        build_db(sys.argv[1])
    else:
        # Download and build
        Path("data/db").mkdir(parents=True, exist_ok=True)
        gz_path = "data/db/variant_summary.txt.gz"
        download_clinvar(gz_path)
        build_db(gz_path)
