#!/usr/bin/env python3
"""Build local ClinGen gene-validity SQLite database."""

import csv
import logging
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

logger = logging.getLogger(__name__)

DEFAULT_DB_PATH = "data/db/clingen.sqlite3"


def build_db(csv_path: str, db_path: str = DEFAULT_DB_PATH) -> str:
    """Build SQLite DB from ClinGen gene-validity CSV export."""
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)

    conn.execute("""CREATE TABLE IF NOT EXISTS gene_validity (
        gene_symbol TEXT NOT NULL,
        hgnc_id TEXT,
        disease TEXT,
        mondo_id TEXT,
        moi TEXT,
        sop TEXT,
        classification TEXT NOT NULL,
        report_url TEXT,
        classification_date TEXT,
        gcep TEXT
    )""")
    conn.execute("DELETE FROM gene_validity")

    with open(csv_path) as f:
        # Skip header lines until we find the column header
        header_cols = None
        for line in f:
            if line.startswith('"GENE SYMBOL"'):
                header_cols = next(csv.reader([line]))
                break
        if header_cols is None:
            logger.error("ClinGen CSV header not found — expected '\"GENE SYMBOL\",...'")
            conn.close()
            return db_path

        # Build column index map so we're resilient to column reordering.
        # Fall back to positional parsing if header names don't match.
        col_map = {name.strip().upper(): idx for idx, name in enumerate(header_cols)}

        def _col(name: str, fallback: int) -> int:
            return col_map.get(name, fallback)

        i_gene = _col("GENE SYMBOL", 0)
        i_hgnc = _col("GENE ID (HGNC)", 1)
        i_disease = _col("DISEASE LABEL", 2)
        i_mondo = _col("DISEASE ID (MONDO)", 3)
        i_moi = _col("MOI", 4)
        i_sop = _col("SOP", 5)
        i_class = _col("CLASSIFICATION", 6)
        i_url = _col("ONLINE REPORT", 7)
        i_date = _col("CLASSIFICATION DATE", 8)
        i_gcep = _col("GCEP", 9)

        reader = csv.reader(f)
        for row in reader:
            if len(row) < 7 or row[0].startswith("+"):
                continue
            conn.execute(
                "INSERT INTO gene_validity VALUES (?,?,?,?,?,?,?,?,?,?)",
                (
                    row[i_gene] if i_gene < len(row) else "",
                    row[i_hgnc] if i_hgnc < len(row) else "",
                    row[i_disease] if i_disease < len(row) else "",
                    row[i_mondo] if i_mondo < len(row) else "",
                    row[i_moi] if i_moi < len(row) else "",
                    row[i_sop] if i_sop < len(row) else "",
                    row[i_class] if i_class < len(row) else "",
                    row[i_url] if i_url < len(row) else "",
                    row[i_date] if i_date < len(row) else "",
                    row[i_gcep] if i_gcep < len(row) else "",
                ),
            )

    conn.execute("CREATE INDEX IF NOT EXISTS idx_cv_gene ON gene_validity(gene_symbol)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_cv_class ON gene_validity(classification)")

    now = datetime.now(timezone.utc).isoformat()
    conn.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (now,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('source', 'ClinGen (clinicalgenome.org)')")
    count = conn.execute("SELECT COUNT(*) FROM gene_validity").fetchone()[0]
    gene_count = conn.execute("SELECT COUNT(DISTINCT gene_symbol) FROM gene_validity").fetchone()[0]
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('row_count', ?)", (str(count),))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('gene_count', ?)", (str(gene_count),))

    conn.commit()
    conn.close()
    logger.info(f"ClinGen DB built: {gene_count} genes, {count} curations → {db_path}")
    return db_path


if __name__ == "__main__":
    import argparse

    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path")
    parser.add_argument("--db-path", default=DEFAULT_DB_PATH)
    args = parser.parse_args()
    build_db(args.csv_path, args.db_path)
