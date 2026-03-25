#!/usr/bin/env python3
"""Build local ClinGen gene-validity SQLite database."""

import csv
import sqlite3
import logging
from pathlib import Path
from datetime import datetime

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
        sop TEXT,
        classification TEXT NOT NULL,
        report_url TEXT,
        classification_date TEXT,
        gcep TEXT
    )""")
    conn.execute("DELETE FROM gene_validity")

    with open(csv_path) as f:
        # Skip header lines until we find the column header
        for line in f:
            if line.startswith('"GENE SYMBOL"'):
                break
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 6:
                continue
            conn.execute(
                "INSERT INTO gene_validity VALUES (?,?,?,?,?,?,?,?,?)",
                (
                    row[0],  # gene_symbol
                    row[1],  # hgnc_id
                    row[2],  # disease
                    row[3],  # mondo_id
                    row[4],  # sop
                    row[5],  # classification
                    row[6] if len(row) > 6 else "",
                    row[7] if len(row) > 7 else "",
                    row[8] if len(row) > 8 else "",
                ),
            )

    conn.execute("CREATE INDEX IF NOT EXISTS idx_cv_gene ON gene_validity(gene_symbol)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_cv_class ON gene_validity(classification)")

    now = datetime.utcnow().isoformat()
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
