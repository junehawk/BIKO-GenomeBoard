#!/usr/bin/env python3
"""Build local HPO gene-phenotype SQLite database from HPO annotation file."""

import logging
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

logger = logging.getLogger(__name__)

DEFAULT_TSV_URL = "https://hpo.jax.org/data/annotations/genes_to_phenotype.txt"
DEFAULT_DB_PATH = "data/db/hpo.sqlite3"


def build_db(tsv_path: str, db_path: str = DEFAULT_DB_PATH) -> str:
    """Build SQLite DB from HPO genes_to_phenotype.txt."""
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)

    conn.execute("""CREATE TABLE IF NOT EXISTS gene_phenotype (
        gene_symbol TEXT NOT NULL,
        hpo_id TEXT NOT NULL,
        hpo_name TEXT,
        ncbi_gene_id TEXT,
        frequency TEXT,
        disease_id TEXT
    )""")

    conn.execute("DELETE FROM gene_phenotype")

    with open(tsv_path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("ncbi_gene_id"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            conn.execute(
                "INSERT INTO gene_phenotype VALUES (?,?,?,?,?,?)",
                (
                    parts[1],  # gene_symbol
                    parts[2],  # hpo_id
                    parts[3],  # hpo_name
                    parts[0],  # ncbi_gene_id
                    parts[4] if len(parts) > 4 else "",
                    parts[5] if len(parts) > 5 else "",
                ),
            )

    conn.execute("CREATE INDEX IF NOT EXISTS idx_gp_gene ON gene_phenotype(gene_symbol)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_gp_hpo ON gene_phenotype(hpo_id)")

    now = datetime.now(timezone.utc).isoformat()
    conn.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (now,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('source', 'HPO (hpo.jax.org)')")
    count = conn.execute("SELECT COUNT(*) FROM gene_phenotype").fetchone()[0]
    gene_count = conn.execute("SELECT COUNT(DISTINCT gene_symbol) FROM gene_phenotype").fetchone()[0]
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('row_count', ?)", (str(count),))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('gene_count', ?)", (str(gene_count),))

    conn.commit()
    conn.close()
    logger.info(f"HPO DB built: {gene_count} genes, {count} associations → {db_path}")
    return db_path


if __name__ == "__main__":
    import argparse

    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument("tsv_path", nargs="?", default="data/db/genes_to_phenotype.txt")
    parser.add_argument("--db-path", default=DEFAULT_DB_PATH)
    args = parser.parse_args()
    build_db(args.tsv_path, args.db_path)
