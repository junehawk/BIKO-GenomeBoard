#!/usr/bin/env python3
"""Build OMIM gene mapping SQLite from mim2gene.txt."""

import sqlite3
import logging
from pathlib import Path
from datetime import datetime, timezone

logger = logging.getLogger(__name__)
DEFAULT_DB_PATH = "data/db/omim_mapping.sqlite3"


def build_db(txt_path: str, db_path: str = DEFAULT_DB_PATH) -> str:
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)

    conn.execute("""CREATE TABLE IF NOT EXISTS mim2gene (
        mim_number TEXT NOT NULL,
        entry_type TEXT,
        entrez_id TEXT,
        gene_symbol TEXT NOT NULL,
        ensembl_id TEXT
    )""")
    conn.execute("DELETE FROM mim2gene")

    with open(txt_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            gene_symbol = parts[3].strip()
            if not gene_symbol:
                continue
            conn.execute(
                "INSERT INTO mim2gene VALUES (?,?,?,?,?)",
                (parts[0], parts[1], parts[2], gene_symbol, parts[4] if len(parts) > 4 else ""),
            )

    conn.execute("CREATE INDEX IF NOT EXISTS idx_omim_gene ON mim2gene(gene_symbol)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_omim_mim ON mim2gene(mim_number)")

    now = datetime.now(timezone.utc).isoformat()
    conn.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (now,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('source', 'OMIM (omim.org)')")
    count = conn.execute("SELECT COUNT(*) FROM mim2gene").fetchone()[0]
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('gene_count', ?)", (str(count),))

    conn.commit()
    conn.close()
    logger.info(f"OMIM mapping DB built: {count} gene entries → {db_path}")
    return db_path


if __name__ == "__main__":
    import argparse

    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument("txt_path", nargs="?", default="data/db/mim2gene.txt")
    parser.add_argument("--db-path", default=DEFAULT_DB_PATH)
    args = parser.parse_args()
    build_db(args.txt_path, args.db_path)
