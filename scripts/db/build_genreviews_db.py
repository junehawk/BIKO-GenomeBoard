#!/usr/bin/env python3
"""Build GeneReviews SQLite from NCBI FTP identifier files."""

import sqlite3
import logging
from pathlib import Path
from datetime import datetime

logger = logging.getLogger(__name__)
DEFAULT_DB_PATH = "data/db/genreviews.sqlite3"


def build_db(genes_path: str, titles_path: str, db_path: str = DEFAULT_DB_PATH) -> str:
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)

    conn.execute("""CREATE TABLE IF NOT EXISTS genreviews (
        gene_symbol TEXT NOT NULL,
        nbk_id TEXT NOT NULL,
        gr_shortname TEXT,
        disease_name TEXT,
        title TEXT,
        pmid TEXT
    )""")
    conn.execute("DELETE FROM genreviews")

    # Load titles (shortname → title + PMID)
    titles = {}
    with open(titles_path, encoding="latin-1") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                titles[parts[0]] = {"title": parts[1], "nbk_id": parts[2], "pmid": parts[3]}

    # Load gene-disease mappings
    with open(genes_path, encoding="latin-1") as f:
        for line in f:
            if line.startswith("#"):
                continue
            # Genes file uses pipe delimiter; titles file uses tab
            parts = line.strip().split("|")
            if len(parts) < 4:
                continue
            shortname, nbk_id, gene, disease = parts[0], parts[1], parts[2], parts[3]
            title_info = titles.get(shortname, {})
            conn.execute(
                "INSERT INTO genreviews VALUES (?,?,?,?,?,?)",
                (gene, nbk_id, shortname, disease,
                 title_info.get("title", disease),
                 title_info.get("pmid", "")),
            )

    conn.execute("CREATE INDEX IF NOT EXISTS idx_gr_gene ON genreviews(gene_symbol)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_gr_nbk ON genreviews(nbk_id)")

    now = datetime.utcnow().isoformat()
    conn.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (now,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('source', 'GeneReviews (NCBI)')")
    gene_count = conn.execute("SELECT COUNT(DISTINCT gene_symbol) FROM genreviews").fetchone()[0]
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('gene_count', ?)", (str(gene_count),))

    conn.commit()
    conn.close()
    logger.info(f"GeneReviews DB built: {gene_count} genes → {db_path}")
    return db_path


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", default="data/db/GRshortname_NBKid_genesymbol_dzname.txt")
    parser.add_argument("--titles", default="data/db/GRtitle_shortname_NBKid.txt")
    parser.add_argument("--db-path", default=DEFAULT_DB_PATH)
    args = parser.parse_args()
    build_db(args.genes, args.titles, args.db_path)
