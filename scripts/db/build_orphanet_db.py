#!/usr/bin/env python3
"""Build Orphanet prevalence SQLite from Product 9 XML."""

import sqlite3
import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime, timezone

logger = logging.getLogger(__name__)
DEFAULT_DB_PATH = "data/db/orphanet.sqlite3"


def build_db(xml_path: str, db_path: str = DEFAULT_DB_PATH) -> str:
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)

    conn.execute("""CREATE TABLE IF NOT EXISTS prevalence (
        orpha_code TEXT,
        disease_name TEXT NOT NULL,
        gene_symbol TEXT,
        gene_name TEXT,
        prevalence_type TEXT,
        prevalence_class TEXT,
        val_moy REAL,
        geographic TEXT
    )""")
    conn.execute("DELETE FROM prevalence")

    tree = ET.parse(xml_path)
    root = tree.getroot()

    for disorder in root.iter("Disorder"):
        orpha = disorder.findtext("OrphaCode", "")
        name = disorder.findtext("Name", "")

        # Extract gene associations
        genes = []
        for ga in disorder.iter("DisorderGeneAssociation"):
            gene_el = ga.find("Gene")
            if gene_el is not None:
                genes.append(
                    {
                        "symbol": gene_el.findtext("Symbol", ""),
                        "name": gene_el.findtext("Name", ""),
                    }
                )

        # Extract prevalence data
        for prev in disorder.iter("Prevalence"):
            prev_type = ""
            pt = prev.find("PrevalenceType")
            if pt is not None:
                prev_type = pt.findtext("Name", "")
            prev_class = ""
            pc = prev.find("PrevalenceClass")
            if pc is not None:
                prev_class = pc.findtext("Name", "")
            val_moy = prev.findtext("ValMoy", "")
            geo = ""
            pg = prev.find("PrevalenceGeographic")
            if pg is not None:
                geo = pg.findtext("Name", "")

            for gene in genes or [{"symbol": "", "name": ""}]:
                conn.execute(
                    "INSERT INTO prevalence VALUES (?,?,?,?,?,?,?,?)",
                    (
                        orpha,
                        name,
                        gene["symbol"],
                        gene["name"],
                        prev_type,
                        prev_class,
                        float(val_moy) if val_moy else None,
                        geo,
                    ),
                )

    conn.execute("CREATE INDEX IF NOT EXISTS idx_prev_gene ON prevalence(gene_symbol)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_prev_orpha ON prevalence(orpha_code)")

    now = datetime.now(timezone.utc).isoformat()
    conn.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (now,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('source', 'Orphanet (orphadata.com)')")
    count = conn.execute("SELECT COUNT(*) FROM prevalence").fetchone()[0]
    gene_count = conn.execute("SELECT COUNT(DISTINCT gene_symbol) FROM prevalence WHERE gene_symbol != ''").fetchone()[
        0
    ]
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('row_count', ?)", (str(count),))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('gene_count', ?)", (str(gene_count),))

    conn.commit()
    conn.close()
    logger.info(f"Orphanet DB built: {gene_count} genes, {count} prevalence entries → {db_path}")
    return db_path


if __name__ == "__main__":
    import argparse

    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument("xml_path", nargs="?", default="data/db/en_product9_prev.xml")
    parser.add_argument("--db-path", default=DEFAULT_DB_PATH)
    args = parser.parse_args()
    build_db(args.xml_path, args.db_path)
