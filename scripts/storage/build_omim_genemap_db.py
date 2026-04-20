#!/usr/bin/env python3
"""Build OMIM genemap2 gene-phenotype-inheritance SQLite database."""

import logging
import re
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

logger = logging.getLogger(__name__)
DEFAULT_DB_PATH = "data/db/omim_genemap.sqlite3"

# Inheritance pattern normalisation map
INHERITANCE_MAP = {
    "autosomal dominant": "AD",
    "autosomal recessive": "AR",
    "x-linked": "XL",
    "x-linked dominant": "XLD",
    "x-linked recessive": "XLR",
    "mitochondrial": "MT",
    "digenic": "DIG",
    "somatic mutation": "SOM",
}


def _normalise_inheritance(raw: str) -> str:
    """Normalise an inheritance string.

    Handles combined patterns like "Autosomal dominant/Autosomal recessive"
    as well as standard single patterns.
    """
    raw = raw.strip()
    if not raw:
        return ""

    # Handle slash-separated dual inheritance (e.g. "Autosomal dominant/Autosomal recessive")
    if "/" in raw:
        parts = [p.strip() for p in raw.split("/")]
        normalised = []
        for part in parts:
            key = part.lower()
            if key in INHERITANCE_MAP:
                normalised.append(INHERITANCE_MAP[key])
        if normalised:
            return "/".join(normalised)

    key = raw.lower()
    return INHERITANCE_MAP.get(key, raw)


def _parse_phenotypes(phenotype_str: str) -> list[dict]:
    """Parse OMIM Phenotypes column into structured records.

    Format: {Phenotype name}, MIM_NUMBER (mapping_key), Inheritance; ...
    Entries are separated by semicolons.  Phenotype name may be in curly
    braces (susceptibility) or square brackets (non-disease) or bare.
    """
    if not phenotype_str or not phenotype_str.strip():
        return []

    results = []
    # Split on '; {' or '; [' or '; ' followed by a capital letter to
    # separate entries, but the most reliable split is on '; ' when followed
    # by '{' or '[' or a capital letter that starts a new phenotype.
    # The standard delimiter is '; ' between entries.
    entries = re.split(r";\s+(?=[\{\[]|[A-Z])", phenotype_str.strip())

    for entry in entries:
        entry = entry.strip().rstrip(";").strip()
        if not entry:
            continue

        phenotype_name = ""
        phenotype_mim = ""
        inheritance = ""

        # Extract phenotype name (may be in braces, brackets, or bare)
        # Pattern: optional {Name} or [Name] or bare Name, then comma,
        #          MIM (mapping_key), then inheritance
        m = re.match(
            r"[\{\[]?([^\}\]]+?)[\}\]]?\s*,\s*(\d{6})\s*\(\d\)\s*,?\s*(.*)",
            entry,
        )
        if m:
            phenotype_name = m.group(1).strip()
            phenotype_mim = m.group(2).strip()
            inheritance = _normalise_inheritance(m.group(3).strip())
        else:
            # Fallback: treat entire entry as phenotype name
            phenotype_name = entry

        if phenotype_name:
            results.append(
                {
                    "phenotype": phenotype_name,
                    "phenotype_mim": phenotype_mim,
                    "inheritance": inheritance,
                }
            )

    return results


def build_db(txt_path: str, db_path: str = DEFAULT_DB_PATH) -> str:
    """Build SQLite DB from OMIM genemap2.txt."""
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)

    conn.execute("""CREATE TABLE IF NOT EXISTS omim_genemap (
        gene TEXT NOT NULL,
        mim_number TEXT,
        phenotype TEXT,
        inheritance TEXT,
        phenotype_mim TEXT
    )""")
    conn.execute("DELETE FROM omim_genemap")

    rows: list[tuple] = []

    with open(txt_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            # genemap2 columns: 0-Chromosome ... 5-MIM Number ... 8-Approved Gene Symbol ... 12-Phenotypes
            if len(parts) < 13:
                continue
            mim_number = parts[5].strip()
            gene = parts[8].strip()
            phenotype_str = parts[12].strip()

            if not gene:
                continue

            phenotypes = _parse_phenotypes(phenotype_str)
            if not phenotypes:
                # Gene exists but has no annotated phenotypes
                rows.append((gene, mim_number, "", "", ""))
            else:
                for p in phenotypes:
                    rows.append(
                        (
                            gene,
                            mim_number,
                            p["phenotype"],
                            p["inheritance"],
                            p["phenotype_mim"],
                        )
                    )

    conn.executemany(
        "INSERT INTO omim_genemap (gene, mim_number, phenotype, inheritance, phenotype_mim) VALUES (?, ?, ?, ?, ?)",
        rows,
    )

    conn.execute("CREATE INDEX IF NOT EXISTS idx_genemap_gene ON omim_genemap(gene)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_genemap_mim ON omim_genemap(mim_number)")

    now = datetime.now(timezone.utc).isoformat()
    conn.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (now,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('source', ?)", (txt_path,))
    count = conn.execute("SELECT COUNT(*) FROM omim_genemap").fetchone()[0]
    gene_count = conn.execute("SELECT COUNT(DISTINCT gene) FROM omim_genemap").fetchone()[0]
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('record_count', ?)", (str(count),))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('gene_count', ?)", (str(gene_count),))

    conn.commit()
    conn.close()
    logger.info(f"OMIM genemap DB built: {gene_count} genes, {count} records → {db_path}")
    return db_path


if __name__ == "__main__":
    import argparse

    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description="Build OMIM genemap2 SQLite database")
    parser.add_argument(
        "txt_path",
        nargs="?",
        default="data/db/genemap2.txt",
        help="Path to genemap2.txt (default: data/db/genemap2.txt)",
    )
    parser.add_argument("--db-path", default=DEFAULT_DB_PATH, help=f"Output SQLite path (default: {DEFAULT_DB_PATH})")
    args = parser.parse_args()
    build_db(args.txt_path, args.db_path)
