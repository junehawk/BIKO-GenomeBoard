#!/usr/bin/env python3
"""Build local CIViC SQLite database from downloaded TSV files."""

import csv
import logging
import re
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

logger = logging.getLogger(__name__)

DEFAULT_DB_PATH = "data/db/civic.sqlite3"


def build_db(civic_dir: str = "data/db/civic", db_path: str = DEFAULT_DB_PATH):
    """Build SQLite DB from CIViC TSV files."""
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)

    # === Genes table ===
    conn.execute("""CREATE TABLE IF NOT EXISTS genes (
        feature_id INTEGER PRIMARY KEY,
        name TEXT NOT NULL,
        description TEXT,
        aliases TEXT,
        entrez_id TEXT
    )""")

    gene_file = Path(civic_dir) / "GeneSummaries.tsv"
    if gene_file.exists():
        with open(gene_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            conn.execute("DELETE FROM genes")
            for row in reader:
                if row.get("feature_type") != "Gene":
                    continue
                conn.execute(
                    "INSERT OR REPLACE INTO genes VALUES (?,?,?,?,?)",
                    (
                        row.get("feature_id"),
                        row.get("name", ""),
                        row.get("description", ""),
                        row.get("feature_aliases", ""),
                        row.get("entrez_id", ""),
                    ),
                )

    # === Variants table ===
    conn.execute("""CREATE TABLE IF NOT EXISTS variants (
        variant_id INTEGER PRIMARY KEY,
        gene TEXT NOT NULL,
        variant_name TEXT NOT NULL,
        variant_types TEXT,
        chrom TEXT,
        start INTEGER,
        stop INTEGER,
        ref TEXT,
        alt TEXT,
        entrez_id TEXT
    )""")

    variant_file = Path(civic_dir) / "VariantSummaries.tsv"
    if variant_file.exists():
        with open(variant_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            conn.execute("DELETE FROM variants")
            for row in reader:
                chrom = row.get("chromosome", "")
                if chrom and not chrom.startswith("chr"):
                    chrom = f"chr{chrom}" if chrom else ""
                conn.execute(
                    "INSERT OR REPLACE INTO variants VALUES (?,?,?,?,?,?,?,?,?,?)",
                    (
                        row.get("variant_id"),
                        row.get("feature_name", ""),
                        row.get("variant", ""),
                        row.get("variant_types", ""),
                        chrom,
                        int(row["start"]) if row.get("start") and row["start"].isdigit() else None,
                        int(row["stop"]) if row.get("stop") and row["stop"].isdigit() else None,
                        row.get("reference_bases", ""),
                        row.get("variant_bases", ""),
                        row.get("entrez_id", ""),
                    ),
                )

    # === Evidence table ===
    conn.execute("""CREATE TABLE IF NOT EXISTS evidence (
        evidence_id INTEGER PRIMARY KEY,
        gene TEXT,
        variant TEXT,
        disease TEXT,
        therapies TEXT,
        therapy_ids TEXT,
        evidence_type TEXT,
        evidence_direction TEXT,
        evidence_level TEXT,
        significance TEXT,
        evidence_statement TEXT,
        citation_id TEXT,
        citation TEXT,
        nct_ids TEXT,
        rating INTEGER,
        variant_origin TEXT,
        molecular_profile TEXT
    )""")
    # Idempotent migration path for pre-v2.2 databases: add therapy_ids if absent.
    existing_evidence_cols = {row[1] for row in conn.execute("PRAGMA table_info(evidence)")}
    if "therapy_ids" not in existing_evidence_cols:
        conn.execute("ALTER TABLE evidence ADD COLUMN therapy_ids TEXT")

    def _extract_therapy_ids(raw_row: dict) -> str:
        """CIViC's ClinicalEvidenceSummaries.tsv exposes therapy identifiers
        under one of several historical column names; try each, comma-joined."""
        for col in ("therapy_ids", "therapy_id", "drug_ids", "drug_id"):
            val = raw_row.get(col)
            if val:
                return val
        return ""

    evidence_file = Path(civic_dir) / "ClinicalEvidenceSummaries.tsv"
    if evidence_file.exists():
        with open(evidence_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            conn.execute("DELETE FROM evidence")
            for row in reader:
                # Extract gene from molecular_profile (e.g., "KRAS G12D" -> "KRAS")
                mp = row.get("molecular_profile", "")
                gene = mp.split()[0] if mp else ""
                variant = " ".join(mp.split()[1:]) if len(mp.split()) > 1 else mp

                conn.execute(
                    "INSERT OR REPLACE INTO evidence "
                    "(evidence_id, gene, variant, disease, therapies, therapy_ids, "
                    "evidence_type, evidence_direction, evidence_level, significance, "
                    "evidence_statement, citation_id, citation, nct_ids, rating, "
                    "variant_origin, molecular_profile) "
                    "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (
                        row.get("evidence_id"),
                        gene,
                        variant,
                        row.get("disease", ""),
                        row.get("therapies", ""),
                        _extract_therapy_ids(row),
                        row.get("evidence_type", ""),
                        row.get("evidence_direction", ""),
                        row.get("evidence_level", ""),
                        row.get("significance", ""),
                        row.get("evidence_statement", ""),
                        row.get("citation_id", ""),
                        row.get("citation", ""),
                        row.get("nct_ids", ""),
                        int(row["rating"]) if row.get("rating") and row["rating"].isdigit() else None,
                        row.get("variant_origin", ""),
                        mp,
                    ),
                )

    # === Hotspots table (extracted from variants with single AA substitution) ===
    conn.execute("""CREATE TABLE IF NOT EXISTS hotspots (
        gene TEXT NOT NULL,
        position INTEGER NOT NULL,
        variants TEXT,
        PRIMARY KEY (gene, position)
    )""")

    conn.execute("DELETE FROM hotspots")
    cursor = conn.execute("SELECT gene, variant_name FROM variants")
    hotspot_map = {}  # (gene, position) -> [variant_names]
    for gene, variant_name in cursor:
        # Match patterns like V600E, G12D, R249M, etc.
        m = re.match(r"^([A-Z])(\d+)([A-Z*]|del|ins|fs)", variant_name)
        if m:
            pos = int(m.group(2))
            key = (gene, pos)
            if key not in hotspot_map:
                hotspot_map[key] = []
            hotspot_map[key].append(variant_name)

    for (gene, pos), variants in hotspot_map.items():
        conn.execute("INSERT INTO hotspots VALUES (?,?,?)", (gene, pos, ",".join(variants)))

    # Indexes
    conn.execute("CREATE INDEX IF NOT EXISTS idx_genes_name ON genes(name)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_variants_gene ON variants(gene)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_evidence_gene ON evidence(gene)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_evidence_variant ON evidence(variant)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_hotspots_gene ON hotspots(gene)")

    # Metadata
    now = datetime.now(timezone.utc).isoformat()
    conn.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (now,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('source', 'CIViC (civicdb.org)')")
    gene_count = conn.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
    variant_count = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    evidence_count = conn.execute("SELECT COUNT(*) FROM evidence").fetchone()[0]
    hotspot_count = conn.execute("SELECT COUNT(*) FROM hotspots").fetchone()[0]
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('gene_count', ?)", (str(gene_count),))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('variant_count', ?)", (str(variant_count),))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('evidence_count', ?)", (str(evidence_count),))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('hotspot_count', ?)", (str(hotspot_count),))

    conn.commit()
    conn.close()
    logger.info(
        f"CIViC DB built: {gene_count} genes, {variant_count} variants, "
        f"{evidence_count} evidence items, {hotspot_count} hotspots → {db_path}"
    )
    return db_path


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_db()
