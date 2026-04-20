#!/usr/bin/env python3
"""Build local gnomAD v4.1 gene-constraint SQLite database.

Source: https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/constraint/
        gnomad.v4.1.constraint_metrics.tsv

Schema rationale
----------------
gnomAD v4.1 ships per-transcript constraint metrics. We keep only the
**MANE Select** transcript per gene (one row per HGNC symbol) so the
``constraint_metrics`` table can be joined on a bare gene symbol — this
matches how the rest of BIKO addresses gene-level evidence (PVS1, PP2,
de novo carve-out admission). Genes whose MANE Select transcript is
absent are skipped; such genes are typically read-through ORFs or
single-exon paralogs without ACMG-relevant constraint signal anyway.

The TSV uses dotted column names (``lof.pLI``, ``mis.z_score``, etc.).
``csv.DictReader`` handles them transparently — the dots are just part
of the header strings.

Idempotent: ``DROP TABLE IF EXISTS`` plus ``CREATE TABLE`` ensures
re-running the build overwrites any prior content. Metadata table
records source URL, build_date, record_count, source_size_bytes for
``version_manager.get_all_db_versions``.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

logger = logging.getLogger(__name__)

DEFAULT_TSV_PATH = "data/db/gnomad_constraint.tsv"
DEFAULT_DB_PATH = "data/db/gnomad_constraint.sqlite3"
SOURCE_URL = (
    "https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
)


def _to_float(raw):
    """Return ``float(raw)`` or ``None`` for empty / NaN / unparseable input."""
    if raw is None:
        return None
    s = str(raw).strip()
    if not s or s.upper() in {"NA", "NAN", "NULL", "."}:
        return None
    try:
        f = float(s)
    except (TypeError, ValueError):
        return None
    # gnomAD encodes missing as NaN sometimes — float("nan") != float("nan")
    if f != f:  # NaN check
        return None
    return f


def _to_bool(raw) -> bool:
    if raw is None:
        return False
    s = str(raw).strip().lower()
    return s in {"true", "1", "t", "yes"}


def build_db(tsv_path: str = DEFAULT_TSV_PATH, db_path: str = DEFAULT_DB_PATH) -> dict:
    """Build SQLite DB from the gnomAD v4.1 constraint metrics TSV.

    Returns a dict with ``records`` (rows inserted) and ``db_path``.
    """
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"gnomAD constraint TSV not found at {tsv_path}. Download from {SOURCE_URL}")

    Path(db_path).parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(db_path)
    try:
        conn.execute("DROP TABLE IF EXISTS constraint_metrics")
        conn.execute(
            """
            CREATE TABLE constraint_metrics (
                gene TEXT PRIMARY KEY,
                transcript TEXT,
                pli REAL,
                loeuf REAL,
                mis_z REAL,
                syn_z REAL,
                oe_lof REAL,
                oe_mis REAL,
                mane_select INTEGER,
                build_date TEXT
            )
            """
        )
        conn.execute("CREATE INDEX IF NOT EXISTS idx_gnomad_gene ON constraint_metrics(gene)")
        conn.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)")

        now = datetime.now(timezone.utc).isoformat()
        rows_to_insert: list[tuple] = []
        seen_genes: set[str] = set()
        scanned = 0
        skipped_no_mane = 0
        skipped_no_gene = 0

        with open(tsv_path, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                scanned += 1
                gene = (row.get("gene") or "").strip()
                if not gene:
                    skipped_no_gene += 1
                    continue
                if not _to_bool(row.get("mane_select")):
                    # We only keep the MANE Select transcript per gene to
                    # collapse to one row per HGNC symbol.
                    skipped_no_mane += 1
                    continue
                if gene in seen_genes:
                    # Defensive: a gene should not have two MANE Select
                    # transcripts, but if it does (data bug) keep the first.
                    continue
                seen_genes.add(gene)

                rows_to_insert.append(
                    (
                        gene,
                        (row.get("transcript") or "").strip() or None,
                        _to_float(row.get("lof.pLI")),
                        _to_float(row.get("lof.oe_ci.upper")),  # LOEUF
                        _to_float(row.get("mis.z_score")),
                        _to_float(row.get("syn.z_score")),
                        _to_float(row.get("lof.oe")),
                        _to_float(row.get("mis.oe")),
                        1,
                        now,
                    )
                )

        conn.executemany(
            "INSERT INTO constraint_metrics VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows_to_insert,
        )

        record_count = len(rows_to_insert)
        try:
            source_size = str(os.path.getsize(tsv_path))
        except OSError:
            source_size = "unknown"

        meta_pairs = [
            ("source", SOURCE_URL),
            ("source_file", os.path.basename(tsv_path)),
            ("source_size_bytes", source_size),
            ("build_date", now),
            ("record_count", str(record_count)),
            ("gnomad_version", "v4.1"),
            ("schema_note", "MANE Select transcript per gene only"),
        ]
        conn.executemany(
            "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
            meta_pairs,
        )
        conn.commit()
    finally:
        conn.close()

    logger.info(
        "gnomAD constraint DB built: %d genes (scanned=%d, no_mane=%d, no_gene=%d) → %s",
        record_count,
        scanned,
        skipped_no_mane,
        skipped_no_gene,
        db_path,
    )
    return {"records": record_count, "db_path": db_path}


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description="Build gnomAD v4.1 constraint SQLite DB from TSV.")
    parser.add_argument(
        "tsv_path",
        nargs="?",
        default=DEFAULT_TSV_PATH,
        help=f"Path to gnomAD constraint TSV (default: {DEFAULT_TSV_PATH})",
    )
    parser.add_argument(
        "--db-path",
        default=DEFAULT_DB_PATH,
        help=f"Output SQLite path (default: {DEFAULT_DB_PATH})",
    )
    args = parser.parse_args()
    result = build_db(args.tsv_path, args.db_path)
    logger.info("Done: %d records → %s", result["records"], result["db_path"])


if __name__ == "__main__":
    main()
