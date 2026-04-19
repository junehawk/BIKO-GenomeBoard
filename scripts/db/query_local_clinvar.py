"""Query local ClinVar SQLite database."""

import logging
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional, Set

from scripts.common.config import get
from scripts.common.hgvs_utils import extract_protein_position
from scripts.common.models import Variant

logger = logging.getLogger(__name__)

_conn = None

# Per-gene memo for get_clinvar_pathogenic_positions. classify_variants() runs
# the same gene many times in a batch (one row per variant) — caching the
# decoded protein-position set per gene avoids re-querying SQLite for every
# variant. Cleared by reset_cache_for_tests() so unit tests can swap DBs.
_PATHOGENIC_POS_CACHE: Dict[str, Set[int]] = {}
# Tri-state availability probe: None = not yet checked, True/False = result.
# Mirrors the lazy-probe pattern used elsewhere so the caller never sees an
# ImportError / FileNotFoundError when the DB is absent.
_HGVSP_AVAILABLE: Optional[bool] = None


def _get_connection() -> Optional[sqlite3.Connection]:
    global _conn
    if _conn is not None:
        return _conn

    db_path = get("paths.clinvar_db", "data/db/clinvar.sqlite3")
    if not Path(db_path).exists():
        logger.warning(f"ClinVar local DB not found: {db_path}")
        return None

    _conn = sqlite3.connect(db_path, check_same_thread=False)
    _conn.row_factory = sqlite3.Row
    return _conn


def get_db_version() -> Dict:
    """Get ClinVar DB metadata (build date, variant count, etc.)"""
    conn = _get_connection()
    if not conn:
        return {"source": "not available", "build_date": "N/A"}

    cursor = conn.execute("SELECT key, value FROM metadata")
    return dict(cursor.fetchall())


def _derive_acmg_codes(significance: str, review_status: str) -> List[str]:
    """Derive ACMG evidence codes from ClinVar local data."""
    codes = []
    sig_lower = significance.lower()
    review_lower = review_status.lower()

    if "pathogenic" in sig_lower and "conflict" not in sig_lower:
        if "expert panel" in review_lower or "practice guideline" in review_lower:
            codes.extend(["PS1", "PP5"])
        elif "multiple submitters" in review_lower:
            codes.append("PS1")
        else:
            codes.append("PP5")

    return codes


def query_local_clinvar(variant: Variant) -> Dict:
    """Query local ClinVar DB for a variant. Returns same structure as API query_clinvar()."""
    conn = _get_connection()

    if conn is None:
        return {
            "agent": "clinical_geneticist",
            "variant": variant.variant_id,
            "gene": variant.gene,
            "clinvar_significance": "Not Found",
            "clinvar_id": None,
            "acmg_codes": [],
            "review_status": None,
            "api_available": False,
        }

    # Strategy 1: Search by rsID
    row = None
    if variant.rsid:
        cursor = conn.execute("SELECT * FROM variants WHERE rsid = ? LIMIT 1", (variant.rsid,))
        row = cursor.fetchone()

    # Strategy 2: Search by position
    if not row:
        cursor = conn.execute(
            "SELECT * FROM variants WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ? LIMIT 1",
            (variant.chrom, variant.pos, variant.ref, variant.alt),
        )
        row = cursor.fetchone()

    # Strategy 3: Search by chrom + pos (relaxed)
    if not row:
        cursor = conn.execute(
            "SELECT * FROM variants WHERE chrom = ? AND pos = ? LIMIT 1", (variant.chrom, variant.pos)
        )
        row = cursor.fetchone()

    if not row:
        return {
            "agent": "clinical_geneticist",
            "variant": variant.variant_id,
            "gene": variant.gene,
            "clinvar_significance": "Not Found",
            "clinvar_id": None,
            "acmg_codes": [],
            "review_status": None,
            "api_available": True,  # DB is available, just no hit
        }

    significance = row["clinical_significance"] or "Unknown"
    review_status = row["review_status"] or ""
    acmg_codes = _derive_acmg_codes(significance, review_status)

    return {
        "agent": "clinical_geneticist",
        "variant": variant.variant_id,
        "gene": variant.gene or row["gene"],
        "clinvar_significance": significance,
        "clinvar_id": f"VCV{row['variation_id'].zfill(12)}" if row["variation_id"] else None,
        "review_status": review_status,
        "acmg_codes": acmg_codes,
        "phenotypes": row["phenotype_list"] or "",
        "api_available": True,
    }


def _probe_hgvsp_availability(conn: sqlite3.Connection) -> bool:
    """Return True if the connected DB has an ``hgvsp`` column.

    Older builds of clinvar.sqlite3 (pre v2.3-T6) lack the column; the
    self-computed PM5 path must short-circuit gracefully on those rather
    than raise. Cached at module level so we only PRAGMA once per process.
    """
    global _HGVSP_AVAILABLE
    if _HGVSP_AVAILABLE is not None:
        return _HGVSP_AVAILABLE
    try:
        cols = {row[1] for row in conn.execute("PRAGMA table_info(variants)")}
    except sqlite3.Error as e:
        logger.warning(f"ClinVar PRAGMA table_info failed: {e}")
        _HGVSP_AVAILABLE = False
        return _HGVSP_AVAILABLE
    _HGVSP_AVAILABLE = "hgvsp" in cols
    if not _HGVSP_AVAILABLE:
        logger.info("ClinVar local DB has no hgvsp column — self-computed PM5 disabled. Rebuild with v2.3-T6 schema.")
    return _HGVSP_AVAILABLE


def get_clinvar_pathogenic_positions(gene: str) -> Set[int]:
    """Return the set of protein positions (int) where ``gene`` has at least
    one ClinVar Pathogenic / Likely Pathogenic entry with a parseable HGVSp.

    Used by ``evidence_collector.collect_additional_evidence`` to fire PM5
    on a novel missense at a residue that already carries a known pathogenic
    variant. Returns an empty set if the DB is unavailable, if ``gene`` has
    no pathogenic entries, or if none of them have an HGVSp that resolves
    to a numeric position.

    Memoises per gene at module scope so the same ``classify_variants()``
    batch does not re-query for every variant.
    """
    if not gene:
        return set()
    cached = _PATHOGENIC_POS_CACHE.get(gene)
    if cached is not None:
        return cached

    conn = _get_connection()
    if conn is None:
        _PATHOGENIC_POS_CACHE[gene] = set()
        return _PATHOGENIC_POS_CACHE[gene]

    if not _probe_hgvsp_availability(conn):
        _PATHOGENIC_POS_CACHE[gene] = set()
        return _PATHOGENIC_POS_CACHE[gene]

    positions: Set[int] = set()
    try:
        cursor = conn.execute(
            """
            SELECT hgvsp FROM variants
            WHERE gene = ?
              AND clinical_significance LIKE '%athogenic%'
              AND clinical_significance NOT LIKE '%conflict%'
              AND hgvsp IS NOT NULL
            """,
            (gene,),
        )
        for row in cursor:
            hgvsp = row["hgvsp"] if isinstance(row, sqlite3.Row) else row[0]
            pos = extract_protein_position(hgvsp)
            if pos is not None:
                positions.add(pos)
    except sqlite3.Error as e:
        logger.warning(f"ClinVar HGVSp lookup failed for {gene}: {e}")

    _PATHOGENIC_POS_CACHE[gene] = positions
    return positions


def reset_cache_for_tests() -> None:
    """Clear the per-gene PM5 cache and the hgvsp-column probe.

    Test helper — production callers should never need this. Tests that
    swap the underlying DB file mid-run rely on it to force the next
    ``get_clinvar_pathogenic_positions`` call to re-query and re-probe.
    """
    global _HGVSP_AVAILABLE
    _PATHOGENIC_POS_CACHE.clear()
    _HGVSP_AVAILABLE = None


def close():
    """Close DB connection."""
    global _conn, _HGVSP_AVAILABLE
    if _conn:
        _conn.close()
        _conn = None
    _PATHOGENIC_POS_CACHE.clear()
    _HGVSP_AVAILABLE = None
