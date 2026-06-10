"""Query local ClinVar SQLite database."""

import logging
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional, Set

from scripts.common.config import get
from scripts.common.hgvs_utils import extract_protein_position, normalize_hgvsp_for_civic
from scripts.common.models import Variant

logger = logging.getLogger(__name__)

_conn = None

# Per-gene memo for get_clinvar_pathogenic_positions. classify_variants() runs
# the same gene many times in a batch (one row per variant) — caching the
# decoded protein-position set per gene avoids re-querying SQLite for every
# variant. Cleared by reset_cache_for_tests() so unit tests can swap DBs.
_PATHOGENIC_POS_CACHE: Dict[str, Set[int]] = {}
# Per-residue pathogenic amino-acid changes (normalised, e.g. {175: {"R175H"}}) per gene,
# for the refined PM5 different-AA gate (T1-07). Same lifecycle as the position cache.
_PATHOGENIC_CHANGES_CACHE: Dict[str, Dict[int, Set[str]]] = {}
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
    """Derive ACMG evidence codes from a ClinVar local-DB hit.

    A generic "ClinVar reports pathogenic" hit is annotated with **PP5 only**
    (scoring-excluded per ClinGen SVI 2018 — see
    ``acmg_engine._count_by_strength``). PS1 is intentionally NOT emitted here.

    PS1 means specifically "same amino-acid change as a previously established
    pathogenic variant" — a claim this gene+position lookup does not verify.
    Emitting PS1 for *any* pathogenic ClinVar hit (a) mislabels the evidence and
    (b) collides with the engine's self-computed PM5 (PS1/PM5 are mutually
    exclusive per ACMG — same vs different AA change), which could manufacture a
    spurious Pathogenic call (ps:1 + pm:n) from a single ClinVar entry. The
    actual classification weight of a high-confidence ClinVar verdict is applied
    independently by ``acmg_engine.apply_clinvar_override`` for exactly the same
    expert-panel / multiple-submitter review tiers, so dropping PS1 here does
    not under-call (v2.7 review CLAS-04). ``review_status`` is retained in the
    signature for callers and a future hgvsp-gated true PS1 path.
    """
    codes = []
    sig_lower = significance.lower()

    if "pathogenic" in sig_lower and "conflict" not in sig_lower:
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
    relaxed = False
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

    # Strategy 3: Search by chrom + pos (relaxed). A hit here is a DIFFERENT allele
    # at the same coordinate (exact ref/alt already missed), so its verdict must not
    # be attributed to this variant — flagged relaxed below.
    if not row:
        cursor = conn.execute(
            "SELECT * FROM variants WHERE chrom = ? AND pos = ? LIMIT 1", (variant.chrom, variant.pos)
        )
        row = cursor.fetchone()
        if row:
            relaxed = True

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

    if relaxed:
        # Coordinate-only match to a different allele. Report it as "Not Found" for
        # classification (so apply_clinvar_override / conflict reconciliation skip it),
        # and surface the neighbouring entry only as informational relaxed_* fields.
        return {
            "agent": "clinical_geneticist",
            "variant": variant.variant_id,
            "gene": variant.gene or row["gene"],
            "clinvar_significance": "Not Found",
            "clinvar_id": None,
            "review_status": None,
            "acmg_codes": [],
            "relaxed_match": True,
            "relaxed_significance": significance,
            "relaxed_review_status": review_status,
            "api_available": True,
        }

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
        "relaxed_match": False,
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


def get_clinvar_pathogenic_changes(gene: str) -> Dict[int, Set[str]]:
    """Return a map of protein position → set of *normalised* pathogenic amino-acid
    changes (e.g. ``{175: {"R175H", "R175C"}}``) for ``gene`` from ClinVar P/LP entries.

    Powers the refined PM5 gate in ``evidence_collector``: PM5 should fire only for a
    *different* change at a residue that already carries a pathogenic missense — an
    identical change is PS1, not PM5 (T1-07). Returns an empty dict when the DB is
    unavailable, the gene has no parseable pathogenic missense, or the hgvsp column is
    absent (legacy build). Memoised per gene at module scope.
    """
    if not gene:
        return {}
    cached = _PATHOGENIC_CHANGES_CACHE.get(gene)
    if cached is not None:
        return cached

    conn = _get_connection()
    if conn is None or not _probe_hgvsp_availability(conn):
        _PATHOGENIC_CHANGES_CACHE[gene] = {}
        return _PATHOGENIC_CHANGES_CACHE[gene]

    changes: Dict[int, Set[str]] = {}
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
            change = normalize_hgvsp_for_civic(hgvsp)
            if pos is not None and change is not None:
                changes.setdefault(pos, set()).add(change)
    except sqlite3.Error as e:
        logger.warning(f"ClinVar HGVSp change lookup failed for {gene}: {e}")

    _PATHOGENIC_CHANGES_CACHE[gene] = changes
    return changes


def reset_cache_for_tests() -> None:
    """Clear the per-gene PM5 caches and the hgvsp-column probe.

    Test helper — production callers should never need this. Tests that
    swap the underlying DB file mid-run rely on it to force the next
    ``get_clinvar_pathogenic_positions`` call to re-query and re-probe.
    """
    global _HGVSP_AVAILABLE
    _PATHOGENIC_POS_CACHE.clear()
    _PATHOGENIC_CHANGES_CACHE.clear()
    _HGVSP_AVAILABLE = None


def close():
    """Close DB connection."""
    global _conn, _HGVSP_AVAILABLE
    if _conn:
        _conn.close()
        _conn = None
    _PATHOGENIC_POS_CACHE.clear()
    _PATHOGENIC_CHANGES_CACHE.clear()
    _HGVSP_AVAILABLE = None
