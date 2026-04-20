"""Query local gnomAD v4.1 gene-constraint SQLite database.

Provides ``get_constraint(gene)`` returning a dict of pLI / LOEUF / missense
Z-score / synonymous Z-score / observed-over-expected ratios, and a
convenience wrapper ``is_constrained(gene)`` that applies the BIKO de novo
carve-out admission threshold (``pLI >= 0.9 AND missense Z >= 3.09``,
Karczewski 2020 / ACMG PVS1+PP2 defaults; PMID 32461654).

Failure modes are handled gracefully without spamming per-variant logs,
following the same "log once then short-circuit" pattern established by
``query_local_clingen.py`` (v2.2 T3 hardening). On a 100-variant rare-disease
run the prior naive approach would emit one warning per query; here the
first probe records a single warning and subsequent calls return silently.

Three failure modes are covered:

1. ``paths.gnomad_constraint_db`` points at a missing file — log once,
   return None / False thereafter.
2. The file exists but lacks a ``constraint_metrics`` table (an empty
   shell from an interrupted ``setup_databases.sh`` run) — log once with
   the build-script reminder, return None / False thereafter.
3. Any other SQLite / OS error — surface a per-call warning (genuinely
   unexpected).
"""

from __future__ import annotations

import logging
import sqlite3
from pathlib import Path
from typing import Dict, Optional

from scripts.common.availability_cache import AvailabilityCache
from scripts.common.config import get

logger = logging.getLogger(__name__)

# Karczewski et al. 2020 (PMID 32461654) — pLI >= 0.9 marks the top decile of
# LoF intolerance; missense Z >= 3.09 marks the top 1% of missense
# intolerance. Both are standard ACMG defaults for PVS1 / PP2 admission.
PLI_THRESHOLD = 0.9
MISSENSE_Z_THRESHOLD = 3.09

# Per-module availability cache — see ``scripts/common/availability_cache.py``.
# Each call-site keeps its own reset surface so tests that rotate through
# multiple DB states do not spill into other modules' memoisation.
_AVAILABILITY = AvailabilityCache()


def _get_db_path(db_path: Optional[str] = None) -> str:
    return db_path or get("paths.gnomad_constraint_db") or "data/db/gnomad_constraint.sqlite3"


def reset_availability_cache() -> None:
    """Clear the memoised availability decision. Tests that build a fresh
    gnomAD constraint DB in tmp_path must call this between runs."""
    _AVAILABILITY.reset()


def _probe(path: str) -> bool:
    """Probe for a usable gnomAD constraint DB at ``path``. Logs a single
    WARNING on any failure mode (missing file, empty shell, SQLite error);
    the :class:`AvailabilityCache` wrapper ensures at most one invocation
    per path per reset cycle.
    """
    if not Path(path).exists():
        logger.warning(
            "gnomAD constraint DB not found at %s — skipping pLI / "
            "missense-Z lookups for this run. To enable, download "
            "gnomad.v4.1.constraint_metrics.tsv from "
            "https://gnomad.broadinstitute.org/downloads#v4-constraint "
            "and run scripts/storage/build_gnomad_constraint_db.py (or "
            "scripts/setup_databases.sh).",
            path,
        )
        return False

    try:
        conn = sqlite3.connect(path)
        try:
            row = conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name='constraint_metrics'"
            ).fetchone()
            if row is None:
                logger.warning(
                    "gnomAD constraint DB at %s has no `constraint_metrics` "
                    "table — it looks like an empty shell. Re-run "
                    "scripts/storage/build_gnomad_constraint_db.py to populate.",
                    path,
                )
                return False
        finally:
            conn.close()
    except sqlite3.Error as e:
        logger.warning("gnomAD constraint DB probe failed at %s: %s", path, e)
        return False

    return True


def _probe_availability(path: str) -> bool:
    """Return True iff ``path`` exposes a usable gnomAD constraint DB.
    Memoised via :class:`AvailabilityCache`; the probe runs at most
    once per path per reset cycle."""
    return _AVAILABILITY.check(path, _probe)


def get_constraint(gene: str, db_path: Optional[str] = None) -> Optional[Dict[str, Optional[float]]]:
    """Return gnomAD v4.1 constraint metrics for a gene.

    Returns dict with keys ``{'pli', 'loeuf', 'missense_z', 'syn_z',
    'oe_lof', 'oe_mis'}`` or None if the gene is not in the DB, or if
    the DB itself is unavailable (missing file or empty shell). Logs a
    single warning per run when the DB is unavailable; per-variant
    failures are silent.
    """
    if not gene:
        return None
    path = _get_db_path(db_path)
    if not _probe_availability(path):
        return None

    try:
        conn = sqlite3.connect(path)
        try:
            row = conn.execute(
                "SELECT pli, loeuf, mis_z, syn_z, oe_lof, oe_mis FROM constraint_metrics WHERE gene = ?",
                (gene,),
            ).fetchone()
        finally:
            conn.close()
    except sqlite3.Error as e:
        logger.warning("gnomAD constraint query failed for %s: %s", gene, e)
        return None

    if row is None:
        return None
    return {
        "pli": row[0],
        "loeuf": row[1],
        "missense_z": row[2],
        "syn_z": row[3],
        "oe_lof": row[4],
        "oe_mis": row[5],
    }


def is_constrained(gene: str, db_path: Optional[str] = None) -> bool:
    """Return True iff ``gene`` has pLI >= 0.9 AND missense Z >= 3.09.

    Convenience wrapper for the variant_selector de novo carve-out
    admission condition. False when the DB is unavailable, the gene is
    missing, or either threshold is not met. Thresholds are the ACMG
    PVS1 / PP2 defaults from Karczewski et al. 2020 (PMID 32461654).
    """
    metrics = get_constraint(gene, db_path)
    if metrics is None:
        return False
    pli = metrics.get("pli")
    mz = metrics.get("missense_z")
    if pli is None or mz is None:
        return False
    try:
        return float(pli) >= PLI_THRESHOLD and float(mz) >= MISSENSE_Z_THRESHOLD
    except (TypeError, ValueError):
        return False


def get_db_version(db_path: Optional[str] = None) -> Dict[str, str]:
    """Read metadata table for version_manager. Returns ``{'source': 'not_available'}``
    when the DB is unavailable (so version_manager can branch cleanly)."""
    path = _get_db_path(db_path)
    if not _probe_availability(path):
        return {"source": "not_available"}

    try:
        conn = sqlite3.connect(path)
        try:
            rows = dict(conn.execute("SELECT key, value FROM metadata").fetchall())
        finally:
            conn.close()
    except sqlite3.Error as e:
        logger.warning("gnomAD constraint metadata read failed at %s: %s", path, e)
        return {"source": "not_available"}

    return {
        "source": "local_db",
        "version": rows.get("gnomad_version", "unknown"),
        "build_date": rows.get("build_date", "unknown"),
        "record_count": rows.get("record_count", "0"),
        "upstream": rows.get("source", ""),
        "path": path,
    }
