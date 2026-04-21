"""Query local KOVA (Korean Variant Archive) allele-frequency SQLite database.

KOVA v7 is the unified Korean population allele-frequency resource that
supersedes the older KRGDB TSV, Korea4K TSV and NARD2 TSV inputs. The
upstream release covers ~43.3M variants with ``kova_AC`` / ``kova_AF`` /
``kova_AN`` / ``KOVA_homozygote_count`` pulled from the KOVA v7 cohort.

This module provides :func:`query_kova`, a variant-level lookup that
returns the Korean allele frequency / count / homozygote count for a
single :class:`scripts.common.models.Variant`. The build lives in
``scripts/storage/build_kova_db.py``; the expected DB path is
``data/db/kova.sqlite3`` (override via ``paths.kova_db``).

Graceful-degradation follows the same "log once then short-circuit"
pattern as :mod:`scripts.storage.query_local_clingen` and
:mod:`scripts.storage.query_gnomad_constraint` — the shared
:class:`~scripts.common.availability_cache.AvailabilityCache` guarantees
a single WARNING per run when the DB is missing or an empty shell, so a
~1000-variant whole-genome rare-disease run does not produce 1000 copies
of the same message.
"""

from __future__ import annotations

import logging
import sqlite3
import threading
from pathlib import Path
from typing import Dict, Optional

from scripts.common.availability_cache import AvailabilityCache
from scripts.common.config import get
from scripts.common.models import Variant

logger = logging.getLogger(__name__)

_DEFAULT_DB_PATH = "data/db/kova.sqlite3"

_AVAILABILITY = AvailabilityCache()

# Thread-safe lazy read-only connection cache keyed by DB path.
_CONN_CACHE: Dict[str, sqlite3.Connection] = {}
_CONN_LOCK = threading.Lock()


def _get_db_path(db_path: Optional[str] = None) -> str:
    return db_path or get("paths.kova_db") or _DEFAULT_DB_PATH


def reset_availability_cache() -> None:
    """Clear the memoised availability decision and close cached
    connections. Tests that rotate through multiple DB states in the
    same process must call this between rotations."""
    _AVAILABILITY.reset()
    with _CONN_LOCK:
        for conn in _CONN_CACHE.values():
            try:
                conn.close()
            except sqlite3.Error:
                pass
        _CONN_CACHE.clear()


def _probe(path: str) -> bool:
    """Probe for a usable KOVA DB at ``path``. Logs a single WARNING on
    any failure mode (missing file, empty shell, SQLite error); the
    :class:`AvailabilityCache` wrapper ensures at most one invocation
    per path per reset cycle.
    """
    if not Path(path).exists():
        logger.warning(
            "KOVA local DB not found at %s — skipping Korean allele-frequency "
            "lookups for this run. To enable, run "
            "scripts/storage/build_kova_db.py --input <1_KOVA.v7.tsv.gz> "
            "--output %s.",
            path,
            path,
        )
        return False

    try:
        conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        try:
            row = conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='kova_af'").fetchone()
            if row is None:
                logger.warning(
                    "KOVA local DB at %s has no `kova_af` table — it looks "
                    "like an empty shell. Re-run scripts/storage/build_kova_db.py "
                    "to populate.",
                    path,
                )
                return False
        finally:
            conn.close()
    except sqlite3.Error as e:
        logger.warning("KOVA local DB probe failed at %s: %s", path, e)
        return False

    return True


def _probe_availability(path: str) -> bool:
    return _AVAILABILITY.check(path, _probe)


def _get_conn(path: str) -> sqlite3.Connection:
    """Return a thread-safe read-only connection cached per-path."""
    with _CONN_LOCK:
        conn = _CONN_CACHE.get(path)
        if conn is not None:
            return conn
        conn = sqlite3.connect(
            f"file:{path}?mode=ro",
            uri=True,
            check_same_thread=False,
        )
        _CONN_CACHE[path] = conn
        return conn


def _normalize_chrom(chrom: str) -> str:
    """KOVA stores chromosomes as ``chr1`` .. ``chrX``; accept either
    form from the caller."""
    if chrom is None:
        return ""
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def query_kova(variant: Variant, db_path: Optional[str] = None) -> Optional[Dict[str, object]]:
    """Return KOVA v7 allele-frequency metadata for ``variant``.

    Returns a dict with keys::

        kova_af          (float | None)   — KOVA v7 allele frequency
        kova_ac          (int   | None)   — allele count
        kova_an          (int   | None)   — total alleles
        kova_homozygote  (int   | None)   — homozygote count
        api_available    (bool)           — True iff the DB was queried

    Returns ``None`` when the DB is unavailable (missing file or empty
    shell — a single WARNING is emitted on the first such call per run
    and subsequent calls are silent) or when the variant is not in the
    DB. Callers that need to distinguish "DB down" from "variant absent"
    should use the sibling :func:`is_available` predicate.
    """
    if variant is None:
        return None
    path = _get_db_path(db_path)
    if not _probe_availability(path):
        return None

    chrom = _normalize_chrom(variant.chrom)
    try:
        conn = _get_conn(path)
        row = conn.execute(
            "SELECT kova_af, kova_ac, kova_an, kova_homozygote "
            "FROM kova_af WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?",
            (chrom, variant.pos, variant.ref, variant.alt),
        ).fetchone()
    except sqlite3.Error as e:
        logger.warning(
            "KOVA local DB query failed for %s:%d %s>%s: %s",
            chrom,
            variant.pos,
            variant.ref,
            variant.alt,
            e,
        )
        return None

    if row is None:
        return None
    return {
        "kova_af": row[0],
        "kova_ac": row[1],
        "kova_an": row[2],
        "kova_homozygote": row[3],
        "api_available": True,
    }


def is_available(db_path: Optional[str] = None) -> bool:
    """Return True iff a usable KOVA DB is reachable at ``db_path``
    (or ``paths.kova_db`` when ``db_path`` is None). Callers that need
    to branch on "DB exists" without firing a query should prefer this
    over :func:`query_kova`."""
    return _probe_availability(_get_db_path(db_path))


def get_db_version(db_path: Optional[str] = None) -> Dict[str, str]:
    """Read the KOVA metadata table for ``version_manager``.

    Returns ``{'source': 'not_available'}`` when the DB is unavailable
    so the caller can branch cleanly."""
    path = _get_db_path(db_path)
    if not _probe_availability(path):
        return {"source": "not_available"}

    try:
        conn = _get_conn(path)
        rows = dict(conn.execute("SELECT key, value FROM metadata").fetchall())
    except sqlite3.Error as e:
        logger.warning("KOVA metadata read failed at %s: %s", path, e)
        return {"source": "not_available"}

    return {
        "source": "local_db",
        "version": rows.get("kova_version", "unknown"),
        "build_date": rows.get("build_date", "unknown"),
        "record_count": rows.get("record_count", "0"),
        "upstream": rows.get("source", ""),
        "assembly": rows.get("assembly", "GRCh38"),
        "path": path,
    }


if __name__ == "__main__":
    import json
    import sys

    if len(sys.argv) < 2:
        print(json.dumps({"error": "Usage: python -m scripts.population.query_kova 'chr17:7675088 C>A' [db_path]"}))
        sys.exit(1)
    v = Variant.from_string(sys.argv[1])
    path_arg = sys.argv[2] if len(sys.argv) > 2 else None
    result = query_kova(v, path_arg)
    print(json.dumps({"variant": v.variant_id, "kova": result}, indent=2, default=str))
