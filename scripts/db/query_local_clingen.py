"""Query local ClinGen gene-validity SQLite database.

Handles three failure modes gracefully without spamming per-variant logs:

1. `paths.clingen_db` points at a missing file — log once, return None thereafter.
2. The file exists but has no `gene_validity` table (setup_databases.sh
   created an empty shell but the manual CSV export was never provided) —
   log once with the manual-download reminder, return None thereafter.
3. Any other SQLite / OS error — log each occurrence (genuinely unexpected).

The "log once then short-circuit" pattern replaces the previous behaviour
where every variant in a pipeline run produced a fresh
``ClinGen local DB query failed for ...: no such table: gene_validity``
warning. On a 100-variant rare-disease run that was 200+ log lines per
run of pure noise.
"""

import logging
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional

from scripts.common.availability_cache import AvailabilityCache
from scripts.common.config import get

logger = logging.getLogger(__name__)

# Classification strength ranking (highest first)
_RANK = {
    "Definitive": 0,
    "Strong": 1,
    "Moderate": 2,
    "Limited": 3,
    "Disputed": 4,
    "Refuted": 5,
}

# Per-module availability cache — gives this module its own reset
# surface (``reset_availability_cache`` below) so tests rotating
# through multiple DB states don't interfere with the shared cache in
# other modules.
_AVAILABILITY = AvailabilityCache()


def _get_db_path(db_path: Optional[str] = None) -> str:
    return db_path or get("paths.clingen_db") or "data/db/clingen.sqlite3"


def reset_availability_cache() -> None:
    """Clear the memoised availability decision. Tests that build a fresh
    ClinGen DB in tmp_path must call this between runs."""
    _AVAILABILITY.reset()


def _probe(path: str) -> bool:
    """Probe for a usable ClinGen DB at ``path``. Logs a single WARNING
    on any failure mode (missing file, empty shell, SQLite error); the
    :class:`AvailabilityCache` wrapper ensures at most one invocation
    per path per reset cycle.
    """
    if not Path(path).exists():
        logger.warning(
            "ClinGen local DB not found at %s — skipping ClinGen gene-validity "
            "lookups for this run. To enable, export the ClinGen gene-validity "
            "CSV from https://search.clinicalgenome.org/kb/gene-validity and "
            "re-run scripts/setup_databases.sh.",
            path,
        )
        return False

    try:
        conn = sqlite3.connect(path)
        try:
            rows = conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='gene_validity'").fetchone()
            if rows is None:
                logger.warning(
                    "ClinGen local DB at %s has no `gene_validity` table — "
                    "it looks like an empty shell. Run the manual ClinGen CSV "
                    "export and re-run scripts/setup_databases.sh to populate.",
                    path,
                )
                return False
        finally:
            conn.close()
    except sqlite3.Error as e:
        logger.warning("ClinGen local DB probe failed at %s: %s", path, e)
        return False

    return True


def _probe_availability(path: str) -> bool:
    """Return True iff ``path`` exposes a usable ClinGen DB. Memoised
    via :class:`AvailabilityCache`; the probe runs at most once per
    path per reset cycle."""
    return _AVAILABILITY.check(path, _probe)


def get_gene_validity_local(gene: str, db_path: Optional[str] = None) -> Optional[str]:
    """Get strongest ClinGen validity classification for a gene.

    Returns ``None`` if the gene is not in the DB, or if the DB itself is
    unavailable (missing file or empty shell). Logs a single warning per
    run when the DB is unavailable; per-variant failures are silent.
    """
    path = _get_db_path(db_path)
    if not _probe_availability(path):
        return None

    try:
        conn = sqlite3.connect(path)
        rows = conn.execute(
            "SELECT classification FROM gene_validity WHERE gene_symbol = ?",
            (gene,),
        ).fetchall()
        conn.close()
    except sqlite3.Error as e:
        logger.warning(f"ClinGen local DB query failed for {gene}: {e}")
        return None

    if not rows:
        return None
    classifications = [r[0] for r in rows]
    return min(classifications, key=lambda c: _RANK.get(c, 99))


def get_gene_disease_pairs(gene: str, db_path: Optional[str] = None) -> List[Dict]:
    """Get all gene-disease validity pairs. Returns [] when the DB is
    unavailable — see get_gene_validity_local for the log-once contract."""
    path = _get_db_path(db_path)
    if not _probe_availability(path):
        return []

    try:
        conn = sqlite3.connect(path)
        rows = conn.execute(
            "SELECT disease, mondo_id, classification, classification_date "
            "FROM gene_validity WHERE gene_symbol = ? ORDER BY classification",
            (gene,),
        ).fetchall()
        conn.close()
    except sqlite3.Error as e:
        logger.warning(f"ClinGen local DB query failed for {gene}: {e}")
        return []

    return [
        {
            "disease": r[0],
            "mondo_id": r[1],
            "classification": r[2],
            "date": r[3],
        }
        for r in rows
    ]
