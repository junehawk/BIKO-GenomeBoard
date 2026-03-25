"""Query local ClinGen gene-validity SQLite database."""

import sqlite3
import logging
from typing import Dict, List, Optional
from scripts.common.config import get

logger = logging.getLogger(__name__)

# Classification strength ranking (highest first)
_RANK = {
    "Definitive": 0, "Strong": 1, "Moderate": 2,
    "Limited": 3, "Disputed": 4, "Refuted": 5,
}


def _get_db_path(db_path: Optional[str] = None) -> str:
    return db_path or get("paths.clingen_db") or "data/db/clingen.sqlite3"


def get_gene_validity_local(
    gene: str, db_path: Optional[str] = None
) -> Optional[str]:
    """Get strongest ClinGen validity classification for a gene."""
    path = _get_db_path(db_path)
    try:
        conn = sqlite3.connect(path)
        rows = conn.execute(
            "SELECT classification FROM gene_validity WHERE gene_symbol = ?",
            (gene,),
        ).fetchall()
        conn.close()
        if not rows:
            return None
        # Return strongest classification
        classifications = [r[0] for r in rows]
        return min(classifications, key=lambda c: _RANK.get(c, 99))
    except Exception as e:
        logger.warning(f"ClinGen local DB query failed for {gene}: {e}")
        return None


def get_gene_disease_pairs(
    gene: str, db_path: Optional[str] = None
) -> List[Dict]:
    """Get all gene-disease validity pairs."""
    path = _get_db_path(db_path)
    try:
        conn = sqlite3.connect(path)
        rows = conn.execute(
            "SELECT disease, mondo_id, classification, classification_date "
            "FROM gene_validity WHERE gene_symbol = ? ORDER BY classification",
            (gene,),
        ).fetchall()
        conn.close()
        return [
            {
                "disease": r[0],
                "mondo_id": r[1],
                "classification": r[2],
                "date": r[3],
            }
            for r in rows
        ]
    except Exception as e:
        logger.warning(f"ClinGen local DB query failed for {gene}: {e}")
        return []
