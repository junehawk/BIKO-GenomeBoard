"""Query Orphanet prevalence SQLite database."""

import sqlite3
import logging
from typing import Dict, List, Optional
from scripts.common.config import get

logger = logging.getLogger(__name__)


def _get_db_path(db_path: Optional[str] = None) -> str:
    return db_path or get("paths.orphanet_db") or "data/db/orphanet.sqlite3"


def get_prevalence_by_gene(gene: str, db_path: Optional[str] = None) -> List[Dict]:
    path = _get_db_path(db_path)
    try:
        conn = sqlite3.connect(path)
        rows = conn.execute(
            "SELECT disease_name, prevalence_type, prevalence_class, val_moy, geographic "
            "FROM prevalence WHERE gene_symbol = ? ORDER BY disease_name",
            (gene,),
        ).fetchall()
        conn.close()
        return [
            {"disease_name": r[0], "prevalence_type": r[1],
             "prevalence_class": r[2], "val_moy": r[3], "geographic": r[4]}
            for r in rows
        ]
    except Exception as e:
        logger.debug(f"Orphanet query failed for {gene}: {e}")
        return []


def get_prevalence_text(gene: str, db_path: Optional[str] = None) -> str:
    """Build a human-readable prevalence summary for a gene."""
    entries = get_prevalence_by_gene(gene, db_path)
    if not entries:
        return ""
    parts = []
    seen = set()
    for e in entries:
        key = (e["disease_name"], e["prevalence_class"])
        if key in seen:
            continue
        seen.add(key)
        geo = f" ({e['geographic']})" if e["geographic"] else ""
        parts.append(f"{e['disease_name']}: {e['prevalence_class']}{geo}")
    return "; ".join(parts[:3])
