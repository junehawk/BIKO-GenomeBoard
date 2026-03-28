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
    """Build a human-readable prevalence summary for a gene.
    Shows one entry per disease (most specific prevalence class), max 3 diseases.
    """
    entries = get_prevalence_by_gene(gene, db_path)
    if not entries:
        return ""
    # Group by disease, keep best (most specific) prevalence per disease
    best = {}
    for e in entries:
        name = e["disease_name"]
        if name not in best or (e["val_moy"] and (not best[name]["val_moy"] or e["val_moy"] > best[name]["val_moy"])):
            best[name] = e
    parts = []
    for name, e in list(best.items())[:3]:
        cls = e["prevalence_class"]
        if cls:
            parts.append(f"{name} (prevalence: {cls})")
        else:
            parts.append(name)
    return ". ".join(parts) + "." if parts else ""


def get_disease_names(gene: str, db_path: Optional[str] = None) -> List[str]:
    """Get unique disease names associated with a gene from Orphanet."""
    entries = get_prevalence_by_gene(gene, db_path)
    seen = set()
    names = []
    for e in entries:
        if e["disease_name"] not in seen:
            seen.add(e["disease_name"])
            names.append(e["disease_name"])
    return names[:5]
