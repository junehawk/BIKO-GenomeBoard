"""Query local HPO SQLite database for gene-phenotype associations."""

import sqlite3
import logging
from typing import Dict, List, Optional
from scripts.common.config import get

logger = logging.getLogger(__name__)


def _get_db_path(db_path: Optional[str] = None) -> str:
    return db_path or get("paths.hpo_db") or "data/db/hpo.sqlite3"


def get_genes_for_hpo(hpo_id: str, db_path: Optional[str] = None) -> List[str]:
    """Get all gene symbols associated with an HPO term."""
    path = _get_db_path(db_path)
    try:
        conn = sqlite3.connect(path)
        rows = conn.execute(
            "SELECT DISTINCT gene_symbol FROM gene_phenotype WHERE hpo_id = ?",
            (hpo_id,),
        ).fetchall()
        conn.close()
        return [r[0] for r in rows]
    except Exception as e:
        logger.warning(f"HPO local DB query failed for {hpo_id}: {e}")
        return []


def get_hpo_for_gene(gene: str, db_path: Optional[str] = None) -> List[Dict]:
    """Get all HPO terms associated with a gene."""
    path = _get_db_path(db_path)
    try:
        conn = sqlite3.connect(path)
        rows = conn.execute(
            "SELECT hpo_id, hpo_name FROM gene_phenotype WHERE gene_symbol = ?",
            (gene,),
        ).fetchall()
        conn.close()
        return [{"hpo_id": r[0], "hpo_name": r[1]} for r in rows]
    except Exception as e:
        logger.warning(f"HPO local DB query failed for {gene}: {e}")
        return []


def resolve_hpo_terms_local(hpo_ids: List[str], db_path: Optional[str] = None) -> List[Dict]:
    """Resolve HPO IDs to names and genes using local DB (no API).
    Same return format as hpo_matcher.resolve_hpo_terms().
    """
    path = _get_db_path(db_path)
    results = []
    try:
        conn = sqlite3.connect(path)
        for hpo_id in hpo_ids:
            hpo_id = hpo_id.strip()
            if not hpo_id.startswith("HP:"):
                continue
            # Get term name
            row = conn.execute(
                "SELECT hpo_name FROM gene_phenotype WHERE hpo_id = ? LIMIT 1",
                (hpo_id,),
            ).fetchone()
            name = row[0] if row else hpo_id
            # Get associated genes
            gene_rows = conn.execute(
                "SELECT DISTINCT gene_symbol FROM gene_phenotype WHERE hpo_id = ?",
                (hpo_id,),
            ).fetchall()
            genes = [r[0] for r in gene_rows]
            results.append({"id": hpo_id, "name": name, "genes": genes})
        conn.close()
    except Exception as e:
        logger.warning(f"HPO local DB resolve failed: {e}")
    return results
