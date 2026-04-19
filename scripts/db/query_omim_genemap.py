"""Query OMIM genemap2 gene-phenotype-inheritance SQLite database."""

import logging
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional

from scripts.common.config import get

logger = logging.getLogger(__name__)


def _get_db_path(db_path: Optional[str] = None) -> str:
    return db_path or get("paths.omim_genemap_db") or "data/db/omim_genemap.sqlite3"


def get_gene_phenotypes(gene: str, db_path: Optional[str] = None) -> Optional[List[Dict]]:
    """Get all phenotype associations for a gene.

    Returns list of dicts with keys: phenotype, mim_number, inheritance,
    phenotype_mim.  Returns None if gene is not found or DB is unavailable.
    """
    path = _get_db_path(db_path)
    if not Path(path).exists():
        logger.debug(f"OMIM genemap DB not found: {path}")
        return None
    try:
        conn = sqlite3.connect(path)
        conn.row_factory = sqlite3.Row
        rows = conn.execute(
            "SELECT gene, mim_number, phenotype, inheritance, phenotype_mim "
            "FROM omim_genemap WHERE gene = ? ORDER BY phenotype",
            (gene,),
        ).fetchall()
        conn.close()
        if not rows:
            return None
        return [
            {
                "phenotype": row["phenotype"],
                "mim_number": row["mim_number"],
                "inheritance": row["inheritance"],
                "phenotype_mim": row["phenotype_mim"],
            }
            for row in rows
        ]
    except Exception as e:
        logger.debug(f"OMIM genemap query failed for {gene}: {e}")
        return None


def get_inheritance_patterns(gene: str, db_path: Optional[str] = None) -> Optional[List[str]]:
    """Get unique inheritance patterns for a gene.

    Returns a sorted list of normalised patterns (e.g. ["AD", "AR"]),
    or None if gene is not found or DB is unavailable.
    """
    path = _get_db_path(db_path)
    if not Path(path).exists():
        logger.debug(f"OMIM genemap DB not found: {path}")
        return None
    try:
        conn = sqlite3.connect(path)
        rows = conn.execute(
            "SELECT DISTINCT inheritance FROM omim_genemap WHERE gene = ? AND inheritance != ''",
            (gene,),
        ).fetchall()
        conn.close()
        if not rows:
            return None
        # Flatten any combined patterns (e.g. "AD/AR") into individual items
        patterns: set[str] = set()
        for row in rows:
            for part in row[0].split("/"):
                part = part.strip()
                if part:
                    patterns.add(part)
        return sorted(patterns) if patterns else None
    except Exception as e:
        logger.debug(f"OMIM genemap inheritance query failed for {gene}: {e}")
        return None
