"""Query GeneReviews local SQLite database."""

import logging
import sqlite3
from typing import Dict, Optional

from scripts.common.config import get

logger = logging.getLogger(__name__)


def _get_db_path(db_path: Optional[str] = None) -> str:
    return db_path or get("paths.genreviews_db") or "data/db/genreviews.sqlite3"


def get_genreviews_for_gene(gene: str, db_path: Optional[str] = None) -> Optional[Dict]:
    """Get GeneReviews entry for a gene (NBK ID, PMID, title, URL)."""
    path = _get_db_path(db_path)
    try:
        conn = sqlite3.connect(path)
        row = conn.execute(
            "SELECT nbk_id, pmid, title, disease_name FROM genreviews WHERE gene_symbol = ? LIMIT 1",
            (gene,),
        ).fetchone()
        conn.close()
        if not row:
            return None
        return {
            "gene": gene,
            "nbk_id": row[0],
            "pmid": row[1],
            "title": row[2],
            "disease_name": row[3],
            "url": f"https://www.ncbi.nlm.nih.gov/books/{row[0]}/",
        }
    except Exception as e:
        logger.debug(f"GeneReviews query failed for {gene}: {e}")
        return None
