"""Query OMIM gene mapping SQLite database."""

import sqlite3
import logging
from typing import Dict, Optional
from scripts.common.config import get

logger = logging.getLogger(__name__)


def _get_db_path(db_path: Optional[str] = None) -> str:
    return db_path or get("paths.omim_mapping_db") or "data/db/omim_mapping.sqlite3"


def get_mim_for_gene(gene: str, db_path: Optional[str] = None) -> Optional[Dict]:
    path = _get_db_path(db_path)
    try:
        conn = sqlite3.connect(path)
        row = conn.execute(
            "SELECT mim_number, entry_type, entrez_id, ensembl_id "
            "FROM mim2gene WHERE gene_symbol = ? LIMIT 1",
            (gene,),
        ).fetchone()
        conn.close()
        if not row:
            return None
        return {
            "gene": gene,
            "mim_number": row[0],
            "entry_type": row[1],
            "entrez_id": row[2],
            "ensembl_id": row[3],
            "url": f"https://omim.org/entry/{row[0]}",
        }
    except Exception as e:
        logger.debug(f"OMIM mapping query failed for {gene}: {e}")
        return None
