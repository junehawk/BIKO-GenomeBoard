"""Query local gnomAD SQLite database."""

import sqlite3
import logging
from typing import Dict, Optional
from pathlib import Path
from scripts.common.config import get
from scripts.common.models import Variant

logger = logging.getLogger(__name__)

_conn = None


def _get_connection() -> Optional[sqlite3.Connection]:
    global _conn
    if _conn is not None:
        return _conn
    db_path = get("paths.gnomad_db", "data/db/gnomad.sqlite3")
    if not Path(db_path).exists():
        logger.warning(f"gnomAD local DB not found: {db_path}")
        return None
    _conn = sqlite3.connect(db_path, check_same_thread=False)
    _conn.row_factory = sqlite3.Row
    return _conn


def get_db_version() -> Dict:
    """Get gnomAD DB metadata (build date, variant count, etc.)"""
    conn = _get_connection()
    if not conn:
        return {"source": "not available", "build_date": "N/A"}
    cursor = conn.execute("SELECT key, value FROM metadata")
    return dict(cursor.fetchall())


def query_local_gnomad(variant: Variant) -> Dict:
    """Query local gnomAD DB. Returns same structure as API query_gnomad()."""
    conn = _get_connection()

    if conn is None:
        return {"gnomad_all": None, "gnomad_eas": None, "api_available": False}

    # Strategy 1: exact match
    cursor = conn.execute(
        "SELECT * FROM variants WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ? LIMIT 1",
        (variant.chrom, variant.pos, variant.ref, variant.alt),
    )
    row = cursor.fetchone()

    # Strategy 2: rsID fallback
    if not row and variant.rsid:
        cursor = conn.execute("SELECT * FROM variants WHERE rsid = ? LIMIT 1", (variant.rsid,))
        row = cursor.fetchone()

    if not row:
        return {"gnomad_all": None, "gnomad_eas": None, "api_available": True}

    return {
        "gnomad_all": row["af_global"],
        "gnomad_eas": row["af_eas"],
        "gnomad_afr": row["af_afr"],
        "gnomad_amr": row["af_amr"],
        "gnomad_nfe": row["af_nfe"],
        "gnomad_sas": row["af_sas"],
        "api_available": True,
    }


def close():
    """Close DB connection."""
    global _conn
    if _conn:
        _conn.close()
        _conn = None
