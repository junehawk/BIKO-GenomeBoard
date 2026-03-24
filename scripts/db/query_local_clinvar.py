"""Query local ClinVar SQLite database."""
import sqlite3
import logging
from typing import Optional, Dict, List
from pathlib import Path
from scripts.common.config import get
from scripts.common.models import Variant

logger = logging.getLogger(__name__)

_conn = None


def _get_connection() -> Optional[sqlite3.Connection]:
    global _conn
    if _conn is not None:
        return _conn

    db_path = get("paths.clinvar_db", "data/db/clinvar.sqlite3")
    if not Path(db_path).exists():
        logger.warning(f"ClinVar local DB not found: {db_path}")
        return None

    _conn = sqlite3.connect(db_path, check_same_thread=False)
    _conn.row_factory = sqlite3.Row
    return _conn


def get_db_version() -> Dict:
    """Get ClinVar DB metadata (build date, variant count, etc.)"""
    conn = _get_connection()
    if not conn:
        return {"source": "not available", "build_date": "N/A"}

    cursor = conn.execute("SELECT key, value FROM metadata")
    return dict(cursor.fetchall())


def _derive_acmg_codes(significance: str, review_status: str) -> List[str]:
    """Derive ACMG evidence codes from ClinVar local data."""
    codes = []
    sig_lower = significance.lower()
    review_lower = review_status.lower()

    if "pathogenic" in sig_lower and "conflict" not in sig_lower:
        if "expert panel" in review_lower or "practice guideline" in review_lower:
            codes.extend(["PS1", "PP5"])
        elif "multiple submitters" in review_lower:
            codes.append("PS1")
        else:
            codes.append("PP5")

    return codes


def query_local_clinvar(variant: Variant) -> Dict:
    """Query local ClinVar DB for a variant. Returns same structure as API query_clinvar()."""
    conn = _get_connection()

    if conn is None:
        return {
            "agent": "clinical_geneticist",
            "variant": variant.variant_id,
            "gene": variant.gene,
            "clinvar_significance": "Not Found",
            "clinvar_id": None,
            "acmg_codes": [],
            "review_status": None,
            "api_available": False,
        }

    # Strategy 1: Search by rsID
    row = None
    if variant.rsid:
        cursor = conn.execute(
            "SELECT * FROM variants WHERE rsid = ? LIMIT 1",
            (variant.rsid,)
        )
        row = cursor.fetchone()

    # Strategy 2: Search by position
    if not row:
        cursor = conn.execute(
            "SELECT * FROM variants WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ? LIMIT 1",
            (variant.chrom, variant.pos, variant.ref, variant.alt)
        )
        row = cursor.fetchone()

    # Strategy 3: Search by chrom + pos (relaxed)
    if not row:
        cursor = conn.execute(
            "SELECT * FROM variants WHERE chrom = ? AND pos = ? LIMIT 1",
            (variant.chrom, variant.pos)
        )
        row = cursor.fetchone()

    if not row:
        return {
            "agent": "clinical_geneticist",
            "variant": variant.variant_id,
            "gene": variant.gene,
            "clinvar_significance": "Not Found",
            "clinvar_id": None,
            "acmg_codes": [],
            "review_status": None,
            "api_available": True,  # DB is available, just no hit
        }

    significance = row["clinical_significance"] or "Unknown"
    review_status = row["review_status"] or ""
    acmg_codes = _derive_acmg_codes(significance, review_status)

    return {
        "agent": "clinical_geneticist",
        "variant": variant.variant_id,
        "gene": variant.gene or row["gene"],
        "clinvar_significance": significance,
        "clinvar_id": f"VCV{row['variation_id'].zfill(12)}" if row["variation_id"] else None,
        "review_status": review_status,
        "acmg_codes": acmg_codes,
        "phenotypes": row["phenotype_list"] or "",
        "api_available": True,
    }


def close():
    """Close DB connection."""
    global _conn
    if _conn:
        _conn.close()
        _conn = None
