"""Query local CIViC database for gene descriptions, variant evidence, and hotspots."""

import logging
import re
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional

from scripts.common.config import get

logger = logging.getLogger(__name__)

_conn = None


def _get_connection():
    global _conn
    if _conn is not None:
        return _conn
    db_path = get("paths.civic_db", "data/db/civic.sqlite3")
    if not Path(db_path).exists():
        return None
    _conn = sqlite3.connect(db_path, check_same_thread=False)
    _conn.row_factory = sqlite3.Row
    return _conn


def reset_civic_connection():
    """Reset the cached connection (useful for testing)."""
    global _conn
    if _conn:
        _conn.close()
        _conn = None


def get_gene_summary(gene: str) -> Optional[Dict]:
    """Get CIViC gene summary (description, aliases)."""
    conn = _get_connection()
    if not conn:
        return None
    row = conn.execute("SELECT * FROM genes WHERE name = ?", (gene,)).fetchone()
    if row:
        return {"gene": row["name"], "description": row["description"], "aliases": row["aliases"]}
    return None


def get_variant_evidence(gene: str, variant_name: str = None) -> List[Dict]:
    """Get CIViC evidence for a gene (optionally filtered by variant).
    Returns list of evidence items sorted by level (A > B > C > D > E).
    """
    conn = _get_connection()
    if not conn:
        return []

    if variant_name:
        cursor = conn.execute(
            "SELECT * FROM evidence WHERE gene = ? AND variant = ? ORDER BY evidence_level, rating DESC",
            (gene, variant_name),
        )
    else:
        cursor = conn.execute("SELECT * FROM evidence WHERE gene = ? ORDER BY evidence_level, rating DESC", (gene,))

    results = []
    for row in cursor:
        results.append(
            {
                "gene": row["gene"],
                "variant": row["variant"],
                "disease": row["disease"],
                "therapies": row["therapies"],
                "evidence_type": row["evidence_type"],
                "evidence_level": row["evidence_level"],
                "significance": row["significance"],
                "statement": row["evidence_statement"],
                "pmid": row["citation_id"],
                "citation": row["citation"],
                "nct_ids": row["nct_ids"],
            }
        )
    return results


def get_treatment_summary(gene: str, variant_name: str = None) -> str:
    """Build a treatment summary string from CIViC evidence for a gene/variant.
    Focuses on Predictive evidence (drug response).
    """
    evidence = get_variant_evidence(gene, variant_name)
    if not evidence:
        evidence = get_variant_evidence(gene)  # Fall back to gene-level

    predictive = [e for e in evidence if e["evidence_type"] == "Predictive" and e["therapies"]]
    if not predictive:
        return ""

    # Group by therapy + significance
    therapies = {}
    for e in predictive:
        key = (e["therapies"], e["significance"])
        if key not in therapies:
            therapies[key] = {"level": e["evidence_level"], "diseases": set(), "pmids": set()}
        therapies[key]["diseases"].add(e["disease"])
        if e["pmid"]:
            therapies[key]["pmids"].add(e["pmid"])

    lines = []
    for (therapy, sig), info in sorted(therapies.items(), key=lambda x: x[1]["level"]):
        diseases = ", ".join(list(info["diseases"])[:3])
        pmid_str = f" (PMID: {', '.join(list(info['pmids'])[:2])})" if info["pmids"] else ""
        lines.append(f"Level {info['level']}: {therapy} — {sig} in {diseases}{pmid_str}")

    return "; ".join(lines[:5])  # Max 5 entries


def is_hotspot(gene: str, protein_position: int) -> bool:
    """Check if a gene + protein position is a known cancer hotspot."""
    conn = _get_connection()
    if not conn:
        return False
    row = conn.execute("SELECT * FROM hotspots WHERE gene = ? AND position = ?", (gene, protein_position)).fetchone()
    return row is not None


def get_hotspot_variants(gene: str, protein_position: int) -> List[str]:
    """Get known hotspot variants at this position."""
    conn = _get_connection()
    if not conn:
        return []
    row = conn.execute(
        "SELECT variants FROM hotspots WHERE gene = ? AND position = ?", (gene, protein_position)
    ).fetchone()
    if row:
        return row["variants"].split(",")
    return []


def extract_protein_position(hgvsp: str) -> Optional[int]:
    """Extract protein position from HGVSp notation.
    e.g., 'p.Arg249Met' -> 249, 'p.Gly12Asp' -> 12, 'p.Val600Glu' -> 600
    Also handles 3-letter and 1-letter codes.
    """
    if not hgvsp:
        return None
    # Match p.Xxx123Yyy (3-letter AA codes)
    m = re.search(r"p\.(?:[A-Z][a-z]{2})(\d+)", hgvsp)
    if m:
        return int(m.group(1))
    # Try 1-letter: p.R249M
    m = re.search(r"p\.[A-Z](\d+)", hgvsp)
    if m:
        return int(m.group(1))
    return None


def get_predictive_evidence_for_tier(gene: str, hgvsp: str, db_path: Optional[str] = None) -> Dict:
    """Get Predictive evidence for tier determination.
    Returns: {"match_level": "variant"|"gene"|"none", "evidence": [...]}
    Only Predictive evidence is returned. match_level distinguishes
    variant-specific from gene-level matches.
    """
    conn = _get_connection()
    if not conn:
        return {"match_level": "none", "evidence": []}

    from scripts.common.hgvs_utils import hgvsp_to_civic_variant

    civic_name = hgvsp_to_civic_variant(hgvsp)

    # Try variant-specific match first
    if civic_name:
        cursor = conn.execute(
            "SELECT * FROM evidence WHERE gene = ? AND variant = ? "
            "AND evidence_type = 'Predictive' ORDER BY evidence_level, rating DESC",
            (gene, civic_name),
        )
        rows = cursor.fetchall()
        if rows:
            evidence = [
                {
                    "gene": r["gene"],
                    "variant": r["variant"],
                    "disease": r["disease"],
                    "therapies": r["therapies"],
                    "evidence_type": r["evidence_type"],
                    "evidence_level": r["evidence_level"],
                    "significance": r["significance"],
                    "statement": r["evidence_statement"],
                    "pmid": r["citation_id"],
                    "citation": r["citation"],
                }
                for r in rows
            ]
            return {"match_level": "variant", "evidence": evidence}

    # Gene-level fallback (for display only, not Tier I elevation)
    cursor = conn.execute(
        "SELECT * FROM evidence WHERE gene = ? AND evidence_type = 'Predictive' ORDER BY evidence_level, rating DESC",
        (gene,),
    )
    rows = cursor.fetchall()
    if rows:
        evidence = [
            {
                "gene": r["gene"],
                "variant": r["variant"],
                "disease": r["disease"],
                "therapies": r["therapies"],
                "evidence_type": r["evidence_type"],
                "evidence_level": r["evidence_level"],
                "significance": r["significance"],
                "statement": r["evidence_statement"],
                "pmid": r["citation_id"],
                "citation": r["citation"],
            }
            for r in rows
        ]
        return {"match_level": "gene", "evidence": evidence}

    return {"match_level": "none", "evidence": []}


def close():
    global _conn
    if _conn:
        _conn.close()
        _conn = None
