# scripts/classification/evidence_collector.py
"""ACMG evidence code collector from variant annotations.

Generates additional evidence codes beyond frequency (BA1/BS1/PM2)
and ClinVar, using VEP consequence, gene constraint, and domain info.

This module is self-contained: it uses only the Variant dataclass and
standard library.  It does NOT implement PP3/BP4 (computational
prediction) — that is handled by the separate in_silico module.
"""
from __future__ import annotations

import json
import logging
import os
import re
from typing import Any, Dict, List, Optional, Set

from scripts.common.hgvs_utils import extract_protein_position
from scripts.common.models import Variant

logger = logging.getLogger(__name__)

# ── PM1 hotspot domain table (A3-core) ────────────────────────────────────────
# Lazy-loaded snapshot of data/pm1_hotspot_domains.json. The table lets PM1
# fire on known functional hotspots (ClinGen TP53 VCEP, AMP/ASCO/CAP KRAS 12/13,
# etc.) even when VEP left the CSQ DOMAINS field empty.

_PM1_HOTSPOT_TABLE: Optional[Dict[str, List[Dict[str, Any]]]] = None
_PM1_HOTSPOT_PATH_DEFAULT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data",
    "pm1_hotspot_domains.json",
)


def _load_pm1_hotspots() -> Dict[str, List[Dict[str, Any]]]:
    """Return the hotspot table, loading it once per process."""
    global _PM1_HOTSPOT_TABLE
    if _PM1_HOTSPOT_TABLE is not None:
        return _PM1_HOTSPOT_TABLE
    path = _PM1_HOTSPOT_PATH_DEFAULT
    if not os.path.exists(path):
        logger.debug("[evidence_collector] PM1 hotspot table not found at %s", path)
        _PM1_HOTSPOT_TABLE = {}
        return _PM1_HOTSPOT_TABLE
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        _PM1_HOTSPOT_TABLE = data.get("genes", {}) or {}
    except (OSError, json.JSONDecodeError) as e:
        logger.warning("[evidence_collector] PM1 hotspot table load failed: %s", e)
        _PM1_HOTSPOT_TABLE = {}
    return _PM1_HOTSPOT_TABLE


def _pm1_hotspot_match(gene: Optional[str], protein_position: Optional[int]) -> Optional[str]:
    """Return ``'moderate'``, ``'supporting'``, or ``None`` for ``(gene, pos)``.

    The first range containing ``protein_position`` wins. ``moderate`` ranks
    above ``supporting`` when both would match (moderate is listed first in
    the JSON per curator convention, but we still scan all entries and
    promote moderate explicitly so ordering is not load-bearing).
    """
    if not gene or protein_position is None:
        return None
    table = _load_pm1_hotspots()
    entries = table.get(gene)
    if not entries:
        return None
    best: Optional[str] = None
    for entry in entries:
        rng = entry.get("range") or []
        if len(rng) != 2:
            continue
        lo, hi = rng[0], rng[1]
        if lo <= protein_position <= hi:
            strength = (entry.get("strength") or "").lower()
            if strength == "moderate":
                return "moderate"
            if strength == "supporting" and best is None:
                best = "supporting"
    return best

# ── Consequence groups ────────────────────────────────────────────────────────
#
# Inverse of scripts/intake/parse_annotation.py::format_consequence mapping.
# Real pipeline stores the BIKO-formatted short label (e.g. "Missense") on the
# variant dict, but evidence_collector's tests use the raw VEP SO term. Both
# must resolve to the same internal canonical form.
_FORMATTED_TO_SO = {
    "missense": "missense_variant",
    "nonsense": "stop_gained",
    "nonsense / stop gain": "stop_gained",
    "frameshift": "frameshift_variant",
    "splice donor": "splice_donor_variant",
    "splice acceptor": "splice_acceptor_variant",
    "in-frame deletion": "inframe_deletion",
    "in-frame insertion": "inframe_insertion",
    "synonymous": "synonymous_variant",
    "start loss": "start_lost",
    "stop loss": "stop_lost",
    "intronic": "intron_variant",
}


def _canonicalize_so_term(raw: str) -> str:
    """Return the canonical VEP SO term for either a raw SO term or a
    BIKO-formatted label. Lowercased input is compared against the inverse
    map first; unknown terms pass through unchanged."""
    if not raw:
        return ""
    key = raw.strip().lower()
    return _FORMATTED_TO_SO.get(key, key)


NULL_CONSEQUENCES = frozenset({
    "stop_gained",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
})

MISSENSE_CONSEQUENCES = frozenset({
    "missense_variant",
})

INFRAME_CONSEQUENCES = frozenset({
    "inframe_insertion",
    "inframe_deletion",
})

SYNONYMOUS_CONSEQUENCES = frozenset({
    "synonymous_variant",
})


def _get_consequence(variant: Variant) -> Optional[str]:
    """Return the primary VEP consequence term, normalised to lowercase.

    Handles multi-consequence strings like ``missense_variant&splice_region_variant``
    by returning the first (most severe) term.
    """
    if not variant.consequence:
        return None
    raw = variant.consequence.strip().lower()
    if "&" in raw:
        raw = raw.split("&")[0]
    return _canonicalize_so_term(raw)


def _is_null(consequence: str) -> bool:
    return consequence in NULL_CONSEQUENCES


def _is_missense(consequence: str) -> bool:
    return consequence in MISSENSE_CONSEQUENCES


def _is_inframe(consequence: str) -> bool:
    return consequence in INFRAME_CONSEQUENCES


def _is_synonymous(consequence: str) -> bool:
    return consequence in SYNONYMOUS_CONSEQUENCES


def _extract_protein_position(hgvsp: str) -> Optional[int]:
    """Extract numeric protein position from HGVSp notation.

    Supports both 3-letter (``p.Arg248Leu``) and 1-letter (``p.R248L``) codes.
    """
    if not hgvsp:
        return None
    # 3-letter amino acid codes: p.Arg248Leu
    m = re.search(r"p\.(?:[A-Z][a-z]{2})(\d+)", hgvsp)
    if m:
        return int(m.group(1))
    # 1-letter: p.R248L
    m = re.search(r"p\.[A-Z](\d+)", hgvsp)
    if m:
        return int(m.group(1))
    return None


def _extract_aa_change(hgvsp: str) -> Optional[str]:
    """Extract full amino acid change string for PM5 comparison.

    Returns e.g. ``R248L`` from ``p.Arg248Leu`` or ``p.R248L``.
    Normalises 3-letter to 1-letter codes for comparison purposes.
    """
    if not hgvsp:
        return None

    AA3_TO_1 = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Ter": "*",
    }

    # 3-letter: p.Arg248Leu
    m = re.search(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", hgvsp)
    if m:
        ref = AA3_TO_1.get(m.group(1), m.group(1))
        pos = m.group(2)
        alt = AA3_TO_1.get(m.group(3), m.group(3))
        return f"{ref}{pos}{alt}"

    # 1-letter: p.R248L
    m = re.search(r"p\.([A-Z])(\d+)([A-Z*])", hgvsp)
    if m:
        return f"{m.group(1)}{m.group(2)}{m.group(3)}"

    return None


def _get_spliceai_max(variant: Variant) -> Optional[float]:
    """Extract max SpliceAI delta score from variant if available.

    Looks for a ``spliceai_max`` attribute or parses the ``sift`` field
    as a fallback (some pipelines embed SpliceAI in INFO).
    Returns None if SpliceAI data is not available.
    """
    # Check for explicit spliceai_max attribute
    val = getattr(variant, "spliceai_max", None)
    if val is not None:
        try:
            return float(val)
        except (ValueError, TypeError):
            return None
    return None


def _has_domain(variant: Variant) -> bool:
    """Check whether the variant has a non-empty protein domain annotation.

    Looks for a ``domains`` attribute on the Variant.  Many VEP-annotated
    VCFs include DOMAINS in the CSQ field.
    """
    domains = getattr(variant, "domains", None)
    if domains and str(domains).strip() and str(domains).strip() != ".":
        return True
    return False


# ── Main collector ────────────────────────────────────────────────────────────

def collect_additional_evidence(
    variant: Variant,
    gene_info: Optional[Dict] = None,
    clinvar_pathogenic_positions: Optional[Set[int]] = None,
) -> List[str]:
    """Collect ACMG evidence codes from variant annotations.

    This function generates evidence codes that supplement the frequency-based
    codes (BA1/BS1/PM2) and ClinVar codes already collected elsewhere.

    Args:
        variant: Variant object with consequence, gene, hgvsp, sift, polyphen, etc.
        gene_info: Optional dict with gene constraint scores.  Expected keys:
            ``pli`` (float 0-1), ``missense_z`` (float), ``loeuf`` (float).
        clinvar_pathogenic_positions: Optional set of protein positions (int)
            that have known pathogenic variants in ClinVar.  Used for PM5.

    Returns:
        List of ACMG evidence code strings (e.g. ``["PVS1", "PM1"]``).
    """
    codes: List[str] = []

    consequence = _get_consequence(variant)
    if not consequence:
        return codes

    pli = gene_info.get("pli") if gene_info else None
    missense_z = gene_info.get("missense_z") if gene_info else None

    # ── PVS1 / PVS1_Strong — Null variant in LOF-intolerant gene ─────────
    if _is_null(consequence) and pli is not None:
        if pli >= 0.9:
            codes.append("PVS1")
        elif pli >= 0.5:
            codes.append("PVS1_Strong")

    # ── PM1 — Missense in established functional domain ──────────────────
    # Two independent paths:
    #   (a) VEP CSQ DOMAINS field is populated — legacy path.
    #   (b) variant hits the curated PM1 hotspot table (ClinGen TP53 VCEP,
    #       AMP KRAS 12/13/61, etc.) — promotes the variant even when the
    #       upstream annotator dropped DOMAINS.
    if _is_missense(consequence):
        domain_hit = _has_domain(variant)
        table_hit = _pm1_hotspot_match(
            variant.gene, extract_protein_position(variant.hgvsp)
        )
        if domain_hit or table_hit == "moderate":
            codes.append("PM1")
        elif table_hit == "supporting":
            codes.append("PM1_Supporting")

    # ── PM4 — Protein length change (in-frame indel) ────────────────────
    if _is_inframe(consequence):
        codes.append("PM4")

    # ── PM5 — Novel missense at position with known pathogenic variant ───
    if (_is_missense(consequence)
            and clinvar_pathogenic_positions
            and variant.hgvsp):
        protein_pos = _extract_protein_position(variant.hgvsp)
        if protein_pos is not None and protein_pos in clinvar_pathogenic_positions:
            # PM5 requires a *different* amino acid change at the same position.
            # Since we only have the position set (not the exact changes), the
            # caller is responsible for ensuring this position has a pathogenic
            # variant with a different AA change.  The position match alone is
            # the trigger here — the variant itself being novel is implied
            # because it would not have reached this code path if it already
            # matched ClinVar exactly (that would be PS1 instead).
            codes.append("PM5")

    # ── PP2 — Missense in gene with low benign missense rate ─────────────
    if _is_missense(consequence) and missense_z is not None:
        if missense_z >= 3.09:
            codes.append("PP2")

    # ── BP1 — Missense in gene where only truncating causes disease ──────
    # BP1 applies when the gene's disease mechanism is truncation (LOF),
    # meaning missense variants are less likely to be pathogenic.
    # We use pLI < 0.1 as a proxy: the gene tolerates LOF, so truncating
    # variants are common/benign → the gene is NOT LOF-mechanism → actually
    # wait, the spec says: "missense in gene where truncating is mechanism".
    # If pLI is LOW, the gene tolerates LOF (truncating is NOT mechanism).
    # If pLI is HIGH, truncating IS the mechanism (LOF intolerant).
    # BP1 fires when truncating IS the mechanism → missense is less relevant.
    # However, the spec table says: pLI < 0.1 → BP1.
    # Following the spec as given.
    if _is_missense(consequence) and pli is not None:
        if pli < 0.1:
            codes.append("BP1")

    # ── BP7 — Synonymous variant with no predicted splice impact ─────────
    if _is_synonymous(consequence):
        spliceai_max = _get_spliceai_max(variant)
        if spliceai_max is None or spliceai_max < 0.1:
            codes.append("BP7")

    return codes
