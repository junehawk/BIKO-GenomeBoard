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

from scripts.annotation.parse_annotation import canonical_consequence
from scripts.classification.in_silico import spliceai_delta_max
from scripts.common.hgvs_utils import extract_protein_position, normalize_hgvsp_for_civic
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
# Consequence canonicalization (raw VEP SO term *or* BIKO short label *or*
# compound "X / Y" → primary SO term) lives in
# scripts/annotation/parse_annotation.py::canonical_consequence — the single
# source of truth shared with the Clinical Board selector and TMB. This module
# previously carried its own divergent copy that split on neither " / " nor
# (in this function) "&", silently mis-canonicalizing compound and
# splice-region labels on production data (v2.7 review CROS-03 / CROS-11).


def _canonicalize_so_term(raw: str) -> str:
    """Return the canonical VEP SO term for a raw term or BIKO label.

    Thin wrapper over :func:`canonical_consequence`; retained for the call
    sites and tests that reference it by name.
    """
    return canonical_consequence(raw)


NULL_CONSEQUENCES = frozenset(
    {
        "stop_gained",
        "frameshift_variant",
        "splice_donor_variant",
        "splice_acceptor_variant",
    }
)

MISSENSE_CONSEQUENCES = frozenset(
    {
        "missense_variant",
    }
)

INFRAME_CONSEQUENCES = frozenset(
    {
        "inframe_insertion",
        "inframe_deletion",
    }
)

SYNONYMOUS_CONSEQUENCES = frozenset(
    {
        "synonymous_variant",
    }
)

# v1 de novo support (PS2/PM6): mirror the selector's B1 consequence gate so
# the classification engine cannot stamp PS2/PM6 on an intronic passenger.
# Rationale — spec Q4 collision point 2: a de novo synonymous variant in a
# constrained gene must not reach LP purely on PS2+PM2+PP3, so the gate runs
# inside the collector.
_DENOVO_PROTEIN_IMPACTING = frozenset(
    {
        "missense_variant",
        "stop_gained",
        "stop_lost",
        "frameshift_variant",
        "splice_donor_variant",
        "splice_acceptor_variant",
        "start_lost",
        "inframe_insertion",
        "inframe_deletion",
    }
)
# Non-coding consequences that can be rescued into the de novo carve-out via
# SpliceAI delta_max >= 0.2 (Tavtigian et al. 2023 PP3-moderate threshold).
_DENOVO_SPLICE_RESCUE = frozenset(
    {
        "synonymous_variant",
        "splice_region_variant",
        "intron_variant",
    }
)
_DENOVO_SPLICEAI_RESCUE_THRESHOLD = 0.2


def _get_consequence(variant: Variant) -> Optional[str]:
    """Return the primary VEP consequence as a canonical SO term, or None.

    Delegates to :func:`canonical_consequence`, which handles raw SO terms,
    BIKO short labels, and compound labels joined by " / " (``format_consequence``)
    or "&" (raw VEP) — returning the first/most-severe component.
    """
    if not variant.consequence:
        return None
    return canonical_consequence(variant.consequence) or None


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
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Gln": "Q",
        "Glu": "E",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
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
    """Extract the max SpliceAI delta score for ``variant``, or None.

    Reads ``variant.in_silico`` — the real-pipeline storage location — via the
    canonical :func:`spliceai_delta_max` helper, which prefers a pre-computed
    ``in_silico['spliceai_max']`` and otherwise derives the max from the four
    raw ``spliceai_pred_ds_*`` keys. Falls back to a legacy ``spliceai_max``
    attribute if present.

    Before v2.7 this read a ``variant.spliceai_max`` attribute that the
    ``Variant`` dataclass never defines, so it always returned None — which made
    the BP7 synonymous-no-splice guard fire on *every* synonymous variant,
    including splice-disrupting ones (a false-benign safety regression).
    """
    in_silico = getattr(variant, "in_silico", None)
    val = spliceai_delta_max(in_silico)
    if val is not None:
        return val
    # Legacy attribute fallback (no production path sets it, kept for safety).
    raw = getattr(variant, "spliceai_max", None)
    if raw is not None:
        try:
            return float(raw)
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
    clinvar_pathogenic_changes: Optional[Dict[int, Set[str]]] = None,
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
        table_hit = _pm1_hotspot_match(variant.gene, extract_protein_position(variant.hgvsp))
        if domain_hit or table_hit == "moderate":
            codes.append("PM1")
        elif table_hit == "supporting":
            codes.append("PM1_Supporting")

    # ── PM4 — Protein length change (in-frame indel) ────────────────────
    if _is_inframe(consequence):
        codes.append("PM4")

    # ── PM5 — Novel missense at position with known pathogenic variant ───
    # PM5 requires a *different* amino-acid change at a residue where another
    # missense is established Pathogenic. When the per-residue pathogenic change
    # set is supplied (clinvar_pathogenic_changes: pos -> {normalised changes}),
    # we fire PM5 only if the residue matches AND the variant's own change is NOT
    # one of the known pathogenic changes — an identical change is PS1, not PM5,
    # and is deliberately not self-asserted here (T1-07). If only the position set
    # is supplied (legacy callers / no AA data) we fall back to the position-only
    # trigger.
    if _is_missense(consequence) and variant.hgvsp:
        protein_pos = _extract_protein_position(variant.hgvsp)
        if protein_pos is not None:
            if clinvar_pathogenic_changes:
                known_changes = clinvar_pathogenic_changes.get(protein_pos)
                if known_changes:
                    variant_change = normalize_hgvsp_for_civic(variant.hgvsp)
                    # Fire only for a parseable, genuinely different change. An
                    # unparseable change (can't confirm it differs) does not fire.
                    if variant_change is not None and variant_change not in known_changes:
                        codes.append("PM5")
            elif clinvar_pathogenic_positions and protein_pos in clinvar_pathogenic_positions:
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


# ── PS2 / PM6 collector (v1 de novo support) ─────────────────────────────────
#
# Spec: _workspace/v23-engineering/00_clinical_denovo_spec.md, Q1 + Q4.
#
# Kept **separate** from :func:`collect_additional_evidence` per the
# implementer checklist — PS2/PM6 come from trio INFO flags on the variant,
# not from gene constraint or VEP annotation, so the two collectors have
# independent inputs and error modes. The gate mirror inside this function
# (protein-impacting + SpliceAI rescue) prevents a de novo intronic
# passenger from collecting PS2 on a purely PS2+PM2+PP3 LP path.


def _denovo_spliceai_rescue(variant: Variant) -> bool:
    """Return True if the variant has SpliceAI delta_max >= 0.2.

    Uses :func:`_get_spliceai_max`, which derives the delta-max from the raw
    ``spliceai_pred_ds_*`` keys when the pre-computed ``spliceai_max`` is
    absent. Before v2.7 this read a never-populated ``in_silico['spliceai_max']``
    key, so the de-novo PS2/PM6 splice rescue was silent dead code on
    production data.
    """
    val = _get_spliceai_max(variant)
    return val is not None and val >= _DENOVO_SPLICEAI_RESCUE_THRESHOLD


def collect_denovo_evidence(variant: Variant) -> List[str]:
    """Collect PS2 / PM6 ACMG codes from trio-aware ``Variant`` flags.

    Rules (spec Q1):
    - ``variant.confirmed_denovo == True`` → ``["PS2_Moderate"]``
      (downgraded from Strong per ClinGen SVI 2018, PMID 30311383, because
      BIKO does not validate phenotype-gene specificity).
    - ``variant.inheritance in {"de_novo", "confirmed_de_novo"}`` (but not
      confirmed) → ``["PM6"]``.
    - No de novo flag → ``[]``.

    Consequence gate (spec Q4 collision point 2) is mirrored here: only fires
    when the variant's consequence is protein-impacting or the variant is
    rescued by SpliceAI delta_max >= 0.2 on a splice-region / synonymous /
    intronic consequence. Otherwise returns ``[]`` even when the de novo
    flag is present — ACMG 2015 PS2/PM6 wording contemplates "disease-causing
    variant", and an intronic passenger without splice signal is not one.
    """
    inheritance = getattr(variant, "inheritance", None)
    confirmed = bool(getattr(variant, "confirmed_denovo", False))
    if inheritance not in ("de_novo", "confirmed_de_novo") and not confirmed:
        return []

    consequence = _get_consequence(variant)
    if not consequence:
        return []

    passes_gate = consequence in _DENOVO_PROTEIN_IMPACTING
    if not passes_gate and consequence in _DENOVO_SPLICE_RESCUE:
        passes_gate = _denovo_spliceai_rescue(variant)

    if not passes_gate:
        return []

    if confirmed:
        return ["PS2_Moderate"]
    return ["PM6"]
