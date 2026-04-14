"""AMP/ASCO/CAP 2017 somatic variant tiering engine.

Supports three strategies:
  A: CIViC evidence priority (variant-level first, OncoKB fallback)
  B: OncoKB + CIViC combined (default) — CIViC can elevate, not lower
  C: OncoKB only (backward compatible)

Reference: Li MM et al. J Mol Diagn. 2017;19(1):4-23.
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional
from scripts.clinical.oncokb import get_cancer_gene_info
from scripts.db.query_civic import is_hotspot, extract_protein_position

logger = logging.getLogger(__name__)

AMP_TIER_LABELS = {
    1: "Tier I — Strong Clinical Significance",
    2: "Tier II — Potential Clinical Significance",
    3: "Tier III — Unknown Clinical Significance",
    4: "Tier IV — Benign or Likely Benign",
}


@dataclass
class TierResult:
    tier: int
    tier_label: str
    evidence_source: str
    civic_match_level: str = "none"
    civic_evidence: List[Dict] = field(default_factory=list)


def _make_result(tier: int, source: str, civic_match: str = "none", civic_ev: List[Dict] = None) -> TierResult:
    return TierResult(
        tier=tier,
        tier_label=AMP_TIER_LABELS.get(tier, "Unknown"),
        evidence_source=source,
        civic_match_level=civic_match,
        civic_evidence=civic_ev or [],
    )


def _best_civic_level(evidence: List[Dict]) -> Optional[str]:
    """Get the best (lowest letter) evidence level from a list."""
    levels = [e.get("evidence_level", "Z") for e in evidence]
    if not levels:
        return None
    return min(levels)


def amp_assign_tier(
    classification: str,
    gene: str,
    hgvsp: str = "",
    strategy: str = "B",
    civic_evidence: Optional[Dict] = None,
) -> TierResult:
    """Assign AMP/ASCO/CAP 2017 tier.

    Args:
        classification: ACMG classification (Pathogenic, LP, VUS, LB, Benign, Drug Response, Risk Factor)
        gene: Gene symbol
        hgvsp: HGVSp notation (e.g., p.Val600Glu)
        strategy: "A" (CIViC priority), "B" (combined, default), "C" (OncoKB only)
        civic_evidence: Pre-fetched CIViC evidence dict with match_level and evidence list.
                       If None, CIViC is not consulted (equivalent to strategy C for this call).
    """
    cls_lower = classification.lower()
    gene_info = get_cancer_gene_info(gene) if gene else None
    oncokb_level = gene_info.get("level", "") if gene_info else ""

    if civic_evidence is None:
        civic_evidence = {"match_level": "none", "evidence": []}

    match_level = civic_evidence.get("match_level", "none")
    evidence_items = civic_evidence.get("evidence", [])
    best_level = _best_civic_level(evidence_items)

    # === Drug Response → Tier I (PGx finding, clinically actionable) ===
    if cls_lower == "drug response":
        return _make_result(1, "pharmacogenomic", match_level, evidence_items)

    # === Risk Factor → Tier IV in cancer context (not a cancer biomarker) ===
    if cls_lower == "risk factor":
        return _make_result(4, "risk-factor", match_level, evidence_items)

    # === Always: Benign/Likely Benign → Tier IV (CIViC cannot override) ===
    if "benign" in cls_lower:
        return _make_result(4, "benign", match_level, evidence_items)

    is_pathogenic = "pathogenic" in cls_lower and "benign" not in cls_lower
    is_vus = cls_lower == "vus"

    # === Strategy A: CIViC priority ===
    if strategy == "A":
        if match_level == "variant" and best_level == "A" and is_pathogenic:
            return _make_result(1, "civic-variant-A", match_level, evidence_items)
        if match_level == "variant" and best_level in ("A", "B"):
            tier = 1 if best_level == "A" and is_pathogenic else 2
            return _make_result(tier, f"civic-variant-{best_level}", match_level, evidence_items)
        # Fall through to OncoKB

    # === Strategy B: Combined (CIViC can elevate) ===
    if strategy == "B":
        # CIViC variant-specific Level A + Pathogenic/LP → Tier I
        if match_level == "variant" and best_level == "A" and is_pathogenic:
            return _make_result(1, "civic-variant-A", match_level, evidence_items)
        # CIViC variant-specific Level B → Tier II
        if match_level == "variant" and best_level == "B":
            return _make_result(2, "civic-variant-B", match_level, evidence_items)

    # === Strategy B: CIViC variant-specific Level C-D + Pathogenic/LP → Tier II ===
    if strategy == "B":
        if match_level == "variant" and best_level in ("C", "D") and is_pathogenic:
            return _make_result(2, f"civic-variant-{best_level}", match_level, evidence_items)

    # === OncoKB gene-level (Strategies A fallback, B, C) ===

    # Pathogenic/LP on high-level gene → Tier I
    if is_pathogenic and gene_info and oncokb_level in ("1", "2"):
        return _make_result(1, f"oncokb-gene-L{oncokb_level}", match_level, evidence_items)

    # Pathogenic/LP on any cancer gene → Tier II
    if is_pathogenic and gene_info:
        return _make_result(2, f"oncokb-gene-L{oncokb_level}", match_level, evidence_items)

    # Pathogenic/LP on non-cancer gene → Tier IV (incidental germline finding, not cancer-relevant)
    if is_pathogenic:
        return _make_result(4, "pathogenic-non-cancer", match_level, evidence_items)

    # VUS on cancer gene — check hotspot
    if is_vus and gene_info:
        protein_pos = extract_protein_position(hgvsp)
        if protein_pos and is_hotspot(gene, protein_pos):
            return _make_result(2, "hotspot", match_level, evidence_items)
        return _make_result(3, "oncokb-gene-vus", match_level, evidence_items)

    # Everything else → Tier IV
    return _make_result(4, "default", match_level, evidence_items)
