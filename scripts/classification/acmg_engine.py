# scripts/classification/acmg_engine.py
from __future__ import annotations

import json
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

from scripts.common.config import get
from scripts.common.models import AcmgEvidence


@dataclass
class ClassificationResult:
    classification: str  # Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign
    evidence_codes: List[str] = field(default_factory=list)
    conflict: bool = False
    matched_rule: str = ""
    clinvar_override_reason: str = ""


_RULES_CACHE: dict = {}
_RULES_LOCK = threading.Lock()


def _load_rules() -> dict:
    with _RULES_LOCK:
        if _RULES_CACHE:
            return _RULES_CACHE
    rules_path = get("paths.acmg_rules") or str(Path(__file__).parent.parent.parent / "data" / "acmg_rules.json")
    with open(rules_path) as f:
        data = json.load(f)
    with _RULES_LOCK:
        _RULES_CACHE.update(data)
    return _RULES_CACHE


def _count_by_strength(evidences: List[AcmgEvidence]) -> dict:
    """Count evidence codes by ACMG strength prefix."""
    counts = {"pvs": 0, "ps": 0, "pm": 0, "pp": 0, "ba": 0, "bs": 0, "bp": 0}
    for e in evidences:
        code_upper = e.code.upper()
        # Handle _Supporting suffix FIRST — downgrade to supporting level
        if "_SUPPORTING" in code_upper:
            if code_upper.startswith(("PVS", "PS", "PM")):
                counts["pp"] += 1
            elif code_upper.startswith(("BA", "BS")):
                counts["bp"] += 1
            continue
        # Then normal counting
        if code_upper.startswith("PVS"):
            counts["pvs"] += 1
        elif code_upper.startswith("PS"):
            counts["ps"] += 1
        elif code_upper.startswith("PM"):
            counts["pm"] += 1
        elif code_upper.startswith("PP"):
            counts["pp"] += 1
        elif code_upper.startswith("BA"):
            counts["ba"] += 1
        elif code_upper.startswith("BS"):
            counts["bs"] += 1
        elif code_upper.startswith("BP"):
            counts["bp"] += 1
    return counts


def _matches_rule(counts: dict, rule: dict) -> bool:
    """Check if evidence counts meet or exceed a rule's requirements."""
    for key, required in rule.items():
        if counts.get(key, 0) < required:
            return False
    return True


PGX_GENES = set(get("pgx.genes", ["CYP2D6", "CYP2C19", "CYP2C9", "HLA-B", "HLA-A", "NUDT15", "TPMT", "DPYD"]))
RISK_FACTOR_GENES = set(get("pgx.risk_factor_genes", ["APOE"]))


def check_clinvar_conflict(classification: str, clinvar_sig: str) -> bool:
    """Check if engine classification conflicts with ClinVar significance."""
    if not clinvar_sig:
        return False

    # Normalize to lowercase for comparison
    engine_lower = classification.lower().strip()
    clinvar_lower = clinvar_sig.lower().strip()

    # Skip non-comparable ClinVar values
    skip_values = {"not found", "no classification for the single variant", "not provided", "other"}
    if clinvar_lower in skip_values:
        return False

    # Drug Response / Risk Factor engine outputs don't conflict with matching ClinVar
    if engine_lower in ("drug response", "risk factor") and clinvar_lower in ("drug response", "risk factor"):
        return False
    # Drug Response engine vs non-standard ClinVar (e.g., "conflicting classifications") — flag it
    if engine_lower in ("drug response", "risk factor") and clinvar_lower not in ("drug response", "risk factor"):
        # Only flag if ClinVar has a strong pathogenic/benign call
        if "pathogenic" in clinvar_lower or "benign" in clinvar_lower:
            return True
        return False

    # Standard ACMG tier comparison (case-insensitive)
    rank = {
        "benign": 0,
        "likely benign": 1,
        "vus": 2,
        "uncertain significance": 2,
        "likely pathogenic": 3,
        "pathogenic": 4,
    }
    engine_rank = rank.get(engine_lower, -1)
    clinvar_rank = rank.get(clinvar_lower, -1)

    # Handle compound ClinVar values like "Pathogenic/Likely pathogenic"
    if clinvar_rank == -1 and "/" in clinvar_lower:
        parts = [p.strip() for p in clinvar_lower.split("/")]
        ranks = [rank.get(p, -1) for p in parts if rank.get(p, -1) >= 0]
        if ranks:
            clinvar_rank = max(ranks)

    # If either is unranked, can't compare
    if engine_rank == -1 or clinvar_rank == -1:
        return False

    # Conflict if difference >= 2 steps (LP vs P is only 1 step, not flagged)
    return abs(engine_rank - clinvar_rank) >= 2


# ── A4: narrow ClinVar-conflict reconciliation ───────────────────────────────
#
# ACMG/AMP 2015 (PMID 25741868) does NOT mandate that a laboratory must defer
# to ClinVar when ClinVar itself is in a *conflicting* state. It mandates
# reconciliation with reasoning documented and reviewable. ClinGen SVI PM1
# refinement (Garrett 2021, PMID 33280026) further contemplates laboratory-
# level application of PM1 at moderate strength in peer-reviewed hotspot
# regions even when ClinVar classifications lag. The override below is the
# narrow reconciliation pattern cleared by clinical-advisor on 2026-04-14
# (see _workspace/v22-phaseA/artifacts/00_clinical_review.md §2).
#
# Fires ONLY when all four conditions hold:
#   1. Engine evidence chain independently reaches P or LP (pre-override).
#   2. ClinVar status category is "Conflicting" (prefix/category match).
#   3. PM1 fires from the protein-domain hotspot table (A3).
#   4. PM5 fires from an adjacent-residue ClinVar P/LP entry.
#
# On fire: preserve engine classification, populate clinvar_override_reason
# with a ≤200-char human-readable citation pointing at PM1 (gene + domain +
# PMID) and PM5. Config gate: acmg.allow_engine_override_on_conflict
# (default: true).

_CONFLICTING_PREFIXES = (
    "conflicting classifications of pathogenicity",
    "conflicting interpretations of pathogenicity",
    "conflicting_pathogenicity",
)

# Short gene-hotspot descriptors for the override_reason template. Keep
# strings compact so the full reason fits under the 200-char cap even when
# the {ADJ_VAR} / {VARIANT} expansions add ~30 chars.
_HOTSPOT_DESCRIPTORS = {
    "TP53": {
        175: ("DBD R175 hotspot", "30224644"),
        245: ("DBD L3 loop", "30224644"),
        246: ("DBD L3 loop", "30224644"),
        247: ("DBD L3 loop", "30224644"),
        248: ("DBD L3 loop", "30224644"),
        249: ("DBD L3 loop", "30224644"),
        273: ("DBD R273 DNA-contact", "30224644"),
        282: ("DBD R282 hotspot", "30224644"),
    },
    "KRAS": {
        12: ("G-domain P-loop", "27993330"),
        13: ("G-domain P-loop", "27993330"),
        61: ("G-domain switch II", "27993330"),
    },
    "NRAS": {
        12: ("G-domain switch I/II", "27993330"),
        13: ("G-domain switch I/II", "27993330"),
        61: ("G-domain switch I/II", "27993330"),
    },
    "HRAS": {
        12: ("G-domain switch I/II", "27993330"),
        13: ("G-domain switch I/II", "27993330"),
        61: ("G-domain switch I/II", "27993330"),
    },
    "BRAF": {600: ("kinase activation loop", "27993330")},
    "EGFR": {
        719: ("ATP-binding P-loop", "27993330"),
        790: ("T790M gatekeeper", "27993330"),
        858: ("L858R activation loop", "27993330"),
    },
    "IDH1": {132: ("active site", "27993330")},
    "IDH2": {140: ("active site", "27993330")},
    "PIK3CA": {
        542: ("helical domain", "27993330"),
        545: ("helical domain", "27993330"),
        1047: ("H1047R kinase domain", "27993330"),
    },
}


def _is_conflicting_status(status) -> bool:
    """Prefix/category match for ClinVar "Conflicting" states.

    Accepts any status string whose lowercase form begins with a known
    conflicting-pathogenicity prefix. Tolerant of whitespace and case.
    """
    if not status:
        return False
    norm = str(status).strip().lower()
    return any(norm.startswith(p) for p in _CONFLICTING_PREFIXES)


def _override_enabled() -> bool:
    """Return the config gate; defaults to True for v2.2."""
    return bool(get("acmg.allow_engine_override_on_conflict", True))


def _hotspot_descriptor(gene: str, hgvsp: str):
    """Return ``(domain_label, pmid)`` for the A4 reason string."""
    from scripts.common.hgvs_utils import extract_protein_position

    pos = extract_protein_position(hgvsp)
    if pos is None:
        return ("hotspot", "30224644")
    gene_map = _HOTSPOT_DESCRIPTORS.get(gene or "", {})
    return gene_map.get(pos, ("hotspot domain", "30224644"))


def _short_hgvsp(hgvsp: str) -> str:
    """Return a compact 1-letter representation like ``R249M`` for the reason."""
    from scripts.common.hgvs_utils import AA3TO1, extract_protein_position

    if not hgvsp:
        return "variant"
    import re as _re

    m = _re.search(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", hgvsp)
    if m:
        ref = AA3TO1.get(m.group(1), m.group(1))
        alt = AA3TO1.get(m.group(3), m.group(3))
        return f"{ref}{m.group(2)}{alt}"
    m = _re.search(r"p\.([A-Z])(\d+)([A-Z*])", hgvsp)
    if m:
        return f"{m.group(1)}{m.group(2)}{m.group(3)}"
    pos = extract_protein_position(hgvsp)
    return f"p.{pos}" if pos is not None else "variant"


def apply_hotspot_conflict_reconciliation(
    result: ClassificationResult,
    clinvar_significance: str,
    gene: str,
    hgvsp: str,
    adjacent_variant: str = "",
    adjacent_class: str = "Pathogenic",
) -> ClassificationResult:
    """Apply the A4 narrow ClinVar-conflict override to ``result``.

    Mutates ``result`` in place (and also returns it) — preserves the engine
    classification and populates ``clinvar_override_reason`` when all four
    gate conditions hold. Returns the same object unchanged otherwise.

    Args:
        result: engine ``ClassificationResult`` (pre-override).
        clinvar_significance: raw ClinVar ``clinical_significance`` string.
        gene: HUGO symbol.
        hgvsp: HGVSp (3-letter or 1-letter) for the current variant.
        adjacent_variant: short HGVSp (e.g. ``R249S``) of the adjacent-residue
            ClinVar P/LP entry that fired PM5. Optional; when blank the
            reason string uses a generic "adjacent PM5" fragment.
        adjacent_class: ``"Pathogenic"`` or ``"Likely Pathogenic"``.
    """
    if not _override_enabled():
        return result

    # Condition 1: engine reached P or LP independently
    engine_class = (result.classification or "").strip()
    engine_upper = engine_class.lower()
    if engine_upper not in ("pathogenic", "likely pathogenic"):
        return result

    # Condition 2: ClinVar status is Conflicting
    if not _is_conflicting_status(clinvar_significance):
        return result

    # Conditions 3+4: PM1 and PM5 fire
    codes = {c.upper() for c in (result.evidence_codes or [])}
    has_pm1 = "PM1" in codes  # PM1_Supporting intentionally NOT accepted
    has_pm5 = "PM5" in codes
    if not (has_pm1 and has_pm5):
        return result

    domain_label, pmid = _hotspot_descriptor(gene, hgvsp)
    short_class = "P" if engine_upper == "pathogenic" else "LP"
    short_var = _short_hgvsp(hgvsp)
    adj_fragment = (
        f"PM5 ({adjacent_variant} ClinVar {adjacent_class})"
        if adjacent_variant
        else f"PM5 (adjacent ClinVar {adjacent_class})"
    )
    reason = (
        f"engine {short_class} override: PM1 ({gene} {domain_label}, "
        f"PMID {pmid}) + {adj_fragment}; ClinVar {short_var} shows "
        f"conflicting submitters"
    )
    if len(reason) > 200:
        reason = reason[:200]
    result.clinvar_override_reason = reason
    return result


def apply_clinvar_override(engine_classification: str, clinvar_significance: str, review_status: str) -> str:
    """Override ACMG engine classification when ClinVar has high-confidence evidence.

    Rules:
    - Expert panel / practice guideline: trust ClinVar completely
    - Multiple submitters, no conflicts: trust for Pathogenic/Benign
    - Single submitter: no override (insufficient confidence)
    - Conflicting: no override

    Returns the final classification (may be same as engine_classification).
    """
    if not clinvar_significance or not review_status:
        return engine_classification

    sig_lower = clinvar_significance.lower().strip()
    review_lower = review_status.lower().strip()

    # Skip non-classifiable ClinVar entries
    skip_sigs = {
        "not found",
        "drug response",
        "risk factor",
        "other",
        "not provided",
        "no classification for the single variant",
        "-",
        "",
    }
    if sig_lower in skip_sigs:
        return engine_classification

    # Determine ClinVar confidence level
    is_expert = "expert panel" in review_lower or "practice guideline" in review_lower
    is_multi = "multiple submitters" in review_lower and "conflicting" not in review_lower

    if not is_expert and not is_multi:
        return engine_classification  # Single submitter or conflicting — don't override

    # Map ClinVar significance to classification
    if "pathogenic" in sig_lower and "likely" not in sig_lower and "conflict" not in sig_lower:
        if is_expert:
            return "Pathogenic"
        elif is_multi:
            return "Likely Pathogenic"
    elif "likely pathogenic" in sig_lower:
        return "Likely Pathogenic"
    elif "pathogenic/likely pathogenic" in sig_lower:
        return "Likely Pathogenic"
    elif "benign" in sig_lower and "likely" not in sig_lower and "conflict" not in sig_lower:
        if is_expert:
            return "Benign"
        elif is_multi:
            return "Likely Benign"
    elif "likely benign" in sig_lower:
        if is_expert or is_multi:
            return "Likely Benign"

    return engine_classification


def classify_variant(evidences: List[AcmgEvidence], gene: str = None) -> ClassificationResult:
    """Classify variant by ACMG rules. If gene is a PGx gene, return 'Drug Response' instead."""
    if gene and gene in PGX_GENES:
        return ClassificationResult(
            classification="Drug Response",
            evidence_codes=[e.code for e in evidences],
            conflict=False,
            matched_rule="pgx_bypass",
        )
    if gene and gene in RISK_FACTOR_GENES:
        return ClassificationResult(
            classification="Risk Factor",
            evidence_codes=[e.code for e in evidences],
            matched_rule="risk_factor_bypass",
        )
    if not evidences:
        return ClassificationResult(classification="VUS")

    rules = _load_rules()
    counts = _count_by_strength(evidences)
    code_list = [e.code for e in evidences]

    has_pathogenic = any(e.direction == "pathogenic" for e in evidences)
    has_benign = any(e.direction == "benign" for e in evidences)

    # Check conflict: both pathogenic and benign evidence at ANY level
    if has_pathogenic and has_benign:
        path_match = any(_matches_rule(counts, r) for r in rules["pathogenic"] + rules["likely_pathogenic"])
        benign_match = any(_matches_rule(counts, r) for r in rules["benign"] + rules["likely_benign"])
        if path_match and benign_match:
            return ClassificationResult(
                classification="VUS", evidence_codes=code_list, conflict=True, matched_rule="conflicting_evidence"
            )

    # Check Benign first (BA1 is stand-alone)
    for i, rule in enumerate(rules["benign"]):
        if _matches_rule(counts, rule):
            return ClassificationResult(classification="Benign", evidence_codes=code_list, matched_rule=f"benign_{i}")

    for i, rule in enumerate(rules["likely_benign"]):
        if _matches_rule(counts, rule):
            return ClassificationResult(
                classification="Likely Benign", evidence_codes=code_list, matched_rule=f"likely_benign_{i}"
            )

    # Check Pathogenic
    for i, rule in enumerate(rules["pathogenic"]):
        if _matches_rule(counts, rule):
            return ClassificationResult(
                classification="Pathogenic", evidence_codes=code_list, matched_rule=f"pathogenic_{i}"
            )

    for i, rule in enumerate(rules["likely_pathogenic"]):
        if _matches_rule(counts, rule):
            return ClassificationResult(
                classification="Likely Pathogenic", evidence_codes=code_list, matched_rule=f"likely_pathogenic_{i}"
            )

    return ClassificationResult(classification="VUS", evidence_codes=code_list)
