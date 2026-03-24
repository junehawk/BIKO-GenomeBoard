# scripts/classification/acmg_engine.py
from __future__ import annotations
import json
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import List
from scripts.common.models import AcmgEvidence
from scripts.common.config import get


@dataclass
class ClassificationResult:
    classification: str  # Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign
    evidence_codes: List[str] = field(default_factory=list)
    conflict: bool = False
    matched_rule: str = ""


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
    skip_sigs = {"not found", "drug response", "risk factor", "other",
                 "not provided", "no classification for the single variant", "-", ""}
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
