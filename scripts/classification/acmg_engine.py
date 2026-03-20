# scripts/classification/acmg_engine.py
from __future__ import annotations
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import List
from scripts.common.models import AcmgEvidence

@dataclass
class ClassificationResult:
    classification: str  # Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign
    evidence_codes: List[str] = field(default_factory=list)
    conflict: bool = False
    matched_rule: str = ""

def _load_rules() -> dict:
    rules_path = Path(__file__).parent.parent.parent / "data" / "acmg_rules.json"
    with open(rules_path) as f:
        return json.load(f)

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

PGX_GENES = {"CYP2D6", "CYP2C19", "CYP2C9", "HLA-B", "HLA-A", "NUDT15", "TPMT", "DPYD"}

def classify_variant(evidences: List[AcmgEvidence], gene: str = None) -> ClassificationResult:
    """Classify variant by ACMG rules. If gene is a PGx gene, return 'Drug Response' instead."""
    if gene and gene in PGX_GENES:
        return ClassificationResult(
            classification="Drug Response",
            evidence_codes=[e.code for e in evidences],
            conflict=False,
            matched_rule="pgx_bypass",
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
        path_match = any(
            _matches_rule(counts, r)
            for r in rules["pathogenic"] + rules["likely_pathogenic"]
        )
        benign_match = any(
            _matches_rule(counts, r)
            for r in rules["benign"] + rules["likely_benign"]
        )
        if path_match and benign_match:
            return ClassificationResult(
                classification="VUS", evidence_codes=code_list, conflict=True,
                matched_rule="conflicting_evidence"
            )

    # Check Benign first (BA1 is stand-alone)
    for i, rule in enumerate(rules["benign"]):
        if _matches_rule(counts, rule):
            return ClassificationResult(
                classification="Benign", evidence_codes=code_list,
                matched_rule=f"benign_{i}"
            )

    for i, rule in enumerate(rules["likely_benign"]):
        if _matches_rule(counts, rule):
            return ClassificationResult(
                classification="Likely Benign", evidence_codes=code_list,
                matched_rule=f"likely_benign_{i}"
            )

    # Check Pathogenic
    for i, rule in enumerate(rules["pathogenic"]):
        if _matches_rule(counts, rule):
            return ClassificationResult(
                classification="Pathogenic", evidence_codes=code_list,
                matched_rule=f"pathogenic_{i}"
            )

    for i, rule in enumerate(rules["likely_pathogenic"]):
        if _matches_rule(counts, rule):
            return ClassificationResult(
                classification="Likely Pathogenic", evidence_codes=code_list,
                matched_rule=f"likely_pathogenic_{i}"
            )

    return ClassificationResult(classification="VUS", evidence_codes=code_list)
