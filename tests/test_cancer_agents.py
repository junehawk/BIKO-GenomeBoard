"""Tests for Cancer-mode AI Board agents and models."""
import json

from scripts.clinical_board.models import BoardOpinion, CancerBoardOpinion


def test_cancer_board_opinion_fields():
    """CancerBoardOpinion has treatment-focused fields."""
    op = CancerBoardOpinion(
        therapeutic_implications="EGFR L858R — TKI sensitive",
        therapeutic_evidence="CIViC Level A",
        treatment_options=[{"drug": "Erlotinib", "evidence_level": "A", "resistance_notes": "T790M"}],
        actionable_findings=["EGFR activating mutation"],
        clinical_actions=["Start erlotinib 150mg daily"],
        immunotherapy_eligibility="TMB-low (3.2 mut/Mb) — not eligible",
    )
    assert op.therapeutic_implications == "EGFR L858R — TKI sensitive"
    assert len(op.treatment_options) == 1
    assert op.treatment_options[0]["resistance_notes"] == "T790M"
    assert op.immunotherapy_eligibility.startswith("TMB-low")


def test_cancer_board_opinion_serializable():
    """CancerBoardOpinion can be JSON serialized (for KB storage)."""
    op = CancerBoardOpinion(
        therapeutic_implications="test",
        therapeutic_evidence="test",
        treatment_options=[],
        actionable_findings=[],
        clinical_actions=[],
        immunotherapy_eligibility="",
    )
    data = json.loads(json.dumps(op.__dict__, default=str))
    assert data["therapeutic_implications"] == "test"


def test_unified_disclaimer_english():
    """Both BoardOpinion and CancerBoardOpinion use the same fixed disclaimer."""
    board = BoardOpinion()
    cancer = CancerBoardOpinion(
        therapeutic_implications="",
        therapeutic_evidence="",
        treatment_options=[],
        actionable_findings=[],
        clinical_actions=[],
        immunotherapy_eligibility="",
    )
    assert "AI-Generated" in board.disclaimer
    assert "Google MedGemma" in board.disclaimer
    assert board.disclaimer == cancer.disclaimer
