"""Tests for Cancer-mode AI Board agents and models."""
import json
from unittest.mock import MagicMock

from scripts.clinical_board.models import BoardOpinion, CancerBoardOpinion


def _mock_client():
    client = MagicMock()
    client.generate_json.return_value = {
        "findings": [{"finding": "test", "evidence": "test", "confidence": "high"}],
        "recommendations": ["test"],
        "concerns": [],
        "references": [],
        "confidence": "high",
    }
    return client


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


# ---------------------------------------------------------------------------
# Cancer agent properties (Task 5)
# ---------------------------------------------------------------------------


def test_therapeutic_target_analyst_properties():
    from scripts.clinical_board.agents.therapeutic_target import TherapeuticTargetAnalyst

    agent = TherapeuticTargetAnalyst(client=_mock_client())
    assert agent.agent_name == "Therapeutic Target Analyst"
    assert agent.domain == "therapeutic_target"
    assert (
        "druggable" in agent.system_prompt.lower()
        or "drug" in agent.system_prompt.lower()
    )


def test_tumor_genomics_specialist_properties():
    from scripts.clinical_board.agents.tumor_genomics import TumorGenomicsSpecialist

    agent = TumorGenomicsSpecialist(client=_mock_client())
    assert agent.agent_name == "Tumor Genomics Specialist"
    assert agent.domain == "tumor_genomics"


def test_clinical_evidence_analyst_properties():
    from scripts.clinical_board.agents.clinical_evidence import ClinicalEvidenceAnalyst

    agent = ClinicalEvidenceAnalyst(client=_mock_client())
    assert agent.agent_name == "Clinical Evidence Analyst"
    assert agent.domain == "clinical_evidence"


def test_cancer_agent_analyze_returns_opinion():
    from scripts.clinical_board.agents.therapeutic_target import TherapeuticTargetAnalyst

    agent = TherapeuticTargetAnalyst(client=_mock_client())
    opinion = agent.analyze("test briefing")
    assert opinion.agent_name == "Therapeutic Target Analyst"
    assert opinion.confidence == "high"
