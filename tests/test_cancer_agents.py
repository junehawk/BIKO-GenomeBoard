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
    assert "druggable" in agent.system_prompt.lower() or "drug" in agent.system_prompt.lower()


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


# ---------------------------------------------------------------------------
# Board Chair mode branching (Task 6)
# ---------------------------------------------------------------------------


def test_board_chair_cancer_mode():
    """Board Chair in cancer mode returns CancerBoardOpinion."""
    from scripts.clinical_board.agents.board_chair import BoardChair

    mock_client = MagicMock()
    mock_client.generate_json.return_value = {
        "therapeutic_implications": "EGFR TKI sensitive",
        "therapeutic_evidence": "Level A evidence",
        "treatment_options": [{"drug": "Erlotinib", "evidence_level": "A"}],
        "actionable_findings": ["EGFR L858R detected"],
        "clinical_actions": ["Start TKI therapy"],
        "immunotherapy_eligibility": "TMB-low",
        "confidence": "high",
        "agent_consensus": "unanimous",
        "dissenting_opinions": [],
        "monitoring_plan": ["CT q3m"],
    }
    chair = BoardChair(client=mock_client)
    result = chair.synthesize("briefing", [], mode="cancer")
    assert isinstance(result, CancerBoardOpinion)
    assert result.therapeutic_implications == "EGFR TKI sensitive"
    assert result.treatment_options[0]["drug"] == "Erlotinib"
    assert result.monitoring_plan == ["CT q3m"]
    assert "AI-Generated" in result.disclaimer


def test_render_cancer_board_opinion():
    """CancerBoardOpinion renders treatment-focused HTML."""
    from scripts.clinical_board.models import AgentOpinion
    from scripts.clinical_board.render import render_board_opinion_html

    opinion = CancerBoardOpinion(
        therapeutic_implications="EGFR L858R — TKI sensitive",
        therapeutic_evidence="CIViC Level A",
        treatment_options=[{"drug": "Erlotinib", "evidence_level": "A", "resistance_notes": "T790M risk"}],
        actionable_findings=["Activating EGFR mutation"],
        clinical_actions=["Start erlotinib"],
        immunotherapy_eligibility="TMB-low",
        agent_opinions=[AgentOpinion(agent_name="Test", domain="test", confidence="high")],
        agent_consensus="unanimous",
    )
    html = render_board_opinion_html(opinion, language="en")
    assert "Therapeutic" in html or "Treatment" in html
    assert "Erlotinib" in html
    assert "T790M" in html
    assert "AI-Generated" in html


def test_render_board_opinion_backward_compatible():
    """Existing BoardOpinion still renders correctly."""
    from scripts.clinical_board.render import render_board_opinion_html

    opinion = BoardOpinion(primary_diagnosis="Li-Fraumeni")
    html = render_board_opinion_html(opinion, language="en")
    assert "Li-Fraumeni" in html


def test_render_includes_selection_metadata_cancer():
    """CancerBoardOpinion with selection_metadata renders the pre-analytic filter caption."""
    from scripts.clinical_board.render import render_board_opinion_html

    opinion = CancerBoardOpinion(
        therapeutic_implications="EGFR L858R — TKI sensitive",
        therapeutic_evidence="CIViC Level A",
        treatment_options=[{"drug": "Erlotinib", "evidence_level": "A", "resistance_notes": ""}],
        selection_metadata={
            "mode": "cancer",
            "total_input": 42,
            "selected": 7,
            "must_included": 3,
            "may_included": 4,
            "excluded": 35,
            "truncated": False,
            "n_dropped": 0,
            "hard_cap_applied": False,
            "empty": False,
            "empty_reason": "",
            "criteria_summary": "Tier I/II + OncoKB 1-2 + ClinVar P/LP",
            "by_selection_reason": {},
        },
    )
    html = render_board_opinion_html(opinion, language="en")
    assert "Pre-analytic filtering" in html
    assert "42" in html
    assert "7" in html
    assert "Tier I/II + OncoKB 1-2 + ClinVar P/LP" in html


def test_render_includes_selection_metadata_rare_disease():
    """BoardOpinion with selection_metadata renders the pre-analytic filter caption,
    including the truncation note when truncated."""
    from scripts.clinical_board.render import render_board_opinion_html

    opinion = BoardOpinion(
        primary_diagnosis="Noonan syndrome",
        selection_metadata={
            "mode": "rare-disease",
            "total_input": 128,
            "selected": 15,
            "must_included": 5,
            "may_included": 10,
            "excluded": 110,
            "truncated": True,
            "n_dropped": 3,
            "hard_cap_applied": False,
            "empty": False,
            "empty_reason": "",
            "criteria_summary": "P/LP + HPO-match + de novo",
            "by_selection_reason": {},
        },
    )
    html = render_board_opinion_html(opinion, language="en")
    assert "Pre-analytic filtering" in html
    assert "128" in html
    assert "15" in html
    assert "P/LP + HPO-match + de novo" in html
    assert "truncated" in html.lower()


def test_agent_panel_placeholder_when_findings_empty():
    """AgentOpinion panels with empty findings render a grey italic placeholder."""
    from scripts.clinical_board.models import AgentOpinion
    from scripts.clinical_board.render import render_board_opinion_html

    empty_agent = AgentOpinion(
        agent_name="Variant Pathologist",
        domain="variant_pathology",
        findings=[],
        confidence="moderate",
    )

    cancer = CancerBoardOpinion(
        therapeutic_implications="Test",
        agent_opinions=[empty_agent],
    )
    cancer_html = render_board_opinion_html(cancer, language="en")
    assert "No specific findings identified for this case." in cancer_html

    rare = BoardOpinion(
        primary_diagnosis="Test",
        agent_opinions=[empty_agent],
    )
    rare_html = render_board_opinion_html(rare, language="en")
    assert "No specific findings identified for this case." in rare_html

    # Korean variant
    rare_ko = render_board_opinion_html(rare, language="ko")
    assert "이 케이스에서 특별한 소견은 확인되지 않았습니다." in rare_ko


def test_agent_panel_normal_when_findings_present():
    """AgentOpinion panels with findings do NOT emit the placeholder message."""
    from scripts.clinical_board.models import AgentOpinion
    from scripts.clinical_board.render import render_board_opinion_html

    agent = AgentOpinion(
        agent_name="Variant Pathologist",
        domain="variant_pathology",
        findings=[{"finding": "BRCA2 truncating"}],
        confidence="high",
    )
    cancer = CancerBoardOpinion(
        therapeutic_implications="Test",
        agent_opinions=[agent],
    )
    html = render_board_opinion_html(cancer, language="en")
    assert "BRCA2 truncating" in html
    assert "No specific findings identified for this case." not in html


def test_render_omits_metadata_when_none():
    """When selection_metadata is None, the caption block is silently skipped."""
    from scripts.clinical_board.render import render_board_opinion_html

    cancer = CancerBoardOpinion(therapeutic_implications="Test")
    rare = BoardOpinion(primary_diagnosis="Test")
    assert cancer.selection_metadata is None
    assert rare.selection_metadata is None
    cancer_html = render_board_opinion_html(cancer, language="en")
    rare_html = render_board_opinion_html(rare, language="en")
    assert "Pre-analytic filtering" not in cancer_html
    assert "Pre-analytic filtering" not in rare_html


def test_orchestrate_format_board_summary_cancer():
    """orchestrate._format_board_summary handles CancerBoardOpinion without AttributeError."""
    from scripts.orchestrate import _format_board_summary

    cancer = CancerBoardOpinion(
        therapeutic_implications="EGFR L858R — TKI sensitive",
        therapeutic_evidence="CIViC Level A",
        treatment_options=[],
        actionable_findings=[],
        clinical_actions=[],
        immunotherapy_eligibility="",
        confidence="moderate",
    )
    line = _format_board_summary(cancer)
    assert "Therapeutic implications" in line
    assert "EGFR L858R" in line
    assert "moderate" in line


def test_orchestrate_format_board_summary_rare_disease():
    """orchestrate._format_board_summary handles BoardOpinion (rare-disease)."""
    from scripts.orchestrate import _format_board_summary

    rd = BoardOpinion(primary_diagnosis="Li-Fraumeni", confidence="high")
    line = _format_board_summary(rd)
    assert "Primary diagnosis" in line
    assert "Li-Fraumeni" in line
    assert "high" in line


def test_cancer_board_opinion_has_optional_headline():
    """CancerBoardOpinion accepts therapeutic_headline; default is empty string."""
    default = CancerBoardOpinion()
    assert default.therapeutic_headline == ""

    populated = CancerBoardOpinion(therapeutic_headline="Stage IV PDAC — KRAS G12D driver")
    assert populated.therapeutic_headline == "Stage IV PDAC — KRAS G12D driver"


def test_board_chair_cancer_emits_headline():
    """Board Chair populates therapeutic_headline from the LLM response."""
    from scripts.clinical_board.agents.board_chair import BoardChair

    mock_client = MagicMock()
    mock_client.generate_json.return_value = {
        "therapeutic_headline": "Stage IV PDAC — KRAS driver",
        "therapeutic_implications": "No approved targeted therapy for KRAS G12D; consider trials.",
        "therapeutic_evidence": "CIViC Level C",
        "treatment_options": [],
        "actionable_findings": [],
        "clinical_actions": [],
        "immunotherapy_eligibility": "",
        "confidence": "moderate",
        "agent_consensus": "majority",
        "dissenting_opinions": [],
        "monitoring_plan": [],
    }
    chair = BoardChair(client=mock_client)
    result = chair.synthesize("briefing", [], mode="cancer")
    assert isinstance(result, CancerBoardOpinion)
    assert result.therapeutic_headline == "Stage IV PDAC — KRAS driver"
    assert result.therapeutic_implications.startswith("No approved")


def test_board_chair_cancer_backward_compat_no_headline():
    """When the LLM omits therapeutic_headline, synthesize still works with empty headline."""
    from scripts.clinical_board.agents.board_chair import BoardChair

    mock_client = MagicMock()
    mock_client.generate_json.return_value = {
        "therapeutic_implications": "Legacy paragraph only.",
        "therapeutic_evidence": "Level B",
        "treatment_options": [],
        "actionable_findings": [],
        "clinical_actions": [],
        "immunotherapy_eligibility": "",
        "confidence": "moderate",
        "agent_consensus": "unanimous",
        "dissenting_opinions": [],
        "monitoring_plan": [],
    }
    chair = BoardChair(client=mock_client)
    result = chair.synthesize("briefing", [], mode="cancer")
    assert isinstance(result, CancerBoardOpinion)
    assert result.therapeutic_headline == ""
    assert result.therapeutic_implications == "Legacy paragraph only."


def test_board_chair_rare_disease_mode_default():
    """Board Chair defaults to rare-disease mode (backward compatible)."""
    from scripts.clinical_board.agents.board_chair import BoardChair

    mock_client = MagicMock()
    mock_client.generate_json.return_value = {
        "primary_diagnosis": "Li-Fraumeni",
        "confidence": "high",
        "agent_consensus": "unanimous",
        "key_findings": [],
        "recommendations": [],
        "differential_diagnoses": [],
        "dissenting_opinions": [],
        "follow_up": [],
    }
    chair = BoardChair(client=mock_client)
    result = chair.synthesize("briefing", [])
    assert isinstance(result, BoardOpinion)
    assert result.primary_diagnosis == "Li-Fraumeni"


def test_render_cancer_headline_and_body():
    """When therapeutic_headline is present, it renders as a bold title and
    therapeutic_implications renders as a separate body paragraph."""
    from scripts.clinical_board.render import render_board_opinion_html

    headline = "Stage IV PDAC — KRAS G12D driver, no standard targeted therapy"
    body = (
        "환자의 종양은 KRAS G12D 변이를 보이며 표준 표적치료제는 부재합니다. "
        "임상시험 등록을 적극 검토할 수 있으며, 면역치료 적응증 평가가 권고됩니다."
    )
    opinion = CancerBoardOpinion(
        therapeutic_headline=headline,
        therapeutic_implications=body,
        therapeutic_evidence="CIViC Level C",
    )
    html = render_board_opinion_html(opinion, language="en")

    # Headline and body must appear as distinct elements.
    assert headline in html
    assert body in html
    headline_idx = html.find(headline)
    body_idx = html.find(body)
    assert headline_idx != -1 and body_idx != -1
    # The body paragraph must come after the headline, not be the same element.
    assert body_idx > headline_idx
    between = html[headline_idx + len(headline) : body_idx]
    assert "</div>" in between, "headline and body should live in separate <div>s"

    # Body must be styled as a readable paragraph, not a 15px bold title.
    body_start = html.rfind("<", 0, body_idx)
    body_tag = html[body_start:body_idx]
    assert "font-size:15px" not in body_tag
    assert "font-weight:700" not in body_tag
    assert "line-height:1.5" in body_tag


def test_render_cancer_fallback_no_headline():
    """When therapeutic_headline is empty, therapeutic_implications renders directly
    as body text — NOT as a 15px bold title."""
    from scripts.clinical_board.render import render_board_opinion_html

    body = (
        "EGFR L858R 활성화 변이가 검출되어 EGFR TKI 표적치료에 민감할 것으로 예측됩니다. "
        "1세대 또는 3세대 TKI 사용을 권고하며, 저항성 모니터링이 필요합니다."
    )
    opinion = CancerBoardOpinion(
        therapeutic_headline="",
        therapeutic_implications=body,
        therapeutic_evidence="CIViC Level A",
    )
    html = render_board_opinion_html(opinion, language="en")

    assert body in html
    body_idx = html.find(body)
    body_start = html.rfind("<", 0, body_idx)
    body_tag = html[body_start:body_idx]
    # Fallback body must not be styled as the old bold title.
    assert "font-size:15px" not in body_tag
    assert "font-weight:700" not in body_tag
    assert "line-height:1.5" in body_tag
    # Section caption stays.
    assert "Therapeutic Implications" in html


def test_render_cancer_evidence_below():
    """therapeutic_evidence still renders as the caption block below the body text."""
    from scripts.clinical_board.render import render_board_opinion_html

    body = "임상시험 등록을 권고합니다."
    evidence = "CIViC Level C — supporting evidence summary"
    opinion = CancerBoardOpinion(
        therapeutic_headline="Stage IV PDAC headline",
        therapeutic_implications=body,
        therapeutic_evidence=evidence,
    )
    html = render_board_opinion_html(opinion, language="en")

    assert evidence in html
    body_idx = html.find(body)
    evidence_idx = html.find(evidence)
    assert evidence_idx > body_idx, "evidence should appear after the body paragraph"

    evidence_start = html.rfind("<", 0, evidence_idx)
    evidence_tag = html[evidence_start:evidence_idx]
    assert "font-size:11px" in evidence_tag
    assert "#0F766E" in evidence_tag
