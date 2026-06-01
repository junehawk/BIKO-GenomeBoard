"""v2.8 Clinical Board agent-prompt + deliberation contract tests.

Pins the v2.8 review fixes so they cannot silently regress:
- a shared scientific-grounding clause reaches EVERY domain agent's prompt
  (both KO and EN), closing the anti-hallucination gap the review found;
- the rare-disease Chair roster is 3 specialists (the v2.5.5 PGx removal),
  not the stale "four specialists ... Pharmacogenomics";
- the PGx agent only asserts a metaboliser phenotype with a PharmCAT call
  (consistency with the v2.7 PGX-02 contract);
- the Round-2 deliberation `revise()` exists and falls back to the Round-1
  opinion on a degenerate/failed revision (never loses signal).

None of these run the LLM.
"""

from __future__ import annotations

import pytest

from scripts.clinical_board.agents.base import GROUNDING_CLAUSE_EN, GROUNDING_CLAUSE_KO
from scripts.clinical_board.agents.clinical_evidence import ClinicalEvidenceAnalyst
from scripts.clinical_board.agents.disease_geneticist import DiseaseGeneticist
from scripts.clinical_board.agents.literature_analyst import LiteratureAnalyst
from scripts.clinical_board.agents.pgx_specialist import PGxSpecialist
from scripts.clinical_board.agents.therapeutic_target import TherapeuticTargetAnalyst
from scripts.clinical_board.agents.tumor_genomics import TumorGenomicsSpecialist
from scripts.clinical_board.agents.variant_pathologist import VariantPathologist
from scripts.clinical_board.models import AgentOpinion

_AGENT_CLASSES = [
    TherapeuticTargetAnalyst,
    TumorGenomicsSpecialist,
    PGxSpecialist,
    ClinicalEvidenceAnalyst,
    VariantPathologist,
    DiseaseGeneticist,
    LiteratureAnalyst,
]


def _agent(cls, language: str):
    a = cls.__new__(cls)
    a.language = language
    return a


def _norm(s: str) -> str:
    """Collapse all whitespace so substring checks ignore prompt line-wrapping."""
    return " ".join(s.split())


# ── Grounding clause reaches every agent's built prompt ───────────────────────


@pytest.mark.parametrize("cls", _AGENT_CLASSES)
def test_grounding_clause_in_built_prompt_en(cls):
    prompt = _norm(_agent(cls, "en")._build_prompt("CASE BRIEFING"))
    assert "Scientific Grounding" in prompt
    assert "not determinable from provided data" in prompt


@pytest.mark.parametrize("cls", _AGENT_CLASSES)
def test_grounding_clause_in_built_prompt_ko(cls):
    prompt = _norm(_agent(cls, "ko")._build_prompt("케이스 정보"))
    assert "과학적 근거" in prompt
    assert "not determinable from provided data" in prompt


def test_grounding_clause_forbids_fabrication_categories():
    for clause in (GROUNDING_CLAUSE_EN, GROUNDING_CLAUSE_KO):
        low = _norm(clause).lower()
        for kw in ("pmid", "pathway", "driver", "allele frequenc", "metabolizer"):
            assert kw in low or "빈도" in clause


# ── Chair roster is 3 (rare-disease), PGx removed ─────────────────────────────


def test_rare_chair_prompt_says_three_not_four():
    from scripts.clinical_board.agents.board_chair import SYSTEM_PROMPT_EN, SYSTEM_PROMPT_KO

    assert "three domain specialists" in SYSTEM_PROMPT_EN
    assert "four domain specialists" not in SYSTEM_PROMPT_EN
    assert "Pharmacogenomics" not in SYSTEM_PROMPT_EN
    assert "3명의 도메인 전문의" in SYSTEM_PROMPT_KO
    assert "4명의" not in SYSTEM_PROMPT_KO
    assert "약물유전체" not in SYSTEM_PROMPT_KO


# ── PGx reconciled with PGX-02 (no phenotype assertion without PharmCAT) ───────


def test_pgx_prompt_gates_phenotype_on_pharmcat():
    en = _agent(PGxSpecialist, "en").system_prompt
    ko = _agent(PGxSpecialist, "ko").system_prompt
    assert "PharmCAT diplotype call" in en
    assert "PharmCAT diplotype" in ko


# ── Round-2 deliberation: revise() exists + never loses signal ────────────────


class _FakeClient:
    def __init__(self, response):
        self._response = response
        self.calls = 0

    def generate_json(self, **kwargs):
        self.calls += 1
        return self._response


def _opinion(name="Variant Pathologist", findings=None):
    return AgentOpinion(
        agent_name=name,
        domain="variant_pathology",
        findings=findings
        if findings is not None
        else [{"finding": "R1 finding", "evidence": "x", "confidence": "moderate"}],
        confidence="moderate",
    )


def test_revise_returns_round1_on_empty_revision():
    """A degenerate ({}) revision must fall back to the Round-1 opinion."""
    agent = VariantPathologist.__new__(VariantPathologist)
    agent.language = "en"
    agent.model = "test"
    agent.client = _FakeClient({})  # empty → no findings/recs
    own = _opinion()
    out = agent.revise("CASE", own, [_opinion("Disease Geneticist")])
    assert out is own  # exact same object — fell back


def test_revise_returns_revised_opinion_when_present():
    agent = VariantPathologist.__new__(VariantPathologist)
    agent.language = "en"
    agent.model = "test"
    agent.client = _FakeClient(
        {
            "findings": [
                {"finding": "revised after peer review", "evidence": "in-silico: REVEL=0.9", "confidence": "high"}
            ],
            "confidence": "high",
        }
    )
    own = _opinion()
    out = agent.revise("CASE", own, [_opinion("Literature Analyst")])
    assert out is not own
    assert out.confidence == "high"
    assert "revised after peer review" in out.findings[0]["finding"]


def test_revision_prompt_contains_peer_block_and_grounding():
    agent = VariantPathologist.__new__(VariantPathologist)
    agent.language = "en"
    prompt = agent._build_revision_prompt(
        "CASE",
        _opinion(),
        [_opinion("Disease Geneticist", [{"finding": "peer says X", "evidence": "HPO", "confidence": "high"}])],
    )
    assert "Peer Round-1 Opinions" in prompt
    assert "peer says X" in prompt
    assert "DEFER" in prompt
    assert "Scientific Grounding" in prompt  # grounding contract carried into R2
