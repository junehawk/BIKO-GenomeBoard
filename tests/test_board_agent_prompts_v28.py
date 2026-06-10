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


# ── Phase 0.3: runner ACTUALLY invokes revise() (the wiring that was lost) ────


def _fake_agent(name, domain):
    class _FA:
        def __init__(self):
            self.agent_name = name
            self.domain = domain
            self.analyze_calls = 0
            self.revise_calls = 0

        def analyze(self, *a, **k):
            self.analyze_calls += 1
            return AgentOpinion(
                agent_name=self.agent_name,
                domain=self.domain,
                findings=[{"finding": "f", "evidence": "e", "confidence": "moderate"}],
                confidence="moderate",
            )

        def revise(self, case_briefing, own_opinion, peer_opinions, **k):
            self.revise_calls += 1
            return own_opinion

    return _FA()


def _run_board(monkeypatch, rounds):
    from unittest.mock import MagicMock

    import scripts.clinical_board.runner as runner_mod
    from scripts.clinical_board.models import CancerBoardOpinion

    fa1 = _fake_agent("Tumor Genomics Specialist", "tumor_genomics")
    fa2 = _fake_agent("Clinical Evidence Analyst", "clinical_evidence")
    agents = [fa1, fa2]

    fake_client = MagicMock()
    fake_client.is_available.return_value = True
    fake_client.has_model.return_value = True

    fake_chair = MagicMock()
    fake_chair.synthesize.return_value = CancerBoardOpinion(
        therapeutic_implications="ok", agent_consensus="majority", confidence="moderate"
    )

    real_get = runner_mod.get

    def fake_get(key, default=None):
        if key == "clinical_board.deliberation_rounds":
            return rounds
        return real_get(key, default)

    monkeypatch.setattr(runner_mod, "OllamaClient", lambda *a, **kw: fake_client, raising=False)
    monkeypatch.setattr(runner_mod, "_load_agents", lambda *a, **kw: agents, raising=False)
    monkeypatch.setattr(runner_mod, "_load_chair", lambda *a, **kw: fake_chair, raising=False)
    monkeypatch.setattr(runner_mod, "curate_treatments", lambda variants, **kw: {}, raising=False)
    monkeypatch.setattr(runner_mod, "get", fake_get, raising=False)

    report_data = {
        "sample_id": "DELIB-TEST",
        "variants": [
            {
                "gene": "TP53",
                "variant": "17:7675088:C:A",
                "chrom": "chr17",
                "pos": 7675088,
                "ref": "C",
                "alt": "A",
                "classification": "Likely Pathogenic",
                "tier": "Tier I",
                "consequence": "Missense",
                "cancer_gene_type": "TSG",
                "oncokb_level": "",
                "variant_type": "SNV",
            },
            {
                "gene": "BRCA2",
                "variant": "13:32356550:C:T",
                "chrom": "chr13",
                "pos": 32356550,
                "ref": "C",
                "alt": "T",
                "classification": "VUS",
                "tier": "Tier III",
                "consequence": "Missense",
                "cancer_gene_type": "TSG",
                "oncokb_level": "",
                "variant_type": "SNV",
            },
        ],
        "summary": {"total": 2},
    }
    runner_mod.run_clinical_board(report_data, mode="cancer")
    return fa1, fa2


def test_runner_invokes_revise_when_deliberation_enabled(monkeypatch):
    """deliberation_rounds=2 → the runner must call each agent's revise() once.

    This is the wiring that was lost in v2.8 (revise() existed but the runner
    never called it, so the feature was dead code). Guards against re-orphaning.
    """
    fa1, fa2 = _run_board(monkeypatch, rounds=2)
    assert (fa1.analyze_calls, fa2.analyze_calls) == (1, 1)
    assert (fa1.revise_calls, fa2.revise_calls) == (1, 1)


def test_runner_skips_revise_when_single_round(monkeypatch):
    """deliberation_rounds=1 → single-round flow, revise() not called."""
    fa1, fa2 = _run_board(monkeypatch, rounds=1)
    assert (fa1.analyze_calls, fa2.analyze_calls) == (1, 1)
    assert (fa1.revise_calls, fa2.revise_calls) == (0, 0)


@pytest.mark.parametrize("lang", ["en", "ko"])
def test_cancer_chair_prompt_does_not_invite_msi(lang):
    """BIKO never computes MSI, so the cancer Chair prompt must not invite MSI reasoning
    (the LLM would otherwise fabricate an MSI status); TMB (which IS computed) stays (T2-01)."""
    from scripts.clinical_board.agents.board_chair import get_system_prompt

    prompt = get_system_prompt("cancer", lang)
    assert "MSI" not in prompt
    assert "TMB" in prompt
