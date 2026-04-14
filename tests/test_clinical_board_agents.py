"""Tests for Clinical Board domain agents and Board Chair."""

from unittest.mock import MagicMock

from scripts.clinical_board.models import AgentOpinion, BoardOpinion


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SAMPLE_CASE_BRIEFING = """
## 환자 정보
Sample ID: SAMPLE-001
분석 모드: rare_disease

## 변이 목록
1. TP53 c.524G>T (p.R175L)
   - Classification: Likely Pathogenic
   - ACMG Evidence: PM1, PM2, PP3, PP5
   - ClinVar: Likely Pathogenic (1 star)
   - REVEL: 0.89, CADD: 28.5
   - gnomAD AF: 0.00001
   - Consequence: missense_variant
   - Zygosity: heterozygous

2. BRCA2 c.7397T>A (p.V2466D)
   - Classification: VUS
   - ACMG Evidence: PM2, PP3
   - ClinVar: Uncertain Significance
   - REVEL: 0.62, CADD: 22.1
   - gnomAD AF: 0.00003
   - Consequence: missense_variant
   - Zygosity: heterozygous

## HPO 매칭
HP:0002664 — Neoplasm (match score: 0.85)

## PGx 결과
CYP2D6 *1/*10 — Intermediate Metabolizer
"""

SAMPLE_CASE_BRIEFING_NO_PGX = """
## 환자 정보
Sample ID: SAMPLE-002
분석 모드: rare_disease

## 변이 목록
1. CFTR c.1521_1523delCTT (p.F508del)
   - Classification: Pathogenic
   - Consequence: inframe_deletion
   - Zygosity: homozygous

## HPO 매칭
HP:0002110 — Bronchiectasis (match score: 0.92)

## PGx 결과
해당 없음
"""


def _mock_agent_response():
    """Standard mock response for domain agents."""
    return {
        "findings": [
            {
                "finding": "TP53 p.R175L은 DNA-binding domain에 위치",
                "evidence": "IARC TP53 database",
                "confidence": "high",
            }
        ],
        "recommendations": ["기능 연구 결과 추적 권고"],
        "concerns": [],
        "references": ["PMID:20522432"],
        "confidence": "high",
    }


def _mock_board_response():
    """Standard mock response for Board Chair."""
    return {
        "primary_diagnosis": "Li-Fraumeni 증후군 의심",
        "primary_diagnosis_evidence": "TP53 p.R175L Likely Pathogenic 변이",
        "differential_diagnoses": [
            {
                "diagnosis": "유전성 유방-난소암 증후군",
                "likelihood": "moderate",
                "evidence": "BRCA2 VUS 존재",
            }
        ],
        "key_findings": [
            "TP53 DNA-binding domain 변이 확인",
            "BRCA2 VUS 추가 평가 필요",
        ],
        "recommendations": [
            "유전 상담 권고",
            "가족력 확인 및 분리분석 권고",
        ],
        "agent_consensus": "majority",
        "dissenting_opinions": ["BRCA2 VUS의 임상적 의미에 대한 이견"],
        "follow_up": ["BRCA2 VUS 재분류 추적", "정기 암 검진"],
        "confidence": "high",
    }


# ---------------------------------------------------------------------------
# Variant Pathologist
# ---------------------------------------------------------------------------


def test_variant_pathologist_analyze():
    mock_client = MagicMock()
    mock_client.generate_json.return_value = _mock_agent_response()

    from scripts.clinical_board.agents.variant_pathologist import VariantPathologist

    agent = VariantPathologist(client=mock_client)
    opinion = agent.analyze(SAMPLE_CASE_BRIEFING)

    assert opinion.agent_name == "Variant Pathologist"
    assert opinion.domain == "variant_pathology"
    assert len(opinion.findings) == 1
    assert opinion.findings[0]["finding"] == "TP53 p.R175L은 DNA-binding domain에 위치"
    assert opinion.confidence == "high"
    assert len(opinion.references) == 1
    assert "PMID:20522432" in opinion.references
    mock_client.generate_json.assert_called_once()


# ---------------------------------------------------------------------------
# Disease Geneticist
# ---------------------------------------------------------------------------


def test_disease_geneticist_analyze():
    mock_client = MagicMock()
    mock_client.generate_json.return_value = {
        "findings": [
            {
                "finding": "TP53 변이와 Li-Fraumeni 증후군의 강한 연관성",
                "evidence": "OMIM #151623",
                "confidence": "high",
            },
            {
                "finding": "BRCA2 VUS는 현재 질환 연관성 불확실",
                "evidence": "ClinVar VUS 분류",
                "confidence": "low",
            },
        ],
        "recommendations": ["가족 내 분리분석 권고", "HPO 재평가 권고"],
        "concerns": ["BRCA2 VUS의 복합 이형접합 가능성 배제 필요"],
        "references": ["PMID:20301425", "OMIM:151623"],
        "confidence": "high",
    }

    from scripts.clinical_board.agents.disease_geneticist import DiseaseGeneticist

    agent = DiseaseGeneticist(client=mock_client)
    opinion = agent.analyze(SAMPLE_CASE_BRIEFING)

    assert opinion.agent_name == "Disease Geneticist"
    assert opinion.domain == "disease_genetics"
    assert len(opinion.findings) == 2
    assert len(opinion.recommendations) == 2
    assert len(opinion.concerns) == 1
    assert opinion.confidence == "high"
    mock_client.generate_json.assert_called_once()


# ---------------------------------------------------------------------------
# PGx Specialist
# ---------------------------------------------------------------------------


def test_pgx_specialist_analyze():
    mock_client = MagicMock()
    mock_client.generate_json.return_value = {
        "findings": [
            {
                "finding": "CYP2D6 *1/*10 — Intermediate Metabolizer로 코데인 등 CYP2D6 기질 약물 대사 저하 예상",
                "evidence": "CPIC guideline for codeine/CYP2D6",
                "confidence": "high",
            }
        ],
        "recommendations": ["코데인 대체 진통제 사용 권고", "타목시펜 용량 조절 고려"],
        "concerns": ["한국인에서 CYP2D6*10 고빈도 (약 40%)"],
        "references": ["PMID:24458010"],
        "confidence": "high",
    }

    from scripts.clinical_board.agents.pgx_specialist import PGxSpecialist

    agent = PGxSpecialist(client=mock_client)
    opinion = agent.analyze(SAMPLE_CASE_BRIEFING)

    assert opinion.agent_name == "PGx Specialist"
    assert opinion.domain == "pharmacogenomics"
    assert len(opinion.findings) == 1
    assert len(opinion.recommendations) == 2
    assert opinion.confidence == "high"
    mock_client.generate_json.assert_called_once()


def test_pgx_specialist_no_pgx_data():
    """When no PGx data is present, the specialist should return minimal output."""
    mock_client = MagicMock()
    mock_client.generate_json.return_value = {
        "findings": [
            {
                "finding": "PGx 관련 변이가 확인되지 않음",
                "evidence": "케이스 데이터에 PGx 결과 없음",
                "confidence": "high",
            }
        ],
        "recommendations": [],
        "concerns": [],
        "references": [],
        "confidence": "high",
    }

    from scripts.clinical_board.agents.pgx_specialist import PGxSpecialist

    agent = PGxSpecialist(client=mock_client)
    opinion = agent.analyze(SAMPLE_CASE_BRIEFING_NO_PGX)

    assert opinion.agent_name == "PGx Specialist"
    assert len(opinion.findings) == 1
    assert "PGx" in opinion.findings[0]["finding"]
    assert len(opinion.recommendations) == 0
    mock_client.generate_json.assert_called_once()


# ---------------------------------------------------------------------------
# Literature Analyst
# ---------------------------------------------------------------------------


def test_literature_analyst_analyze():
    mock_client = MagicMock()
    mock_client.generate_json.return_value = {
        "findings": [
            {
                "finding": "TP53 p.R175L은 다수의 기능 연구에서 dominant-negative 효과 확인",
                "evidence": "Kato et al. 2003 — yeast-based functional assay",
                "confidence": "high",
            },
            {
                "finding": "BRCA2 p.V2466D에 대한 기능 연구는 제한적",
                "evidence": "소수 증례 보고만 존재",
                "confidence": "low",
            },
        ],
        "recommendations": ["BRCA2 VUS에 대한 기능 연구 추적 필요"],
        "concerns": [],
        "references": ["PMID:12826609", "PMID:15146165"],
        "confidence": "moderate",
    }

    from scripts.clinical_board.agents.literature_analyst import LiteratureAnalyst

    agent = LiteratureAnalyst(client=mock_client)
    opinion = agent.analyze(SAMPLE_CASE_BRIEFING)

    assert opinion.agent_name == "Literature Analyst"
    assert opinion.domain == "literature_evidence"
    assert len(opinion.findings) == 2
    assert len(opinion.references) == 2
    assert any("PMID" in ref for ref in opinion.references)
    assert opinion.confidence == "moderate"
    mock_client.generate_json.assert_called_once()


# ---------------------------------------------------------------------------
# Board Chair
# ---------------------------------------------------------------------------


def test_board_chair_synthesize():
    mock_client = MagicMock()
    mock_client.generate_json.return_value = _mock_board_response()

    agent_opinions = [
        AgentOpinion(
            agent_name="Variant Pathologist",
            domain="variant_pathology",
            findings=[{"finding": "TP53 DNA-binding domain 변이", "evidence": "IARC DB", "confidence": "high"}],
            recommendations=["기능 연구 추적"],
            confidence="high",
        ),
        AgentOpinion(
            agent_name="Disease Geneticist",
            domain="disease_genetics",
            findings=[{"finding": "Li-Fraumeni 증후군 가능성", "evidence": "OMIM", "confidence": "high"}],
            recommendations=["가족 분리분석"],
            confidence="high",
        ),
        AgentOpinion(
            agent_name="PGx Specialist",
            domain="pharmacogenomics",
            findings=[{"finding": "CYP2D6 IM 확인", "evidence": "CPIC", "confidence": "high"}],
            recommendations=["코데인 대체"],
            confidence="high",
        ),
        AgentOpinion(
            agent_name="Literature Analyst",
            domain="literature_evidence",
            findings=[{"finding": "TP53 R175L 기능 연구 다수", "evidence": "Kato 2003", "confidence": "high"}],
            references=["PMID:12826609"],
            confidence="high",
        ),
    ]

    from scripts.clinical_board.agents.board_chair import BoardChair

    chair = BoardChair(client=mock_client)
    result = chair.synthesize(SAMPLE_CASE_BRIEFING, agent_opinions)

    assert isinstance(result, BoardOpinion)
    assert result.primary_diagnosis == "Li-Fraumeni 증후군 의심"
    assert result.primary_diagnosis_evidence != ""
    assert len(result.differential_diagnoses) >= 1
    assert len(result.key_findings) >= 1
    assert len(result.recommendations) >= 1
    assert result.agent_consensus == "majority"
    assert len(result.agent_opinions) == 4
    assert len(result.dissenting_opinions) >= 1
    assert len(result.follow_up) >= 1
    assert result.confidence == "high"
    assert "AI-Generated" in result.disclaimer
    mock_client.generate_json.assert_called_once()


# ---------------------------------------------------------------------------
# Error Handling
# ---------------------------------------------------------------------------


def test_agent_handles_error():
    """When the LLM client raises an exception, the agent returns an error AgentOpinion."""
    mock_client = MagicMock()
    mock_client.generate_json.side_effect = ConnectionError("Ollama is offline")

    from scripts.clinical_board.agents.variant_pathologist import VariantPathologist

    agent = VariantPathologist(client=mock_client)
    opinion = agent.analyze(SAMPLE_CASE_BRIEFING)

    assert opinion.agent_name == "Variant Pathologist"
    assert opinion.confidence == "low"
    assert len(opinion.findings) == 1
    assert "failed" in opinion.findings[0]["finding"].lower()


def test_agent_handles_malformed_json():
    """When the LLM returns a non-JSON string, the agent parses gracefully."""
    mock_client = MagicMock()
    mock_client.generate_json.return_value = "이것은 JSON이 아닌 텍스트 응답입니다."

    from scripts.clinical_board.agents.disease_geneticist import DiseaseGeneticist

    agent = DiseaseGeneticist(client=mock_client)
    opinion = agent.analyze(SAMPLE_CASE_BRIEFING)

    assert opinion.agent_name == "Disease Geneticist"
    assert opinion.confidence == "low"
    assert len(opinion.findings) == 1
    assert opinion.raw_response != ""


def test_board_chair_handles_error():
    """When the Board Chair LLM call fails, it returns a fallback BoardOpinion."""
    mock_client = MagicMock()
    mock_client.generate_json.side_effect = ConnectionError("Ollama is offline")

    agent_opinions = [
        AgentOpinion(agent_name="Variant Pathologist", domain="variant_pathology", confidence="high"),
    ]

    from scripts.clinical_board.agents.board_chair import BoardChair

    chair = BoardChair(client=mock_client)
    result = chair.synthesize(SAMPLE_CASE_BRIEFING, agent_opinions)

    assert isinstance(result, BoardOpinion)
    assert "실패" in result.primary_diagnosis
    assert result.confidence == "low"
    assert len(result.agent_opinions) == 1
    assert "AI-Generated" in result.disclaimer


# ---------------------------------------------------------------------------
# System Prompt Validation
# ---------------------------------------------------------------------------


def test_all_agents_have_korean_system_prompt():
    """All domain agents must have system prompts containing Korean text."""
    mock_client = MagicMock()

    from scripts.clinical_board.agents.variant_pathologist import VariantPathologist
    from scripts.clinical_board.agents.disease_geneticist import DiseaseGeneticist
    from scripts.clinical_board.agents.pgx_specialist import PGxSpecialist
    from scripts.clinical_board.agents.literature_analyst import LiteratureAnalyst

    agents = [
        VariantPathologist(client=mock_client),
        DiseaseGeneticist(client=mock_client),
        PGxSpecialist(client=mock_client),
        LiteratureAnalyst(client=mock_client),
    ]

    for agent in agents:
        prompt = agent.system_prompt
        # Check that the system prompt contains Korean characters
        has_korean = any("\uac00" <= ch <= "\ud7a3" for ch in prompt)
        assert has_korean, f"{agent.agent_name} system prompt does not contain Korean text"
        # Check the disclaimer about not modifying classification engine results
        assert "분류 엔진" in prompt, (
            f"{agent.agent_name} system prompt must include disclaimer about deterministic classification engine"
        )
        # Check it instructs to respond in Korean
        assert "한국어" in prompt, f"{agent.agent_name} system prompt must instruct Korean response"


def test_board_chair_has_korean_system_prompt():
    """Board Chair must have a system prompt containing Korean text."""
    from scripts.clinical_board.agents.board_chair import SYSTEM_PROMPT

    has_korean = any("\uac00" <= ch <= "\ud7a3" for ch in SYSTEM_PROMPT)
    assert has_korean, "Board Chair system prompt does not contain Korean text"
    assert "분류 엔진" in SYSTEM_PROMPT
    assert "한국어" in SYSTEM_PROMPT


# ---------------------------------------------------------------------------
# Agent Properties
# ---------------------------------------------------------------------------


def test_agent_names_and_domains():
    """Verify all agents have correct names and domains."""
    mock_client = MagicMock()

    from scripts.clinical_board.agents.variant_pathologist import VariantPathologist
    from scripts.clinical_board.agents.disease_geneticist import DiseaseGeneticist
    from scripts.clinical_board.agents.pgx_specialist import PGxSpecialist
    from scripts.clinical_board.agents.literature_analyst import LiteratureAnalyst

    expected = [
        (VariantPathologist, "Variant Pathologist", "variant_pathology"),
        (DiseaseGeneticist, "Disease Geneticist", "disease_genetics"),
        (PGxSpecialist, "PGx Specialist", "pharmacogenomics"),
        (LiteratureAnalyst, "Literature Analyst", "literature_evidence"),
    ]

    for cls, name, domain in expected:
        agent = cls(client=mock_client)
        assert agent.agent_name == name
        assert agent.domain == domain


def test_agent_empty_response():
    """When the LLM returns an empty dict, the agent handles gracefully."""
    mock_client = MagicMock()
    mock_client.generate_json.return_value = {}

    from scripts.clinical_board.agents.variant_pathologist import VariantPathologist

    agent = VariantPathologist(client=mock_client)
    opinion = agent.analyze(SAMPLE_CASE_BRIEFING)

    assert opinion.agent_name == "Variant Pathologist"
    assert opinion.confidence == "moderate"  # default
    assert opinion.findings == []
    assert opinion.recommendations == []


# ---------------------------------------------------------------------------
# Domain Sheet + Prior Knowledge (AI Board v2 — Task 2)
# ---------------------------------------------------------------------------


def test_build_prompt_with_domain_sheet():
    """Domain sheet is included in the prompt when provided."""
    mock_client = MagicMock()
    from scripts.clinical_board.agents.variant_pathologist import VariantPathologist

    agent = VariantPathologist(client=mock_client)
    prompt = agent._build_prompt(
        "case briefing text",
        domain_sheet="== DOMAIN DATA ==\nClinVar: Pathogenic",
    )
    assert "DOMAIN DATA" in prompt
    assert "ClinVar: Pathogenic" in prompt


def test_build_prompt_without_domain_sheet():
    """Prompt works without domain sheet (backward compatible)."""
    mock_client = MagicMock()
    from scripts.clinical_board.agents.variant_pathologist import VariantPathologist

    agent = VariantPathologist(client=mock_client)
    prompt = agent._build_prompt("case briefing text")
    assert "case briefing text" in prompt


def test_build_prompt_with_prior_knowledge():
    """Prior knowledge is included as separate section."""
    mock_client = MagicMock()
    from scripts.clinical_board.agents.variant_pathologist import VariantPathologist

    agent = VariantPathologist(client=mock_client)
    prompt = agent._build_prompt(
        "briefing",
        prior_knowledge="== PRIOR BOARD KNOWLEDGE ==\n3 prior cases",
    )
    assert "PRIOR BOARD KNOWLEDGE" in prompt
    assert "3 prior cases" in prompt


# ---------------------------------------------------------------------------
# Runner mode-based agent selection (Task 7)
# ---------------------------------------------------------------------------


def test_runner_loads_cancer_agents():
    """Cancer mode loads therapeutic target, tumor genomics, pgx, clinical evidence."""
    from scripts.clinical_board.runner import _load_agents

    client = MagicMock()
    agents = _load_agents(client, "test-model", "en", mode="cancer")
    names = [a.agent_name for a in agents]
    assert "Therapeutic Target Analyst" in names
    assert "Tumor Genomics Specialist" in names
    assert "Clinical Evidence Analyst" in names
    assert len(agents) == 4


def test_runner_loads_rare_disease_agents_default():
    """Default mode loads rare disease agents (backward compatible)."""
    from scripts.clinical_board.runner import _load_agents

    client = MagicMock()
    agents = _load_agents(client, "test-model", "en", mode="rare-disease")
    names = [a.agent_name for a in agents]
    assert "Variant Pathologist" in names
    assert "Disease Geneticist" in names
    assert len(agents) == 4


# ---------------------------------------------------------------------------
# Runner selector wiring (Task 4)
# ---------------------------------------------------------------------------


def test_runner_uses_selector_and_records_metadata(monkeypatch):
    """run_clinical_board runs the selector, pre-populates report_data, and
    propagates selection_metadata onto the returned BoardOpinion."""
    from scripts.clinical_board import runner as runner_mod
    from scripts.clinical_board.models import (
        AgentOpinion,
        BoardOpinion,
    )

    captured: dict = {}

    # --- Mock OllamaClient: always-available, always-model-loaded ---
    class _FakeClient:
        def __init__(self, *a, **kw):
            pass

        def is_available(self):
            return True

        def has_model(self, m):
            return True

    monkeypatch.setattr(runner_mod, "OllamaClient", _FakeClient)

    # --- Mock domain agents: one fake agent, records the variant list it sees ---
    fake_agent = MagicMock()
    fake_agent.agent_name = "FakeAgent"
    fake_agent.domain = "fake_domain"
    fake_agent.analyze.return_value = AgentOpinion(agent_name="FakeAgent", domain="fake_domain", confidence="high")
    monkeypatch.setattr(runner_mod, "_load_agents", lambda *a, **kw: [fake_agent])

    def _fake_build_domain_sheet(domain, mode, variants, report_data):
        captured["domain_variants"] = list(variants)
        return "fake sheet"

    monkeypatch.setattr(runner_mod, "build_domain_sheet", _fake_build_domain_sheet)

    # --- Mock chair: returns a minimal BoardOpinion ---
    fake_chair = MagicMock()
    fake_chair.synthesize.return_value = BoardOpinion(
        primary_diagnosis="test diagnosis",
        confidence="moderate",
    )
    monkeypatch.setattr(runner_mod, "_load_chair", lambda *a, **kw: fake_chair)

    # --- KB disabled to keep the test hermetic ---
    def _fake_get(key, default=None):
        if key.startswith("knowledge_base"):
            return False if key.endswith(".enabled") else ""
        return default

    monkeypatch.setattr(runner_mod, "get", _fake_get)

    # Mix of Tier I (MUST), Tier III passenger (EXCLUDED), Benign (EXCLUDED)
    report_data = {
        "sample_id": "TEST-SEL",
        "variants": [
            {
                "gene": "EGFR",
                "variant": "chr7:1",
                "classification": "VUS",
                "tier": "Tier I",
                "hgvsp": "p.Leu858Arg",
                "consequence": "missense_variant",
                "cancer_gene_type": "",
                "oncokb_level": "",
                "hpo_score": 0,
                "gnomad_af": None,
                "variant_type": "SNV",
            },
            {
                "gene": "RANDOM",
                "variant": "chr1:1",
                "classification": "VUS",
                "tier": "Tier IV",
                "hgvsp": "",
                "consequence": "missense_variant",
                "cancer_gene_type": "",
                "oncokb_level": "",
                "hpo_score": 0,
                "gnomad_af": None,
                "variant_type": "SNV",
            },
            {
                "gene": "NOISE",
                "variant": "chr2:1",
                "classification": "Benign",
                "tier": "Tier IV",
                "hgvsp": "",
                "consequence": "missense_variant",
                "cancer_gene_type": "",
                "oncokb_level": "",
                "hpo_score": 0,
                "gnomad_af": None,
                "variant_type": "SNV",
            },
        ],
        "summary": {"total": 3},
    }

    opinion = runner_mod.run_clinical_board(report_data, mode="cancer")

    assert opinion is not None
    # Selector populated the report_data side-channel
    assert "_board_variants" in report_data
    assert "_board_selection_metadata" in report_data
    selected_genes = {v["gene"] for v in report_data["_board_variants"]}
    assert selected_genes == {"EGFR"}
    meta = report_data["_board_selection_metadata"]
    assert meta["total_input"] == 3
    assert meta["selected"] == 1
    assert meta["must_included"] == 1
    # Domain sheet saw the board list, not the raw one
    assert [v["gene"] for v in captured["domain_variants"]] == ["EGFR"]
    # Metadata propagated onto the opinion for render.py
    assert opinion.selection_metadata is not None
    assert opinion.selection_metadata["selected"] == 1
