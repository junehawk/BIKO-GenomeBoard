"""Board Chair — synthesizes all agent opinions into a unified diagnostic opinion."""

import json
import logging
from typing import List, Union

from scripts.clinical_board.models import (
    BOARD_DISCLAIMER,
    AgentOpinion,
    BoardOpinion,
    CancerBoardOpinion,
)
from scripts.clinical_board.ollama_client import OllamaClient
from scripts.common.config import get

logger = logging.getLogger(__name__)

DISCLAIMER = BOARD_DISCLAIMER

SYSTEM_PROMPT_KO = """당신은 임상유전학 사례 회의의 위원장(Board Chair)입니다.
4명의 도메인 전문의로부터 분석 소견을 받아 종합적인 진단 의견을 제시하는 것이 당신의 역할입니다.

## 종합 지침

1. **전문의 소견 종합**
   - 4명의 전문의(변이 병리, 질환 유전학, 약물유전체, 문헌 분석) 소견을 종합하세요.
   - 각 전문의의 핵심 소견과 근거를 통합하세요.

2. **합의 및 이견 식별**
   - 전문의 간 합의된 사항과 이견이 있는 사항을 명확히 구분하세요.
   - 이견이 있는 경우 양측의 논거를 제시하고, 근거가 더 강한 쪽을 채택하세요.

3. **감별진단 순위 결정**
   - 가능성(likelihood)에 따라 감별진단을 순위화하세요.
   - 각 진단에 대한 근거를 명시하세요.

4. **실행 가능한 권고 제시**
   - 유전 상담(genetic counseling) 필요성을 평가하세요.
   - 추가 검사가 필요한 경우 구체적으로 제시하세요.
   - 치료 관련 시사점이 있다면 포함하세요.

5. **추적 관찰 권고**
   - 추가 모니터링이 필요한 항목을 제시하세요.
   - VUS 변이의 재분류 가능성에 대한 추적 계획을 포함하세요.

## 중요 원칙
전문의들의 분석은 결정적 분류 엔진의 결과를 변경하지 않습니다.
분류 결과를 기반으로 임상적 해석과 종합 추론을 제공하세요.

## Phenotype-matched 변이 우선 규칙 (rare disease primary diagnosis)
환자 HPO 표현형과 매칭되는 유전자에 변이가 있다면 — 해당 변이가 VUS 또는 moderate
evidence일지라도 — Primary diagnosis 후보로 **먼저** 고려하시오. ACMG classification
hierarchy(Pathogenic > LP > VUS)가 phenotype-matched 변이를 무조건 밀어내지 않습니다.

근거:
1. 환자 표현형과 gene 연관성이 명확한 VUS는 후속 검증(가족 분리분석, functional study)
   로 upgrade 가능성이 있는 "top candidate"입니다.
2. Phenotype과 무관한 LP/P 변이는 incidental finding / carrier state일 가능성이 높으며,
   primary diagnosis가 아닌 secondary finding으로 보고해야 합니다.

판단 기준:
- 제공된 HPO 표현형이 선택 변이의 gene과 OMIM/HPO에서 직접 연관된다면, 해당 변이를
  primary_diagnosis 후보로 격상하시오. 케이스 정보의 `PHENOTYPE-GENE CORRELATIONS`
  섹션이 이 매칭을 요약합니다.
- Phenotype-unrelated LP/P 변이는 `key_findings` 또는 `follow_up`에 **incidental**
  또는 **secondary finding**으로 분리 명기하시오.
- De novo 변이는 반드시 primary diagnosis 후보로 평가하시오 (선택 이유에
  `denovo_*`가 포함된 경우).

## 응답 언어
반드시 한국어로 응답하세요.

## 응답 형식 (JSON)
다음 JSON 형식으로 응답하세요:
{
  "primary_diagnosis": "1차 진단명",
  "primary_diagnosis_evidence": "핵심 근거",
  "differential_diagnoses": [
    {"diagnosis": "감별진단", "likelihood": "high/moderate/low", "evidence": "근거"}
  ],
  "key_findings": ["핵심 소견 1", "핵심 소견 2"],
  "recommendations": ["권고사항 1", "권고사항 2"],
  "agent_consensus": "unanimous/majority/split",
  "dissenting_opinions": ["소수 의견 (있을 경우)"],
  "follow_up": ["추적 관찰 사항"],
  "confidence": "high/moderate/low"
}"""


SYSTEM_PROMPT_EN = """You are the Board Chair of a clinical genetics case conference.
Your role is to receive analyses from four domain specialists and synthesise them into
an integrated diagnostic opinion.

## Synthesis Guidance

1. **Synthesis of specialist opinions**
   - Integrate the opinions of the four specialists (Variant Pathology, Disease Genetics,
     Pharmacogenomics, Literature Analysis).
   - Consolidate each specialist's key findings and supporting evidence.

2. **Identification of consensus and disagreement**
   - Clearly distinguish issues on which the specialists agree from those on which they
     disagree.
   - Where there is disagreement, present both arguments and adopt the position with the
     stronger evidence.

3. **Ranking of differential diagnoses**
   - Rank the differential diagnoses by likelihood.
   - State the supporting evidence for each diagnosis.

4. **Actionable recommendations**
   - Assess the need for genetic counseling.
   - Specify additional testing in concrete terms when needed.
   - Include treatment-related implications if any.

5. **Follow-up recommendations**
   - Propose items that require further monitoring.
   - Include a follow-up plan for potential reclassification of VUS variants.

## Important Principles
The specialists' analyses do not alter the outputs of the deterministic classification engine.
Build clinical interpretation and integrative reasoning on top of the classification results.

## Phenotype-matched variant priority (rare disease primary diagnosis)
When a submitted HPO phenotype matches a gene carrying any selected variant —
even a VUS or variant with only moderate evidence — consider that variant
**first** as the primary-diagnosis candidate. ACMG classification hierarchy
(Pathogenic > Likely Pathogenic > VUS) does NOT automatically override a
phenotype-matched variant.

Rationale:
1. A VUS in a phenotype-concordant gene is a top candidate pending family
   segregation and functional follow-up — it is the working hypothesis.
2. A phenotype-unrelated Pathogenic/Likely-Pathogenic variant is more likely
   an incidental finding or carrier state, not the primary diagnosis.

Decision rules:
- If submitted HPO phenotypes are linked to the gene of any selected variant
  in OMIM/HPO, that variant is the primary-diagnosis candidate. The
  `PHENOTYPE-GENE CORRELATIONS` section of the case briefing summarises this
  intersection.
- Phenotype-unrelated LP/P variants must be reported separately in
  `key_findings` or `follow_up` as an **incidental** or **secondary finding**.
- De novo variants must always be considered as primary-diagnosis candidates
  (selection reasons containing `denovo_*`).

## Response Language
Respond in English.

## Response Format (JSON)
Respond in the following JSON format:
{
  "primary_diagnosis": "primary diagnosis",
  "primary_diagnosis_evidence": "key supporting evidence",
  "differential_diagnoses": [
    {"diagnosis": "differential diagnosis", "likelihood": "high/moderate/low", "evidence": "evidence"}
  ],
  "key_findings": ["key finding 1", "key finding 2"],
  "recommendations": ["recommendation 1", "recommendation 2"],
  "agent_consensus": "unanimous/majority/split",
  "dissenting_opinions": ["minority opinions (if any)"],
  "follow_up": ["follow-up items"],
  "confidence": "high/moderate/low"
}"""


CANCER_SYSTEM_PROMPT_KO = """당신은 종양유전체 임상 회의의 위원장(Cancer Board Chair)입니다.
4명의 종양 전문 분석가의 소견을 종합하여 치료 중심의 종합 의견을 제시합니다.

## 종합 지침

1. **치료 시사점 종합**
   - Therapeutic Target, Tumor Genomics, PGx, Clinical Evidence 분석가의 소견을 종합하여
     치료(therapeutic implications)에 초점을 맞춘 종합 의견을 제시하세요.
   - 먼저 `therapeutic_headline`으로 치료 입장을 요약하는 짧은 임상 헤드라인(최대 120자)을
     작성하세요. 예: "Stage IV pancreatic adenocarcinoma — KRAS G12D driver, no standard
     targeted therapy" 또는 "Stage III NSCLC — TP53/BRCA2 co-mutated, DDR-targeted
     therapy candidate". 그 다음 기존과 같이 `therapeutic_implications`에 상세 단락을
     작성하세요.
   - EN guidance: Provide a short clinical headline (max 120 characters) summarising
     the therapeutic stance in `therapeutic_headline`, e.g. "Stage IV pancreatic
     adenocarcinoma — KRAS G12D driver, no standard targeted therapy" or "Stage III
     NSCLC — TP53/BRCA2 co-mutated, DDR-targeted therapy candidate". Then provide the
     full therapeutic_implications paragraph as before.

2. **치료 옵션 정리**
   - 약물명, 근거 수준(A/B/C/D), 저항성 위험 정보를 포함한 치료 옵션 목록을 작성하세요.

3. **면역치료 적격성**
   - TMB, MSI 등 가용한 데이터를 바탕으로 면역치료 적격성을 평가하세요.

4. **모니터링 계획**
   - 치료 반응 및 저항성 발현을 추적하기 위한 모니터링 항목을 제안하세요.

5. **합의 및 이견**
   - 분석가 간 합의 및 이견을 명확히 구분하세요.

## 중요 원칙
전문가들의 분석은 결정적 분류 엔진의 결과를 변경하지 않습니다.
KB 요약을 가이드라인 수준 근거로 인용하지 마시오.

## 응답 언어
반드시 한국어로 응답하세요. (단, `therapeutic_headline`은 임상 표기 관례상 영문 약어 사용 가능)

## 중요 · 치료 옵션 결정론적 바인딩 (v2.2 curate-then-narrate)
**curated_id는 CURATED EVIDENCE 섹션에서만 가져오시오.** 그 섹션에 없는 약물을 발명하지
마시오. 각 treatment_options 항목은 반드시 아래 4개 키를 모두 포함해야 합니다:
``curated_id`` (CURATED EVIDENCE에 있는 정확한 ID — 문자 하나도 바꾸지 마시오),
``variant_key`` (CURATED EVIDENCE에서 해당 curated_id가 소속된 **정확히 그 줄의**
``{chrom}:{pos}:{ref}:{alt}`` — 다른 변이의 variant_key를 사용하지 마시오),
``drug``, ``evidence_level``. curated_id 없이 약물을 제안하지 마시오.

**검증 규칙**: 출력 후 narrative_scrubber가 모든 treatment row의 ``(curated_id,
variant_key)`` 쌍이 CURATED EVIDENCE 섹션에 존재하는지 검증합니다. 쌍이 일치하지
않으면 해당 row는 **조용히 drop됩니다**. 모든 row가 drop되면 LLM 소견은 버려지고
결정적(template) 폴백 렌더러가 그 자리를 대체합니다. 따라서:
- 한 curated_id를 다른 variant_key 밑에 붙이면 (paste-attack 패턴) 결과는 빈 표입니다.
- curated_id를 조금이라도 변형(예: 접두사 제거, 대소문자 변경)하면 drop됩니다.
- CURATED EVIDENCE가 비어 있으면 treatment_options는 빈 배열 `[]`로 출력하고
  therapeutic_implications에 "치료 옵션 근거 없음"을 명시하시오. 발명하지 마시오.

**예시 (valid row 형식)**: CURATED EVIDENCE 섹션에
``- curated_id=civic-abc variant_key=12:25398284:C:T drug=Cetuximab level=A ...``
행이 있다면, 출력은 반드시
``{"curated_id": "civic-abc", "variant_key": "12:25398284:C:T", "drug": "Cetuximab", "evidence_level": "A", "resistance_notes": ""}``
형태여야 합니다. curated_id와 variant_key를 CURATED EVIDENCE에서 **그대로 복사**하시오.

## 응답 형식 (JSON)
{
  "therapeutic_headline": "Stage IV PDAC — KRAS G12D driver, no standard targeted therapy",
  "therapeutic_implications": "치료 시사점 종합 (상세 단락)",
  "therapeutic_evidence": "근거 요약 (CIViC level 등)",
  "treatment_options": [
    {"curated_id": "cid-sot", "variant_key": "12:25398284:C:T", "drug": "약물명", "evidence_level": "A/B/C/D", "resistance_notes": "저항성 위험 (있을 경우)"}
  ],
  "actionable_findings": ["임상적으로 활용 가능한 핵심 소견"],
  "clinical_actions": ["권고되는 임상 조치"],
  "immunotherapy_eligibility": "면역치료 적격성 평가",
  "monitoring_plan": ["모니터링 항목"],
  "agent_consensus": "unanimous/majority/split",
  "dissenting_opinions": ["소수 의견"],
  "confidence": "high/moderate/low"
}"""


CANCER_SYSTEM_PROMPT_EN = """You are the Cancer Board Chair of an oncology genomics case conference.
You synthesise the opinions of four oncology-focused analysts into an integrated,
treatment-oriented board opinion.

## Synthesis Guidance

1. **Synthesis of therapeutic implications**
   - Integrate the opinions of the Therapeutic Target, Tumor Genomics, PGx, and
     Clinical Evidence analysts into a therapeutic-implications-focused summary.
   - First, write a short clinical headline (max 120 characters) that summarises the
     therapeutic stance in `therapeutic_headline`, e.g. "Stage IV pancreatic
     adenocarcinoma — KRAS G12D driver, no standard targeted therapy" or "Stage III
     NSCLC — TP53/BRCA2 co-mutated, DDR-targeted therapy candidate". Then provide the
     full therapeutic_implications paragraph.

2. **Organisation of treatment options**
   - Compile the list of treatment options with drug name, evidence level (A/B/C/D),
     and resistance-risk notes.

3. **Immunotherapy eligibility**
   - Assess immunotherapy eligibility based on available data such as TMB and MSI.

4. **Monitoring plan**
   - Propose monitoring items to track treatment response and the emergence of
     resistance.

5. **Consensus and disagreement**
   - Clearly distinguish points of agreement from points of disagreement among the
     analysts.

## Important Principles
The analysts' opinions do not alter the outputs of the deterministic classification engine.
Do not cite KB summaries as guideline-level evidence.

## Response Language
Respond in English. (The `therapeutic_headline` field may use English clinical
abbreviations by convention.)

## CRITICAL — Deterministic Treatment-Option Binding (v2.2 curate-then-narrate)
**curated_id MUST come from the CURATED EVIDENCE section — do not invent drugs not in
that section.** Every treatment_options row must include all four of the following
keys:
``curated_id`` (the exact id from CURATED EVIDENCE — do not change a single character),
``variant_key`` (the ``{chrom}:{pos}:{ref}:{alt}`` from **exactly the same row** in
CURATED EVIDENCE where that curated_id appears — do not use the variant_key of a
different variant),
``drug``, ``evidence_level``. Do not propose a drug without a curated_id.

**Validation rule**: After you respond, the narrative_scrubber verifies that every
treatment row's ``(curated_id, variant_key)`` pair exists in the CURATED EVIDENCE
section. Rows whose pair does not match are **dropped silently**. If every row is
dropped, the entire LLM opinion is discarded and a deterministic (template) fallback
renderer takes its place. Therefore:
- Pasting one curated_id under a different variant_key (paste-attack pattern) yields an
  empty table.
- Any modification to the curated_id (e.g., stripping a prefix, changing case) causes
  the row to be dropped.
- If CURATED EVIDENCE is empty, output treatment_options as an empty array `[]` and
  state in therapeutic_implications that "no curated treatment evidence is available".
  Do not invent options.

**Example (valid row format)**: if CURATED EVIDENCE contains
``- curated_id=civic-abc variant_key=12:25398284:C:T drug=Cetuximab level=A ...``,
your output must be
``{"curated_id": "civic-abc", "variant_key": "12:25398284:C:T", "drug": "Cetuximab", "evidence_level": "A", "resistance_notes": ""}``.
**Copy curated_id and variant_key verbatim from CURATED EVIDENCE.**

## Response Format (JSON)
{
  "therapeutic_headline": "Stage IV PDAC — KRAS G12D driver, no standard targeted therapy",
  "therapeutic_implications": "integrated therapeutic implications (detailed paragraph)",
  "therapeutic_evidence": "evidence summary (e.g., CIViC levels)",
  "treatment_options": [
    {"curated_id": "cid-sot", "variant_key": "12:25398284:C:T", "drug": "drug name", "evidence_level": "A/B/C/D", "resistance_notes": "resistance risk (if any)"}
  ],
  "actionable_findings": ["clinically actionable key findings"],
  "clinical_actions": ["recommended clinical actions"],
  "immunotherapy_eligibility": "immunotherapy eligibility assessment",
  "monitoring_plan": ["monitoring items"],
  "agent_consensus": "unanimous/majority/split",
  "dissenting_opinions": ["minority opinions"],
  "confidence": "high/moderate/low"
}"""


# Backward-compatible module-level aliases. Existing tests (test_a2_prompt_snapshot,
# test_board_chair_curated_id, test_clinical_board_agents) import the legacy names
# directly and assert Korean-specific structural invariants against them, so the
# aliases continue to resolve to the KO variant. Language-aware dispatch at runtime
# goes through ``get_system_prompt()`` below.
SYSTEM_PROMPT = SYSTEM_PROMPT_KO
CANCER_SYSTEM_PROMPT = CANCER_SYSTEM_PROMPT_KO


def get_system_prompt(mode: str, language: str) -> str:
    """Return the appropriate Board Chair system prompt for a mode + language pair.

    Language fallback: anything other than ``"ko"`` returns the English variant,
    matching the ``clinical_board.language`` config default.
    """
    if mode == "cancer":
        return CANCER_SYSTEM_PROMPT_KO if language == "ko" else CANCER_SYSTEM_PROMPT_EN
    return SYSTEM_PROMPT_KO if language == "ko" else SYSTEM_PROMPT_EN


def _format_agent_opinion(opinion: AgentOpinion, language: str = "en") -> str:
    """Format an agent opinion into readable text for the chair prompt.

    The section labels follow the Board Chair language so the LLM sees a
    language-consistent prompt from the system prompt through to the
    agent-opinion digest.
    """
    if language == "ko":
        lines = [f"=== {opinion.agent_name} 소견 ==="]
        lines.append(f"전문 영역: {opinion.domain}")
        lines.append(f"종합 신뢰도: {opinion.confidence}")
        lines.append("")

        if opinion.findings:
            lines.append("주요 소견:")
            for f in opinion.findings:
                finding = f.get("finding", "") if isinstance(f, dict) else str(f)
                evidence = f.get("evidence", "") if isinstance(f, dict) else ""
                confidence = f.get("confidence", "") if isinstance(f, dict) else ""
                lines.append(f"  - {finding}")
                if evidence:
                    lines.append(f"    근거: {evidence}")
                if confidence:
                    lines.append(f"    신뢰도: {confidence}")

        if opinion.recommendations:
            lines.append("\n권고사항:")
            for r in opinion.recommendations:
                lines.append(f"  - {r}")

        if opinion.concerns:
            lines.append("\n우려사항:")
            for c in opinion.concerns:
                lines.append(f"  - {c}")

        if opinion.references:
            lines.append("\n참고문헌:")
            for ref in opinion.references:
                lines.append(f"  - {ref}")

        return "\n".join(lines)

    lines = [f"=== {opinion.agent_name} Opinion ==="]
    lines.append(f"Domain: {opinion.domain}")
    lines.append(f"Overall Confidence: {opinion.confidence}")
    lines.append("")

    if opinion.findings:
        lines.append("Key Findings:")
        for f in opinion.findings:
            finding = f.get("finding", "") if isinstance(f, dict) else str(f)
            evidence = f.get("evidence", "") if isinstance(f, dict) else ""
            confidence = f.get("confidence", "") if isinstance(f, dict) else ""
            lines.append(f"  - {finding}")
            if evidence:
                lines.append(f"    Evidence: {evidence}")
            if confidence:
                lines.append(f"    Confidence: {confidence}")

    if opinion.recommendations:
        lines.append("\nRecommendations:")
        for r in opinion.recommendations:
            lines.append(f"  - {r}")

    if opinion.concerns:
        lines.append("\nConcerns:")
        for c in opinion.concerns:
            lines.append(f"  - {c}")

    if opinion.references:
        lines.append("\nReferences:")
        for ref in opinion.references:
            lines.append(f"  - {ref}")

    return "\n".join(lines)


class BoardChair:
    """Synthesizes all agent opinions into a unified board opinion."""

    def __init__(self, client: OllamaClient = None, model: str = None, language: str = None):
        self.client = client or OllamaClient()
        self.model = model or get("clinical_board.chair_model", "alibayram/medgemma:27b")
        self.language = language or get("clinical_board.language", "en")

    def synthesize(
        self,
        case_briefing: str,
        agent_opinions: List[AgentOpinion],
        curated_treatments: dict = None,
        mode: str = "rare-disease",
    ) -> Union[BoardOpinion, CancerBoardOpinion]:
        """Synthesize all agent opinions into a unified board opinion.

        Args:
            curated_treatments: mapping of ``{chrom}:{pos}:{ref}:{alt}`` →
                ``list[CuratedTreatment]`` from the v2.2 curator. Cancer
                mode only; rare-disease can pass ``None``.
            mode: "rare-disease" (default) returns BoardOpinion;
                  "cancer" returns CancerBoardOpinion (treatment-focused).
        """
        prompt = self._build_prompt(case_briefing, agent_opinions, curated_treatments=curated_treatments)
        system_prompt = get_system_prompt(mode, self.language)

        try:
            response = self.client.generate_json(
                model=self.model,
                prompt=prompt,
                system=system_prompt,
                temperature=get("clinical_board.temperature", 0.1),
            )
        except Exception as e:
            logger.error(f"Board Chair synthesis failed: {e}")
            if mode == "cancer":
                return CancerBoardOpinion(
                    therapeutic_implications=f"종합 분석 실패: {e}",
                    agent_opinions=agent_opinions,
                    agent_consensus="split",
                    confidence="low",
                    disclaimer=DISCLAIMER,
                )
            return BoardOpinion(
                primary_diagnosis=f"종합 분석 실패: {e}",
                primary_diagnosis_evidence="",
                key_findings=[f"오류: {e}"],
                agent_opinions=agent_opinions,
                agent_consensus="split",
                confidence="low",
                disclaimer=DISCLAIMER,
            )

        if mode == "cancer":
            return self._parse_cancer_response(response, agent_opinions)
        return self._parse_response(response, agent_opinions)

    def _format_curated_evidence(self, curated_treatments: dict) -> str:
        """Render the curator output as a prompt-friendly, LLM-citable block.

        Every row lists ``curated_id``, ``variant_key``, drug, evidence level,
        disease context and PMIDs so the LLM has exactly one authoritative
        source to cite from. Empty or missing curator output is rendered as a
        "no curated rows available" placeholder so the prompt stays stable.
        """
        if not curated_treatments:
            return "(no curated rows available — do not invent drug/target pairings)"
        lines = []
        for variant_key, rows in curated_treatments.items():
            if not rows:
                lines.append(f"- variant_key={variant_key}: (no curated rows)")
                continue
            for row in rows:
                drug = getattr(row, "drug", "") or (row.get("drug") if isinstance(row, dict) else "")
                cid = getattr(row, "curated_id", "") or (row.get("curated_id") if isinstance(row, dict) else "")
                level = getattr(row, "evidence_level", "") or (
                    row.get("evidence_level") if isinstance(row, dict) else ""
                )
                disease = getattr(row, "disease_context", "") or (
                    row.get("disease_context") if isinstance(row, dict) else ""
                )
                pmids = getattr(row, "pmids", []) or (row.get("pmids") if isinstance(row, dict) else [])
                source = getattr(row, "source", "") or (row.get("source") if isinstance(row, dict) else "")
                pmid_str = ",".join(str(p) for p in (pmids or []))
                lines.append(
                    f"- curated_id={cid} variant_key={variant_key} drug={drug} "
                    f"level={level} source={source} disease={disease} pmids={pmid_str}"
                )
        return "\n".join(lines)

    def _build_prompt(
        self,
        case_briefing: str,
        agent_opinions: List[AgentOpinion],
        curated_treatments: dict = None,
    ) -> str:
        """Build the synthesis prompt with case briefing and all agent opinions."""
        opinions_text = "\n\n".join(_format_agent_opinion(op, self.language) for op in agent_opinions)
        curated_block = self._format_curated_evidence(curated_treatments)

        if self.language == "ko":
            return f"""다음은 임상유전학 사례 회의 자료입니다.
케이스 정보와 4명의 전문의 소견을 종합하여 최종 진단 의견을 제시하세요.
반드시 한국어로 응답하세요.

## 케이스 정보
{case_briefing}

## CURATED EVIDENCE (authoritative — 치료 옵션의 유일한 출처)
{curated_block}

## 전문의 소견

{opinions_text}

## 요청
위 케이스 정보와 전문의 소견을 종합하여, 지정된 JSON 형식으로 최종 진단 의견을 제시하세요.
**치료 옵션은 반드시 CURATED EVIDENCE 섹션의 curated_id만 사용하여 작성하시오.**"""
        else:
            return f"""The following is a clinical genetics case conference.
Synthesize the case information and 4 domain specialists' opinions into a final diagnostic opinion.
Respond in English.

## Case Information
{case_briefing}

## CURATED EVIDENCE (authoritative — the only source for treatment options)
{curated_block}

## Specialist Opinions

{opinions_text}

## Request
Synthesize the case information and specialist opinions above into a final diagnostic opinion in the specified JSON format.
**Treatment options MUST cite only curated_id values from the CURATED EVIDENCE section above.**"""

    def _parse_response(self, response, agent_opinions: List[AgentOpinion]) -> BoardOpinion:
        """Parse LLM response into BoardOpinion."""
        if isinstance(response, str):
            try:
                response = json.loads(response)
            except json.JSONDecodeError:
                return BoardOpinion(
                    primary_diagnosis=response[:200],
                    primary_diagnosis_evidence="",
                    key_findings=[response],
                    agent_opinions=agent_opinions,
                    agent_consensus="split",
                    confidence="low",
                    disclaimer=DISCLAIMER,
                    raw_response=response,
                )

        return BoardOpinion(
            primary_diagnosis=response.get("primary_diagnosis", ""),
            primary_diagnosis_evidence=response.get("primary_diagnosis_evidence", ""),
            differential_diagnoses=response.get("differential_diagnoses", []),
            key_findings=response.get("key_findings", []),
            recommendations=response.get("recommendations", []),
            agent_opinions=agent_opinions,
            agent_consensus=response.get("agent_consensus", "split"),
            dissenting_opinions=response.get("dissenting_opinions", []),
            follow_up=response.get("follow_up", []),
            confidence=response.get("confidence", "moderate"),
            disclaimer=DISCLAIMER,
            raw_response=(json.dumps(response, ensure_ascii=False) if isinstance(response, dict) else str(response)),
        )

    def _parse_cancer_response(self, response, agent_opinions: List[AgentOpinion]) -> CancerBoardOpinion:
        """Parse LLM response into CancerBoardOpinion (treatment-focused)."""
        if isinstance(response, str):
            try:
                response = json.loads(response)
            except json.JSONDecodeError:
                return CancerBoardOpinion(
                    therapeutic_implications=response[:200],
                    agent_opinions=agent_opinions,
                    agent_consensus="split",
                    confidence="low",
                    disclaimer=DISCLAIMER,
                    raw_response=response,
                )

        # Preserve curated_id + variant_key on every treatment row so the
        # narrative_scrubber downstream can enforce the (cid, variant_key)
        # pair binding.
        treatment_options = []
        for opt in response.get("treatment_options", []) or []:
            if isinstance(opt, dict):
                treatment_options.append(
                    {
                        "drug": opt.get("drug", ""),
                        "curated_id": opt.get("curated_id", ""),
                        "variant_key": opt.get("variant_key", ""),
                        "evidence_level": opt.get("evidence_level", ""),
                        "resistance_notes": opt.get("resistance_notes", ""),
                    }
                )

        return CancerBoardOpinion(
            therapeutic_headline=response.get("therapeutic_headline", ""),
            therapeutic_implications=response.get("therapeutic_implications", ""),
            therapeutic_evidence=response.get("therapeutic_evidence", ""),
            treatment_options=treatment_options,
            actionable_findings=response.get("actionable_findings", []),
            clinical_actions=response.get("clinical_actions", []),
            immunotherapy_eligibility=response.get("immunotherapy_eligibility", ""),
            agent_opinions=agent_opinions,
            agent_consensus=response.get("agent_consensus", "split"),
            dissenting_opinions=response.get("dissenting_opinions", []),
            monitoring_plan=response.get("monitoring_plan", []),
            confidence=response.get("confidence", "moderate"),
            disclaimer=DISCLAIMER,
            raw_response=(json.dumps(response, ensure_ascii=False) if isinstance(response, dict) else str(response)),
        )
