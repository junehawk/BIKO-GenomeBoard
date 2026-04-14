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

SYSTEM_PROMPT = """당신은 임상유전학 사례 회의의 위원장(Board Chair)입니다.
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


CANCER_SYSTEM_PROMPT = """당신은 종양유전체 임상 회의의 위원장(Cancer Board Chair)입니다.
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

## 응답 형식 (JSON)
{
  "therapeutic_headline": "Stage IV PDAC — KRAS G12D driver, no standard targeted therapy",
  "therapeutic_implications": "치료 시사점 종합 (상세 단락)",
  "therapeutic_evidence": "근거 요약 (CIViC level 등)",
  "treatment_options": [
    {"drug": "약물명", "evidence_level": "A/B/C/D", "resistance_notes": "저항성 위험 (있을 경우)"}
  ],
  "actionable_findings": ["임상적으로 활용 가능한 핵심 소견"],
  "clinical_actions": ["권고되는 임상 조치"],
  "immunotherapy_eligibility": "면역치료 적격성 평가",
  "monitoring_plan": ["모니터링 항목"],
  "agent_consensus": "unanimous/majority/split",
  "dissenting_opinions": ["소수 의견"],
  "confidence": "high/moderate/low"
}"""


def _format_agent_opinion(opinion: AgentOpinion) -> str:
    """Format an agent opinion into readable text for the chair prompt."""
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
        mode: str = "rare-disease",
    ) -> Union[BoardOpinion, CancerBoardOpinion]:
        """Synthesize all agent opinions into a unified board opinion.

        Args:
            mode: "rare-disease" (default) returns BoardOpinion;
                  "cancer" returns CancerBoardOpinion (treatment-focused).
        """
        prompt = self._build_prompt(case_briefing, agent_opinions)
        system_prompt = CANCER_SYSTEM_PROMPT if mode == "cancer" else SYSTEM_PROMPT

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

    def _build_prompt(
        self, case_briefing: str, agent_opinions: List[AgentOpinion]
    ) -> str:
        """Build the synthesis prompt with case briefing and all agent opinions."""
        opinions_text = "\n\n".join(
            _format_agent_opinion(op) for op in agent_opinions
        )

        if self.language == "ko":
            return f"""다음은 임상유전학 사례 회의 자료입니다.
케이스 정보와 4명의 전문의 소견을 종합하여 최종 진단 의견을 제시하세요.
반드시 한국어로 응답하세요.

## 케이스 정보
{case_briefing}

## 전문의 소견

{opinions_text}

## 요청
위 케이스 정보와 전문의 소견을 종합하여, 지정된 JSON 형식으로 최종 진단 의견을 제시하세요."""
        else:
            return f"""The following is a clinical genetics case conference.
Synthesize the case information and 4 domain specialists' opinions into a final diagnostic opinion.
Respond in English.

## Case Information
{case_briefing}

## Specialist Opinions

{opinions_text}

## Request
Synthesize the case information and specialist opinions above into a final diagnostic opinion in the specified JSON format."""

    def _parse_response(
        self, response, agent_opinions: List[AgentOpinion]
    ) -> BoardOpinion:
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
            primary_diagnosis_evidence=response.get(
                "primary_diagnosis_evidence", ""
            ),
            differential_diagnoses=response.get("differential_diagnoses", []),
            key_findings=response.get("key_findings", []),
            recommendations=response.get("recommendations", []),
            agent_opinions=agent_opinions,
            agent_consensus=response.get("agent_consensus", "split"),
            dissenting_opinions=response.get("dissenting_opinions", []),
            follow_up=response.get("follow_up", []),
            confidence=response.get("confidence", "moderate"),
            disclaimer=DISCLAIMER,
            raw_response=(
                json.dumps(response, ensure_ascii=False)
                if isinstance(response, dict)
                else str(response)
            ),
        )

    def _parse_cancer_response(
        self, response, agent_opinions: List[AgentOpinion]
    ) -> CancerBoardOpinion:
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

        return CancerBoardOpinion(
            therapeutic_headline=response.get("therapeutic_headline", ""),
            therapeutic_implications=response.get("therapeutic_implications", ""),
            therapeutic_evidence=response.get("therapeutic_evidence", ""),
            treatment_options=response.get("treatment_options", []),
            actionable_findings=response.get("actionable_findings", []),
            clinical_actions=response.get("clinical_actions", []),
            immunotherapy_eligibility=response.get("immunotherapy_eligibility", ""),
            agent_opinions=agent_opinions,
            agent_consensus=response.get("agent_consensus", "split"),
            dissenting_opinions=response.get("dissenting_opinions", []),
            monitoring_plan=response.get("monitoring_plan", []),
            confidence=response.get("confidence", "moderate"),
            disclaimer=DISCLAIMER,
            raw_response=(
                json.dumps(response, ensure_ascii=False)
                if isinstance(response, dict)
                else str(response)
            ),
        )
