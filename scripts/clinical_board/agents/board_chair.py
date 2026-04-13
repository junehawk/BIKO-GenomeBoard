"""Board Chair — synthesizes all agent opinions into a unified diagnostic opinion."""

import json
import logging
from typing import List

from scripts.clinical_board.models import AgentOpinion, BoardOpinion
from scripts.clinical_board.ollama_client import OllamaClient
from scripts.common.config import get

logger = logging.getLogger(__name__)

DISCLAIMER = (
    "[AI-Generated] 본 소견은 로컬 LLM 기반 다전문가 시스템에 의해 생성되었습니다. "
    "최종 임상 판단은 반드시 담당 임상의가 수행해야 합니다. "
    "변이 분류(ACMG/AMP)는 결정적 분류 엔진에 의해 수행되었으며, "
    "본 소견은 해당 분류를 변경하지 않습니다."
)

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
        self, case_briefing: str, agent_opinions: List[AgentOpinion]
    ) -> BoardOpinion:
        """Synthesize all agent opinions into a unified board opinion."""
        prompt = self._build_prompt(case_briefing, agent_opinions)

        try:
            response = self.client.generate_json(
                model=self.model,
                prompt=prompt,
                system=SYSTEM_PROMPT,
                temperature=0.3,
            )
        except Exception as e:
            logger.error(f"Board Chair synthesis failed: {e}")
            return BoardOpinion(
                primary_diagnosis=f"종합 분석 실패: {e}",
                primary_diagnosis_evidence="",
                key_findings=[f"오류: {e}"],
                agent_opinions=agent_opinions,
                agent_consensus="split",
                confidence="low",
                disclaimer=DISCLAIMER,
            )

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
