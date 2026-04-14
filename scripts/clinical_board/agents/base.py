"""Base class for Clinical Board domain agents."""

import json
import logging
from abc import ABC, abstractmethod

from scripts.clinical_board.models import AgentOpinion
from scripts.clinical_board.ollama_client import OllamaClient
from scripts.common.config import get

logger = logging.getLogger(__name__)


class BaseAgent(ABC):
    """Abstract base for domain-specific clinical agents."""

    def __init__(self, client: OllamaClient = None, model: str = None, language: str = None):
        self.client = client or OllamaClient()
        self.model = model or get("clinical_board.agent_model", "medgemma:27b")
        self.language = language or get("clinical_board.language", "en")

    @property
    @abstractmethod
    def agent_name(self) -> str: ...

    @property
    @abstractmethod
    def domain(self) -> str: ...

    @property
    @abstractmethod
    def system_prompt(self) -> str: ...

    def analyze(
        self,
        case_briefing: str,
        domain_sheet: str = "",
        prior_knowledge: str = "",
    ) -> AgentOpinion:
        """Run analysis on a case briefing. Returns structured opinion."""
        prompt = self._build_prompt(case_briefing, domain_sheet, prior_knowledge)
        temperature = get("clinical_board.temperature", 0.1)

        try:
            response = self.client.generate_json(
                model=self.model,
                prompt=prompt,
                system=self.system_prompt,
                temperature=temperature,
            )
        except Exception as e:
            logger.error(f"{self.agent_name} analysis failed: {e}")
            return AgentOpinion(
                agent_name=self.agent_name,
                domain=self.domain,
                findings=[
                    {
                        "finding": f"Analysis failed: {e}",
                        "evidence": "",
                        "confidence": "low",
                    }
                ],
                confidence="low",
            )

        return self._parse_response(response)

    def _build_prompt(
        self,
        case_briefing: str,
        domain_sheet: str = "",
        prior_knowledge: str = "",
    ) -> str:
        """Build the analysis prompt with case briefing and output format instructions."""
        domain_section = f"\n\n## Domain-Specific Data\n{domain_sheet}" if domain_sheet else ""
        prior_section = f"\n\n## Prior Board Knowledge\n{prior_knowledge}" if prior_knowledge else ""

        if self.language == "ko":
            return f"""다음 환자의 유전체 분석 결과를 검토하고, 당신의 전문 분야 관점에서 소견을 제시하세요.

## 케이스 정보
{case_briefing}{domain_section}{prior_section}

Note: Clinical notes may be provided in Korean or English. Interpret accordingly.

## 응답 형식 (JSON)
반드시 한국어로 응답하세요. 다음 JSON 형식으로 응답하세요:
{{
  "findings": [
    {{"finding": "소견 내용", "evidence": "근거", "confidence": "high/moderate/low"}}
  ],
  "recommendations": ["권고사항 1", "권고사항 2"],
  "concerns": ["우려사항 (있을 경우)"],
  "references": ["PMID:xxxxx 등 참고문헌"],
  "confidence": "high/moderate/low"
}}"""
        else:
            return f"""Review the following genomic analysis results and provide your expert opinion from your domain perspective.

## Case Information
{case_briefing}{domain_section}{prior_section}

Note: Clinical notes may be provided in Korean or English. Interpret accordingly.

## Response Format (JSON)
Respond in English. Use the following JSON format:
{{
  "findings": [
    {{"finding": "finding description", "evidence": "supporting evidence", "confidence": "high/moderate/low"}}
  ],
  "recommendations": ["recommendation 1", "recommendation 2"],
  "concerns": ["concerns if any"],
  "references": ["PMID:xxxxx etc."],
  "confidence": "high/moderate/low"
}}"""

    def _parse_response(self, response) -> AgentOpinion:
        """Parse LLM response (dict or str) into AgentOpinion."""
        if isinstance(response, str):
            try:
                response = json.loads(response)
            except json.JSONDecodeError:
                return AgentOpinion(
                    agent_name=self.agent_name,
                    domain=self.domain,
                    findings=[
                        {"finding": response, "evidence": "", "confidence": "low"}
                    ],
                    confidence="low",
                    raw_response=response,
                )

        return AgentOpinion(
            agent_name=self.agent_name,
            domain=self.domain,
            findings=response.get("findings", []),
            recommendations=response.get("recommendations", []),
            concerns=response.get("concerns", []),
            references=response.get("references", []),
            confidence=response.get("confidence", "moderate"),
            raw_response=(
                json.dumps(response, ensure_ascii=False)
                if isinstance(response, dict)
                else str(response)
            ),
        )
