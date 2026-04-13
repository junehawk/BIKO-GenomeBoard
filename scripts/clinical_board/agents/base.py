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

    def __init__(self, client: OllamaClient = None, model: str = None):
        self.client = client or OllamaClient()
        self.model = model or get("clinical_board.agent_model", "medgemma:27b")

    @property
    @abstractmethod
    def agent_name(self) -> str: ...

    @property
    @abstractmethod
    def domain(self) -> str: ...

    @property
    @abstractmethod
    def system_prompt(self) -> str: ...

    def analyze(self, case_briefing: str) -> AgentOpinion:
        """Run analysis on a case briefing. Returns structured opinion."""
        prompt = self._build_prompt(case_briefing)

        try:
            response = self.client.generate_json(
                model=self.model,
                prompt=prompt,
                system=self.system_prompt,
                temperature=0.3,
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

    def _build_prompt(self, case_briefing: str) -> str:
        """Build the analysis prompt with case briefing and output format instructions."""
        return f"""다음 환자의 유전체 분석 결과를 검토하고, 당신의 전문 분야 관점에서 소견을 제시하세요.

## 케이스 정보
{case_briefing}

## 응답 형식 (JSON)
다음 JSON 형식으로 응답하세요:
{{
  "findings": [
    {{"finding": "소견 내용", "evidence": "근거", "confidence": "high/moderate/low"}}
  ],
  "recommendations": ["권고사항 1", "권고사항 2"],
  "concerns": ["우려사항 (있을 경우)"],
  "references": ["PMID:xxxxx 등 참고문헌"],
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
