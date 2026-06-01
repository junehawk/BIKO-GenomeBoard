"""Base class for Clinical Board domain agents."""

import json
import logging
from abc import ABC, abstractmethod

from scripts.clinical_board.models import AgentOpinion
from scripts.clinical_board.ollama_client import OllamaClient
from scripts.common.config import get

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Shared scientific-grounding / anti-hallucination clause (v2.8).
#
# Injected into EVERY domain agent's user prompt by ``_build_prompt`` so the
# guardrail can never drift between agents or between the KO/EN variants. Before
# v2.8 only the Literature Analyst (and, for drugs, the Clinical Evidence
# Analyst) had an anti-fabrication clause; the genomics/genetics/pathology
# agents could recall pathway membership, driver status, gene-disease strength,
# allele frequencies, PMIDs, and trial names from training memory with no fence.
# Keep the two language variants structurally parallel — a test asserts this.
# ─────────────────────────────────────────────────────────────────────────────

GROUNDING_CLAUSE_EN = """## Scientific Grounding (MANDATORY — applies to every finding, recommendation, and concern)
Base every factual claim ONLY on data explicitly present in the case information above
(the variant table, in-silico scores, HPO terms, the CURATED EVIDENCE block, the CIViC
evidence section, and the gene_knowledge section). You may interpret and reason over that
data, but you may NOT introduce a fact that is not in the briefing. This specifically
forbids inventing or recalling from training memory: pathway/network membership,
driver-vs-passenger status, gene-disease association strength or inheritance mode, allele
frequencies, PMIDs or citations, clinical-trial / NCT identifiers, drug names, and
metabolizer phenotypes. If a claim is not supported by a specific briefing item, either
omit it or write exactly "not determinable from provided data". The "evidence" field of
each finding must name the briefing element it came from (e.g. curated_id, CIViC EID,
"in-silico: REVEL=…", "HPO: …", gene_knowledge); a finding whose evidence cannot name a
source must be dropped. Do not restate, re-derive, or contradict the deterministic
ACMG/AMP classification or tiering — it is fixed; you may only interpret it."""

GROUNDING_CLAUSE_KO = """## 과학적 근거 (필수 — 모든 finding·recommendation·concern에 적용)
모든 사실 주장은 위 케이스 정보(변이표, in-silico 점수, HPO 용어, CURATED EVIDENCE 블록,
CIViC 근거 섹션, gene_knowledge 섹션)에 명시된 데이터에만 근거하시오. 그 데이터를 해석·추론할
수는 있으나, briefing에 없는 사실을 새로 만들어내지 마시오. 특히 학습 기억에서 다음을 지어내거나
회상하지 마시오: pathway/네트워크 소속, driver-vs-passenger 여부, gene-disease 연관 강도 또는
유전 양식, 대립유전자 빈도, PMID/인용, 임상시험/NCT 식별자, 약물명, metabolizer 표현형.
특정 briefing 항목으로 뒷받침되지 않는 주장은 생략하거나 정확히 "not determinable from
provided data"라고 적으시오. 각 finding의 "evidence" 필드는 출처가 된 briefing 요소를
명시해야 합니다 (예: curated_id, CIViC EID, "in-silico: REVEL=…", "HPO: …", gene_knowledge).
출처를 댈 수 없는 finding은 버리시오. 결정적 ACMG/AMP 분류·tiering을 재진술·재도출·반박하지
마시오 — 고정된 값이며, 해석만 가능합니다."""


def grounding_clause(language: str) -> str:
    """Return the scientific-grounding clause for the given language."""
    return GROUNDING_CLAUSE_KO if language == "ko" else GROUNDING_CLAUSE_EN


def _format_opinion_for_peer(opinion: AgentOpinion, language: str = "en") -> str:
    """Compact, readable digest of one agent opinion for the deliberation round."""
    conf = getattr(opinion, "confidence", "") or ""
    if language == "ko":
        lines = [f"=== {opinion.agent_name} ({opinion.domain}) · 신뢰도: {conf} ==="]
        label_f, label_r, label_c = "소견", "권고", "우려"
    else:
        lines = [f"=== {opinion.agent_name} ({opinion.domain}) · confidence: {conf} ==="]
        label_f, label_r, label_c = "Findings", "Recommendations", "Concerns"
    findings = getattr(opinion, "findings", None) or []
    if findings:
        lines.append(f"{label_f}:")
        for f in findings:
            if isinstance(f, dict):
                txt = f.get("finding", "")
                ev = f.get("evidence", "")
                fc = f.get("confidence", "")
                suffix = f" [{fc}]" if fc else ""
                lines.append(f"  - {txt}{suffix}" + (f" (evidence: {ev})" if ev else ""))
            else:
                lines.append(f"  - {f}")
    recs = getattr(opinion, "recommendations", None) or []
    if recs:
        lines.append(f"{label_r}:")
        for r in recs:
            lines.append(f"  - {r}")
    concerns = getattr(opinion, "concerns", None) or []
    if concerns:
        lines.append(f"{label_c}:")
        for c in concerns:
            lines.append(f"  - {c}")
    return "\n".join(lines)


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

{grounding_clause("ko")}

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

{grounding_clause("en")}

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

    # ── Round 2: deliberation / peer cross-review (v2.8) ─────────────────────
    def revise(
        self,
        case_briefing: str,
        own_opinion: AgentOpinion,
        peer_opinions: list,
        domain_sheet: str = "",
        prior_knowledge: str = "",
    ) -> AgentOpinion:
        """Re-examine this agent's Round-1 opinion in light of the peers' opinions.

        The agent may: DEFER where a peer owns the domain, RETRACT a claim a peer
        credibly refuted, RAISE/LOWER its confidence on convergence/divergence, and
        FLAG genuine contradictions in ``concerns``. The scientific-grounding
        contract still holds — no new fabricated facts may be introduced during
        revision. On any failure (or a degenerate empty revision) the original
        Round-1 opinion is returned unchanged, so the deliberation round can never
        lose signal relative to a single round.
        """
        prompt = self._build_revision_prompt(case_briefing, own_opinion, peer_opinions, domain_sheet, prior_knowledge)
        temperature = get("clinical_board.temperature", 0.1)
        try:
            response = self.client.generate_json(
                model=self.model,
                prompt=prompt,
                system=self.system_prompt,
                temperature=temperature,
            )
        except Exception as e:
            logger.error(f"{self.agent_name} revision failed: {e}")
            return own_opinion

        revised = self._parse_response(response)
        # Guard against a degenerate revision wiping out a useful Round-1 result.
        if not revised.findings and not revised.recommendations:
            return own_opinion
        return revised

    def _build_revision_prompt(
        self,
        case_briefing: str,
        own_opinion: AgentOpinion,
        peer_opinions: list,
        domain_sheet: str = "",
        prior_knowledge: str = "",
    ) -> str:
        """Build the Round-2 revision prompt (own opinion + peer digest + rules)."""
        own = _format_opinion_for_peer(own_opinion, self.language)
        peers = "\n\n".join(_format_opinion_for_peer(op, self.language) for op in peer_opinions) or "(none)"
        domain_section = f"\n\n## Domain-Specific Data\n{domain_sheet}" if domain_sheet else ""
        prior_section = f"\n\n## Prior Board Knowledge\n{prior_knowledge}" if prior_knowledge else ""

        if self.language == "ko":
            return f"""당신은 1차 소견을 이미 제출했습니다. 이제 동료 전문의들의 1차 소견을 읽고
당신의 소견을 **재검토**하여 최종 소견을 제출하세요.

## 케이스 정보
{case_briefing}{domain_section}{prior_section}

## 당신의 1차 소견
{own}

## 동료 전문의 1차 소견
{peers}

## 재검토 지침
- 동료가 더 전문성 있는 영역(예: 당신의 영역이 아닌 부분)은 그 동료에게 **양보(defer)**하고
  중복 소견을 제거하시오.
- 동료가 신뢰할 만하게 반박한 주장은 **철회(retract)**하시오.
- 의견이 수렴하면 confidence를 높이고, 갈리면 낮추시오.
- 실제 모순이 있으면 `concerns`에 명시하시오.
- **새로운 사실을 만들어내지 마시오.** 위 케이스 정보(briefing)에 있는 데이터로만 근거하시오.
  (1차 때와 동일한 grounding 계약)
- 변경이 없다면 1차 소견을 그대로 다시 제출해도 됩니다.

{grounding_clause("ko")}

## 응답 형식 (JSON)
반드시 한국어로, 1차와 동일한 JSON 형식(findings/recommendations/concerns/references/confidence)으로
최종(수정) 소견을 제출하세요."""
        return f"""You have already submitted a Round-1 opinion. Now read your peers' Round-1
opinions and **revise** your own into a final opinion.

## Case Information
{case_briefing}{domain_section}{prior_section}

## Your Round-1 Opinion
{own}

## Peer Round-1 Opinions
{peers}

## Revision Guidance
- DEFER to a peer who owns a domain outside your own and remove duplicate findings.
- RETRACT any claim a peer credibly refuted.
- RAISE confidence where opinions converge; LOWER it where they diverge.
- FLAG genuine contradictions in `concerns`.
- Do NOT introduce any new fact — ground every claim only in the case information above
  (same grounding contract as Round 1).
- If nothing changes, you may resubmit your Round-1 opinion unchanged.

{grounding_clause("en")}

## Response Format (JSON)
Respond in English in the SAME JSON format as Round 1
(findings/recommendations/concerns/references/confidence) with your final (revised) opinion."""

    def _parse_response(self, response) -> AgentOpinion:
        """Parse LLM response (dict or str) into AgentOpinion."""
        if isinstance(response, str):
            try:
                response = json.loads(response)
            except json.JSONDecodeError:
                return AgentOpinion(
                    agent_name=self.agent_name,
                    domain=self.domain,
                    findings=[{"finding": response, "evidence": "", "confidence": "low"}],
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
            raw_response=(json.dumps(response, ensure_ascii=False) if isinstance(response, dict) else str(response)),
        )
