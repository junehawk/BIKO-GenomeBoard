"""Data models for Clinical Board system."""
from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class AgentOpinion:
    """Single domain agent's analysis opinion."""
    agent_name: str           # e.g., "Variant Pathologist"
    domain: str               # e.g., "variant_pathology"
    findings: List[dict] = field(default_factory=list)      # [{finding: str, evidence: str, confidence: str}]
    recommendations: List[str] = field(default_factory=list)
    concerns: List[str] = field(default_factory=list)
    references: List[str] = field(default_factory=list)     # PMIDs or citations
    confidence: str = "moderate"  # high/moderate/low
    raw_response: str = ""        # Original LLM response for debugging


@dataclass
class BoardOpinion:
    """Synthesized opinion from the full Clinical Board."""
    primary_diagnosis: str = ""
    primary_diagnosis_evidence: str = ""
    differential_diagnoses: List[dict] = field(default_factory=list)  # [{diagnosis, likelihood, evidence}]
    key_findings: List[str] = field(default_factory=list)
    recommendations: List[str] = field(default_factory=list)
    agent_opinions: List[AgentOpinion] = field(default_factory=list)
    agent_consensus: str = "unknown"  # unanimous/majority/split
    dissenting_opinions: List[str] = field(default_factory=list)
    follow_up: List[str] = field(default_factory=list)
    confidence: str = "moderate"
    disclaimer: str = (
        "[AI-Generated] 본 진단 종합 의견은 로컬 LLM에 의해 생성되었으며, "
        "임상의의 검토와 확인이 반드시 필요합니다. 이 의견은 최종 진단이 아니며, "
        "임상적 의사결정의 보조 자료로만 활용되어야 합니다."
    )
    raw_response: str = ""
