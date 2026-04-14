"""Template-renderer chair (AI Clinical Board v2.2 · A2 fallback).

Deterministic fallback that emits a ``CancerBoardOpinion`` directly from
curator rows when the Board Chair LLM response contains zero valid
``(curated_id, variant_key)`` pairs. Keeps the clinician-facing therapy
table populated even when the LLM synthesiser mis-formats its JSON or
hallucinates the whole block.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from scripts.clinical_board.models import (
    BOARD_DISCLAIMER,
    AgentOpinion,
    CancerBoardOpinion,
)


def _row_get(row: Any, name: str, default: Any = "") -> Any:
    if isinstance(row, dict):
        return row.get(name, default)
    return getattr(row, name, default)


def render_from_curated(
    curated_by_variant: Dict[str, List[Any]],
    agent_opinions: Optional[List[AgentOpinion]] = None,
) -> CancerBoardOpinion:
    """Build a ``CancerBoardOpinion`` from curator rows without any LLM call.

    The emitted rows carry a ``_hydrated`` sub-dict so ``render.py`` can
    display drug/target/PMIDs/disease_context without re-walking the curator
    output. Confidence is set to ``low`` and consensus to ``fallback`` so
    downstream readers can tell the opinion came from the deterministic
    path.
    """
    treatment_rows: List[Dict[str, Any]] = []
    for variant_key, rows in (curated_by_variant or {}).items():
        for row in rows or []:
            drug = _row_get(row, "drug", "")
            if not drug:
                continue
            treatment_rows.append(
                {
                    "drug": drug,
                    "curated_id": _row_get(row, "curated_id", ""),
                    "variant_key": variant_key,
                    "evidence_level": _row_get(row, "evidence_level", "D"),
                    "resistance_notes": "",
                    "_hydrated": {
                        "drug": drug,
                        "target": _row_get(row, "target", ""),
                        "pmids": list(_row_get(row, "pmids", []) or []),
                        "disease_context": _row_get(row, "disease_context", ""),
                        "source": _row_get(row, "source", ""),
                        "significance": _row_get(row, "significance", "sensitivity"),
                    },
                }
            )

    # Preserve AMP level ranking (A > B > C > D) in the fallback table
    _rank = {"A": 0, "B": 1, "C": 2, "D": 3}
    treatment_rows.sort(key=lambda r: (_rank.get(r.get("evidence_level", "D"), 9), r.get("drug", "").lower()))

    if treatment_rows:
        headline = f"{len(treatment_rows)} curated therapy option(s) — LLM synthesis unavailable"
        implications = (
            "Deterministic fallback: the Board Chair LLM response did not "
            "produce any valid (curated_id, variant_key) pair. The treatment "
            "options below are drawn directly from the OncoKB/CIViC curator "
            "output. No LLM narrative is included."
        )
    else:
        headline = "No curated therapy options"
        implications = (
            "No curated therapy rows are available for the selected variants "
            "and the Board Chair synthesis failed validation."
        )

    return CancerBoardOpinion(
        therapeutic_headline=headline,
        therapeutic_implications=implications,
        therapeutic_evidence="Curated rows rendered without LLM narrative.",
        treatment_options=treatment_rows,
        actionable_findings=[],
        clinical_actions=[],
        immunotherapy_eligibility="",
        agent_opinions=list(agent_opinions or []),
        agent_consensus="fallback",
        dissenting_opinions=[],
        monitoring_plan=[],
        confidence="low",
        disclaimer=BOARD_DISCLAIMER,
    )


__all__ = ["render_from_curated"]
