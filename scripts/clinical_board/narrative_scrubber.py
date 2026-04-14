"""Narrative scrubber (AI Clinical Board v2.2 · A2).

Post-processor that walks the full ``CancerBoardOpinion`` dataclass and:

1. Drops any ``treatment_options`` row whose ``(curated_id, variant_key)``
   pair was not emitted by the curator (catches copy-paste attacks where an
   LLM reuses a valid curated_id under the wrong variant).
2. Scrubs banned drug tokens out of every prose field — headline, body,
   evidence paragraph, action lists, monitoring plan, agent-opinion findings
   and recommendations. "Banned" = a drug name that appeared in a dropped
   treatment row AND is not present in the curator's allowed-drugs set.

The scrubber is the patient-safety gate enforced *after* the LLM response is
parsed and *before* it is serialised to ``raw_opinion_json`` or rendered.
"""
from __future__ import annotations

import logging
import re
from typing import Any, Dict, Iterable, List, Set, Tuple

logger = logging.getLogger(__name__)


def _allowed_drug_set(curated_by_variant: Dict[str, list]) -> Set[str]:
    """Lowercased union of every drug name the curator emitted."""
    allowed: Set[str] = set()
    for rows in (curated_by_variant or {}).values():
        for row in rows or []:
            drug = _row_attr(row, "drug")
            if drug:
                allowed.add(drug.strip().lower())
    return allowed


def _valid_pair_set(curated_by_variant: Dict[str, list]) -> Set[Tuple[str, str]]:
    """Set of ``(curated_id, variant_key)`` pairs the curator emitted."""
    pairs: Set[Tuple[str, str]] = set()
    for variant_key, rows in (curated_by_variant or {}).items():
        for row in rows or []:
            cid = _row_attr(row, "curated_id")
            if cid:
                pairs.add((cid, variant_key))
    return pairs


def _row_attr(row: Any, name: str) -> str:
    """Read an attribute from a ``CuratedTreatment`` dataclass *or* a dict."""
    if isinstance(row, dict):
        return str(row.get(name, "") or "")
    return str(getattr(row, name, "") or "")


def _scrub_text(text: str, banned: Iterable[str]) -> str:
    """Remove every banned token (case-insensitive, word-boundary)."""
    if not text or not banned:
        return text
    out = text
    for token in banned:
        if not token:
            continue
        pattern = re.compile(rf"\b{re.escape(token)}\b", re.IGNORECASE)
        out = pattern.sub("[REDACTED-DRUG]", out)
    out = re.sub(r"\s{2,}", " ", out).strip()
    return out


def validate_treatment_option(
    option: Dict[str, Any], curated_by_variant: Dict[str, list]
) -> bool:
    """Return True iff the row's ``(curated_id, variant_key)`` pair is valid.

    A row is valid iff:
    * it carries a non-empty ``curated_id`` AND a non-empty ``variant_key``
    * that exact pair is present in the curator output (the curator is
      authoritative — if it didn't emit the pair, the row is hallucinated)
    """
    cid = str(option.get("curated_id", "") or "")
    vk = str(option.get("variant_key", "") or "")
    if not cid or not vk:
        return False
    return (cid, vk) in _valid_pair_set(curated_by_variant)


_PROSE_STR_FIELDS = (
    "therapeutic_headline",
    "therapeutic_implications",
    "therapeutic_evidence",
    "immunotherapy_eligibility",
    "primary_diagnosis",
    "primary_diagnosis_evidence",
    "raw_response",
)

_PROSE_LIST_FIELDS = (
    "actionable_findings",
    "clinical_actions",
    "monitoring_plan",
    "dissenting_opinions",
    "key_findings",
    "recommendations",
    "follow_up",
)


def scrub_opinion(opinion: Any, curated_by_variant: Dict[str, list]) -> Dict[str, int]:
    """Walk ``opinion`` in place: drop bad rows, scrub banned tokens.

    Args:
        opinion: a ``CancerBoardOpinion`` (or ``BoardOpinion``) instance —
            mutated in place.
        curated_by_variant: output of ``curate_treatments()``.

    Returns:
        Stats dict ``{"kept": int, "dropped": int, "banned_terms": int}``
        for caller-side logging.
    """
    allowed_drugs = _allowed_drug_set(curated_by_variant)
    valid_pairs = _valid_pair_set(curated_by_variant)

    banned: Set[str] = set()
    kept: List[Dict[str, Any]] = []
    dropped = 0

    rows = list(getattr(opinion, "treatment_options", []) or [])
    for row in rows:
        if not isinstance(row, dict):
            continue
        cid = str(row.get("curated_id", "") or "")
        vk = str(row.get("variant_key", "") or "")
        drug = str(row.get("drug", "") or "").strip()
        if (cid, vk) in valid_pairs and drug.lower() in allowed_drugs:
            kept.append(row)
            continue
        dropped += 1
        if drug and drug.lower() not in allowed_drugs:
            banned.add(drug.lower())
        logger.warning(
            "[narrative_scrubber] dropped treatment row drug=%r curated_id=%r variant_key=%r",
            drug, cid, vk,
        )

    if hasattr(opinion, "treatment_options"):
        opinion.treatment_options = kept

    # Any dropped row means the raw LLM response is now a patient-safety
    # liability — it may contain the hallucinated drug mention (fixture 1)
    # or the paste-attacked curated_id (fixture 2) that no longer survives
    # in the validated opinion. Clear raw_response so it cannot leak back
    # into ``json.dumps(dataclasses.asdict(opinion))``.
    if dropped and hasattr(opinion, "raw_response"):
        opinion.raw_response = ""

    # Scrub prose fields
    if banned:
        for field in _PROSE_STR_FIELDS:
            val = getattr(opinion, field, None)
            if isinstance(val, str):
                setattr(opinion, field, _scrub_text(val, banned))
        for field in _PROSE_LIST_FIELDS:
            val = getattr(opinion, field, None)
            if isinstance(val, list):
                new_list = [
                    _scrub_text(item, banned) if isinstance(item, str) else item
                    for item in val
                ]
                setattr(opinion, field, new_list)

        # Agent opinions (nested dataclasses)
        for agent_op in getattr(opinion, "agent_opinions", None) or []:
            new_findings = []
            for f in getattr(agent_op, "findings", None) or []:
                if isinstance(f, dict):
                    new_findings.append({
                        k: (_scrub_text(v, banned) if isinstance(v, str) else v)
                        for k, v in f.items()
                    })
                else:
                    new_findings.append(f)
            agent_op.findings = new_findings
            agent_op.recommendations = [
                _scrub_text(r, banned) if isinstance(r, str) else r
                for r in (getattr(agent_op, "recommendations", None) or [])
            ]
            agent_op.concerns = [
                _scrub_text(c, banned) if isinstance(c, str) else c
                for c in (getattr(agent_op, "concerns", None) or [])
            ]
            if isinstance(getattr(agent_op, "raw_response", None), str):
                agent_op.raw_response = _scrub_text(agent_op.raw_response, banned)

    return {
        "kept": len(kept),
        "dropped": dropped,
        "banned_terms": len(banned),
    }


__all__ = ["scrub_opinion", "validate_treatment_option"]
