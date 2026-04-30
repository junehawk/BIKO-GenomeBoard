"""Narrative scrubber (AI Clinical Board v2.2 · A2, v2.5.5 lookup backfill).

Post-processor that walks the full ``CancerBoardOpinion`` dataclass and:

1. Validates each ``treatment_options`` row against the curator and drops any
   row whose ``curated_id`` is not in the curator output. The ``variant_key``
   is **deterministically backfilled** from the curator dict — the scrubber
   does not trust the LLM's variant_key field. This pins each row to the
   variant the curator originally associated with that curated_id, closing
   the v2.2 paste-attack threat model: pasting a valid curated_id under a
   different variant in the LLM narrative is a no-op (the row is reassigned
   to the original variant).
2. Scrubs banned drug tokens out of every prose field — headline, body,
   evidence paragraph, action lists, monitoring plan, agent-opinion findings
   and recommendations. "Banned" = a drug name that appeared in a dropped
   treatment row AND is not present in the curator's allowed-drugs set.

The scrubber is the patient-safety gate enforced *after* the LLM response is
parsed and *before* it is serialised to ``raw_opinion_json`` or rendered.

v2.5.5 (2026-04-30): switched from ``(curated_id, variant_key)`` pair-matching
to ``curated_id`` lookup with variant_key backfill. Motivation: supergemma4-31b
miscopied variant_key on long Cancer Chair prompts (e.g., ``chr13:32356550:C:T``
→ ``chr13:32356 own:C:T``), causing 25 % drop rate on real demo runs. Pair
matching is mathematically equivalent to lookup-then-backfill, but the latter
is robust to BPE/SentencePiece-level token confusion in the LLM's variant_key
copy. Paste-attack defence is preserved (and arguably strengthened): a pasted
curated_id is reassigned to its original variant_key, never to a different
variant.
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


def _curated_lookup(curated_by_variant: Dict[str, list]) -> Dict[str, Tuple[str, Any]]:
    """Map ``curated_id`` → ``(variant_key, original_curated_row)``.

    The curator deduplicates curated_ids per variant; cross-variant collisions
    are not expected (curated_ids are high-entropy hex). On a hypothetical
    collision, last-wins.
    """
    lookup: Dict[str, Tuple[str, Any]] = {}
    for variant_key, rows in (curated_by_variant or {}).items():
        for row in rows or []:
            cid = _row_attr(row, "curated_id")
            if cid:
                lookup[cid] = (variant_key, row)
    return lookup


def _valid_pair_set(curated_by_variant: Dict[str, list]) -> Set[Tuple[str, str]]:
    """Set of ``(curated_id, variant_key)`` pairs the curator emitted.

    Retained for backward compatibility with callers that still want the
    pair-set view. The scrubber itself uses ``_curated_lookup`` so that
    variant_key is treated as derived data, not a verification key.
    """
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


def validate_treatment_option(option: Dict[str, Any], curated_by_variant: Dict[str, list]) -> bool:
    """Return True iff the row's ``curated_id`` is present in the curator output.

    A row is valid iff it carries a non-empty ``curated_id`` and that id is
    present in the curator's lookup. ``variant_key`` is ignored at validation
    time — :func:`scrub_opinion` overwrites it with the curator's authoritative
    value before keeping the row.
    """
    cid = str(option.get("curated_id", "") or "").strip()
    if not cid:
        return False
    return cid in _curated_lookup(curated_by_variant)


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
    """Walk ``opinion`` in place: drop bad rows, backfill variant_key, scrub banned tokens.

    Args:
        opinion: a ``CancerBoardOpinion`` (or ``BoardOpinion``) instance —
            mutated in place.
        curated_by_variant: output of ``curate_treatments()``.

    Returns:
        Stats dict ``{"kept": int, "dropped": int, "banned_terms": int}``
        for caller-side logging.
    """
    allowed_drugs = _allowed_drug_set(curated_by_variant)
    lookup = _curated_lookup(curated_by_variant)

    banned: Set[str] = set()
    kept: List[Dict[str, Any]] = []
    dropped = 0

    rows = list(getattr(opinion, "treatment_options", []) or [])
    for row in rows:
        if not isinstance(row, dict):
            continue
        cid = str(row.get("curated_id", "") or "").strip()
        drug = str(row.get("drug", "") or "").strip()
        llm_vk = str(row.get("variant_key", "") or "").strip()

        if cid in lookup and drug.lower() in allowed_drugs:
            true_vk, _ = lookup[cid]
            if llm_vk and llm_vk != true_vk:
                logger.info(
                    "[narrative_scrubber] backfilled variant_key for curated_id=%r: "
                    "LLM emitted %r, curator says %r — using curator value",
                    cid,
                    llm_vk,
                    true_vk,
                )
            row["variant_key"] = true_vk
            kept.append(row)
            continue

        dropped += 1
        if drug and drug.lower() not in allowed_drugs:
            banned.add(drug.lower())
        logger.warning(
            "[narrative_scrubber] dropped treatment row drug=%r curated_id=%r variant_key=%r",
            drug,
            cid,
            llm_vk,
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
                new_list = [_scrub_text(item, banned) if isinstance(item, str) else item for item in val]
                setattr(opinion, field, new_list)

        # Agent opinions (nested dataclasses)
        for agent_op in getattr(opinion, "agent_opinions", None) or []:
            new_findings = []
            for f in getattr(agent_op, "findings", None) or []:
                if isinstance(f, dict):
                    new_findings.append(
                        {k: (_scrub_text(v, banned) if isinstance(v, str) else v) for k, v in f.items()}
                    )
                else:
                    new_findings.append(f)
            agent_op.findings = new_findings
            agent_op.recommendations = [
                _scrub_text(r, banned) if isinstance(r, str) else r
                for r in (getattr(agent_op, "recommendations", None) or [])
            ]
            agent_op.concerns = [
                _scrub_text(c, banned) if isinstance(c, str) else c for c in (getattr(agent_op, "concerns", None) or [])
            ]
            if isinstance(getattr(agent_op, "raw_response", None), str):
                agent_op.raw_response = _scrub_text(agent_op.raw_response, banned)

    return {
        "kept": len(kept),
        "dropped": dropped,
        "banned_terms": len(banned),
    }


__all__ = ["scrub_opinion", "validate_treatment_option"]
