"""OncoKB free-tier variant API client (curate-then-narrate infrastructure).

Thin wrapper over ``/annotate/mutations/byProteinChange``. Responses are
cached via ``scripts.common.cache.get_cached_ns``/``set_cached_ns`` in the
``oncokb`` namespace so repeated calls for the same (gene, alteration)
don't re-hit the upstream API.

Any transient network failure (HTTP 429, 5xx, connection-refused, timeout)
raises :class:`OncoKBUnavailable` — the curator catches this and degrades
to CIViC-only rather than letting a network blip fail the whole board run.
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

import requests

from scripts.common.cache import get_cached_ns, set_cached_ns
from scripts.common.config import get

logger = logging.getLogger(__name__)

_NAMESPACE = "oncokb"


class OncoKBUnavailable(Exception):
    """Raised when the OncoKB free-tier API is unreachable or rate-limited."""


def _base_url() -> str:
    return get("clinical_board.curated_treatments.oncokb_base_url",
               "https://www.oncokb.org/api/v1")


def _timeout_s() -> float:
    return float(get("clinical_board.curated_treatments.http_timeout_s", 5.0))


def _normalise_level(raw: Optional[str]) -> str:
    """OncoKB returns ``LEVEL_1``/``LEVEL_2A``/``LEVEL_R1`` etc. Normalise
    to AMP-style single letters (``A``/``B``/``C``/``D``) for the curator's
    ranking layer. Unknown levels fall through as 'D' (weakest)."""
    if not raw:
        return "D"
    upper = str(raw).upper()
    if "LEVEL_1" in upper or upper == "1":
        return "A"
    if "LEVEL_2" in upper or upper.startswith("2"):
        return "B"
    if "LEVEL_3" in upper or upper.startswith("3"):
        return "C"
    if "LEVEL_4" in upper or upper.startswith("4"):
        return "D"
    if "LEVEL_R1" in upper or upper == "R1":
        return "A"  # FDA-recognised resistance
    if "LEVEL_R2" in upper or upper == "R2":
        return "C"
    return "D"


def _parse_response(payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Flatten OncoKB response to a list of dicts the curator can consume.

    Each dict carries: ``drug``, ``level``, ``pmids``, ``disease``,
    ``significance``, ``therapy_ids`` (empty string — OncoKB free tier
    doesn't expose therapy IDs), ``raw_row`` (the original payload chunk
    for audit).
    """
    out: List[Dict[str, Any]] = []
    treatments = payload.get("treatments") or []
    for t in treatments:
        drugs = t.get("drugs") or []
        drug_names = ", ".join(
            d.get("drugName") or d.get("name") or ""
            for d in drugs if isinstance(d, dict)
        ).strip(", ")
        if not drug_names and isinstance(t.get("drug"), str):
            drug_names = t["drug"]  # Some free-tier shims emit flat {drug: ...}
        if not drug_names:
            continue

        pmids_raw = t.get("pmids") or []
        if isinstance(pmids_raw, str):
            pmids = [p.strip() for p in pmids_raw.split(",") if p.strip()]
        else:
            pmids = [str(p) for p in pmids_raw if p]

        level = _normalise_level(t.get("level") or t.get("levelOfEvidence"))
        out.append({
            "drug": drug_names,
            "level": level,
            "pmids": pmids,
            "disease": t.get("indication") or t.get("tumorType") or "",
            "significance": "sensitivity" if "R" not in str(t.get("level", "")).upper() else "resistance",
            "therapy_ids": "",  # free tier doesn't expose stable therapy IDs
            "raw_row": t,
        })
    return out


def annotate_protein_change(
    gene: str,
    alteration: str,
    offline_mode: bool = False,
) -> List[Dict[str, Any]]:
    """Return OncoKB curated treatment rows for ``(gene, alteration)``.

    Args:
        gene: HUGO symbol (e.g. ``"KRAS"``).
        alteration: protein change in either CIViC-style (``"G12D"``) or
            HGVSp short form. The upstream endpoint accepts both.
        offline_mode: when ``True``, return ``[]`` without touching the
            network. Used by the curator's config-gated offline path.

    Returns:
        A list of dicts as produced by :func:`_parse_response`. Empty list
        on a successful query that returned no hits.

    Raises:
        OncoKBUnavailable: on 429 / 5xx / ConnectionError / Timeout.
    """
    if offline_mode:
        return []
    if not gene or not alteration:
        return []

    cache_key = f"{gene}::{alteration}"
    cached = get_cached_ns(_NAMESPACE, cache_key)
    if cached is not None:
        return cached

    url = f"{_base_url().rstrip('/')}/annotate/mutations/byProteinChange"
    params = {"hugoSymbol": gene, "alteration": alteration}
    try:
        response = requests.get(url, params=params, timeout=_timeout_s())
    except requests.exceptions.Timeout as e:
        raise OncoKBUnavailable(f"timeout querying OncoKB for {gene} {alteration}: {e}") from e
    except requests.exceptions.ConnectionError as e:
        raise OncoKBUnavailable(f"connection refused by OncoKB: {e}") from e
    except requests.exceptions.RequestException as e:
        raise OncoKBUnavailable(f"OncoKB request failed: {e}") from e

    status = getattr(response, "status_code", None)
    if status == 429:
        raise OncoKBUnavailable(f"OncoKB rate-limit (HTTP 429) for {gene} {alteration}")
    if status is not None and status >= 500:
        raise OncoKBUnavailable(f"OncoKB HTTP {status} for {gene} {alteration}")
    if status is not None and status >= 400:
        # 4xx other than 429 — treat as "no hit" rather than crash the run
        logger.warning(
            "[oncokb_client] HTTP %s for %s %s — treating as empty", status, gene, alteration
        )
        set_cached_ns(_NAMESPACE, cache_key, [])
        return []

    try:
        payload = response.json() if callable(getattr(response, "json", None)) else {}
    except Exception as e:  # pragma: no cover — defensive
        raise OncoKBUnavailable(f"OncoKB response not JSON: {e}") from e

    parsed = _parse_response(payload or {})
    set_cached_ns(_NAMESPACE, cache_key, parsed)
    return parsed
