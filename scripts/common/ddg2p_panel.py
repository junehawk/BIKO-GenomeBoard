"""DDG2P neurodevelopmental gene panel loader.

Reads ``data/ddg2p_neurodev_genes.json`` (curated subset in v1, full DDG2P in
v1.1 via ``scripts/tools/build_ddg2p_table.py``) and exposes two helpers used
by the variant selector's rare-disease de novo carve-out:

- :func:`is_admitted_neurodev_gene` — is this gene in the admission set
  (definitive/strong/moderate confidence per spec Q2)?
- :func:`get_neurodev_info` — return the full record or ``None``.

Load failures (file missing, malformed JSON) are logged **once per process**
and the module degrades to returning ``False`` / ``None`` for every query.
This mirrors the v2.2 T3 ClinGen hardening pattern so a fresh CI runner without
the optional data file does not break the rare-disease mode — it just silently
disables the de novo neurodev carve-out.
"""

from __future__ import annotations

import json
import logging
import os
import threading
from typing import Any, Dict, Optional, Set

from scripts.common.availability_cache import AvailabilityCache

logger = logging.getLogger(__name__)

_PANEL_PATH_DEFAULT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data",
    "ddg2p_neurodev_genes.json",
)

_PANEL_CACHE: Optional[Dict[str, Any]] = None
_ADMITTED_CACHE: Optional[Set[str]] = None
_LOAD_LOCK = threading.Lock()

# Log-once memoisation for the "panel missing / malformed" warnings.
# Distinct from _PANEL_CACHE (which caches the parsed data) so that
# _reset_for_tests can clear the warning state without rebuilding the
# panel on unrelated tests.
_LOAD_AVAILABILITY = AvailabilityCache()


def _load_panel() -> Dict[str, Any]:
    """Return the panel dict, loading it once per process (thread-safe)."""
    global _PANEL_CACHE, _ADMITTED_CACHE

    if _PANEL_CACHE is not None:
        return _PANEL_CACHE

    with _LOAD_LOCK:
        if _PANEL_CACHE is not None:
            return _PANEL_CACHE

        path = _PANEL_PATH_DEFAULT

        def _probe(p: str) -> bool:
            if not os.path.exists(p):
                logger.warning(
                    "[ddg2p_panel] panel file not found at %s — de novo neurodev "
                    "carve-out will be disabled for this process",
                    p,
                )
                return False
            return True

        if not _LOAD_AVAILABILITY.check(path, _probe):
            _PANEL_CACHE = {"genes": {}, "admission_confidences": []}
            _ADMITTED_CACHE = set()
            return _PANEL_CACHE

        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
        except (OSError, json.JSONDecodeError) as exc:
            logger.warning("[ddg2p_panel] panel load failed: %s", exc)
            _PANEL_CACHE = {"genes": {}, "admission_confidences": []}
            _ADMITTED_CACHE = set()
            return _PANEL_CACHE

        _PANEL_CACHE = data
        admitted_confs = {c.lower() for c in (data.get("admission_confidences") or [])}
        genes = data.get("genes") or {}
        _ADMITTED_CACHE = {
            sym
            for sym, rec in genes.items()
            if isinstance(rec, dict) and (rec.get("confidence") or "").lower() in admitted_confs
        }
        return _PANEL_CACHE


def is_admitted_neurodev_gene(gene: Optional[str]) -> bool:
    """Return True if ``gene`` is in the DDG2P admission set (defin/strong/mod).

    Returns False for ``None``, empty strings, and any gene outside the
    admission set. Load failures degrade to False for every gene.
    """
    if not gene:
        return False
    _load_panel()
    return gene in (_ADMITTED_CACHE or set())


def get_neurodev_info(gene: Optional[str]) -> Optional[Dict[str, Any]]:
    """Return the per-gene DDG2P record (``confidence``, ``disease``,
    ``allelic_requirement``) or ``None`` if the gene is not in the panel."""
    if not gene:
        return None
    panel = _load_panel()
    genes = panel.get("genes") or {}
    rec = genes.get(gene)
    if isinstance(rec, dict):
        return rec
    return None


def _reset_for_tests() -> None:
    """Test helper: clear the module-level cache so tests can force a reload."""
    global _PANEL_CACHE, _ADMITTED_CACHE
    with _LOAD_LOCK:
        _PANEL_CACHE = None
        _ADMITTED_CACHE = None
        _LOAD_AVAILABILITY.reset()
