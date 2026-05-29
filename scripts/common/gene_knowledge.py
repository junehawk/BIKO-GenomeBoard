import json
import logging
from pathlib import Path
from typing import Dict, Optional

from scripts.common.config import get

logger = logging.getLogger(__name__)

_KNOWLEDGE = None


def _load_knowledge() -> dict:
    global _KNOWLEDGE
    if _KNOWLEDGE is None:
        path = get("paths.gene_knowledge") or str(Path(__file__).parent.parent.parent / "data" / "gene_knowledge.json")
        try:
            with open(path) as f:
                data = json.load(f)
            # Degrade gracefully on a structurally-valid-but-wrong file (missing
            # "genes" key, or entries lacking a "gene" key) instead of raising a
            # KeyError that would crash the whole pipeline (v2.7 review ENRI-04).
            _KNOWLEDGE = {g["gene"]: g for g in data.get("genes", []) if isinstance(g, dict) and g.get("gene")}
        except (FileNotFoundError, json.JSONDecodeError, TypeError, AttributeError) as exc:
            logger.warning("gene_knowledge load failed (%s): %s", path, exc)
            _KNOWLEDGE = {}
    return _KNOWLEDGE


def get_gene_info(gene: str) -> Optional[Dict]:
    """Get clinical knowledge for a gene."""
    return _load_knowledge().get(gene)
