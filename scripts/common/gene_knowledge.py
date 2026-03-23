import json
from pathlib import Path
from typing import Optional, Dict
from scripts.common.config import get

_KNOWLEDGE = None


def _load_knowledge() -> dict:
    global _KNOWLEDGE
    if _KNOWLEDGE is None:
        path = get("paths.gene_knowledge") or str(Path(__file__).parent.parent.parent / "data" / "gene_knowledge.json")
        try:
            with open(path) as f:
                data = json.load(f)
            _KNOWLEDGE = {g["gene"]: g for g in data["genes"]}
        except (FileNotFoundError, json.JSONDecodeError):
            _KNOWLEDGE = {}
    return _KNOWLEDGE


def get_gene_info(gene: str) -> Optional[Dict]:
    """Get clinical knowledge for a gene."""
    return _load_knowledge().get(gene)
