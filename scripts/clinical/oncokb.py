"""OncoKB cancer gene lookup and variant tiering."""

import json
import logging
from pathlib import Path
from typing import Dict, Optional

from scripts.common.config import get

logger = logging.getLogger(__name__)

_ONCOKB_DATA = None


def _load_oncokb() -> dict:
    global _ONCOKB_DATA
    if _ONCOKB_DATA is None:
        path = get(
            "paths.oncokb_genes",
            str(Path(__file__).parent.parent.parent / "data" / "oncokb_cancer_genes.json"),
        )
        try:
            with open(path) as f:
                data = json.load(f)
            _ONCOKB_DATA = data.get("genes", {})
        except (FileNotFoundError, json.JSONDecodeError):
            logger.warning(f"OncoKB gene list not found: {path}")
            _ONCOKB_DATA = {}
    return _ONCOKB_DATA


def reset_oncokb_cache() -> None:
    """Reset cached OncoKB data (useful for testing)."""
    global _ONCOKB_DATA
    _ONCOKB_DATA = None


def is_cancer_gene(gene: str) -> bool:
    """Check if gene is in OncoKB cancer gene list."""
    return gene in _load_oncokb()


def get_cancer_gene_info(gene: str) -> Optional[Dict]:
    """Get OncoKB info for a cancer gene. Returns None if not a cancer gene."""
    data = _load_oncokb()
    if gene in data:
        return {"gene": gene, **data[gene]}
    return None


def assign_tier(classification: str, gene: str, clinvar_significance: str = "", hgvsp: str = "") -> int:
    """DEPRECATED: Use scripts.somatic.amp_tiering.amp_assign_tier() instead.
    Kept for backward compatibility — delegates to amp_assign_tier with strategy C.
    """
    from scripts.somatic.amp_tiering import amp_assign_tier

    result = amp_assign_tier(classification, gene, hgvsp=hgvsp, strategy="C")
    return result.tier


def get_tier_label(tier: int) -> str:
    """Human-readable tier label (AMP/ASCO/CAP 2017)."""
    from scripts.somatic.amp_tiering import AMP_TIER_LABELS

    return AMP_TIER_LABELS.get(tier, "Unknown")
