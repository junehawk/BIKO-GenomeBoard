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


def assign_tier(classification: str, gene: str, clinvar_significance: str = "") -> int:
    """Assign reporting tier based on classification + OncoKB gene status.

    Tier 1: Pathogenic/LP with therapeutic implications (OncoKB Level 1-2)
    Tier 2: Pathogenic/LP on any cancer gene (OncoKB Level 3+)
    Tier 3: VUS/Benign but on OncoKB cancer gene → abbreviated report
    Tier 4: VUS/Benign on non-cancer gene → count only, no report

    Returns: 1, 2, 3, or 4
    """
    cls_lower = classification.lower()
    gene_info = get_cancer_gene_info(gene) if gene else None

    # Drug Response and Risk Factor always reported
    if cls_lower in ("drug response", "risk factor"):
        return 1

    # Pathogenic / Likely Pathogenic
    if "pathogenic" in cls_lower and "benign" not in cls_lower:
        if gene_info and gene_info.get("level") in ("1", "2"):
            return 1  # High-level therapeutic target
        return 2  # Pathogenic on any gene

    # VUS on cancer gene
    if cls_lower == "vus" and gene_info:
        return 3  # Cancer gene VUS → abbreviated report

    # Benign on cancer gene (still worth noting)
    if "benign" in cls_lower and gene_info:
        return 4  # Cancer gene but benign → count only

    # Everything else
    return 4  # No report detail


def get_tier_label(tier: int) -> str:
    """Human-readable tier label."""
    labels = {
        1: "Tier 1 — Therapeutic Target",
        2: "Tier 2 — Clinically Significant",
        3: "Tier 3 — Cancer Gene (VUS)",
        4: "Tier 4 — No Clinical Significance",
    }
    return labels.get(tier, "Unknown")
