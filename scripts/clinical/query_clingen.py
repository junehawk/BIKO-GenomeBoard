"""ClinGen gene-disease validity lookup — local DB with static fallback."""
import logging
from typing import Optional

logger = logging.getLogger(__name__)

# Minimal static fallback (used only when local DB unavailable)
_STATIC_FALLBACK = {
    "TP53": "Definitive",
    "BRCA2": "Definitive",
    "CFTR": "Definitive",
    "ATM": "Definitive",
    "MUTYH": "Definitive",
    "PALB2": "Definitive",
    "PTPN11": "Definitive",
}


def get_gene_validity(gene: str) -> Optional[str]:
    """Get ClinGen gene-disease validity classification.
    Strategy: local DB first → static fallback.
    """
    try:
        from scripts.db.query_local_clingen import get_gene_validity_local
        result = get_gene_validity_local(gene)
        if result:
            return result
    except Exception as e:
        logger.debug(f"ClinGen local DB not available: {e}")

    return _STATIC_FALLBACK.get(gene)
