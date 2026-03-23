"""ClinGen gene-disease validity lookup."""
import logging
from typing import Optional
from scripts.common.api_utils import fetch_with_retry

logger = logging.getLogger(__name__)

CLINGEN_API = "https://search.clinicalgenome.org/kb/gene-validity"

# Static data for demo genes (ClinGen search API can be unreliable)
CLINGEN_VALIDITY = {
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
    Returns: Definitive, Strong, Moderate, Limited, Disputed, Refuted, or None
    """
    return CLINGEN_VALIDITY.get(gene)
