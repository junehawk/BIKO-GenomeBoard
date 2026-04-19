"""NCBI Gene E-utilities — fetch gene summaries from NCBI Gene database."""

import logging
import time
from typing import Dict, Optional

from scripts.common.api_utils import fetch_with_retry
from scripts.common.config import get

logger = logging.getLogger(__name__)

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def _get_api_key() -> str:
    return get("api.ncbi_api_key", "") or ""


def fetch_gene_summary(gene_symbol: str) -> Optional[Dict]:
    """Fetch gene summary from NCBI Gene database.

    Returns: {
        "gene": "TP53",
        "entrez_id": "7157",
        "full_name": "tumor protein p53",
        "summary": "This gene encodes...",
        "aliases": "p53, LFS1",
    } or None on failure.
    """
    api_key = _get_api_key()
    key_param = f"&api_key={api_key}" if api_key else ""

    # Step 1: Search for gene ID
    search_url = f"{ESEARCH}?db=gene&term={gene_symbol}[Gene Name]+AND+Homo+sapiens[Organism]&retmode=json{key_param}"
    search_data = fetch_with_retry(search_url)
    if not search_data:
        return None

    id_list = search_data.get("esearchresult", {}).get("idlist", [])
    if not id_list:
        return None

    gene_id = id_list[0]

    # Rate limit: 3/sec without key, 10/sec with key
    time.sleep(0.35 if not api_key else 0.1)

    # Step 2: Fetch summary
    summary_url = f"{ESUMMARY}?db=gene&id={gene_id}&retmode=json{key_param}"
    summary_data = fetch_with_retry(summary_url)
    if not summary_data:
        return None

    result = summary_data.get("result", {}).get(gene_id)
    if not result:
        return None

    return {
        "gene": gene_symbol,
        "entrez_id": gene_id,
        "full_name": result.get("description", ""),
        "summary": result.get("summary", ""),
        "aliases": result.get("otheraliases", ""),
    }
