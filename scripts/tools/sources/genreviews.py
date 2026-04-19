"""GeneReviews PMID lookup via PubMed E-utilities."""

import logging
import time
from typing import Dict, Optional

from scripts.common.api_utils import fetch_with_retry
from scripts.common.config import get

logger = logging.getLogger(__name__)

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def fetch_genreviews_info(gene_symbol: str) -> Optional[Dict]:
    """Search PubMed for GeneReviews entry for a gene.

    Returns: {
        "gene": "TP53",
        "pmid": "20301371",
        "title": "Li-Fraumeni Syndrome",
        "nbk_id": "NBK1311",
        "url": "https://www.ncbi.nlm.nih.gov/books/NBK1311/",
    } or None.
    """
    api_key = get("api.ncbi_api_key", "") or ""
    key_param = f"&api_key={api_key}" if api_key else ""

    # Search PubMed for GeneReviews entries mentioning this gene
    search_url = f"{ESEARCH}?db=pubmed&term={gene_symbol}[Title]+AND+GeneReviews[Book]&retmode=json&retmax=1{key_param}"
    search_data = fetch_with_retry(search_url)
    if not search_data:
        return None

    id_list = search_data.get("esearchresult", {}).get("idlist", [])
    if not id_list:
        return None

    pmid = id_list[0]

    time.sleep(0.35 if not api_key else 0.1)

    # Fetch article summary
    summary_url = f"{ESUMMARY}?db=pubmed&id={pmid}&retmode=json{key_param}"
    summary_data = fetch_with_retry(summary_url)
    if not summary_data:
        return {"gene": gene_symbol, "pmid": pmid, "title": "", "nbk_id": "", "url": ""}

    article = summary_data.get("result", {}).get(pmid, {})
    title = article.get("title", "")
    nbk_id = article.get("bookshelfaccession", "")

    # Try to extract NBK from elocationid if bookshelfaccession not present
    if not nbk_id:
        eloc = article.get("elocationid", "")
        if "NBK" in eloc:
            nbk_id = eloc.split("NBK")[-1].split(".")[0]
            nbk_id = f"NBK{nbk_id}"

    url = f"https://www.ncbi.nlm.nih.gov/books/{nbk_id}/" if nbk_id else ""

    return {
        "gene": gene_symbol,
        "pmid": pmid,
        "title": title,
        "nbk_id": nbk_id,
        "url": url,
    }
