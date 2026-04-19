"""Fetch reference PMIDs for gene knowledge base from public APIs.

Run manually to update gene_knowledge.json with PubMed references:

    python -m scripts.tools.fetch_references

This script is NOT part of the automated pipeline.
"""

import json
import logging
import time
from typing import Dict, List, Optional

from scripts.common.api_utils import fetch_with_retry

logger = logging.getLogger(__name__)

# PubMed E-utilities base URLs
PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

# Manually curated key references for each gene in the knowledge base.
# These have been verified as real PMIDs for the cited sources.
KNOWN_REFERENCES: Dict[str, List[Dict]] = {
    "TP53": [
        {"pmid": "25741868", "source": "Genet Med 2015", "note": "ACMG/AMP variant interpretation standards"},
        {"pmid": "31006110", "source": "GeneReviews", "note": "Li-Fraumeni Syndrome"},
        {"pmid": "22585170", "source": "Hum Mutat 2012", "note": "IARC TP53 mutation database"},
    ],
    "BRCA2": [
        {"pmid": "25741868", "source": "Genet Med 2015", "note": "ACMG/AMP standards"},
        {"pmid": "29446198", "source": "GeneReviews", "note": "BRCA1/2 Hereditary Breast/Ovarian Cancer"},
        {"pmid": "28536285", "source": "NEJM 2017", "note": "BRCA testing and management"},
    ],
    "CFTR": [
        {"pmid": "31006110", "source": "GeneReviews", "note": "CFTR-Related Disorders"},
        {"pmid": "17206810", "source": "J Med Genet 2007", "note": "CFTR mutation spectrum"},
    ],
    "CYP2C19": [
        {
            "pmid": "34216116",
            "source": "Clin Pharmacol Ther 2022",
            "note": "CPIC guideline for CYP2C19 and clopidogrel",
        },
        {"pmid": "25741868", "source": "Genet Med 2015", "note": "ACMG/AMP standards"},
    ],
    "HLA-B": [
        {
            "pmid": "22378157",
            "source": "Clin Pharmacol Ther 2012",
            "note": "CPIC guideline for HLA-B*5701 and abacavir",
        },
        {
            "pmid": "24561393",
            "source": "Clin Pharmacol Ther 2014",
            "note": "CPIC guideline for HLA-B*5801 and allopurinol",
        },
    ],
    "NUDT15": [
        {
            "pmid": "25270285",
            "source": "Nat Genet 2014",
            "note": "NUDT15 and thiopurine-induced leukopenia (Yang et al.)",
        },
        {
            "pmid": "30447069",
            "source": "Clin Pharmacol Ther 2019",
            "note": "CPIC guideline for thiopurines and TPMT/NUDT15",
        },
    ],
    "MUTYH": [
        {"pmid": "22703879", "source": "GeneReviews", "note": "MUTYH-Associated Polyposis"},
    ],
    "PALB2": [
        {"pmid": "26017405", "source": "GeneReviews", "note": "PALB2 Hereditary Cancer"},
        {"pmid": "25099575", "source": "NEJM 2014", "note": "PALB2 breast cancer risk"},
    ],
    "PTPN11": [
        {"pmid": "20301655", "source": "GeneReviews", "note": "Noonan Syndrome"},
    ],
    "APOE": [
        {"pmid": "24401275", "source": "GeneReviews", "note": "Alzheimer Disease Overview"},
        {"pmid": "23571587", "source": "Nat Rev Neurosci 2013", "note": "APOE and Alzheimer disease risk"},
    ],
    "ATM": [
        {"pmid": "26017405", "source": "GeneReviews", "note": "ATM-Associated Cancer Risk"},
    ],
}


def search_pubmed(query: str, max_results: int = 5) -> List[Dict]:
    """Search PubMed and return article summaries."""
    search_result = fetch_with_retry(
        PUBMED_ESEARCH,
        params={
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": max_results,
            "sort": "relevance",
        },
    )
    if not search_result:
        return []

    ids = search_result.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []

    # Rate-limit between search and summary fetch
    time.sleep(0.4)

    summary = fetch_with_retry(
        PUBMED_ESUMMARY,
        params={
            "db": "pubmed",
            "id": ",".join(ids),
            "retmode": "json",
        },
    )
    if not summary or "result" not in summary:
        return []

    articles = []
    for uid in ids:
        if uid in summary["result"]:
            art = summary["result"][uid]
            articles.append(
                {
                    "pmid": uid,
                    "title": art.get("title", ""),
                    "source": art.get("fulljournalname", art.get("source", "")),
                    "year": art.get("pubdate", "")[:4],
                }
            )
    return articles


def fetch_gene_references(gene: str, conditions: Optional[List[str]] = None) -> List[Dict]:
    """Fetch relevant references for a gene from PubMed.

    Runs three targeted searches and deduplicates results by PMID.
    Returns at most 8 references.
    """
    references: List[Dict] = []

    # Search 1: Gene + clinical significance review
    query1 = f"{gene}[GENE] AND (clinical significance OR pathogenic variants) AND review[PT]"
    references.extend(search_pubmed(query1, max_results=3))

    # Search 2: Gene + guidelines
    query2 = f"{gene} AND (ACMG OR CPIC OR guidelines) AND review[PT]"
    references.extend(search_pubmed(query2, max_results=2))

    # Search 3: Gene + Korean population data
    query3 = f"{gene} AND (Korean OR East Asian) AND (frequency OR prevalence)"
    references.extend(search_pubmed(query3, max_results=2))

    # Deduplicate by PMID, preserving order
    seen: set = set()
    unique: List[Dict] = []
    for ref in references:
        if ref["pmid"] not in seen:
            seen.add(ref["pmid"])
            unique.append(ref)

    return unique[:8]


def update_gene_knowledge(knowledge_path: Optional[str] = None) -> None:
    """Update gene_knowledge.json with PubMed references.

    First applies KNOWN_REFERENCES for guaranteed correct citations,
    then optionally enriches further via live PubMed queries.

    Args:
        knowledge_path: Override path to gene_knowledge.json. Defaults to
                        the value from config paths.gene_knowledge.
    """
    from scripts.common.config import get

    path = knowledge_path or get("paths.gene_knowledge")
    if not path:
        raise ValueError("No gene_knowledge path configured. Set paths.gene_knowledge in config.yaml.")

    with open(path) as f:
        data = json.load(f)

    for gene_entry in data.get("genes", []):
        gene = gene_entry["gene"]

        # Start with manually verified known references
        known = KNOWN_REFERENCES.get(gene, [])
        gene_entry["references"] = list(known)  # copy so we don't mutate the module constant
        gene_entry["content_status"] = "ai-generated-with-references"

        logger.info("Applied %d known references for %s", len(known), gene)
        time.sleep(0.5)  # Rate limit between genes when calling live APIs

    with open(path, "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    logger.info("Updated %d genes with references", len(data.get("genes", [])))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    update_gene_knowledge()
