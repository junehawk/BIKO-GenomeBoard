"""HPO phenotype matching for rare disease variant prioritization."""
import logging
from typing import Dict, List, Optional
from scripts.common.api_utils import fetch_with_retry

logger = logging.getLogger(__name__)

HPO_API = "https://ontology.jax.org/api/hp"


def resolve_hpo_terms(hpo_ids: List[str]) -> List[Dict]:
    """Resolve HPO IDs to names and associated genes.
    Returns: [{"id": "HP:0001250", "name": "Seizures", "genes": ["SCN1A", "KCNQ2", ...]}]
    """
    results = []
    for hpo_id in hpo_ids:
        hpo_id = hpo_id.strip()
        if not hpo_id.startswith("HP:"):
            continue

        # Get term name
        term_data = fetch_with_retry(f"{HPO_API}/terms/{hpo_id}")
        name = ""
        if term_data:
            name = term_data.get("name", hpo_id)

        # Get associated genes
        genes_data = fetch_with_retry(f"{HPO_API}/terms/{hpo_id}/genes")
        genes = []
        if genes_data and isinstance(genes_data, dict):
            for gene_entry in genes_data.get("genes", []):
                gene_symbol = gene_entry.get("symbol") or gene_entry.get("name", "")
                if gene_symbol:
                    genes.append(gene_symbol)
        elif genes_data and isinstance(genes_data, list):
            for gene_entry in genes_data:
                gene_symbol = gene_entry.get("symbol") or gene_entry.get("name", "")
                if gene_symbol:
                    genes.append(gene_symbol)

        results.append({"id": hpo_id, "name": name or hpo_id, "genes": genes})

    return results


def calculate_hpo_score(gene: str, hpo_results: List[Dict]) -> int:
    """Calculate HPO overlap score for a gene. +1 per HPO term whose gene list includes this gene."""
    score = 0
    for hpo in hpo_results:
        if gene and gene.upper() in [g.upper() for g in hpo.get("genes", [])]:
            score += 1
    return score


def get_matching_hpo_terms(gene: str, hpo_results: List[Dict]) -> List[str]:
    """Get list of HPO term names that match a gene."""
    matching = []
    for hpo in hpo_results:
        if gene and gene.upper() in [g.upper() for g in hpo.get("genes", [])]:
            matching.append(f"{hpo['id']} ({hpo['name']})")
    return matching
