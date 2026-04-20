"""HPO phenotype matching for rare disease variant prioritization."""

import logging
from typing import Dict, List

from scripts.common.api_utils import fetch_with_retry

logger = logging.getLogger(__name__)

HPO_API = "https://ontology.jax.org/api/hp"


def resolve_hpo_terms(hpo_ids: List[str]) -> List[Dict]:
    """Resolve HPO IDs to names and associated genes.
    Strategy: API first → local DB fallback if API fails or returns no genes.
    """
    results = []
    for hpo_id in hpo_ids:
        hpo_id = hpo_id.strip()
        if not hpo_id.startswith("HP:"):
            continue

        name = ""
        genes = []

        # Try API first
        term_data = fetch_with_retry(f"{HPO_API}/terms/{hpo_id}")
        if term_data:
            name = term_data.get("name", "")

        genes_data = fetch_with_retry(f"{HPO_API}/terms/{hpo_id}/genes")
        if genes_data and isinstance(genes_data, dict):
            for ge in genes_data.get("genes", []):
                sym = ge.get("symbol") or ge.get("name", "")
                if sym:
                    genes.append(sym)
        elif genes_data and isinstance(genes_data, list):
            for ge in genes_data:
                sym = ge.get("symbol") or ge.get("name", "")
                if sym:
                    genes.append(sym)

        # Fallback to local DB if API returned no genes
        if not genes:
            try:
                from scripts.storage.query_local_hpo import resolve_hpo_terms_local

                local = resolve_hpo_terms_local([hpo_id])
                if local:
                    if not name:
                        name = local[0].get("name", hpo_id)
                    genes = local[0].get("genes", [])
                    if genes:
                        logger.info(f"HPO local fallback: {hpo_id} → {len(genes)} genes")
            except Exception as e:
                logger.debug(f"HPO local DB not available: {e}")

        results.append({"id": hpo_id, "name": name or hpo_id, "genes": genes})

    return results


def calculate_hpo_score(gene: str, hpo_results: List[Dict]) -> int:
    """Calculate HPO overlap score for a gene."""
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
