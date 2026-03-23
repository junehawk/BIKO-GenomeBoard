# scripts/clinical/query_clinvar.py
import logging
import os
from typing import Optional, Dict, List
from scripts.common.models import Variant
from scripts.common.api_utils import fetch_with_retry
from scripts.common.config import get

CLINVAR_ESEARCH = get("api.clinvar_esearch", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi")
CLINVAR_ESUMMARY = get("api.clinvar_esummary", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi")

logger = logging.getLogger(__name__)

def _fetch_summary(uid: str, api_key: str = "") -> Optional[dict]:
    """Fetch ClinVar summary for a given UID."""
    summary_params = {"db": "clinvar", "id": uid, "retmode": "json"}
    if api_key:
        summary_params["api_key"] = api_key
    summary = fetch_with_retry(CLINVAR_ESUMMARY, params=summary_params)
    if summary and "result" in summary and uid in summary["result"]:
        return summary["result"][uid]
    return None

def _search_clinvar_variant(variant: Variant) -> Optional[dict]:
    """Search ClinVar for a variant and return summary."""
    api_key = get("api.ncbi_api_key", "") or os.environ.get("NCBI_API_KEY", "")

    # Strategy 1: Search by rsID (most reliable)
    if variant.rsid:
        params = {"db": "clinvar", "term": variant.rsid, "retmode": "json"}
        if api_key:
            params["api_key"] = api_key
        search_result = fetch_with_retry(CLINVAR_ESEARCH, params=params)
        if search_result and search_result.get("esearchresult", {}).get("idlist"):
            uid = search_result["esearchresult"]["idlist"][0]
            logger.info(f"ClinVar hit via rsID {variant.rsid}: uid={uid}")
            return _fetch_summary(uid, api_key)

    # Strategy 2: Search by gene + position
    if variant.gene:
        params = {"db": "clinvar", "term": f"{variant.gene}[GENE] AND {variant.pos}[CHRPOS]", "retmode": "json"}
        if api_key:
            params["api_key"] = api_key
        search_result = fetch_with_retry(CLINVAR_ESEARCH, params=params)
        if search_result and search_result.get("esearchresult", {}).get("idlist"):
            uid = search_result["esearchresult"]["idlist"][0]
            logger.info(f"ClinVar hit via gene+pos {variant.gene}/{variant.pos}: uid={uid}")
            return _fetch_summary(uid, api_key)

    return None

def _extract_significance(clinvar_data: dict) -> tuple:
    """Extract clinical significance and review status from ClinVar data.
    Handles both old format (clinical_significance) and new format (germline_classification).
    """
    # New format (2024+): germline_classification
    gc = clinvar_data.get("germline_classification", {})
    if gc and gc.get("description"):
        return gc["description"], gc.get("review_status", "")

    # Old format: clinical_significance
    cs = clinvar_data.get("clinical_significance", {})
    if cs and cs.get("description"):
        return cs["description"], clinvar_data.get("review_status", "")

    return "Unknown", ""

def _derive_acmg_codes(clinvar_data: dict) -> List[str]:
    """Derive ACMG evidence codes from ClinVar data."""
    codes = []
    sig, review = _extract_significance(clinvar_data)
    sig_lower = sig.lower()

    if "pathogenic" in sig_lower and "conflict" not in sig_lower:
        if "expert panel" in review.lower() or "practice guideline" in review.lower():
            codes.append("PS1")
            codes.append("PP5")
        elif "multiple submitters" in review.lower():
            codes.append("PS1")
        else:
            codes.append("PP5")

    return codes

def query_clinvar(variant: Variant) -> Dict:
    """Query ClinVar for variant and return structured result."""
    clinvar_data = _search_clinvar_variant(variant)

    if clinvar_data is None:
        return {
            "agent": "clinical_geneticist",
            "variant": variant.variant_id,
            "gene": variant.gene,
            "clinvar_significance": "Not Found",
            "clinvar_id": None,
            "acmg_codes": [],
            "review_status": None,
            "api_available": False,
        }

    sig, review = _extract_significance(clinvar_data)
    acmg_codes = _derive_acmg_codes(clinvar_data)
    accession = clinvar_data.get("accession", clinvar_data.get("variation_id"))

    return {
        "agent": "clinical_geneticist",
        "variant": variant.variant_id,
        "gene": variant.gene or clinvar_data.get("gene", {}).get("symbol"),
        "clinvar_significance": sig,
        "clinvar_id": accession,
        "review_status": review,
        "acmg_codes": acmg_codes,
        "api_available": True,
    }


if __name__ == "__main__":
    import sys, json
    if len(sys.argv) < 2:
        print(json.dumps({"error": "Usage: python -m scripts.clinical.query_clinvar 'chr17:7577120 G>A'"}))
        sys.exit(1)
    v = Variant.from_string(sys.argv[1])
    print(json.dumps(query_clinvar(v), indent=2))
