# scripts/clinical/query_clinvar.py
import os
from typing import Optional, Dict, List
from scripts.common.models import Variant
from scripts.common.api_utils import fetch_with_retry

CLINVAR_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
CLINVAR_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

def _search_clinvar_variant(variant: Variant) -> Optional[dict]:
    """Search ClinVar for a variant and return summary."""
    chrom_num = variant.chrom.replace("chr", "")
    query = f"{chrom_num}[CHR] AND {variant.pos}[CHRPOS] AND {variant.ref}>{variant.alt}[VARNAME]"
    api_key = os.environ.get("NCBI_API_KEY", "")
    params = {"db": "clinvar", "term": query, "retmode": "json"}
    if api_key:
        params["api_key"] = api_key

    search_result = fetch_with_retry(CLINVAR_ESEARCH, params=params)
    if not search_result or not search_result.get("esearchresult", {}).get("idlist"):
        return None

    uid = search_result["esearchresult"]["idlist"][0]
    summary_params = {"db": "clinvar", "id": uid, "retmode": "json"}
    if api_key:
        summary_params["api_key"] = api_key

    summary = fetch_with_retry(CLINVAR_ESUMMARY, params=summary_params)
    if summary and "result" in summary and uid in summary["result"]:
        return summary["result"][uid]
    return None

def _derive_acmg_codes(clinvar_data: dict) -> List[str]:
    """Derive ACMG evidence codes from ClinVar data."""
    codes = []
    sig = clinvar_data.get("clinical_significance", {}).get("description", "").lower()
    review = clinvar_data.get("review_status", "")

    if "pathogenic" in sig and "conflict" not in sig:
        if "multiple submitters" in review:
            codes.append("PS1")
        else:
            codes.append("PP3")

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
        }

    sig = clinvar_data.get("clinical_significance", {}).get("description", "Unknown")
    acmg_codes = _derive_acmg_codes(clinvar_data)

    return {
        "agent": "clinical_geneticist",
        "variant": variant.variant_id,
        "gene": variant.gene or clinvar_data.get("gene", {}).get("symbol"),
        "clinvar_significance": sig,
        "clinvar_id": clinvar_data.get("variation_id"),
        "review_status": clinvar_data.get("review_status"),
        "acmg_codes": acmg_codes,
    }
