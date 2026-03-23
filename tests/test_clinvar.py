# tests/test_clinvar.py
import pytest
from unittest.mock import call
from scripts.clinical.query_clinvar import query_clinvar, _search_clinvar_variant
from scripts.common.models import Variant

SAMPLE_CLINVAR_RESPONSE = {
    "result": {
        "12375": {
            "clinical_significance": {"description": "Pathogenic"},
            "gene": {"symbol": "TP53"},
            "trait_set": [{"trait_name": "Li-Fraumeni syndrome"}],
            "review_status": "criteria provided, multiple submitters, no conflicts",
            "variation_id": "12375"
        }
    }
}

def test_query_clinvar_pathogenic(mocker):
    mocker.patch(
        "scripts.clinical.query_clinvar._search_clinvar_variant",
        return_value=SAMPLE_CLINVAR_RESPONSE["result"]["12375"]
    )
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")
    result = query_clinvar(variant)
    assert result is not None
    assert result["clinvar_significance"] == "Pathogenic"
    assert "PVS1" in result["acmg_codes"] or "PS1" in result["acmg_codes"]

def test_query_clinvar_not_found(mocker):
    mocker.patch(
        "scripts.clinical.query_clinvar._search_clinvar_variant",
        return_value=None
    )
    variant = Variant(chrom="chr1", pos=12345, ref="A", alt="T")
    result = query_clinvar(variant)
    assert result is not None
    assert result["clinvar_significance"] == "Not Found"
    assert result["acmg_codes"] == []

# I-5c: api_available tracking
def test_query_clinvar_api_available_true(mocker):
    mocker.patch(
        "scripts.clinical.query_clinvar._search_clinvar_variant",
        return_value=SAMPLE_CLINVAR_RESPONSE["result"]["12375"]
    )
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")
    result = query_clinvar(variant)
    assert result["api_available"] is True

def test_query_clinvar_api_available_false(mocker):
    mocker.patch(
        "scripts.clinical.query_clinvar._search_clinvar_variant",
        return_value=None
    )
    variant = Variant(chrom="chr1", pos=12345, ref="A", alt="T")
    result = query_clinvar(variant)
    assert result["api_available"] is False

def test_search_clinvar_uses_rsid_first(mocker):
    """_search_clinvar_variant should use rsID as first search strategy."""
    esearch_response = {"esearchresult": {"idlist": ["12375"]}}
    esummary_response = {
        "result": {
            "12375": SAMPLE_CLINVAR_RESPONSE["result"]["12375"]
        }
    }
    mock_fetch = mocker.patch(
        "scripts.clinical.query_clinvar.fetch_with_retry",
        side_effect=[esearch_response, esummary_response]
    )
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53", rsid="rs28934578")
    result = _search_clinvar_variant(variant)
    # First call should search by rsID term
    first_call_params = mock_fetch.call_args_list[0][1]["params"]
    assert first_call_params["term"] == "rs28934578"
    assert result is not None
    assert result["clinical_significance"]["description"] == "Pathogenic"

def test_search_clinvar_falls_back_to_gene_pos(mocker):
    """When no rsID, fall back to gene+position search."""
    esearch_response = {"esearchresult": {"idlist": ["12375"]}}
    esummary_response = {
        "result": {
            "12375": SAMPLE_CLINVAR_RESPONSE["result"]["12375"]
        }
    }
    mock_fetch = mocker.patch(
        "scripts.clinical.query_clinvar.fetch_with_retry",
        side_effect=[esearch_response, esummary_response]
    )
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53", rsid=None)
    result = _search_clinvar_variant(variant)
    first_call_params = mock_fetch.call_args_list[0][1]["params"]
    assert "TP53[GENE]" in first_call_params["term"]
    assert "7577120[CHRPOS]" in first_call_params["term"]
    assert result is not None
