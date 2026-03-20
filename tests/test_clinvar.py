# tests/test_clinvar.py
import pytest
from scripts.clinical.query_clinvar import query_clinvar
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
