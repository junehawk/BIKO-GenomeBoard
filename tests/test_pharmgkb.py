# tests/test_pharmgkb.py
from scripts.pharma.query_pharmgkb import query_pharmgkb
from scripts.common.models import Variant


def test_query_pharmgkb_found(mocker):
    mocker.patch(
        "scripts.pharma.query_pharmgkb._fetch_pharmgkb",
        return_value={"name": "clopidogrel", "recommendation": "Use alternative antiplatelet"},
    )
    variant = Variant(chrom="chr10", pos=96541616, ref="G", alt="A", gene="CYP2C19")
    result = query_pharmgkb(variant)
    assert result is not None
    assert "clopidogrel" in result["drug_name"]
