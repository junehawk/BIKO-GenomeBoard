from scripts.korean_pop.query_gnomad import query_gnomad
from scripts.common.models import Variant

def test_query_gnomad_returns_frequencies(mocker):
    mock_response = {
        "data": {
            "variant": {
                "genome": {
                    "af": 0.0002,
                    "populations": [
                        {"id": "eas", "af": 0.0003}
                    ]
                }
            }
        }
    }
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value=mock_response)
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A")
    result = query_gnomad(variant)
    assert result["gnomad_all"] == 0.0002
    assert result["gnomad_eas"] == 0.0003

def test_query_gnomad_not_found(mocker):
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value={"data": {"variant": None}})
    variant = Variant(chrom="chr1", pos=99999999, ref="A", alt="T")
    result = query_gnomad(variant)
    assert result["gnomad_all"] is None
    assert result["gnomad_eas"] is None

# I-5c: api_available tracking
def test_query_gnomad_api_available_true(mocker):
    mock_response = {
        "data": {
            "variant": {
                "genome": {"af": 0.0002, "populations": [{"id": "eas", "af": 0.0003}]}
            }
        }
    }
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value=mock_response)
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A")
    result = query_gnomad(variant)
    assert result["api_available"] is True

def test_query_gnomad_api_available_false(mocker):
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value={"data": {"variant": None}})
    variant = Variant(chrom="chr1", pos=99999999, ref="A", alt="T")
    result = query_gnomad(variant)
    assert result["api_available"] is False
