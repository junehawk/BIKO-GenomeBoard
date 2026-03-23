from scripts.korean_pop.query_gnomad import query_gnomad
from scripts.common.models import Variant


# Mock response using ac/an schema (gnomAD v4 format)
def _make_genome_response(ac, an, eas_ac, eas_an):
    return {
        "data": {
            "variant": {
                "variant_id": "17-7577120-G-A",
                "genome": {"ac": ac, "an": an, "populations": [{"id": "eas", "ac": eas_ac, "an": eas_an}]},
                "exome": None,
            }
        }
    }


def test_query_gnomad_returns_frequencies(mocker):
    # ac=200, an=1000000 → af=0.0002; eas ac=3, an=10000 → af=0.0003
    mock_response = _make_genome_response(ac=200, an=1000000, eas_ac=3, eas_an=10000)
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value=mock_response)
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A")
    result = query_gnomad(variant)
    assert abs(result["gnomad_all"] - 0.0002) < 1e-9
    assert abs(result["gnomad_eas"] - 0.0003) < 1e-9


def test_query_gnomad_not_found(mocker):
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value={"data": {"variant": None}})
    variant = Variant(chrom="chr1", pos=99999999, ref="A", alt="T")
    result = query_gnomad(variant)
    assert result["gnomad_all"] is None
    assert result["gnomad_eas"] is None


def test_query_gnomad_fallback_to_exome(mocker):
    """When genome has no data, fall back to exome."""
    mock_response = {
        "data": {
            "variant": {
                "variant_id": "17-7577120-G-A",
                "genome": {"ac": 0, "an": 0, "populations": []},
                "exome": {"ac": 100, "an": 500000, "populations": [{"id": "eas", "ac": 5, "an": 10000}]},
            }
        }
    }
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value=mock_response)
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A")
    result = query_gnomad(variant)
    assert abs(result["gnomad_all"] - 0.0002) < 1e-9


def test_query_gnomad_multi_dataset_fallback(mocker):
    """First dataset returns no variant; second dataset returns data."""
    hit_response = _make_genome_response(ac=200, an=1000000, eas_ac=3, eas_an=10000)
    mocker.patch(
        "scripts.korean_pop.query_gnomad._graphql_query", side_effect=[{"data": {"variant": None}}, hit_response]
    )
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A")
    result = query_gnomad(variant)
    assert result["gnomad_all"] is not None
    assert result["api_available"] is True


# I-5c: api_available tracking
def test_query_gnomad_api_available_true(mocker):
    mock_response = _make_genome_response(ac=200, an=1000000, eas_ac=3, eas_an=10000)
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value=mock_response)
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A")
    result = query_gnomad(variant)
    assert result["api_available"] is True


def test_query_gnomad_api_available_false(mocker):
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value={"data": {"variant": None}})
    variant = Variant(chrom="chr1", pos=99999999, ref="A", alt="T")
    result = query_gnomad(variant)
    assert result["api_available"] is False
