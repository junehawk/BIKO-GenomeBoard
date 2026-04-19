# tests/conftest.py
import pytest

from scripts.common.models import FrequencyData, Variant


@pytest.fixture
def sample_tp53_variant():
    return Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")


@pytest.fixture
def sample_brca2_variant():
    return Variant(chrom="chr13", pos=32337326, ref="C", alt="T", gene="BRCA2")


@pytest.fixture
def sample_cyp2c19_variant():
    return Variant(chrom="chr10", pos=96541616, ref="G", alt="A", gene="CYP2C19")


@pytest.fixture
def mock_clinvar_pathogenic():
    return {
        "clinical_significance": {"description": "Pathogenic"},
        "gene": {"symbol": "TP53"},
        "trait_set": [{"trait_name": "Li-Fraumeni syndrome"}],
        "review_status": "criteria provided, multiple submitters, no conflicts",
        "variation_id": "12375",
        "molecular_consequence": "nonsense",
    }


@pytest.fixture
def mock_clinvar_not_found():
    return None


@pytest.fixture
def mock_gnomad_response():
    return {"data": {"variant": {"genome": {"af": 0.0002, "populations": [{"id": "eas", "af": 0.0003}]}}}}


@pytest.fixture
def mock_gnomad_not_found():
    return {"data": {"variant": None}}


@pytest.fixture
def sample_frequency_data():
    return FrequencyData(krgdb=0.0001, gnomad_eas=0.0003, gnomad_all=0.0002)
