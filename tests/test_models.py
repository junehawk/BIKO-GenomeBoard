# tests/test_models.py
import pytest
from scripts.common.models import Variant, AcmgEvidence, FrequencyData, PgxResult


def test_variant_from_string():
    v = Variant.from_string("chr17:7577120 G>A")
    assert v.chrom == "chr17"
    assert v.pos == 7577120
    assert v.ref == "G"
    assert v.alt == "A"


def test_variant_from_string_with_chr_prefix():
    v = Variant.from_string("17:7577120 G>A")
    assert v.chrom == "chr17"  # auto chr prefix


def test_variant_id():
    v = Variant(chrom="chr17", pos=7577120, ref="G", alt="A")
    assert v.variant_id == "chr17:7577120:G>A"


def test_acmg_evidence_creation():
    e = AcmgEvidence(code="PVS1", source="clinical_geneticist", description="null variant")
    assert e.strength == "very_strong"  # PVS1 → very_strong auto mapping


def test_acmg_evidence_strength_mapping():
    assert AcmgEvidence(code="PS1", source="test", description="").strength == "strong"
    assert AcmgEvidence(code="PM2_Supporting", source="test", description="").strength == "supporting"
    assert AcmgEvidence(code="BA1", source="test", description="").strength == "stand_alone"
    assert AcmgEvidence(code="BS1", source="test", description="").strength == "strong"


def test_frequency_data():
    f = FrequencyData(krgdb=0.001, gnomad_eas=0.002, gnomad_all=0.003)
    assert f.korean_vs_global_ratio() == pytest.approx(0.333, rel=0.01)


def test_pgx_result():
    p = PgxResult(
        gene="CYP2C19",
        star_allele="*2",
        phenotype="Poor Metabolizer",
        cpic_level="A",
        korean_prevalence=0.15,
        western_prevalence=0.03,
    )
    assert p.korean_flag is True  # >= 5x difference
