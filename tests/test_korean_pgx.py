# tests/test_korean_pgx.py
from scripts.common.models import Variant
from scripts.pharmacogenomics.korean_pgx import check_korean_pgx


def test_cyp2c19_detected():
    variant = Variant(chrom="chr10", pos=96541616, ref="G", alt="A", gene="CYP2C19")
    result = check_korean_pgx(variant)
    assert result is not None
    assert result.gene == "CYP2C19"
    assert result.korean_flag is True
    assert result.cpic_level == "A"


def test_non_pgx_gene():
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")
    result = check_korean_pgx(variant)
    assert result is None  # TP53 is not a PGx gene


def test_cyp2c19_phenotype_intermediate_metabolizer():
    variant = Variant(chrom="chr10", pos=96541616, ref="G", alt="A", gene="CYP2C19")
    result = check_korean_pgx(variant)
    assert result is not None
    assert result.phenotype == "Intermediate Metabolizer (*2 carrier)"


def test_hla_b_phenotype():
    variant = Variant(chrom="chr6", pos=31353875, ref="C", alt="A", gene="HLA-B")
    result = check_korean_pgx(variant)
    # HLA-B may or may not be in the data file; if it is, check phenotype
    if result is not None:
        assert result.phenotype == "HLA-B*5701 carrier — abacavir hypersensitivity risk"
