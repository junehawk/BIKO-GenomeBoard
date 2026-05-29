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


def test_cyp2c19_builtin_does_not_assert_star_allele():
    """v2.7 PGX-02: the gene-only builtin path must NOT assert a star-allele
    metabolizer phenotype — it has no locus to verify the genotype."""
    variant = Variant(chrom="chr10", pos=96541616, ref="G", alt="A", gene="CYP2C19")
    result = check_korean_pgx(variant)
    assert result is not None
    assert result.star_allele == ""
    assert "carrier" not in result.phenotype.lower()
    assert "metabolizer" not in result.phenotype.lower()
    assert "PharmCAT" in result.phenotype
    assert "CYP2C19" in result.phenotype
    assert result.cpic_level == "A"  # gene-level fact preserved


def test_arbitrary_variant_in_pgx_gene_not_called_a_carrier():
    """An unrelated variant in a PGx gene must not be reported as a *2 carrier."""
    variant = Variant(chrom="chr10", pos=1, ref="T", alt="C", gene="CYP2C19")
    result = check_korean_pgx(variant)
    assert result is not None
    assert "*2" not in result.phenotype
    assert result.star_allele == ""


def test_hla_b_builtin_no_carrier_assertion():
    variant = Variant(chrom="chr6", pos=31353875, ref="C", alt="A", gene="HLA-B")
    result = check_korean_pgx(variant)
    if result is not None:
        assert result.star_allele == ""
        assert "carrier" not in result.phenotype.lower()


def test_korean_flag_unified_across_paths():
    """The builtin PgxResult.korean_flag and the shared helper must agree
    (v2.7 ENRI-06 — the PharmCAT path no longer uses a divergent formula)."""
    from scripts.common.models import PgxResult, pgx_korean_flag

    for korean, western in [(0.15, 0.02), (0.03, 0.03), (0.0, 0.05), (0.10, 0.0)]:
        r = PgxResult(
            gene="X",
            star_allele="",
            phenotype="",
            cpic_level="",
            korean_prevalence=korean,
            western_prevalence=western,
        )
        assert r.korean_flag == pgx_korean_flag(korean, western)
