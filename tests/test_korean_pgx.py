# tests/test_korean_pgx.py
from scripts.pharma.korean_pgx import check_korean_pgx
from scripts.common.models import Variant

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
