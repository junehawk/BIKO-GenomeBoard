"""TMB calculation tests."""
import pytest
from dataclasses import dataclass
from typing import Optional


@dataclass
class MockVariant:
    consequence: Optional[str] = None
    gene: Optional[str] = None


def test_tmb_high():
    from scripts.somatic.tmb import calculate_tmb
    variants = [MockVariant(consequence="missense_variant", gene=f"G{i}") for i in range(330)]
    result = calculate_tmb(variants, panel_size_mb=33.0)
    assert result.score == 10.0
    assert result.level == "High"
    assert result.variant_count == 330


def test_tmb_intermediate():
    from scripts.somatic.tmb import calculate_tmb
    variants = [MockVariant(consequence="missense_variant") for _ in range(231)]
    result = calculate_tmb(variants, panel_size_mb=33.0)
    assert result.level == "Intermediate"
    assert 6.0 <= result.score < 10.0


def test_tmb_low():
    from scripts.somatic.tmb import calculate_tmb
    variants = [MockVariant(consequence="missense_variant") for _ in range(100)]
    result = calculate_tmb(variants, panel_size_mb=33.0)
    assert result.level == "Low"
    assert result.score < 6.0


def test_tmb_empty():
    from scripts.somatic.tmb import calculate_tmb
    result = calculate_tmb([], panel_size_mb=33.0)
    assert result.score == 0.0
    assert result.level == "Low"
    assert result.variant_count == 0


def test_tmb_filters_synonymous():
    from scripts.somatic.tmb import calculate_tmb
    variants = [
        MockVariant(consequence="missense_variant"),
        MockVariant(consequence="missense_variant"),
        MockVariant(consequence="synonymous_variant"),
        MockVariant(consequence="3_prime_UTR_variant"),
        MockVariant(consequence="intron_variant"),
    ]
    result = calculate_tmb(variants, panel_size_mb=1.0)
    assert result.variant_count == 2
    assert result.total_variants == 5


def test_tmb_counts_all_nonsynonymous():
    from scripts.somatic.tmb import calculate_tmb
    variants = [
        MockVariant(consequence="missense_variant"),
        MockVariant(consequence="stop_gained"),
        MockVariant(consequence="frameshift_variant"),
        MockVariant(consequence="splice_donor_variant"),
        MockVariant(consequence="splice_acceptor_variant"),
        MockVariant(consequence="inframe_deletion"),
        MockVariant(consequence="inframe_insertion"),
    ]
    result = calculate_tmb(variants, panel_size_mb=1.0)
    assert result.variant_count == 7
    assert result.score == 7.0


def test_tmb_custom_panel_size():
    from scripts.somatic.tmb import calculate_tmb
    variants = [MockVariant(consequence="missense_variant") for _ in range(11)]
    result = calculate_tmb(variants, panel_size_mb=1.1)
    assert result.score == pytest.approx(10.0, abs=0.1)
    assert result.panel_size_mb == 1.1


def test_tmb_custom_thresholds():
    from scripts.somatic.tmb import calculate_tmb
    variants = [MockVariant(consequence="missense_variant") for _ in range(264)]
    result = calculate_tmb(variants, panel_size_mb=33.0, high_threshold=8.0, intermediate_threshold=4.0)
    assert result.level == "High"


def test_tmb_result_fields():
    from scripts.somatic.tmb import calculate_tmb
    variants = [MockVariant(consequence="missense_variant")]
    result = calculate_tmb(variants, panel_size_mb=33.0)
    assert hasattr(result, "score")
    assert hasattr(result, "level")
    assert hasattr(result, "variant_count")
    assert hasattr(result, "total_variants")
    assert hasattr(result, "panel_size_mb")
    assert hasattr(result, "counted_consequences")
    assert isinstance(result.counted_consequences, list)


def test_tmb_bed_file(tmp_path):
    from scripts.somatic.tmb import calculate_panel_size_from_bed
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t100\t1100\nchr1\t2000\t3000\nchr2\t500\t1500\n")
    size = calculate_panel_size_from_bed(str(bed))
    assert size == pytest.approx(0.003, abs=0.0001)
