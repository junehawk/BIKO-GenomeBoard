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


# ── Integration tests ──


def test_cancer_pipeline_tmb(tmp_path):
    """Cancer pipeline에서 TMB 자동 계산."""
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path="data/sample_vcf/demo_variants_grch38_annotated.vcf",
        output_path=str(tmp_path / "report.html"),
        skip_api=True,
        mode="cancer",
    )
    assert result is not None
    assert result["tmb"] is not None
    assert result["tmb"]["score"] >= 0
    assert result["tmb"]["level"] in ("High", "Intermediate", "Low")
    assert result["tmb"]["panel_size_mb"] == 33.0


def test_rare_disease_no_tmb(tmp_path):
    """Rare disease에서는 TMB 미계산."""
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path="data/sample_vcf/rare_disease_demo.vcf",
        output_path=str(tmp_path / "report.html"),
        skip_api=True,
        mode="rare-disease",
        hpo_ids=["HP:0001250"],
    )
    assert result is not None
    assert result.get("tmb") is None


def test_tmb_in_html_report(tmp_path):
    """HTML 리포트에 TMB 뱃지가 표시됨."""
    from scripts.orchestrate import run_pipeline

    run_pipeline(
        vcf_path="data/sample_vcf/codegen-Tumor_WB.mutect.passed.vep.vcf",
        output_path=str(tmp_path / "report.html"),
        skip_api=True,
        mode="cancer",
    )
    html = (tmp_path / "report.html").read_text()
    assert "Tumor Mutational Burden" in html
    assert "mut/Mb" in html


def test_tmb_in_methodology(tmp_path):
    """Methodology 섹션에 TMB 계산 방법이 기록됨."""
    from scripts.orchestrate import run_pipeline

    run_pipeline(
        vcf_path="data/sample_vcf/demo_variants_grch38_annotated.vcf",
        output_path=str(tmp_path / "report.html"),
        skip_api=True,
        mode="cancer",
    )
    html = (tmp_path / "report.html").read_text()
    assert "nonsynonymous coding variants" in html
    assert "coding region" in html
