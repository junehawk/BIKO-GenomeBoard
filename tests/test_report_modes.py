# tests/test_report_modes.py
"""Tests for report mode abstraction (Task 1.3)."""

from pathlib import Path

from scripts.counselor.generate_pdf import generate_report_html

# ── Shared fixture data ──────────────────────────────────────────────────────

MINIMAL_REPORT = {
    "sample_id": "TEST_001",
    "date": "2026-03-20",
    "variants": [
        {
            "variant": "chr17:7577120:G>A",
            "gene": "TP53",
            "classification": "Pathogenic",
            "acmg_codes": ["PVS1", "PS1"],
            "clinvar_significance": "Pathogenic",
            "agents": {
                "clinical": {"clinvar_significance": "Pathogenic"},
                "korean_pop": {"korean_flag": "Rare variant (Korean population)"},
            },
            "conflict": False,
        }
    ],
    "pgx_results": [],
    "summary": {"total": 1, "pathogenic": 1, "vus": 0, "benign": 0},
    "db_versions": {"clinvar": "2026-03-15", "gnomad": "4.0"},
}


# ── Mode selection tests ─────────────────────────────────────────────────────


def test_cancer_mode_default():
    """generate_report_html(data) works with no mode arg — backward compatible."""
    html = generate_report_html(MINIMAL_REPORT)
    assert "BIKO GenomeBoard" in html
    assert "TP53" in html
    assert "Research Use Only" in html


def test_cancer_mode_explicit():
    """generate_report_html(data, mode='cancer') produces same output as default."""
    html_default = generate_report_html(MINIMAL_REPORT)
    html_explicit = generate_report_html(MINIMAL_REPORT, mode="cancer")
    assert html_default == html_explicit


def test_rare_disease_mode():
    """generate_report_html(data, mode='rare-disease') produces HTML with 'Rare Disease' in it."""
    html = generate_report_html(MINIMAL_REPORT, mode="rare-disease")
    assert "Rare Disease" in html
    assert "BIKO GenomeBoard" in html


def test_rare_disease_has_variant_table():
    """Rare disease report shows the variant table."""
    html = generate_report_html(MINIMAL_REPORT, mode="rare-disease")
    assert "TP53" in html
    assert "Pathogenic" in html
    assert "PVS1" in html


def test_unknown_mode_falls_back_to_cancer():
    """generate_report_html(data, mode='nonexistent') falls back to cancer template."""
    html = generate_report_html(MINIMAL_REPORT, mode="nonexistent")
    # Should render cancer template content (not raise an error)
    assert "BIKO GenomeBoard" in html
    assert "TP53" in html
    assert "Research Use Only" in html


# ── Orchestrate mode flag tests ──────────────────────────────────────────────


def test_orchestrate_mode_flag(tmp_path):
    """run_pipeline(..., mode='rare-disease') sets mode in report_data."""
    from scripts.orchestrate import run_pipeline

    demo_vcf = str(Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants.vcf")
    output_path = tmp_path / "report.html"

    result = run_pipeline(
        vcf_path=demo_vcf,
        output_path=str(output_path),
        skip_api=True,
        mode="rare-disease",
    )

    assert result is not None
    assert result["mode"] == "rare-disease"


def test_orchestrate_mode_default(tmp_path):
    """run_pipeline(...) without mode uses 'cancer' default."""
    from scripts.orchestrate import run_pipeline

    demo_vcf = str(Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants.vcf")
    output_path = tmp_path / "report.html"

    result = run_pipeline(
        vcf_path=demo_vcf,
        output_path=str(output_path),
        skip_api=True,
    )

    assert result is not None
    assert result["mode"] == "cancer"


def test_orchestrate_rare_disease_report_content(tmp_path):
    """run_pipeline with rare-disease mode writes Rare Disease HTML."""
    from scripts.orchestrate import run_pipeline

    demo_vcf = str(Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants.vcf")
    output_path = tmp_path / "report.html"

    result = run_pipeline(
        vcf_path=demo_vcf,
        output_path=str(output_path),
        skip_api=True,
        mode="rare-disease",
    )

    assert result is not None
    html = output_path.read_text()
    assert "Rare Disease" in html


def test_rare_disease_mode_no_amp_tier(tmp_path):
    """Rare disease mode는 AMP tier 필드 없이 기존 방식 유지."""
    from pathlib import Path

    from scripts.orchestrate import run_pipeline

    rare_vcf = str(Path(__file__).parent.parent / "data" / "sample_vcf" / "rare_disease_demo.vcf")
    result = run_pipeline(
        vcf_path=rare_vcf,
        output_path=str(tmp_path / "report.html"),
        skip_api=True,
        mode="rare-disease",
        hpo_ids=["HP:0001250"],
    )
    assert result is not None
    # Rare disease: classification is primary, hpo_score present
    for v in result["variants"]:
        assert "classification" in v
        assert "hpo_score" in v
