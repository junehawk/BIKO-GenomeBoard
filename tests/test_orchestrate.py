# tests/test_orchestrate.py
import json
from pathlib import Path

DEMO_VCF = str(Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants.vcf")


def test_orchestrate_demo_vcf(mocker, tmp_path):
    """End-to-end: orchestrate.py produces report from demo VCF"""
    mocker.patch("scripts.clinical.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value=None)

    output_path = tmp_path / "report.html"
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
    )

    assert result is not None
    assert output_path.exists()
    assert "GenomeBoard" in output_path.read_text()


def test_orchestrate_skip_api(tmp_path):
    """Offline mode skips API calls and still produces a report"""
    output_path = tmp_path / "report.html"
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
        skip_api=True,
    )

    assert result is not None
    assert output_path.exists()
    assert "GenomeBoard" in output_path.read_text()


def test_orchestrate_json_output(mocker, tmp_path):
    """JSON data written alongside HTML when json_output is specified"""
    mocker.patch("scripts.clinical.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value=None)

    output_path = tmp_path / "report.html"
    json_path = tmp_path / "report.json"
    from scripts.orchestrate import run_pipeline

    run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
        json_output=str(json_path),
    )

    assert json_path.exists()
    data = json.loads(json_path.read_text())
    assert data["summary"]["total"] == 10


def test_orchestrate_missing_vcf(tmp_path):
    """Missing VCF file returns None without crashing"""
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path=str(tmp_path / "nonexistent.vcf"),
        output_path=str(tmp_path / "report.html"),
        skip_api=True,
    )

    assert result is None


def test_orchestrate_sample_id(mocker, tmp_path):
    """Custom sample_id is reflected in report data"""
    mocker.patch("scripts.clinical.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value=None)

    output_path = tmp_path / "report.html"
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
        sample_id="PATIENT_XYZ",
    )

    assert result is not None
    assert result["sample_id"] == "PATIENT_XYZ"


def test_orchestrate_summary_counts(mocker, tmp_path):
    """Summary counts in report data are consistent with variant records"""
    mocker.patch("scripts.clinical.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value=None)

    output_path = tmp_path / "report.html"
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
    )

    assert result is not None
    summary = result["summary"]
    variants = result["variants"]

    assert summary["total"] == len(variants)
    assert summary["pathogenic"] == sum(1 for v in variants if v["classification"] == "Pathogenic")
    assert summary["drug_response"] == sum(1 for v in variants if v["classification"] == "Drug Response")


def test_orchestrate_pgx_results(mocker, tmp_path):
    """PGx hits are captured in report data for known PGx genes"""
    mocker.patch("scripts.clinical.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value=None)

    output_path = tmp_path / "report.html"
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
    )

    assert result is not None
    # demo VCF has CYP2C19 and HLA-B which are PGx genes
    pgx_genes = {r["gene"] for r in result["pgx_results"]}
    assert "CYP2C19" in pgx_genes or "HLA-B" in pgx_genes


def test_orchestrate_krgdb_missing(mocker, tmp_path):
    """Missing KRGDB file logs a warning but pipeline continues"""
    mocker.patch("scripts.clinical.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value=None)

    output_path = tmp_path / "report.html"
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
        krgdb_path=str(tmp_path / "no_such_krgdb.tsv"),
    )

    assert result is not None
    # All KRGDB frequencies should be None
    for v in result["variants"]:
        assert v["krgdb_freq"] is None


def test_orchestrate_report_structure(mocker, tmp_path):
    """Report data contains all required top-level keys"""
    mocker.patch("scripts.clinical.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.korean_pop.query_gnomad._graphql_query", return_value=None)

    output_path = tmp_path / "report.html"
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
    )

    assert result is not None
    for key in ("sample_id", "date", "variants", "pgx_results", "summary", "db_versions"):
        assert key in result, f"Missing key: {key}"
    for v in result["variants"]:
        for vkey in ("variant", "gene", "classification", "acmg_codes", "conflict"):
            assert vkey in v, f"Variant record missing key: {vkey}"


def test_cancer_pipeline_produces_amp_tier_fields(tmp_path):
    """Cancer pipeline output에 AMP tier 필드가 포함되어야 함."""
    from scripts.orchestrate import run_pipeline

    annotated_vcf = str(Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants_grch38_annotated.vcf")
    result = run_pipeline(
        vcf_path=annotated_vcf,
        output_path=str(tmp_path / "report.html"),
        skip_api=True,
        mode="cancer",
    )
    assert result is not None
    for v in result["variants"]:
        assert "tier" in v
        assert "tier_label" in v
        assert "tier_evidence_source" in v
        assert v["tier"] in (1, 2, 3, 4)
        assert "Tier" in v["tier_label"]
