# tests/test_orchestrate.py
import json
from pathlib import Path

DEMO_VCF = str(Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants.vcf")


def test_orchestrate_demo_vcf(mocker, tmp_path):
    """End-to-end: orchestrate.py produces report from demo VCF"""
    mocker.patch("scripts.enrichment.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.population.query_gnomad._graphql_query", return_value=None)

    output_path = tmp_path / "report.html"
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
    )

    assert result is not None
    assert output_path.exists()
    assert "BIKO GenomeBoard" in output_path.read_text()


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
    assert "BIKO GenomeBoard" in output_path.read_text()


def test_orchestrate_json_output(mocker, tmp_path):
    """JSON data written alongside HTML when json_output is specified"""
    mocker.patch("scripts.enrichment.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.population.query_gnomad._graphql_query", return_value=None)

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
    mocker.patch("scripts.enrichment.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.population.query_gnomad._graphql_query", return_value=None)

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
    mocker.patch("scripts.enrichment.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.population.query_gnomad._graphql_query", return_value=None)

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
    mocker.patch("scripts.enrichment.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.population.query_gnomad._graphql_query", return_value=None)

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


def test_orchestrate_kova_missing(mocker, tmp_path, monkeypatch):
    """Missing KOVA DB logs a single warning but pipeline continues"""
    mocker.patch("scripts.enrichment.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.population.query_gnomad._graphql_query", return_value=None)

    # Point query_kova at a non-existent DB path; availability cache should
    # fire exactly one WARNING and every subsequent query short-circuits to None.
    from scripts.population import query_kova as qk

    monkeypatch.setattr(qk, "_DEFAULT_DB_PATH", str(tmp_path / "no_such_kova.sqlite3"))
    mocker.patch.object(qk, "get", return_value=str(tmp_path / "no_such_kova.sqlite3"))
    qk.reset_availability_cache()

    output_path = tmp_path / "report.html"
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
    )

    assert result is not None
    # All KOVA frequencies should be None
    for v in result["variants"]:
        assert v.get("kova_freq") is None


def test_orchestrate_report_structure(mocker, tmp_path):
    """Report data contains all required top-level keys"""
    mocker.patch("scripts.enrichment.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.population.query_gnomad._graphql_query", return_value=None)

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


def test_pipeline_works_with_empty_gene_knowledge(tmp_path):
    """gene_knowledge.json이 비어있어도 파이프라인 정상 동작."""
    from scripts.orchestrate import run_pipeline

    annotated_vcf = str(Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants_grch38_annotated.vcf")
    result = run_pipeline(
        vcf_path=annotated_vcf,
        output_path=str(tmp_path / "report.html"),
        skip_api=True,
        mode="cancer",
    )
    assert result is not None
    assert len(result["variants"]) > 0
    # CIViC enrichment still works at runtime for cancer
    for v in result["variants"]:
        assert "classification" in v
        assert "tier" in v


def test_run_pgx_forwards_proband_sample_id(monkeypatch):
    """_run_pgx must forward the resolved proband sample_id to get_pgx_results, so
    PharmCAT's multi-sample report selection can pick the proband instead of silently
    defaulting to the alphabetically-first sample report (PGX-01)."""
    import scripts.pharmacogenomics.korean_pgx as kpgx
    from scripts.orchestration.canonical import _run_pgx

    captured = {}

    def fake_get_pgx_results(variants, germline_vcf=None, sample_id=None, config=None):
        captured["sample_id"] = sample_id
        return {
            "pgx_source": "pharmcat",
            "pgx_hits": [],
            "germline_provided": True,
            "pharmcat_version": "3.2.0",
            "warnings": [],
        }

    monkeypatch.setattr(kpgx, "get_pgx_results", fake_get_pgx_results)

    _run_pgx([], germline_vcf="/tmp/germline.vcf", db_results={}, sample_id="PROBAND_01")
    assert captured["sample_id"] == "PROBAND_01"


def test_linkify_pmids_escapes_surrounding_markup():
    """CIViC / gene_knowledge prose passed to _linkify_pmids is rendered via {{ ...|safe }};
    markup in the source text must be escaped so only the linkifier's own anchors are raw HTML
    (XSEC-02)."""
    from scripts.orchestration.canonical import _linkify_pmids

    out = _linkify_pmids("see <script>alert(1)</script> in PMID:12345 done")
    assert "<script>" not in out
    assert "&lt;script&gt;" in out
    assert 'href="https://pubmed.ncbi.nlm.nih.gov/12345/"' in out


def test_pgx_warnings_surfaced_in_report_data(mocker, tmp_path):
    """PGx fallback warnings (why PGx degraded to the builtin table) must reach report_data
    so the reviewer can see them, instead of being computed and silently dropped (T4-04/ORCH-03)."""
    mocker.patch("scripts.enrichment.query_clinvar._search_clinvar_variant", return_value=None)
    mocker.patch("scripts.population.query_gnomad._graphql_query", return_value=None)
    mocker.patch(
        "scripts.pharmacogenomics.korean_pgx.get_pgx_results",
        return_value={
            "pgx_source": "builtin",
            "pgx_hits": [],
            "germline_provided": False,
            "pharmcat_version": "",
            "warnings": ["PharmCAT not available (Java 17+ or JAR missing); falling back to builtin PGx"],
        },
    )
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(vcf_path=DEMO_VCF, output_path=str(tmp_path / "r.html"))
    assert any("PharmCAT not available" in w for w in result.get("pgx_warnings", []))
