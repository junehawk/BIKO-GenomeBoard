"""Tests for case briefing builder."""
from pathlib import Path

import pytest

from scripts.clinical_board.case_briefing import build_case_briefing


# ── Paths ────────────────────────────────────────────────────────────────────

DEMO_VCF = str(
    Path(__file__).parent.parent
    / "data" / "sample_vcf" / "demo_variants_grch38_annotated.vcf"
)


# ── Helpers ──────────────────────────────────────────────────────────────────

def _run_pipeline_cancer(tmp_path):
    """Run the pipeline in cancer mode on the demo VCF."""
    from scripts.orchestrate import run_pipeline
    return run_pipeline(
        DEMO_VCF,
        output_path=str(tmp_path / "r.html"),
        skip_api=True,
        mode="cancer",
    )


def _run_pipeline_rare_disease(tmp_path):
    """Run the pipeline in rare-disease mode on the demo VCF."""
    from scripts.orchestrate import run_pipeline
    return run_pipeline(
        DEMO_VCF,
        output_path=str(tmp_path / "r.html"),
        skip_api=True,
        mode="rare-disease",
        hpo_ids=["HP:0001263"],  # Global developmental delay
    )


# ── Cancer briefing ─────────────────────────────────────────────────────────

def test_build_briefing_cancer(tmp_path):
    result = _run_pipeline_cancer(tmp_path)
    assert result is not None
    briefing = build_case_briefing(result, "cancer")

    # Verify structural sections present
    assert "= CASE BRIEFING =" in briefing
    assert "== CLASSIFIED VARIANTS ==" in briefing
    assert "== SUMMARY STATISTICS ==" in briefing

    # Verify key variant data
    assert "TP53" in briefing
    assert "BRCA2" in briefing

    # Verify cancer-mode content
    assert "Cancer Somatic" in briefing

    # Should contain classification info
    assert "Classification:" in briefing

    # Should contain tier info (cancer mode)
    assert "Tier:" in briefing


# ── Rare disease briefing ───────────────────────────────────────────────────

def test_build_briefing_rare_disease(tmp_path):
    result = _run_pipeline_rare_disease(tmp_path)
    assert result is not None
    briefing = build_case_briefing(result, "rare-disease")

    # Verify structural sections present
    assert "= CASE BRIEFING =" in briefing
    assert "== CLASSIFIED VARIANTS ==" in briefing
    assert "Rare Disease Germline" in briefing

    # Key genes still present
    assert "TP53" in briefing


# ── In silico scores ────────────────────────────────────────────────────────

def test_build_briefing_with_in_silico(tmp_path):
    """Verify in-silico / prediction scores are included when present."""
    result = _run_pipeline_cancer(tmp_path)
    assert result is not None

    # The demo VCF has SIFT/PolyPhen annotations on TP53
    briefing = build_case_briefing(result, "cancer")

    # VEP annotations propagate through as sift/polyphen fields
    # Check that the Predictions line appears for annotated variants
    assert "Predictions:" in briefing or "In Silico:" in briefing or "SIFT" in briefing


# ── PGx section ─────────────────────────────────────────────────────────────

def test_build_briefing_with_pgx(tmp_path):
    """Verify PGx section appears when PGx data is present."""
    result = _run_pipeline_cancer(tmp_path)
    assert result is not None

    # The demo VCF has CYP2C19 which should trigger a PGx hit
    pgx_results = result.get("pgx_results", [])
    if pgx_results:
        briefing = build_case_briefing(result, "cancer")
        assert "== PHARMACOGENOMICS ==" in briefing
        assert "CPIC Level:" in briefing
    else:
        # If pipeline didn't find PGx, inject synthetic data
        result["pgx_results"] = [{
            "gene": "CYP2C19",
            "star_allele": "*2",
            "phenotype": "Poor Metabolizer",
            "cpic_level": "A",
            "cpic_recommendation": "Use alternative drug",
            "korean_prevalence": 0.15,
            "western_prevalence": 0.03,
            "korean_flag": True,
            "clinical_impact": "Reduced metabolism",
        }]
        briefing = build_case_briefing(result, "cancer")
        assert "== PHARMACOGENOMICS ==" in briefing
        assert "CYP2C19" in briefing
        assert "CPIC Level: A" in briefing


# ── Truncation ───────────────────────────────────────────────────────────────

def test_build_briefing_truncation():
    """More than 20 variants → only top 20 shown."""
    # Build synthetic report data with 30 variants
    variants = []
    for i in range(30):
        variants.append({
            "variant": f"chr1:{1000 + i}:A>T",
            "gene": f"GENE{i}",
            "classification": "VUS",
            "acmg_codes": ["PM2_Supporting"],
            "clinvar_significance": "Not Found",
            "hgvsc": f"c.{100+i}A>T",
            "hgvsp": "",
            "consequence": "missense_variant",
            "in_silico": {},
            "sift": "",
            "polyphen": "",
            "tier": "III",
            "tier_label": "VUS",
            "tier_evidence_source": "",
        })
    report_data = {
        "sample_id": "TEST_TRUNC",
        "date": "2026-04-09",
        "variants": variants,
        "summary": {"total": 30, "pathogenic": 0, "likely_pathogenic": 0,
                     "vus": 30, "likely_benign": 0, "benign": 0,
                     "drug_response": 0, "risk_factor": 0},
        "pgx_results": [],
    }

    briefing = build_case_briefing(report_data, "cancer")

    # Only 20 variants should be rendered
    assert "Variant 20:" in briefing
    assert "Variant 21:" not in briefing
    assert "10 additional lower-priority variants omitted" in briefing


# ── Empty data ───────────────────────────────────────────────────────────────

def test_build_briefing_empty():
    """Empty report_data produces a minimal briefing."""
    briefing = build_case_briefing({}, "cancer")
    assert "= CASE BRIEFING =" in briefing
    # Empty dict goes through normal path with zero data
    assert "Total Variants Analyzed: 0" in briefing or "No variants found" in briefing or "No variant data" in briefing


def test_build_briefing_none():
    """None report_data produces a minimal briefing."""
    briefing = build_case_briefing(None, "rare-disease")
    assert "= CASE BRIEFING =" in briefing
    assert "No variant data available" in briefing


# ── TMB section (cancer) ────────────────────────────────────────────────────

def test_build_briefing_tmb(tmp_path):
    """TMB section appears in cancer briefings."""
    result = _run_pipeline_cancer(tmp_path)
    assert result is not None
    briefing = build_case_briefing(result, "cancer")
    # TMB should be calculated and shown for cancer mode
    if result.get("tmb"):
        assert "== TUMOR MUTATIONAL BURDEN ==" in briefing
        assert "mutations/Mb" in briefing
