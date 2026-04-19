"""Tests for case briefing builder."""

from pathlib import Path

import pytest

from scripts.clinical_board.case_briefing import build_case_briefing

# ── Paths ────────────────────────────────────────────────────────────────────

DEMO_VCF = str(Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants_grch38_annotated.vcf")

# Tests that run the full pipeline and assert specific variants appear in the
# case briefing depend on the civic.sqlite3 hotspot lookup. CI runs without
# that DB provisioned, so is_hotspot() returns False and TP53 R249M (a hotspot
# VUS the tests expect to see) gets dropped from the selector's MAY arm. Skip
# the full-pipeline briefing tests when the DB is absent — the selector is
# covered by tests/test_variant_selector_v22.py with mocked is_hotspot.
_CIVIC_DB = Path(__file__).parent.parent / "data" / "db" / "civic.sqlite3"
_requires_civic_db = pytest.mark.skipif(
    not _CIVIC_DB.exists(),
    reason="requires data/db/civic.sqlite3 for is_hotspot() (not provisioned on CI)",
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


@_requires_civic_db
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


@_requires_civic_db
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
        result["pgx_results"] = [
            {
                "gene": "CYP2C19",
                "star_allele": "*2",
                "phenotype": "Poor Metabolizer",
                "cpic_level": "A",
                "cpic_recommendation": "Use alternative drug",
                "korean_prevalence": 0.15,
                "western_prevalence": 0.03,
                "korean_flag": True,
                "clinical_impact": "Reduced metabolism",
            }
        ]
        briefing = build_case_briefing(result, "cancer")
        assert "== PHARMACOGENOMICS ==" in briefing
        assert "CYP2C19" in briefing
        assert "CPIC Level: A" in briefing


# ── Board variant selection ─────────────────────────────────────────────────


def _passenger_vus_report(n: int = 30) -> dict:
    variants = []
    for i in range(n):
        variants.append(
            {
                "variant": f"chr1:{1000 + i}:A>T",
                "gene": f"GENE{i}",
                "classification": "VUS",
                "acmg_codes": ["PM2_Supporting"],
                "clinvar_significance": "Not Found",
                "hgvsc": f"c.{100 + i}A>T",
                "hgvsp": "",
                "consequence": "missense_variant",
                "in_silico": {},
                "sift": "",
                "polyphen": "",
                "tier": "Tier IV",
                "tier_label": "VUS",
                "tier_evidence_source": "",
                "hpo_score": 0,
                "gnomad_af": None,
            }
        )
    return {
        "sample_id": "TEST_SELECT",
        "date": "2026-04-14",
        "variants": variants,
        "summary": {
            "total": n,
            "pathogenic": 0,
            "likely_pathogenic": 0,
            "vus": n,
            "likely_benign": 0,
            "benign": 0,
            "drug_response": 0,
            "risk_factor": 0,
        },
        "pgx_results": [],
    }


def test_build_briefing_selection_section_present():
    """Briefing contains the BOARD VARIANT SELECTION audit section."""
    report_data = _passenger_vus_report(5)
    briefing = build_case_briefing(report_data, "cancer")
    assert "== BOARD VARIANT SELECTION ==" in briefing
    assert "Criteria:" in briefing
    assert "Input variants: 5" in briefing


def test_build_briefing_passes_fail_selector():
    """Pure passenger VUS (no hotspot, no TSG LoF) are excluded — empty board."""
    report_data = _passenger_vus_report(30)
    briefing = build_case_briefing(report_data, "cancer")
    # None of the passenger GENEn should appear in the board variant list
    assert "--- Variant 1:" not in briefing
    assert "No variants selected for the board." in briefing
    # Empty-reason surfaces in the selection section
    assert "Empty selection:" in briefing


def test_build_briefing_mixed_pass_fail_selection():
    """Mixed fixture: only P/LP + Tier I/II survive; passengers are dropped."""
    variants = [
        # Passes: P/LP
        {
            "variant": "chr17:7675088:C>A",
            "gene": "TP53",
            "classification": "Pathogenic",
            "clinvar_significance": "Pathogenic",
            "hgvsc": "c.524G>T",
            "hgvsp": "p.Arg175Leu",
            "consequence": "missense_variant",
            "acmg_codes": ["PVS1", "PM2"],
            "tier": "Tier I",
            "tier_label": "Strong",
            "tier_evidence_source": "Tier I (AMP)",
            "in_silico": {},
        },
        # Passes: Tier II
        {
            "variant": "chr7:55191822:T>A",
            "gene": "EGFR",
            "classification": "VUS",
            "clinvar_significance": "Not Found",
            "hgvsc": "c.2369C>T",
            "hgvsp": "p.Thr790Met",
            "consequence": "missense_variant",
            "acmg_codes": [],
            "tier": "Tier II",
            "tier_label": "Potential",
            "tier_evidence_source": "Tier II (AMP)",
            "in_silico": {},
        },
        # Fails: passenger VUS
        {
            "variant": "chr1:100:A>T",
            "gene": "PASSENGER",
            "classification": "VUS",
            "clinvar_significance": "Not Found",
            "hgvsc": "c.50A>T",
            "hgvsp": "",
            "consequence": "missense_variant",
            "acmg_codes": [],
            "tier": "Tier IV",
            "tier_label": "VUS",
            "tier_evidence_source": "",
            "in_silico": {},
        },
    ]
    report_data = {
        "sample_id": "MIXED",
        "date": "2026-04-14",
        "variants": variants,
        "summary": {
            "total": 3,
            "pathogenic": 1,
            "likely_pathogenic": 0,
            "vus": 2,
            "likely_benign": 0,
            "benign": 0,
            "drug_response": 0,
            "risk_factor": 0,
        },
        "pgx_results": [],
    }
    briefing = build_case_briefing(report_data, "cancer")
    assert "TP53" in briefing
    assert "EGFR" in briefing
    # PASSENGER is in the raw list, but the selector should drop it before
    # the CLASSIFIED VARIANTS section is rendered.
    assert "--- Variant 3: PASSENGER ---" not in briefing
    assert "Selection Reason: P_LP" in briefing
    assert "Selection Reason: Tier_II" in briefing


def test_build_briefing_consumes_precomputed_selection():
    """If runner attaches _board_variants, the selector is NOT re-run."""
    pre_filtered = [
        {
            "variant": "chr17:1:G>T",
            "gene": "PRECOMPUTED_GENE",
            "classification": "Pathogenic",
            "hgvsc": "c.1G>T",
            "hgvsp": "",
            "consequence": "missense_variant",
            "acmg_codes": ["PVS1"],
            "tier": "Tier I",
            "tier_label": "Strong",
            "tier_evidence_source": "Tier I (AMP)",
            "in_silico": {},
            "selection_reason": "P_LP",
        }
    ]
    meta = {
        "mode": "cancer",
        "total_input": 100,  # runner claims 100 → briefing must echo it
        "selected": 1,
        "must_included": 1,
        "may_included": 0,
        "excluded": 99,
        "truncated": False,
        "n_dropped": 0,
        "hard_cap_applied": False,
        "empty": False,
        "empty_reason": "",
        "tmb_high_footnote": False,
        "criteria_summary": "precomputed test criteria",
        "by_selection_reason": {"P_LP": 1},
    }
    report_data = {
        "sample_id": "PRE",
        "date": "2026-04-14",
        "variants": [{"gene": "SHOULD_NOT_APPEAR", "classification": "Pathogenic"}],
        "_board_variants": pre_filtered,
        "_board_selection_metadata": meta,
        "summary": {"total": 1},
        "pgx_results": [],
    }
    briefing = build_case_briefing(report_data, "cancer")
    assert "PRECOMPUTED_GENE" in briefing
    assert "SHOULD_NOT_APPEAR" not in briefing
    assert "precomputed test criteria" in briefing
    assert "Input variants: 100" in briefing


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
