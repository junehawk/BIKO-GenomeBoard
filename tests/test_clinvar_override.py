# tests/test_clinvar_override.py
"""Tests for ClinVar direct classification override and VUS filtering."""

import pytest
from pathlib import Path

from scripts.classification.acmg_engine import apply_clinvar_override


# ── apply_clinvar_override unit tests ─────────────────────────────────────────

def test_override_expert_panel_pathogenic():
    """Expert panel + Pathogenic → Pathogenic override."""
    result = apply_clinvar_override("VUS", "Pathogenic", "reviewed by expert panel")
    assert result == "Pathogenic"


def test_override_multi_submitter_pathogenic():
    """Multiple submitters (no conflict) + Pathogenic → Likely Pathogenic."""
    result = apply_clinvar_override(
        "VUS", "Pathogenic", "criteria provided, multiple submitters, no conflicts"
    )
    assert result == "Likely Pathogenic"


def test_override_single_submitter_no_change():
    """Single submitter — no override regardless of ClinVar significance."""
    result = apply_clinvar_override(
        "VUS", "Pathogenic", "criteria provided, single submitter"
    )
    assert result == "VUS"


def test_override_conflicting_no_change():
    """Conflicting interpretations — no override."""
    result = apply_clinvar_override(
        "VUS", "Conflicting classifications of pathogenicity",
        "criteria provided, multiple submitters, conflicting interpretations"
    )
    assert result == "VUS"


def test_override_expert_benign():
    """Expert panel + Benign → Benign override."""
    result = apply_clinvar_override("VUS", "Benign", "reviewed by expert panel")
    assert result == "Benign"


def test_override_multi_benign():
    """Multiple submitters + Benign → Likely Benign."""
    result = apply_clinvar_override(
        "Likely Pathogenic", "Benign", "criteria provided, multiple submitters, no conflicts"
    )
    assert result == "Likely Benign"


def test_override_drug_response_skipped():
    """Drug response ClinVar significance is skipped — no override."""
    result = apply_clinvar_override("VUS", "drug response", "reviewed by expert panel")
    assert result == "VUS"


def test_override_not_provided_skipped():
    """'Not provided' ClinVar is skipped — no override."""
    result = apply_clinvar_override("VUS", "not provided", "reviewed by expert panel")
    assert result == "VUS"


def test_override_empty_clinvar_no_change():
    """Empty ClinVar significance — no override."""
    result = apply_clinvar_override("Pathogenic", "", "reviewed by expert panel")
    assert result == "Pathogenic"


def test_override_empty_review_status_no_change():
    """Empty review status — no override."""
    result = apply_clinvar_override("VUS", "Pathogenic", "")
    assert result == "VUS"


def test_override_likely_pathogenic_multi():
    """Likely pathogenic with multi submitters → Likely Pathogenic."""
    result = apply_clinvar_override(
        "VUS", "Likely pathogenic", "criteria provided, multiple submitters, no conflicts"
    )
    assert result == "Likely Pathogenic"


def test_override_likely_pathogenic_expert():
    """Likely pathogenic with expert panel → Likely Pathogenic."""
    result = apply_clinvar_override("VUS", "Likely pathogenic", "reviewed by expert panel")
    assert result == "Likely Pathogenic"


def test_override_pathogenic_likely_pathogenic_compound():
    """Compound 'Pathogenic/Likely pathogenic' → Likely Pathogenic."""
    result = apply_clinvar_override(
        "VUS", "Pathogenic/Likely pathogenic", "criteria provided, multiple submitters, no conflicts"
    )
    assert result == "Likely Pathogenic"


def test_override_likely_benign_expert():
    """Likely benign with expert panel → Likely Benign."""
    result = apply_clinvar_override("VUS", "Likely benign", "reviewed by expert panel")
    assert result == "Likely Benign"


def test_override_practice_guideline():
    """Practice guideline counts as expert-level confidence."""
    result = apply_clinvar_override("VUS", "Pathogenic", "practice guideline")
    assert result == "Pathogenic"


# ── hide_vus integration tests ────────────────────────────────────────────────

DEMO_VCF = str(Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants.vcf")


def test_hide_vus_splits_variants(tmp_path):
    """With hide_vus=True, detailed_variants excludes VUS/Benign, omitted_variants contains them."""
    from scripts.orchestrate import run_pipeline

    output_path = tmp_path / "report.html"
    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
        skip_api=True,
        hide_vus=True,
    )

    assert result is not None
    vus_classes = {"VUS", "Benign", "Likely Benign"}

    # detailed_variants must not contain VUS/Benign/Likely Benign
    for v in result["detailed_variants"]:
        assert v["classification"] not in vus_classes, (
            f"detailed_variants contains {v['classification']} for {v['variant']}"
        )

    # omitted_variants must only contain VUS/Benign/Likely Benign
    for v in result["omitted_variants"]:
        assert v["classification"] in vus_classes, (
            f"omitted_variants contains unexpected {v['classification']} for {v['variant']}"
        )

    # Together they account for all variants
    total = len(result["detailed_variants"]) + len(result["omitted_variants"])
    assert total == len(result["variants"])

    # hide_vus flag is reflected in report data
    assert result["hide_vus"] is True


def test_hide_vus_default_false(tmp_path):
    """With hide_vus=False (default), detailed_variants equals all variants and omitted is empty."""
    from scripts.orchestrate import run_pipeline

    output_path = tmp_path / "report.html"
    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
        skip_api=True,
        hide_vus=False,
    )

    assert result is not None
    assert result["hide_vus"] is False
    assert result["omitted_variants"] == []
    assert len(result["detailed_variants"]) == len(result["variants"])
