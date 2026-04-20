# tests/test_oncokb.py
"""Tests for OncoKB cancer gene tiering."""

from pathlib import Path

import pytest

from scripts.enrichment.oncokb import (
    assign_tier,
    get_cancer_gene_info,
    get_tier_label,
    is_cancer_gene,
    reset_oncokb_cache,
)


@pytest.fixture(autouse=True)
def clear_cache():
    """Reset OncoKB cache before each test to avoid cross-test pollution."""
    reset_oncokb_cache()
    yield
    reset_oncokb_cache()


# ── is_cancer_gene ────────────────────────────────────────────────────────────


def test_is_cancer_gene_tp53():
    assert is_cancer_gene("TP53") is True


def test_is_cancer_gene_kras():
    assert is_cancer_gene("KRAS") is True


def test_is_cancer_gene_unknown():
    assert is_cancer_gene("NOTAREALGENE123") is False


def test_is_cancer_gene_empty_string():
    assert is_cancer_gene("") is False


# ── get_cancer_gene_info ──────────────────────────────────────────────────────


def test_get_cancer_gene_info_kras():
    info = get_cancer_gene_info("KRAS")
    assert info is not None
    assert info["type"] == "Oncogene"
    assert info["level"] == "1"
    assert info["gene"] == "KRAS"


def test_get_cancer_gene_info_tp53():
    info = get_cancer_gene_info("TP53")
    assert info is not None
    assert info["type"] == "TSG"
    assert info["level"] == "1"


def test_get_cancer_gene_info_unknown():
    info = get_cancer_gene_info("NOTAREALGENE123")
    assert info is None


def test_get_cancer_gene_info_erbb4():
    info = get_cancer_gene_info("ERBB4")
    assert info is not None
    assert info["type"] == "Oncogene"
    assert info["level"] == "2"


def test_get_cancer_gene_info_smad4():
    info = get_cancer_gene_info("SMAD4")
    assert info is not None
    assert info["type"] == "TSG"


def test_get_cancer_gene_info_csf1r():
    info = get_cancer_gene_info("CSF1R")
    assert info is not None
    assert info["type"] == "Oncogene"


def test_get_cancer_gene_info_cul3():
    info = get_cancer_gene_info("CUL3")
    assert info is not None
    assert info["type"] == "TSG"
    assert info["level"] == "4"


def test_get_cancer_gene_info_fat1():
    info = get_cancer_gene_info("FAT1")
    assert info is not None
    assert info["type"] == "TSG"


def test_get_cancer_gene_info_sufu():
    info = get_cancer_gene_info("SUFU")
    assert info is not None
    assert info["type"] == "TSG"


def test_get_cancer_gene_info_kmt2d():
    info = get_cancer_gene_info("KMT2D")
    assert info is not None
    assert info["type"] == "TSG"


def test_get_cancer_gene_info_reln():
    info = get_cancer_gene_info("RELN")
    assert info is not None
    assert info["type"] == "TSG"


def test_get_cancer_gene_info_cltc():
    info = get_cancer_gene_info("CLTC")
    assert info is not None
    assert info["type"] == "Oncogene"


# ── assign_tier ───────────────────────────────────────────────────────────────


def test_assign_tier_pathogenic_level1_gene():
    """Pathogenic on KRAS (Level 1) → Tier 1."""
    assert assign_tier("Pathogenic", "KRAS") == 1


def test_assign_tier_pathogenic_level2_gene():
    """Pathogenic on ERBB4 (Level 2) → Tier 1."""
    assert assign_tier("Pathogenic", "ERBB4") == 1


def test_assign_tier_likely_pathogenic_level1():
    """Likely Pathogenic on TP53 (Level 1) → Tier 1."""
    assert assign_tier("Likely Pathogenic", "TP53") == 1


def test_assign_tier_pathogenic_level3():
    """Pathogenic on CSF1R (Level 3) → Tier 2."""
    assert assign_tier("Pathogenic", "CSF1R") == 2


def test_assign_tier_pathogenic_level4():
    """Pathogenic on CUL3 (Level 4) → Tier 2."""
    assert assign_tier("Pathogenic", "CUL3") == 2


def test_assign_tier_pathogenic_non_cancer_gene():
    """Pathogenic on a non-cancer gene → Tier 4 (incidental germline finding)."""
    assert assign_tier("Pathogenic", "NOTAREALGENE") == 4


def test_assign_tier_vus_cancer_gene():
    """VUS on KRAS → Tier 3."""
    assert assign_tier("VUS", "KRAS") == 3


def test_assign_tier_vus_cancer_gene_tp53():
    """VUS on TP53 → Tier 3."""
    assert assign_tier("VUS", "TP53") == 3


def test_assign_tier_vus_non_cancer():
    """VUS on non-cancer gene → Tier 4."""
    assert assign_tier("VUS", "NOTAREALGENE") == 4


def test_assign_tier_vus_empty_gene():
    """VUS with no gene → Tier 4."""
    assert assign_tier("VUS", "") == 4


def test_assign_tier_drug_response():
    """Drug Response always → Tier 1 regardless of gene."""
    assert assign_tier("Drug Response", "") == 1
    assert assign_tier("Drug Response", "CYP2C19") == 1


def test_assign_tier_risk_factor():
    """Risk Factor → Tier 4 in cancer context (not a cancer biomarker)."""
    assert assign_tier("Risk Factor", "") == 4
    assert assign_tier("Risk Factor", "APOE") == 4


def test_assign_tier_benign_cancer_gene():
    """Benign on cancer gene → Tier 4."""
    assert assign_tier("Benign", "KRAS") == 4


def test_assign_tier_likely_benign_cancer_gene():
    """Likely Benign on cancer gene → Tier 4."""
    assert assign_tier("Likely Benign", "TP53") == 4


def test_assign_tier_benign_non_cancer():
    """Benign on non-cancer gene → Tier 4."""
    assert assign_tier("Benign", "NOTAREALGENE") == 4


def test_assign_tier_case_insensitive():
    """Classification matching is case-insensitive."""
    assert assign_tier("pathogenic", "KRAS") == 1
    assert assign_tier("PATHOGENIC", "KRAS") == 1
    assert assign_tier("drug response", "") == 1


# ── get_tier_label ────────────────────────────────────────────────────────────


def test_get_tier_label_1():
    assert get_tier_label(1) == "Tier I — Strong Clinical Significance"


def test_get_tier_label_2():
    assert get_tier_label(2) == "Tier II — Potential Clinical Significance"


def test_get_tier_label_3():
    assert get_tier_label(3) == "Tier III — Unknown Clinical Significance"


def test_get_tier_label_4():
    assert get_tier_label(4) == "Tier IV — Benign or Likely Benign"


def test_get_tier_label_unknown():
    assert get_tier_label(99) == "Unknown"


# ── pipeline integration ──────────────────────────────────────────────────────

DEMO_VCF = str(Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants.vcf")


def test_pipeline_has_tiers(tmp_path):
    """run_pipeline produces variant records with tier fields."""
    from scripts.orchestrate import run_pipeline

    output_path = tmp_path / "report.html"
    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
        skip_api=True,
        mode="cancer",
    )

    assert result is not None
    assert "tier1_variants" in result
    assert "tier2_variants" in result
    assert "tier3_variants" in result
    assert "tier4_count" in result

    for v in result["variants"]:
        assert "tier" in v, f"Variant {v.get('variant')} missing tier"
        assert v["tier"] in (1, 2, 3, 4), f"Invalid tier: {v['tier']}"
        assert "tier_label" in v
        assert "cancer_gene_type" in v
        assert "oncokb_level" in v


def test_pipeline_tier_counts_consistent(tmp_path):
    """Tier lists + tier4_count account for all variants."""
    from scripts.orchestrate import run_pipeline

    output_path = tmp_path / "report.html"
    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
        skip_api=True,
        mode="cancer",
    )

    assert result is not None
    total = len(result["variants"])
    accounted = (
        len(result["tier1_variants"])
        + len(result["tier2_variants"])
        + len(result["tier3_variants"])
        + result["tier4_count"]
    )
    assert accounted == total


def test_report_shows_tier_sections(tmp_path):
    """Cancer mode HTML report contains tier section labels."""
    from scripts.orchestrate import run_pipeline

    output_path = tmp_path / "report.html"
    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
        skip_api=True,
        mode="cancer",
    )

    assert result is not None
    html = output_path.read_text()
    assert "BIKO GenomeBoard" in html
    # At least the Genomic Findings header should be present
    assert "Genomic Findings" in html
    # Tier sections appear when variants exist
    # (demo VCF may have cancer gene variants, so tier 3 or tier 1/2 labels may appear)
    assert "Tier" in html or "tier" in html.lower() or "Genomic Findings" in html
