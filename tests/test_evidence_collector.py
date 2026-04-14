# tests/test_evidence_collector.py
"""Tests for ACMG evidence collector module."""

from scripts.common.models import Variant
from scripts.classification.evidence_collector import collect_additional_evidence


def _make_variant(**kwargs) -> Variant:
    """Create a Variant with sensible defaults; override with kwargs."""
    defaults = dict(chrom="chr17", pos=7675088, ref="C", alt="A")
    defaults.update(kwargs)
    return Variant(**defaults)


# ── PVS1 ─────────────────────────────────────────────────────────────────────


def test_pvs1_stop_gained_lof_intolerant():
    """stop_gained + pLI >= 0.9 -> PVS1"""
    v = _make_variant(consequence="stop_gained", gene="TP53")
    codes = collect_additional_evidence(v, gene_info={"pli": 0.98})
    assert "PVS1" in codes
    assert "PVS1_Strong" not in codes


def test_pvs1_frameshift_lof_intolerant():
    """frameshift_variant + pLI >= 0.9 -> PVS1"""
    v = _make_variant(consequence="frameshift_variant", gene="BRCA1")
    codes = collect_additional_evidence(v, gene_info={"pli": 0.95})
    assert "PVS1" in codes


def test_pvs1_strong_moderate_lof():
    """stop_gained + 0.5 <= pLI < 0.9 -> PVS1_Strong"""
    v = _make_variant(consequence="stop_gained", gene="ATM")
    codes = collect_additional_evidence(v, gene_info={"pli": 0.72})
    assert "PVS1_Strong" in codes
    assert "PVS1" not in codes


def test_pvs1_not_applicable_tolerant():
    """stop_gained + pLI < 0.5 -> no PVS1 at any strength"""
    v = _make_variant(consequence="stop_gained", gene="FAKE")
    codes = collect_additional_evidence(v, gene_info={"pli": 0.3})
    assert "PVS1" not in codes
    assert "PVS1_Strong" not in codes


def test_pvs1_not_applicable_missense():
    """missense_variant -> no PVS1 regardless of pLI (only null variants)"""
    v = _make_variant(consequence="missense_variant", gene="TP53")
    codes = collect_additional_evidence(v, gene_info={"pli": 0.99})
    assert "PVS1" not in codes
    assert "PVS1_Strong" not in codes


def test_pvs1_boundary_0_9():
    """pLI exactly 0.9 -> PVS1 (>= threshold)"""
    v = _make_variant(consequence="stop_gained", gene="TP53")
    codes = collect_additional_evidence(v, gene_info={"pli": 0.9})
    assert "PVS1" in codes


def test_pvs1_boundary_0_5():
    """pLI exactly 0.5 -> PVS1_Strong (>= 0.5 threshold)"""
    v = _make_variant(consequence="stop_gained", gene="ATM")
    codes = collect_additional_evidence(v, gene_info={"pli": 0.5})
    assert "PVS1_Strong" in codes


# ── PM1 ──────────────────────────────────────────────────────────────────────


def test_pm1_domain_missense():
    """missense + DOMAINS present -> PM1"""
    v = _make_variant(consequence="missense_variant", gene="TP53")
    # Add domains attribute dynamically (simulating VEP annotation)
    v.domains = "Pfam:PF00870"
    codes = collect_additional_evidence(v)
    assert "PM1" in codes


def test_pm1_no_domain():
    """missense without DOMAINS -> no PM1"""
    v = _make_variant(consequence="missense_variant", gene="TP53")
    codes = collect_additional_evidence(v)
    assert "PM1" not in codes


def test_pm1_not_missense():
    """Non-missense with DOMAINS -> no PM1"""
    v = _make_variant(consequence="stop_gained", gene="TP53")
    v.domains = "Pfam:PF00870"
    codes = collect_additional_evidence(v, gene_info={"pli": 0.1})
    assert "PM1" not in codes


# ── PM4 ──────────────────────────────────────────────────────────────────────


def test_pm4_inframe_deletion():
    """inframe_deletion -> PM4"""
    v = _make_variant(consequence="inframe_deletion", gene="ERBB2")
    codes = collect_additional_evidence(v)
    assert "PM4" in codes


def test_pm4_inframe_insertion():
    """inframe_insertion -> PM4"""
    v = _make_variant(consequence="inframe_insertion", gene="EGFR")
    codes = collect_additional_evidence(v)
    assert "PM4" in codes


def test_pm4_not_frameshift():
    """frameshift_variant -> NOT PM4 (that is PVS1)"""
    v = _make_variant(consequence="frameshift_variant", gene="BRCA1")
    codes = collect_additional_evidence(v)
    assert "PM4" not in codes


# ── PM5 ──────────────────────────────────────────────────────────────────────


def test_pm5_novel_at_pathogenic_position():
    """Different AA at same pos as ClinVar pathogenic -> PM5"""
    v = _make_variant(
        consequence="missense_variant",
        gene="TP53",
        hgvsp="p.Arg248Leu",  # Novel change at position 248
    )
    # Position 248 has known pathogenic variant (e.g., R248W)
    codes = collect_additional_evidence(v, clinvar_pathogenic_positions={248, 273})
    assert "PM5" in codes


def test_pm5_no_clinvar_positions():
    """No ClinVar position data -> no PM5"""
    v = _make_variant(
        consequence="missense_variant",
        gene="TP53",
        hgvsp="p.Arg248Leu",
    )
    codes = collect_additional_evidence(v)
    assert "PM5" not in codes


def test_pm5_position_not_in_clinvar():
    """Variant position not in ClinVar pathogenic set -> no PM5"""
    v = _make_variant(
        consequence="missense_variant",
        gene="TP53",
        hgvsp="p.Arg100Leu",
    )
    codes = collect_additional_evidence(v, clinvar_pathogenic_positions={248, 273})
    assert "PM5" not in codes


# ── PP2 ──────────────────────────────────────────────────────────────────────


def test_pp2_missense_constrained():
    """missense + missense_z >= 3.09 -> PP2"""
    v = _make_variant(consequence="missense_variant", gene="TP53")
    codes = collect_additional_evidence(v, gene_info={"missense_z": 4.5})
    assert "PP2" in codes


def test_pp2_missense_not_constrained():
    """missense + missense_z < 3.09 -> no PP2"""
    v = _make_variant(consequence="missense_variant", gene="BARD1")
    codes = collect_additional_evidence(v, gene_info={"missense_z": 1.5})
    assert "PP2" not in codes


def test_pp2_boundary():
    """missense_z exactly 3.09 -> PP2 (>= threshold)"""
    v = _make_variant(consequence="missense_variant", gene="TP53")
    codes = collect_additional_evidence(v, gene_info={"missense_z": 3.09})
    assert "PP2" in codes


# ── BP1 ──────────────────────────────────────────────────────────────────────


def test_bp1_missense_lof_tolerant():
    """missense + pLI < 0.1 -> BP1"""
    v = _make_variant(consequence="missense_variant", gene="FAKE")
    codes = collect_additional_evidence(v, gene_info={"pli": 0.05})
    assert "BP1" in codes


def test_bp1_not_triggered_high_pli():
    """missense + pLI >= 0.1 -> no BP1"""
    v = _make_variant(consequence="missense_variant", gene="TP53")
    codes = collect_additional_evidence(v, gene_info={"pli": 0.95})
    assert "BP1" not in codes


# ── BP7 ──────────────────────────────────────────────────────────────────────


def test_bp7_synonymous_no_splice():
    """synonymous + no SpliceAI data -> BP7"""
    v = _make_variant(consequence="synonymous_variant", gene="CDH1")
    codes = collect_additional_evidence(v)
    assert "BP7" in codes


def test_bp7_synonymous_low_splice():
    """synonymous + spliceai_max < 0.1 -> BP7"""
    v = _make_variant(consequence="synonymous_variant", gene="CDH1")
    v.spliceai_max = 0.02
    codes = collect_additional_evidence(v)
    assert "BP7" in codes


def test_bp7_synonymous_high_splice():
    """synonymous + spliceai_max >= 0.1 -> no BP7 (splice impact)"""
    v = _make_variant(consequence="synonymous_variant", gene="CDH1")
    v.spliceai_max = 0.85
    codes = collect_additional_evidence(v)
    assert "BP7" not in codes


def test_bp7_not_synonymous():
    """missense -> no BP7 (only for synonymous)"""
    v = _make_variant(consequence="missense_variant", gene="CDH1")
    codes = collect_additional_evidence(v)
    assert "BP7" not in codes


# ── Edge cases ───────────────────────────────────────────────────────────────


def test_no_gene_info():
    """gene_info=None -> only consequence-based codes (PM4, BP7 etc.)"""
    v = _make_variant(consequence="inframe_deletion", gene="ERBB2")
    codes = collect_additional_evidence(v, gene_info=None)
    assert "PM4" in codes
    assert "PVS1" not in codes
    assert "PP2" not in codes
    assert "BP1" not in codes


def test_no_consequence():
    """Variant without consequence -> empty list"""
    v = _make_variant(gene="TP53")
    codes = collect_additional_evidence(v, gene_info={"pli": 0.99})
    assert codes == []


def test_no_consequence_none():
    """Variant with consequence=None -> empty list"""
    v = _make_variant(consequence=None, gene="TP53")
    codes = collect_additional_evidence(v, gene_info={"pli": 0.99})
    assert codes == []


def test_combined_evidence():
    """Variant matching multiple rules -> all codes returned."""
    v = _make_variant(
        consequence="missense_variant",
        gene="TP53",
        hgvsp="p.Arg248Leu",
    )
    v.domains = "Pfam:PF00870"
    gene_info = {"pli": 0.98, "missense_z": 4.5}
    codes = collect_additional_evidence(
        v,
        gene_info=gene_info,
        clinvar_pathogenic_positions={248},
    )
    assert "PM1" in codes
    assert "PM5" in codes
    assert "PP2" in codes
    # PVS1 should NOT be present (missense, not null)
    assert "PVS1" not in codes
    # BP1 should NOT be present (high pLI)
    assert "BP1" not in codes


def test_multi_consequence_string():
    """VEP multi-consequence like 'missense_variant&splice_region_variant' uses first term."""
    v = _make_variant(
        consequence="missense_variant&splice_region_variant",
        gene="TP53",
    )
    v.domains = "Pfam:PF00870"
    codes = collect_additional_evidence(v)
    assert "PM1" in codes  # missense with domain


def test_empty_consequence_string():
    """Empty string consequence -> empty list"""
    v = _make_variant(consequence="", gene="TP53")
    codes = collect_additional_evidence(v, gene_info={"pli": 0.99})
    assert codes == []


def test_gene_info_missing_keys():
    """gene_info dict with missing keys should not error."""
    v = _make_variant(consequence="stop_gained", gene="TP53")
    # gene_info without 'pli' key
    codes = collect_additional_evidence(v, gene_info={"missense_z": 4.0})
    # PVS1 needs pli, should skip gracefully
    assert "PVS1" not in codes
    assert "PVS1_Strong" not in codes
