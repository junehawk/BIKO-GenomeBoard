"""v2.2 Phase B1 regression — protein-impacting consequence gate.

Scope:
- B1: Tier III (and Tier I/II) consequence gate rejects non-coding VUS.
- B1: P_LP branch bypasses the gate unconditionally.
- B1: SpliceAI >= 0.2 rescues synonymous / splice_region consequences.

All tests mock ``is_hotspot`` / ``is_cancer_gene`` to isolate the
selector logic from CIViC DB availability. B2 MMR/Lynch carve-out
tests live in the sibling commit that adds the carve-out.
"""
from __future__ import annotations

from unittest.mock import patch

import pytest


def _mk_snv(
    gene: str,
    classification: str = "VUS",
    tier: str = "Tier IV",
    hgvsp: str = "",
    consequence: str = "missense_variant",
    cancer_gene_type: str = "",
    oncokb_level: str = "",
    hpo_score: int = 0,
    gnomad_af=None,
    in_silico: dict | None = None,
    variant_type: str = "SNV",
) -> dict:
    return {
        "gene": gene,
        "variant": f"{gene}_var",
        "classification": classification,
        "tier": tier,
        "hgvsp": hgvsp,
        "consequence": consequence,
        "cancer_gene_type": cancer_gene_type,
        "oncokb_level": oncokb_level,
        "hpo_score": hpo_score,
        "gnomad_af": gnomad_af,
        "in_silico": in_silico or {},
        "variant_type": variant_type,
    }


def _patch_clinical(*, hotspot_positions=None, cancer_genes=None):
    hotspot_positions = hotspot_positions or set()
    cancer_genes = cancer_genes or set()

    def _hs(gene: str, pos: int) -> bool:
        return (gene, pos) in hotspot_positions

    def _cg(gene: str) -> bool:
        return gene in cancer_genes

    return patch.multiple(
        "scripts.clinical_board.variant_selector",
        is_hotspot=_hs,
        is_cancer_gene=_cg,
    )


# ---------------------------------------------------------------------------
# B1 — consequence gate
# ---------------------------------------------------------------------------


def test_notch2_intronic_tier3_rejected():
    """Tier III VUS with intron_variant consequence must NOT be admitted,
    even if the gene is in the OncoKB cancer gene list. Closes the F1 leak
    for NOTCH2 intronic VUS diluting the clinician attention.
    """
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "NOTCH2",
            classification="VUS",
            tier="Tier III",
            hgvsp="",
            consequence="intron_variant",
            oncokb_level="1",
        ),
    ]
    with _patch_clinical(cancer_genes={"NOTCH2"}):
        selected, meta = select_board_variants(variants, mode="cancer")

    assert selected == []
    assert meta["empty"] is True


def test_kras_g12d_missense_tier3_admitted():
    """Tier III missense VUS at a known hotspot position IS admitted.
    Sanity check that the consequence gate does not over-reject the
    clinically valuable cases.
    """
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "KRAS",
            classification="VUS",
            tier="Tier III",
            hgvsp="p.Gly12Asp",
            consequence="missense_variant",
        ),
    ]
    with _patch_clinical(hotspot_positions={("KRAS", 12)}):
        selected, meta = select_board_variants(variants, mode="cancer")

    assert len(selected) == 1
    assert selected[0]["gene"] == "KRAS"
    assert selected[0]["selection_reason"] == "Tier_III_hotspot"


def test_deep_intronic_clinvar_p_passes():
    """A deep-intronic ClinVar-Pathogenic splice variant MUST still be
    admitted via the unconditional P_LP bypass — the gate applies only
    to lower-confidence classifications.
    """
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "ATM",
            classification="Pathogenic",
            tier="Tier IV",
            hgvsp="",
            consequence="intron_variant",
        ),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="cancer")

    assert len(selected) == 1
    assert selected[0]["selection_reason"] == "P_LP"
    assert selected[0]["gene"] == "ATM"


def test_synonymous_spliceai_025_admitted():
    """A synonymous VUS with SpliceAI delta_max=0.25 is rescued by the
    splice-adjacent path and admitted via an existing must-reason
    (Tier III hotspot in this case).
    """
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "TP53",
            classification="VUS",
            tier="Tier III",
            hgvsp="p.Arg248Arg",
            consequence="synonymous_variant",
            in_silico={"spliceai_max": 0.25},
        ),
    ]
    with _patch_clinical(hotspot_positions={("TP53", 248)}):
        selected, meta = select_board_variants(variants, mode="cancer")

    assert len(selected) == 1
    assert selected[0]["selection_reason"] == "Tier_III_hotspot"


def test_synonymous_spliceai_015_rejected():
    """A synonymous VUS with SpliceAI delta_max=0.15 is BELOW the 0.20
    rescue threshold and must be rejected by the consequence gate.
    """
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "TP53",
            classification="VUS",
            tier="Tier III",
            hgvsp="p.Arg248Arg",
            consequence="synonymous_variant",
            in_silico={"spliceai_max": 0.15},
        ),
    ]
    with _patch_clinical(hotspot_positions={("TP53", 248)}):
        selected, meta = select_board_variants(variants, mode="cancer")

    assert selected == []
    assert meta["empty"] is True


def test_splice_region_rescued_at_exact_threshold():
    """Boundary: SpliceAI delta_max = 0.20 exactly IS admitted (>= threshold)."""
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "BRAF",
            classification="VUS",
            tier="Tier III",
            hgvsp="p.Val600Val",
            consequence="splice_region_variant",
            in_silico={"spliceai_max": 0.20},
        ),
    ]
    with _patch_clinical(hotspot_positions={("BRAF", 600)}):
        selected, _ = select_board_variants(variants, mode="cancer")

    assert len(selected) == 1
    assert selected[0]["selection_reason"] == "Tier_III_hotspot"


def test_cancer_may_arm_also_gated():
    """The MAY arm (VUS_hotspot, etc.) is gated by consequence too.
    A VUS at a hotspot position with intron_variant consequence must
    NOT be admitted via VUS_hotspot.
    """
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "BRAF",
            classification="VUS",
            tier="Tier IV",
            hgvsp="p.Val600Glu",
            consequence="intron_variant",
        ),
    ]
    with _patch_clinical(hotspot_positions={("BRAF", 600)}):
        selected, meta = select_board_variants(variants, mode="cancer")

    assert selected == []


