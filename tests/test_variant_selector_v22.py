"""v2.2 Phase B1/B2 regression — consequence gate + MMR/Lynch carve-out.

Scope:
- B1: Tier III (and Tier I/II) consequence gate rejects non-coding VUS.
- B1: P_LP branch bypasses the gate unconditionally.
- B1: SpliceAI >= 0.2 rescues synonymous / splice_region consequences.
- B2: MMR/Lynch VUS carve-out admits protein-impacting VUS in
  MLH1/MSH2/MSH6/PMS2/EPCAM regardless of hotspot / TSG / HPO.

All tests mock ``is_hotspot`` / ``is_cancer_gene`` to isolate the
selector logic from CIViC DB availability.
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


def test_gate_accepts_biko_formatted_consequence_labels():
    """v2.2 regression — the real pipeline stores consequence as the
    BIKO-formatted short label (e.g. "Missense") via
    scripts/intake/parse_annotation.py::format_consequence, not the raw VEP
    SO term. The B1 gate and downstream reason detection must match both
    forms or every Tier III MAY-arm VUS silently disappears from the Board.

    Root cause of the v2.2 post-Phase-B showcase regression — TP53 R249M
    was dropped from the real-VCF run because consequence="Missense"
    failed set membership against {"missense_variant", ...}.
    """
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "KRAS",
            classification="VUS",
            tier="Tier III",
            hgvsp="p.Gly12Asp",
            consequence="Missense",  # BIKO formatted form, not SO term
        ),
        _mk_snv(
            "TP53",
            classification="VUS",
            tier="Tier II",
            hgvsp="ENSP00000269305.4:p.Arg249Met",
            consequence="Missense",
        ),
    ]
    with _patch_clinical(
        hotspot_positions={("KRAS", 12), ("TP53", 249)},
    ):
        selected, meta = select_board_variants(variants, mode="cancer")

    genes = {v["gene"] for v in selected}
    assert "KRAS" in genes
    assert "TP53" in genes
    assert meta["selected"] == 2


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


# ---------------------------------------------------------------------------
# B2 — MMR / Lynch carve-out
# ---------------------------------------------------------------------------


def test_pms2_r563q_missense_vus_admitted_as_mmr_lynch():
    """PMS2 R563Q — missense VUS, no hotspot, no TSG_LoF — must be
    admitted with reason ``VUS_MMR_Lynch``. This is the F1-parity case
    where BIKO previously dropped the row.
    """
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "PMS2",
            classification="VUS",
            tier="Tier IV",
            hgvsp="p.Arg563Gln",
            consequence="missense_variant",
        ),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="cancer")

    assert len(selected) == 1
    assert selected[0]["gene"] == "PMS2"
    assert selected[0]["selection_reason"] == "VUS_MMR_Lynch"


def test_mlh1_intronic_vus_not_admitted_b1_b2_intersection():
    """B1 consequence gate still applies inside the MMR carve-out —
    an MLH1 intronic VUS is NOT admitted. There is no Lynch exception
    to the consequence gate.
    """
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "MLH1",
            classification="VUS",
            tier="Tier IV",
            hgvsp="",
            consequence="intron_variant",
        ),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="cancer")

    assert selected == []


def test_all_mmr_lynch_genes_covered():
    """All five MMR/Lynch panel genes admit protein-impacting VUS."""
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(g, classification="VUS", tier="Tier IV", hgvsp="p.Arg100Lys")
        for g in ("MLH1", "MSH2", "MSH6", "PMS2", "EPCAM")
    ]
    with _patch_clinical():
        selected, _ = select_board_variants(variants, mode="cancer")

    assert {v["gene"] for v in selected} == {
        "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM",
    }
    assert all(v["selection_reason"] == "VUS_MMR_Lynch" for v in selected)


def test_mmr_lynch_priority_between_hotspot_and_tsg_lof():
    """Ordering invariant: VUS_hotspot < VUS_MMR_Lynch < VUS_TSG_LoF."""
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "TP53",
            classification="VUS",
            tier="Tier IV",
            consequence="frameshift_variant",
            cancer_gene_type="TSG",
        ),
        _mk_snv(
            "PMS2",
            classification="VUS",
            tier="Tier IV",
            hgvsp="p.Arg563Gln",
            consequence="missense_variant",
        ),
        _mk_snv(
            "BRAF",
            classification="VUS",
            tier="Tier IV",
            hgvsp="p.Val600Glu",
            consequence="missense_variant",
        ),
    ]
    with _patch_clinical(hotspot_positions={("BRAF", 600)}):
        selected, _ = select_board_variants(variants, mode="cancer")

    order = [v["gene"] for v in selected]
    assert order.index("BRAF") < order.index("PMS2") < order.index("TP53")


def test_mmr_lynch_non_protein_consequence_rejected():
    """A PMS2 synonymous VUS (no SpliceAI signal) is NOT admitted —
    Lynch carve-out still respects the consequence gate.
    """
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "PMS2",
            classification="VUS",
            tier="Tier IV",
            hgvsp="p.Arg563Arg",
            consequence="synonymous_variant",
        ),
    ]
    with _patch_clinical():
        selected, _ = select_board_variants(variants, mode="cancer")

    assert selected == []


def test_mmr_lynch_hotspot_takes_higher_priority_reason():
    """A hotspot-positive MMR variant takes the higher-priority
    VUS_hotspot reason, not VUS_MMR_Lynch — MMR carve-out is below
    VUS_hotspot in the MAY branch order.
    """
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "MSH2",
            classification="VUS",
            tier="Tier IV",
            hgvsp="p.Arg200Gln",
            consequence="missense_variant",
        ),
    ]
    with _patch_clinical(hotspot_positions={("MSH2", 200)}):
        selected, _ = select_board_variants(variants, mode="cancer")

    assert len(selected) == 1
    assert selected[0]["selection_reason"] == "VUS_hotspot"


