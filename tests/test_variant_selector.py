"""Tests for scripts.clinical_board.variant_selector.

Implements the 25-test spec from Task 1. Validates the 5 clinical-review
corrections from _workspace/variant-selector/00_clinical_review.md:

1. Tier III hotspot + OncoKB-gene MUST-include
2. "protein-impacting alone" rejected as VUS admission arm
3. No top-3 fallback — empty selection is legitimate
4. Soft caps — never truncate MUST items
5. ACMG SF v3.2 VUS silently excluded in rare-disease mode
"""

from __future__ import annotations

from unittest.mock import patch


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


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
    variant_type: str = "SNV",
    tier_evidence_source: str = "",
) -> dict:
    return {
        "gene": gene,
        "variant": f"{gene}_var",
        "classification": classification,
        "tier": tier,
        "tier_evidence_source": tier_evidence_source,
        "hgvsp": hgvsp,
        "consequence": consequence,
        "cancer_gene_type": cancer_gene_type,
        "oncokb_level": oncokb_level,
        "hpo_score": hpo_score,
        "gnomad_af": gnomad_af,
        "variant_type": variant_type,
    }


def _patch_clinical(*, hotspot_positions=None, cancer_genes=None):
    """Patch is_hotspot and is_cancer_gene used inside variant_selector."""
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
# Cancer mode — MUST include
# ---------------------------------------------------------------------------


def test_cancer_selects_all_p_lp():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv("TP53", classification="Pathogenic"),
        _mk_snv("BRCA1", classification="Likely Pathogenic"),
        _mk_snv("FOO", classification="VUS"),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="cancer")

    genes = [v["gene"] for v in selected]
    assert "TP53" in genes
    assert "BRCA1" in genes
    reasons = {v["gene"]: v["selection_reason"] for v in selected}
    assert reasons["TP53"] == "P_LP"
    assert reasons["BRCA1"] == "P_LP"


def test_cancer_selects_all_tier1_tier2():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv("EGFR", classification="VUS", tier="Tier I"),
        _mk_snv("KRAS", classification="VUS", tier="Tier II"),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="cancer")

    assert {v["gene"] for v in selected} == {"EGFR", "KRAS"}
    reasons = {v["gene"]: v["selection_reason"] for v in selected}
    assert reasons["EGFR"] == "Tier_I"
    assert reasons["KRAS"] == "Tier_II"


def test_cancer_tier3_hotspot_must_include():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "BRAF",
            classification="VUS",
            tier="Tier III",
            hgvsp="p.Val600Glu",
        ),
    ]
    with _patch_clinical(hotspot_positions={("BRAF", 600)}):
        selected, meta = select_board_variants(variants, mode="cancer")

    assert len(selected) == 1
    assert selected[0]["selection_reason"] == "Tier_III_hotspot"
    assert meta["must_included"] == 1


def test_cancer_tier3_oncokb_gene_must_include():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "KRAS",
            classification="VUS",
            tier="Tier III",
            hgvsp="p.Gly13Asp",
            oncokb_level="1",
        ),
    ]
    with _patch_clinical(cancer_genes={"KRAS"}):
        selected, meta = select_board_variants(variants, mode="cancer")

    assert len(selected) == 1
    assert selected[0]["selection_reason"] == "Tier_III_oncokb_gene"
    assert meta["must_included"] == 1


# ---------------------------------------------------------------------------
# Cancer mode — MAY include VUS
# ---------------------------------------------------------------------------


def test_cancer_may_vus_hotspot():
    from scripts.clinical_board.variant_selector import select_board_variants

    # 15 hotspot VUS → cap at 10
    variants = [
        _mk_snv(
            f"GENE{i}",
            classification="VUS",
            tier="Tier IV",
            hgvsp=f"p.Arg{100 + i}Lys",
        )
        for i in range(15)
    ]
    hotspots = {(f"GENE{i}", 100 + i) for i in range(15)}
    with _patch_clinical(hotspot_positions=hotspots):
        selected, meta = select_board_variants(variants, mode="cancer")

    may = [v for v in selected if v["selection_reason"] == "VUS_hotspot"]
    assert len(may) == 10
    assert meta["may_included"] == 10
    assert meta["truncated"] is True
    assert meta["n_dropped"] == 5


def test_cancer_may_vus_tsg_lof():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "TP53",
            classification="VUS",
            tier="Tier IV",
            consequence="frameshift_variant",
            cancer_gene_type="TSG",
        ),
        # inframe — not LoF, should NOT match TSG_LoF arm
        _mk_snv(
            "APC",
            classification="VUS",
            tier="Tier IV",
            consequence="inframe_deletion",
            cancer_gene_type="TSG",
        ),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="cancer")

    genes = {v["gene"] for v in selected}
    assert "TP53" in genes
    assert "APC" not in genes
    tp53 = next(v for v in selected if v["gene"] == "TP53")
    assert tp53["selection_reason"] == "VUS_TSG_LoF"


def test_cancer_rejects_protein_impacting_alone():
    from scripts.clinical_board.variant_selector import select_board_variants

    # missense VUS in random non-driver gene — no hotspot, no TSG — excluded
    variants = [
        _mk_snv(
            "RANDOM1",
            classification="VUS",
            tier="Tier IV",
            consequence="missense_variant",
            hgvsp="p.Arg42Lys",
        ),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="cancer")

    assert selected == []
    assert meta["empty"] is True


def test_cancer_soft_cap_never_truncates_must():
    from scripts.clinical_board.variant_selector import select_board_variants

    # 35 Tier I variants → soft cap 30, but MUST never truncated
    variants = [_mk_snv(f"GENE{i}", classification="VUS", tier="Tier I") for i in range(35)]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="cancer")

    assert len(selected) == 35
    assert meta["must_included"] == 35
    assert meta["hard_cap_applied"] is False


def test_cancer_excludes_tier4_benign_lb():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv("G1", classification="Benign", tier="Tier IV"),
        _mk_snv("G2", classification="Likely Benign", tier="Tier IV"),
        _mk_snv("G3", classification="VUS", tier="Tier IV"),  # no driver signal
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="cancer")

    assert selected == []
    assert meta["excluded"] >= 3


def test_cancer_tmb_high_footnote():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [_mk_snv("TP53", classification="Pathogenic")]
    report_data = {"tmb": {"score": 15.0}}
    with _patch_clinical():
        _, meta = select_board_variants(variants, mode="cancer", report_data=report_data)

    assert meta["tmb_high_footnote"] is True


def test_cancer_cnv_sv_class4_must_include():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        {
            "gene": "MYCN",
            "variant": "chr2:16080543",
            "classification": "Pathogenic",
            "tier": "Tier IV",
            "hgvsp": "",
            "consequence": "",
            "cancer_gene_type": "",
            "oncokb_level": "",
            "hpo_score": 0,
            "gnomad_af": None,
            "variant_type": "CNV",
            "tier_evidence_source": "",
        },
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="cancer")

    assert len(selected) == 1
    assert selected[0]["variant_type"] == "CNV"
    assert selected[0]["selection_reason"] == "P_LP"
    assert meta["must_included"] == 1


def test_cancer_pgx_excluded():
    from scripts.clinical_board.variant_selector import select_board_variants

    # CYP2C19 with no clinical signal — should be excluded even if consequence
    # is missense (no cancer criterion matches)
    variants = [
        _mk_snv(
            "CYP2C19",
            classification="Drug Response",
            tier="Tier IV",
            consequence="missense_variant",
        ),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="cancer")

    genes = {v["gene"] for v in selected}
    assert "CYP2C19" not in genes


def test_cancer_empty_selection_sets_empty_reason():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv("G1", classification="VUS", tier="Tier IV"),
        _mk_snv("G2", classification="Benign"),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="cancer")

    assert selected == []
    assert meta["empty"] is True
    assert "Tier" in meta["empty_reason"] or "hotspot" in meta["empty_reason"]


# ---------------------------------------------------------------------------
# Rare-disease mode
# ---------------------------------------------------------------------------


def test_rare_all_p_lp_included():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv("G1", classification="Pathogenic"),
        _mk_snv("G2", classification="Likely Pathogenic"),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="rare-disease")

    assert {v["gene"] for v in selected} == {"G1", "G2"}
    assert all(v["selection_reason"] == "P_LP" for v in selected)


def test_rare_vus_requires_hpo_match_and_rare_af():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        # HPO 0 → excluded
        _mk_snv("G1", classification="VUS", hpo_score=0, gnomad_af=0.0001),
        # HPO≥1 + rare AF → included
        _mk_snv(
            "G2",
            classification="VUS",
            hpo_score=2,
            gnomad_af=0.0001,
            consequence="missense_variant",
        ),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="rare-disease")

    genes = {v["gene"] for v in selected}
    assert "G1" not in genes
    assert "G2" in genes


def test_rare_vus_excludes_common_af():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "G1",
            classification="VUS",
            hpo_score=3,
            gnomad_af=0.05,
            consequence="missense_variant",
        ),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="rare-disease")

    assert selected == []


def test_rare_vus_excludes_sf_gene():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "BRCA1",
            classification="VUS",
            hpo_score=5,
            gnomad_af=0.0001,
            consequence="missense_variant",
        ),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="rare-disease")

    genes = {v["gene"] for v in selected}
    assert "BRCA1" not in genes


def test_rare_may_capped_at_10():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            f"GENE{i}",
            classification="VUS",
            hpo_score=1,
            gnomad_af=0.0001,
            consequence="missense_variant",
        )
        for i in range(15)
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="rare-disease")

    assert len([v for v in selected if v["selection_reason"] == "VUS_HPO_match"]) == 10
    assert meta["truncated"] is True
    assert meta["n_dropped"] == 5


def test_rare_soft_cap_never_truncates_p_lp():
    from scripts.clinical_board.variant_selector import select_board_variants

    # 25 P/LP → exceeds default cap 20, but MUST never truncated
    variants = [_mk_snv(f"G{i}", classification="Pathogenic") for i in range(25)]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="rare-disease")

    assert len(selected) == 25
    assert meta["must_included"] == 25


def test_rare_empty_sets_empty_reason():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv("G1", classification="Benign"),
        _mk_snv("G2", classification="VUS", hpo_score=0),
    ]
    with _patch_clinical():
        selected, meta = select_board_variants(variants, mode="rare-disease")

    assert selected == []
    assert meta["empty"] is True
    assert "HPO" in meta["empty_reason"] or "P/LP" in meta["empty_reason"]


# ---------------------------------------------------------------------------
# Cross-cutting
# ---------------------------------------------------------------------------


def test_metadata_fields_present():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [_mk_snv("TP53", classification="Pathogenic")]
    with _patch_clinical():
        _, meta = select_board_variants(variants, mode="cancer")

    required = {
        "mode",
        "total_input",
        "selected",
        "must_included",
        "may_included",
        "excluded",
        "truncated",
        "n_dropped",
        "hard_cap_applied",
        "empty",
        "empty_reason",
        "tmb_high_footnote",
        "criteria_summary",
        "by_selection_reason",
    }
    assert required.issubset(meta.keys())


def test_selection_reason_propagated():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv("TP53", classification="Pathogenic"),
        _mk_snv("EGFR", classification="VUS", tier="Tier I"),
        _mk_snv(
            "KRAS",
            classification="VUS",
            tier="Tier IV",
            hgvsp="p.Gly12Asp",
        ),
    ]
    with _patch_clinical(hotspot_positions={("KRAS", 12)}):
        selected, _ = select_board_variants(variants, mode="cancer")

    assert all("selection_reason" in v for v in selected)


def test_ordering_follows_amp_table2():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            "HOTSPOT_VUS",
            classification="VUS",
            tier="Tier IV",
            hgvsp="p.Arg42Lys",
        ),
        _mk_snv("TIER2", classification="VUS", tier="Tier II"),
        _mk_snv("PATH", classification="Pathogenic"),
        _mk_snv("TIER1", classification="VUS", tier="Tier I"),
    ]
    with _patch_clinical(hotspot_positions={("HOTSPOT_VUS", 42)}):
        selected, _ = select_board_variants(variants, mode="cancer")

    order = [v["gene"] for v in selected]
    assert order.index("PATH") < order.index("TIER1")
    assert order.index("TIER1") < order.index("TIER2")
    assert order.index("TIER2") < order.index("HOTSPOT_VUS")


def test_config_overrides():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [
        _mk_snv(
            f"GENE{i}",
            classification="VUS",
            tier="Tier IV",
            hgvsp=f"p.Arg{100 + i}Lys",
        )
        for i in range(5)
    ]
    hotspots = {(f"GENE{i}", 100 + i) for i in range(5)}
    overrides = {"max_cancer_may_include": 2}
    with _patch_clinical(hotspot_positions=hotspots):
        selected, meta = select_board_variants(variants, mode="cancer", config_overrides=overrides)

    assert len([v for v in selected if v["selection_reason"] == "VUS_hotspot"]) == 2
    assert meta["n_dropped"] == 3


def test_reads_config_overrides(monkeypatch):
    """Config keys under clinical_board.variant_selection propagate to the selector."""
    from scripts.clinical_board import variant_selector as vs

    calls: list[tuple] = []

    def _fake_get(key, default=None):
        calls.append((key, default))
        if key == "clinical_board.variant_selection.max_cancer_may_include":
            return 3
        if key == "clinical_board.variant_selection.max_cancer_board_variants":
            return 30
        return default

    monkeypatch.setattr(vs, "get", _fake_get)

    variants = [
        _mk_snv(
            f"GENE{i}",
            classification="VUS",
            tier="Tier IV",
            hgvsp=f"p.Arg{100 + i}Lys",
        )
        for i in range(8)
    ]
    hotspots = {(f"GENE{i}", 100 + i) for i in range(8)}
    with _patch_clinical(hotspot_positions=hotspots):
        selected, meta = vs.select_board_variants(variants, mode="cancer")

    may = [v for v in selected if v["selection_reason"] == "VUS_hotspot"]
    assert len(may) == 3  # cap from fake config
    assert meta["n_dropped"] == 5
    # At least one call asked for the relevant config key
    keys_asked = {c[0] for c in calls}
    assert "clinical_board.variant_selection.max_cancer_may_include" in keys_asked


def test_fallback_top3_NOT_implemented():
    from scripts.clinical_board.variant_selector import select_board_variants

    # Empty input → empty output — no fabricated top-3
    with _patch_clinical():
        selected, meta = select_board_variants([], mode="cancer")

    assert selected == []
    assert meta["empty"] is True
    assert meta["total_input"] == 0
