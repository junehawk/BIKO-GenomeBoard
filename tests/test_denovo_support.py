"""Tests for v1 de novo variant support (PS2/PM6 + rare-disease carve-out).

Spec: ``_workspace/v23-engineering/00_clinical_denovo_spec.md``
Implementer checklist covers models, parse_vcf, evidence_collector, ddg2p
panel, variant_selector, and the narrow ClinVar-override sanity check.

Fixture notes
-------------
- Uses the BIKO short-label consequence form (``"Missense"``, ``"Intronic"``)
  where the real pipeline would surface it via ``format_consequence``. Tests
  that pass raw SO terms directly would be silent dead code on real data
  (CLAUDE.md "consequence-form normalization" pitfall — commit 564f5da).
- The test variants stay under ``chr17:7675088`` (an arbitrary non-hotspot
  TP53 position) unless a specific gene signal is required.
"""

from __future__ import annotations

import pytest

from scripts.classification.evidence_collector import (
    collect_additional_evidence,
    collect_denovo_evidence,
)
from scripts.classification.acmg_engine import apply_hotspot_conflict_reconciliation, ClassificationResult
from scripts.common.ddg2p_panel import is_admitted_neurodev_gene, get_neurodev_info
from scripts.common.models import Variant


# ---------------------------------------------------------------------------
# 1–3. parse_vcf INFO flag plumbing
# ---------------------------------------------------------------------------


_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    '##INFO=<ID=DN,Number=0,Type=Flag,Description="De novo">\n'
    '##INFO=<ID=CONFIRMED_DN,Number=0,Type=Flag,Description="Trio-confirmed de novo">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _write_vcf(tmp_path, info: str):
    body = f"chr2\t166000000\t.\tA\tT\t50\tPASS\t{info}\n"
    path = tmp_path / "sample.vcf"
    path.write_text(_VCF_HEADER + body)
    return str(path)


def test_parse_vcf_denovo_flag_sets_inheritance(tmp_path):
    from scripts.intake.parse_vcf import parse_vcf

    variants = parse_vcf(_write_vcf(tmp_path, "Gene=SCN1A;DN=1"))
    assert len(variants) == 1
    v = variants[0]
    assert v.inheritance == "de_novo"
    assert v.confirmed_denovo is False


def test_parse_vcf_confirmed_denovo_flag(tmp_path):
    from scripts.intake.parse_vcf import parse_vcf

    variants = parse_vcf(_write_vcf(tmp_path, "Gene=SCN1A;CONFIRMED_DN=1"))
    v = variants[0]
    assert v.confirmed_denovo is True
    assert v.inheritance == "confirmed_de_novo"


def test_parse_vcf_ps2_flag_also_confirms(tmp_path):
    """Upstream tool tagging with PS2=1 is honoured as confirmed."""
    from scripts.intake.parse_vcf import parse_vcf

    variants = parse_vcf(_write_vcf(tmp_path, "Gene=SCN1A;PS2=1"))
    v = variants[0]
    assert v.confirmed_denovo is True


def test_parse_vcf_no_denovo_flag_stays_none(tmp_path):
    from scripts.intake.parse_vcf import parse_vcf

    variants = parse_vcf(_write_vcf(tmp_path, "Gene=SCN1A"))
    v = variants[0]
    assert v.inheritance is None
    assert v.confirmed_denovo is False


def test_parse_vcf_regression_basic_info_still_parses(tmp_path):
    """Regression: adding DN parsing must not break existing Gene= parsing."""
    from scripts.intake.parse_vcf import parse_vcf

    variants = parse_vcf(_write_vcf(tmp_path, "Gene=TP53;AF=0.01"))
    assert variants[0].gene == "TP53"
    assert variants[0].inheritance is None


# ---------------------------------------------------------------------------
# 4–7. evidence_collector de novo path
# ---------------------------------------------------------------------------


def _denovo_variant(**kwargs) -> Variant:
    defaults = dict(
        chrom="chr2",
        pos=166000000,
        ref="A",
        alt="T",
        gene="SCN1A",
        consequence="Missense",
        hgvsp="p.Arg101Leu",
        inheritance="de_novo",
    )
    defaults.update(kwargs)
    return Variant(**defaults)


def test_collect_denovo_evidence_pm6():
    v = _denovo_variant()
    assert collect_denovo_evidence(v) == ["PM6"]


def test_collect_denovo_evidence_ps2_moderate():
    v = _denovo_variant(confirmed_denovo=True, inheritance="confirmed_de_novo")
    assert collect_denovo_evidence(v) == ["PS2_Moderate"]


def test_collect_denovo_evidence_consequence_gate_intronic():
    """Spec Q4 collision point 2: intronic de novo without SpliceAI → no code."""
    v = _denovo_variant(consequence="Intronic", hgvsp=None)
    assert collect_denovo_evidence(v) == []


def test_collect_denovo_evidence_consequence_gate_synonymous():
    v = _denovo_variant(consequence="Synonymous", hgvsp=None)
    assert collect_denovo_evidence(v) == []


def test_collect_denovo_evidence_splice_rescue_intronic():
    """Intronic de novo rescued by SpliceAI >= 0.2 fires PM6."""
    v = _denovo_variant(consequence="Intronic", hgvsp=None)
    v.in_silico = {"spliceai_max": 0.3}
    assert collect_denovo_evidence(v) == ["PM6"]


def test_collect_denovo_evidence_splice_rescue_synonymous():
    v = _denovo_variant(consequence="Synonymous", hgvsp=None)
    v.in_silico = {"spliceai_max": 0.25}
    assert collect_denovo_evidence(v) == ["PM6"]


def test_collect_denovo_evidence_returns_empty_without_flag():
    v = _denovo_variant(inheritance=None)
    assert collect_denovo_evidence(v) == []


def test_collect_denovo_evidence_frameshift_confirmed():
    v = _denovo_variant(consequence="Frameshift", confirmed_denovo=True, inheritance="confirmed_de_novo")
    assert collect_denovo_evidence(v) == ["PS2_Moderate"]


# ---------------------------------------------------------------------------
# 8–9. DDG2P panel loader
# ---------------------------------------------------------------------------


def test_ddg2p_panel_admits_scn1a():
    assert is_admitted_neurodev_gene("SCN1A") is True
    info = get_neurodev_info("SCN1A")
    assert info is not None
    assert info["confidence"] in ("definitive", "strong", "moderate")


def test_ddg2p_panel_rejects_unknown():
    assert is_admitted_neurodev_gene("DOES_NOT_EXIST_FAKE_GENE") is False
    assert get_neurodev_info("DOES_NOT_EXIST_FAKE_GENE") is None


def test_ddg2p_panel_handles_none_gene():
    assert is_admitted_neurodev_gene(None) is False
    assert is_admitted_neurodev_gene("") is False


# ---------------------------------------------------------------------------
# 10–12. variant_selector rare-disease carve-out
# ---------------------------------------------------------------------------


def _rare_variant(**kwargs) -> dict:
    defaults = {
        "gene": "SCN1A",
        "classification": "VUS",
        "consequence": "Missense",
        "hgvsp": "p.Arg101Leu",
        "variant_inheritance": "de_novo",
        "confirmed_denovo": False,
        "hpo_score": 0,
        "gnomad_af": 0.0,
        "in_silico": {},
        "tier": "",
        "cancer_gene_type": "",
        "oncokb_level": "",
    }
    defaults.update(kwargs)
    return defaults


def test_selector_admits_denovo_neurodev_missense_vus():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [_rare_variant()]
    selected, meta = select_board_variants(variants, mode="rare-disease")
    assert len(selected) == 1
    assert selected[0]["selection_reason"] == "VUS_denovo_neurodev"
    assert meta["selected"] == 1


def test_selector_rejects_benign_denovo_synonymous():
    """B1 gate: de novo synonymous without SpliceAI is rejected."""
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [_rare_variant(consequence="Synonymous")]
    selected, meta = select_board_variants(variants, mode="rare-disease")
    assert len(selected) == 0


def test_selector_rescues_denovo_splice_with_spliceai():
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [_rare_variant(consequence="Synonymous", in_silico={"spliceai_max": 0.3})]
    selected, _ = select_board_variants(variants, mode="rare-disease")
    assert len(selected) == 1
    assert selected[0]["selection_reason"] == "VUS_denovo_splice"


def test_selector_rescues_denovo_intronic_with_spliceai():
    """Spec Q3.2: intronic de novo + SpliceAI >= 0.2 opens narrow window."""
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [_rare_variant(consequence="Intronic", in_silico={"spliceai_max": 0.25})]
    selected, _ = select_board_variants(variants, mode="rare-disease")
    assert len(selected) == 1
    assert selected[0]["selection_reason"] == "VUS_denovo_splice"


def test_selector_rejects_denovo_unknown_gene():
    """De novo in a non-DDG2P, unconstrained gene is rejected."""
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [_rare_variant(gene="FAKE_UNLISTED_GENE")]
    selected, _ = select_board_variants(variants, mode="rare-disease")
    # Falls through to HPO path, which requires hpo_score >= 1 → excluded.
    assert len(selected) == 0


def test_selector_denovo_with_hpo_records_both_reasons():
    """Collision point 3: when HPO also matches, both reasons are recorded."""
    from scripts.clinical_board.variant_selector import select_board_variants

    variants = [_rare_variant(hpo_score=3)]
    selected, _ = select_board_variants(variants, mode="rare-disease")
    assert len(selected) == 1
    assert selected[0]["selection_reason"] == "VUS_denovo_neurodev"
    # Auxiliary HPO signal is preserved for the board renderer.
    assert selected[0].get("selection_reason_list") == ["VUS_denovo_neurodev", "VUS_HPO_match"]


# ---------------------------------------------------------------------------
# 13. ClinVar override sanity — PS2/PM6 must NOT participate
# ---------------------------------------------------------------------------


def test_no_ps2_in_clinvar_override_conditions():
    """Spec Q4 collision point 4: PS2/PM6 alone must not trip the narrow
    ClinVar conflict override. The override requires PM1 + PM5; a LP result
    built on PS2_Moderate + PM2 + PP3 must leave clinvar_override_reason
    empty."""
    result = ClassificationResult(
        classification="Likely Pathogenic",
        evidence_codes=["PS2_Moderate", "PM2_Supporting", "PP3"],
    )
    apply_hotspot_conflict_reconciliation(
        result,
        clinvar_significance="Conflicting interpretations of pathogenicity",
        gene="SCN1A",
        hgvsp="p.Arg101Leu",
    )
    assert result.clinvar_override_reason == ""


# ---------------------------------------------------------------------------
# 14–15. Coexistence + pass-through sanity
# ---------------------------------------------------------------------------


def test_denovo_fires_pm2_supporting_coexist():
    """Spec Q4 collision point 1: PM6 + PM2_Supporting co-firing allowed.

    collect_additional_evidence and collect_denovo_evidence are independent
    — neither short-circuits the other, so a de novo missense with a
    missense-Z-scored constraint gene gets PM6 + PP2 stamped simultaneously.
    """
    v = _denovo_variant()
    denovo_codes = collect_denovo_evidence(v)
    extra_codes = collect_additional_evidence(v, gene_info={"missense_z": 3.5})
    all_codes = denovo_codes + extra_codes
    assert "PM6" in all_codes
    assert "PP2" in all_codes


def test_inherited_variant_unchanged(tmp_path):
    """Variant without any de novo flag is parsed + classified exactly as
    before the v1 de novo support landed."""
    from scripts.intake.parse_vcf import parse_vcf

    variants = parse_vcf(_write_vcf(tmp_path, "Gene=TP53"))
    v = variants[0]
    assert v.inheritance is None
    assert v.confirmed_denovo is False
    # Trio arm is dormant → evidence collector returns no de novo codes.
    v.consequence = "Missense"
    assert collect_denovo_evidence(v) == []


@pytest.mark.parametrize(
    "confidence,expected",
    [
        ("SCN1A", True),
        ("MECP2", True),
        ("STXBP1", True),
        ("TSC2", True),
    ],
)
def test_ddg2p_well_known_neurodev_genes_admitted(confidence, expected):
    """Sanity: the v1 curated panel covers the well-known NDD genes."""
    assert is_admitted_neurodev_gene(confidence) is expected
