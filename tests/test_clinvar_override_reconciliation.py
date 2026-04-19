"""A4-core: narrow ClinVar-conflict override when PM1 + PM5 both fire.

Covers the five clinically blessed fixtures from
``_workspace/v22-phaseA/artifacts/00_clinical_review.md`` §3. Exactly one of
the five fires the override — the other four exercise the four negative
paths around the gate.
"""

from __future__ import annotations

import pytest

from scripts.classification.acmg_engine import (
    ClassificationResult,
    _is_conflicting_status,
    apply_hotspot_conflict_reconciliation,
)


def _lp_result(codes):
    return ClassificationResult(
        classification="Likely Pathogenic",
        evidence_codes=list(codes),
        matched_rule="likely_pathogenic_0",
    )


def _p_result(codes):
    return ClassificationResult(
        classification="Pathogenic",
        evidence_codes=list(codes),
        matched_rule="pathogenic_0",
    )


def _vus_result(codes):
    return ClassificationResult(
        classification="VUS",
        evidence_codes=list(codes),
        matched_rule="",
    )


# ── _is_conflicting_status prefix/category match ────────────────────────────


@pytest.mark.parametrize(
    "status",
    [
        "Conflicting classifications of pathogenicity",
        "Conflicting interpretations of pathogenicity",
        "conflicting_pathogenicity",
        "  Conflicting classifications of pathogenicity  ",
        "CONFLICTING classifications of pathogenicity",
    ],
)
def test_is_conflicting_status_matches(status):
    assert _is_conflicting_status(status) is True


@pytest.mark.parametrize(
    "status",
    ["", None, "Pathogenic", "Benign", "not provided", "Likely pathogenic"],
)
def test_is_conflicting_status_rejects(status):
    assert _is_conflicting_status(status) is False


# ── Fixture 1 — TP53 R249M + conflicting + PM1 + PM5 → override fires ──────


def test_fixture_1_tp53_r249m_override_fires():
    result = _lp_result(["PM1", "PM5", "PM2_Supporting", "PP3"])
    out = apply_hotspot_conflict_reconciliation(
        result,
        clinvar_significance="Conflicting classifications of pathogenicity",
        gene="TP53",
        hgvsp="p.Arg249Met",
    )
    assert out.classification == "Likely Pathogenic"
    assert out.clinvar_override_reason, "override_reason must be populated"
    assert "PM1" in out.clinvar_override_reason
    assert "PM5" in out.clinvar_override_reason
    assert "TP53" in out.clinvar_override_reason
    assert "PMID" in out.clinvar_override_reason
    assert len(out.clinvar_override_reason) <= 200


def test_fixture_1_preserves_p_when_engine_reaches_p():
    result = _p_result(["PM1", "PM5", "PS3", "PP3"])
    out = apply_hotspot_conflict_reconciliation(
        result,
        clinvar_significance="Conflicting classifications of pathogenicity",
        gene="TP53",
        hgvsp="p.Arg249Met",
    )
    assert out.classification == "Pathogenic"
    assert "engine P override" in out.clinvar_override_reason


# ── Fixture 2 — TP53 R100Q + conflicting + no PM1/PM5 → VUS, no override ───


def test_fixture_2_tp53_r100q_no_override():
    result = _vus_result(["PM2_Supporting"])
    out = apply_hotspot_conflict_reconciliation(
        result,
        clinvar_significance="Conflicting classifications of pathogenicity",
        gene="TP53",
        hgvsp="p.Arg100Gln",
    )
    assert out.classification == "VUS"
    assert out.clinvar_override_reason == ""


# ── Fixture 3 — BRCA1 VUS + ClinVar Benign → no A4 override ────────────────


def test_fixture_3_brca1_benign_clinvar_wins():
    result = _vus_result([])
    out = apply_hotspot_conflict_reconciliation(
        result,
        clinvar_significance="Benign",
        gene="BRCA1",
        hgvsp="p.Glu1134Lys",
    )
    assert out.classification == "VUS"
    assert out.clinvar_override_reason == ""


# ── Fixture 4 — KRAS G12D + ClinVar Pathogenic → no A4 override ────────────


def test_fixture_4_kras_g12d_clinvar_pathogenic_no_override():
    result = _lp_result(["PM1", "PM2_Supporting", "PP3"])
    out = apply_hotspot_conflict_reconciliation(
        result,
        clinvar_significance="Pathogenic",
        gene="KRAS",
        hgvsp="p.Gly12Asp",
    )
    assert out.classification == "Likely Pathogenic"
    assert out.clinvar_override_reason == ""


# ── Fixture 5 — TP53 R175H + ClinVar Pathogenic → no A4 override ───────────


def test_fixture_5_tp53_r175h_clinvar_pathogenic_no_override():
    result = _lp_result(["PM1", "PM5", "PM2_Supporting", "PP3", "PS3_Supporting"])
    out = apply_hotspot_conflict_reconciliation(
        result,
        clinvar_significance="Pathogenic",
        gene="TP53",
        hgvsp="p.Arg175His",
    )
    assert out.classification == "Likely Pathogenic"
    assert out.clinvar_override_reason == ""


# ── Config gate ─────────────────────────────────────────────────────────────


def test_config_gate_disables_override(monkeypatch):
    from scripts.classification import acmg_engine

    monkeypatch.setattr(
        acmg_engine,
        "_override_enabled",
        lambda: False,
    )
    result = _lp_result(["PM1", "PM5"])
    out = apply_hotspot_conflict_reconciliation(
        result,
        clinvar_significance="Conflicting classifications of pathogenicity",
        gene="TP53",
        hgvsp="p.Arg249Met",
    )
    assert out.clinvar_override_reason == ""


# ── Schema delta ────────────────────────────────────────────────────────────


def test_classification_result_has_override_reason_field():
    import dataclasses as dc

    field_names = {f.name for f in dc.fields(ClassificationResult)}
    assert "clinvar_override_reason" in field_names
    result = ClassificationResult(classification="VUS")
    assert result.clinvar_override_reason == ""


# ── Engine must have reached P/LP before override considers firing ─────────


def test_engine_vus_does_not_get_upgraded_even_with_pm1_pm5():
    """Condition 1: engine must independently reach P/LP pre-override."""
    result = _vus_result(["PM1", "PM5"])
    out = apply_hotspot_conflict_reconciliation(
        result,
        clinvar_significance="Conflicting classifications of pathogenicity",
        gene="TP53",
        hgvsp="p.Arg249Met",
    )
    assert out.classification == "VUS"
    assert out.clinvar_override_reason == ""


def test_missing_pm1_blocks_override():
    result = _lp_result(["PM5", "PM2_Supporting", "PP3", "PS3_Supporting"])
    out = apply_hotspot_conflict_reconciliation(
        result,
        clinvar_significance="Conflicting classifications of pathogenicity",
        gene="TP53",
        hgvsp="p.Arg249Met",
    )
    assert out.clinvar_override_reason == ""


def test_missing_pm5_blocks_override():
    result = _lp_result(["PM1", "PM2_Supporting", "PP3", "PS3_Supporting"])
    out = apply_hotspot_conflict_reconciliation(
        result,
        clinvar_significance="Conflicting classifications of pathogenicity",
        gene="TP53",
        hgvsp="p.Arg249Met",
    )
    assert out.clinvar_override_reason == ""
