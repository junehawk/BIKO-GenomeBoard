# tests/test_acmg_logic.py
from scripts.classification.acmg_engine import (
    _count_by_strength,
    check_clinvar_conflict,
    classify_variant,
)
from scripts.common.models import AcmgEvidence


def _ev(code, src="test"):
    return AcmgEvidence(code=code, source=src, description="test")


def test_pathogenic_pvs1_ps1():
    codes = [_ev("PVS1"), _ev("PS1")]
    result = classify_variant(codes)
    assert result.classification == "Pathogenic"


def test_pathogenic_two_ps():
    codes = [_ev("PS1"), _ev("PS3")]
    result = classify_variant(codes)
    assert result.classification == "Pathogenic"


def test_pathogenic_pvs1_two_pm():
    codes = [_ev("PVS1"), _ev("PM1"), _ev("PM4")]
    result = classify_variant(codes)
    assert result.classification == "Pathogenic"


def test_likely_pathogenic_pvs1_pm():
    codes = [_ev("PVS1"), _ev("PM1")]
    result = classify_variant(codes)
    assert result.classification == "Likely Pathogenic"


def test_likely_pathogenic_three_pm():
    codes = [_ev("PM1"), _ev("PM4"), _ev("PM5")]
    result = classify_variant(codes)
    assert result.classification == "Likely Pathogenic"


def test_benign_ba1():
    codes = [_ev("BA1")]
    result = classify_variant(codes)
    assert result.classification == "Benign"


def test_benign_two_bs():
    codes = [_ev("BS1"), _ev("BS2")]
    result = classify_variant(codes)
    assert result.classification == "Benign"


def test_likely_benign():
    codes = [_ev("BS1"), _ev("BP1")]
    result = classify_variant(codes)
    assert result.classification == "Likely Benign"


def test_vus_insufficient():
    codes = [_ev("PM1")]
    result = classify_variant(codes)
    assert result.classification == "VUS"


def test_vus_conflicting():
    """Pathogenic + Benign evidence both present → VUS"""
    codes = [_ev("PVS1"), _ev("PS1"), _ev("BS1"), _ev("BS2")]
    result = classify_variant(codes)
    assert result.classification == "VUS"
    assert result.conflict is True


def test_empty_codes():
    result = classify_variant([])
    assert result.classification == "VUS"


def test_result_has_evidence_summary():
    codes = [_ev("PVS1"), _ev("PS1")]
    result = classify_variant(codes)
    assert "PVS1" in result.evidence_codes
    assert "PS1" in result.evidence_codes


def test_pgx_bypass():
    """PGx genes should return 'Drug Response', not an ACMG classification."""
    codes = [_ev("PVS1"), _ev("PS1")]
    result = classify_variant(codes, gene="CYP2D6")
    assert result.classification == "Drug Response"
    assert result.matched_rule == "pgx_bypass"
    assert result.conflict is False


def test_pgx_bypass_cyp2c19():
    codes = [_ev("PM2")]
    result = classify_variant(codes, gene="CYP2C19")
    assert result.classification == "Drug Response"


def test_pgx_bypass_hla_b():
    codes = []
    result = classify_variant(codes, gene="HLA-B")
    assert result.classification == "Drug Response"


def test_non_pgx_gene_not_bypassed():
    """Non-PGx gene should still go through normal ACMG classification."""
    codes = [_ev("PVS1"), _ev("PS1")]
    result = classify_variant(codes, gene="TP53")
    assert result.classification == "Pathogenic"


def test_pm2_supporting_counts_as_pp():
    """PM2_Supporting should be counted as PP (supporting), not PM (moderate)."""
    codes = [_ev("PM2_Supporting")]
    result = classify_variant(codes)
    # One PP alone is VUS, not Likely Pathogenic (which would require PM)
    assert result.classification == "VUS"


def test_pm2_supporting_not_pm():
    """PM2_Supporting must not inflate PM count — three PM2_Supporting alone
    should NOT reach Likely Pathogenic (which requires >=3 PM)."""
    codes = [_ev("PM2_Supporting"), _ev("PM2_Supporting"), _ev("PM2_Supporting")]
    result = classify_variant(codes)
    # Three PP is not enough for LP under standard ACMG rules
    assert result.classification != "Likely Pathogenic"


# I-5: Risk Factor bypass
def test_apoe_risk_factor_bypass():
    """APOE should return 'Risk Factor' regardless of evidence."""
    codes = [_ev("PM1"), _ev("PP3")]
    result = classify_variant(codes, gene="APOE")
    assert result.classification == "Risk Factor"
    assert result.matched_rule == "risk_factor_bypass"


def test_apoe_risk_factor_no_evidence():
    codes = []
    result = classify_variant(codes, gene="APOE")
    assert result.classification == "Risk Factor"


def test_non_risk_factor_gene_not_bypassed():
    codes = [_ev("PM1")]
    result = classify_variant(codes, gene="BRCA1")
    assert result.classification != "Risk Factor"


# I-4: ClinVar conflict detection
def test_check_clinvar_conflict_lp_vs_pathogenic():
    """LP (rank 3) vs Pathogenic (rank 4) — diff=1 → NOT conflict (clinically acceptable)."""
    assert check_clinvar_conflict("Likely Pathogenic", "Pathogenic") is False


def test_check_clinvar_conflict_benign_vs_pathogenic():
    """Benign (rank 0) vs Pathogenic (rank 4) — diff=4 → conflict."""
    assert check_clinvar_conflict("Benign", "Pathogenic") is True


def test_check_clinvar_no_conflict_same():
    assert check_clinvar_conflict("Pathogenic", "Pathogenic") is False


def test_check_clinvar_no_conflict_not_found():
    assert check_clinvar_conflict("Likely Pathogenic", "Not Found") is False


def test_check_clinvar_no_conflict_drug_response():
    assert check_clinvar_conflict("Drug Response", "Drug Response") is False


def test_check_clinvar_no_conflict_empty():
    assert check_clinvar_conflict("Pathogenic", "") is False


# ─── _count_by_strength: strength modifier suffix handling ──────────────────
# Regression tests for Pejaver 2022 3-tier REVEL support. Plain codes count
# at their prefix strength; _Strong/_Moderate upgrade PP→ps/pm; _Supporting
# downgrades PVS/PS/PM→pp (ClinGen SVI convention).


def test_count_pp3_strong_upgrades_to_ps():
    """PP3_Strong (REVEL >= 0.932) counts as 1 ps, not 1 pp."""
    counts = _count_by_strength([_ev("PP3_Strong")])
    assert counts["ps"] == 1
    assert counts["pp"] == 0


def test_count_pp3_moderate_upgrades_to_pm():
    """PP3_Moderate (REVEL >= 0.773) counts as 1 pm, not 1 pp."""
    counts = _count_by_strength([_ev("PP3_Moderate")])
    assert counts["pm"] == 1
    assert counts["pp"] == 0


def test_count_plain_pp3_counts_as_pp():
    """Plain PP3 (REVEL [0.644, 0.773)) counts as 1 pp."""
    counts = _count_by_strength([_ev("PP3")])
    assert counts["pp"] == 1
    assert counts["ps"] == 0
    assert counts["pm"] == 0


def test_count_bp4_strong_upgrades_to_bs():
    """BP4_Strong counts as 1 bs."""
    counts = _count_by_strength([_ev("BP4_Strong")])
    assert counts["bs"] == 1
    assert counts["bp"] == 0


def test_count_bp4_moderate_counts_as_bs():
    """BP4_Moderate counts as 1 bs (benign upgrade path)."""
    counts = _count_by_strength([_ev("BP4_Moderate")])
    assert counts["bs"] == 1
    assert counts["bp"] == 0


def test_count_plain_bp4_counts_as_bp():
    """Plain BP4 counts as 1 bp (natural prefix strength)."""
    counts = _count_by_strength([_ev("BP4")])
    assert counts["bp"] == 1
    assert counts["bs"] == 0


def test_count_pvs1_strong_counts_as_ps():
    """PVS1_Strong (pLI 0.5-0.9 downgrade) counts as 1 ps, not 1 pp."""
    counts = _count_by_strength([_ev("PVS1_Strong")])
    assert counts["ps"] == 1
    assert counts["pvs"] == 0
    assert counts["pp"] == 0


def test_count_ps2_moderate_stays_pm():
    """PS2_Moderate (ClinGen SVI 2018 de novo downgrade) counts as 1 pm."""
    counts = _count_by_strength([_ev("PS2_Moderate")])
    assert counts["pm"] == 1
    assert counts["ps"] == 0


def test_count_pm2_supporting_stays_pp():
    """PM2_Supporting (existing downgrade) still counts as 1 pp (unchanged)."""
    counts = _count_by_strength([_ev("PM2_Supporting")])
    assert counts["pp"] == 1
    assert counts["pm"] == 0


def test_count_pm6_supporting_stays_pp():
    """PM6_Supporting (existing downgrade) still counts as 1 pp (unchanged)."""
    counts = _count_by_strength([_ev("PM6_Supporting")])
    assert counts["pp"] == 1
    assert counts["pm"] == 0


def test_count_pm1_supporting_stays_pp():
    """PM1_Supporting (hotspot supporting-tier) still counts as 1 pp."""
    counts = _count_by_strength([_ev("PM1_Supporting")])
    assert counts["pp"] == 1
    assert counts["pm"] == 0


def test_classify_pm2_plus_pp3_strong_reaches_lp():
    """PM2 + PP3_Strong → LP via ps:1 + pm:1 (under-call fix verification).

    Before the fix, PP3_Strong silently counted as pp, giving pm:1 + pp:1
    which is only VUS. With the fix, PP3_Strong counts as ps:1, triggering
    the 'ps:1 + pm:1' LP rule.
    """
    codes = [_ev("PM2"), _ev("PP3_Strong")]
    result = classify_variant(codes)
    assert result.classification == "Likely Pathogenic"


def test_classify_pm2_plus_pp3_moderate_reaches_lp():
    """PM2 + PP3_Moderate → LP via pm:2 (three-PM rule edge).

    PP3_Moderate counts as pm:1, combined with PM2 yields pm:2, which is
    NOT by itself LP (LP needs pm:3, ps:1+pm:1, ps:1+pm:2, ps:1+pp:2,
    pvs:1+pm:1, pm:2+pp:2, or pm:1+pp:4). pm:2 alone stays VUS.
    """
    codes = [_ev("PM2"), _ev("PP3_Moderate")]
    result = classify_variant(codes)
    # pm:2 alone does not reach LP (needs pm:3, or pm:2+pp:2 etc.)
    assert result.classification == "VUS"


def test_classify_plain_pp3_does_not_inflate_pm():
    """Plain PP3 must count as pp, not pm — regression for naive prefix match."""
    # PM1 + PP3 (plain) → pm:1 + pp:1 → VUS (not LP; LP needs pm:1+pp:4 or ps:1+pm:1 etc.)
    codes = [_ev("PM1"), _ev("PP3")]
    result = classify_variant(codes)
    assert result.classification == "VUS"


# ─── v2.5.3 ClinGen SVI 2018 PP5 / BP6 deprecation ──────────────────────────
# PP5 / BP6 double-count the same ClinVar evidence already captured by PS1 /
# PM5 / population / function / segregation codes. SVI 2018 therefore retired
# them from scoring. BIKO v2.5.3 drops PP5 / BP6 from _count_by_strength by
# default while keeping emission (query_clinvar / query_local_clinvar /
# parse_intervar) unchanged, so the codes still land in the report evidence
# list — only the counter is affected.


def test_count_by_strength_pp5_excluded_by_default():
    """PP5 must NOT contribute to counts['pp'] under default config (SVI 2018)."""
    counts = _count_by_strength([_ev("PS1"), _ev("PP5"), _ev("PM2_Supporting")])
    # PM2_Supporting downgrades to pp → counts['pp'] == 1 (PM2_Supporting only).
    # PP5 is silently skipped from scoring.
    assert counts["pp"] == 1
    assert counts["ps"] == 1
    assert counts["pm"] == 0


def test_count_by_strength_bp6_excluded_by_default():
    """BP6 must NOT contribute to counts['bp'] under default config (SVI 2018)."""
    counts = _count_by_strength([_ev("BA1"), _ev("BP6"), _ev("BP4")])
    # BP4 alone should be counted; BP6 skipped.
    assert counts["bp"] == 1
    assert counts["ba"] == 1


def test_count_by_strength_pp5_bp6_restored_via_config(monkeypatch):
    """acmg.use_pp5_bp6 = true restores pre-v2.5.3 scoring for back-compat."""
    from scripts.common import config as _cfg

    original_get = _cfg.get

    def _patched_get(key, default=None):
        if key == "acmg.use_pp5_bp6":
            return True
        return original_get(key, default)

    monkeypatch.setattr("scripts.classification.acmg_engine.get", _patched_get)

    counts = _count_by_strength([_ev("PS1"), _ev("PP5"), _ev("BP6")])
    # Legacy behaviour: PP5 counted as pp, BP6 as bp.
    assert counts["pp"] == 1
    assert counts["bp"] == 1
    assert counts["ps"] == 1


def test_classify_clinvar_pathogenic_no_double_count():
    """ClinVar Pathogenic emits both PS1 and PP5 (expert-panel path). Under
    SVI 2018 scoring only PS1 should drive classification — i.e. 1 PS alone
    does NOT reach LP (requires 1 PS + 1 PM or 2 PS).
    """
    codes = [_ev("PS1"), _ev("PP5")]  # as emitted by query_local_clinvar for expert-panel Pathogenic
    result = classify_variant(codes)
    # With PP5 ignored: counts = {ps: 1}. 1 PS alone is VUS, not LP.
    # If PP5 were still counted: {ps: 1, pp: 1} → still VUS under Richards 2015
    # (LP needs 1 PS + 1 PM or 1 PS + 2 PP), so this test mainly anchors the
    # scoring floor: PS1 + PP5 never reaches LP on its own.
    assert result.classification == "VUS"
    # Annotation still carries both codes — scoring and annotation are independent.
    assert "PS1" in result.evidence_codes
    assert "PP5" in result.evidence_codes


def test_classify_pp5_alone_is_vus():
    """A bare PP5 (no other codes) must not lift a variant above VUS."""
    codes = [_ev("PP5")]
    result = classify_variant(codes)
    assert result.classification == "VUS"
    # Annotation retained
    assert "PP5" in result.evidence_codes
