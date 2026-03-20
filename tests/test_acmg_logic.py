# tests/test_acmg_logic.py
from scripts.classification.acmg_engine import classify_variant, check_clinvar_conflict
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
    """LP (rank 3) vs Pathogenic (rank 4) — diff=1 → conflict."""
    assert check_clinvar_conflict("Likely Pathogenic", "Pathogenic") is True

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
