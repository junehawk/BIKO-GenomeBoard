# tests/test_acmg_logic.py
from scripts.classification.acmg_engine import classify_variant
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
