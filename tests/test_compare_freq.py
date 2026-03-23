from scripts.korean_pop.compare_freq import compare_frequencies
from scripts.common.models import FrequencyData


def test_korean_rare_pm2():
    freq = FrequencyData(krgdb=0.0001, gnomad_eas=0.0003, gnomad_all=0.0002)
    result = compare_frequencies(freq)
    assert "PM2_Supporting" in result["acmg_codes"]


def test_korean_common_bs1():
    freq = FrequencyData(krgdb=0.05, gnomad_eas=0.04, gnomad_all=0.03)
    result = compare_frequencies(freq)
    assert "BS1" in result["acmg_codes"]


def test_very_common_ba1():
    freq = FrequencyData(krgdb=0.06, gnomad_eas=0.05, gnomad_all=0.05)
    result = compare_frequencies(freq)
    assert "BA1" in result["acmg_codes"]


def test_no_data():
    freq = FrequencyData(krgdb=None, gnomad_eas=None, gnomad_all=None)
    result = compare_frequencies(freq)
    assert result["acmg_codes"] == []
    assert "No frequency data available" in result["korean_flag"]


def test_korean_specific_variant():
    freq = FrequencyData(krgdb=0.001, gnomad_eas=None, gnomad_all=None)
    result = compare_frequencies(freq)
    assert "Korean-specific variant" in result["korean_flag"]


# I-5b: PM2 moderate for freq between PM2_THRESHOLD and BS1_THRESHOLD
def test_pm2_moderate_freq_mid_range():
    """Freq 0.005 is between 0.001 and 0.01 → PM2 (moderate)."""
    freq = FrequencyData(krgdb=0.005, gnomad_eas=None, gnomad_all=None)
    result = compare_frequencies(freq)
    assert "PM2" in result["acmg_codes"]
    assert "PM2_Supporting" not in result["acmg_codes"]


def test_pm2_moderate_at_low_end():
    """Freq 0.0011 just above PM2_Supporting threshold → PM2."""
    freq = FrequencyData(krgdb=0.0011, gnomad_eas=None, gnomad_all=None)
    result = compare_frequencies(freq)
    assert "PM2" in result["acmg_codes"]


def test_pm2_supporting_at_threshold():
    """Freq exactly 0.001 → PM2_Supporting, not PM2."""
    freq = FrequencyData(krgdb=0.001, gnomad_eas=None, gnomad_all=None)
    result = compare_frequencies(freq)
    assert "PM2_Supporting" in result["acmg_codes"]
    assert "PM2" not in result["acmg_codes"]
