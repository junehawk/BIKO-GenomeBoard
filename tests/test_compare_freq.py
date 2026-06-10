from scripts.common.models import FrequencyData
from scripts.population.compare_freq import compare_frequencies


def test_korean_rare_pm2():
    freq = FrequencyData(kova=0.0001, gnomad_eas=0.0003, gnomad_all=0.0002)
    result = compare_frequencies(freq)
    assert "PM2_Supporting" in result["acmg_codes"]


def test_korean_common_bs1():
    freq = FrequencyData(kova=0.05, gnomad_eas=0.04, gnomad_all=0.03)
    result = compare_frequencies(freq)
    assert "BS1" in result["acmg_codes"]


def test_very_common_ba1():
    freq = FrequencyData(kova=0.06, gnomad_eas=0.05, gnomad_all=0.05)
    result = compare_frequencies(freq)
    assert "BA1" in result["acmg_codes"]


def test_absent_from_all_dbs_fires_pm2_supporting():
    """Absent from every queried frequency DB (gnomAD ALL/EAS + KOVA) is the canonical
    ACMG PM2 case ('absent from controls'), emitted at Supporting strength per ClinGen
    SVI 2020 — not 'no data / no code' (T1-03)."""
    freq = FrequencyData(kova=None, gnomad_eas=None, gnomad_all=None)
    result = compare_frequencies(freq)
    assert "PM2_Supporting" in result["acmg_codes"]
    assert "Absent" in result["korean_flag"]


def test_korean_specific_variant():
    freq = FrequencyData(kova=0.001, gnomad_eas=None, gnomad_all=None)
    result = compare_frequencies(freq)
    assert "Korean-specific variant (KOVA only)" in result["korean_flag"]


def test_korean_enrichment_5x():
    """KOVA AF 5× higher than gnomAD EAS → enrichment flag."""
    freq = FrequencyData(kova=0.005, gnomad_eas=0.0005, gnomad_all=0.001)
    result = compare_frequencies(freq)
    assert "Korean frequency 5x+ higher than East Asian" in result["korean_flag"]


def test_korean_depletion_vs_eas():
    """KOVA AF much lower than gnomAD EAS → depletion flag."""
    freq = FrequencyData(kova=0.0001, gnomad_eas=0.002, gnomad_all=0.001)
    result = compare_frequencies(freq)
    assert "Korean frequency much lower than East Asian" in result["korean_flag"]


# I-5b: PM2 fires at Supporting strength across both rare bands.
# ClinGen SVI 2020 recommends PM2 at Supporting by default; the prior code's
# full-strength PM2 for the moderately-rare band was an evidence-strength
# inversion (rarer variants got *weaker* PM2_Supporting) — v2.7 CLAS-06.
def test_pm2_supporting_freq_mid_range():
    """Freq 0.005 (between 0.001 and 0.01) → PM2_Supporting, not full PM2."""
    freq = FrequencyData(kova=0.005, gnomad_eas=None, gnomad_all=None)
    result = compare_frequencies(freq)
    assert "PM2_Supporting" in result["acmg_codes"]
    assert "PM2" not in result["acmg_codes"]


def test_pm2_supporting_at_low_end():
    """Freq 0.0011 just above the PM2_Supporting threshold → PM2_Supporting."""
    freq = FrequencyData(kova=0.0011, gnomad_eas=None, gnomad_all=None)
    result = compare_frequencies(freq)
    assert "PM2_Supporting" in result["acmg_codes"]
    assert "PM2" not in result["acmg_codes"]


def test_pm2_supporting_at_threshold():
    """Freq exactly 0.001 → PM2_Supporting, not PM2."""
    freq = FrequencyData(kova=0.001, gnomad_eas=None, gnomad_all=None)
    result = compare_frequencies(freq)
    assert "PM2_Supporting" in result["acmg_codes"]
    assert "PM2" not in result["acmg_codes"]


def test_kova_homozygote_preserved():
    """KOVA homozygote count should be preserved on the FrequencyData record."""
    freq = FrequencyData(kova=0.02, gnomad_eas=0.01, gnomad_all=0.01, kova_homozygote=3)
    result = compare_frequencies(freq)
    assert result["frequencies"].kova_homozygote == 3
