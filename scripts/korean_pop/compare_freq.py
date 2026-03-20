from typing import Dict, List
from scripts.common.models import FrequencyData

BA1_THRESHOLD = 0.05
BS1_THRESHOLD = 0.01
PM2_THRESHOLD = 0.001  # for PM2_Supporting

def compare_frequencies(freq: FrequencyData) -> Dict:
    acmg_codes: List[str] = []
    flags: List[str] = []

    max_freq = max(f for f in [freq.krgdb, freq.gnomad_eas, freq.gnomad_all] if f is not None) \
        if any(f is not None for f in [freq.krgdb, freq.gnomad_eas, freq.gnomad_all]) else None

    if max_freq is None:
        return {
            "acmg_codes": [],
            "korean_flag": "빈도 데이터 없음",
            "frequencies": freq,
        }

    # BA1: stand-alone benign
    if max_freq > BA1_THRESHOLD:
        acmg_codes.append("BA1")
        flags.append("매우 흔한 변이")
    # BS1: strong benign
    elif max_freq >= BS1_THRESHOLD:
        acmg_codes.append("BS1")
        flags.append("흔한 변이")
    # PM2_Supporting: rare
    elif max_freq <= PM2_THRESHOLD:
        acmg_codes.append("PM2_Supporting")
        flags.append("희귀 변이")

    # Korean-specific flags
    if freq.krgdb is not None and freq.gnomad_all is not None and freq.gnomad_all > 0:
        ratio = freq.krgdb / freq.gnomad_all
        if ratio >= 5:
            flags.append("한국인 빈도 글로벌 대비 5배 이상 높음")
        elif ratio <= 0.2:
            flags.append("한국인 빈도 글로벌 대비 매우 낮음")
    elif freq.krgdb is not None and freq.gnomad_eas is None and freq.gnomad_all is None:
        flags.append("한국인 특이 변이 (KRGDB only)")

    return {
        "acmg_codes": acmg_codes,
        "korean_flag": " | ".join(flags) if flags else "특이사항 없음",
        "frequencies": freq,
    }
