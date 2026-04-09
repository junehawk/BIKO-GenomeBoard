from typing import Dict, List
from scripts.common.models import FrequencyData
from scripts.common.config import get

BA1_THRESHOLD = get("thresholds.ba1", 0.05)
BS1_THRESHOLD = get("thresholds.bs1", 0.01)
PM2_THRESHOLD = get("thresholds.pm2", 0.001)  # for PM2_Supporting


def compare_frequencies(freq: FrequencyData) -> Dict:
    acmg_codes: List[str] = []
    flags: List[str] = []

    all_freqs = [freq.krgdb, freq.korea4k, freq.nard2, freq.gnomad_eas, freq.gnomad_all]
    max_freq = (
        max(f for f in all_freqs if f is not None)
        if any(f is not None for f in all_freqs)
        else None
    )

    if max_freq is None:
        return {
            "acmg_codes": [],
            "korean_flag": "No frequency data available",
            "frequencies": freq,
        }

    # BA1: stand-alone benign
    if max_freq > BA1_THRESHOLD:
        acmg_codes.append("BA1")
        flags.append("Very common variant")
    # BS1: strong benign
    elif max_freq >= BS1_THRESHOLD:
        acmg_codes.append("BS1")
        flags.append("Common variant")
    # PM2_Supporting: rare
    elif max_freq <= PM2_THRESHOLD:
        acmg_codes.append("PM2_Supporting")
        flags.append("Rare variant")
    # PM2: moderately rare (between PM2_Supporting and BS1 thresholds)
    elif max_freq < BS1_THRESHOLD:
        acmg_codes.append("PM2")
        flags.append("Low frequency variant")

    # Korean-specific flags
    korean_max = freq.korean_max
    if korean_max is not None and freq.gnomad_all is not None and freq.gnomad_all > 0:
        ratio = korean_max / freq.gnomad_all
        if ratio >= 5:
            flags.append("Korean frequency 5x+ higher than global")
        elif ratio <= 0.2:
            flags.append("Korean frequency much lower than global")
    elif korean_max is not None and freq.gnomad_eas is None and freq.gnomad_all is None:
        korean_sources = []
        if freq.krgdb is not None:
            korean_sources.append("KRGDB")
        if freq.korea4k is not None:
            korean_sources.append("Korea4K")
        if freq.nard2 is not None:
            korean_sources.append("NARD2")
        flags.append(f"Korean-specific variant ({'/'.join(korean_sources)} only)")

    return {
        "acmg_codes": acmg_codes,
        "korean_flag": " | ".join(flags) if flags else "No notable findings",
        "frequencies": freq,
    }


if __name__ == "__main__":
    import sys
    import json

    krgdb = float(sys.argv[1]) if len(sys.argv) > 1 else None
    eas = float(sys.argv[2]) if len(sys.argv) > 2 else None
    gnomad_all = float(sys.argv[3]) if len(sys.argv) > 3 else None
    freq = FrequencyData(krgdb=krgdb, gnomad_eas=eas, gnomad_all=gnomad_all)
    result = compare_frequencies(freq)
    result["frequencies"] = {"krgdb": freq.krgdb, "korea4k": freq.korea4k, "nard2": freq.nard2, "gnomad_eas": freq.gnomad_eas, "gnomad_all": freq.gnomad_all}
    print(json.dumps(result, indent=2, ensure_ascii=False))
