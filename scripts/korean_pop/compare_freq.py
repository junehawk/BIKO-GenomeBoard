from typing import Dict, List
from scripts.common.models import FrequencyData
from scripts.common.config import get

BA1_THRESHOLD = get("thresholds.ba1", 0.05)
BS1_THRESHOLD = get("thresholds.bs1", 0.01)
PM2_THRESHOLD = get("thresholds.pm2", 0.001)  # for PM2_Supporting


def compare_frequencies(freq: FrequencyData) -> Dict:
    acmg_codes: List[str] = []
    flags: List[str] = []

    max_freq = (
        max(f for f in [freq.krgdb, freq.gnomad_eas, freq.gnomad_all] if f is not None)
        if any(f is not None for f in [freq.krgdb, freq.gnomad_eas, freq.gnomad_all])
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
    if freq.krgdb is not None and freq.gnomad_all is not None and freq.gnomad_all > 0:
        ratio = freq.krgdb / freq.gnomad_all
        if ratio >= 5:
            flags.append("Korean frequency 5x+ higher than global")
        elif ratio <= 0.2:
            flags.append("Korean frequency much lower than global")
    elif freq.krgdb is not None and freq.gnomad_eas is None and freq.gnomad_all is None:
        flags.append("Korean-specific variant (KRGDB only)")

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
    result["frequencies"] = {"krgdb": freq.krgdb, "gnomad_eas": freq.gnomad_eas, "gnomad_all": freq.gnomad_all}
    print(json.dumps(result, indent=2, ensure_ascii=False))
