"""3-tier allele-frequency comparator (KOVA + gnomAD EAS + gnomAD ALL).

KOVA v7 replaced the legacy KRGDB / Korea4K / NARD2 TSV sources in v2.4.
The Korean-enrichment ratio is defined as ``kova / gnomad_eas`` to compare
Korean allele frequencies against the broader East Asian reference.
"""

from typing import Dict, List

from scripts.common.config import get
from scripts.common.models import FrequencyData

BA1_THRESHOLD = get("thresholds.ba1", 0.05)
BS1_THRESHOLD = get("thresholds.bs1", 0.01)
PM2_THRESHOLD = get("thresholds.pm2", 0.001)  # for PM2_Supporting


def compare_frequencies(freq: FrequencyData) -> Dict:
    """Compare a variant's 3-tier frequency record and derive ACMG/korean flags."""
    acmg_codes: List[str] = []
    flags: List[str] = []

    all_freqs = [freq.kova, freq.gnomad_eas, freq.gnomad_all]
    max_freq = max(f for f in all_freqs if f is not None) if any(f is not None for f in all_freqs) else None

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

    # Korean-specific flags — KOVA vs gnomAD EAS enrichment
    kova_af = freq.kova
    if kova_af is not None and freq.gnomad_eas is not None and freq.gnomad_eas > 0:
        ratio = kova_af / freq.gnomad_eas
        if ratio >= 5:
            flags.append("Korean frequency 5x+ higher than East Asian")
        elif ratio <= 0.2:
            flags.append("Korean frequency much lower than East Asian")
    elif kova_af is not None and freq.gnomad_eas is None and freq.gnomad_all is None:
        flags.append("Korean-specific variant (KOVA only)")

    return {
        "acmg_codes": acmg_codes,
        "korean_flag": " | ".join(flags) if flags else "No notable findings",
        "frequencies": freq,
    }


if __name__ == "__main__":
    import json
    import sys

    kova_af = float(sys.argv[1]) if len(sys.argv) > 1 else None
    eas = float(sys.argv[2]) if len(sys.argv) > 2 else None
    gnomad_all = float(sys.argv[3]) if len(sys.argv) > 3 else None
    freq = FrequencyData(kova=kova_af, gnomad_eas=eas, gnomad_all=gnomad_all)
    result = compare_frequencies(freq)
    result["frequencies"] = {
        "kova": freq.kova,
        "kova_homozygote": freq.kova_homozygote,
        "gnomad_eas": freq.gnomad_eas,
        "gnomad_all": freq.gnomad_all,
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))
