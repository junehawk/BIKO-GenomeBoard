"""3-tier allele-frequency comparator (KOVA + gnomAD EAS + gnomAD ALL).

KOVA v7 is the sole Korean-cohort allele-frequency source. The Korean-
enrichment ratio is defined as ``kova / gnomad_eas`` to compare Korean
allele frequencies against the broader East Asian reference.
"""

from typing import Dict, List, Optional

from scripts.common.config import get
from scripts.common.models import FrequencyData

BA1_THRESHOLD = get("thresholds.ba1", 0.05)
BS1_THRESHOLD = get("thresholds.bs1", 0.01)
PM2_THRESHOLD = get("thresholds.pm2", 0.001)  # for PM2_Supporting
# BS2: minimum homozygous observations in the healthy KOVA cohort to call a variant
# benign for a recessive condition (ACMG BS2). Configurable; conservative default 2.
BS2_HOMOZYGOTE_COUNT = get("thresholds.bs2_homozygote_count", 2)


def compare_frequencies(freq: FrequencyData, inheritance: Optional[str] = None) -> Dict:
    """Compare a variant's 3-tier frequency record and derive ACMG/korean flags.

    ``inheritance`` (the gene's OMIM mode, e.g. "AR" / "AD/AR") gates BS2: a variant
    observed homozygous in healthy KOVA controls for a recessive condition is benign
    (ACMG BS2). Omitted → BS2 never fires (T2-07).
    """
    acmg_codes: List[str] = []
    flags: List[str] = []

    all_freqs = [freq.kova, freq.gnomad_eas, freq.gnomad_all]
    max_freq = max(f for f in all_freqs if f is not None) if any(f is not None for f in all_freqs) else None

    if max_freq is None:
        # Absent from every queried frequency DB (gnomAD ALL/EAS + KOVA) → ACMG PM2
        # ("absent from controls"), emitted at Supporting strength per ClinGen SVI 2020.
        # NB: assumes the lookups ran and genuinely found nothing. If a frequency DB is
        # unavailable the availability cache logs a WARNING upstream; a degraded run can
        # therefore over-apply PM2_Supporting — surfaced provenance, not a silent error.
        return {
            "acmg_codes": ["PM2_Supporting"],
            "korean_flag": "Absent from population databases",
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
    # PM2_Supporting: rare. ClinGen SVI 2020 recommends PM2 at Supporting
    # strength by default. Both rare bands therefore fire PM2_Supporting — the
    # prior code emitted full-strength PM2 for the *moderately*-rare band while
    # the *rarer* band got PM2_Supporting, an evidence-strength inversion that
    # over-called borderline-frequency variants (v2.7 review CLAS-06).
    elif max_freq <= PM2_THRESHOLD:
        acmg_codes.append("PM2_Supporting")
        flags.append("Rare variant")
    # Moderately rare (between PM2_Supporting and BS1 thresholds) — still only
    # Supporting strength, but flagged distinctly for the reader.
    elif max_freq < BS1_THRESHOLD:
        acmg_codes.append("PM2_Supporting")
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

    # BS2 — observed homozygous in healthy KOVA controls for a recessive condition.
    # Only for AR/recessive genes (a homozygous benign observation does not inform a
    # dominant disorder), and only above the configured homozygote count (T2-07).
    if (
        inheritance
        and "AR" in inheritance.upper()
        and freq.kova_homozygote is not None
        and freq.kova_homozygote >= BS2_HOMOZYGOTE_COUNT
    ):
        acmg_codes.append("BS2")
        flags.append(f"Homozygous in {freq.kova_homozygote} healthy KOVA control(s) — BS2 (recessive)")

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
