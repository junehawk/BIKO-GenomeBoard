"""Shared HGVSp conversion utilities.

Extended in v2.2 (A1 curated treatments) with the CIViC normaliser and a
generic protein-position extractor that tolerates 1-letter / 3-letter forms,
start-lost (``p.Met1?``), frameshift (``p.Arg249fs``), empty, and whitespace-
wrapped strings.
"""

from __future__ import annotations

import re
from typing import Optional

AA3TO1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Glu": "E",
    "Gln": "Q",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
}


def _strip(hgvsp: Optional[str]) -> str:
    if not hgvsp:
        return ""
    return str(hgvsp).strip()


def hgvsp_to_civic_variant(hgvsp: Optional[str]) -> Optional[str]:
    """Convert HGVSp to CIViC variant name format (``p.Gly12Asp`` → ``G12D``).

    Legacy entry point — preserved for ``scripts/storage/query_civic.py`` callers.
    Only recognises clean single-residue substitutions. For broader tolerance
    use :func:`normalize_hgvsp_for_civic`.
    """
    stripped = _strip(hgvsp)
    if not stripped:
        return None
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", stripped)
    if m:
        aa1 = AA3TO1.get(m.group(1), "?")
        pos = m.group(2)
        aa2 = AA3TO1.get(m.group(3), "?")
        return f"{aa1}{pos}{aa2}"
    return None


def normalize_hgvsp_for_civic(hgvsp: Optional[str]) -> Optional[str]:
    """Normalise an HGVSp string to a CIViC variant name.

    Handles both 3-letter (``p.Gly12Asp``) and 1-letter (``p.G12D``) forms
    plus whitespace wrapping. Returns ``None`` for unparseable input, empty
    strings, and non-substitution notations (``p.Met1?``, ``p.Arg249fs``)
    where there is no clean CIViC variant-name mapping — the curator falls
    through to gene-level in those cases.
    """
    stripped = _strip(hgvsp)
    if not stripped:
        return None

    # Three-letter AA substitution: p.Gly12Asp
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", stripped)
    if m:
        ref = AA3TO1.get(m.group(1))
        alt = AA3TO1.get(m.group(3))
        if ref and alt:
            return f"{ref}{m.group(2)}{alt}"
        return None

    # One-letter AA substitution: p.G12D
    m = re.match(r"p\.([A-Z])(\d+)([A-Z*])$", stripped)
    if m:
        return f"{m.group(1)}{m.group(2)}{m.group(3)}"

    return None


def extract_protein_position(hgvsp: Optional[str]) -> Optional[int]:
    """Extract the numeric residue position from an HGVSp string.

    Works across 1-letter, 3-letter, start-lost (``p.Met1?``), frameshift
    (``p.Arg249fs``), and whitespace-wrapped forms. Returns ``None`` when no
    residue number can be recovered.
    """
    stripped = _strip(hgvsp)
    if not stripped:
        return None

    # 3-letter prefix: p.Arg249Met / p.Met1? / p.Arg249fs
    m = re.search(r"p\.(?:[A-Z][a-z]{2})(\d+)", stripped)
    if m:
        return int(m.group(1))

    # 1-letter prefix: p.R249M / p.M1? / p.R249fs
    m = re.search(r"p\.[A-Z](\d+)", stripped)
    if m:
        return int(m.group(1))

    return None
