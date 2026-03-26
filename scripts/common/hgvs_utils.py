"""Shared HGVSp conversion utilities."""
import re
from typing import Optional

AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}


def hgvsp_to_civic_variant(hgvsp: Optional[str]) -> Optional[str]:
    """Convert HGVSp to CIViC variant name format. p.Gly12Asp -> G12D"""
    if not hgvsp:
        return None
    m = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvsp)
    if m:
        aa1 = AA3TO1.get(m.group(1), "?")
        pos = m.group(2)
        aa2 = AA3TO1.get(m.group(3), "?")
        return f"{aa1}{pos}{aa2}"
    return None
