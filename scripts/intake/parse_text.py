# scripts/intake/parse_text.py
import re
from typing import List
from scripts.common.models import Variant


def parse_text(text: str) -> List[Variant]:
    """Parse comma- or newline-separated variant strings.
    Supports: chr17:7577120 G>A format.
    """
    entries = re.split(r"[,\n]+", text.strip())
    variants = []
    for entry in entries:
        entry = entry.strip()
        if not entry:
            continue
        variants.append(Variant.from_string(entry))
    return variants
