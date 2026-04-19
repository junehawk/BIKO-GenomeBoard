"""Tumor Mutational Burden (TMB) calculation.

TMB = nonsynonymous coding variant count / coding region size (Mb)

Reference: Chalmers et al. Genome Med. 2017;9(1):34.
FDA-approved FoundationOne CDx methodology.
"""

import logging
from dataclasses import dataclass, field
from typing import List, Optional

from scripts.common.config import get

logger = logging.getLogger(__name__)

DEFAULT_PANEL_SIZE_MB = 33.0
DEFAULT_HIGH_THRESHOLD = 10.0
DEFAULT_INTERMEDIATE_THRESHOLD = 6.0
DEFAULT_COUNTED_CONSEQUENCES = [
    "missense_variant",
    "stop_gained",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "inframe_deletion",
    "inframe_insertion",
]


@dataclass
class TmbResult:
    score: float
    level: str
    variant_count: int
    total_variants: int
    panel_size_mb: float
    counted_consequences: List[str] = field(default_factory=list)


def calculate_tmb(
    variants,
    panel_size_mb: Optional[float] = None,
    high_threshold: Optional[float] = None,
    intermediate_threshold: Optional[float] = None,
    counted_consequences: Optional[List[str]] = None,
) -> TmbResult:
    """Calculate TMB from a list of variant objects.

    Args:
        variants: List of objects with a .consequence attribute (str)
        panel_size_mb: Coding region size in Mb (default: 33.0 for WGS exome)
        high_threshold: TMB-High cutoff (default: 10.0 mut/Mb)
        intermediate_threshold: TMB-Intermediate cutoff (default: 6.0 mut/Mb)
        counted_consequences: Which consequence types to count
    """
    panel_size = panel_size_mb or float(
        get("somatic.tmb_default_panel_size_mb", DEFAULT_PANEL_SIZE_MB) or DEFAULT_PANEL_SIZE_MB
    )
    high_thresh = high_threshold or float(
        get("somatic.tmb_high_threshold", DEFAULT_HIGH_THRESHOLD) or DEFAULT_HIGH_THRESHOLD
    )
    int_thresh = intermediate_threshold or float(
        get("somatic.tmb_intermediate_threshold", DEFAULT_INTERMEDIATE_THRESHOLD) or DEFAULT_INTERMEDIATE_THRESHOLD
    )
    consequences = counted_consequences or get("somatic.tmb_counted_consequences") or DEFAULT_COUNTED_CONSEQUENCES

    # Build normalized lookup set (lowercase) for matching both
    # VEP raw (missense_variant) and formatted (Missense) consequence names
    _FORMATTED_TO_RAW = {
        "missense": "missense_variant",
        "nonsense": "stop_gained",
        "nonsense / stop gain": "stop_gained",
        "stop gain": "stop_gained",
        "frameshift": "frameshift_variant",
        "splice donor": "splice_donor_variant",
        "splice acceptor": "splice_acceptor_variant",
        "inframe deletion": "inframe_deletion",
        "inframe insertion": "inframe_insertion",
        "splice region variant": "splice_region_variant",
    }
    csq_set = set(consequences)
    csq_lower_set = {c.lower() for c in csq_set}

    total = len(variants)
    counted = 0
    for v in variants:
        csq = getattr(v, "consequence", None) or ""
        # Direct match (VEP raw format)
        if csq in csq_set:
            counted += 1
            continue
        # Normalized match (formatted consequence)
        csq_lower = csq.lower()
        if csq_lower in csq_lower_set:
            counted += 1
            continue
        # Check formatted → raw mapping (e.g., "Missense" → "missense_variant")
        # Also handle compound consequences like "Missense / Splice Region Variant"
        for part in csq_lower.split(" / "):
            part = part.strip()
            raw = _FORMATTED_TO_RAW.get(part, "")
            if raw in csq_set:
                counted += 1
                break

    score = round(counted / panel_size, 1) if panel_size > 0 else 0.0

    if score >= high_thresh:
        level = "High"
    elif score >= int_thresh:
        level = "Intermediate"
    else:
        level = "Low"

    return TmbResult(
        score=score,
        level=level,
        variant_count=counted,
        total_variants=total,
        panel_size_mb=panel_size,
        counted_consequences=list(consequences),
    )


def calculate_panel_size_from_bed(bed_path: str) -> float:
    """Calculate total panel size in Mb from a BED file."""
    total_bp = 0
    try:
        with open(bed_path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    start = int(parts[1])
                    end = int(parts[2])
                    total_bp += end - start
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Error reading BED file {bed_path}: {e}")
        return DEFAULT_PANEL_SIZE_MB

    return total_bp / 1_000_000
