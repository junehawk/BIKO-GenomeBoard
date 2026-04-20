# scripts/annotation/parse_intervar.py
"""Parse InterVar TSV output and extract ACMG evidence codes per variant.

InterVar produces a tab-separated file with per-variant rows.  The first
columns describe the variant (Chr, Start, End, Ref, Alt, ...) and later
columns hold 0/1 flags for each ACMG evidence code (PVS1 through BP7).

This module provides:
  - parse_intervar(path)          → dict mapping variant_key to evidence list
  - get_intervar_evidence(v, d)   → evidence list for a single Variant object
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List

from scripts.common.models import Variant

logger = logging.getLogger(__name__)

# Ordered ACMG evidence codes expected in InterVar output columns.
EVIDENCE_CODES = [
    "PVS1",
    "PS1",
    "PS2",
    "PS3",
    "PS4",
    "PM1",
    "PM2",
    "PM3",
    "PM4",
    "PM5",
    "PM6",
    "PP1",
    "PP2",
    "PP3",
    "PP4",
    "PP5",
    "BA1",
    "BS1",
    "BS2",
    "BS3",
    "BS4",
    "BP1",
    "BP2",
    "BP3",
    "BP4",
    "BP5",
    "BP6",
    "BP7",
]


def _strip_bom(text: str) -> str:
    """Remove UTF-8 BOM if present."""
    return text.lstrip("\ufeff")


def _normalise_chrom(chrom: str) -> str:
    """Ensure chromosome has 'chr' prefix for consistency."""
    chrom = chrom.strip()
    if not chrom.lower().startswith("chr"):
        return f"chr{chrom}"
    return chrom


def _build_variant_key(chrom: str, pos: str, ref: str, alt: str) -> str:
    """Build BIKO GenomeBoard-style variant key: chr17:7675088:C>A"""
    return f"{_normalise_chrom(chrom)}:{pos.strip()}:{ref.strip()}>{alt.strip()}"


def parse_intervar(path: str) -> Dict[str, List[str]]:
    """Parse an InterVar TSV file and return triggered evidence codes per variant.

    Args:
        path: Path to InterVar output TSV file.

    Returns:
        Dict mapping variant_key (``chr:pos:ref>alt``) to a list of evidence
        code strings that were triggered (value == "1").

    Raises:
        FileNotFoundError: If *path* does not exist.
    """
    filepath = Path(path)
    if not filepath.exists():
        raise FileNotFoundError(f"InterVar file not found: {path}")

    result: Dict[str, List[str]] = {}

    with open(filepath, encoding="utf-8-sig") as fh:
        # Peek at the first line to detect the header.
        first_line = fh.readline()
        if not first_line.strip():
            return result  # empty file

        first_line = _strip_bom(first_line)

        # Build column name → index map from the header.
        headers = first_line.rstrip("\n\r").split("\t")
        col_idx = {h.strip(): i for i, h in enumerate(headers)}

        # Identify positional columns.  InterVar uses "#Chr" or "Chr".
        chr_col = col_idx.get("#Chr", col_idx.get("Chr"))
        start_col = col_idx.get("Start")
        ref_col = col_idx.get("Ref")
        alt_col = col_idx.get("Alt")

        if chr_col is None or start_col is None or ref_col is None or alt_col is None:
            logger.warning("InterVar header missing required columns (#Chr/Start/Ref/Alt)")
            return result

        # Map evidence code names to column indices.
        code_indices: List[tuple] = []
        for code in EVIDENCE_CODES:
            idx = col_idx.get(code)
            if idx is not None:
                code_indices.append((code, idx))

        if not code_indices:
            logger.warning("InterVar file has no evidence code columns")
            return result

        # Parse data rows.
        for line_no, line in enumerate(fh, start=2):
            line = line.rstrip("\n\r")
            if not line or line.startswith("##"):
                continue

            fields = line.split("\t")
            # Guard against truncated lines.
            max_needed = max(chr_col, start_col, ref_col, alt_col, max(idx for _, idx in code_indices)) + 1
            if len(fields) < max_needed:
                logger.debug("Skipping truncated line %d (%d fields)", line_no, len(fields))
                continue

            try:
                chrom = fields[chr_col]
                pos = fields[start_col]
                ref = fields[ref_col]
                alt = fields[alt_col]
            except IndexError:
                continue

            if not chrom or not pos or not ref or not alt:
                continue

            key = _build_variant_key(chrom, pos, ref, alt)

            triggered: List[str] = []
            for code, idx in code_indices:
                val = fields[idx].strip()
                if val == "1":
                    triggered.append(code)

            # Only store variants that have at least one triggered code
            # (but always create entry so caller can distinguish "all zero" from "missing").
            result[key] = triggered

    logger.info("Parsed %d variants from InterVar file %s", len(result), path)
    return result


def get_intervar_evidence(variant: Variant, intervar_data: Dict[str, List[str]]) -> List[str]:
    """Look up InterVar evidence codes for a Variant object.

    Args:
        variant: A Variant with chrom, pos, ref, alt fields.
        intervar_data: Dict returned by :func:`parse_intervar`.

    Returns:
        List of triggered evidence code strings, or empty list if the variant
        is not present in *intervar_data*.
    """
    if not intervar_data or variant is None:
        return []

    key = f"{variant.chrom}:{variant.pos}:{variant.ref}>{variant.alt}"
    return intervar_data.get(key, [])
