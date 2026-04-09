"""In silico prediction scoring and PP3/BP4 evidence generation.

Based on ClinGen SVI 2022 recommendations for computational evidence.

References:
  - Pejaver et al. (2022) "Calibration of computational tools for missense
    variant pathogenicity classification and ClinGen recommendations for
    PP3/BP4 criteria" AJHG 109(12):2163-2177
  - Tavtigian et al. (2023) SpliceAI thresholds for PP3/BP4
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional


# ---------------------------------------------------------------------------
# Default ClinGen SVI 2022 thresholds
# ---------------------------------------------------------------------------

DEFAULT_THRESHOLDS: Dict[str, Any] = {
    # REVEL (missense primary predictor)
    "revel_pp3_strong": 0.932,
    "revel_pp3_moderate": 0.644,
    "revel_bp4_strong": 0.016,
    "revel_bp4_moderate": 0.183,
    # SpliceAI (splice primary predictor)
    "spliceai_pp3_strong": 0.5,
    "spliceai_pp3_moderate": 0.2,
    "spliceai_bp4": 0.1,
    # Fallback: CADD + AlphaMissense (when REVEL unavailable)
    "cadd_pathogenic": 25.0,
}

# CSQ field name mappings — the left side is the lowercased CSQ field name,
# the right side is the InSilicoScores attribute it maps to.
_CSQ_FIELD_MAP: Dict[str, str] = {
    "revel_score": "revel",
    "revel": "revel",
    "cadd_phred": "cadd_phred",
    "cadd_raw_rankscore": None,          # ignored
    "am_class": "alphamissense_class",
    "am_pathogenicity": "alphamissense_score",
    "spliceai_pred_ds_ag": "spliceai_ds_ag",
    "spliceai_pred_ds_al": "spliceai_ds_al",
    "spliceai_pred_ds_dg": "spliceai_ds_dg",
    "spliceai_pred_ds_dl": "spliceai_ds_dl",
    "sift": "sift",
    "polyphen": "polyphen",
}

# VEP missing-value sentinels
_MISSING_VALUES = {"", ".", "-", "NA", "None"}


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------

@dataclass
class InSilicoScores:
    """In silico prediction scores parsed from VEP CSQ."""

    revel: Optional[float] = None              # 0-1, higher = more pathogenic
    cadd_phred: Optional[float] = None         # phred-scaled, >20 = top 1%
    alphamissense_score: Optional[float] = None  # 0-1
    alphamissense_class: Optional[str] = None  # likely_pathogenic / ambiguous / likely_benign
    spliceai_max: Optional[float] = None       # max of 4 delta scores
    spliceai_ds_ag: Optional[float] = None     # acceptor gain
    spliceai_ds_al: Optional[float] = None     # acceptor loss
    spliceai_ds_dg: Optional[float] = None     # donor gain
    spliceai_ds_dl: Optional[float] = None     # donor loss
    sift: Optional[str] = None                 # e.g. deleterious(0.01)
    polyphen: Optional[str] = None             # e.g. probably_damaging(0.998)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _safe_float(value: str) -> Optional[float]:
    """Convert a string to float, returning None for missing/invalid values."""
    if value is None or str(value).strip() in _MISSING_VALUES:
        return None
    try:
        f = float(value)
        # Reject NaN / Inf
        if f != f:  # NaN check
            return None
        return f
    except (ValueError, TypeError):
        return None


def _safe_str(value: str) -> Optional[str]:
    """Return a cleaned string, or None for VEP missing sentinels."""
    if value is None:
        return None
    s = str(value).strip()
    if s in _MISSING_VALUES:
        return None
    return s


def parse_in_silico_from_csq(csq_fields: Dict[str, str]) -> InSilicoScores:
    """Extract in silico scores from a parsed CSQ field dict.

    Parameters
    ----------
    csq_fields : dict
        Keys are **lowercased** CSQ field names (as produced by
        ``parse_annotation.parse_csq_value`` internals), values are raw
        strings from the VCF.

    Returns
    -------
    InSilicoScores
        Populated scores; missing values are None.
    """
    scores = InSilicoScores()

    for csq_key, attr_name in _CSQ_FIELD_MAP.items():
        if attr_name is None:
            continue
        raw = csq_fields.get(csq_key)
        if raw is None:
            continue

        # String fields
        if attr_name in ("sift", "polyphen", "alphamissense_class"):
            setattr(scores, attr_name, _safe_str(raw))
        else:
            setattr(scores, attr_name, _safe_float(raw))

    # Compute spliceai_max from the four delta scores
    scores.spliceai_max = _compute_spliceai_max(scores)

    return scores


def _compute_spliceai_max(scores: InSilicoScores) -> Optional[float]:
    """Return max of the four SpliceAI delta scores, or None if all absent."""
    vals = [
        v for v in (
            scores.spliceai_ds_ag,
            scores.spliceai_ds_al,
            scores.spliceai_ds_dg,
            scores.spliceai_ds_dl,
        )
        if v is not None
    ]
    return max(vals) if vals else None


# ---------------------------------------------------------------------------
# PP3 / BP4 evidence generation
# ---------------------------------------------------------------------------

def generate_pp3_bp4(
    scores: InSilicoScores,
    thresholds: Optional[Dict[str, Any]] = None,
) -> List[str]:
    """Generate PP3 or BP4 evidence codes based on in silico scores.

    Implements ClinGen SVI 2022 recommendations:
    - Missense: REVEL is primary; fall back to CADD + AlphaMissense
    - Splice: SpliceAI is primary
    - PP3 and BP4 are mutually exclusive

    Parameters
    ----------
    scores : InSilicoScores
        Parsed prediction scores.
    thresholds : dict, optional
        Override default thresholds. Keys mirror ``DEFAULT_THRESHOLDS``.

    Returns
    -------
    list[str]
        Zero or one evidence code, e.g. ``["PP3_Strong"]``, ``["BP4_Moderate"]``,
        or ``[]``.
    """
    if thresholds is None:
        t = DEFAULT_THRESHOLDS.copy()
    else:
        t = {**DEFAULT_THRESHOLDS, **thresholds}

    # --- Splice evidence (SpliceAI) ---
    splice_evidence = _evaluate_spliceai(scores, t)

    # --- Missense evidence (REVEL primary, CADD+AM fallback) ---
    missense_evidence = _evaluate_missense(scores, t)

    # Combine: if both present, pick the stronger one (higher rank).
    # PP3 and BP4 are mutually exclusive — never return both.
    candidates = []
    if splice_evidence:
        candidates.append(splice_evidence)
    if missense_evidence:
        candidates.append(missense_evidence)

    if not candidates:
        return []

    # If we have both splice and missense, pick the strongest.
    # If they conflict in direction (one PP3, one BP4), the PP3 wins per
    # ClinGen guidance (pathogenic evidence takes precedence when both
    # are applicable).
    pp3_candidates = [c for c in candidates if c.startswith("PP3")]
    bp4_candidates = [c for c in candidates if c.startswith("BP4")]

    if pp3_candidates and bp4_candidates:
        # Conflict — go with PP3 (pathogenic evidence)
        return [_strongest(pp3_candidates)]
    if pp3_candidates:
        return [_strongest(pp3_candidates)]
    if bp4_candidates:
        return [_strongest(bp4_candidates)]

    return []


_STRENGTH_RANK = {
    "Strong": 3,
    "Moderate": 2,
    "Supporting": 1,
}


def _strength_of(code: str) -> int:
    """Return numeric rank of an evidence code's strength suffix."""
    for suffix, rank in _STRENGTH_RANK.items():
        if code.endswith(suffix):
            return rank
    return 0


def _strongest(codes: List[str]) -> str:
    """Return the code with the highest strength from a list."""
    return max(codes, key=_strength_of)


def _evaluate_spliceai(
    scores: InSilicoScores,
    t: Dict[str, Any],
) -> Optional[str]:
    """Evaluate SpliceAI score against thresholds."""
    if scores.spliceai_max is None:
        return None

    val = scores.spliceai_max

    if val >= t["spliceai_pp3_strong"]:
        return "PP3_Strong"
    if val >= t["spliceai_pp3_moderate"]:
        return "PP3_Moderate"
    if val < t["spliceai_bp4"]:
        return "BP4_Supporting"

    # Between moderate and BP4 threshold — indeterminate
    return None


def _evaluate_missense(
    scores: InSilicoScores,
    t: Dict[str, Any],
) -> Optional[str]:
    """Evaluate missense predictors: REVEL primary, CADD+AM fallback."""
    # --- REVEL (primary) ---
    if scores.revel is not None:
        return _revel_evidence(scores.revel, t)

    # --- Fallback: CADD >= threshold AND AlphaMissense pathogenic ---
    return _fallback_cadd_am(scores, t)


def _revel_evidence(revel: float, t: Dict[str, Any]) -> Optional[str]:
    """Map a REVEL score to an evidence code."""
    if revel >= t["revel_pp3_strong"]:
        return "PP3_Strong"
    if revel >= t["revel_pp3_moderate"]:
        return "PP3_Moderate"
    if revel <= t["revel_bp4_strong"]:
        return "BP4_Strong"
    if revel <= t["revel_bp4_moderate"]:
        return "BP4_Moderate"
    # In the indeterminate zone — no evidence
    return None


def _fallback_cadd_am(
    scores: InSilicoScores,
    t: Dict[str, Any],
) -> Optional[str]:
    """Fallback when REVEL is unavailable: CADD + AlphaMissense."""
    cadd_ok = (
        scores.cadd_phred is not None
        and scores.cadd_phred >= t["cadd_pathogenic"]
    )
    am_pathogenic = (
        scores.alphamissense_class is not None
        and scores.alphamissense_class.lower() in ("likely_pathogenic", "pathogenic")
    )

    if cadd_ok and am_pathogenic:
        return "PP3_Supporting"

    # Check for benign fallback: low CADD + benign AlphaMissense
    cadd_benign = (
        scores.cadd_phred is not None
        and scores.cadd_phred < 15.0
    )
    am_benign = (
        scores.alphamissense_class is not None
        and scores.alphamissense_class.lower() in ("likely_benign", "benign")
    )

    if cadd_benign and am_benign:
        return "BP4_Supporting"

    return None


# ---------------------------------------------------------------------------
# Display formatting
# ---------------------------------------------------------------------------

def format_scores_for_display(scores: InSilicoScores) -> Dict[str, str]:
    """Format in silico scores for report display.

    Returns a dict with human-readable keys and formatted values.
    Only includes scores that are not None.

    Example::

        {
            "REVEL": "0.890",
            "CADD": "32.0",
            "AlphaMissense": "likely_pathogenic (0.950)",
            "SpliceAI": "0.85 (AL)",
            "SIFT": "deleterious(0.01)",
            "PolyPhen": "probably_damaging(0.998)",
        }
    """
    result: Dict[str, str] = {}

    if scores.revel is not None:
        result["REVEL"] = f"{scores.revel:.3f}"

    if scores.cadd_phred is not None:
        result["CADD"] = f"{scores.cadd_phred:.1f}"

    if scores.alphamissense_class is not None or scores.alphamissense_score is not None:
        parts = []
        if scores.alphamissense_class is not None:
            parts.append(scores.alphamissense_class)
        if scores.alphamissense_score is not None:
            parts.append(f"({scores.alphamissense_score:.3f})")
        if parts:
            result["AlphaMissense"] = " ".join(parts)

    if scores.spliceai_max is not None:
        # Identify which delta score is the max
        label = _spliceai_max_label(scores)
        result["SpliceAI"] = f"{scores.spliceai_max:.2f} ({label})"

    if scores.sift is not None:
        result["SIFT"] = scores.sift

    if scores.polyphen is not None:
        result["PolyPhen"] = scores.polyphen

    return result


def _spliceai_max_label(scores: InSilicoScores) -> str:
    """Return a short label for which SpliceAI delta score is the max."""
    mapping = {
        "AG": scores.spliceai_ds_ag,
        "AL": scores.spliceai_ds_al,
        "DG": scores.spliceai_ds_dg,
        "DL": scores.spliceai_ds_dl,
    }
    best_label = "AG"
    best_val = -1.0
    for label, val in mapping.items():
        if val is not None and val > best_val:
            best_val = val
            best_label = label
    return best_label
