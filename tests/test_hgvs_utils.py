"""HGVSp normalisation utilities — fuzz + edge cases for A1 curator.

Covers ``scripts.common.hgvs_utils`` extensions:
- ``normalize_hgvsp_for_civic`` — CIViC variant-name mapping (``p.Gly12Asp`` → ``G12D``)
- ``extract_protein_position`` — residue extraction across 1-letter / 3-letter / edge forms
"""
from __future__ import annotations

import pytest


def test_normalize_hgvsp_three_letter_missense():
    from scripts.common.hgvs_utils import normalize_hgvsp_for_civic
    assert normalize_hgvsp_for_civic("p.Gly12Asp") == "G12D"
    assert normalize_hgvsp_for_civic("p.Val600Glu") == "V600E"
    assert normalize_hgvsp_for_civic("p.Arg249Met") == "R249M"


def test_normalize_hgvsp_one_letter_missense():
    from scripts.common.hgvs_utils import normalize_hgvsp_for_civic
    assert normalize_hgvsp_for_civic("p.G12D") == "G12D"
    assert normalize_hgvsp_for_civic("p.V600E") == "V600E"
    assert normalize_hgvsp_for_civic("p.R249M") == "R249M"


def test_normalize_hgvsp_whitespace_wrapped():
    from scripts.common.hgvs_utils import normalize_hgvsp_for_civic
    assert normalize_hgvsp_for_civic("  p.Gly12Asp  ") == "G12D"
    assert normalize_hgvsp_for_civic("\tp.V600E\n") == "V600E"


def test_normalize_hgvsp_unparseable_returns_none():
    from scripts.common.hgvs_utils import normalize_hgvsp_for_civic
    assert normalize_hgvsp_for_civic("") is None
    assert normalize_hgvsp_for_civic(None) is None
    assert normalize_hgvsp_for_civic("not hgvsp at all") is None
    # p.Met1? is a start-lost notation — not a simple substitution
    assert normalize_hgvsp_for_civic("p.Met1?") is None
    # frameshift — no simple AA-to-AA mapping
    assert normalize_hgvsp_for_civic("p.Arg249fs") is None


@pytest.mark.parametrize(
    "hgvsp,expected",
    [
        ("p.Arg249Met", 249),
        ("p.R249M", 249),
        ("p.Gly12Asp", 12),
        ("p.G12D", 12),
        ("p.Val600Glu", 600),
        ("  p.Arg249Met  ", 249),
        ("p.Met1?", 1),
        ("p.Arg249fs", 249),
        ("", None),
        (None, None),
        ("not hgvsp", None),
    ],
)
def test_extract_protein_position_fuzz(hgvsp, expected):
    from scripts.common.hgvs_utils import extract_protein_position
    assert extract_protein_position(hgvsp) == expected


def test_legacy_hgvsp_to_civic_variant_still_works():
    """Backward-compat: existing callers (query_civic) keep working."""
    from scripts.common.hgvs_utils import hgvsp_to_civic_variant
    assert hgvsp_to_civic_variant("p.Gly12Asp") == "G12D"
    assert hgvsp_to_civic_variant(None) is None
