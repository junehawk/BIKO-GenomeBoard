"""PM1 protein-domain hotspot table lookup (A3-core).

Tests that :func:`collect_additional_evidence` consults
``data/pm1_hotspot_domains.json`` so PM1 can fire even when VEP left the
``domains`` field empty (ClinGen TP53 VCEP / AMP hotspot-table pathway).
"""
from __future__ import annotations

import pytest

from scripts.classification.evidence_collector import collect_additional_evidence
from scripts.common.models import Variant


def _missense(gene: str, hgvsp: str, *, domains: str = "") -> Variant:
    v = Variant(
        chrom="1",
        pos=1,
        ref="A",
        alt="T",
        gene=gene,
        consequence="missense_variant",
        hgvsp=hgvsp,
    )
    v.domains = domains
    return v


def test_tp53_r249m_fires_pm1_via_table():
    """TP53 R249M sits in the L3 loop hotspot range 245-249."""
    variant = _missense("TP53", "p.Arg249Met", domains="")
    codes = collect_additional_evidence(variant)
    assert "PM1" in codes


def test_tp53_r100q_no_pm1():
    """R100 is outside every TP53 hotspot range and has no VEP domain."""
    variant = _missense("TP53", "p.Arg100Gln", domains="")
    codes = collect_additional_evidence(variant)
    assert "PM1" not in codes
    assert "PM1_Supporting" not in codes


def test_kras_g12d_fires_pm1_via_table():
    """KRAS G12D is the canonical AMP Tier I example."""
    variant = _missense("KRAS", "p.Gly12Asp", domains="")
    codes = collect_additional_evidence(variant)
    assert "PM1" in codes


def test_kras_a146t_fires_pm1_supporting():
    """KRAS 146 is marked as supporting in the hotspot table."""
    variant = _missense("KRAS", "p.Ala146Thr", domains="")
    codes = collect_additional_evidence(variant)
    assert "PM1_Supporting" in codes
    assert "PM1" not in codes


def test_pm1_fires_with_biko_formatted_consequence_label():
    """v2.2 regression — real pipeline stores consequence as the BIKO-formatted
    short label (e.g. "Missense") via scripts/intake/parse_annotation.py::
    format_consequence, not the raw VEP SO term. ``_get_consequence`` must
    canonicalise both forms or PM1 silently never fires on real data.

    This is the root cause of TP53 R249M staying VUS with empty acmg_codes
    in the post-Phase-B showcase regeneration even though A3 was meant to
    fix it.
    """
    variant = Variant(
        chrom="17",
        pos=7674217,
        ref="C",
        alt="A",
        gene="TP53",
        consequence="Missense",  # BIKO formatted form, not SO term
        hgvsp="ENSP00000269305.4:p.Arg249Met",
    )
    variant.domains = ""
    codes = collect_additional_evidence(variant)
    assert "PM1" in codes


def test_existing_domain_path_still_fires_pm1():
    """Regression: the legacy VEP DOMAINS-based path still works."""
    variant = _missense("FOOBAR", "p.Val100Leu", domains="PFAM:PF00001")
    codes = collect_additional_evidence(variant)
    assert "PM1" in codes


def test_non_hotspot_gene_without_domain_no_pm1():
    variant = _missense("FOOBAR", "p.Val100Leu", domains="")
    codes = collect_additional_evidence(variant)
    assert "PM1" not in codes
    assert "PM1_Supporting" not in codes


@pytest.mark.parametrize(
    "hgvsp",
    [
        "p.Arg249Met",  # 3-letter substitution
        "p.R249M",      # 1-letter substitution
        "p.Arg249fs",   # frameshift — position only
        " p.Arg249Met ",  # whitespace
    ],
)
def test_hgvsp_fuzz_r249_fires(hgvsp):
    variant = _missense("TP53", hgvsp, domains="")
    codes = collect_additional_evidence(variant)
    assert "PM1" in codes, f"expected PM1 for hgvsp={hgvsp!r}"


def test_hgvsp_empty_no_crash():
    variant = _missense("TP53", "", domains="")
    codes = collect_additional_evidence(variant)
    assert "PM1" not in codes


def test_hgvsp_met1_start_lost_tp53_no_pm1():
    """p.Met1? is a start-lost — position 1 is not in any TP53 hotspot."""
    variant = _missense("TP53", "p.Met1?", domains="")
    codes = collect_additional_evidence(variant)
    assert "PM1" not in codes


def test_supporting_does_not_upgrade_with_domain_hit():
    """Domain-only hit should still fire the full PM1 even in a gene whose
    table entry for the position is only ``supporting`` — the two paths are
    additive, domain hit wins."""
    variant = _missense("KRAS", "p.Ala146Thr", domains="PFAM:PF00071")
    codes = collect_additional_evidence(variant)
    assert "PM1" in codes
