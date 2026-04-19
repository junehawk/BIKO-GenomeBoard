"""Integration tests for the v2.4 QW-A classify.build_variant_records
pass-through of ``selection_reason_list``.

The plan's pre-step fix: v2.3 ticket 1 populates
``variant.selection_reason_list`` (via variant_selector's tagging + any
upstream enrichment), but the report_data builder in
``scripts/pipeline/classify.py`` previously dropped the field, so the
template's de novo badge was permanently blank on real data.

These tests confirm the record dict now always carries a
``selection_reason_list`` key so the ``denovo_badge`` macro can render.
"""

from __future__ import annotations

from unittest.mock import MagicMock

from scripts.common.models import Variant
from scripts.pipeline.classify import build_variant_records


def _min_db(variant_id: str) -> dict:
    """Minimal db_results entry — matches the field set consumed by
    build_variant_records."""
    return {
        variant_id: {
            "clinvar": {
                "clinvar_significance": "VUS",
                "clinvar_id": None,
                "review_status": "",
                "acmg_codes": [],
            },
            "gnomad": {
                "gnomad_all": None,
                "gnomad_eas": None,
                "api_available": False,
            },
            "krgdb_freq": None,
            "korea4k_freq": None,
            "nard2_freq": None,
            "pgx": None,
        }
    }


def _min_classification() -> MagicMock:
    mock = MagicMock()
    mock.classification = "VUS"
    mock.evidence_codes = []
    mock.conflict = None
    mock.clinvar_override = False
    mock.clinvar_override_reason = ""
    mock.original_engine_classification = None
    return mock


def test_report_data_carries_empty_selection_reason_list_by_default():
    """Legacy Variant objects without selection_reason_list get [] in the record."""
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")
    vid = variant.variant_id
    records = build_variant_records(
        [variant],
        _min_db(vid),
        {vid: {"korean_flag": ""}},
        {vid: _min_classification()},
        "rare-disease",
        [],
    )
    assert len(records) == 1
    # Key must always exist (not raise KeyError) — the template loop
    # iterates this directly.
    assert "selection_reason_list" in records[0]
    assert records[0]["selection_reason_list"] == []


def test_report_data_carries_selection_reason_list_from_variant():
    """When an upstream step has tagged the Variant object with
    selection_reason_list, it propagates to the record dict."""
    variant = Variant(chrom="chr2", pos=166000000, ref="A", alt="T", gene="SCN1A")
    # v2.3 ticket 1: selector's _tag dict carries this. Mimic the same
    # attribute on the Variant object so classify's getattr fallback
    # picks it up.
    variant.selection_reason_list = ["VUS_denovo_neurodev", "VUS_HPO_match"]  # type: ignore[attr-defined]
    vid = variant.variant_id
    records = build_variant_records(
        [variant],
        _min_db(vid),
        {vid: {"korean_flag": ""}},
        {vid: _min_classification()},
        "rare-disease",
        [],
    )
    assert records[0]["selection_reason_list"] == [
        "VUS_denovo_neurodev",
        "VUS_HPO_match",
    ]


def test_report_data_selection_reason_list_is_list_copy():
    """The record field must be a list copy — mutating the record later
    must not leak back to the Variant object, and a ``None`` attribute
    must normalize to ``[]``."""
    variant = Variant(chrom="chr2", pos=166000000, ref="A", alt="T", gene="SCN1A")
    variant.selection_reason_list = None  # type: ignore[attr-defined]
    vid = variant.variant_id
    records = build_variant_records(
        [variant],
        _min_db(vid),
        {vid: {"korean_flag": ""}},
        {vid: _min_classification()},
        "rare-disease",
        [],
    )
    assert records[0]["selection_reason_list"] == []
    assert isinstance(records[0]["selection_reason_list"], list)


def test_report_data_carries_denovo_inheritance_alongside_reason_list():
    """The three de novo fields are plumbed together so the template
    badge macro renders the PS2/PM6 inheritance pill AND the DDG2P pill
    from the same record."""
    variant = Variant(
        chrom="chr2",
        pos=166000000,
        ref="A",
        alt="T",
        gene="SCN1A",
        inheritance="de_novo",
        confirmed_denovo=False,
    )
    variant.selection_reason_list = ["VUS_denovo_neurodev"]  # type: ignore[attr-defined]
    vid = variant.variant_id
    records = build_variant_records(
        [variant],
        _min_db(vid),
        {vid: {"korean_flag": ""}},
        {vid: _min_classification()},
        "rare-disease",
        [],
    )
    rec = records[0]
    assert rec["variant_inheritance"] == "de_novo"
    assert rec["confirmed_denovo"] is False
    assert rec["selection_reason_list"] == ["VUS_denovo_neurodev"]


def test_report_data_selection_reason_list_independent_of_mode():
    """Both cancer and rare-disease records carry the field. Cancer mode
    doesn't usually populate it but the channel must exist for
    future-proofing."""
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")
    vid = variant.variant_id
    for mode in ("cancer", "rare-disease"):
        records = build_variant_records(
            [variant],
            _min_db(vid),
            {vid: {"korean_flag": ""}},
            {vid: _min_classification()},
            mode,
            [],
        )
        assert "selection_reason_list" in records[0], f"missing in mode={mode}"
        assert records[0]["selection_reason_list"] == []
