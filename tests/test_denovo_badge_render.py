"""Macro-level unit tests for the ``denovo_badge`` Jinja macro (v2.4 QW-A).

Exercises the rare-disease report template's variant-level de novo badge in
isolation so we don't depend on the full report fixture. We extract the
macro source text out of ``templates/rare-disease/report.html`` by its
``{% macro denovo_badge(variant) %}`` delimiter and compile it via
``Environment.from_string`` with ``StrictUndefined`` — any forgotten
``default('')`` guard then raises loudly rather than silently rendering
``None``. This matches the QA reviewer's "macro unit isolation" guidance.

Spec: ``_workspace/v24-engineering/abc_plan_v2.md`` section A.3.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from jinja2 import Environment, StrictUndefined


_TEMPLATE_PATH = Path(__file__).resolve().parents[1] / "templates" / "rare-disease" / "report.html"


def _extract_macro_source(template_text: str, macro_name: str) -> str:
    """Slice the named macro definition out of a larger template.

    Relies on the ``{% macro NAME(...) %}…{% endmacro %}`` syntax. We match
    the first opening delimiter and its matching ``{% endmacro %}``; nested
    macros of the same name are not supported (and the real template does
    not do that).
    """
    pattern = re.compile(
        r"(\{%\s*macro\s+" + re.escape(macro_name) + r"\s*\([^)]*\)\s*%\}.*?\{%\s*endmacro\s*%\})",
        re.DOTALL,
    )
    match = pattern.search(template_text)
    if not match:
        raise AssertionError(f"macro {macro_name!r} not found in {_TEMPLATE_PATH}")
    return match.group(1)


@pytest.fixture(scope="module")
def denovo_badge():
    """Return the ``denovo_badge`` macro bound to a minimal Jinja env.

    ``StrictUndefined`` is intentional — the plan's final test
    (``test_legacy_variant_without_denovo_attrs``) must still pass when a
    variant dict lacks both ``variant_inheritance`` and
    ``selection_reason_list``. If that passes under StrictUndefined it
    proves the ``default('')`` / ``default([])`` guards inside the macro
    are doing real work.
    """
    text = _TEMPLATE_PATH.read_text()
    macro_src = _extract_macro_source(text, "denovo_badge")
    # Append a tiny harness that calls the macro so ``from_string`` has
    # somewhere to bind ``denovo_badge`` at render time.
    env = Environment(autoescape=True, undefined=StrictUndefined)
    tmpl = env.from_string(macro_src + "\n{{ denovo_badge(variant) }}")

    def _render(variant: dict) -> str:
        return tmpl.render(variant=variant)

    return _render


# ---------------------------------------------------------------------------
# variant_inheritance axis (PS2 / PM6 badges)
# ---------------------------------------------------------------------------


def test_denovo_badge_renders_for_de_novo_inheritance(denovo_badge):
    variant = {
        "variant_inheritance": "de_novo",
        "confirmed_denovo": False,
        "selection_reason_list": [],
    }
    html = denovo_badge(variant)
    assert "De novo (PM6, assumed)" in html
    # Should NOT render the confirmed variant badge.
    assert "Confirmed de novo" not in html


def test_confirmed_denovo_badge_ps2(denovo_badge):
    variant = {
        "variant_inheritance": "confirmed_de_novo",
        "confirmed_denovo": True,
        "selection_reason_list": [],
    }
    html = denovo_badge(variant)
    assert "Confirmed de novo (PS2)" in html
    # The elif branch must not also fire.
    assert "PM6, assumed" not in html


def test_denovo_badge_absent_when_inheritance_none(denovo_badge):
    variant = {
        "variant_inheritance": None,
        "confirmed_denovo": False,
        "selection_reason_list": [],
    }
    html = denovo_badge(variant).strip()
    # No inheritance + no reasons → empty output (ignoring whitespace).
    assert "De novo" not in html
    assert "Confirmed" not in html
    assert "Neurodev" not in html
    assert "Splice" not in html


# ---------------------------------------------------------------------------
# selection_reason_list axis (Neurodev / Splice badges)
# ---------------------------------------------------------------------------


def test_selection_reason_empty_list(denovo_badge):
    variant = {
        "variant_inheritance": None,
        "selection_reason_list": [],
    }
    html = denovo_badge(variant).strip()
    assert "Neurodev" not in html
    assert "Splice" not in html


def test_selection_reason_ddg2p_neurodev_badge(denovo_badge):
    variant = {
        "variant_inheritance": None,
        "selection_reason_list": ["VUS_denovo_neurodev"],
    }
    html = denovo_badge(variant)
    assert "Neurodev (DDG2P)" in html
    # Splice badge must not appear for a neurodev-only reason.
    assert "De novo Splice" not in html


def test_selection_reason_denovo_splice_badge(denovo_badge):
    variant = {
        "variant_inheritance": None,
        "selection_reason_list": ["VUS_denovo_splice"],
    }
    html = denovo_badge(variant)
    assert "De novo Splice (SpliceAI)" in html
    # Neurodev badge must not appear for a splice-only reason.
    assert "Neurodev (DDG2P)" not in html


def test_selection_reason_unknown_codes_ignored(denovo_badge):
    """Reasons we don't render (VUS_HPO_match, VUS_hotspot, ...) stay silent."""
    variant = {
        "variant_inheritance": None,
        "selection_reason_list": ["VUS_HPO_match", "VUS_hotspot", "Tier_I"],
    }
    html = denovo_badge(variant).strip()
    assert "Neurodev" not in html
    assert "Splice" not in html
    assert "De novo" not in html


# ---------------------------------------------------------------------------
# Combination axis — both inheritance AND selection reasons
# ---------------------------------------------------------------------------


def test_both_inheritance_and_selection_reason(denovo_badge):
    """A DDG2P-admitted de novo neurodev missense should render BOTH badges."""
    variant = {
        "variant_inheritance": "de_novo",
        "confirmed_denovo": False,
        "selection_reason_list": ["VUS_denovo_neurodev"],
    }
    html = denovo_badge(variant)
    assert "De novo (PM6, assumed)" in html
    assert "Neurodev (DDG2P)" in html


def test_confirmed_denovo_plus_splice_rescue(denovo_badge):
    """Confirmed de novo + SpliceAI splice rescue → PS2 + Splice badge both."""
    variant = {
        "variant_inheritance": "confirmed_de_novo",
        "confirmed_denovo": True,
        "selection_reason_list": ["VUS_denovo_splice"],
    }
    html = denovo_badge(variant)
    assert "Confirmed de novo (PS2)" in html
    assert "De novo Splice (SpliceAI)" in html


def test_denovo_with_hpo_coreason_records_both(denovo_badge):
    """variant_selector collision-point 3 output: neurodev + HPO → neurodev badge only."""
    variant = {
        "variant_inheritance": "de_novo",
        "selection_reason_list": ["VUS_denovo_neurodev", "VUS_HPO_match"],
    }
    html = denovo_badge(variant)
    assert "De novo (PM6, assumed)" in html
    assert "Neurodev (DDG2P)" in html
    # VUS_HPO_match is audit-only — no badge for it.
    assert "HPO match" not in html


# ---------------------------------------------------------------------------
# Legacy / defensive paths — default('') guards
# ---------------------------------------------------------------------------


def test_legacy_variant_without_denovo_attrs(denovo_badge):
    """A variant dict missing both ``variant_inheritance`` and
    ``selection_reason_list`` must render empty output, not raise."""
    variant: dict = {}  # worst-case legacy record
    html = denovo_badge(variant).strip()
    assert "De novo" not in html
    assert "Neurodev" not in html
    assert "Splice" not in html


def test_selection_reason_list_missing_but_inheritance_present(denovo_badge):
    variant = {"variant_inheritance": "de_novo"}  # no selection_reason_list key
    html = denovo_badge(variant)
    assert "De novo (PM6, assumed)" in html


def test_inheritance_missing_but_selection_reason_present(denovo_badge):
    variant = {"selection_reason_list": ["VUS_denovo_neurodev"]}
    html = denovo_badge(variant)
    assert "Neurodev (DDG2P)" in html


def test_empty_string_inheritance_ignored(denovo_badge):
    """An empty-string ``variant_inheritance`` must not match any branch."""
    variant = {"variant_inheritance": "", "selection_reason_list": []}
    html = denovo_badge(variant).strip()
    assert "De novo" not in html
    assert "Confirmed" not in html
