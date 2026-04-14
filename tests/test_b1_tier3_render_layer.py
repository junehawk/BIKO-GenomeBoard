"""Render-layer regression for v2.2 B1 (qa-engineer fixture 4).

Pipes a mixed variant set through ``runner.run_clinical_board()`` and the
Clinical Board Jinja render, then asserts that the consequence-gated
selector result is what actually reaches the final HTML.

- KRAS G12D missense Tier III (hotspot) MUST appear in the rendered HTML
  because the selector admits it.
- NOTCH2 intronic Tier III MUST NOT appear because the B1 consequence
  gate rejects intron_variant consequences outside the P/LP bypass.

This closes the Foundation One leak end-to-end (not just at the selector
boundary) — if the gate ever regresses, the rendered Board HTML will
start carrying NOTCH2 again.
"""
from __future__ import annotations

from unittest.mock import MagicMock

import pytest


def _build_mixed_cancer_variants() -> list[dict]:
    """Build the minimum mixed set the render-layer test needs."""
    return [
        {
            "gene": "KRAS",
            "variant": "12:25398284:C:T",
            "chrom": "12",
            "pos": 25398284,
            "ref": "C",
            "alt": "T",
            "classification": "VUS",
            "tier": "Tier III",
            "hgvsp": "p.Gly12Asp",
            "hgvsc": "c.35G>A",
            "consequence": "missense_variant",
            "cancer_gene_type": "Oncogene",
            "oncokb_level": "1",
            "variant_type": "SNV",
        },
        {
            "gene": "NOTCH2",
            "variant": "1:120478000:A:G",
            "chrom": "1",
            "pos": 120478000,
            "ref": "A",
            "alt": "G",
            "classification": "VUS",
            "tier": "Tier III",
            "hgvsp": "",
            "hgvsc": "c.123+500A>G",
            "consequence": "intron_variant",
            "cancer_gene_type": "",
            "oncokb_level": "1",
            "variant_type": "SNV",
        },
    ]


def test_b1_render_layer_notch2_intronic_excluded_kras_present(monkeypatch):
    runner_mod = pytest.importorskip("scripts.clinical_board.runner")
    render_mod = pytest.importorskip("scripts.clinical_board.render")
    selector_mod = pytest.importorskip("scripts.clinical_board.variant_selector")

    # Route the hotspot / cancer-gene checks through deterministic fakes so
    # the test is hermetic: KRAS G12 hotspot, both genes in the cancer gene
    # list. Without this, is_hotspot / is_cancer_gene would hit the CIViC DB.
    monkeypatch.setattr(
        selector_mod,
        "is_hotspot",
        lambda gene, pos: (gene, pos) == ("KRAS", 12),
    )
    monkeypatch.setattr(
        selector_mod,
        "is_cancer_gene",
        lambda gene: gene in {"KRAS", "NOTCH2"},
    )

    # No domain agents, no curator rows — we only want to exercise the
    # selector → briefing → chair → render path.
    monkeypatch.setattr(runner_mod, "_load_agents", lambda *a, **kw: [], raising=False)
    monkeypatch.setattr(
        runner_mod, "curate_treatments", lambda variants, **kw: {}, raising=False
    )

    # Fake OllamaClient: the Board Chair calls .generate_json(prompt=...)
    # with a prompt that embeds the case briefing (which is built from the
    # already-gated variant list). Echo the briefing back into
    # therapeutic_implications so the render layer surfaces whichever gene
    # names actually made it past the selector.
    def fake_generate_json(*, model, prompt, system, temperature, **_kw):
        return {
            "therapeutic_headline": "B1 render-layer probe",
            "therapeutic_implications": prompt,
            "therapeutic_evidence": "",
            "treatment_options": [],
            "actionable_findings": [],
            "clinical_actions": [],
            "immunotherapy_eligibility": "",
            "confidence": "moderate",
            "agent_consensus": "majority",
            "dissenting_opinions": [],
            "monitoring_plan": [],
        }

    fake_client = MagicMock()
    fake_client.is_available.return_value = True
    fake_client.has_model.return_value = True
    fake_client.generate_json.side_effect = fake_generate_json
    monkeypatch.setattr(
        runner_mod, "OllamaClient", lambda *a, **kw: fake_client, raising=False
    )

    report_data = {
        "sample_id": "B1-RENDER-FIXTURE-4",
        "variants": _build_mixed_cancer_variants(),
        "summary": {"total": 2},
    }

    opinion = runner_mod.run_clinical_board(report_data, mode="cancer")
    assert opinion is not None, "runner returned None — check fake client wiring"

    # Selector audit trail: KRAS admitted (Tier_III_hotspot), NOTCH2 excluded.
    meta = report_data.get("_board_selection_metadata") or {}
    assert meta.get("selected") == 1
    admitted_genes = {v.get("gene") for v in (report_data.get("_board_variants") or [])}
    assert admitted_genes == {"KRAS"}

    # Render the cancer board opinion to HTML and assert on the final artifact.
    html = render_mod.render_board_opinion_html(opinion, language="en")
    assert "KRAS" in html, (
        "KRAS missense Tier III should reach the rendered Clinical Board HTML — "
        "the consequence gate is over-rejecting."
    )
    assert "NOTCH2" not in html, (
        "NOTCH2 intronic Tier III leaked into the rendered Clinical Board HTML "
        "— B1 consequence gate regressed. Foundation One audit regression."
    )
