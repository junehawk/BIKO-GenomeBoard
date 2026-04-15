"""Template renderer chair — deterministic fallback when LLM synthesis fails validation.

When the Board Chair returns zero valid ``(curated_id, variant_key)`` pairs,
``runner.run_clinical_board`` falls back to the template renderer, which
emits a ``CancerBoardOpinion`` directly from the curator output so the
clinician still sees the validated therapy table.
"""

from __future__ import annotations

from scripts.clinical_board.models import AgentOpinion, CancerBoardOpinion


def _stub_row(drug, variant_key, curated_id, level="A", pmids=None, disease=""):
    class _R:
        pass

    r = _R()
    r.drug = drug
    r.curated_id = curated_id
    r.variant_key = variant_key
    r.target = ""
    r.evidence_level = level
    r.source = "oncokb"
    r.pmids = pmids or ["12345"]
    r.disease_context = disease
    r.significance = "sensitivity"
    r.therapy_ids = ""
    r.raw_row = {}
    return r


def test_fallback_produces_nonempty_treatment_options():
    from scripts.clinical_board.template_renderer_chair import render_from_curated

    curated = {
        "12:25:C:T": [
            _stub_row("Sotorasib", "12:25:C:T", "cid-sot", "A", ["32955176"], "NSCLC"),
            _stub_row("Adagrasib", "12:25:C:T", "cid-adag", "A", ["35658005"], "NSCLC"),
        ]
    }
    opinion = render_from_curated(curated, agent_opinions=[])
    assert isinstance(opinion, CancerBoardOpinion)
    assert len(opinion.treatment_options) == 2
    drugs = [r["drug"] for r in opinion.treatment_options]
    assert "Sotorasib" in drugs
    assert "Adagrasib" in drugs


def test_fallback_confidence_is_low():
    from scripts.clinical_board.template_renderer_chair import render_from_curated

    curated = {"12:25:C:T": [_stub_row("Sotorasib", "12:25:C:T", "cid-sot")]}
    opinion = render_from_curated(curated)
    assert opinion.confidence == "low"
    assert opinion.agent_consensus == "fallback"


def test_fallback_preserves_pmids_in_hydrated():
    from scripts.clinical_board.template_renderer_chair import render_from_curated

    curated = {"12:25:C:T": [_stub_row("Sotorasib", "12:25:C:T", "cid-sot", "A", ["32955176", "35110016"], "NSCLC")]}
    opinion = render_from_curated(curated)
    row = opinion.treatment_options[0]
    assert row["curated_id"] == "cid-sot"
    assert row["variant_key"] == "12:25:C:T"
    hydrated = row.get("_hydrated", {})
    assert "32955176" in hydrated.get("pmids", [])
    assert "35110016" in hydrated.get("pmids", [])
    assert hydrated.get("disease_context") == "NSCLC"


def test_fallback_empty_curated_returns_no_findings_opinion():
    from scripts.clinical_board.template_renderer_chair import render_from_curated

    opinion = render_from_curated({})
    assert opinion.treatment_options == []
    assert opinion.confidence == "low"


def test_fallback_headline_does_not_imply_recommended_therapies():
    """Regression: the fallback headline must not read like 'N treatment
    options are recommended'. An early reviewer mis-read '211 curated
    therapy option(s) — LLM synthesis unavailable' as 'we have 211
    recommended treatments'. The headline must explicitly disclaim that
    the table is an evidence library, not a list of recommended
    treatments."""
    from scripts.clinical_board.template_renderer_chair import render_from_curated

    curated = {"12:25:C:T": [_stub_row(f"Drug{i}", "12:25:C:T", f"cid-{i}", "A") for i in range(50)]}
    opinion = render_from_curated(curated)

    # Forbidden phrasing from the pre-fix headline
    assert "curated therapy option" not in opinion.therapeutic_headline.lower()
    assert "treatment option" not in opinion.therapeutic_headline.lower()

    # Required framing
    assert "no variant-specific treatment recommendation" in opinion.therapeutic_headline.lower()
    assert "evidence library" in opinion.therapeutic_headline.lower()

    # Body must tell the reader these rows are NOT recommended treatments
    body = opinion.therapeutic_implications.lower()
    assert "not" in body and ("recommended treatment" in body or "research reference" in body)
    assert "qualified oncologist" in body or "clinical interpretation" in body

    # Evidence field must also carry the disclaimer
    assert "not" in opinion.therapeutic_evidence.lower()
    assert (
        "research reference" in opinion.therapeutic_evidence.lower()
        or "evidence" in opinion.therapeutic_evidence.lower()
    )


def test_fallback_empty_curated_has_no_matched_messaging():
    from scripts.clinical_board.template_renderer_chair import render_from_curated

    opinion = render_from_curated({})
    assert opinion.treatment_options == []
    assert "no matched" in opinion.therapeutic_headline.lower()
    assert "clinical interpretation" in opinion.therapeutic_implications.lower()


def test_fallback_carries_agent_opinions_through():
    from scripts.clinical_board.template_renderer_chair import render_from_curated

    agents = [AgentOpinion(agent_name="Clinical Evidence Analyst", domain="clinical_evidence")]
    curated = {"12:25:C:T": [_stub_row("Sotorasib", "12:25:C:T", "cid-sot")]}
    opinion = render_from_curated(curated, agent_opinions=agents)
    assert opinion.agent_opinions == agents
