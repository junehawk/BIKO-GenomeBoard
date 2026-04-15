"""CI smoke — exercises A1/A2 curate-then-narrate → fallback with all external dependencies stubbed.

This test fills the end-to-end regression gap between the Phase A unit tests
and the real Ollama-backed Clinical Board run. CI (ubuntu-latest) has no
Ollama and no CIViC DB, so we stub:

    * ``scripts.clinical.oncokb_client.annotate_protein_change`` — returns a
      canned KRAS G12D row set (sotorasib + adagrasib, NO futibatinib) so the
      real ``curate_treatments`` deterministic merge path runs end-to-end.
    * ``scripts.clinical_board.agents.board_chair.BoardChair.synthesize`` —
      returns a ``CancerBoardOpinion`` with ``treatment_options=[]`` so the
      runner's fallback condition (``curated_nonempty and pre_scrub_count == 0``)
      fires and ``template_renderer_chair.render_from_curated`` is exercised.
    * ``runner.OllamaClient`` / ``runner._load_agents`` — bypass live LLM
      dependency so the whole suite runs in <5s with zero network.

Guarded invariants (patient safety):

    1. ``futibatinib`` MUST NOT appear anywhere in the serialised opinion.
    2. Fallback headline ``"No variant-specific treatment recommendations …"``
       is rendered, confirming the deterministic path ran.
    3. ``confidence == "low"`` and ``agent_consensus == "fallback"`` — the
       fallback tombstones that tell downstream render.py this opinion came
       from the safe deterministic path, not the LLM.
    4. Curated rows (sotorasib, adagrasib) survive into
       ``treatment_options`` via the fallback renderer.
    5. The returned opinion is a ``CancerBoardOpinion`` instance — no Nones,
       no silent early-return from the runner.
"""

from __future__ import annotations

import dataclasses
import json
import time
from unittest.mock import MagicMock


def _kras_g12d_variant() -> dict:
    """Minimal but complete KRAS G12D variant dict for the board runner."""
    return {
        "gene": "KRAS",
        "hgvsp": "p.Gly12Asp",
        "chrom": "12",
        "pos": 25398284,
        "ref": "C",
        "alt": "T",
        "variant": "12:25398284:C:T",
        "classification": "Likely Pathogenic",
        "tier": "Tier I",
        "consequence": "missense_variant",
        "cancer_gene_type": "Oncogene",
        "oncokb_level": "1",
        "variant_type": "SNV",
    }


def _oncokb_stub_sotorasib_adagrasib(gene: str, alteration: str, **kwargs) -> list[dict]:
    """Canned OncoKB response — sotorasib + adagrasib, explicitly NO futibatinib.

    The task description singles out futibatinib as the audit-trail failure
    case (Foundation One fixture 1); its absence from this stub is a
    deliberate guard that the hallucination guardrail isn't trivially
    circumvented by the stub leaking the banned drug.
    """
    return [
        {
            "drug": "Sotorasib",
            "level": "A",
            "pmids": ["32955176"],
            "disease": "Non-Small Cell Lung Cancer",
            "significance": "sensitivity",
            "therapy_ids": "sotorasib",
            "raw_row": {"gene": gene, "alteration": alteration},
        },
        {
            "drug": "Adagrasib",
            "level": "A",
            "pmids": ["35658005"],
            "disease": "Non-Small Cell Lung Cancer",
            "significance": "sensitivity",
            "therapy_ids": "adagrasib",
            "raw_row": {"gene": gene, "alteration": alteration},
        },
    ]


def test_ci_ai_board_curate_then_narrate_fallback(monkeypatch):
    """CI smoke — runner → curator → chair → scrubber → fallback, no Ollama."""
    from scripts.clinical import oncokb_client
    from scripts.clinical_board import runner as runner_mod
    from scripts.clinical_board.agents import board_chair as chair_mod
    from scripts.clinical_board.models import BOARD_DISCLAIMER, CancerBoardOpinion

    # 1) Stub OncoKB HTTP — real curate_treatments() will consume these rows.
    monkeypatch.setattr(
        oncokb_client,
        "annotate_protein_change",
        _oncokb_stub_sotorasib_adagrasib,
        raising=True,
    )

    # 2) Stub Ollama client — is_available=True, has_model=True, generate_json
    #    should never actually be called in this test because _load_chair is
    #    also stubbed, but we set a safe default just in case.
    fake_client = MagicMock()
    fake_client.is_available.return_value = True
    fake_client.has_model.return_value = True
    fake_client.generate_json.return_value = {"treatment_options": []}
    monkeypatch.setattr(runner_mod, "OllamaClient", lambda *a, **kw: fake_client, raising=True)

    # 3) Zero domain agents — skip all TherapeuticTarget/TumorGenomics/etc.
    #    LLM calls entirely.
    monkeypatch.setattr(runner_mod, "_load_agents", lambda *a, **kw: [], raising=True)

    # 4) Stub BoardChair.synthesize — return a CancerBoardOpinion with an
    #    EMPTY treatment_options list. This forces ``pre_scrub_count == 0``
    #    in runner.py:208, which combined with ``curated_nonempty == True``
    #    (from the OncoKB stub above) triggers the
    #    template_renderer_chair.render_from_curated fallback path.
    def _stub_synthesize(self, briefing, opinions, curated_treatments=None, mode="rare-disease"):
        return CancerBoardOpinion(
            therapeutic_headline="stub — LLM produced nothing usable",
            therapeutic_implications="stub body",
            therapeutic_evidence="stub evidence",
            treatment_options=[],
            agent_opinions=list(opinions or []),
            agent_consensus="split",
            confidence="low",
            disclaimer=BOARD_DISCLAIMER,
        )

    monkeypatch.setattr(chair_mod.BoardChair, "synthesize", _stub_synthesize, raising=True)

    # 5) Run the board end-to-end.
    report_data = {
        "sample_id": "CI-SMOKE-KRAS-G12D",
        "variants": [_kras_g12d_variant()],
        "summary": {"total": 1},
    }

    t0 = time.time()
    opinion = runner_mod.run_clinical_board(report_data, mode="cancer")
    elapsed = time.time() - t0

    # ── Invariant 0: runner returned a valid CancerBoardOpinion ─────────
    assert opinion is not None, "run_clinical_board returned None — Ollama stub mis-wired"
    assert isinstance(opinion, CancerBoardOpinion), f"expected CancerBoardOpinion, got {type(opinion).__name__}"

    # ── Invariant 1: futibatinib MUST NOT appear anywhere ──────────────
    serialised = json.dumps(dataclasses.asdict(opinion), ensure_ascii=False, default=str).lower()
    assert "futibatinib" not in serialised, (
        "Patient-safety invariant failed: 'futibatinib' leaked into the "
        "cancer board opinion. The OncoKB stub does not return futibatinib "
        "and the LLM stub returned nothing — any appearance here is a "
        "hallucination the scrubber+fallback failed to catch."
    )

    # ── Invariant 2: fallback headline is rendered ─────────────────────
    assert "no variant-specific treatment recommendations" in opinion.therapeutic_headline.lower(), (
        f"Fallback path did not fire — therapeutic_headline={opinion.therapeutic_headline!r}. "
        "Expected template_renderer_chair.render_from_curated output."
    )

    # ── Invariant 3: confidence + consensus tombstones ─────────────────
    assert opinion.confidence == "low", f"fallback opinion must have confidence='low', got {opinion.confidence!r}"
    assert opinion.agent_consensus == "fallback", (
        f"fallback opinion must have agent_consensus='fallback', got {opinion.agent_consensus!r}"
    )

    # ── Invariant 4: curated rows survive via fallback renderer ────────
    assert len(opinion.treatment_options) >= 1, "fallback produced zero treatment_options despite curated rows"
    drugs = {(row.get("drug") or "").lower() for row in opinion.treatment_options}
    assert drugs & {"sotorasib", "adagrasib"}, f"fallback treatment_options missing curated drugs — got drugs={drugs!r}"

    # ── Invariant 5: runtime budget for CI ────────────────────────────
    assert elapsed < 5.0, f"CI smoke exceeded 5s budget: {elapsed:.2f}s"
