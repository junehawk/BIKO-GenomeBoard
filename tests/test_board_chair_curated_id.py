"""Board Chair — curated evidence prompt block and (curated_id, variant_key) binding.

Covers the A2 signature change: ``BoardChair.synthesize(briefing, opinions,
curated_treatments, mode=...)``. The prompt must inject a ``CURATED EVIDENCE``
section the LLM can cite by ``curated_id``, and the runner-level narrative
scrubber must reject any row whose ``(curated_id, variant_key)`` pair was not
emitted by the curator.
"""

from __future__ import annotations

import inspect
from unittest.mock import MagicMock

from scripts.clinical_board.agents.board_chair import BoardChair


def _stub(drug, variant_key, cid):
    class _R:
        pass

    r = _R()
    r.drug = drug
    r.curated_id = cid
    r.variant_key = variant_key
    r.target = ""
    r.evidence_level = "A"
    r.source = "oncokb"
    r.pmids = ["1"]
    r.disease_context = "NSCLC"
    r.significance = "sensitivity"
    r.therapy_ids = ""
    r.raw_row = {}
    return r


def test_synthesize_signature_accepts_curated_treatments():
    sig = inspect.signature(BoardChair.synthesize)
    assert "curated_treatments" in sig.parameters, (
        "A2 signature change missing: BoardChair.synthesize must accept curated_treatments"
    )


def test_prompt_includes_curated_evidence_block():
    chair = BoardChair(client=MagicMock(), model="stub", language="en")
    curated = {"12:25:C:T": [_stub("Sotorasib", "12:25:C:T", "cid-sot")]}
    prompt = chair._build_prompt("briefing", [], curated_treatments=curated)
    assert "CURATED EVIDENCE" in prompt.upper()
    assert "cid-sot" in prompt
    assert "Sotorasib" in prompt
    assert "12:25:C:T" in prompt


def test_prompt_requires_curated_id_in_schema():
    """The CANCER system prompt JSON schema must mention curated_id + variant_key."""
    from scripts.clinical_board.agents.board_chair import CANCER_SYSTEM_PROMPT

    assert "curated_id" in CANCER_SYSTEM_PROMPT
    assert "variant_key" in CANCER_SYSTEM_PROMPT


def test_cross_variant_binding_rejected_by_scrubber():
    """Fixture-2 regression at the scrubber layer (not end-to-end).

    Even when the LLM pastes a valid EGFR curated_id under a TP53 row,
    ``scrub_opinion`` drops the row because the (cid, variant_key) pair
    was never emitted by the curator.
    """
    from scripts.clinical_board.models import CancerBoardOpinion
    from scripts.clinical_board.narrative_scrubber import scrub_opinion

    curated = {
        "7:55:T:G": [_stub("Osimertinib", "7:55:T:G", "cid-osi")],
        "17:76:G:A": [],
    }
    opinion = CancerBoardOpinion(
        therapeutic_headline="TP53 R249M + EGFR L858R",
        treatment_options=[
            {"drug": "Osimertinib", "curated_id": "cid-osi", "variant_key": "17:76:G:A", "evidence_level": "A"},
        ],
    )
    scrub_opinion(opinion, curated)
    assert opinion.treatment_options == [], "cross-variant paste attack not caught by narrative_scrubber"
