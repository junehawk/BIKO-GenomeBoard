"""Narrative scrubber — strips drug hallucinations + invalid (curated_id, variant_key) rows.

Unit tests for the A2 patient-safety post-processor. The scrubber walks the
full ``CancerBoardOpinion`` dataclass and drops any treatment row whose
``(curated_id, variant_key)`` pair did not come from the curator output,
then scrubs banned drug tokens out of every prose field.
"""

from __future__ import annotations

import dataclasses
import json

from scripts.clinical_board.models import AgentOpinion, CancerBoardOpinion


def _stub_curated(drug: str, variant_key: str, curated_id: str):
    """Build a tiny stand-in for a CuratedTreatment without importing it."""

    class _Row:
        pass

    r = _Row()
    r.drug = drug
    r.curated_id = curated_id
    r.variant_key = variant_key
    r.target = ""
    r.evidence_level = "A"
    r.source = "oncokb"
    r.pmids = ["1"]
    r.disease_context = ""
    r.significance = "sensitivity"
    r.therapy_ids = ""
    r.raw_row = {}
    return r


def _opinion_with_rows(rows):
    return CancerBoardOpinion(
        therapeutic_headline="headline",
        therapeutic_implications="implications",
        therapeutic_evidence="evidence",
        treatment_options=rows,
        actionable_findings=["finding 1"],
        clinical_actions=["action 1"],
    )


def test_validate_pair_match_ok():
    from scripts.clinical_board.narrative_scrubber import validate_treatment_option

    curated = {"12:25:C:T": [_stub_curated("Sotorasib", "12:25:C:T", "cid-sot")]}
    row = {"drug": "Sotorasib", "curated_id": "cid-sot", "variant_key": "12:25:C:T"}
    assert validate_treatment_option(row, curated) is True


def test_validate_pair_rejects_missing_curated_id():
    from scripts.clinical_board.narrative_scrubber import validate_treatment_option

    curated = {"12:25:C:T": [_stub_curated("Sotorasib", "12:25:C:T", "cid-sot")]}
    row = {"drug": "Sotorasib", "curated_id": "", "variant_key": "12:25:C:T"}
    assert validate_treatment_option(row, curated) is False


def test_validate_accepts_known_curated_id_regardless_of_variant_key():
    """v2.5.5: validation looks up curated_id only.

    A wrong ``variant_key`` no longer fails validation — the scrubber
    backfills the curator's authoritative value.  ``validate_treatment_option``
    therefore returns True so long as the ``curated_id`` is real.
    """
    from scripts.clinical_board.narrative_scrubber import validate_treatment_option

    curated = {
        "7:55:T:G": [_stub_curated("Osimertinib", "7:55:T:G", "cid-osi")],
        "17:76:G:A": [],
    }
    bad = {"drug": "Osimertinib", "curated_id": "cid-osi", "variant_key": "17:76:G:A"}
    assert validate_treatment_option(bad, curated) is True


def test_validate_rejects_unknown_curated_id():
    from scripts.clinical_board.narrative_scrubber import validate_treatment_option

    curated = {"7:55:T:G": [_stub_curated("Osimertinib", "7:55:T:G", "cid-osi")]}
    bad = {"drug": "Osimertinib", "curated_id": "cid-fabricated", "variant_key": "7:55:T:G"}
    assert validate_treatment_option(bad, curated) is False


def test_scrub_drops_row_with_fabricated_curated_id():
    from scripts.clinical_board.narrative_scrubber import scrub_opinion

    curated = {"12:25:C:T": [_stub_curated("Sotorasib", "12:25:C:T", "cid-sot")]}
    rows = [
        {
            "drug": "Sotorasib",
            "curated_id": "cid-sot",
            "variant_key": "12:25:C:T",
            "evidence_level": "A",
            "resistance_notes": "",
        },
        {
            "drug": "Futibatinib",
            "curated_id": "fabricated",
            "variant_key": "12:25:C:T",
            "evidence_level": "A",
            "resistance_notes": "",
        },
    ]
    op = _opinion_with_rows(rows)
    scrub_opinion(op, curated)
    drugs = [r["drug"].lower() for r in op.treatment_options]
    assert "sotorasib" in drugs
    assert "futibatinib" not in drugs


def test_scrub_backfills_variant_key_when_llm_emits_wrong_one():
    """v2.5.5: scrubber overwrites the LLM's ``variant_key`` with the curator value.

    Pre-v2.5.5 this row was dropped because the (cid, variant_key) pair did not
    match the curator. The new behaviour reassigns the row to the variant the
    curator originally bound to that curated_id, neutralising both LLM
    token-confusion miscopies (the v2.5.5 motivating bug, where supergemma4
    rewrote ``chr13:32356550:C:T`` as ``chr13:32356 own:C:T``) and paste-attack
    attempts.
    """
    from scripts.clinical_board.narrative_scrubber import scrub_opinion

    curated = {
        "7:55:T:G": [_stub_curated("Osimertinib", "7:55:T:G", "cid-osi")],
        "17:76:G:A": [],
    }
    rows = [
        {
            "drug": "Osimertinib",
            "curated_id": "cid-osi",
            "variant_key": "17:76:G:A",  # WRONG — curator says 7:55:T:G
            "evidence_level": "A",
            "resistance_notes": "",
        },
    ]
    op = _opinion_with_rows(rows)
    stats = scrub_opinion(op, curated)
    assert stats["dropped"] == 0
    assert stats["kept"] == 1
    assert len(op.treatment_options) == 1
    # variant_key has been overwritten with the curator's authoritative value
    assert op.treatment_options[0]["variant_key"] == "7:55:T:G"
    assert op.treatment_options[0]["drug"] == "Osimertinib"


def test_scrub_backfills_missing_variant_key():
    """v2.5.5: when the LLM emits curated_id only, scrubber fills variant_key.

    The v2.5.5 Cancer Chair prompt instructs the LLM to omit variant_key
    entirely.  This test exercises that contract end-to-end at the scrubber
    layer.
    """
    from scripts.clinical_board.narrative_scrubber import scrub_opinion

    curated = {"12:25:C:T": [_stub_curated("Sotorasib", "12:25:C:T", "cid-sot")]}
    rows = [
        {
            "drug": "Sotorasib",
            "curated_id": "cid-sot",
            # NO variant_key — Chair prompt says do not emit
            "evidence_level": "A",
            "resistance_notes": "",
        },
    ]
    op = _opinion_with_rows(rows)
    stats = scrub_opinion(op, curated)
    assert stats["kept"] == 1
    assert stats["dropped"] == 0
    assert op.treatment_options[0]["variant_key"] == "12:25:C:T"


def test_scrub_drops_row_with_unknown_curated_id():
    """Hallucinated curated_id is the only valid drop reason in v2.5.5+."""
    from scripts.clinical_board.narrative_scrubber import scrub_opinion

    curated = {"12:25:C:T": [_stub_curated("Sotorasib", "12:25:C:T", "cid-sot")]}
    rows = [
        {
            "drug": "Sotorasib",
            "curated_id": "cid-sot",
            "variant_key": "12:25:C:T",
            "evidence_level": "A",
            "resistance_notes": "",
        },
        {
            "drug": "Futibatinib",
            "curated_id": "cid-fabricated",
            "variant_key": "12:25:C:T",
            "evidence_level": "A",
            "resistance_notes": "",
        },
    ]
    op = _opinion_with_rows(rows)
    stats = scrub_opinion(op, curated)
    assert stats["kept"] == 1
    assert stats["dropped"] == 1
    drugs = [r["drug"].lower() for r in op.treatment_options]
    assert "sotorasib" in drugs
    assert "futibatinib" not in drugs


def test_scrub_strips_banned_drug_from_prose():
    from scripts.clinical_board.narrative_scrubber import scrub_opinion

    curated = {"12:25:C:T": [_stub_curated("Sotorasib", "12:25:C:T", "cid-sot")]}
    rows = [
        {
            "drug": "Futibatinib",
            "curated_id": "fake",
            "variant_key": "12:25:C:T",
            "evidence_level": "A",
            "resistance_notes": "",
        },
    ]
    op = CancerBoardOpinion(
        therapeutic_headline="futibatinib is our headline",
        therapeutic_implications="Consider futibatinib as primary therapy.",
        therapeutic_evidence="futibatinib PMID 123",
        treatment_options=rows,
        clinical_actions=["Start Futibatinib 20mg"],
    )
    scrub_opinion(op, curated)
    payload = json.dumps(dataclasses.asdict(op), ensure_ascii=False).lower()
    assert "futibatinib" not in payload


def test_scrub_preserves_allowed_drug_mentioned_in_prose():
    from scripts.clinical_board.narrative_scrubber import scrub_opinion

    curated = {"12:25:C:T": [_stub_curated("Sotorasib", "12:25:C:T", "cid-sot")]}
    rows = [
        {
            "drug": "Sotorasib",
            "curated_id": "cid-sot",
            "variant_key": "12:25:C:T",
            "evidence_level": "A",
            "resistance_notes": "",
        },
    ]
    op = CancerBoardOpinion(
        therapeutic_headline="Sotorasib candidate",
        therapeutic_implications="Sotorasib is supported by curated evidence.",
        treatment_options=rows,
    )
    scrub_opinion(op, curated)
    assert "Sotorasib" in op.therapeutic_headline
    assert "sotorasib" in op.therapeutic_implications.lower()


def test_scrub_walks_agent_opinions():
    from scripts.clinical_board.narrative_scrubber import scrub_opinion

    curated = {"12:25:C:T": [_stub_curated("Sotorasib", "12:25:C:T", "cid-sot")]}
    agent_op = AgentOpinion(
        agent_name="Clinical Evidence Analyst",
        domain="clinical_evidence",
        findings=[{"finding": "Futibatinib may be effective", "evidence": "bogus", "confidence": "low"}],
        recommendations=["Try Futibatinib"],
        concerns=["Futibatinib side effects"],
    )
    rows = [
        {
            "drug": "Futibatinib",
            "curated_id": "fake",
            "variant_key": "12:25:C:T",
            "evidence_level": "A",
            "resistance_notes": "",
        },
    ]
    op = CancerBoardOpinion(treatment_options=rows, agent_opinions=[agent_op])
    scrub_opinion(op, curated)
    payload = json.dumps(dataclasses.asdict(op), ensure_ascii=False).lower()
    assert "futibatinib" not in payload


def test_scrub_empty_curated_drops_all_rows():
    from scripts.clinical_board.narrative_scrubber import scrub_opinion

    curated: dict = {}
    rows = [
        {
            "drug": "Osimertinib",
            "curated_id": "anything",
            "variant_key": "X:1:A:T",
            "evidence_level": "A",
            "resistance_notes": "",
        },
    ]
    op = _opinion_with_rows(rows)
    scrub_opinion(op, curated)
    assert op.treatment_options == []


def test_scrub_returns_stats():
    from scripts.clinical_board.narrative_scrubber import scrub_opinion

    curated = {"12:25:C:T": [_stub_curated("Sotorasib", "12:25:C:T", "cid-sot")]}
    rows = [
        {"drug": "Sotorasib", "curated_id": "cid-sot", "variant_key": "12:25:C:T"},
        {"drug": "Futibatinib", "curated_id": "fake", "variant_key": "12:25:C:T"},
    ]
    op = _opinion_with_rows(rows)
    stats = scrub_opinion(op, curated)
    assert stats["kept"] == 1
    assert stats["dropped"] == 1
