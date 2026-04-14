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


def test_validate_pair_rejects_cross_variant_paste():
    from scripts.clinical_board.narrative_scrubber import validate_treatment_option

    curated = {
        "7:55:T:G": [_stub_curated("Osimertinib", "7:55:T:G", "cid-osi")],
        "17:76:G:A": [],
    }
    bad = {"drug": "Osimertinib", "curated_id": "cid-osi", "variant_key": "17:76:G:A"}
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


def test_scrub_drops_row_with_mismatched_variant_key():
    from scripts.clinical_board.narrative_scrubber import scrub_opinion

    curated = {
        "7:55:T:G": [_stub_curated("Osimertinib", "7:55:T:G", "cid-osi")],
        "17:76:G:A": [],
    }
    rows = [
        {
            "drug": "Osimertinib",
            "curated_id": "cid-osi",
            "variant_key": "17:76:G:A",
            "evidence_level": "A",
            "resistance_notes": "",
        },
    ]
    op = _opinion_with_rows(rows)
    scrub_opinion(op, curated)
    assert op.treatment_options == []


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
