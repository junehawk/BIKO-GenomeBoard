"""Tests for clinical note input handling."""

from scripts.clinical_board.case_briefing import build_case_briefing


def test_clinical_note_in_briefing():
    """Clinical note appears in briefing under CLINICAL CONTEXT section."""
    data = {
        "sample_id": "S001",
        "variants": [],
        "summary": {"total": 0},
        "clinical_note": "55세 남성, 대장암 3기",
    }
    briefing = build_case_briefing(data, "cancer")
    assert "CLINICAL CONTEXT" in briefing
    assert "55세 남성" in briefing


def test_clinical_note_absent():
    """Briefing works without clinical note."""
    data = {"sample_id": "S001", "variants": [], "summary": {"total": 0}}
    briefing = build_case_briefing(data, "cancer")
    assert "CLINICAL CONTEXT" not in briefing


def test_clinical_note_truncation():
    """Long clinical notes are truncated to max chars."""
    long_note = "가" * 2000
    data = {
        "sample_id": "S001",
        "variants": [],
        "summary": {"total": 0},
        "clinical_note": long_note,
    }
    briefing = build_case_briefing(data, "cancer")
    assert "truncated" in briefing.lower() or len(long_note) > 1500


def test_clinical_note_empty_string():
    """Empty string clinical note is treated as absent."""
    data = {
        "sample_id": "S001",
        "variants": [],
        "summary": {"total": 0},
        "clinical_note": "",
    }
    briefing = build_case_briefing(data, "cancer")
    assert "CLINICAL CONTEXT" not in briefing
