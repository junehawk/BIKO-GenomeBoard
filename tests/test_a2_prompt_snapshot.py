"""A2 prompt structure snapshot — regression guard for the narrate-only rewrite.

We don't snapshot the literal bytes (that would flake on any copy-edit) but
assert the structural invariants qa-engineer's review requires:

1. Cancer chair prompt injects a ``CURATED EVIDENCE`` block
2. Cancer chair JSON schema REQUIRES both ``curated_id`` and ``variant_key``
   on every ``treatment_options`` row
3. Clinical Evidence Analyst prompt carries the 금지사항 block mentioning
   ``curated_id`` and the forbidden verb ("invent"/"발명")
4. Forbidden block preserves the 반드시 한국어로 응답하세요 constraint
"""
from __future__ import annotations

from unittest.mock import MagicMock


def test_cancer_chair_prompt_has_curated_evidence_section():
    from scripts.clinical_board.agents.board_chair import BoardChair

    chair = BoardChair(client=MagicMock(), model="stub", language="ko")

    class _R:
        pass

    r = _R()
    r.drug = "Sotorasib"
    r.curated_id = "cid-sot"
    r.variant_key = "12:25:C:T"
    r.target = ""
    r.evidence_level = "A"
    r.source = "oncokb"
    r.pmids = ["32955176"]
    r.disease_context = "NSCLC"
    r.significance = "sensitivity"
    r.therapy_ids = ""
    r.raw_row = {}

    prompt = chair._build_prompt(
        "case briefing",
        [],
        curated_treatments={"12:25:C:T": [r]},
    )
    up = prompt.upper()
    assert "CURATED EVIDENCE" in up
    assert "cid-sot" in prompt
    assert "Sotorasib" in prompt


def test_cancer_chair_system_prompt_schema_requires_curated_id_and_variant_key():
    from scripts.clinical_board.agents.board_chair import CANCER_SYSTEM_PROMPT

    assert "curated_id" in CANCER_SYSTEM_PROMPT
    assert "variant_key" in CANCER_SYSTEM_PROMPT
    assert "treatment_options" in CANCER_SYSTEM_PROMPT


def test_clinical_evidence_forbidden_block_present():
    from scripts.clinical_board.agents.clinical_evidence import ClinicalEvidenceAnalyst

    agent = ClinicalEvidenceAnalyst(client=MagicMock(), model="stub", language="ko")
    text = agent.system_prompt
    assert "금지사항" in text
    assert "curated_id" in text
    assert "한국어" in text  # preserved language constraint
