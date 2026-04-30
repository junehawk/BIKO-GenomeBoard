"""Therapeutic Target Analyst prompt — v2.6 F2 differentiation.

Empirically (n=5 reproducibility, 2026-04-30) the prior Therapeutic Target
prompt produced an empty ``{}`` response on every run of the demo Cancer
VCF. The clinical-advisor consult attributed this to medgemma:27b seeing
the curator's 88-row treatment list in the case briefing and concluding
there was nothing to add.

The v2.6 prompt rewrites the agent's role: it explicitly forbids restating
curated drugs and orients the agent toward sequencing, resistance,
off-label / investigational evidence, druggability gaps, and
sensitivity / resistance co-occurrence. It also installs a
"never emit an empty response" floor.

These tests pin the new prompt's contract — they do not run the LLM.
"""

from __future__ import annotations

from scripts.clinical_board.agents.therapeutic_target import TherapeuticTargetAnalyst


def _en():
    agent = TherapeuticTargetAnalyst.__new__(TherapeuticTargetAnalyst)
    agent.language = "en"
    return agent.system_prompt


def _ko():
    agent = TherapeuticTargetAnalyst.__new__(TherapeuticTargetAnalyst)
    agent.language = "ko"
    return agent.system_prompt


def test_en_prompt_acknowledges_curator_input():
    """The English prompt must tell the agent the curator's list is in the briefing."""
    p = _en()
    assert "deterministic curator" in p.lower() or "curator" in p.lower()
    assert "case briefing" in p.lower() or "in the case briefing" in p.lower()


def test_en_prompt_forbids_restating_curated_drugs():
    p = _en()
    assert "do not restate" in p.lower()


def test_en_prompt_anti_empty_floor():
    p = _en()
    # "never emit an empty response" or equivalent strong phrase
    assert "never emit an empty response" in p.lower() or "do not emit an empty" in p.lower()


def test_en_prompt_lists_added_value_categories():
    p = _en().lower()
    # Sequencing, resistance, off-label/investigational, druggability gaps,
    # co-occurrence — at least four of these five must be named.
    categories = [
        "sequencing",
        "resistance",
        "off-label",
        "investigational",
        "druggability",
        "co-occurrence",
    ]
    hits = sum(1 for c in categories if c in p)
    assert hits >= 4, f"only {hits} added-value categories found in EN prompt: {categories}"


def test_ko_prompt_acknowledges_curator_input():
    p = _ko()
    assert "큐레이터" in p


def test_ko_prompt_forbids_restating_curated_drugs():
    p = _ko()
    assert "다시 나열" in p or "다시 읊" in p


def test_ko_prompt_anti_empty_floor():
    p = _ko()
    assert "빈 응답을 절대 emit하지" in p or "빈 응답을 절대" in p


def test_ko_prompt_lists_added_value_categories():
    p = _ko()
    categories = ["순서", "저항", "off-label", "investigational", "Druggability", "동반 변이"]
    hits = sum(1 for c in categories if c in p)
    assert hits >= 4, f"only {hits} added-value categories found in KO prompt"


def test_prompt_does_not_mention_v2_5_5_chair_keys():
    """The agent prompt must not leak Chair-side key contracts."""
    p = _en() + _ko()
    # variant_key / curated_id are Chair-only concerns
    assert "variant_key" not in p
    assert "curated_id" not in p
