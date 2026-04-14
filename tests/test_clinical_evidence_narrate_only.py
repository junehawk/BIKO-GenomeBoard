"""Clinical Evidence Analyst — narrate-only prompt rewrite (A2).

The domain agent prompt is rewritten so the LLM is instructed that curated
evidence is authoritative and it MUST NOT invent drug/target pairings. The
scrubber runs as the downstream gate, but the prompt-level instruction
prevents the common-case failure where the model confidently asserts a
plausible-but-wrong drug association.
"""

from __future__ import annotations

import dataclasses
import json
from unittest.mock import MagicMock


def _stub_row(drug, vk, cid):
    class _R:
        pass

    r = _R()
    r.drug = drug
    r.curated_id = cid
    r.variant_key = vk
    r.target = ""
    r.evidence_level = "A"
    r.source = "oncokb"
    r.pmids = ["32955176"]
    r.disease_context = "NSCLC"
    r.significance = "sensitivity"
    r.therapy_ids = ""
    r.raw_row = {}
    return r


def test_clinical_evidence_prompt_mentions_curated_evidence():
    from scripts.clinical_board.agents.clinical_evidence import ClinicalEvidenceAnalyst

    agent = ClinicalEvidenceAnalyst(client=MagicMock(), model="stub", language="ko")
    text = agent.system_prompt
    assert "CURATED EVIDENCE" in text.upper() or "curated" in text.lower()


def test_clinical_evidence_prompt_contains_forbidden_keywords():
    from scripts.clinical_board.agents.clinical_evidence import ClinicalEvidenceAnalyst

    agent = ClinicalEvidenceAnalyst(client=MagicMock(), model="stub", language="ko")
    text = agent.system_prompt
    assert "금지사항" in text
    assert "curated_id" in text
    assert "발명" in text or "invent" in text.lower()


def test_clinical_evidence_retains_korean_constraint():
    from scripts.clinical_board.agents.clinical_evidence import ClinicalEvidenceAnalyst

    agent = ClinicalEvidenceAnalyst(client=MagicMock(), model="stub", language="ko")
    assert "한국어" in agent.system_prompt


def test_no_futibatinib_for_kras_g12d_end_to_end(monkeypatch):
    """End-to-end scrub: stub LLM emits 'futibatinib' in a fabricated treatment
    row; after narrative_scrubber runs, 'futibatinib' must not appear anywhere
    in ``json.dumps(dataclasses.asdict(opinion))``.
    """
    from scripts.clinical_board import runner as runner_mod

    variant_key = "12:25398284:C:T"
    curated = {variant_key: [_stub_row("Sotorasib", variant_key, "cid-sot")]}
    monkeypatch.setattr(
        runner_mod,
        "curate_treatments",
        lambda variants, **kw: curated,
        raising=False,
    )

    fake_chair = {
        "therapeutic_headline": "KRAS G12D headline",
        "therapeutic_implications": "limited direct options",
        "therapeutic_evidence": "curated rows only",
        "treatment_options": [
            {
                "drug": "Sotorasib",
                "curated_id": "cid-sot",
                "variant_key": variant_key,
                "evidence_level": "A",
                "resistance_notes": "",
            },
            {
                "drug": "Futibatinib",
                "curated_id": "bogus",
                "variant_key": variant_key,
                "evidence_level": "A",
                "resistance_notes": "",
            },
        ],
        "actionable_findings": [],
        "clinical_actions": [],
        "immunotherapy_eligibility": "",
        "confidence": "moderate",
        "agent_consensus": "majority",
        "dissenting_opinions": [],
        "monitoring_plan": [],
    }
    fake_client = MagicMock()
    fake_client.generate_json.return_value = fake_chair
    fake_client.is_available.return_value = True
    fake_client.has_model.return_value = True
    monkeypatch.setattr(runner_mod, "OllamaClient", lambda *a, **kw: fake_client, raising=False)
    monkeypatch.setattr(runner_mod, "_load_agents", lambda *a, **kw: [], raising=False)

    report_data = {
        "sample_id": "CE-KRAS",
        "variants": [
            {
                "gene": "KRAS",
                "hgvsp": "p.Gly12Asp",
                "chrom": "12",
                "pos": 25398284,
                "ref": "C",
                "alt": "T",
                "variant": variant_key,
                "classification": "Likely Pathogenic",
                "tier": "Tier I",
                "consequence": "missense_variant",
                "variant_type": "SNV",
            },
        ],
        "summary": {"total": 1},
    }

    opinion = runner_mod.run_clinical_board(report_data, mode="cancer")
    assert opinion is not None
    payload = json.dumps(dataclasses.asdict(opinion), ensure_ascii=False).lower()
    assert "futibatinib" not in payload
    assert "sotorasib" in payload
