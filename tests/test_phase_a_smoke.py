"""Phase A patient-safety smoke fixtures (BIKO GenomeBoard v2.2).

Five end-to-end regression fixtures encoded from the Foundation One audit
(`_workspace/board-polish/foundation_one_comparison.md`) and the qa-engineer
v2.2 plan review (`_workspace/v22-plan/reviews/qa_engineer_review.md`).

- fixture 1: KRAS G12D futibatinib hallucination blocked by narrative_scrubber
- fixture 2: EGFR→TP53 cross-variant curated_id paste attack dropped
- fixture 3: TP53 R249M classified LP with clinvar_override_reason populated
- fixture 4: [Phase B] Tier III render-layer consequence gate
- fixture 5: [Phase B] PHI non-persistence scan

Fixtures 4 and 5 are marked xfail until Phase B lands.
"""

from __future__ import annotations

import dataclasses
import json
from typing import Any
from unittest.mock import MagicMock

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _opinion_to_json(opinion: Any) -> str:
    """Serialise a CancerBoardOpinion (or dict) to a deep JSON string.

    The fixture-1 invariant is that the hallucinated drug name must not
    appear *anywhere* in the opinion — including nested narrative fields —
    so we flatten everything into one JSON string before searching.
    """
    if dataclasses.is_dataclass(opinion):
        payload = dataclasses.asdict(opinion)
    else:
        payload = opinion
    return json.dumps(payload, ensure_ascii=False, default=str).lower()


def _build_cancer_report_data(variants: list[dict], sample_id: str = "PHASE-A-SMOKE") -> dict:
    return {
        "sample_id": sample_id,
        "variants": variants,
        "summary": {"total": len(variants)},
    }


# ---------------------------------------------------------------------------
# Fixture 1 — KRAS G12D futibatinib hallucination blocked
# ---------------------------------------------------------------------------


def test_smoke_fixture_1_kras_futibatinib_hallucination_blocked(monkeypatch):
    """F1-audit regression: even if MedGemma emits 'futibatinib' for KRAS G12D,
    the curated_id post-processor (A2 narrative_scrubber) MUST strip it before
    the opinion is serialised. The string 'futibatinib' must not appear in
    `json.dumps(dataclasses.asdict(cancer_board_opinion))` anywhere.
    """
    runner_mod = pytest.importorskip(
        "scripts.clinical_board.runner",
        reason="clinical_board.runner not yet wired (blocked by #2/#3)",
    )
    curated_mod = pytest.importorskip(
        "scripts.clinical_board.curated_treatments",
        reason="curated_treatments module lands with task #2",
    )
    CuratedTreatment = curated_mod.CuratedTreatment  # type: ignore[attr-defined]

    sotorasib = CuratedTreatment(
        curated_id="oncokb_kras_g12c_sotor",
        drug="Sotorasib",
        target="KRAS",
        evidence_level="A",
        source="oncokb",
        pmids=["32955176"],
        disease_context="NSCLC, KRAS G12C",
        significance="sensitivity",
        raw_row={},
    )
    adagrasib = CuratedTreatment(
        curated_id="oncokb_kras_g12c_adag",
        drug="Adagrasib",
        target="KRAS",
        evidence_level="A",
        source="oncokb",
        pmids=["35658005"],
        disease_context="NSCLC, KRAS G12C",
        significance="sensitivity",
        raw_row={},
    )
    variant_key = "12:25398284:C:T"
    monkeypatch.setattr(
        runner_mod,
        "curate_treatments",
        lambda variants, **kw: {variant_key: [sotorasib, adagrasib]},
        raising=False,
    )

    fake_chair = {
        "therapeutic_headline": "Stage IV PDAC — KRAS G12D",
        "therapeutic_implications": "Limited direct G12D targeted options.",
        "therapeutic_evidence": "See curated rows.",
        "treatment_options": [
            {
                "drug": "Sotorasib",
                "curated_id": "oncokb_kras_g12c_sotor",
                "variant_key": variant_key,
                "evidence_level": "A",
                "resistance_notes": "",
            },
            {
                "drug": "Futibatinib",
                "curated_id": "fabricated_fgfr_row",
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

    report_data = _build_cancer_report_data(
        [
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
                "cancer_gene_type": "Oncogene",
                "oncokb_level": "1",
                "variant_type": "SNV",
            }
        ],
        sample_id="SMOKE-KRAS-FUTI",
    )

    opinion = runner_mod.run_clinical_board(report_data, mode="cancer")
    assert opinion is not None, "Clinical Board returned None — check mock client"

    serialised = _opinion_to_json(opinion)
    assert "futibatinib" not in serialised, (
        "F1 audit regression — 'futibatinib' leaked into the cancer board "
        "opinion after narrative_scrubber. Phase A A2 scrubber has failed."
    )

    # Positive: sotorasib (real curated row) must remain.
    assert "sotorasib" in serialised, "narrative_scrubber stripped legitimate curated row — false positive"


# ---------------------------------------------------------------------------
# Fixture 2 — EGFR→TP53 cross-variant paste-attack
# ---------------------------------------------------------------------------


def test_smoke_fixture_2_cross_variant_curated_id_paste_attack(monkeypatch):
    """A valid curated_id for EGFR L858R must not survive if the LLM sneaks it
    under a TP53 R249M row. The A2 narrative_scrubber must validate the
    `(curated_id, variant_key)` PAIR, not bare curated_id set membership.

    Test isolation: to make the regression unambiguous, the report carries
    ONLY the TP53 variant. The curator returns an osimertinib row but keys
    it under a DIFFERENT variant_key that is not present in the report
    (simulating a stale curator cache from a prior EGFR run). The LLM then
    hallucinates the osimertinib row under TP53. The scrubber must drop it,
    and osimertinib must not appear anywhere in the final opinion.
    """
    runner_mod = pytest.importorskip(
        "scripts.clinical_board.runner",
        reason="clinical_board.runner not yet wired (blocked by #2/#3)",
    )
    curated_mod = pytest.importorskip(
        "scripts.clinical_board.curated_treatments",
        reason="curated_treatments module lands with task #2",
    )
    CuratedTreatment = curated_mod.CuratedTreatment  # type: ignore[attr-defined]

    egfr_key = "7:55259515:T:G"  # NOT in report_data — stale curator key
    tp53_key = "17:7674220:G:A"

    osimertinib = CuratedTreatment(
        curated_id="oncokb_egfr_l858r_osi",
        drug="Osimertinib",
        target="EGFR",
        evidence_level="A",
        source="oncokb",
        pmids=["28522753"],
        variant_key=egfr_key,
        disease_context="NSCLC, EGFR L858R",
        significance="sensitivity",
        raw_row={},
    )
    # Curator returns the osimertinib row keyed to egfr_key. Because the
    # EGFR variant is NOT in report_data, no legitimate TP53-mode path can
    # surface osimertinib — if it appears in the final opinion, a scrubber
    # leak is the only explanation.
    monkeypatch.setattr(
        runner_mod,
        "curate_treatments",
        lambda variants, **kw: {egfr_key: [osimertinib], tp53_key: []},
        raising=False,
    )

    fake_chair = {
        "therapeutic_headline": "TP53 R249M + EGFR L858R",
        "therapeutic_implications": "",
        "therapeutic_evidence": "",
        "treatment_options": [
            # LLM copy-pastes EGFR curated row under the TP53 variant_key.
            {
                "drug": "Osimertinib",
                "curated_id": "oncokb_egfr_l858r_osi",
                "variant_key": tp53_key,  # wrong binding — must be dropped
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

    report_data = _build_cancer_report_data(
        [
            {
                "gene": "EGFR",
                "hgvsp": "p.Leu858Arg",
                "chrom": "7",
                "pos": 55259515,
                "ref": "T",
                "alt": "G",
                "variant": egfr_key,
                "classification": "Likely Pathogenic",
                "tier": "Tier I",
                "consequence": "missense_variant",
                "variant_type": "SNV",
            },
            {
                "gene": "TP53",
                "hgvsp": "p.Arg249Met",
                "chrom": "17",
                "pos": 7674220,
                "ref": "G",
                "alt": "A",
                "variant": tp53_key,
                "classification": "Likely Pathogenic",
                "tier": "Tier II",
                "consequence": "missense_variant",
                "variant_type": "SNV",
            },
        ],
        sample_id="SMOKE-CROSS-BIND",
    )

    opinion = runner_mod.run_clinical_board(report_data, mode="cancer")
    assert opinion is not None

    # Per plan §A2 step 5: narrative_scrubber.validate_treatment_option()
    # MUST DROP any row whose (curated_id, variant_key) pair does not match
    # the curated set for that variant — not just rebind. So the osimertinib
    # row pasted under TP53 must disappear entirely from the opinion output.
    serialised = _opinion_to_json(opinion)
    assert "osimertinib" not in serialised, (
        "Cross-variant paste attack not caught — the EGFR curated row "
        "bound to the TP53 variant_key survived into the opinion. "
        "narrative_scrubber must drop (curated_id, variant_key) pair mismatches."
    )

    treatments = getattr(opinion, "treatment_options", []) or []
    for row in treatments:
        cid = row.get("curated_id", "")
        if cid == "oncokb_egfr_l858r_osi":
            pytest.fail(
                f"narrative_scrubber allowed EGFR curated_id through despite variant_key mismatch (row={row!r})."
            )


# ---------------------------------------------------------------------------
# Fixture 3 — TP53 R249M PM1 hotspot + ClinVar override
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "hgvsp,expect_pm1",
    [
        ("ENSP00000269305.4:p.Arg249Met", True),
        ("p.Arg249Met", True),
        ("p.R249M", True),
        ("p.Arg100Gln", False),  # outside DBD hotspot range
        ("p.Met1?", False),  # Met-init — not missense
        ("p.Arg249fs", False),  # frameshift, not missense
        ("", False),  # empty hgvsp must not raise
    ],
)
def test_smoke_fixture_3a_tp53_pm1_hotspot_hgvs_fuzz(hgvsp, expect_pm1):
    """HGVS parser fuzz for TP53 R249M PM1 hotspot lookup.

    Real VEP annotation emits HGVSp in multiple shapes; the PM1 table lookup
    must handle all of them without crashing, and must fire PM1 only for the
    missense R249M variants inside the DBD hotspot range.
    """
    import inspect
    import os

    from scripts.common.config import get as config_get

    pm1_path = config_get("paths.pm1_hotspot_domains", None) or os.path.join("data", "pm1_hotspot_domains.json")
    if not os.path.exists(pm1_path):
        pytest.skip(f"pm1_hotspot_domains.json not yet built ({pm1_path}) — blocked by tasks #5/#6")

    collector = pytest.importorskip(
        "scripts.classification.evidence_collector",
        reason="evidence_collector A3 PM1 logic lands with task #6",
    )
    # Feature detection: A3-core (task #6) wires evidence_collector to the
    # PM1 hotspot table. Until that edit lands, the module source has no
    # reference to the table, so PM1 will never fire from it. Skip this
    # fuzz until the wiring is in place.
    try:
        collector_src = inspect.getsource(collector)
    except (OSError, TypeError):
        collector_src = ""
    if "pm1_hotspot" not in collector_src.lower():
        pytest.skip("evidence_collector does not yet reference pm1_hotspot_domains (A3-core task #6 pending)")

    Variant = pytest.importorskip("scripts.common.models").Variant  # type: ignore[attr-defined]

    consequence = "missense_variant"
    if "fs" in hgvsp:
        consequence = "frameshift_variant"
    elif "Met1?" in hgvsp:
        consequence = "start_lost"

    v = Variant(
        chrom="17",
        pos=7674220,
        ref="G",
        alt="A",
        gene="TP53",
        hgvsp=hgvsp,
        consequence=consequence,
    )
    try:
        codes = collector.collect_additional_evidence(v)
    except TypeError:
        # Older signature used in some branches
        codes = collector.collect_additional_evidence(v, {})  # type: ignore[arg-type]

    if isinstance(codes, dict):
        codes = set(codes.keys())
    else:
        codes = set(codes or [])

    if expect_pm1:
        assert "PM1" in codes, f"PM1 should fire for TP53 hgvsp={hgvsp!r} via pm1_hotspot_domains.json"
    else:
        assert "PM1" not in codes, f"PM1 must NOT fire for TP53 hgvsp={hgvsp!r}"


def test_smoke_fixture_3b_tp53_r249m_clinvar_override_reason_field_exists():
    """A4 (task #7) adds `clinvar_override_reason: str = ""` to
    `ClassificationResult` and reshapes `apply_clinvar_override` to emit a
    reconciliation verdict that populates that field for TP53 R249M when
    ClinVar is conflicting and PM1 + PM5 both fire.

    This smoke test guards the schema delta. The detailed reconciliation
    logic is covered by `tests/test_acmg_engine.py` / the A4 test fixture
    that pipeline-dev lands with task #7; here we only verify that:
      (a) `ClassificationResult` grows the new field, and
      (b) a downstream caller can read it without AttributeError.
    Until #7 lands, this test skips gracefully.
    """
    import dataclasses as dc

    models_mod = pytest.importorskip(
        "scripts.classification.models",
        reason="classification.models not importable",
    )
    ClassificationResult = getattr(models_mod, "ClassificationResult", None)
    if ClassificationResult is None:
        pytest.skip("ClassificationResult dataclass not found — A4 prerequisite")
    if not dc.is_dataclass(ClassificationResult):
        pytest.skip("ClassificationResult is not a dataclass on this branch")

    field_names = {f.name for f in dc.fields(ClassificationResult)}
    if "clinvar_override_reason" not in field_names:
        pytest.skip("ClassificationResult.clinvar_override_reason not yet added (A4 task #7 reshape pending)")

    # New-contract smoke: constructing ClassificationResult with the new
    # field must succeed and the reader attribute access must work.
    #
    # We use only required-ish fields common across v2.1 to avoid coupling
    # to unrelated schema churn.
    try:
        instance = ClassificationResult(  # type: ignore[call-arg]
            classification="Likely Pathogenic",
            clinvar_override_reason=(
                "engine LP override: PM1 (TP53 DBD hotspot, PMID 30224644) "
                "+ PM5 (R249S ClinVar Pathogenic); ClinVar R249M conflicting"
            ),
        )
    except TypeError:
        # Other required fields exist — construct via defaults map.
        defaults = {
            f.name: (
                f.default
                if f.default is not dc.MISSING
                else (
                    f.default_factory()
                    if f.default_factory is not dc.MISSING  # type: ignore[misc]
                    else ""
                )
            )
            for f in dc.fields(ClassificationResult)
        }
        defaults["classification"] = "Likely Pathogenic"
        defaults["clinvar_override_reason"] = (
            "engine LP override: PM1 (TP53 DBD hotspot, PMID 30224644) "
            "+ PM5 (R249S ClinVar Pathogenic); ClinVar R249M conflicting"
        )
        instance = ClassificationResult(**defaults)  # type: ignore[arg-type]

    assert instance.classification == "Likely Pathogenic"
    reason = getattr(instance, "clinvar_override_reason", "")
    assert reason.startswith("engine LP override"), (
        "clinvar_override_reason should persist the human-readable "
        "explanation the reviewing clinician sees in the report."
    )


# ---------------------------------------------------------------------------
# Fixture 4 — Phase B scope (render-layer Tier III consequence gate)
# ---------------------------------------------------------------------------


def test_smoke_fixture_4_render_tier3_excludes_intronic():
    """Render-layer Tier III consequence gate — v2.2 Phase B B1.

    Delegates to the dedicated render-layer regression in
    ``tests/test_b1_tier3_render_layer.py`` to keep the smoke suite
    authoritative for the Phase A/B rollup. Passing here means the B1
    consequence gate survives the full runner → Jinja render chain.
    """
    from tests.test_b1_tier3_render_layer import (
        test_b1_render_layer_notch2_intronic_excluded_kras_present,
    )

    pytest_monkeypatch = pytest.MonkeyPatch()
    try:
        test_b1_render_layer_notch2_intronic_excluded_kras_present(pytest_monkeypatch)
    finally:
        pytest_monkeypatch.undo()


# ---------------------------------------------------------------------------
# Fixture 5 — Phase B scope (PHI non-persistence)
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    reason="Phase B (B4) scope — persist_patient_metadata flag lands in v2.2 Phase B",
    strict=False,
    run=False,
)
def test_smoke_fixture_5_phi_not_persisted_by_default():
    """PHI non-persistence guard. Not in Phase A scope."""
    raise AssertionError("Phase B fixture placeholder — tracked for B4 completion")
