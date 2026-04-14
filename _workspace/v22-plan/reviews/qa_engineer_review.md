# QA-Engineer Review — v2.2 Plan

**Reviewer:** qa-engineer
**Date:** 2026-04-14
**Plan reviewed:** `docs/superpowers/plans/2026-04-14-ai-board-v2.2.md` (600 lines, draft)
**Baseline:** `pytest tests/ -q` → **728 passed, 0 failed, 41 warnings, 50.10 s** (run 2026-04-14)
**Reference audit:** `_workspace/board-polish/foundation_one_comparison.md`, `_workspace/variant-selector/00_clinical_review.md`, `_workspace/variant-selector/qa_report.md`, `_workspace/ai-board-v2/final_report.md`

---

## Verdict

**APPROVED with 8 required test additions + 1 baseline correction.**

The plan is well-structured and the TDD steps are concrete. Phase A correctly
encodes the futibatinib/KRAS and TP53 R249M lessons from the F1 audit. However,
several test fixtures stop one step short of catching the failure modes that
matter most for a clinical pilot:

- the curated_id post-processor is only tested for "non-existent ID", not for
  "valid ID that belongs to a different variant" (the most likely shape of an
  LLM hallucination once curated lists exist);
- the Tier III consequence gate is tested at the selector level but not at the
  render-layer path that the F1 audit identified as the actual leak site;
- there is no performance budget for the new OncoKB HTTP loop, which will fire
  per selected variant (up to 30 in cancer mode);
- the PM1 hotspot table tests cover the happy path but not the HGVS parser
  edge cases (`p.R249M`, `p.Arg249Met`, missing `hgvsp`) that real VEP output
  exhibits.

None of the gaps are blockers — they are bounded test additions, all S effort.
Concrete fixtures are listed in §6.

---

## Baseline correction (cover-note item)

The plan's task-5 brief states "current pytest pass count (~722 post
board-polish)". The actual baseline on `main @ e5de9f7` is **728 passed**.
The +6 delta is the cancer-board headline split + render tests
(`test_render_cancer_headline_and_body`, `test_render_cancer_fallback_no_headline`,
`test_render_cancer_evidence_below`, `test_cancer_board_opinion_has_optional_headline`,
`test_board_chair_cancer_emits_headline`, `test_board_chair_cancer_backward_compat_no_headline`).
**Update the cover note before kicking off implementation** so the
green-bar gate is calibrated to the right number; otherwise pipeline-dev will
either silently pass at <728 or fail to recognise an actual regression.

---

## Section-by-section critique

### Phase A — Patient safety

#### A1 — Curated treatments module

- TDD steps are concrete and the empty-hit / offline-mode tests are correctly
  scoped. ✅
- **Gap 1 — partial-failure modes not tested.** OncoKB free-tier returns three
  failure shapes the curator will see in production:
  1. HTTP 200 with an empty `[]` body (variant not annotated).
  2. HTTP 429 / 5xx (rate limit or transient outage).
  3. Connection refused / DNS failure (offline appliance).
  The plan only tests #1 implicitly via the empty list. **Add explicit tests
  for #2 and #3** — the curator must degrade to "CIViC-only" rather than
  raising into the agent loop.
- **Gap 2 — drug-name normaliser fragility.** Plan says "lowercase + strip
  brand suffixes". Real CIViC has rows like `"Trastuzumab (Herceptin)"`,
  `"5-Fluorouracil"`, `"trametinib + dabrafenib"` (combination therapy). Add
  a parameterised test with at least 6 alias pairs (e.g. `vemurafenib` vs
  `Zelboraf`, `5-FU` vs `5-Fluorouracil`).
- **Gap 3 — performance budget.** A cancer-mode case with 30 selected variants
  triggers up to 30 HTTP round-trips. With OncoKB free-tier rate limits
  (~3 req/s un-keyed), worst case is **10 s of OncoKB-only latency** added
  to every cancer board run. The current pipeline already takes 358 s
  (per `qa_report.md` cancer smoke). Define a hard budget — e.g.
  `curate_treatments(...)` total wall time ≤ 30 s for ≤ 30 variants when
  network-degraded — and add a test that times the loop with a stubbed
  100 ms-per-call OncoKB. **No budget = no regression detection.**

#### A2 — Narrate-only Clinical Evidence Analyst + Cancer Board Chair

- **Gap 4 — curated_id binding is not validated.** Plan step 4 says
  *"drops any row whose curated_id is not in report_data['_curated_treatments']"*.
  This catches a fabricated ID but **does not catch the more likely
  hallucination shape**: a valid `curated_id` from a *different* variant
  copy-pasted into the wrong row (e.g. an EGFR L858R curated row appearing
  under a TP53 R249M heading). The post-processor must validate
  `(curated_id, variant_key)` pairs, not just `curated_id` set membership.
  See fixture `test_adversarial_curated_id_cross_variant_binding` in §6.
- **Gap 5 — escape-hatch fallback is documented but not tested.** Plan says
  *"if the LLM refuses to emit curated_ids, fall back to a template-renderer
  Chair"*. This is the fail-safe the pilot will lean on if MedGemma is flaky.
  **Add a test** that simulates an LLM response with zero `curated_id` fields
  and asserts the template-renderer path produces a non-empty Treatment
  Options table (drugs sourced directly from `_curated_treatments`). Without
  this test, the escape hatch is dead code.
- **Gap 6 — prompt-rewrite snapshot regression.** Plan acknowledges "prompt
  rewrites can destabilise existing cancer-mode snapshots" (Risk H). The
  mitigation says "run the existing v2 showcase fixture" — that's a manual
  smoke, not a regression gate. **Add a test** that asserts the existing
  `tests/test_cancer_agents.py::test_board_chair_cancer_mode` output schema
  is unchanged (specifically: `treatment_options[0].drug` still extracts to
  `"Erlotinib"` after the prompt rewrite).

#### A3 — PM1 protein-domain hotspot table

- TDD list correctly tests TP53 R249M positive + TP53 R100Q negative + KRAS
  G12D positive. ✅
- **Gap 7 — HGVS parser edge cases.** `_extract_protein_position` (called
  from `evidence_collector.py:204` for PM5) takes `variant.hgvsp` as input.
  Real VEP `HGVSp` strings vary:
  - 3-letter: `ENSP00000269305.4:p.Arg249Met`
  - 1-letter: `p.R249M`
  - synonymous Met-init: `p.Met1?` (should not fire)
  - missing entirely: `""` (should not fire, must not raise)
  - frameshift suffix: `p.Arg249fs` (not missense — must not fire)
  **Add a parameterised test** that runs all five variants through the
  PM1 lookup and asserts the expected fire/no-fire result. Without this,
  TP53 R249M may still be missed in production if the local VEP run emits
  3-letter HGVS while the parser expects 1-letter (or vice versa).
- **Gap 8 — JSON schema validation.** Plan says the table includes a
  `$schema` and `version` field. Add a unit test that loads
  `data/pm1_hotspot_domains.json`, validates the structure against a small
  Pydantic / dataclass schema, and asserts every entry has a non-empty
  `source` containing a PMID substring (`r"PMID\s+\d+"`). This is a
  one-line fixture that prevents future careless edits.
- The choice to make `KRAS 117` `supporting` rather than `moderate` is
  defensible against Garrett 2021 PMID 33280026 — keep it.

#### A4 — ClinVar conflict reconciliation

- Four required test cases (R249M LP, R100Q VUS, BRCA1 Benign, KRAS P) are
  the right scope. ✅
- **Add fifth case (recommended, not blocker):** TP53 R175H + ClinVar
  **Pathogenic** + PM1 + PM5 → engine emits Pathogenic via ClinVar (no
  override fires, `clinvar_override_reason == ""`). This guards against
  the override accidentally firing on already-pathogenic ClinVar entries
  and re-classifying them.
- **Gap 9 — override audit log.** The plan says the override is *recorded
  and displayed in the report*. Add an explicit assertion that
  `clinvar_override_reason` propagates through `_build_variant_records`
  (in orchestrate.py) into the rendered HTML, not just into the
  `ClassificationResult` dataclass. Without this the override is invisible
  to the reviewing clinician — defeats the whole "narrow auditable override"
  framing.

### Phase B — Clinical credibility

#### B1 — Selector consequence gate

- **Critical observation.** The current `_cancer_must_reason` (verified at
  `scripts/clinical_board/variant_selector.py:251-270`) does **NOT** apply
  any consequence filter — Tier I/II/III variants pass purely on
  `(tier, classification, gene)`. The current `_cancer_may_reason`
  (`:273-301`) applies partial gating only. The rare-disease selector
  applies the full gate (line 369). **The plan's wording "already in the
  set; just enforce it" is misleading — it is in the set as a constant
  but not invoked in cancer mode at all.** Update step 2's wording so
  pipeline-dev does not assume this is a no-op.
- **Gap 10 — render-layer leak is the actual bug site.** The F1 audit
  found intronic NOTCH2/ERBB4 in the rendered Tier III table. Plan step 3
  says *"Audit render.py for any path that bypasses the selector"*. This
  is correct, but **the test fixture in step 1 only verifies the selector
  output, not the rendered HTML**. Add a fixture that:
  1. Builds a `report_data` with one intronic NOTCH2 (Tier III) and one
     missense KRAS (Tier III).
  2. Calls `render_full_report(report_data)` (or whichever function builds
     the Tier III VUS table).
  3. Asserts the rendered HTML contains `KRAS` but NOT `NOTCH2`.
  This is the only way to actually prove the F1 leak is closed.
- Splice-region / SpliceAI ≥ 0.2 rescue is a good design decision. Add one
  test for SpliceAI = 0.21 admitted, SpliceAI = 0.19 rejected.

#### B2 — MMR/Lynch-panel VUS carve-out

- Acceptance criteria are clean.
- **Gap 11 — `_REASON_PRIORITY` shuffle regression risk.** Inserting
  `VUS_MMR_Lynch` at priority 5.5 (or renumbering existing values) will
  affect the sort order in `combined.sort(key=...)` at line 246. The
  existing test `test_ordering_follows_amp_table2` (variant_selector.py
  test, line 515) asserts a specific order that does **not** include
  MMR. **Add a new ordering test** with a Lynch-panel VUS in the mix to
  pin the new priority; rerun the existing ordering test to confirm no
  drift.
- The existing `test_metadata_fields_present` checks
  `by_selection_reason` — verify the new `VUS_MMR_Lynch` reason key
  appears in the metadata after the carve-out.

#### B3 — PMID references

- **Gap 12 — PMID-format parser fuzz.** Plan acknowledges
  `"PMID: 30224644"`, `"PMID 30224644"`, `"pmid:30224644"`, bare numbers.
  Add a parameterised test with at least these forms plus the malformed
  cases the LLM will actually emit:
  - `"Pubmed: 30224644"` (wrong prefix)
  - `"PMID30224644"` (no separator)
  - `"PMID: 3022 4644"` (whitespace in number)
  - `"30224644."` (trailing punctuation)
  - `"30224644, 30224645"` (multiple in one string)
  Test both the normaliser (rejects malformed) and the renderer (handles
  multi-PMID strings without crashing).
- The "unverified source" footnote for unmatched PMIDs is the right call.
  Add a test that the unverified footnote uses a distinct CSS class so
  the report-dev review can verify it's visually distinguishable.

#### B4 — Patient demographics header

- CLI argument list is complete. ✅
- **Gap 13 — PHI persistence test.** Plan says
  `reporting.persist_patient_metadata: false` gates downstream persistence.
  **Add an explicit test**: invoke `run_pipeline(..., patient={...})` with
  the flag false (default), verify that `data/knowledge_base/kb.sqlite3`
  contains zero patient name strings after the run. Without this, a future
  refactor could quietly start persisting PHI without any test catching it.
- DOB validation: plan says "DOB parsable, sex in {M,F,U,Other}". Spell
  out the parser format(s) accepted (ISO-8601 only? Korean YYYY.MM.DD?)
  and add a negative test for the unaccepted format.

### Phase C — Data / ops

- C1 / C2 are documentation-only, no test surface. ✅
- Recommend C2 also documents the new pytest baseline (728 → expected
  ~745 after v2.2 lands, given the new test files in §6).

---

## Punch list (concrete, ordered by clinical impact)

1. **[A1/perf]** Define and test a wall-time budget for `curate_treatments()` —
   ≤ 30 s for ≤ 30 variants under degraded network. Pipeline currently runs
   358 s for cancer; another uncapped HTTP loop is unacceptable for on-prem.
2. **[A2/binding]** Validate `(curated_id, variant_key)` pair, not just
   `curated_id` set membership. Add `test_adversarial_curated_id_cross_variant_binding`
   (fixture in §6).
3. **[B1/render]** Add a render-layer test for the Tier III leak. Selector-only
   test does not prove the F1 bug is closed.
4. **[A1/network]** Add explicit tests for OncoKB HTTP 429 and connection-refused
   degradation paths. The curator must degrade to CIViC-only, never raise.
5. **[A3/HGVS]** Parameterise PM1 lookup test across 1-letter / 3-letter HGVS
   and missing-hgvsp cases. TP53 R249M may still be missed in production
   without this.
6. **[A2/escape]** Add a test for the template-renderer fallback when LLM
   emits no `curated_id` fields. The fallback is currently undocumented dead
   code in the plan.
7. **[B4/PHI]** Add a regression test that asserts patient name strings
   never reach the KB SQLite when `persist_patient_metadata=false`.
8. **[B2/priority]** Add a new ordering test with `VUS_MMR_Lynch` in the
   mix; update the existing `test_ordering_follows_amp_table2` if it drifts.
9. **[meta/baseline]** Update the task-5 cover note: pytest baseline is **728**,
   not 722. Set the post-v2.2 expected count (probably ~745) before
   pipeline-dev starts.

---

## Adversarial test fixtures (concrete, paste-ready)

### Fixture 1 — futibatinib hallucination (A2 regression)

`tests/test_clinical_evidence_narrate_only.py`:

```python
"""Adversarial test: even if MedGemma hallucinates futibatinib for KRAS G12D,
the curated_id post-processor MUST drop the row. This regression fixture
encodes the exact failure mode found in the FoundationOne CDx audit
(_workspace/board-polish/foundation_one_comparison.md §3 "KRAS inhibitor
family errors").
"""
from unittest.mock import MagicMock
import pytest


def test_adversarial_kras_futibatinib_hallucination_blocked(monkeypatch):
    from scripts.clinical_board import runner as runner_mod
    from scripts.clinical_board.curated_treatments import CuratedTreatment
    from scripts.clinical_board.models import CancerBoardOpinion

    sotorasib = CuratedTreatment(
        curated_id="abc123sotor",
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
        curated_id="def456adag",
        drug="Adagrasib",
        target="KRAS",
        evidence_level="A",
        source="oncokb",
        pmids=["35658005"],
        disease_context="NSCLC, KRAS G12C",
        significance="sensitivity",
        raw_row={},
    )
    monkeypatch.setattr(
        runner_mod,
        "curate_treatments",
        lambda variants, **kw: {"12:25398284:C:T": [sotorasib, adagrasib]},
    )

    # MedGemma simulated to hallucinate futibatinib AND fabricate a curated_id.
    fake_chair_response = {
        "therapeutic_headline": "Stage IV PDAC — KRAS G12D",
        "therapeutic_implications": "...",
        "therapeutic_evidence": "Curated rows below.",
        "treatment_options": [
            {"drug": "Sotorasib", "curated_id": "abc123sotor",
             "evidence_level": "A", "resistance_notes": ""},
            {"drug": "Futibatinib", "curated_id": "ghi789fabricated",
             "evidence_level": "A", "resistance_notes": ""},  # hallucination
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
    fake_client.generate_json.return_value = fake_chair_response
    fake_client.is_available.return_value = True
    fake_client.has_model.return_value = True

    monkeypatch.setattr(runner_mod, "OllamaClient", lambda *a, **kw: fake_client)
    monkeypatch.setattr(runner_mod, "_load_agents", lambda *a, **kw: [])

    report_data = {
        "sample_id": "ADV-FUTI-001",
        "variants": [{
            "gene": "KRAS", "hgvsp": "p.Gly12Asp",
            "chrom": "12", "pos": 25398284, "ref": "C", "alt": "T",
            "classification": "Likely Pathogenic", "tier": "Tier I",
            "consequence": "missense_variant",
            "cancer_gene_type": "Oncogene", "oncokb_level": "1",
            "hpo_score": 0, "gnomad_af": None, "variant_type": "SNV",
        }],
        "summary": {"total": 1},
    }

    opinion = runner_mod.run_clinical_board(report_data, mode="cancer")

    drugs = [r["drug"].lower() for r in opinion.treatment_options]
    assert "futibatinib" not in drugs, (
        "Curated_id post-processor failed: futibatinib reached the rendered "
        "Treatment Options. This is the exact F1-audit regression."
    )
    assert all(
        r["curated_id"] in {"abc123sotor", "def456adag"}
        for r in opinion.treatment_options
    ), "Fabricated curated_id was not stripped"
```

### Fixture 2 — cross-variant curated_id binding (A2 binding)

```python
def test_adversarial_curated_id_cross_variant_binding(monkeypatch):
    """The post-processor must reject a curated_id that exists for a DIFFERENT
    variant than the one being discussed. Without (curated_id, variant_key)
    pair validation, a real-world LLM can copy-paste an EGFR row into a TP53
    section and pass the existing 'set membership' check.
    """
    from scripts.clinical_board import runner as runner_mod
    from scripts.clinical_board.curated_treatments import CuratedTreatment
    from unittest.mock import MagicMock

    osimertinib = CuratedTreatment(
        curated_id="egfr_osi_001",
        drug="Osimertinib", target="EGFR", evidence_level="A",
        source="oncokb", pmids=["28522753"],
        disease_context="NSCLC EGFR L858R",
        significance="sensitivity", raw_row={},
    )
    monkeypatch.setattr(
        runner_mod, "curate_treatments",
        lambda variants, **kw: {
            "7:55259515:T:G": [osimertinib],     # EGFR L858R
            "17:7674220:G:A": [],                 # TP53 R249M — empty
        },
    )

    # LLM sneaks the EGFR row into the TP53 Treatment Options section.
    fake_response = {
        "therapeutic_headline": "TP53 R249M",
        "therapeutic_implications": "...",
        "therapeutic_evidence": "",
        "treatment_options": [
            # Wrong: this row's curated_id binds to the EGFR variant, not TP53.
            {"drug": "Osimertinib", "curated_id": "egfr_osi_001",
             "variant_key": "17:7674220:G:A",  # claims to be TP53!
             "evidence_level": "A", "resistance_notes": ""},
        ],
        "actionable_findings": [], "clinical_actions": [],
        "immunotherapy_eligibility": "", "confidence": "moderate",
        "agent_consensus": "majority", "dissenting_opinions": [],
        "monitoring_plan": [],
    }
    fake_client = MagicMock()
    fake_client.generate_json.return_value = fake_response
    fake_client.is_available.return_value = True
    fake_client.has_model.return_value = True
    monkeypatch.setattr(runner_mod, "OllamaClient", lambda *a, **kw: fake_client)
    monkeypatch.setattr(runner_mod, "_load_agents", lambda *a, **kw: [])

    report_data = {
        "sample_id": "ADV-BIND-001",
        "variants": [
            {"gene": "EGFR", "hgvsp": "p.Leu858Arg", "chrom": "7",
             "pos": 55259515, "ref": "T", "alt": "G",
             "classification": "Likely Pathogenic", "tier": "Tier I",
             "consequence": "missense_variant", "variant_type": "SNV"},
            {"gene": "TP53", "hgvsp": "p.Arg249Met", "chrom": "17",
             "pos": 7674220, "ref": "G", "alt": "A",
             "classification": "Likely Pathogenic", "tier": "Tier II",
             "consequence": "missense_variant", "variant_type": "SNV"},
        ],
        "summary": {"total": 2},
    }

    opinion = runner_mod.run_clinical_board(report_data, mode="cancer")

    drugs = [r["drug"].lower() for r in opinion.treatment_options]
    assert "osimertinib" not in drugs, (
        "Cross-variant curated_id rebinding was not caught — the post-processor "
        "must validate (curated_id, variant_key) PAIRS, not bare curated_id."
    )
```

### Fixture 3 — TP53 R249M ACMG PM1 + override regression

`tests/test_pm1_hotspot_domains.py`:

```python
import pytest
from scripts.classification.evidence_collector import collect_additional_evidence
from scripts.common.models import Variant


@pytest.mark.parametrize("hgvsp,expect_pm1", [
    ("ENSP00000269305.4:p.Arg249Met", True),
    ("p.Arg249Met",                   True),
    ("p.R249M",                       True),
    ("p.Arg100Gln",                   False),  # not in DBD hotspot range
    ("p.Met1?",                       False),  # synonymous Met-init
    ("p.Arg249fs",                    False),  # frameshift, not missense
    ("",                              False),  # empty hgvsp must not raise
])
def test_pm1_tp53_dbd_hotspot_range(hgvsp, expect_pm1):
    v = Variant(
        chrom="17", pos=7674220, ref="G", alt="A",
        gene="TP53", hgvsp=hgvsp,
        consequence="missense_variant" if "fs" not in hgvsp and "?" not in hgvsp else "frameshift_variant",
        domains="",  # empty — table must rescue
    )
    codes = collect_additional_evidence(v)
    if expect_pm1:
        assert "PM1" in codes, f"{hgvsp} should fire PM1 via DBD table"
    else:
        assert "PM1" not in codes, f"{hgvsp} should NOT fire PM1"
```

### Fixture 4 — Render-layer Tier III consequence gate (B1)

`tests/test_variant_selector_v22.py`:

```python
def test_render_tier3_excludes_intronic_after_v22():
    """The F1 audit found NOTCH2/ERBB4/RELN intronic variants in BIKO's Tier III
    table. The selector AND the render path must both apply the consequence gate.
    """
    from scripts.clinical_board.render import render_full_report  # or equivalent

    report_data = {
        "sample_id": "TEST-T3",
        "variants": [
            {"gene": "NOTCH2", "tier": "Tier III", "classification": "VUS",
             "consequence": "intron_variant", "hgvsp": "",
             "variant": "chr1:120498000C>T", "variant_type": "SNV"},
            {"gene": "KRAS",   "tier": "Tier III", "classification": "VUS",
             "consequence": "missense_variant", "hgvsp": "p.Gly12Asp",
             "variant": "chr12:25398284C>T", "variant_type": "SNV"},
        ],
        "summary": {"total": 2},
    }
    html = render_full_report(report_data)
    assert "KRAS" in html
    assert "NOTCH2" not in html, (
        "Render-layer Tier III leak: intronic NOTCH2 was emitted into the VUS "
        "table. v2.2 B1 step 3 (render.py audit) is incomplete."
    )
```

### Fixture 5 — PHI non-persistence (B4)

`tests/test_patient_header.py`:

```python
def test_patient_metadata_not_persisted_by_default(tmp_path, monkeypatch):
    """When persist_patient_metadata=false (default), patient name strings
    MUST NOT appear in kb.sqlite3 after a run.
    """
    import sqlite3
    from scripts.orchestrate import run_pipeline

    kb_path = tmp_path / "kb.sqlite3"
    monkeypatch.setenv("BIKO_KB_PATH", str(kb_path))

    run_pipeline(
        vcf_path="data/sample_vcf/demo_with_in_silico.vcf",
        output_path=str(tmp_path / "out.html"),
        skip_api=True, mode="cancer",
        patient={"name": "TestPatient_Sentinel_XYZ", "dob": "1970-01-01", "sex": "M"},
    )
    if not kb_path.exists():
        pytest.skip("KB not created by this run")

    conn = sqlite3.connect(kb_path)
    rows = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'"
    ).fetchall()
    for (table,) in rows:
        try:
            cur = conn.execute(f"SELECT * FROM {table}")
            for row in cur:
                assert "TestPatient_Sentinel_XYZ" not in str(row), (
                    f"PHI leak: patient name found in KB table {table}"
                )
        except sqlite3.OperationalError:
            continue
    conn.close()
```

---

## Performance gates (summary)

| Path | Current | v2.2 budget | Test |
|---|---|---|---|
| Cancer board (Ollama + curator) | 358 s (qa_report.md) | ≤ 400 s | smoke (manual) |
| `curate_treatments(30 vars)` happy | n/a | ≤ 5 s w/ network, ≤ 1 s offline | new |
| `curate_treatments(30 vars)` degraded | n/a | ≤ 30 s (timeout-bounded) | new |
| Pytest full suite | 50 s | ≤ 70 s | CI |

The 30 s degraded budget assumes per-call timeout ≤ 1 s + retry once. Define
this in `clinical_board.curated_treatments.http_timeout_s` and test the
config plumbing.

---

## Clinical acceptance criteria (alignment with plan §Patient-safety gating)

The plan correctly states *"A pilot cannot begin until Phase A is merged and
A4's override reason has been clinically reviewed by the attending oncologist
on the pilot team."* Operationalise this with:

1. A `_workspace/v22-plan/clinical_signoff.md` template containing:
   - Reviewer name, date, ACMG/AMP citation count.
   - Three explicit fixture cases the reviewer signs off on (TP53 R249M,
     KRAS G12D futibatinib block, PMS2 R563Q surfacing).
   - A "rollback trigger" line: any clinical incident at the pilot site
     reverts to v2.1 within 24 h.
2. A pre-pilot smoke that runs all five §6 fixtures plus the existing
   integration tests (`test_integration.py`, `test_orchestrate.py`,
   `test_render_rerender.py`) and produces a one-page summary report.
3. The Foundation-One comparison findings 1, 2, 3, 6, 7 from
   `_workspace/board-polish/foundation_one_comparison.md` §7 are explicitly
   re-verified against the v2.2 codegen showcase regen.

---

## Final word

The plan is solid. The TDD discipline is consistent, the file map is
accurate, the dependency graph is correct, and Phase A targets exactly the
right risks. The 8 test additions in §6 close the gaps that would otherwise
let the F1-audit regressions reappear in a different shape. Adopt them and
this plan is pilot-ready from a verification standpoint.

**Sign-off:** APPROVED with required test additions.
