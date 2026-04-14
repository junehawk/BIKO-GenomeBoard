# v2.2 Plan — Draft Cover Note

**Author:** clinical-advisor
**Date:** 2026-04-14
**Plan:** `docs/superpowers/plans/2026-04-14-ai-board-v2.2.md`
**Status:** DRAFT — pending multi-agent review (pipeline-dev, db-dev, report-dev, qa-engineer)

## What I chose to include and why

I framed this plan around the four "blocker" findings from `_workspace/board-polish/foundation_one_comparison.md`, ranked by patient-safety impact:

1. **Curate-then-narrate (A1 + A2)** — the futibatinib/KRAS hallucination in the codegen showcase is the single most dangerous artifact in the current AI Board. A clinician reviewing that row and trusting it could refer a patient to an irrelevant drug. The curator module (OncoKB free API + local CIViC, merged) is the structural fix: the LLM cannot hallucinate drugs it was never shown. I gave this two tasks because the infrastructure (A1) and the prompt rewrite (A2) are different kinds of work, and A2's risk profile (prompt destabilisation) deserves its own commit and review.

2. **TP53 R249M LP classification (A3 + A4)** — this is the classic "ACMG engine vs ClinVar conflict" problem. The fix has two parts: (a) a PM1 hotspot table keyed by codon-range (not VEP `DOMAINS`), (b) a narrow ClinVar-conflict override that fires only when engine evidence is strong enough. I deliberately designed A4 as a *narrow* override with an explicit reason field so that every override is auditable in the rendered report — this was a clinical-judgment call, not a purely technical one.

3. **Selector polish (B1 + B2)** — these are easy wins but I kept them in Phase B, not Phase A, because they do not directly endanger a patient. Mis-ranked VUS are a quality issue; mis-cited drugs and mis-classified hotspots are a safety issue.

4. **PMID provenance (B3)** — A1 and A3 naturally surface PMIDs, so B3 just extends the plumbing. I put it in Phase B because it is credibility polish, not safety — a rendered report without citations is not *dangerous*, it is *unacceptable for pilot*.

5. **Patient demographics header (B4)** — pure regulatory baseline. Easy implementation, gated on `reporting.persist_patient_metadata: false` by default so that PHI does not leak into the KB.

## What I deferred and why

**Clinical trials matching (ClinicalTrials.gov + CRIS).** The team-lead task description explicitly excluded this for on-prem Korean operation. I agree with the exclusion for v2.2 but flag it as the highest-value v2.3 item — stage IV patients depend on trial enrollment.

**MSI / HRD / fusions / mutational signatures.** All gated on data sources or reference DBs we do not yet have. Belongs in a separate v2.3+ roadmap plan.

**Tumor-type-scoped narrative** ("in pancreatic cancer, G12D is observed in 91-95%"). F1 does this well; BIKO's current narrative is tumor-agnostic. This is high-value readability but not patient-safety — deferred to v2.3 to avoid scope creep.

**OncoKB commercial license.** The free-tier API is sufficient for variant-level queries. If a clinical partner later demands the commercial API (for full therapy-level actionability), it is a drop-in config change to `oncokb_api.py`.

## Where I applied clinical judgment

- **A4 ClinVar override narrowness.** A blanket "engine wins over ClinVar" rule would be dangerous. The override I specified fires only when PM1 (from the hotspot table, with a PMID source) AND PM5 (from an adjacent ClinVar P/LP entry) both fire AND ClinVar is specifically `"Conflicting"` (not Pathogenic or Benign). This is the clinically defensible window — we are not disagreeing with ClinVar's *consensus*, we are resolving a case where ClinVar has no consensus.
- **PM1_Supporting strength.** I recommended `PM1_Supporting` (not full `PM1`) for less-established hotspots like KRAS 117. This avoids over-weighting marginal residues while still flowing into a PP3+PM2_Supporting combination that can reach LP when other evidence is present.
- **MMR carve-out scope.** I restricted it to the canonical Lynch panel (MLH1/MSH2/MSH6/PMS2/EPCAM) rather than a broader DNA-repair gene list (ATM/BRCA/PALB2 etc.) because the Lynch panel is the clinically actionable screen (pembrolizumab tumour-agnostic, reflex germline counseling). Broader DDR-gene VUS admission is a separate question for v2.3.
- **No LLM in the Treatment Options loop, as a fallback.** I flagged in A2's risk section that if the narrate-only prompt destabilises too aggressively, the escape hatch is a *template-renderer* Chair — zero LLM for treatment_options, pure curated-row rendering. This is not my first-choice design (narrative quality suffers), but I wanted the plan to name the fallback so implementation is not blocked if the prompt rewrite proves fragile.
- **Phase A gating.** I made explicit that a pilot cannot begin until A4's override reason has been reviewed by an attending oncologist on the pilot team. This is a human-in-the-loop gate, not a code gate, and it needs to be on the record.

## Open questions for reviewers

### For **pipeline-dev**

1. Does the curator module live correctly under `scripts/clinical_board/` (as "Board infrastructure") or under `scripts/somatic/` (as "somatic evidence layer")? I chose clinical_board/ because the consumer is the Board; push back if this conflicts with module boundary convention.
2. Is there a cleaner place to stash `_curated_treatments` than on `report_data`? The `_` prefix matches `_board_variants` / `_board_selection_metadata`, but if the convention is evolving, please flag.
3. For A4's ClinVar override: should the override be *computed inline* in `classify()`, or should we add a post-classification reconciliation pass so the override logic is isolated and auditable? I leaned toward inline for simplicity; the post-pass option is cleaner but more plumbing.

### For **db-dev**

1. Does the PM1 hotspot JSON table (`data/pm1_hotspot_domains.json`) need DB-version tracking alongside ClinVar/gnomAD version manifests? I drafted it as a static JSON with a top-level `version` field, but if the project convention is to track all curated tables via `version_manager.py`, please adjust.
2. The OncoKB free-tier API requires an internet-optional mode for on-prem deployments. I specified a local JSON cache at `data/cache/oncokb/`. Does the on-prem deployment have a standing cache-sync job, or do we need to add one?
3. Should the `curate_treatments` module write its results to a persistent table for audit (traceability for every curated row that entered a report)? I did not specify this — it is a downstream v2.3 "report-DB" concern but worth confirming.

### For **report-dev**

1. Is the existing `templates/cancer/report.html` layout flexible enough for a running patient header on every page, or does adding a header require a Jinja layout overhaul? The plan assumes this is a small change; please correct the effort estimate if it is not.
2. PMID superscript rendering (B3) — is there an existing reference-numbering utility in the codebase, or do we need a new one? I did not find one.
3. Should the curated Treatment Options table replace the current "Treatment Options" section entirely, or live alongside it (with a "curated" label) during a transition period? I favored replacement (cleaner) but a side-by-side A/B phase is defensible if the pilot team wants it.
4. The `clinvar_override_reason` footnote (A4) — where should it render? I suggested "under the variant card with the override applied", but if the current variant-card renderer is rigid, we may need a separate "Override Notes" block at the bottom of the Classifications page.

### For **qa-engineer**

1. The regression test that matters most to me is **"no 'futibatinib' in the cancer board output for any KRAS variant"** — a trivial string check that would have caught the incident that motivated this plan. Please add this as a CI assertion, not just an ad-hoc test.
2. For A3 (PM1 table), I listed a gene set (TP53, KRAS, NRAS, HRAS, BRAF, PIK3CA, EGFR, IDH1, IDH2). I want the test suite to have one positive case and one negative case per gene — please confirm that's a reasonable coverage target, or propose a leaner one.
3. For A4 (ClinVar override), a negative case matters more than the positive: I want a test that proves the override does NOT fire for a non-hotspot residue even when engine LP is reached through unrelated evidence. Please design that test carefully — it is the safety net against the override expanding silently.
4. The patient header test (B4) needs a golden HTML snapshot. Is the snapshot approach already in the test suite, or should I use a string-contains assertion?

## Risk summary

| Risk | Phase | Severity | Mitigation |
|---|---|---|---|
| Prompt rewrite destabilises cancer Board snapshots | A2 | H | Run v2 showcase fixture through new prompt; fallback to template-renderer Chair if needed |
| OncoKB rate-limit / offline | A1 | M | `offline_mode` config flag; local cache |
| ClinVar override expands silently | A4 | M | Narrow gate + override reason field + negative test (qa question #3) |
| CIViC drug-name merge fragile | A1 | M | Normaliser + unmerged-near-match log |
| Patient PHI leak into KB | B4 | H | `reporting.persist_patient_metadata: false` by default; code review |

## Final thought

The v2.2 scope is intentionally narrow. The temptation is to bundle clinical trials, MSI, and tumor-type narrative into the same plan because they all showed up in the F1 comparison. I resisted. The two most dangerous things BIKO does right now are (a) hallucinating drug/target pairings and (b) mis-classifying a canonical DBD hotspot as VUS. Fix those, add the credibility polish (citations, demographics), and ship. Everything else is v2.3.

— clinical-advisor
