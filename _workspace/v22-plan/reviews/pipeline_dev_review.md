# Pipeline-dev Review — v2.2 Plan

**Reviewer:** pipeline-dev (team v22-plan)
**Date:** 2026-04-14
**Plan:** `docs/superpowers/plans/2026-04-14-ai-board-v2.2.md`
**Cover note:** `_workspace/v22-plan/draft_cover.md`
**Verdict:** **approve-with-changes** — core architecture is sound, but three integration gaps must be closed before execution starts.

---

## TL;DR

The curate-then-narrate pattern is the right fix and the ACMG PM1 refinement is well-motivated. Three concrete blockers, though:

1. **The LLM-hallucination guard is leaky.** The plan's post-filter strips `treatment_options[]` rows without a `curated_id`, but it does NOT scrub free-text fields (`therapeutic_implications`, `therapeutic_evidence`, `clinical_actions`, `actionable_findings`). Those are *exactly* where "based on general knowledge, futibatinib is a tyrosine kinase inhibitor…" sentences would land today.
2. **`AgentOpinion.references` is already `List[str]` on `models.py:23`.** Replacing it with `List[Reference]` is a silently-breaking ABI change for the KB `raw_opinion_json` round-trip and for every cached showcase JSON under `docs/showcase/*.json`.
3. **The Board Chair never sees `_curated_treatments`.** Today the Chair synthesizes from agent-opinion *text* only (`board_chair._build_prompt` at `board_chair.py:217`). If the Chair is asked to enforce `curated_id`, it must receive the curated rows directly — but the plan only wires them into the Clinical Evidence Analyst's domain sheet.

Fix those three, tighten the effort numbers, and this is green for execution.

---

## Section-by-section critique

### Architecture / Curate-then-narrate (A1 + A2) — **concerns**

**Module placement (answers cover-note Q1).** `scripts/clinical_board/curated_treatments.py` is the right home. The consumer is the Board; none of the existing `scripts/somatic/` modules (`amp_tiering.py`, TMB) import from `clinical_board`, and the curator is a Board-only artifact. `scripts/clinical/oncokb_api.py` as a thin client is also fine — note that the existing `scripts/clinical/oncokb.py` is purely a cancer-gene-list JSON loader (`is_cancer_gene`, `get_cancer_gene_info`) and does **not** collide with the proposed new file. Consider calling the new one `oncokb_client.py` for intention-revealing symmetry with `ollama_client.py`, but this is a nit, not a blocker.

**`report_data["_curated_treatments"]` naming (Q2).** Correct convention. It matches `_board_variants` and `_board_selection_metadata` exactly (`runner.py:123-124`). Keep the underscore prefix.

**HGVSp → CIViC variant-name normalization — missing.** The CIViC schema uses its own variant name column (`SELECT * FROM evidence WHERE gene = ? AND variant = ?` at `query_civic.py:53-56`). CIViC names are things like `"G12D"`, `"V600E"`, `"AMPLIFICATION"`, `"EXON 19 DELETION"` — *not* HGVSp. The plan invokes `scripts.db.query_civic.get_variant_evidence(gene, variant_name)` but never explains how `variant_name` is derived from the variant dict's `hgvsp`. The project already has `scripts/common/hgvs_utils.py` for HGVSp → CIViC conversion (per `genomeboard-conventions`). The plan **must** name that module and budget the integration, otherwise `curate_treatments()` will silently miss variant-specific hits and fall through to gene-level — effectively reverting to the current (hallucinogenic) behaviour for CIViC.

**Drug-name merge fragility.** The cover note already flags this (A1 risk M). A concrete hardening step: every unmerged near-match should be logged to `data/cache/oncokb/unmerged_near_matches.log` with `(variant_key, oncokb_drug, civic_drug, distance)` so the curator can be tuned in-place without code changes. This is already implicit in the risk note; make it a task step.

**Clinical Evidence Analyst prompt (A2) — the real hallucination surface.** The current `ClinicalEvidenceAnalyst.system_prompt` at `scripts/clinical_board/agents/clinical_evidence.py:18-47` is a free-form CIViC narrative prompt with no curated-list constraint at all. The plan's added "CURATED EVIDENCE (authoritative source)" + "금지사항 — 목록에 없는 약물을 발명하지 마시오" block is the right rewrite.

**But: the Cancer Board Chair is where synthesis happens, not the Analyst.** Inspect `agents/board_chair.py:217-252` — the Chair's prompt is literally `case_briefing + formatted agent opinions`. The Chair never sees the curated rows directly; it only sees whatever the Analyst wrote about them. So:

- The Analyst paraphrases a curated row into Korean narrative.
- The Chair reads that narrative + the three other agents' text.
- The Chair synthesizes a `treatment_options[]` JSON list.

The plan wants the Chair to emit `curated_id` per row. **How does the Chair know which `curated_id` applies?** The curated list is on `report_data["_curated_treatments"]` but the prompt builder at `board_chair._build_prompt` never reads `report_data`. The plan must:
(a) extend `BoardChair.synthesize()` signature to take the curated list (or read it from `report_data` if the method gets that), AND
(b) inject a "CURATED EVIDENCE (authoritative, cite by curated_id)" block into the Chair prompt alongside the agent opinions.

Without (a) + (b), the `curated_id` post-filter drops *every* row because the Chair never emits a valid `curated_id` in the first place.

**Narrative hallucination — the critical gap.** A `CancerBoardOpinion` has:
```
therapeutic_headline, therapeutic_implications, therapeutic_evidence,
treatment_options[], actionable_findings[], clinical_actions[],
immunotherapy_eligibility, monitoring_plan[]
```
The plan's post-filter only guards `treatment_options[]` (rows without `curated_id` are dropped). Every other field is free-text. A clinician reads `therapeutic_implications` first — that's the executive summary rendered as the section body in `render.py`. If the LLM writes "…futibatinib has shown activity in KRAS-mutant cohorts…" into `therapeutic_implications`, the post-filter does nothing. **This is the original futibatinib failure mode and the plan does not close it.**

Required mitigations (pick at least one):
1. **Scrub narrative** — after the Chair responds, extract all drug names from `therapeutic_implications | therapeutic_evidence | clinical_actions | actionable_findings` (a simple noun-phrase scanner against a drug vocabulary; use the curated list's drug names as the positive allow-list) and reject the response if any out-of-list drug appears. This is the closest analogue to how F1-CDx isolates curated claims.
2. **Structured-narrative Chair** — force the Chair to emit `therapeutic_implications` as a list of `{claim, curated_ids: [...]}` objects, and have `render.py` format them into prose. Then the post-filter can drop unsupported claims.
3. **Template-renderer Chair fallback** (cover note's escape hatch) — zero LLM for `treatment_options`, pure curated-row rendering. This guarantees no hallucination in that field but leaves `therapeutic_implications` still free-text. **Only an acceptable fallback if combined with option 1.**

Recommendation: add option 1 as the baseline guard in A2. It is cheap, auditable, and would have caught the futibatinib incident.

---

### ACMG PM1 + ClinVar override (A3 + A4) — **approve with one nit**

**A3 table design.** Sound. The evidence collector's current PM1 rule is at `evidence_collector.py:193-194`:
```python
if _is_missense(consequence) and _has_domain(variant):
    codes.append("PM1")
```
`_has_domain` just checks VEP `DOMAINS` is non-empty. The plan's extension (consult `data/pm1_hotspot_domains.json` when `_is_missense` is true, fire PM1 if either the table or `_has_domain` matches) is a clean, minimal change that preserves existing behaviour.

`_extract_protein_position` already handles both 3-letter (`p.Arg249Met`) and 1-letter (`p.R249M`) HGVSp forms (`evidence_collector.py:74-89`), so R249M will be matched. Multi-residue mutations (rare) are not handled; the plan should acknowledge they fall through to the `_has_domain` path.

**Combination edge case in the pseudo-code.** The plan's snippet:
```python
if table_hit == "moderate" or domain_hit:
    codes.append("PM1")
elif table_hit == "supporting":
    codes.append("PM1_Supporting")
```
is correct. Double-check: a variant with `table_hit="supporting"` AND `domain_hit=True` falls into the first branch (PM1, not PM1_Supporting), which is what we want — if VEP's domain annotation agrees, upgrade to full PM1. Fine.

**Over-triggering risk (cover-note question implicit).** The plan's TP53 range `245–249` is wider than I would defend on clinical-advisor's own ClinGen SVI framing. Giacomelli 2018 identifies R175, R248, R273 as the top-3; 245–249 as a block is a compromise. This is a **clinical judgement call, not a pipeline concern**, but the plan should ask clinical-advisor whether `245–249` should be `PM1_Supporting` (safer, combines with PP3+PM2 to reach LP) or full `PM1` (moderate). Flag for clinical re-review in Task 6.

**Regression on 722-test baseline.** The plan says "existing tests that rely on VEP DOMAINS-based PM1 still pass." Confirmed by reading the collector — the new table is an OR, not a replacement. Good.

**A4 override — narrow and defensible.** The gate is:
`engine ∈ {LP, P}` AND `clinvar == "Conflicting"` AND `PM1 (table)` AND `PM5`.
All four conditions are independently verifiable. I prefer this to an "inline in `classify()`" flag because `acmg_engine.classify_variant` currently takes only `(evidences, gene)` and has no access to the ClinVar significance. The ClinVar override path is currently at `acmg_engine.apply_clinvar_override` (line 129, called by `orchestrate._classify_variants`, not `classify_variant` itself) — so A4's reconciliation should live in `apply_clinvar_override` or a new sibling function, NOT in `classify()`. The plan says "inline in classify()" in the cover note Q3 — **this is wrong**, the override must live at the post-classify reconciliation seam where ClinVar data is already in scope.

**Answer to Q3:** post-classification reconciliation pass. Not inline. The existing `apply_clinvar_override` is already that seam; add `apply_hotspot_conflict_reconciliation` alongside it or extend `apply_clinvar_override` to accept the engine's evidence codes as an argument.

**`clinvar_override_reason` field.** Must be plumbed from `classify()` (or the reconciliation function) → `orchestrate._classify_variants` → variant record dict → `render.py`. Currently `ClassificationResult` has only `classification`, `evidence_codes`, `conflict`, `matched_rule`. Adding `clinvar_override_reason` is one-line + one plumbing site. Effort estimate holds.

---

### Selector polish (B1 + B2) — **approve**

**B1 consequence gate.** Correct read. `_cancer_must_reason` at `variant_selector.py:251-270` today does not gate on consequence at all — any Tier I/II/III hits `return reason` before checking the VEP term. So a Tier III intronic variant in an OncoKB gene currently passes. The plan's fix is to apply `_PROTEIN_IMPACTING_CONSEQUENCES` as a universal gate at the top of both `_cancer_must_reason` and `_cancer_may_reason`. Straightforward.

One subtlety: a **P/LP** variant with non-protein-impacting consequence exists in practice (e.g. a deep-intronic ClinVar Pathogenic that was called P by ClinVar itself). Applying the gate to `_cancer_must_reason` would drop it silently from the board. I'd **exclude the `P_LP` branch from the consequence gate** — if ClinVar or the engine said it's P/LP, we respect that regardless of consequence. Apply the gate only to `Tier_I`, `Tier_II`, `Tier_III_hotspot`, `Tier_III_oncokb_gene`, and the MAY-path reasons.

**SpliceAI rescue threshold (≥0.2).** Reasonable. The field to read is `in_silico["spliceai_max"]` from the variant dict, which is already populated by the intake path. The rescue only matters for variants whose primary consequence is `splice_region_variant`, `intron_variant`, or `synonymous_variant` — add that branch explicitly.

**B2 MMR carve-out.** Safe. `_ACMG_SF_V32` at `variant_selector.py:39-48` includes the Lynch panel, but the SF exclusion is only applied in `_select_rare_disease`, not in `_select_cancer`. So adding `_MMR_LYNCH_GENES = frozenset({"MLH1","MSH2","MSH6","PMS2","EPCAM"})` and a MAY branch for cancer mode does not conflict with SF handling. Good.

**Separate concern — not raised in the plan:** PMS2/MLH1/MSH2/MSH6/EPCAM VUS are *silently excluded* from the **rare-disease MAY list** today (`_select_rare_disease` drops SF genes before HPO scoring). This is a cross-mode inconsistency: cancer mode will now surface Lynch VUS unconditionally, but rare-disease mode still drops them. Either this is intentional (SF opt-in flow, per the docstring) or it's a second bug. Flag for Task 6 clinical re-review.

**Priority renumbering (5.5 vs integer shift).** Renumber to integers. `_REASON_PRIORITY` at `variant_selector.py:79-89` is consumed by tuple sort keys; a float works but a renumber keeps the file tidy.

---

### PMID references ABI (B3) — **concerns**

**Already-existing field — this is a breaking change.** `AgentOpinion.references: List[str]` already exists on `models.py:23` as a plain list of strings (PMIDs or citation snippets). The plan's proposal is to replace this with `List[Reference]` where `Reference` is a new nested dataclass. Consequences:

1. **KB round-trip break.** `runner.py:206` calls `kb.save_decision(..., raw_opinion_json=json.dumps(board_opinion.__dict__, default=str))`. With `Reference` dataclasses nested inside, `default=str` will serialize each Reference as its `repr(...)` string — losing all structure. **Fix:** use `dataclasses.asdict(board_opinion)` before `json.dumps`. The plan must call this out as a modification to `runner.py`, not just `models.py`.

2. **Existing-JSON replay break.** `scripts/tools/rerender_report.py:25-32` uses `_filter_dataclass_kwargs` to tolerate missing keys, and `references` has `default_factory=list`, so absent → `[]` is fine. **But cached JSONs that already contain `references=["30224644"]` (string form) will load into a field typed `List[Reference]`** — and then `render.py` will call `ref.pmid` on a bare string and crash.

   **Fix options:**
   - (a) Type the field as `List[Union[str, Reference]]` with a renderer-side coercion (`ref.pmid if isinstance(ref, Reference) else str(ref)`).
   - (b) Add a `__post_init__` on `AgentOpinion` that upgrades string references to `Reference(pmid=s, quote="", context="", source="legacy")`.
   - (c) Rename the new field to `references_v2` and keep the old `references: List[str]` untouched. Cleanest, but the plan would need to rename in all downstream touch points.
   
   Pick (b). It's invisible to callers, works with the filter in rerender, and requires one method on the dataclass.

3. **The cached showcases under `docs/showcase/sample_codegen_777_report_v2.json` et al. will all hit (2).** I'd audit the set before committing the change — grep for `"references"` in those JSONs.

**Effort estimate is too low.** Plan says M (1.5–2 days). With (1) and (2), realistic is **M → 2.5–3 days**. The code change is small but the serialization cleanup + backfill for cached JSON is hidden cost.

**Schema flexibility for LLM output.** The plan already mentions LLMs emitting inconsistent PMID formats ("PMID: X", "pmid:X", bare numbers). Good — normalize in `_parse_cancer_response` at `board_chair.py:294-330` with a single regex.

---

### Patient demographics header (B4) — **approve**

CLI-args → `report_data["patient"]` → render is the right flow. The PHI-safety config flag `reporting.persist_patient_metadata: false` is correct.

Two minor plumbing notes:
- `run_clinical_board()` never touches CLI args directly; the entrypoint is `orchestrate.run_pipeline()`. The new CLI args live on `orchestrate.py` (argparse) and the `report_data["patient"]` is already a dict so no model change needed. Clean.
- The `_summarize_context` helper at `runner.py:62-75` is called for KB storage and currently includes `sample_id` and `clinical_note`. It should **not** be extended to include patient identifiers — that would defeat the `persist_patient_metadata` gate. Explicitly call this out as a task acceptance criterion: "KB context summary does not include patient name/DOB/MRN/facility even when the header is rendered."

---

### Showcase regen + docs (C1 + C2) — **approve**

No pipeline-dev concerns. Read-only documentation.

---

## Effort estimate review

| Task | Plan estimate | My estimate | Reason |
|---|---|---|---|
| A1 | M (2–3d) | **M (3d)** | HGVSp→CIViC variant-name mapping is unaccounted in the plan |
| A2 | M (2–3d) | **M (3d)** | Must thread curated list into Chair prompt AND add narrative scrubber |
| A3 | S (1d) | **S (1d)** | Fine |
| A4 | S (1d) | **S (1d)** | Fine — but the override must live at the post-classify seam, not inline in `classify()` |
| B1 | S (<1d) | **S (1d)** | P/LP consequence-gate exclusion adds a branch |
| B2 | S (<0.5d) | **S (<0.5d)** | Fine |
| B3 | M (1.5–2d) | **M (2.5–3d)** | Backward-compat coercion for old `List[str]` dumps is non-trivial |
| B4 | S (1d) | **S (1d)** | Fine |
| C1, C2 | <1d each | <1d each | Fine |
| **Total** | **10–13 days** | **12–15 days** | |

With two agents parallel (pipeline-dev on A1+A2+A3+A4, db-dev+report-dev on B-phase), wall-clock slips from ~6–8 days to **~8–10 days**. Still fits a 2-sprint window.

---

## Missing items

1. **HGVSp → CIViC variant-name normalization** — already discussed above.
2. **Board Chair curated-list injection** — already discussed above.
3. **Narrative scrubber for free-text opinion fields** — already discussed above.
4. **`raw_opinion_json` serialization with nested Reference dataclasses** — already discussed above.
5. **Rare-disease Lynch-panel handling** — cross-mode inconsistency, flag for Task 6.
6. **Mode-gating of `curate_treatments()`** — plan does not say "cancer-mode only". Rare-disease `run_clinical_board()` should not call the curator. Add as a guard in `runner.py` after `mode == "cancer"` check.
7. **Regression test for `raw_opinion_json` KB round-trip** — qa-engineer should own this, but the plan's B3 task must list it.

---

## Prioritized punch list (in order)

1. **[BLOCKER] Add narrative scrubber to A2.** Post-processing must reject any Chair response whose `therapeutic_implications / therapeutic_evidence / clinical_actions / actionable_findings` text mentions a drug not in the curated list. Without this, the futibatinib incident is not actually prevented — only relocated from `treatment_options[].drug` to `therapeutic_implications` prose. The regression test must grep the full `CancerBoardOpinion.__dict__` serialization, not just `treatment_options`.

2. **[BLOCKER] Inject `_curated_treatments` into the Board Chair prompt.** `BoardChair.synthesize()` currently has no access to `report_data`. Extend the signature to `synthesize(briefing, opinions, mode, curated_treatments=None)` and inject a `## CURATED EVIDENCE` block into `_build_prompt` alongside the agent-opinion text. Without this, the Chair cannot emit a valid `curated_id` and the A2 post-filter will drop every `treatment_options[]` row.

3. **[BLOCKER] Fix `AgentOpinion.references` ABI break.** The field already exists as `List[str]`. Add a `__post_init__` coercion that upgrades bare strings to `Reference(pmid=s, quote="", context="", source="legacy")`, and switch `runner.py:206` from `json.dumps(__dict__, default=str)` to `json.dumps(dataclasses.asdict(board_opinion))`. Without this, cached showcase JSONs crash on re-render and KB persistence silently loses structure.

4. **[HIGH] Resolve HGVSp → CIViC name mapping path.** Name `scripts/common/hgvs_utils.py` (or the equivalent) as a dependency in A1 and add a test case for a variant whose HGVSp doesn't cleanly map to a CIViC name (expect fall-through to gene-level + logged near-match).

5. **[HIGH] Move A4 reconciliation out of `classify()` into `apply_clinvar_override` or a sibling function.** `classify_variant` has no access to ClinVar data. The cover note Q3 answer is "post-classification reconciliation pass."

6. **[MED] Exclude P/LP from the B1 consequence gate.** `_cancer_must_reason` should preserve the existing "P/LP always included" invariant even when a variant is non-protein-impacting.

7. **[MED] Update B3 effort estimate to M (2.5–3d).** Serialization cleanup and legacy-JSON coercion are not captured in the 1.5–2d estimate.

8. **[MED] Mode-gate `curate_treatments()` to cancer mode only.** Add an explicit `if mode == "cancer":` guard in `runner.py` so rare-disease paths are unaffected.

9. **[LOW] Add unmerged-drug-name log to A1 for ongoing tuning.** One-line `logger.warning(...)` in the curator merge loop.

10. **[LOW] Rename `oncokb_api.py` → `oncokb_client.py`** for intention-revealing parallelism with `ollama_client.py`. Nit only.

---

## Verdict

**approve-with-changes.**

The plan is clinically and architecturally sound — curate-then-narrate is the right pattern, the PM1 table is well-motivated, and the selector tightening is overdue. The three blockers above are all implementation-integration gaps, not design flaws. If the author folds items 1–3 into the plan (tighten A2's guard, thread curated list through the Chair, handle the `references` ABI break), this plan is ready to execute.

Blocking gate for my sign-off:
- A2 must grow a narrative scrubber step (item 1).
- A2 must wire curated rows into `BoardChair.synthesize()` (item 2).
- B3 must include the `__post_init__` coercion and the `asdict` swap (item 3).

Without those three, executing the plan as written would leave the futibatinib failure mode intact, break the KB JSON round-trip, and produce a Chair that can't satisfy the curated_id contract. With them, the plan is a clean path to the clinical pilot gate.

— pipeline-dev
