# v2.2 Plan — DB-dev Review

**Reviewer:** db-dev
**Date:** 2026-04-14
**Plan:** `docs/superpowers/plans/2026-04-14-ai-board-v2.2.md`
**Scope of this review:** data sources, SQLite schemas, build/query modules, caching, licensing, version tracking, on-prem operability.

---

## 0. Executive verdict

**Approve with required revisions.** The curate-then-narrate architecture is the correct structural fix and maps cleanly onto the existing CIViC/ClinVar/OncoKB layer. None of the four Phase A tasks need a new SQLite database to be built — this is good news for the pilot timeline. However, five DB/data-layer issues must be resolved before A1 can be implemented as written:

1. **OncoKB free-tier licensing is under-specified** and on-prem-unfriendly. The plan calls it "free academic/non-commercial variant API" but does not address the registered-token requirement, offline-token rotation, or the redistribution clause that matters for a Korean clinical partner.
2. **Existing assets are not acknowledged.** `data/db/hotspots/hotspots_v2_single.tsv` (cancerhotspots.org) already exists and is unused; CIViC already has a `hotspots` table; neither is referenced in A3. The plan risks creating redundant curation.
3. **`pm1_hotspot_domains.json` is not integrated with `version_manager.py`** — the project convention is that every curated data table appears in `db_versions` so reports show provenance. The plan's "top-level `version` field" is necessary but not sufficient.
4. **The OncoKB/CIViC drug-name merge key is fragile** and not aligned with any controlled vocabulary. The plan's "lowercase + strip brand suffixes" will mis-merge or miss-merge on common cases (Herceptin/trastuzumab, binimetinib/MEK162, etc.).
5. **Cache layout and TTL are not aligned with existing `scripts/common/cache.py`** conventions. A new `data/cache/oncokb/` tree is fine but must go through the same idempotency/concurrency path the rest of the pipeline uses, or it will drift.

None of these are blockers for the clinical direction. All five can be resolved with a ~1-paragraph edit each to the affected task sections; together they add roughly 0.5–1 engineering day to Phase A.

---

## 1. Section-by-section critique

### 1.1 Task A1 — Curated Treatments Module

**What is correct.**

- The module boundary (`scripts/clinical_board/curated_treatments.py`) and the thin wrapper (`scripts/clinical/oncokb_api.py`) match the existing convention in `scripts/clinical/` (cf. `oncokb.py`, `query_clinvar.py`, `hpo_matcher.py`). A1 can reuse `scripts.common.config.get` for paths and `scripts.common.cache` for response caching — both patterns are in place.
- `query_civic.get_variant_evidence(gene, variant_name)` and `get_predictive_evidence_for_tier(gene, hgvsp)` already return `pmid`/`citation` fields and are sorted by evidence level. The curator can wrap these without any schema change to `civic.sqlite3`. No CIViC rebuild required.
- Offline mode is explicit. Good — on-prem deployments may run fully disconnected and the plan names the flag (`clinical_board.curated_treatments.offline_mode`) and the degradation path (CIViC-only).

**Where the plan is wrong or incomplete.**

**(1a) OncoKB free-tier licensing is under-specified.** The plan says "OncoKB free academic/non-commercial variant API". In practice the OncoKB public API requires:

- a registered API token (`https://www.oncokb.org/account/register`), even for the variant-annotation endpoint;
- acceptance of the OncoKB Terms of Use, which restrict **commercial and clinical** redistribution of OncoKB-derived content;
- an explicit academic-license agreement for non-academic institutional use.

For a Korean *clinical pilot* (memory: `feedback_onprem_local_db.md` — on-prem operation), the free-tier ToS are marginal at best. The plan should either (a) scope OncoKB to research-only pilot showcases and *not* into clinical reports, or (b) commit to obtaining a commercial/clinical OncoKB license before the pilot begins. Silently embedding the free-tier API into a Korean clinical pipeline is a licensing landmine. Open question #1 in the cover note partly addresses this ("commercial license is a drop-in config change") but the default path in A1 is the free API — the legal posture needs to be explicit in the plan itself, not deferred to a config flag.

**(1b) Cached OncoKB responses must go through `scripts/common/cache.py`, not a parallel cache tree.** The existing pipeline has a thread-safe SQLite cache (`data/cache/variant_cache.sqlite3`, 7-day TTL) used by every API wrapper (`scripts/common/api_utils.py`, HPO, ClinVar, gnomAD). Creating a new `data/cache/oncokb/` JSON-file tree duplicates invalidation logic, loses the TTL, and breaks concurrency guarantees. Recommendation: extend `variant_cache.sqlite3` with a `namespace` column (or an `oncokb_*` key prefix) and let the existing cache handle TTL, locking, and cleanup. This also means one cache-sync job (not two) for on-prem deployments, answering open question #2.

**(1c) Offline graceful-degradation must be surfaced to the reader, not silent.** If `offline_mode=True` and OncoKB is skipped, the rendered Treatment Options table will be CIViC-only — potentially missing rows OncoKB would have returned. The report must display a "curation source: CIViC-only (OncoKB offline)" footnote. This is a 2-line rendering change but belongs in A1's acceptance criteria so report-dev picks it up.

**(1d) Drug-name merge key is fragile.** The plan proposes "lowercase + strip brand suffixes" for the OncoKB ↔ CIViC merge. Known failure cases:

- `Herceptin` (brand) vs `trastuzumab` (INN) — these are the same drug.
- `Binimetinib` vs `MEK162` (code name) — same drug.
- `Trametinib + hydroxychloroquine` (F1's KRAS row) — a combo; simple string match will miss either half.
- Salt/form suffixes: `Imatinib mesylate` vs `imatinib`.

The plan should either (i) use the CIViC `therapy_id` (if CIViC exposes a stable therapy identifier in its TSV — worth checking `scripts/db/build_civic_db.py:102-130`; the current schema only stores the `therapies` free-text column, so a schema extension to store `therapy_ids` from ClinicalEvidenceSummaries.tsv would be a one-commit change), or (ii) accept the ambiguity and gate merges through a hand-curated synonym list checked into `data/drug_synonyms.json`. Option (i) is cleaner and should be a sub-task of A1.

**(1e) `variant_key` uniqueness assumption.** The plan keys curated rows by `{chrom}:{pos}:{ref}:{alt}`. This is fine for genomic variants but breaks for HGVSp-only inputs (the curator receives variants that may have been selected by consequence, not by coordinate). The KRAS G12D unit test in Step 1 of A1 uses both a genomic key *and* an HGVSp — curator inputs should be normalized to genomic-key-required. Add this to acceptance criteria: "curator rejects variants lacking `chrom/pos/ref/alt` with a clear ValueError".

**Db-dev punch list for A1:**

- [ ] **A1-db-1** Add a "Licensing" paragraph: explicit ToS posture for OncoKB free-tier in a clinical pilot; recommend a formal academic-license or commercial-license step before pilot start.
- [ ] **A1-db-2** Route OncoKB caching through `scripts/common/cache.py` (add a `namespace` or key-prefix column to `variant_cache.sqlite3`). Delete the proposed `data/cache/oncokb/` path. Update config to reference `cache.path` (already in config.yaml).
- [ ] **A1-db-3** Add acceptance criterion: "report must display curation-source footnote when `offline_mode` causes CIViC-only degradation".
- [ ] **A1-db-4** Extend `civic.sqlite3.evidence` schema with a `therapy_ids` column populated from `ClinicalEvidenceSummaries.tsv` (one-line build-script change + re-run); use it as the primary merge key. Fall back to string match only when `therapy_id` is missing.
- [ ] **A1-db-5** Add explicit input contract: curator requires genomic coordinates on every input variant.

### 1.2 Task A2 — Narrate-only prompts

No DB/data-layer concerns. A2 is prompt engineering and runner plumbing. Two minor observations:

- The `_curated_treatments` stash key on `report_data` is fine and matches existing `_board_variants` / `_board_selection_metadata` conventions (answers open question #2 for pipeline-dev).
- The fallback "template-renderer Chair" (cover-note bullet) is a sound escape hatch. If it fires, the curator's rows become the sole source for the Treatment Options table, which places *even more* weight on A1's correctness — reinforcing punch list items A1-db-2/3/4.

### 1.3 Task A3 — PM1 Protein-Domain Hotspot Table

**What is correct.**

- Keying by `(gene, residue_range)` with `strength: moderate|supporting` is the right abstraction. It also matches Garrett 2021 (ClinGen SVI PM1 refinement) and is portable across VCEPs.
- Storing as JSON (not SQLite) is appropriate for a small, read-only, git-tracked table. No DB build step needed.
- The gene list (TP53, KRAS, NRAS, HRAS, BRAF, PIK3CA, EGFR, IDH1, IDH2) covers the highest-yield cancer hotspots.

**Where the plan is wrong or incomplete.**

**(3a) Existing assets are not acknowledged.** Two hotspot data sources already exist in the repo and are not mentioned:

1. `data/db/hotspots/hotspots_v2_single.tsv` — this is cancerhotspots.org v2 single-residue hotspots (MSK). ~400 gene-position pairs with statistical evidence. It is downloaded and staged but currently unused by `evidence_collector.py`.
2. `civic.sqlite3.hotspots` — a CIViC-derived hotspot table already built by `scripts/db/build_civic_db.py:132-154`. Spot-check: `KRAS 12/13/61/146` ✓, `TP53 175/248/273` ✓, but `TP53 249` is **missing** (because CIViC's single-variant naming convention for R249M didn't match the regex on line 145, or there's no R249 SNV in CIViC variants). This is exactly the BIKO bug the plan is trying to fix.

The plan should explicitly state: "`data/pm1_hotspot_domains.json` is the **ACMG PM1 authoritative** table. It supersedes CIViC's derived hotspots and cancerhotspots.org for PM1 logic — cancerhotspots/CIViC may still be used for Tier-III selector hotspot signalling, which is a separate concern." Without this clarification, a future maintainer will not know which of three tables `evidence_collector._pm1_hotspot_match` should consult, and the existing two sources will silently drift from the new JSON.

**(3b) Version tracking.** The plan says the JSON has a top-level `version` field. That is necessary but `scripts/db/version_manager.py:get_all_db_versions()` must also surface it so the rendered report shows "PM1 hotspot table v1.0 / 2026-04-14". Otherwise a reader cannot audit which hotspot set was active when the report was generated. Add a `versions["PM1_Hotspots"]` block to `version_manager.py` (2-line change) and list this as an A3 modified-file.

**(3c) Bootstrapping the JSON from existing sources.** The plan writes the JSON by hand from ClinGen VCEP publications. Faster and more defensible: a one-off `scripts/db/build_pm1_hotspots.py` that merges (a) cancerhotspots_v2_single.tsv, (b) ClinGen SVI VCEP tables for TP53/BRCA1/BRCA2/PTEN/MLH1/MSH2, and (c) clinical-advisor-supplied canonical entries (TP53 DBD L3 loop 245–249, KRAS G-domain 12–13/61). Output is the same JSON. This turns a hand-curated artifact into a reproducible build, and the resulting JSON still lives in `data/` and is git-tracked. It also answers open question #1 for db-dev (the build script is the traceability story; version_manager surfaces it).

**(3d) JSON schema.** The `$schema` and `version` fields should live alongside an explicit `build_date`, `source_refs: [{pmid, note}, ...]`, and `source_hash` (sha1 of the payload) so `version_manager.py` can detect tampering. This is parallel to the `metadata` table convention for SQLite DBs.

**Db-dev punch list for A3:**

- [ ] **A3-db-1** Plan must acknowledge `data/db/hotspots/hotspots_v2_single.tsv` and `civic.sqlite3.hotspots` and clarify precedence: `pm1_hotspot_domains.json` is ACMG-authoritative for PM1; the others serve Tier-III hotspot signalling.
- [ ] **A3-db-2** Add a new helper `scripts/db/build_pm1_hotspots.py` that merges cancerhotspots + ClinGen VCEP + clinical-advisor overrides → JSON. Keep it as part of A3 (same day; it's a dozen lines of Python + a merge dict).
- [ ] **A3-db-3** Register `versions["PM1_Hotspots"]` in `scripts/db/version_manager.py` and include `build_date`, `source_refs`, `record_count`.
- [ ] **A3-db-4** JSON schema: top-level `version`, `build_date`, `source_refs`, `source_hash`, then `genes` payload.

### 1.4 Task A4 — ClinVar Conflict Reconciliation

No new data source or schema. One small observation:

**(4a) The override is keyed off `clinvar_class == "Conflicting"`.** The current `scripts/db/query_local_clinvar.py` normalises ClinVar's `ClinicalSignificance` column — but the exact string "Conflicting classifications of pathogenicity" has shifted historically ("Conflicting interpretations of pathogenicity" pre-2023, now "Conflicting classifications of pathogenicity"). The override gate should match on a prefix (`startswith("Conflicting")`) or on the normalised category (`VARIANT_CATEGORY_CONFLICTING`), not the exact string. One-line hardening. Mention in A4 acceptance criteria: "gate fires on any ClinVar status whose normalised category is Conflicting, regardless of wording revision".

**Db-dev punch list for A4:**

- [ ] **A4-db-1** Normalise the ClinVar "Conflicting" match — prefix/category, not exact string.

### 1.5 Task B1 — Variant Selector consequence gate

No DB/data-layer change. Approve as written.

### 1.6 Task B2 — MMR/Lynch carve-out

No DB/data-layer change. The MMR gene set is hard-coded in Python which is appropriate (5 genes, stable). Approve.

### 1.7 Task B3 — PMID references

**What is correct.**

- All the PMID data already exists: `civic.sqlite3.evidence.citation_id` is populated (spot-checked: KRAS has 159 Predictive evidence rows, each with a PMID), and `pm1_hotspot_domains.json` will carry PMIDs natively. No new DB or source.

**Where the plan is incomplete.**

**(b3a) PMID normalisation.** The plan's `Reference` dataclass has `pmid: str`. Good. But CIViC stores `citation_id` which may be a bare PMID (`"30224644"`), a multi-id comma-joined string (`"30224644,33280026"`), or empty. The curator (A1) must split/normalise these before constructing `Reference` objects. Otherwise B3's test fixture (`references=[Reference(pmid="30224644", ...)]`) will not match what the curator produces for evidence rows that cite two papers. Add a helper `scripts/clinical_board/references.py::normalize_pmids(citation_id: str) -> list[str]` (6 lines) and use it from both A1 and B3.

**(b3b) Numbered-list dedup.** If the Board chair cites PMID 30224644 three times in the same opinion, B3 should render one numbered reference (`[1]`) cited from three superscripts, not three separate entries. The plan's acceptance criteria don't state this. This is a rendering detail but also a data-shape decision — the `References` list should be de-duped at the `CancerBoardOpinion` level, not at the render stage, so that the de-duped list is part of the stable data model.

**Db-dev punch list for B3:**

- [ ] **B3-db-1** Add `scripts/clinical_board/references.py::normalize_pmids` helper for CIViC/OncoKB citation parsing; use from A1 and B3.
- [ ] **B3-db-2** Dedup `CancerBoardOpinion.references` by PMID at the model/runner level, not at render.

### 1.8 Task B4 — Patient demographics header

**What is correct.**

- The `reporting.persist_patient_metadata: false` default is the correct PHI posture for on-prem deployments.

**Where the plan is incomplete — PHI and the future Report-DB (v2.3).**

**(b4a) No PHI sink audit.** When `persist_patient_metadata: false`, the demographics must not leak into:
- `variant_cache.sqlite3` (cache keys must not include patient name/MRN)
- `_workspace/` artifacts (the existing codegen showcase writes JSON dumps)
- `gene_knowledge.json` re-builds
- any `version_manager` metadata

The plan says "gates downstream persistence" but does not enumerate the sinks. I recommend a small audit checklist in B4: "grep for `report_data['patient']` usage after this diff; confirm the only read sites are the render path and the case_briefing assembler". This is a 15-minute check, but it should be in the plan so it doesn't get skipped.

**(b4b) v2.3 Report-DB forward compatibility.** Memory: `project_v2_report_db.md` — v2.0 is about persisting report results to a DB. If B4 lands with `persist_patient_metadata: false` and v2.3 later needs patient lookups, the schema decision should already be made now: patient metadata lives in a *separate* table with a per-institution retention policy, not inline in the report row. The plan should note this constraint so the v2.3 schema does not inherit a PHI-inline design by accident. One sentence in B4's clinical-rationale paragraph.

**Db-dev punch list for B4:**

- [ ] **B4-db-1** Add PHI sink audit checklist: enumerate the read/write sites for `report_data["patient"]` and confirm the cache, `_workspace/`, and `gene_knowledge` paths cannot touch it.
- [ ] **B4-db-2** Add a forward-compatibility note for v2.3 report-DB: patient metadata must be in a separate table with its own retention policy.

### 1.9 Task C1 / C2 — Showcase regen / Docs

No DB concerns. Approve.

---

## 2. Cross-cutting DB-layer observations

### 2.1 Schema and volume — no new SQLite databases needed

Phase A adds zero new SQLite files. The only new data artifact is `data/pm1_hotspot_domains.json` (~10 KB) and a handful of cached OncoKB responses (tiny — ~1 KB per variant × maybe 50 variants per report = ~50 KB per report). The existing `civic.sqlite3` (~200 MB) and `variant_cache.sqlite3` (~200 MB with TTL) absorb all of A1's traffic. This is a good outcome for the clinical pilot — no schema migration, no version-manager rework, no re-indexing — as long as the A1-db-2 caching route is adopted.

### 2.2 CIViC therapy field extension (A1-db-4)

The one schema extension worth doing is a `therapy_ids` column on `civic.sqlite3.evidence`. This is a one-commit change:

```python
# scripts/db/build_civic_db.py — add to CREATE TABLE evidence
therapy_ids TEXT,          # comma-joined list from ClinicalEvidenceSummaries.tsv

# add to the INSERT
row.get("therapy_ids", ""),
```

`ClinicalEvidenceSummaries.tsv` from CIViC does expose a `therapy_ids` column (CIViC's internal therapy IDs). Using them as the merge key for the OncoKB↔CIViC join eliminates most of the drug-name-string fragility and reduces A1's risk from **Medium** to **Low**. The cost is running `python scripts/db/build_civic_db.py` once after the column is added.

### 2.3 On-prem cache-sync posture (open question #2)

The cover note asks if the on-prem deployment has a cache-sync job. Answer: **there is no standing cache-sync job today.** The existing `variant_cache.sqlite3` is populated lazily on first query and expires at 7 days. For an air-gapped on-prem pilot, this means:

- Every report against a previously-unseen variant triggers a first-query delay (or fails if air-gapped).
- A cache-warm step (`scripts/tools/warm_variant_cache.py` — not yet written) should be added as a pre-flight for the pilot.

This is out of scope for v2.2 but worth flagging to team-lead. Recommend tracking as a v2.2.1 operational task, not blocking the plan.

### 2.4 Version tracking summary

After v2.2, `scripts/db/version_manager.py::get_all_db_versions()` should also surface:

- `CIViC` (currently absent! — `civic.sqlite3` has a `metadata` table populated by the build script, but version_manager does not read it)
- `PM1_Hotspots` (new, from A3)
- `cancerhotspots_v2_single` (currently absent, for completeness — existing file, not used by PM1 but used by selector)

The CIViC omission is a pre-existing bug that v2.2 should fix as a bycatch — it is trivially two functions added to `version_manager.py`, and the rendered report will stop under-reporting its data provenance. Adding it to the plan as a 0.5-day bycatch item under Task C2 (documentation) is reasonable; alternatively, it can be its own v2.2.1 cleanup commit.

### 2.5 Licensing summary

| Data source | Licensing posture | v2.2 risk | Recommended action |
|---|---|---|---|
| **OncoKB free-tier API** (A1 new) | Registered token required; ToS restrict clinical/commercial redistribution | **HIGH** — direct clinical use ambiguous | Obtain formal academic/commercial license before pilot; explicit ToS paragraph in A1 |
| **CIViC** (existing, reused) | CC0 public domain | None | ✓ |
| **cancerhotspots.org** (existing, unused) | Open academic; citation required | None | ✓ |
| **ClinGen VCEP publications** (PMID-cited in A3) | Free to cite; tables are factual | None | ✓ (cite PMIDs in `source_refs`) |
| **ClinVar** (existing, reused) | Open | None | ✓ |
| **PubMed PMIDs** (B3, reused via CIViC) | PMIDs themselves are identifiers, not content; no licensing issue | None | ✓ |

The only licensing landmine is OncoKB. Everything else is clean.

---

## 3. Consolidated punch list

**Blockers (must be resolved before A1 implementation begins):**

1. **A1-db-1** — Explicit OncoKB licensing paragraph in A1. Confirm academic or commercial license is obtained before pilot.
2. **A1-db-2** — Route OncoKB caching through `scripts/common/cache.py` (existing `variant_cache.sqlite3`). Delete the proposed `data/cache/oncokb/` tree.
3. **A3-db-1** — Acknowledge `cancerhotspots_v2_single.tsv` and `civic.sqlite3.hotspots`; declare `pm1_hotspot_domains.json` as the ACMG-authoritative PM1 source.

**Required (must be resolved before A1/A3 merge):**

4. **A1-db-3** — Offline-mode degradation footnote in the rendered report.
5. **A1-db-4** — Extend `civic.sqlite3.evidence` with a `therapy_ids` column; use as primary merge key.
6. **A1-db-5** — Require genomic coordinates on curator input; reject with ValueError if missing.
7. **A3-db-2** — Add `scripts/db/build_pm1_hotspots.py` helper that merges cancerhotspots + ClinGen VCEP + clinical overrides into the JSON.
8. **A3-db-3** — Register `PM1_Hotspots` in `version_manager.py`.
9. **A3-db-4** — Rich JSON schema: `version`, `build_date`, `source_refs`, `source_hash`, `genes`.
10. **B3-db-1** — `normalize_pmids` helper for comma-joined `citation_id` strings.
11. **B3-db-2** — Dedup `references` by PMID at model/runner level.
12. **B4-db-1** — PHI sink audit checklist in B4.

**Recommended (fold into v2.2 if possible, otherwise v2.2.1):**

13. **A4-db-1** — Normalise ClinVar "Conflicting" match (prefix/category, not exact string).
14. **B4-db-2** — Forward-compatibility note for v2.3 report-DB patient-metadata separation.
15. **C2-db-1** — Fix the CIViC omission in `version_manager.py` (it already has a metadata table).

**Operational (out of scope for v2.2, track separately):**

16. **ops-db-1** — Cache-warm pre-flight script for air-gapped pilot.

---

## 4. Answers to the cover note's open questions (for db-dev)

> **Q1.** Does `pm1_hotspot_domains.json` need DB-version tracking alongside ClinVar/gnomAD version manifests?

**Yes.** A top-level `version` field in the JSON is necessary but not sufficient — `scripts/db/version_manager.py::get_all_db_versions()` must also surface it so the report renders provenance. See A3-db-3. In addition, CIViC is *currently* missing from version_manager despite having a metadata table; v2.2 should fix this as a bycatch (C2-db-1).

> **Q2.** Does the on-prem deployment have a standing cache-sync job, or do we need to add one?

**There is no standing cache-sync job today.** The existing `variant_cache.sqlite3` is lazy-populated with a 7-day TTL. For air-gapped pilot operation, a cache-warm pre-flight (`scripts/tools/warm_variant_cache.py`, not yet written) is required. This is out of scope for v2.2 but should be tracked as a v2.2.1 operational task. In the meantime, A1 should route OncoKB caching through `scripts/common/cache.py` (A1-db-2) so when the warm script is eventually written, it uses one path, not two.

> **Q3.** Should `curate_treatments` write results to a persistent table for audit?

**Not in v2.2. Yes in v2.3 (report-DB).** For v2.2, `_curated_treatments` lives on `report_data` in-memory per run and is rendered into the HTML. Traceability for v2.2 is: the rendered table shows `source` (OncoKB/CIViC/both), `curated_id`, and `pmids[]`, which lets a reviewer reconstruct the curation post-hoc. For v2.3 (report-DB), the curated rows should persist to a `curated_treatments` table keyed by `(report_id, variant_key, curated_id)` with full raw-row storage — but that is a v2.3 schema decision, not a v2.2 one. Blocking v2.2 on this would be scope creep. Add one sentence to A1 acknowledging the v2.3 forward-compatibility.

---

## 5. Verdict

**Approve with the 12 required revisions (punch list items 1–12).** The curate-then-narrate architecture is sound and the Phase A/B work is proportionate to the patient-safety risk. No new SQLite DB needs to be built. The CIViC `therapy_ids` column extension (A1-db-4) is the single highest-leverage technical improvement — it takes the drug-name merge from fragile to robust with one schema column and a build-script rerun. The OncoKB licensing posture (A1-db-1) is the single highest-leverage non-technical improvement — it prevents a legal surprise at pilot go-live.

Phase A effort revised with db-dev additions: **+0.5 day** on A1 (cache routing + `therapy_ids` column + licensing paragraph), **+0.5 day** on A3 (build helper + version_manager integration). New Phase A total: **~7–9 engineering days** (up from 6–8). Phase B unchanged. Recommend proceeding to Task 6 re-review with these revisions applied.

— db-dev
