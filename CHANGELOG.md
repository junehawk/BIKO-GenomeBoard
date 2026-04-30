# Changelog

All notable changes to BIKO GenomeBoard are documented here.

This project adheres to [Keep a Changelog 1.1](https://keepachangelog.com/en/1.1.0/)
and uses [Semantic Versioning](https://semver.org/).

BIKO GenomeBoard is a **research reference tool**, not a clinical
decision instrument. Every release carries the same disclaimer: output
is intended for independent review by a researcher or clinician.

## [Unreleased]

### Changed — v2.5.5 (Chair JSON simplification + rare-disease PGx skip)
- **Cancer Chair `variant_key` is now backfilled deterministically.** Before
  v2.5.5 the Cancer Board Chair LLM was required to emit
  `{curated_id, variant_key, drug, evidence_level, resistance_notes}` and the
  `narrative_scrubber` validated each row against a `(curated_id, variant_key)`
  pair set. Empirically (`data/sample_vcf/demo_variants_grch38_annotated.vcf`,
  NSCLC, 2026-04-30) supergemma4-31b miscopied `variant_key` on the long
  Chair prompt (23K chars), rewriting `chr13:32356550:C:T` as
  `chr13:32356 own:C:T` — a BPE/SentencePiece-level token confusion the
  prompt-side `"copy character-by-character"` instruction could not prevent.
  This dropped 25 % of LLM-emitted treatment rows on the demo and silently
  removed the BRCA2 platinum chemotherapy line (CIViC EID:879, OncoKB Level
  3B; NCCN NSCLC v.3.2024 first-line standard) from the report.
  The Chair prompt now omits `variant_key` from the required keys, and
  `scrub_opinion` looks up the curator's authoritative `variant_key` via
  `curated_id` and writes it onto the row before keeping it. `validate_treatment_option`
  validates `curated_id` only.
- **v2.2 Phase A "curate-then-narrate" contract strengthened, not weakened.**
  Pasting a valid `curated_id` into a narrative paragraph describing a
  different variant no longer silently drops the row — it is reassigned to
  the variant the curator originally bound to that `curated_id`. The
  EGFR-curated_id-pasted-under-TP53 attack still cannot put an EGFR drug
  under a TP53 row in the rendered table; instead the row appears under
  EGFR. `tests/test_board_chair_curated_id.py::test_cross_variant_paste_pinned_to_origin_by_scrubber`
  pins this new semantic.
- **Rare-disease mode no longer runs the PGx Specialist agent.**
  `scripts/clinical_board/runner.py::_load_agents` returns three agents
  (Variant Pathologist, Disease Geneticist, Literature Analyst) for
  `mode == "rare-disease"`. The PGx Specialist always returned 0 findings
  on rare-disease cases (no curated PGx data path), wasting ~15 s of GPU
  time per run. Cancer mode still runs four agents including PGx.
- **Non-goal preserved.** The `template_renderer_chair` fallback trigger
  (`post_scrub_count == 0 AND curator_nonempty`) is unchanged. The fallback
  still fires when the LLM emits zero treatment rows (or every row has a
  hallucinated `curated_id`); the only change is that variant_key
  miscopies no longer reach the trigger. Phase A fixture #1 (fabricated
  `curated_id`) regression remains green.
- **Therapeutic Target Analyst empty response (deferred).** The Cancer
  measurement also showed Therapeutic Target Analyst returning `{}` on
  the demo VCF (n=1). Hypothesis: medgemma:27b sees the curator's 88
  authoritative rows and concludes there is nothing to add. v2.5.5 does
  not modify this — the deterministic curator continues to produce its
  rows, the Chair narrative continues to integrate the other three
  agents, and the visible report is preserved. Treated as a separate
  agent-prompt-fit issue to be measured (n≥10) and addressed in v2.6.
- Rollback tag: `pre-v2.5.5`.

### Changed — v2.5.4 (batch-canonical refactor + report assembly hygiene)
- **Single source of truth for per-sample assembly.** New
  `scripts/orchestration/canonical.py` holds the canonical
  single-sample pipeline (`build_sample_report`). Both
  `scripts/orchestrate.py::run_pipeline` (CLI wrapper) and
  `scripts/orchestration/batch.py::run_batch_pipeline` (per-sample loop)
  route through it. The pre-v2.5.4 `batch.py` had its own ad-hoc
  pipeline that diverged from single-mode on half a dozen behaviours;
  that divergence is gone.
- **H1 — Variant reuse across samples fixed.** Pre-v2.5.4 batch mode
  keyed a shared `unique_variants` dict on `chrom:pos:ref:alt` and
  re-used the **first sample's Variant object** for every later sample.
  Variant carries per-sample annotation state (canonical transcript,
  HGVSp, trio genotype, selection_reason_list, source, matching_hpo,
  …) so later samples' reports silently inherited earlier samples'
  annotations — a correctness bug. The new batch calls
  `build_sample_report` once per VCF; Variant identity never crosses a
  sample boundary.
- **M2 — CIViC priority fix.** The pre-v2.5.4 renderer chained
  `setdefault()` calls on `finding_summary`, `treatment_strategies`,
  `references`, and `content_status`, filling from gene_knowledge
  first and CIViC second. Because `setdefault` cannot overwrite, the
  order meant **gene_knowledge always beat CIViC** — the opposite of
  what the comment claimed and of the v2.2 curate-then-narrate
  contract. Enrichment now lives in canonical (assembly time), and
  CIViC values win over gene_knowledge for the overlapping clinical
  fields. Caller-supplied values (e.g. a Board scrubber setting
  `content_status="ai-generated"`) are still preserved verbatim.
- **L7 — `generate_report_html` is pure.** The renderer used to mutate
  `report_data["variants"][*]` / `report_data["pgx_results"]` via
  `setdefault`. It now deep-copies those lists into a render-scoped
  view model so the caller's dict is never touched. Re-running the
  same report through the renderer is idempotent.
- **M3 — Real `skip_api` in batch provenance.** Pre-v2.5.4 batch wrote
  `pipeline.skip_api = False` and `db_versions = get_all_db_versions(skip_api=False)`
  regardless of the `--skip-api` CLI flag. The per-sample `report_data`
  now reflects the caller's actual value; research reproducibility is
  no longer dependent on the runner remembering to ignore these fields.
- **M4 — `sample_id` consistency across modes.** Single mode used
  `Path(vcf_path).stem.upper()` (gives `SAMPLE.VCF` for
  `sample.vcf.gz` — the inner `.vcf` survived and everything was
  uppercased). Batch mode stripped extensions properly but kept case.
  The two disagreed on the same file. The new `normalize_sample_id`
  helper — used by both modes — strips `.vcf.bgz` / `.vcf.gz` / `.vcf`
  in that order and preserves case. **Breaking**: single-mode sample
  id changes from `SAMPLE.VCF` to `sample`. Update downstream consumers
  that parse sample id from the report filename.
- **M5 — `.vcf.bgz` batch discovery.** Directory-mode discovery now
  globs `.vcf` + `.vcf.gz` + `.vcf.bgz`. The CLI + VCF parser
  already supported `.vcf.bgz`; batch discovery silently skipped them
  before.
- **M6 — Per-sample feature parity in batch.** `BatchSample` dataclass
  carries optional `germline_vcf`, `ped_path`, `hpo_ids`, `sv_path`,
  `intervar_path`, `clinical_note`, and `panel_size_mb`. Manifest CSV
  can opt into these columns; the legacy two-column
  (`sample_id,vcf_path`) manifest keeps working. HPO IDs merge (batch
  default + per-sample override). When `--clinical-board` is enabled
  batch assembly falls back to serial because the LLM layer is not
  thread-safe.
- **Side effect — `content_status="curated-civic"` frequency rises.**
  Because CIViC now wins over gene_knowledge on cancer reports, the
  curated-civic watermark appears on more rows than before. Downstream
  filters that tracked this value should re-baseline. (Annotated-VCF
  smoke: 0 → 6 rows carrying `curated-civic` on the demo fixture.)
- **Non-goals held.** Classification (ACMG / AMP tiering), ClinVar /
  gnomAD / KOVA query logic, Clinical Board prompts / LLM choice, and
  the Jinja template engine are untouched.
- Rollback tag: `pre-batch-refactor`.

### Changed — v2.5 (KOVA v7 Korean frequency migration)
- **Korean frequency source consolidated to KOVA v7.** The earlier
  multi-cohort Korean-frequency strategy (several parallel Korean-only
  cohorts) is replaced by a single KOGO / gene2korea Korean Variant
  Archive (KOVA) v7 source. KOVA v7 ships 43 M variants with per-variant
  AF **and homozygote counts**; the homozygote count is the feature
  that drives BS2 candidate flagging for autosomal-recessive review and
  is the main clinical advantage over gnomAD for Korean patients.
- `scripts/population/query_kova.py` is the single query module for
  the Korean cohort. `compare_freq.py` is rewritten to a 3-tier
  comparison (KOVA v7 / gnomAD EAS / gnomAD ALL) and surfaces KOVA
  homozygote counts in the `FrequencyData` payload for downstream BS2
  evaluation.
- `config.yaml` registers `paths.kova` (default
  `data/kova_freq.tsv`); the earlier per-cohort path entries are
  removed. `docs/KOREAN_STRATEGY.md`, `docs/SETUP.md`,
  `docs/TIERING_PRINCIPLES.md`, and `docs/ARCHITECTURE.md` are updated
  in lock-step.
- Historical design specs under `docs/superpowers/specs/` that refer
  to the prior multi-cohort plan carry a "superseded by KOVA v7
  (2026-04-21)" banner at the top. The specs themselves are
  preserved for audit continuity — they are no longer the active
  implementation target.
- User-facing docs (README.md / README.en.md / showcase copy) are
  rewritten around the 3-tier KOVA v7 comparison and call out the
  homozygote-count advantage explicitly.

### Added — v2.4 (germline integration & curator expansion)
- `--germline <vcf>` CLI flag enables somatic + germline dual-input
  pipeline. In rare-disease mode `extract_inherited_variants` intersects
  the germline VCF against a ClinVar P/LP **point-level** target BED
  (5000-variant cap) and yields inherited-variant rows that flow through
  the same ACMG path as the proband variants. Inherited variants gene-key
  off the BED's 4th column with a fallback to the VCF annotation when
  the column is empty.
- PharmCAT 3.2.0 integration via `scripts/pharma/pharmcat_runner.py`.
  Java 17 is auto-installed on first invocation if missing; the PharmCAT
  JAR is auto-downloaded and cached. The runner pre-filters the germline
  VCF with `tabix` against ~1000 PharmCAT-relevant positions to stay
  inside the 120 s harness ceiling. `pgx_source` in report metadata
  distinguishes `"pharmcat"` from `"builtin"` fallback.
- 12-gene PGx gap fill: DPYD, TPMT, HLA-A, CYP2B6, CYP4F2, ABCG2, NAT2,
  CACNA1S, CFTR, CYP3A4, MT-RNR1, RYR1 added to
  `data/korean_pgx_table.json` (24 genes total). Each entry now carries
  a `default_phenotype` field; the prior hardcoded `elif` chain in
  `scripts/pharma/korean_pgx.py` was replaced by data-driven lookup.
- `--ped <ped>` CLI flag with strict-mode trio / quartet resolution.
  Any unresolved role raises immediately. Filename heuristic
  (`*_proband.vcf`, `*_father.vcf`, ...) remains the fallback when
  `--ped` is omitted.
- De novo evidence codes in rare-disease mode: PS2 (confirmed) and PM6
  (assumed) fire when `parse_vcf` reads trio FORMAT/GT and detects a
  de-novo proband call. DDG2P neurodevelopmental panel ingest (2,201
  admitted genes) gives a separate selector carve-out. SpliceAI
  `delta_max >= 0.2` rescue applies to splice-region and synonymous
  variants in rare-disease mode.
- gnomAD v4.1 constraint OR-branch activated for de-novo admission.
- ClinGen ingest now parses the public CSV export by **header name**
  rather than column index, so future column reorders no longer
  silently mis-load fields.
- Literature Analyst (Rare-Disease Board) is now grounded against the
  local CIViC build (4,811 evidence_statements). PMIDs that cannot be
  matched against a known CIViC evidence statement are dropped from
  emitted prose.
- Two-model hybrid Clinical Board: domain agents run on
  **MedGemma 27B**, Board Chair runs on **SuperGemma4 31B**. Hybrid
  merge preserves the Chair narrative when the deterministic curator
  and the LLM Chair disagree on row count.
- De-novo badges on the rare-disease variant card, driven by
  `selection_reason_list` written back to `report_data` after Board
  selection.
- Clinical-priority sort applied to the report after Board selection:
  `(board_admitted, classification_rank, hpo_score, has_gene, variant_id)`
  — Board-admitted variants float to the top regardless of classification
  rank. Board-admitted VUS are auto-promoted from the
  `detailed_variants` block into the headline variant table.

### Changed — v2.4
- Default Clinical Board model wiring split: agent calls go to
  `alibayram/medgemma:27b`, Chair calls go to the SuperGemma4 31B
  endpoint.
- Variant-selector consequence gate (v2.2 B1) now also applies to the
  rare-disease selector path. P/LP variants and SpliceAI ≥ 0.2 splice
  rescues bypass the gate.
- README.md / README.en.md / CLAUDE.md document the v2.3 / v2.4
  features, system-tool dependencies (Java 17, bcftools, tabix,
  pysam), and the `--germline` / `--ped` usage examples.

### Tests — v2.4
- 988 → 1053+ tests across the v2.3 / v2.4 series. Trio FORMAT/GT
  parsing, PED resolution, PharmCAT runner, inherited-variant
  extraction, DDG2P carve-out, and de-novo PS2/PM6 each carry their
  own test module.

## [2.2.1] — 2026-04-15

Privacy cleanup and public-release baseline. No feature changes —
purely the hygiene work required before the repository transitioned
from private to public on 2026-04-15.

### Added
- Bilingual README (Korean primary `README.md` + English alternate
  `README.en.md`) with a language switcher at the top of both files.
- Mermaid pipeline flow diagram with per-mode Clinical Board specialist
  rosters (Cancer: Therapeutic Target Analyst / Tumor Genomics / PGx /
  Clinical Evidence; Rare Disease: Variant Pathologist / Disease
  Geneticist / PGx / Literature Analyst).
- GitHub Pages landing at `docs/index.html`, published at
  <https://junehawk.github.io/BIKO-GenomeBoard/>.
- Mandatory pre-push `ruff check` + `ruff format` checklist added to
  `CLAUDE.md` to prevent CI format-gate drift.

### Changed
- Dark-mode-friendly diagram colors (saturated fills + explicit white
  text) so Mermaid boxes remain legible in both themes.
- `template_renderer_chair` fallback headline rewritten from
  `"{N} curated therapy option(s) — LLM synthesis unavailable"` to
  `"No variant-specific treatment recommendations — evidence library
  only"`, with an explicit disclaimer body.
- Showcase reports regenerated end-to-end: `sample_cancer_report.{html,json}`
  and `sample_rare_disease_report.{html,json}`.

### Fixed
- 404 MedGemma links replaced with verified `deepmind.google` and
  `ai.google.dev` URLs.
- Legacy `gnomad.sqlite3` fallback removed from the query path (was
  emitting ~200 `WARNING` lines per cancer run).
- Consequence-form normalization (`fix(v22)`, commit `564f5da`):
  `evidence_collector` and `variant_selector` accept both raw VEP SO
  terms and BIKO-formatted short labels, restoring A3 PM1 hotspot
  lookup and B1 consequence gate on real pipeline data.
- CI green: ruff format drift eliminated across 127 files; CI-only test
  failures made DB-state-independent; CI passes on ubuntu-latest /
  Python 3.10, 3.11, 3.12.

### Removed
- `docs/superpowers/plans/` — contained plan material with real
  patient name, DOB, specimen ID, hospital, and ordering physician.
- `docs/API_KEYS.md` — deprecated external-API key documentation from
  a pre-v2.2 architecture BIKO no longer uses (Ollama + local DBs
  only).
- `codegen-Tumor_WB.*` real-patient VCFs and the derived
  `sample_codegen_777_report.{html,json}` showcase.
- All `codegen` / "777 variants" / patient-specific references from
  `README.md`, the Pages landing, the project intro HTML, spec files,
  and tests that previously used the codegen VCF as input.

### Security
- PHI removed from the tracked tree at HEAD. **History is not
  rewritten**; earlier commits still contain the codegen content. A
  full `filter-repo` history rewrite is tracked as a follow-up. This
  tag marks HEAD-clean, not history-clean.

### Tests
- **901 passed, 1 xfailed** (fixture #5 PHI non-persistence, deferred
  to v2.3+).

## [2.2.0] — 2026-04-14

Clinical-credibility release. Closes four patient-safety gaps from the
FoundationOne CDx head-to-head audit and refreshes the showcase folder
with regenerated sample reports.

### Added
- **A1/A2 — Curate-then-narrate architecture.** Treatment options are
  now deterministically curated from OncoKB and the local CIViC build.
  The Board Chair LLM may only narrate curated rows by
  `(curated_id, variant_key)` pair. `narrative_scrubber` strips any
  drug mention outside the curated set, and `template_renderer_chair`
  provides a deterministic Jinja fallback when the LLM cannot produce
  valid pairs — so no hallucinated drugs can reach the final report.
- **A3 — PM1 hotspot table.** `data/pm1_hotspot_domains.json` lets the
  ACMG engine fire PM1 via a curated gene/residue lookup
  (TP53/KRAS/NRAS/BRAF/PIK3CA/EGFR/IDH1-2) even when VEP `DOMAINS` is
  empty. Every entry cites a PMID.
- **A4 — ClinVar conflict reconciliation.** Narrow four-condition
  override (engine ≥ LP + ClinVar "Conflicting" + PM1 + PM5) emits
  `clinvar_override_reason` at the post-classify seam. Config-gated
  via `acmg.allow_engine_override_on_conflict`.
- **B1 — Variant selector consequence gate.** Tier III VUS admission
  rejects non-coding variants (intronic, UTR, upstream, downstream)
  with a P/LP bypass and SpliceAI ≥ 0.2 splice-rescue.
- **B2 — MMR/Lynch carve-out.** Protein-impacting VUS in
  `MLH1/MSH2/MSH6/PMS2/EPCAM` are admitted as `VUS_MMR_Lynch`
  regardless of hotspot status.
- `CLAUDE.md`, `docs/ARCHITECTURE.md`, `docs/TIERING_PRINCIPLES.md`
  v2.2 architecture notes.

### Changed
- README.md: new **v2.2 Highlights** section, Variant Selector section
  documents B1/B2, test count 673 → 899.
- Intro HTML (`BIKO_GenomeBoard_소개.html`) refreshed: overflow cards
  trimmed, new **Curate-then-Narrate (v2.2 A1/A2)** card, stats
  updated to the 899-test state.
- Showcase folder scoped to intro HTML + assets + three sample
  reports, all regenerated end-to-end with the current pipeline.

### Fixed
- Consequence-form normalization (`fix(v22)`, commit `564f5da`):
  `evidence_collector` and `variant_selector` now accept both raw VEP
  SO terms and the BIKO-formatted short labels (`"Missense"` etc.)
  produced by `parse_annotation.format_consequence`. Without this,
  A3/B1 were silent dead code on real data despite unit tests passing
  against raw-SO-term fixtures.
- gnomAD cleanup (`6712472`): dropped the legacy `gnomad.sqlite3`
  fallback from the query path. Tabix/pysam direct query is now the
  sole local gnomAD source.
- CI green (`eed3549`): ruff check + ruff format clean across 127
  files; four CI-only test failures made DB-state-independent
  (`civic.sqlite3` / `hpo.sqlite3` absence no longer breaks briefing
  and HPO fallback tests).

### Tests
- **899 passed, 1 xfailed** (`test_phase_a_smoke.py::fixture_5` PHI
  non-persistence, deferred to v2.3+).
- Pre-v2.2 baseline: 728 → post-v2.2: **+173 tests**.
- CI green on ubuntu-latest / Python 3.10, 3.11, 3.12.

## [pre-v2.2-phaseA] — 2026-04-14

Git tag marking the last v2.1 state before the Phase A curate-then-narrate
inversion. Preserved as a rollback checkpoint only; not a supported
release.

```bash
git reset --hard pre-v2.2-phaseA
```

[Unreleased]: https://github.com/junehawk/BIKO-GenomeBoard/compare/v2.2.1...HEAD
[2.2.1]: https://github.com/junehawk/BIKO-GenomeBoard/compare/v2.2.0...v2.2.1
[2.2.0]: https://github.com/junehawk/BIKO-GenomeBoard/compare/pre-v2.2-phaseA...v2.2.0
[pre-v2.2-phaseA]: https://github.com/junehawk/BIKO-GenomeBoard/releases/tag/pre-v2.2-phaseA
