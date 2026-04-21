# Changelog

All notable changes to BIKO GenomeBoard are documented here.

This project adheres to [Keep a Changelog 1.1](https://keepachangelog.com/en/1.1.0/)
and uses [Semantic Versioning](https://semver.org/).

BIKO GenomeBoard is a **research reference tool**, not a clinical
decision instrument. Every release carries the same disclaimer: output
is intended for independent review by a researcher or clinician.

## [Unreleased]

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
