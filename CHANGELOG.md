# Changelog

All notable changes to BIKO GenomeBoard are documented here.

This project adheres to [Keep a Changelog 1.1](https://keepachangelog.com/en/1.1.0/)
and uses [Semantic Versioning](https://semver.org/).

BIKO GenomeBoard is a **research reference tool**, not a clinical
decision instrument. Every release carries the same disclaimer: output
is intended for independent review by a researcher or clinician.

## [Unreleased]

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
