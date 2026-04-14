# CLAUDE.md — BIKO GenomeBoard agent guide

Architectural notes for AI coding agents working on this repo. See
`README.md` for the user-facing project description and
`.claude/skills/genomeboard-conventions` for code style conventions.

## Framing — research reference, not clinical instrument

BIKO produces a **research reference document**. It is **not** a clinical
decision instrument. Do not add CLIA-baseline requirements, attending
oncologist sign-off hooks, or "patient-ready" gates that assume the
output will be handed directly to a patient or physician. The report is
consumed by a reviewing researcher or clinician who makes their own
independent judgement.

Patient-safety guardrails (narrative scrubber, narrow classification
override, curator-only treatment rows) stay in place — they protect
against LLM hallucination regardless of downstream use. CLIA-specific
features like full patient header metadata persistence (v2.2 task B4)
and PMID citation superscripts (v2.2 task B3) are deferred to v2.3 and
beyond. The `reporting.persist_patient_metadata` config flag is the
forward-compat stopgap gate — keep it default-false.

## v2.2 Phase A — curate-then-narrate architecture

v2.2 inverted the Clinical Board data flow. Before v2.2, LLM agents
generated treatment recommendations directly from a free-form prompt,
then the output was post-filtered. Starting with Phase A (commits
`370edeb..388ec69`, rollback tag **`pre-v2.2-phaseA`**), every
treatment row must originate from a deterministic curator; the LLM
may only narrate curated rows.

### Curator (`scripts/clinical_board/curated_treatments.py`)

`curate_treatments(variants, offline_mode=...)` returns
`{variant_key: list[CuratedTreatment]}` by querying OncoKB's free-tier
variant-level API plus the local CIViC build. Each `CuratedTreatment`
carries `curated_id`, `drug`, `evidence_level`, `source`, `pmids`,
`disease_context`, `significance`, and the original `raw_row`.

`runner.run_clinical_board()` calls `curate_treatments` after the
selector runs, stores the result at `report_data["_curated_treatments"]`,
and renders the curated block into the Board Chair prompt as an
authoritative CURATED EVIDENCE section. The Chair is instructed to cite
only `curated_id` values from that block.

After synthesis, `narrative_scrubber.scrub_opinion()` validates every
emitted treatment row against the curated set by matching the
`(curated_id, variant_key)` **pair** — a bare curated_id is not
enough, which closes the EGFR→TP53 paste-attack path (Phase A fixture
2). Any row whose pair does not match is dropped. A banned-drug
denylist additionally strips hallucinated names out of every prose
field before rendering.

If the Chair emits zero treatment rows despite a non-empty curator
output, `template_renderer_chair.render_from_curated()` falls back to
a deterministic Jinja renderer over the curated rows. This guarantees
at least the curated rows reach the report, but never adds rows the
LLM did not endorse.

### ACMG PM1 hotspot domains (`data/pm1_hotspot_domains.json`)

`pm1_hotspot_domains.json` is the ACMG-authoritative PM1 source for
protein-domain hotspots. It is built by
`scripts/tools/build_pm1_hotspot_table.py` and registered through
`scripts/db/version_manager.py`. `evidence_collector.collect_additional_evidence`
consults this table (not the legacy `cancerhotspots_v2_single` build)
when firing PM1 for domain hotspots.

`cancerhotspots_v2_single` stays in the repo for **Tier III cancer
signalling only** — the variant selector's `is_hotspot()` check and
AMP tiering hotspot branch. Do not use it for ACMG PM1.

### Narrow ClinVar conflict override (`acmg_engine.py`)

`apply_hotspot_conflict_reconciliation` implements a four-condition
gate for overriding a ClinVar "conflicting interpretations" entry with
the engine's Likely Pathogenic verdict:

1. ClinVar review status is "conflicting".
2. The BIKO engine independently reaches LP or P.
3. PM1 fires from `pm1_hotspot_domains.json` (protein-domain hotspot).
4. PM5 fires (novel missense at a residue where another missense has
   an established ClinVar Pathogenic entry).

Only when all four hold does the override engage. The override reason
(including the PMID of the hotspot source and the PM5 peer residue) is
recorded in `ClassificationResult.clinvar_override_reason` and rendered
in the report as an `.override-notice` macro, so the reviewing
researcher sees *why* BIKO disagreed with ClinVar. This is the sole
circumstance in which BIKO overrides a curated database — the narrow
scope is deliberate.

## v2.2 Phase B — quality and production polish

Phase B is the production-polish pass on top of Phase A. Scope landed
on `main`:

- **B1 — selector consequence gate.** `scripts/clinical_board/variant_selector.py`
  applies a `_PROTEIN_IMPACTING_CONSEQUENCES` gate to `_cancer_must_reason`
  and `_cancer_may_reason`. Non-coding VUS (intronic, UTR, upstream,
  downstream, synonymous) are rejected from Tier I/II/III and the MAY
  arm. The `P_LP` must-reason branch bypasses the gate unconditionally
  so a deep-intronic ClinVar-Pathogenic splice variant still reaches
  the Clinical Board. Splice-region and synonymous variants are
  rescued at SpliceAI `delta_max >= 0.2` (Tavtigian et al. 2023
  PP3-moderate threshold), read from `v["in_silico"]["spliceai_max"]`.
- **B2 — MMR/Lynch panel carve-out.** `_MMR_LYNCH_GENES` admits any
  protein-impacting VUS in MLH1/MSH2/MSH6/PMS2/EPCAM with reason
  `VUS_MMR_Lynch` (priority 6, between `VUS_hotspot` and
  `VUS_TSG_LoF`), regardless of hotspot or TSG-LoF status. The B1
  consequence gate still applies inside the carve-out.

Deferred to v2.3 (long-term, per user 2026-04-14):

- **B3 — PMID references on AgentOpinion / CancerBoardOpinion.**
  Not needed for the research-reference use case; citation infrastructure
  will land when clinical-grade output is actually required.
- **B4 — patient header metadata.** `--patient-name`, `--patient-dob`,
  `--patient-mrn`, `--patient-ordering-physician`, `--patient-facility`,
  `--patient-specimen-id`, `--patient-collected-date`,
  `--patient-report-date` CLI flags are **not** wired in v2.2. The
  `reporting.persist_patient_metadata` config flag exists as a
  forward-compat gate for v2.3 — leave it default-false. Do not add
  PHI persistence paths without an explicit v2.3 ticket.

## Rollback

The rollback checkpoint for the entire v2.2 Phase A + B stack is git
tag **`pre-v2.2-phaseA`**. `git checkout pre-v2.2-phaseA` returns the
repo to the last v2.1 state.

## Pre-push checklist — run these BEFORE `git push`

GitHub Actions CI (`.github/workflows/ci.yml`) enforces **ruff
`0.15.10`** for lint and format, plus `pytest tests/ -k "not
integration"` across Python 3.10 / 3.11 / 3.12. New or edited files
that pass pytest locally can still fail CI on a single missed
`ruff format` pass — this has happened twice in the v2.2 release
history. Run the full sequence locally before every push:

```bash
# 1. Lint — autofix what ruff can, review the rest
ruff check --fix scripts/ tests/

# 2. Format — must match CI's `ruff format --check` exactly
ruff format scripts/ tests/

# 3. Verify both gates are green (what CI actually runs)
ruff check scripts/ tests/          # must say "All checks passed!"
ruff format --check scripts/ tests/  # must say "N files already formatted"

# 4. Full test suite (CI runs without the integration marker)
python -m pytest tests/ -q
```

Every edit that touches `scripts/` or `tests/` — including a small
one-line assertion change, a new regression test, or a generated
file — goes through this sequence. Committing a file without running
`ruff format` on it is the single most common cause of CI red on
this repo: unit tests pass locally, `ruff check` passes locally, but
`ruff format --check` fails on CI because the editor-saved formatting
differs from ruff's canonical output. The four-step sequence above
eliminates that failure mode.

Special notes:

- **Tests that depend on optional local DBs** (`civic.sqlite3`,
  `hpo.sqlite3`, etc.) must guard on actual query results, not on
  `os.path.exists`. CI runners do not provision these DBs and an
  empty/stub file can trip `os.path.exists` while queries return
  nothing — see the v2.2 fix in `tests/test_rare_disease.py` and
  `tests/test_case_briefing.py` for the pattern.
- **Consequence-form normalization** — the real pipeline stores
  `variant.consequence` as a BIKO-formatted short label
  (`"Missense"`, not `"missense_variant"`) via
  `scripts/intake/parse_annotation.py::format_consequence`. Any new
  code that filters on VEP SO terms must either call the
  `_canonicalize_so_term` / `_canonical_consequence` helpers or
  include both forms in its frozenset, or it will be silent dead
  code on real data while still passing unit tests that use raw SO
  terms in fixtures. This bit us in v2.2 A3/B1; see commit
  `564f5da`.

## Things to never do

- Do **not** put the classification or tiering engine behind an LLM
  call. Classification must stay deterministic Python.
- Do **not** let the Board Chair or any domain agent invent a drug or
  `curated_id` that is not in the CURATED EVIDENCE block — the
  narrative scrubber exists to catch this, but defence-in-depth means
  prompts must also instruct agents to cite only curated rows.
- Do **not** widen the ClinVar conflict override beyond the
  four-condition gate without a clinical signoff; it is the sole
  database-override path BIKO uses.
- Do **not** persist patient metadata by default; the
  `reporting.persist_patient_metadata` gate must stay off until v2.3
  B4 lands.
- Do **not** use `--no-verify` on commits or push directly to
  `origin/main` from agent workflows.
- Do **not** `git push` without running the Pre-push checklist
  (above) first. `ruff format` missed on a single edited file is
  the #1 cause of CI red on this repo — cheaper to run the local
  commands than to chase a CI failure notification.
