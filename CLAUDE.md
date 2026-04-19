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

## v2.3 / v2.4 — germline integration & curator expansion

v2.3 and v2.4 extended BIKO from a somatic-VCF tool into a dual-input
pipeline (somatic + germline) and filled the deterministic-curator
side of the rare-disease and PGx paths. The Phase A inversion from
v2.2 still holds: every treatment / drug row originates from a curator,
never from free-form LLM generation.

### --germline VCF + inherited variant extraction

`orchestrate.py` accepts an optional `--germline <vcf>` alongside the
primary somatic / proband VCF. In `--mode rare-disease`,
`scripts/pipeline/extract_germline.py::extract_inherited_variants`
intersects the germline VCF against a target BED and yields
inherited-variant rows that flow through the same ACMG path as the
proband variants. The somatic VCF still drives Tier I–IV (cancer mode)
or proband classification (rare-disease mode); the germline VCF
populates the inherited block and the PGx path. When `--germline` is
omitted in cancer or rare-disease mode the pipeline emits a single
warning explaining that PGx accuracy degrades without a germline call
set, then continues with the somatic input only.

The target BED for germline extraction is **point-level** (one row per
ClinVar P/LP position) with a 5000-variant cap. Earlier gene-span
BEDs ran tabix over hundreds of kilobases per gene, which made the
extraction run for tens of minutes on whole-genome germline VCFs.
Do not revert to gene-span BEDs — see commit `9c5e9bb` for the
optimisation history.

The extracted inherited variants carry their gene assignment from the
BED's 4th column (gene symbol). When the BED 4th column is empty the
extractor falls back to the VCF's annotation gene; do not remove this
fallback, it is what allows older annotated germline VCFs (no gene in
INFO) to still produce gene-keyed rows in the report (commit `65c6057`).

### PharmCAT 3.2.0 integration

`scripts/pharma/pharmcat_runner.py` invokes Oracle PharmCAT 3.2.0 over
the germline VCF for the curator-grade PGx path. Setup is automated:
the runner installs Java 17 on first use if `java -version` reports
< 17, and downloads the PharmCAT JAR on first run. `pgx_source` in
the report metadata distinguishes `"pharmcat"` (full PharmCAT pipeline
ran) from `"builtin"` (fallback to the in-repo `korean_pgx_table.json`
when PharmCAT is unavailable or the germline VCF is missing).

PharmCAT requires VCF positions to be limited or it walks the entire
genome and times out at the harness's 120 s ceiling. Before invoking
PharmCAT the runner pre-filters the germline VCF with `tabix` against
~1000 PharmCAT-relevant positions. **Do not remove the pre-filter** —
without it the runner times out on whole-genome germline VCFs and the
PGx block silently falls back to builtin. See `9c5e9bb`.

When PharmCAT runs, the runner v3 parser merges the PharmCAT phenotype
output with the builtin metadata table so Korean-vs-Western prevalence
and clinical-impact text from `korean_pgx_table.json` survive in the
report alongside the diplotype call.

### 12-gene PGx gap fill

`data/korean_pgx_table.json` was extended in v2.4 to cover the
PharmCAT gene set: DPYD, TPMT, HLA-A, CYP2B6, CYP4F2, ABCG2, NAT2,
CACNA1S, CFTR, CYP3A4, MT-RNR1, RYR1 — bringing the total to 24
genes. Each entry now carries a `default_phenotype` field used when
PharmCAT is unavailable. The previous hardcoded `elif` chain in
`scripts/pharma/korean_pgx.py` was replaced by data-driven lookup
against this field (commit `decbec3`); new genes are added by editing
JSON, not Python.

### --ped CLI + strict-mode trio / quartet

`scripts/intake/parse_ped.py` resolves trio (proband + 2 parents) and
quartet (+ 1 sibling) topologies from a standard PED file. The
`--ped` CLI flag is **strict mode**: when `--ped` is supplied the
parser raises if any role cannot be resolved unambiguously from the
file. When `--ped` is omitted, the parser falls back to filename
heuristics (`*_proband.vcf`, `*_father.vcf`, ...) for backward
compatibility with the pre-v2.4 sample layout.

Do **not** remove the strict mode and route everything through the
filename fallback — strict mode is the only code path that gives a
loud failure when a PED file disagrees with the trio expectation.
See `44027af` for the rationale.

### De novo evidence codes (PS2 / PM6 / DDG2P / SpliceAI)

In rare-disease mode with a resolved trio, `parse_vcf` reads the
trio FORMAT/GT cells and computes a de-novo flag per variant. The
ACMG engine then fires:

- **PS2** — confirmed de novo (both parents 0/0, proband het) in a
  gene with a maternally-or-paternally-confirmed disease association
  and an established gene–disease link.
- **PM6** — assumed de novo (parental phenotype absent but parental
  genotype not biologically confirmed) in the same gene set.
- **DDG2P carve-out** — variants in genes on the DDG2P
  neurodevelopmental panel (2,201 admitted genes, ingested in
  `3f34ec0`) are admitted to the rare-disease selector even when
  proband-only evidence would otherwise leave them at VUS, provided
  de novo (PS2/PM6) fires or the gnomAD v4.1 constraint OR-branch is
  satisfied (`693cc23`).
- **SpliceAI rescue** — splice-region and synonymous variants are
  rescued at `delta_max >= 0.2` (Tavtigian et al. 2023 PP3-moderate
  threshold) from `v["in_silico"]["spliceai_max"]`. This is the same
  rescue introduced in v2.2 B1 but now also applies to the rare-disease
  selector path.

The report renders these as **de-novo badges** on the variant card
(`4592758`). The badges are driven by `selection_reason_list`, which
the Board selector writes back to `report_data` after admitting a
variant (`bb9433e`); without this write-back the report would show
zero badges even when the selector admitted a de-novo variant.

### Clinical priority sort + board-admitted VUS promotion

After the Clinical Board selection completes, `report_data` is
re-sorted by clinical priority before render:

```
(board_admitted, classification_rank, hpo_score, has_gene, variant_id)
```

`board_admitted` is the primary key so any variant the Board admits
floats to the top of the report regardless of classification rank
(commit `a0240c9`). Board-admitted VUS are also auto-promoted from
the `detailed_variants` block into the headline variant table
(`65c6057`) — without the promotion the report would surface a
Board recommendation for a variant the reader could not find on
the first page.

### MedGemma + SuperGemma4 hybrid Board config

The Clinical Board now runs a **two-model hybrid**: domain agents
(Therapeutic Target Analyst, Tumor Genomics Specialist, PGx
Specialist, Clinical Evidence Analyst, Variant Pathologist, Disease
Geneticist, Literature Analyst) all run on **MedGemma 27B**; the
Board Chair runs on **SuperGemma4 31B**. The split was introduced in
`7290586` because SuperGemma4 produced more coherent narrative
synthesis than MedGemma at the Chair seam, while MedGemma kept its
domain-grounding edge for the agent layer.

The hybrid merge step (`0a21386`) preserves the Chair narrative when
the deterministic curator and the LLM Chair disagree on row count —
the Chair narrative is kept, but the curator rows backstop missing
treatment lines. Do not let merge logic drop the Chair narrative when
falling back to curated rows; the deterministic Jinja renderer
(`template_renderer_chair.render_from_curated()`) is only the last
resort, used when the Chair emits zero rows.

### Literature Analyst CIViC grounding

The Rare Disease Board's Literature Analyst is now grounded against
the local CIViC build (4,811 evidence_statements). Free-form prose
that names a PMID is checked against CIViC's evidence corpus before
emission (`bb8ba8d`). Literature claims that cannot be matched to a
known evidence statement are dropped, mirroring the curate-then-narrate
contract that Phase A established for treatments.

### ClinGen column-map CSV parser

The ClinGen ingest now parses the public CSV export by **header
name** rather than column index (`eb4c760`). When ClinGen reorders
columns in a future export, the build no longer silently mis-loads
fields — the loader looks up `MOI`, `Disease`, `GENE_SYMBOL`,
etc. from the header. `setup_databases.sh` was updated to pass the
CSV path to the new parser; the legacy positional parser is gone.

### Patient-metadata flags (still deferred)

v2.4 does **not** wire `--patient-name`, `--patient-dob`, `--patient-mrn`,
etc. The `reporting.persist_patient_metadata` config flag remains
default-false. The `--ped` flag is the closest v2.4 came to patient
metadata, and it deliberately stays scoped to PED relationship
resolution — it does not introduce a PHI-persistence path.

## Rollback

The rollback checkpoint for the entire v2.2 Phase A + B stack is git
tag **`pre-v2.2-phaseA`**. `git checkout pre-v2.2-phaseA` returns the
repo to the last v2.1 state.

The v2.4 work landed incrementally on `main` between commits `8ec13eb`
(2026-04-15, PharmCAT integration) and `a0240c9` (2026-04-19,
clinical priority sort). There is no consolidated `pre-v2.4` tag
yet — for a clean rollback, target `pre-v2.2-phaseA` and replay the
desired v2.3/v2.4 commits selectively.

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
- Do **not** remove the PharmCAT pre-filter (tabix-narrowed germline
  VCF) — without it, PharmCAT walks the whole genome and trips the
  120 s timeout, which silently degrades the PGx block to the builtin
  fallback. See commit `9c5e9bb`.
- Do **not** revert the ClinVar P/LP target BED back to gene-span
  ranges. The point-level (one row per P/LP position) BED with the
  5000-variant cap is what makes whole-genome germline extraction
  finish in seconds instead of tens of minutes.
- Do **not** drop the `selection_reason_list` write-back from the
  Board selector to `report_data`. Without it, the de-novo /
  promoted-VUS badges in the rare-disease report render empty even
  when the selector admitted a variant for those reasons (commit
  `bb9433e`).
- Do **not** route the trio-resolution path solely through the
  filename heuristic and remove `--ped` strict mode. Strict mode is
  the only path that fails loudly when a PED file disagrees with the
  trio expectation; the heuristic is the fallback for legacy sample
  layouts only.
