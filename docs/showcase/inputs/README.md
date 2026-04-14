# Showcase inputs — multi-caller VCF staging

This directory is the staging location for the **input files** required to
regenerate the cancer showcase reports under `docs/showcase/`. The files
themselves are not committed (often patient-derived, large, or upstream-tool
outputs) — this README documents what the regen process expects to find here.

## Why a dedicated inputs directory

The cancer showcase needs more than a single SNV VCF — it needs the full
multi-caller picture (SNV + CNV + SV + ACMG evidence + panel BED). Scattering
those across `data/sample_vcf/`, `data/sample_sv/`, etc. makes regen scripts
brittle. Co-locating them under `docs/showcase/inputs/` keeps the showcase
self-contained and the regen command short.

## Expected files (cancer showcase: `codegen_777`)

| File | Caller / source | Required by | Notes |
|---|---|---|---|
| `codegen_777.mutect.passed.vep.vcf` | MuTect2 → VEP | `orchestrate.py` (positional `vcf_path`) | Somatic SNV/indel calls, PASS-filtered, VEP-annotated. May be replaced/augmented with Strelka. |
| `codegen_777.annotsv.tsv` | Canvas (CNV) + Manta (SV) → AnnotSV | `--sv` | AnnotSV merges Canvas CNV calls and Manta SV calls into a single annotated TSV. Without this, CNV/SV findings (e.g. CDKN2A/B loss, ERBB2 amplification) are not ingested. |
| `codegen_777.intervar.tsv` | InterVar | `--intervar` | Pre-computed ACMG evidence codes for germline interpretation. Optional but recommended for parity with the production pipeline. |
| `codegen_777.panel.bed` | exome/panel design | `--bed` | Defines the TMB denominator. If absent, fall back to `--panel-size <Mb>`. |

For a synthetic cancer demo (`sample_cancer`) the equivalent files live in
`data/sample_vcf/` and `data/sample_sv/` — see those directories for the
existing fixtures.

## Caller compatibility notes

- **SNV/indel** — `parse_vcf` accepts MuTect2, Strelka, or any VEP-annotated
  VCF. The current real-data fixtures under `data/sample_vcf/` use the naming
  pattern `<sample>.<caller>.passed.vep.vcf` (see
  `codegen-Tumor_WB.mutect.passed.vep.vcf` and
  `codegen-Tumor_WB.strelka.passed.vep.vcf`).
- **CNV** — Canvas calls are expected to be merged into AnnotSV upstream;
  `orchestrate.py` does not parse Canvas VCFs directly. The integration point
  is `scripts/intake/parse_annotsv.py`.
- **SV** — Manta calls follow the same pattern: merge into AnnotSV upstream.
- **TMB** — uses the SNV VCF + panel BED (or `--panel-size` Mb fallback).

## What does NOT belong here

- Any patient-identifying metadata (PHI). Showcase inputs must be either
  synthetic or fully de-identified.
- Raw FASTQ / BAM files — those are upstream of BIKO's scope.
- Output reports — those land in `docs/showcase/` (parent directory).

## Regeneration entry point

See `docs/showcase/README.md` § "Regenerating the showcase" for the exact
`orchestrate.py` invocation. A `make showcase-regen` Make target is planned
but not yet implemented; once present, it will read inputs from this directory
and write rendered HTML/JSON back to `docs/showcase/`.
