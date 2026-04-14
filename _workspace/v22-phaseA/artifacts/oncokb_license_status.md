# OncoKB Licensing Status — v2.2 Phase A

**Date:** 2026-04-14
**Author:** db-dev
**Scope:** A1 curated treatments runner (`scripts/clinical_board/curated_treatments.py`)

## 1. Default posture

The v2.2 `curated_treatments` module uses the **OncoKB public variant-annotation
endpoint** as its default upstream:

```
POST {oncokb_base_url}/annotate/mutations/byProteinChange
```

- Base URL configured under `clinical_board.curated_treatments.oncokb_base_url`
  (default `https://www.oncokb.org/api/v1`).
- HTTP timeout: `clinical_board.curated_treatments.http_timeout_s` (5 s default).
- Responses are cached in the shared `variant_cache.sqlite3` via
  `set_cached_ns("oncokb", key, value)` — **no parallel `data/cache/oncokb/` tree**.
- 7-day TTL inherited from `cache.ttl_seconds`; single concurrency-safe
  SQLite path for the whole pipeline.

## 2. Licensing posture for a Korean clinical pilot

OncoKB's public API is gated by:

1. A registered API token (`https://www.oncokb.org/account/register`), even for
   the variant-annotation endpoint.
2. Terms of Use that restrict **commercial and clinical redistribution** of
   OncoKB-derived content.
3. An explicit academic-license agreement for non-academic institutional use.

**Risk:** silently embedding the free-tier API into a Korean clinical pipeline
is a licensing landmine. For the pilot, clinical deployment must obtain a
formal academic or commercial OncoKB license **before go-live**.

## 3. Degradation path — `offline_mode`

Toggle `clinical_board.curated_treatments.offline_mode: true` to skip OncoKB
entirely. When set:

- `curated_treatments.curate(...)` does not hit the OncoKB API.
- The curator returns CIViC-only rows. Merge-key fallback is string-based on
  `drug name`.
- The rendered Treatment Options table must display a
  **"curation source: CIViC-only (OncoKB offline)"** footnote (owned by
  report-dev in A1-db-3).

This provides a clean kill-switch for air-gapped deployments and for any
period during which the licensing question is still open.

## 4. Merge-key strategy

To minimise drug-name string fragility the CIViC `evidence` table now carries
a `therapy_ids` column (A1-db-4, this task). Merge order is:

1. **CIViC `therapy_ids` ↔ OncoKB therapy identifier** — deterministic.
2. **Normalised drug name** — fallback when `therapy_ids` is empty on either
   side. Case-folded, brand/salt suffix stripped.

Pure string fallback is retained only as a safety net; curator rows should
always prefer identifier-based merges.

## 5. Open items tracked for v2.2.1

- Cache-warm preflight script for air-gapped pilots
  (`scripts/tools/warm_variant_cache.py`) — not yet written.
- Formal OncoKB license agreement signed with the pilot institution.
- Add `data/drug_synonyms.json` if the `therapy_ids` merge has gaps in
  real-world CIViC/OncoKB pairs.

— db-dev
