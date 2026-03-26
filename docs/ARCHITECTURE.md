# GenomeBoard Architecture

## Design Principle

> All variant classification and database queries run as deterministic Python scripts. No LLM is involved in classification logic.

This ensures reproducibility: the same VCF always produces the same ACMG codes, tier, and classification.

---

## Pipeline Overview

```
VCF Input
    │
    ▼
[1] VCF Parser (parse_vcf.py)
    - cyvcf2-based parsing
    - CSQ/ANN field extraction (VEP/SnpEff pre-annotated VCFs)
    - rsID extraction from Existing_variation
    - HGVS, consequence, transcript, SIFT/PolyPhen parsing
    │
    ▼ (parallel per variant, ThreadPoolExecutor)
[2] Database Queries
    ├── ClinVar      → query_local_clinvar.py / query_clinvar.py (API fallback)
    ├── gnomAD       → query_tabix_gnomad.py → query_local_gnomad.py → query_gnomad.py
    ├── KRGDB        → query_krgdb.py (local TSV)
    └── PGx          → korean_pgx.py (local lookup)
    │
    ▼
[3] Frequency Comparison (compare_freq.py)
    - KRGDB → gnomAD EAS → gnomAD ALL (3-tier Korean enrichment)
    - Produces BA1, BS1, PM2_Supporting ACMG codes
    │
    ▼
[4] ACMG Classification (acmg_engine.py)
    - Aggregates codes from ClinVar + frequency comparison
    - Applies ACMG/AMP 2015 rules from acmg_rules.json
    - PGx genes → "Drug Response"; APOE → "Risk Factor"
    - ClinVar conflict detection
    - ClinVar override (expert panel / multi-submitter)
    │
    ▼
[5] AMP/ASCO/CAP 2017 Tiering (amp_tiering.py)
    - Assigns AMP Tier I–IV to every variant (cancer mode)
    - Strategy A/B/C selectable via config.yaml somatic.tiering_strategy
    - CIViC variant-level evidence (get_predictive_evidence_for_tier)
    - OncoKB gene-level evidence (oncokb.py — gene data provider)
    - Hotspot detection via CIViC variant coordinates
    │
    ▼ (rare-disease mode only)
[5b] Rare Disease Enrichment
    - HPO phenotype matching (hpo_matcher.py)
    - OMIM gene-disease lookup (query_omim.py)
    - ClinGen validity (query_clingen.py)
    - Candidate ranking: classification rank, then HPO score
    │
    ▼
[6] Report Generation (generate_pdf.py)
    - Jinja2 templates → HTML
    - Optional WeasyPrint PDF conversion
    - VUS filtering when --hide-vus is set
```

---

## AMP/ASCO/CAP 2017 Somatic Tiering

Cancer mode uses the AMP/ASCO/CAP 2017 guidelines (Li MM et al., *J Mol Diagn* 2017). See [docs/TIERING_PRINCIPLES.md](TIERING_PRINCIPLES.md) for full design rationale.

`scripts/somatic/amp_tiering.py` is the tiering engine. It combines OncoKB gene-level evidence with CIViC variant-level evidence under **Modified Approach B** (default):

| Tier | Label | Criteria |
|------|-------|----------|
| I | Strong Clinical Significance | FDA-approved therapy / professional guideline biomarker; Drug Response; Risk Factor |
| II | Potential Clinical Significance | Clinical trial evidence; emerging clinical significance; VUS at known hotspot |
| III | Unknown Clinical Significance | VUS on cancer gene (non-hotspot) |
| IV | Benign or Likely Benign | Benign/Likely Benign; non-cancer gene VUS |

### Tiering Strategy (config.yaml)

```yaml
somatic:
  tiering_strategy: "B"  # "A" (CIViC priority), "B" (combined, default), "C" (OncoKB only)
```

| Strategy | Behaviour |
|----------|-----------|
| **A** | CIViC variant-level evidence takes priority; OncoKB as fallback |
| **B** | OncoKB + CIViC combined — CIViC can elevate tier, never lower (default) |
| **C** | OncoKB only — backward compatible with legacy `assign_tier()` |

### Priority Order (Strategy B)

1. CIViC variant-specific Level A + Pathogenic/LP → Tier I
2. CIViC variant-specific Level B → Tier II
3. OncoKB Level 1–2 + Pathogenic/LP → Tier I
4. OncoKB Level 1–2 + VUS + hotspot → Tier II
5. CIViC variant-specific Level C–D + Pathogenic/LP → Tier II
6. Pathogenic/LP on any cancer gene → Tier II
7. VUS on cancer gene → Tier III
8. All other → Tier IV

The OncoKB cancer gene list is stored in `data/oncokb_cancer_genes.json` (gene name → `{type, level}`).

`scripts/clinical/oncokb.py::assign_tier()` is a **deprecated wrapper** that delegates to `amp_assign_tier(strategy="C")` for backward compatibility.

---

## CIViC Integration

`scripts/db/query_civic.py` provides:

- **`get_gene_summary(gene)`** — gene description and aliases from CIViC
- **`get_variant_evidence(gene, variant_name)`** — treatment evidence items sorted by evidence level (A > B > C > D > E), with disease, drugs, clinical significance, and PMID citations
- **`get_predictive_evidence_for_tier(gene, hgvsp)`** — returns `{match_level, evidence}` for tiering; tries variant-specific match first (using `scripts/common/hgvs_utils.py` for HGVSp→CIViC name conversion), falls back to gene-level. Only Predictive evidence type is returned.
- **`is_hotspot(gene, protein_position)`** — checks whether a protein position falls within any CIViC variant coordinate range for that gene

The local SQLite database (`data/db/civic.sqlite3`) is built from CIViC's public TSV exports via `scripts/db/build_civic_db.py`. Tables: `genes`, `variants`, `evidence`.

HGVSp conversion utilities (e.g., `p.Val600Glu` → `V600E`) are shared via `scripts/common/hgvs_utils.py`.

---

## Hotspot Detection

Hotspot detection uses CIViC variant coordinates:

1. `extract_protein_position(hgvsp)` parses the amino-acid position from an HGVSp string (e.g., `p.Val600Glu` → `600`).
2. `is_hotspot(gene, position)` queries the CIViC `variants` table for any record with `start <= position <= stop` for that gene.
3. A VUS at a hotspot position is promoted to Tier II (`amp_assign_tier` in `scripts/somatic/amp_tiering.py`).

---

## ClinVar Override Logic

`apply_clinvar_override()` in `scripts/classification/acmg_engine.py` overrides the ACMG engine's classification when ClinVar provides high-confidence evidence:

| ClinVar review status | Pathogenic | Likely Pathogenic | Benign | Likely Benign |
|-----------------------|-----------|-------------------|--------|---------------|
| Expert panel / practice guideline | → Pathogenic | → Likely Pathogenic | → Benign | → Likely Benign |
| Multiple submitters, no conflicts | → Likely Pathogenic | → Likely Pathogenic | → Likely Benign | → Likely Benign |
| Single submitter | no override | no override | no override | no override |
| Conflicting submissions | no override | no override | no override | no override |

When an override occurs, `classification.clinvar_override = True` and `original_engine_classification` is preserved for display in the report.

`check_clinvar_conflict()` independently flags cases where the engine and ClinVar differ by 2+ steps on the 5-tier scale (e.g., engine says VUS, ClinVar says Pathogenic).

---

## ACMG Classification Engine

`scripts/classification/acmg_engine.py` is a pure deterministic rule engine:

1. Evidence codes from ClinVar and frequency comparison are collected as `AcmgEvidence` objects.
2. `_count_by_strength()` buckets codes into PVS/PS/PM/PP/BA/BS/BP strength tiers. `_Supporting` suffixes downgrade by one tier.
3. Rules from `data/acmg_rules.json` define minimum counts per strength tier for each classification (Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign).
4. Conflicting evidence (both pathogenic and benign codes meeting their respective rules) → VUS with `conflict=True`.
5. PGx genes and APOE bypass the rule engine and receive fixed classifications.

---

## Batch Pipeline Flow

```
[1] Discover samples (directory glob or manifest CSV)
        │
        ▼
[2] Parse all VCFs → collect unique variants
    (variant key = chrom:pos:ref:alt)
        │
        ▼
[3] Annotate unique variants in parallel (ThreadPoolExecutor, up to 10 workers)
    Each unique variant annotated once regardless of how many samples share it
        │
        ▼
[4] HPO resolution (rare-disease mode only)
        │
        ▼
[5] Per-sample assembly (reuse cached annotations)
    → ACMG classification
    → OncoKB tiering
    → Rare disease enrichment (if applicable)
    → Build report_data dict
        │
        ▼
[6] Parallel report generation (ProcessPoolExecutor)
    One HTML report per sample
```

Variant deduplication is the key performance optimization: a cohort of 100 samples sharing common variants annotates each variant once rather than 100 times.

---

## Data Source Priority

`annotation.source` in `config.yaml` controls the fallback chain:

| Mode | ClinVar | gnomAD |
|------|---------|--------|
| `auto` | Local SQLite → ClinVar E-utilities API | tabix VCF → local SQLite → GraphQL API |
| `local` | Local SQLite only | tabix VCF → local SQLite |
| `api` | ClinVar E-utilities API | gnomAD GraphQL API |
| `--skip-api` | Local SQLite only (forced) | tabix VCF → local SQLite (forced) |

---

## File Structure

```
gb/
├── scripts/
│   ├── orchestrate.py              # CLI + pipeline orchestration
│   ├── classification/
│   │   └── acmg_engine.py          # ACMG/AMP 2015 rule engine
│   ├── somatic/
│   │   └── amp_tiering.py          # AMP/ASCO/CAP 2017 tiering engine (TierResult, strategy A/B/C)
│   ├── clinical/
│   │   ├── oncokb.py               # OncoKB gene data provider; assign_tier() deprecated wrapper
│   │   ├── query_clinvar.py        # ClinVar E-utilities API
│   │   ├── hpo_matcher.py          # HPO phenotype scoring
│   │   ├── query_omim.py           # OMIM gene-disease
│   │   └── query_clingen.py        # ClinGen gene validity
│   ├── db/
│   │   ├── build_clinvar_db.py     # ClinVar SQLite builder
│   │   ├── build_civic_db.py       # CIViC SQLite builder
│   │   ├── build_gnomad_db.py      # gnomAD SQLite builder (legacy)
│   │   ├── query_local_clinvar.py  # ClinVar SQLite queries
│   │   ├── query_tabix_gnomad.py   # gnomAD tabix queries (primary)
│   │   ├── query_local_gnomad.py   # gnomAD SQLite queries (fallback)
│   │   ├── query_civic.py          # CIViC evidence + hotspot queries
│   │   └── version_manager.py      # DB version metadata
│   ├── korean_pop/
│   │   ├── query_krgdb.py          # KRGDB local TSV
│   │   ├── query_gnomad.py         # gnomAD GraphQL API
│   │   └── compare_freq.py         # Frequency comparison + ACMG codes
│   ├── pharma/
│   │   └── korean_pgx.py           # PGx gene patterns + Korean prevalence
│   ├── counselor/
│   │   └── generate_pdf.py         # Jinja2 HTML + WeasyPrint PDF
│   └── common/
│       ├── models.py               # Variant, AcmgEvidence, FrequencyData
│       ├── cache.py                # SQLite response cache
│       ├── config.py              # config.yaml loader
│       └── hgvs_utils.py          # HGVSp↔CIViC variant name conversion (shared)
├── data/
│   ├── krgdb_freq.tsv              # Korean population frequencies
│   ├── oncokb_cancer_genes.json    # OncoKB gene list
│   ├── acmg_rules.json             # ACMG classification rules
│   ├── gene_knowledge.json         # Gene metadata
│   └── db/                         # Built databases (gitignored)
├── templates/                      # Jinja2 report templates
├── tests/                          # pytest suite
├── config.yaml
├── Dockerfile
└── docker-compose.yml
```
