# BIKO GenomeBoard

Korean Population-Aware Genomic Variant Interpretation Platform

> **Research reference only** — BIKO GenomeBoard produces research reference documents for clinicians to consult. It is not a clinical decision instrument and is not intended to make treatment recommendations.

---

## v2.2 Highlights (2026-04-14)

Patient-safety hot-fix release closing four clinical-credibility gaps surfaced by the FoundationOne CDx head-to-head audit:

- **Curate-then-narrate (A1/A2)** — Treatment Options are now deterministically curated from OncoKB + local CIViC. The Board Chair LLM may only *narrate* curated rows. A `narrative_scrubber` drops any drug mention outside the curated set, a `template_renderer_chair` deterministic fallback fires when the LLM cannot produce valid `(curated_id, variant_key)` pairs, and futibatinib-as-KRAS-inhibitor-style hallucinations are structurally impossible.
- **PM1 hotspot table (A3)** — `data/pm1_hotspot_domains.json` lets the ACMG engine fire PM1 via a curated gene/residue lookup (TP53, KRAS, NRAS, BRAF, PIK3CA, EGFR, IDH1/2, …) even when VEP `DOMAINS` is empty. Every entry cites a PMID.
- **ClinVar conflict reconciliation (A4)** — `apply_hotspot_conflict_reconciliation` at the post-classify seam emits `clinvar_override_reason` when the engine independently reaches LP, ClinVar is "Conflicting", and PM1 + PM5 both fire. Narrow, auditable, config-gated.
- **Selector tightening (B1) + MMR/Lynch carve-out (B2)** — Tier III VUS selector now enforces a protein-impacting consequence gate (with P/LP bypass and SpliceAI ≥ 0.2 splice-rescue) and admits all protein-impacting MMR-gene VUS regardless of hotspot status.

Rollback: `git reset --hard pre-v2.2-phaseA` (tag at the pre-Phase-A checkpoint).

Full v2.2 plan: [`docs/superpowers/plans/2026-04-14-ai-board-v2.2.md`](docs/superpowers/plans/2026-04-14-ai-board-v2.2.md)

---

## Features

### Variant Classification
- **ACMG/AMP 2015** deterministic rule engine (germline)
- **AMP/ASCO/CAP 2017** somatic tiering (Tier 1–4) with OncoKB + CIViC
- **ClinGen SVI 2022** in silico thresholds (REVEL, CADD, AlphaMissense, SpliceAI → PP3/BP4)
- Additional ACMG evidence collection (PVS1, PM1, PM4, PM5, PP2, BP1, BP7)
- InterVar integration for external ACMG evidence
- ClinVar expert panel / multi-submitter direct override

### Korean Population Analysis
- **5-tier frequency comparison**: KRGDB, Korea4K, NARD2, gnomAD EAS, gnomAD ALL
- 12-gene PGx screening with Korean vs Western prevalence comparison
- CPIC Level A/B drug-gene interaction recommendations

### AI Clinical Board (v2.2)
- Local multi-agent diagnostic synthesis powered by **MedGemma 27B** (Ollama)
- **Mode-specific agent sets** — rare-disease (Variant Pathologist, Disease Geneticist, PGx Specialist, Literature Analyst) and cancer (Therapeutic Target Analyst, Tumor Genomics Specialist, PGx Specialist, Clinical Evidence Analyst)
- **Grounded prompting** — per-agent domain sheets inject ClinVar / OncoKB / CIViC / HPO / OMIM / CPIC context instead of free-form reasoning
- **Curate-then-narrate (v2.2 A1/A2)** — treatment options come from a deterministic OncoKB + local CIViC curator (`scripts/clinical_board/curated_treatments.py`); the LLM may only narrate curated rows by `(curated_id, variant_key)`. A `narrative_scrubber` enforces that no drug outside the curated set appears in the final opinion, and a `template_renderer_chair` deterministic fallback fires when the LLM cannot produce valid curated pairs.
- **`CancerBoardOpinion`** treatment-focused output schema (therapeutic implications, treatment options with PMID, immunotherapy eligibility, monitoring plan)
- **Clinical note input** — `--clinical-note` / `--clinical-note-file` feeds free-text history (KR/EN, 1500 char cap) into briefings without altering deterministic classification
- **Knowledge Base (hybrid SQLite + Wiki)** — prior board decisions stored in `data/knowledge_base/kb.sqlite3`, cross-case variant stats injected as Prior Knowledge (with anti-anchoring guardrail)
- Bilingual output (English / Korean), temperature 0.1 for consistency
- Board Chair synthesis with consensus opinion, differential diagnoses (rare) or therapeutic strategy (cancer)
- Deterministic classification is never altered by AI — AI provides interpretive synthesis only

### Variant Selector (AI Board Input Filter)
Clinical tiered filter that decides which variants reach the AI Board, implemented per AMP/ASCO/CAP 2017 (PMID 27993330) and ACMG/AMP 2015 (PMID 25741868).

- **MUST-include (cancer)** — ClinVar P/LP, AMP Tier I/II, and Tier III variants at Cancer Hotspots residues or OncoKB oncogenic annotations in COSMIC CGC Tier 1 genes
- **MUST-include (rare disease)** — ClinVar P/LP, ACMG P/LP, and CNV/SV Class 4–5
- **MAY-include (cancer VUS)** — Cancer Hotspots residue, OncoKB oncogenic/likely-oncogenic, or truncating LoF in CGC Tier 1 tumor suppressors. Capped at 10.
- **v2.2 B1 consequence gate** — Tier III VUS selector rejects variants whose primary VEP consequence is non-coding (intronic, UTR, upstream, downstream). `P/LP` must-reason bypasses the gate unconditionally so deep-intronic ClinVar-Pathogenic splice variants still pass. Splice-region / synonymous variants are rescued when SpliceAI `delta_max >= 0.2`.
- **v2.2 B2 MMR/Lynch carve-out** — protein-impacting VUS in any MMR gene (MLH1, MSH2, MSH6, PMS2, EPCAM) is admitted as `VUS_MMR_Lynch` (priority between `VUS_hotspot` and `VUS_TSG_LoF`), independent of hotspot / TSG / HPO status. The B1 consequence gate still applies — no Lynch exception for intronic variants.
- **MAY-include (rare disease VUS)** — requires HPO phenotype match above threshold. ACMG SF v3.2 genes are excluded from silent inclusion (opt-in path not yet modeled).
- **Soft caps** — 30 for cancer, 20 for rare disease; MUST-include items are never truncated. When MAY-include overflows, selector emits `truncated=True` and `n_dropped=<n>` metadata.
- **No top-k fallback** — when nothing qualifies, the selector emits an empty selection with `empty_reason`, and the report renders "No reportable somatic alterations" / "No candidate variants meeting diagnostic criteria" instead of fabricating low-evidence picks. QC, coverage, TMB, and methodology are still reported.
- **Audit trail** — every selected variant carries `selection_reason` ∈ {`P/LP`, `Tier_I`, `Tier_II`, `Tier_III_hotspot`, `VUS_hotspot`, `VUS_MMR_Lynch`, `VUS_TSG_LoF`, `VUS_HPO_match`, …}
- **Pre-analytic filtering caption** — every AI Board HTML section shows "Pre-analytic filtering: N variants → M presented to Board / Criteria: …" so clinicians can see exactly how much the selector pruned.
- **No-findings placeholder** — empty `AgentOpinion` panels render "No specific findings identified for this case." (EN) / "이 케이스에서 특별한 소견은 확인되지 않았습니다." (KO)

### Report Regeneration Tooling
- **`orchestrate.py`** now dumps `clinical_board` as a dict (asdict) in the run JSON, so reports can be regenerated from cached results
- **`scripts/rerender_report.py`** — regenerate HTML from the cached orchestrate JSON without re-running Ollama. Useful for template tweaks and showcase rebuilds.

### Cancer Report
- FoundationOne CDx-style tiered layout
- OncoKB-based Tier 1–4 with CIViC treatment evidence and PMID citations
- Tumor Mutational Burden (TMB) — nonsynonymous variants per Mb
- Cancer hotspot detection (CIViC-derived)

### Rare Disease Report
- HPO phenotype-driven candidate gene ranking
- OMIM genemap2 gene-disease associations with inheritance patterns (AD/AR/XL/XLD/XLR/MT)
- ClinGen gene validity scores
- Orphanet prevalence data
- GeneReviews cross-references

### Structural Variants
- AnnotSV TSV parsing for CNV/SV
- ACMG CNV 2020 classification display (Class 1–5)
- Cytoband, size, phenotype annotations

### Pipeline
- VEP pre-annotated VCF parsing (CSQ fields including in silico scores)
- Batch processing with variant deduplication across samples
- On-premise deployment (Docker) — no external API dependency required
- 100% offline analysis with local databases

---

## Quick Start

### 1. Install Dependencies
```bash
pip install -r requirements.txt
```

### 2. Set Up Databases
```bash
# Downloads and builds all public reference databases
bash scripts/setup_databases.sh

# For manual downloads (OMIM genemap2, ClinGen), see script output
```

### 3. Run Analysis
```bash
# Cancer report (default)
python scripts/orchestrate.py sample.vcf -o report.html --skip-api --hide-vus

# Rare disease with HPO terms
python scripts/orchestrate.py patient.vcf --mode rare-disease \
  --hpo HP:0001250,HP:0001263 -o report.html

# With AI Clinical Board
python scripts/orchestrate.py sample.vcf --clinical-board --board-lang en -o report.html

# Batch processing
python scripts/orchestrate.py --batch vcf_dir/ --output-dir reports/ --workers 8
```

### Docker
```bash
docker build -t biko-genomeboard .
docker run -v ./data/db:/app/data/db -v ./input:/app/input -v ./output:/app/output \
  biko-genomeboard /app/input/sample.vcf -o /app/output/report.html --skip-api --hide-vus
```

---

## Data Sources

| Source | Type | Size | Setup |
|--------|------|------|-------|
| ClinVar | SQLite (NCBI TSV) | ~1.5 GB | `setup_databases.sh` |
| gnomAD v4.1 | Tabix VCF (direct query) | ~700 GB | Manual download |
| CIViC | SQLite (civicdb.org TSV) | ~5 MB | `setup_databases.sh` |
| HPO | SQLite (gene-phenotype) | ~5 MB | `setup_databases.sh` |
| OMIM mim2gene | SQLite (MIM mapping) | ~1 MB | `setup_databases.sh` |
| OMIM genemap2 | SQLite (gene-disease) | ~2 MB | Manual (OMIM account) |
| Orphanet | SQLite (prevalence XML) | ~2 MB | `setup_databases.sh` |
| GeneReviews | SQLite (NCBI FTP) | ~1 MB | `setup_databases.sh` |
| ClinGen | SQLite (gene validity CSV) | ~1 MB | Manual (web export) |
| KRGDB | TSV (Korean freq) | ~1 MB | Pre-included |
| Korea4K | TSV (Korean freq) | ~1 MB | Pre-included |
| NARD2 | TSV (Korean freq) | ~1 MB | Pre-included |

---

## Configuration

All settings in `config.yaml`:

| Key | Default | Description |
|-----|---------|-------------|
| `annotation.source` | `auto` | `local` / `api` / `auto` |
| `thresholds.ba1` | `0.05` | Stand-alone benign AF |
| `thresholds.bs1` | `0.01` | Strong benign AF |
| `thresholds.pm2` | `0.001` | PM2_Supporting AF cutoff |
| `in_silico.revel_pp3_strong` | `0.932` | REVEL PP3_Strong threshold |
| `in_silico.revel_pp3_moderate` | `0.644` | REVEL PP3_Moderate threshold |
| `in_silico.revel_bp4_strong` | `0.016` | REVEL BP4_Strong threshold |
| `in_silico.spliceai_pp3_strong` | `0.5` | SpliceAI PP3_Strong threshold |
| `clinical_board.enabled` | `false` | Enable AI Clinical Board |
| `clinical_board.agent_model` | `alibayram/medgemma:27b` | Ollama model for specialists |
| `pgx.genes` | 12 genes | PGx gene list |

---

## Project Structure

```
gb/
├── scripts/
│   ├── orchestrate.py              # Main CLI entry point
│   ├── setup_databases.sh          # Database download & build helper
│   ├── pipeline/
│   │   ├── query.py                # Variant database queries
│   │   ├── classify.py             # Classification & record building
│   │   └── batch.py                # Batch processing engine
│   ├── classification/
│   │   ├── acmg_engine.py          # ACMG/AMP 2015 classifier
│   │   ├── in_silico.py            # REVEL/CADD/AlphaMissense/SpliceAI → PP3/BP4
│   │   └── evidence_collector.py   # PVS1, PM1, PM4, PM5, PP2, BP7
│   ├── clinical/
│   │   ├── oncokb.py               # OncoKB tiering + hotspot
│   │   ├── hpo_matcher.py          # HPO phenotype scoring
│   │   ├── query_omim.py           # OMIM gene-disease lookup
│   │   └── query_clingen.py        # ClinGen validity
│   ├── clinical_board/
│   │   ├── agents/                 # 4 specialists + Board Chair
│   │   ├── ollama_client.py        # Ollama REST API client
│   │   ├── case_briefing.py        # Pipeline results → LLM briefing
│   │   ├── runner.py               # Sequential agent execution
│   │   └── render.py               # Board opinion → HTML
│   ├── somatic/
│   │   └── amp_tiering.py          # AMP/ASCO/CAP 2017 tiering
│   ├── db/                         # Database builders (build_*.py)
│   ├── korean_pop/                 # KRGDB, Korea4K, NARD2, gnomAD, freq comparison
│   ├── pharma/                     # PGx analysis
│   ├── intake/                     # VCF/AnnotSV/InterVar parsers
│   ├── counselor/                  # Report generation (Jinja2 → HTML/PDF)
│   └── common/                     # Config, models, cache
├── templates/
│   ├── cancer/report.html          # Cancer report template
│   └── rare-disease/report.html    # Rare disease report template
├── data/
│   ├── sample_vcf/                 # Demo VCF files
│   └── db/                         # Local databases (built by setup script)
├── tests/                          # 899 pytest tests
├── docs/
│   ├── showcase/                   # Sample reports & intro document
│   ├── SETUP.md
│   ├── ARCHITECTURE.md
│   ├── KOREAN_STRATEGY.md
│   └── TIERING_PRINCIPLES.md
├── config.yaml
├── Dockerfile
└── requirements.txt
```

---

## Testing

```bash
pip install -r requirements-dev.txt
python -m pytest tests/ -v
```

899 tests covering ACMG classification, in silico PP3/BP4 thresholds, ClinVar override + v2.2 hotspot conflict reconciliation, CIViC integration, OncoKB curator + 401-degradation, PM1 hotspot table lookup, TMB calculation, CNV/SV parsing, HPO matching, Korean frequency comparison, PGx analysis, AI Clinical Board agents (rare-disease + cancer), curate-then-narrate curator, narrative scrubber and template-renderer fallback, variant selector v2.2 B1/B2 consequence gate and MMR Lynch carve-out, knowledge base CRUD, batch deduplication, config validation, and report generation. CI green on ubuntu-latest / Python 3.10, 3.11, 3.12.

---

## AI Clinical Board Attribution

The AI Clinical Board feature uses [Google MedGemma](https://ai.google.dev/gemma/docs/medgemma) via Ollama for local inference. MedGemma is provided under the [MedGemma Terms of Use](https://ai.google.dev/gemma/docs/medgemma/terms). This feature is for research assistance only, not for diagnostic purposes.

---

## Documentation

| Doc | Contents |
|-----|----------|
| [docs/SETUP.md](docs/SETUP.md) | Installation and database setup |
| [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) | System architecture and data flow |
| [docs/KOREAN_STRATEGY.md](docs/KOREAN_STRATEGY.md) | Korean population analysis strategy |
| [docs/TIERING_PRINCIPLES.md](docs/TIERING_PRINCIPLES.md) | Variant tiering principles |

---

## License

MIT — see [LICENSE](LICENSE).
