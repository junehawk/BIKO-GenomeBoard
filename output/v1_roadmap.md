# GenomeBoard v1.0 Roadmap

**Date:** 2026-03-23
**Scope:** v0.3 (current) to v1.0 (production-ready)
**Audience:** Development team

---

## Current State Assessment

GenomeBoard is a Korean population-aware genomic variant interpretation system. The standalone CLI pipeline (`scripts/orchestrate.py`) parses VCF files, queries ClinVar/gnomAD/KRGDB, runs ACMG classification with PGx bypass and Risk Factor bypass, and produces a FoundationOne-style HTML/PDF report with Pretendard font.

**What works today:**
- 10-variant demo VCF end-to-end (cancer somatic style)
- ClinVar E-utilities + gnomAD GraphQL + KRGDB TSV frequency pipeline
- ACMG/AMP 2015 classification engine with PGx gene bypass (CYP2D6, CYP2C19, CYP2C9, HLA-B, HLA-A, NUDT15, TPMT, DPYD) and Risk Factor bypass (APOE)
- Korean enrichment/depletion flagging (KRGDB vs gnomAD ratio)
- 5 PGx genes with Korean vs Western prevalence comparison
- Thread-safe caching (ACMG rules, PGx data)
- 127 tests passing
- GRCh38 support (gnomAD v4 primary, v2.1 fallback)
- WeasyPrint graceful degradation (HTML fallback when unavailable)
- Jinja2 autoescape enabled
- 4xx retry exclusion in `api_utils.py`

**What is missing for v1.0:**
- Only one report mode (cancer somatic); no rare disease / germline support
- No web UI; CLI only
- Gene knowledge is LLM-generated (`data/gene_knowledge.json`), no PMID citations
- No VCF annotation (no HGVS, no consequence type from VEP/SnpEff)
- No Docker packaging
- No real VCF stress testing (demo has 10 variants)
- Only 5 PGx genes (missing CYP3A5, CYP1A2, UGT1A1, SLCO1B1, G6PD, VKORC1, etc.)
- No configuration file; paths and thresholds are hardcoded or scattered across modules
- No CI/CD pipeline
- No liftover support (GRCh37 input requires manual conversion)

---

## Phase 1: Foundation (v0.5) — 1-2 weeks

Core infrastructure that every later feature depends on. Focus on configuration management, annotation pipeline, and the report mode abstraction that Phase 2 builds on.

---

### 1.1 Configuration Management (YAML config)

**Description:** Replace all hardcoded paths, thresholds, and settings with a single `config.yaml` file loaded at pipeline startup. This unblocks Docker packaging, Web UI settings, and multi-environment deployments.

**Effort:** M
**Priority:** Must-have
**Dependencies:** None
**Technical approach:**
- Create `config.yaml` at project root with sections: `paths` (krgdb, gene_knowledge, pgx_table, acmg_rules, templates), `thresholds` (BA1, BS1, PM2), `api` (gnomad_url, clinvar_url, ncbi_api_key, timeouts, max_retries), `report` (default_mode, default_genome_build, font_url), `pgx` (gene list)
- Add `scripts/common/config.py` module using `pyyaml` with environment variable override support (`GB_CONFIG_PATH`, `NCBI_API_KEY`)
- Refactor all modules that hardcode paths (`gene_knowledge.py`, `korean_pgx.py`, `query_krgdb.py`, `acmg_engine.py`, `generate_pdf.py`) to accept config or use the singleton config loader
- Keep backward compatibility: if no config file exists, use current defaults
- Add `pyyaml` to `requirements.txt`

**Acceptance criteria:**
- `config.yaml` exists with all current hardcoded values documented
- `python scripts/orchestrate.py demo.vcf` works with and without config file
- Every module path reference traces back to config, not a hardcoded string
- Tests pass with `GB_CONFIG_PATH` override pointing to a test config

---

### 1.2 VCF Annotation Pipeline (VEP/SnpEff integration)

**Description:** Add automated VCF annotation to produce HGVS nomenclature (c./p.), consequence type (missense, frameshift, etc.), and transcript ID. Currently `gene_knowledge.json` provides static HGVS per gene, but real variants need per-variant annotation.

**Effort:** L
**Priority:** Must-have
**Dependencies:** 1.1 (config for tool paths and reference genome settings)
**Technical approach:**
- Create `scripts/intake/annotate_vcf.py` module
- Primary strategy: Ensembl VEP REST API (`https://rest.ensembl.org/vep/human/region/`) for zero-install annotation. POST batch endpoint accepts up to 200 variants per request.
- Fallback strategy: local SnpEff (`java -jar snpEff.jar GRCh38.105`) when configured in `config.yaml`
- Parse VEP JSON response to extract: `transcript_id`, `hgvsc` (c. notation), `hgvsp` (p. notation), `consequence_terms` (missense_variant, frameshift_variant, etc.), `impact` (HIGH/MODERATE/LOW/MODIFIER)
- Extend `Variant` dataclass in `models.py`: add `hgvsc`, `hgvsp`, `consequence`, `transcript_id`, `impact` fields (all Optional[str])
- Wire into `orchestrate.py` as step 1.5 (after VCF parse, before database queries): batch-annotate all variants, attach results to Variant objects
- Update report template to show real HGVS instead of static gene_knowledge HGVS
- Config section: `annotation.method` (vep_rest | snpeff_local | none), `annotation.genome_build` (GRCh38 | GRCh37), `annotation.snpeff_jar_path`

**Acceptance criteria:**
- Demo VCF produces correct HGVS for all 10 variants (verify TP53 c.524G>T, BRCA2 c.7878C>T against Ensembl)
- Report shows per-variant HGVS c./p. notation and consequence type
- Pipeline still runs when VEP API is unreachable (graceful fallback: HGVS fields stay empty, warning logged)
- Annotation step adds < 10 seconds for 10 variants (VEP REST batch)

---

### 1.3 Report Mode Abstraction (cancer vs rare-disease CLI flag)

**Description:** Introduce `--mode cancer` (default, current behavior) and `--mode rare-disease` as a CLI flag and config option. Phase 1 only adds the plumbing (mode parameter threading, template selection logic). Phase 2 fills in the rare disease template content.

**Effort:** S
**Priority:** Must-have
**Dependencies:** 1.1 (config for default mode)
**Technical approach:**
- Add `--mode` argument to `orchestrate.py` argparse: choices=["cancer", "rare-disease"], default from config
- Thread `mode` parameter through `run_pipeline()` into `generate_report_html(report_data, mode=mode)`
- Create template directory structure: `templates/cancer/report.html` (move current template), `templates/rare-disease/report.html` (placeholder in Phase 1, built out in Phase 2), `templates/shared/` (CSS, components, ACMG appendix)
- Refactor `report.html` (1769 lines) into Jinja2 includes: `{% include "shared/head.html" %}`, `{% include "shared/masthead.html" %}`, `{% include "shared/methodology.html" %}`, `{% include "shared/appendix.html" %}` — mode-specific middle sections
- `generate_pdf.py` selects template based on mode parameter
- Add `mode` field to `report_data` JSON output

**Acceptance criteria:**
- `python scripts/orchestrate.py demo.vcf --mode cancer` produces identical output to current behavior
- `python scripts/orchestrate.py demo.vcf --mode rare-disease` produces a valid HTML report with a placeholder body ("Rare disease report mode — under construction") and shared masthead/appendix/CSS
- Report JSON includes `"mode": "cancer"` or `"mode": "rare-disease"`
- Existing tests pass without modification

---

### 1.4 CI/CD Pipeline (GitHub Actions)

**Description:** Add automated linting, testing, and build verification on every push and PR.

**Effort:** S
**Priority:** Should-have
**Dependencies:** None
**Technical approach:**
- Create `.github/workflows/ci.yml`
- Jobs: (1) `lint` — ruff check + ruff format --check, (2) `test` — pytest with `--skip-api` flag (no external API calls in CI), (3) `build` — verify `python scripts/orchestrate.py --help` exits cleanly
- Python matrix: 3.10, 3.11, 3.12
- Cache pip dependencies
- Add `ruff` to `requirements-dev.txt` (new file for dev dependencies)
- Add `pyproject.toml` with ruff config (line-length=120, target-version=py310)

**Acceptance criteria:**
- CI runs on push to `main` and on all PRs
- All 127+ tests pass in CI (with `--skip-api`)
- Ruff lint passes with zero errors
- CI completes in under 3 minutes

---

## Phase 2: Clinical Expansion (v0.7) — 2-3 weeks

Rare disease mode, gene knowledge sourcing from real databases, and expanded PGx coverage.

---

### 2.1 Rare Disease Report Mode — Full Implementation

**Description:** Build the germline/rare disease report template and pipeline extensions. This mode targets clinical genetics workflows: phenotype-driven candidate gene ranking, inheritance pattern display, OMIM disease associations, and ClinGen gene-disease validity scores.

**Effort:** L
**Priority:** Must-have
**Dependencies:** 1.2 (VEP annotation for consequence type), 1.3 (report mode abstraction)
**Technical approach:**

**A. Phenotype input (HPO terms):**
- Add `--hpo` CLI argument accepting comma-separated HPO IDs (e.g., `--hpo HP:0001250,HP:0001263`)
- Create `scripts/clinical/hpo_matcher.py`: query HPO API (`https://ontology.jax.org/api/hp/`) to resolve HPO IDs to names, fetch associated genes
- Rank variants by HPO-gene overlap score: variant gene in HPO-associated gene list = +1 per matching HPO term

**B. OMIM integration:**
- Create `scripts/clinical/query_omim.py`: use OMIM API (`https://api.omim.org/api/`) with `OMIM_API_KEY` from config
- For each variant gene, fetch: MIM number, phenotype description, inheritance pattern (AD/AR/XL/Mito)
- Extend `report_data.variants[]` with: `omim_mim`, `omim_phenotype`, `inheritance_pattern`

**C. ClinGen gene-disease validity:**
- Create `scripts/clinical/query_clingen.py`: query ClinGen Gene Validity API (`https://search.clinicalgenome.org/kb/gene-validity/`)
- Fetch validity classification per gene: Definitive, Strong, Moderate, Limited, Disputed, Refuted
- Display in report alongside OMIM data

**D. Report template (`templates/rare-disease/report.html`):**
- Patient phenotype summary section (HPO terms listed)
- Candidate gene ranking table (sorted by HPO overlap score + ACMG classification)
- Inheritance pattern column in variant table
- OMIM disease association per variant
- ClinGen validity badge per gene
- De novo indicator placeholder (for trio VCF support in future)
- Shared sections: Korean population comparison, PGx, methodology, appendix

**E. ACMG criteria adjustments for germline:**
- PP4 (phenotype specificity): if HPO terms match known gene-phenotype association, add PP4
- PS4 (prevalence in affected): if ClinVar review >= 2 stars pathogenic, weight higher
- No changes to PGx bypass or Risk Factor bypass (shared logic)

**Acceptance criteria:**
- `python scripts/orchestrate.py rare_disease_demo.vcf --mode rare-disease --hpo HP:0001250,HP:0001263` produces a complete report
- Rare disease report shows: HPO phenotypes listed, candidate gene ranking, inheritance patterns, OMIM associations, ClinGen validity
- Cancer mode is unaffected (regression test)
- At least 3 genes in demo VCF have OMIM data resolved

---

### 2.2 Gene Knowledge from Real Sources (PMID citations)

**Description:** Replace LLM-generated gene knowledge text in `data/gene_knowledge.json` with content sourced from public databases, each statement backed by PMID or database accession. This is the single most important credibility improvement.

**Effort:** L
**Priority:** Must-have
**Dependencies:** 1.1 (config for API keys), 1.2 (VEP for consequence data)
**Technical approach:**

**A. Data source integration (new modules):**
- `scripts/common/sources/genereviews.py`: NCBI E-utilities to fetch GeneReviews entries by gene symbol. Extract: gene function summary, clinical significance, management recommendations. Each entry has a GeneReviews NBK ID and associated PMIDs.
- `scripts/common/sources/cpic_api.py`: CPIC API (`https://api.cpicpgx.org/v1/`) for PGx gene-drug pairs. Extract: guideline recommendations, evidence level, PMID of guideline publication.
- `scripts/common/sources/clinvar_knowledge.py`: ClinVar Variation Services API for variant-level evidence summaries, submitter assertions, and review status.
- `scripts/common/sources/pubmed.py`: PubMed E-utilities (`esearch` + `esummary`) to resolve PMIDs to citation strings (Author et al., Journal Year).

**B. Knowledge schema overhaul:**
- Add `references` field to every text block in `gene_knowledge.json`:
  ```json
  {
    "treatment_strategies": "...",
    "treatment_references": ["PMID:25741868", "GeneReviews:NBK1247"],
    "frequency_prognosis": "...",
    "frequency_references": ["PMID:30311383", "KRGDB:2026-03"]
  }
  ```
- Add `source_type` field: "database" (verified), "ai_generated" (legacy), "curated" (manually reviewed)

**C. Report template updates:**
- Superscript citation numbers next to each statement (FoundationOne style)
- References section at bottom of each gene detail page listing full citations
- AI-generated text watermark: "AI-generated summary, not peer-reviewed" badge when `source_type == "ai_generated"`

**D. Build script for knowledge base refresh:**
- `scripts/tools/refresh_knowledge.py`: batch-fetch from GeneReviews, CPIC, ClinVar for all genes in `gene_knowledge.json`; merge with existing data; flag any entries that could not be sourced

**Acceptance criteria:**
- Every `treatment_strategies` and `frequency_prognosis` field in gene_knowledge.json has at least one PMID or database accession in its references array
- Report renders citation numbers (superscript) linking to reference list
- AI-generated content clearly labeled when no database source is available
- `python scripts/tools/refresh_knowledge.py` runs without error and updates at least TP53, BRCA2, CYP2C19 entries
- PGx gene entries cite CPIC guideline PMIDs

---

### 2.3 Expanded PGx Gene Panel

**Description:** Expand from 5 PGx genes to 12+, covering all CPIC Level A gene-drug pairs relevant to the Korean population.

**Effort:** M
**Priority:** Should-have
**Dependencies:** 1.1 (config for PGx gene list), 2.2 (CPIC API integration for recommendations)
**Technical approach:**
- Add to `PGX_GENES` set (in both `acmg_engine.py` and `korean_pgx.py`): CYP3A5, CYP1A2, UGT1A1, SLCO1B1, VKORC1, G6PD, IFNL3
- Extend `data/korean_pgx_table.json` with Korean vs Western prevalence for each new gene (source: CPIC + published Korean PGx studies)
- Add star allele / diplotype-to-phenotype mapping per gene in `korean_pgx.py` (currently only CYP2C19, HLA-B, NUDT15 have phenotype strings)
- Add CPIC recommendation text per gene-drug pair (from CPIC API, with PMID)
- Update report template PGx section to handle 12+ rows cleanly (pagination if > 8)

**New genes with Korean significance:**
| Gene | Key Drug | Korean Prevalence Note |
|------|----------|----------------------|
| CYP3A5 | Tacrolimus | ~75% expressors in Koreans vs ~15% in Caucasians |
| UGT1A1 | Irinotecan | *6 allele ~13% in Koreans (East Asian specific) |
| SLCO1B1 | Statins | *15 ~16% in Koreans |
| VKORC1 | Warfarin | -1639G>A ~90% in Koreans vs ~37% in Caucasians |
| CYP1A2 | Caffeine, Clozapine | *1C allele higher in East Asians |
| G6PD | Rasburicase | Deficiency ~2-5% in Korean males |
| IFNL3 | PEG-IFN/Ribavirin | rs12979860 CC genotype ~85% in Koreans |

**Acceptance criteria:**
- `korean_pgx_table.json` has 12+ gene entries with Korean/Western frequencies and CPIC source
- Demo VCF with CYP3A5, UGT1A1 variants produces correct PGx phenotype and recommendation
- PGx gene list is configurable in `config.yaml` (not only hardcoded frozenset)
- Existing 5-gene tests still pass

---

### 2.4 Cancer Report Enhancements (TMB/MSI placeholders)

**Description:** Add structured placeholders and data model fields for tumor-specific metrics that a real cancer panel would include: Tumor Mutational Burden (TMB), Microsatellite Instability (MSI) status, and a therapeutic options section.

**Effort:** S
**Priority:** Should-have
**Dependencies:** 1.3 (report mode abstraction)
**Technical approach:**
- Add to `report_data` schema: `tmb` (float, mutations/Mb, default null), `msi_status` (string: MSS/MSI-L/MSI-H, default null), `therapeutic_options` (list of dicts with `drug`, `evidence_level`, `reference`)
- Add CLI arguments: `--tmb 10.5`, `--msi MSI-H` (optional manual input; automated calculation requires whole-exome data beyond current scope)
- Cancer report template: TMB/MSI summary box in header section, therapeutic options table per variant (conditional: only shown if data present)
- Therapeutic options populated from CIVIC/OncoKB lookups (Phase 4 stretch goal) or manual JSON input via `--therapeutic-json path.json`

**Acceptance criteria:**
- `--mode cancer --tmb 12.3 --msi MSI-H` renders TMB and MSI in report header
- Without `--tmb`/`--msi` flags, those sections are gracefully hidden (not "null" displayed)
- Therapeutic options table renders when JSON data is provided

---

## Phase 3: User Experience (v0.9) — 2-3 weeks

Web UI, Docker packaging, and polish.

---

### 3.1 Web UI (FastAPI + HTMX)

**Description:** A lightweight web interface for VCF upload, pipeline execution with progress feedback, and report download. No heavy frontend framework; server-rendered HTML with HTMX for interactivity.

**Effort:** L
**Priority:** Should-have
**Dependencies:** 1.1 (config), 1.3 (report mode abstraction)
**Technical approach:**

**A. Backend (`scripts/web/app.py`):**
- FastAPI application with endpoints:
  - `GET /` — upload form (VCF file, mode selector, HPO input, sample ID)
  - `POST /analyze` — accepts multipart VCF upload, validates file, queues analysis
  - `GET /status/{job_id}` — returns job status (queued/running/step/complete/error) as JSON
  - `GET /report/{job_id}` — serves generated HTML report
  - `GET /report/{job_id}/pdf` — serves PDF (if WeasyPrint available)
  - `GET /report/{job_id}/json` — serves raw JSON data
  - `GET /history` — lists past analyses with date, sample ID, variant count, status
- Async pipeline execution: `run_pipeline()` called via `asyncio.to_thread()` in a background task
- Progress tracking: pipeline `_progress()` function writes to a shared progress dict keyed by job_id; `/status` endpoint reads it
- Job storage: SQLite database (`data/jobs.db`) with columns: job_id, sample_id, mode, status, created_at, completed_at, report_path, variant_count
- File storage: uploaded VCFs in `data/uploads/{job_id}/`, reports in `output/{job_id}/`

**B. Frontend (`templates/web/`):**
- `index.html` — upload form with HTMX-powered mode switcher (shows HPO field when rare-disease selected)
- `status.html` — progress page with HTMX polling (`hx-get="/status/{job_id}" hx-trigger="every 2s"`), shows pipeline step progress (1/6 Parsing... 2/6 ClinVar...)
- `history.html` — table of past analyses
- `report_viewer.html` — embedded report with download buttons (HTML, PDF, JSON)
- Shared layout with Pretendard font, matching report color palette
- Mobile-responsive (Tailwind CSS via CDN, minimal)

**C. Dependencies:**
- Add `fastapi`, `uvicorn`, `python-multipart`, `aiosqlite` to requirements
- No Node.js build step; HTMX via CDN (`<script src="https://unpkg.com/htmx.org@1.9"`)

**Acceptance criteria:**
- `uvicorn scripts.web.app:app --port 8000` starts the server
- Upload a VCF, see real-time progress, download HTML/PDF report — full flow works
- History page shows past analyses
- No JavaScript build step required
- Works in Docker container (Phase 3.2)

---

### 3.2 Docker Packaging

**Description:** Single Dockerfile that packages the entire GenomeBoard stack: Python runtime, WeasyPrint system dependencies (cairo, pango, gdk-pixbuf), Pretendard font, and all data files. One `docker run` for both CLI and Web UI.

**Effort:** M
**Priority:** Must-have
**Dependencies:** 1.1 (config), 3.1 (web UI)
**Technical approach:**
- Base image: `python:3.11-slim-bookworm` (Debian Bookworm for WeasyPrint compat)
- System deps: `libcairo2 libpango-1.0-0 libgdk-pixbuf-2.0-0 libffi-dev fonts-noto-cjk`
- Install Pretendard font from GitHub release into `/usr/share/fonts/`
- Multi-stage build: stage 1 builds wheels, stage 2 copies only runtime deps
- Entrypoints:
  - CLI: `docker run genomeboard:latest cli input.vcf -o report.html --mode cancer`
  - Web: `docker run -p 8000:8000 genomeboard:latest web`
  - Entrypoint script (`docker-entrypoint.sh`) dispatches to `orchestrate.py` or `uvicorn`
- Volume mounts: `/data` (for custom KRGDB, config), `/output` (for reports)
- Environment variables: `NCBI_API_KEY`, `OMIM_API_KEY`, `GB_CONFIG_PATH`
- `docker-compose.yml` for convenience with volume mappings pre-configured
- Image size target: < 500MB

**Acceptance criteria:**
- `docker build -t genomeboard .` completes without error
- `docker run genomeboard cli data/sample_vcf/demo_variants.vcf -o /output/report.pdf` produces a PDF with Pretendard font
- `docker run -p 8000:8000 genomeboard web` serves the web UI
- Image size < 500MB
- `docker-compose up` starts web UI with persistent data/output volumes

---

### 3.3 Report Polish and PDF Quality

**Description:** Improve report visual quality, fix PDF rendering edge cases, and add report customization options.

**Effort:** S
**Priority:** Should-have
**Dependencies:** 1.3 (template refactoring), 2.1 (rare disease template)
**Technical approach:**
- Fix WeasyPrint page-break issues: add `page-break-inside: avoid` to variant detail blocks, ensure tables don't split across pages mid-row
- Add print-optimized CSS media query (`@media print`) for browser-based printing
- Cover page with GenomeBoard logo, patient info summary, date, report ID
- Table of contents with page numbers (WeasyPrint supports `target-counter()`)
- Configurable report header/footer text (institution name, lab director) via config.yaml
- Korean language option for report text (`--lang ko` flag; template conditional blocks)

**Acceptance criteria:**
- PDF report has no split-table rendering artifacts
- Cover page renders with configurable institution name
- Report prints cleanly from browser (Chrome, Safari)

---

### 3.4 Liftover Support (GRCh37 auto-detection)

**Description:** Automatically detect VCF reference genome build and perform liftover from GRCh37 to GRCh38 when needed, since gnomAD v4 and current KRGDB coordinates are GRCh38.

**Effort:** M
**Priority:** Should-have
**Dependencies:** 1.1 (config for chain file path), 1.2 (annotation needs correct coordinates)
**Technical approach:**
- Create `scripts/intake/liftover.py`
- Auto-detection: check VCF header for `##reference=` line, check contig lengths (chr1 length differs between GRCh37 and GRCh38), check if `chr` prefix present
- Use `pyliftover` library (pure Python, no external deps) for coordinate conversion
- Chain file: `data/hg19ToHg38.over.chain.gz` (bundled or downloaded on first use)
- Wire into pipeline: after VCF parse, before annotation, if build == GRCh37: liftover all variants, log conversions, flag any that fail to lift
- Add `--genome-build` CLI flag: auto (default), GRCh37, GRCh38 — overrides auto-detection
- Config: `liftover.chain_file`, `liftover.enabled` (default true)

**Acceptance criteria:**
- GRCh37 VCF input produces correct GRCh38 coordinates for all 10 demo variants
- Auto-detection correctly identifies `demo_variants.vcf` (no `chr` prefix = ambiguous, warns) and `demo_variants_grch38.vcf` (GRCh38)
- `--genome-build GRCh37` forces liftover even if auto-detection says GRCh38
- Failed liftover variants are flagged in report, not silently dropped
- `pyliftover` added to requirements.txt

---

## Phase 4: Production (v1.0) — 2-4 weeks

Stress testing, validation, performance, and documentation for production deployment.

---

### 4.1 Real VCF Stress Testing

**Description:** Validate the pipeline against real-world VCF files with 100 to 10,000+ variants. Identify performance bottlenecks, API rate limits, and classification accuracy at scale.

**Effort:** L
**Priority:** Must-have
**Dependencies:** 1.2 (annotation), all Phase 2 features
**Technical approach:**

**A. Test VCF corpus:**
- Curate 5 test VCFs of increasing complexity:
  1. `small_panel.vcf` — 50 variants (targeted gene panel)
  2. `medium_panel.vcf` — 200 variants (comprehensive cancer panel, e.g., MSK-IMPACT-like)
  3. `exome_filtered.vcf` — 1,000 variants (filtered whole exome)
  4. `exome_full.vcf` — 5,000 variants (less-filtered whole exome)
  5. `genome_sampled.vcf` — 10,000 variants (sampled whole genome)
- Sources: ClinVar FTP (downloadable variant VCFs), gnomAD truth sets, synthetic VCFs from `vcf_generator` tools
- Each VCF annotated with expected classifications for a validation subset

**B. Performance profiling:**
- Add `--profile` CLI flag that outputs timing per pipeline step
- Identify bottleneck: likely ClinVar E-utilities (rate limit: 3 req/sec without API key, 10 req/sec with key)
- Mitigation: batch ClinVar queries where possible (efetch supports multiple UIDs), connection pooling with `requests.Session`, configurable concurrency in ThreadPoolExecutor

**C. API rate limit handling:**
- ClinVar: respect NCBI 429 responses, implement exponential backoff with jitter
- gnomAD: GraphQL endpoint has undocumented rate limits; add configurable delay between requests
- Add `--max-concurrent-api` CLI flag (default 4) to control ThreadPoolExecutor workers

**D. Classification accuracy validation:**
- Compare GenomeBoard ACMG classification against ClinVar expert panel classification for 100+ variants
- Compute concordance rate (target: > 90% agreement within 1 tier)
- Document all discordant cases with explanation

**E. Stress test harness:**
- `scripts/tools/stress_test.py`: runs pipeline on all test VCFs, collects timing and accuracy metrics, outputs markdown summary

**Acceptance criteria:**
- 50-variant VCF completes in < 2 minutes
- 1,000-variant VCF completes in < 15 minutes (with API key)
- No crashes or unhandled exceptions on any test VCF
- Classification concordance with ClinVar expert panel > 90% on 100-variant validation set
- Performance report generated as markdown

---

### 4.2 gnomAD Coverage Improvement

**Description:** Improve gnomAD hit rate, which is currently limited by the variant ID format matching and dataset availability. Many clinically relevant variants have gnomAD data but are missed by the current query strategy.

**Effort:** M
**Priority:** Should-have
**Dependencies:** 1.2 (VEP annotation provides rsID and normalized variant IDs)
**Technical approach:**
- Current limitation: `query_gnomad.py` constructs `{chrom}-{pos}-{ref}-{alt}` and queries GraphQL. Multi-allelic sites, indels with left-alignment differences, and VCF normalization mismatches cause missed lookups.
- Solution 1: Use VEP-annotated `colocated_variants` field which includes gnomAD frequencies directly (eliminates separate gnomAD query for VEP-annotated variants)
- Solution 2: For non-VEP runs, add rsID-based gnomAD lookup as fallback (`variant(rsid: "rs28934578")` GraphQL query)
- Solution 3: Normalize variant representation before query (left-align, trim common prefix/suffix) using `scripts/intake/normalize_variant.py`
- Add gnomAD v4.1 dataset support (update `VARIANT_QUERY` dataset parameter)

**Acceptance criteria:**
- gnomAD hit rate on demo VCF improves from current baseline (measure before/after)
- Indel variants (e.g., CFTR c.1521_1523del) successfully resolve gnomAD frequency
- rsID fallback works when position-based query fails

---

### 4.3 Test Suite Expansion

**Description:** Expand test coverage to cover all critical classification paths, edge cases, and integration scenarios. Current 127 tests cover basic functionality but miss boundary conditions in ACMG logic and PGx phenotyping.

**Effort:** M
**Priority:** Must-have
**Dependencies:** All Phase 2 features
**Technical approach:**
- **ACMG engine tests** (`tests/test_acmg_logic.py`): add tests for every rule combination in `acmg_rules.json` (8 pathogenic + 7 likely_pathogenic + 2 benign + 2 likely_benign = 19 rules), conflict resolution, PGx bypass, Risk Factor bypass
- **PGx phenotype tests**: verify diplotype-to-phenotype mapping for all 12+ genes, especially boundary cases (heterozygous vs homozygous)
- **Frequency comparison tests**: test all threshold boundaries (0.001, 0.01, 0.05), Korean enrichment ratio edge cases (zero denominators)
- **VEP annotation tests**: mock VEP API responses, verify HGVS extraction
- **Rare disease mode tests**: HPO matching, OMIM integration, ClinGen validity
- **Integration tests**: end-to-end pipeline with mock API responses for reproducibility
- **Report rendering tests**: verify Jinja2 template renders without errors for all classification types, empty PGx results, missing fields
- Target: 250+ tests, > 85% line coverage on `scripts/` directory

**Acceptance criteria:**
- Every ACMG rule in `acmg_rules.json` has at least one dedicated test
- PGx phenotype correct for all star allele combinations in `korean_pgx_table.json`
- `pytest --cov=scripts` shows > 85% line coverage
- All tests run in < 30 seconds (no real API calls in test suite)

---

### 4.4 Documentation

**Description:** Developer documentation, API documentation (for web UI), and clinical validation documentation.

**Effort:** M
**Priority:** Should-have
**Dependencies:** All phases
**Technical approach:**
- `README.md` overhaul: quick start (CLI + Docker + Web), architecture diagram, configuration reference
- `docs/clinical_validation.md`: ACMG rule implementation details, ClinVar concordance results, PGx accuracy validation, known limitations
- `docs/api.md`: Web UI REST API documentation (FastAPI auto-generates OpenAPI spec at `/docs`)
- `docs/deployment.md`: Docker deployment guide, environment variable reference, NCBI/OMIM API key setup
- `docs/development.md`: local setup, running tests, adding new genes, adding new PGx entries
- Inline docstrings for all public functions (enforce via ruff rule D100)
- `CHANGELOG.md` maintained per version

**Acceptance criteria:**
- New developer can set up local environment from README in < 10 minutes
- Clinical validation document covers all ACMG rules with test evidence
- FastAPI `/docs` endpoint serves interactive API documentation
- All public functions have docstrings

---

### 4.5 Paperclip Agent Architecture Cleanup

**Description:** Rationalize the multi-agent Paperclip integration. Currently 7 agents where 4 are unnecessary LLM wrappers around deterministic Python functions. Reduce to 3 LLM agents + direct Python function calls.

**Effort:** L
**Priority:** Nice-to-have
**Dependencies:** None (independent track, but benefits from stable pipeline)
**Technical approach:**
- Remove LLM agents for: krgdb-agent, vcf-parser-agent, pdf-generator-agent, clinvar-agent (replace with direct Python function calls from the orchestrating agent)
- Keep LLM agents for: counselor-agent (natural language interpretation), cto-agent (orchestration), ceo-agent (clinical judgment)
- Downgrade cto-agent model from `claude-opus-4-6` to `claude-sonnet-4-6` (per roadmap N-2)
- Fix hardcoded CWD in agent configs (use `GB_PROJECT_ROOT` env var)
- Implement Worker Preamble Protocol for remaining LLM agents (replace `CRITICAL: YOU MUST USE` prompt patterns)
- Update `paperclip.manifest.json` to reflect 3-agent architecture

**Acceptance criteria:**
- Agent count reduced from 7 to 3
- End-to-end Paperclip pipeline produces identical report to standalone CLI
- No hardcoded absolute paths in agent configs
- CTO agent uses Sonnet model
- Pipeline response time improves (target: < 3 minutes worst case, down from 10 minutes)

---

## Dependency Graph

```
Phase 1 (Foundation)
  1.1 Config ─────────────┬──────────────────────────────────────┐
  1.2 VEP Annotation ─────┤ (depends on 1.1)                    │
  1.3 Report Mode ─────────┤ (depends on 1.1)                    │
  1.4 CI/CD ───────────────┘ (independent)                       │
                                                                  │
Phase 2 (Clinical)                                                │
  2.1 Rare Disease Mode ────── (depends on 1.2, 1.3)             │
  2.2 Gene Knowledge Sources ── (depends on 1.1, 1.2)            │
  2.3 Expanded PGx ──────────── (depends on 1.1, 2.2)            │
  2.4 Cancer Enhancements ───── (depends on 1.3)                 │
                                                                  │
Phase 3 (UX)                                                      │
  3.1 Web UI ──────────────────── (depends on 1.1, 1.3)          │
  3.2 Docker ──────────────────── (depends on 1.1, 3.1)          │
  3.3 Report Polish ───────────── (depends on 1.3, 2.1)          │
  3.4 Liftover ────────────────── (depends on 1.1, 1.2)          │
                                                                  │
Phase 4 (Production)                                              │
  4.1 Stress Testing ──────────── (depends on 1.2, Phase 2 all)  │
  4.2 gnomAD Coverage ─────────── (depends on 1.2)               │
  4.3 Test Suite ──────────────── (depends on Phase 2 all)       │
  4.4 Documentation ───────────── (depends on all)               │
  4.5 Paperclip Cleanup ────────── (independent)                 │
```

---

## Effort Summary

| Phase | Item | Effort | Priority |
|-------|------|--------|----------|
| **1 Foundation** | 1.1 Config Management | M | Must-have |
| | 1.2 VEP/SnpEff Annotation | L | Must-have |
| | 1.3 Report Mode Abstraction | S | Must-have |
| | 1.4 CI/CD Pipeline | S | Should-have |
| **2 Clinical** | 2.1 Rare Disease Mode | L | Must-have |
| | 2.2 Gene Knowledge Sources | L | Must-have |
| | 2.3 Expanded PGx Panel | M | Should-have |
| | 2.4 Cancer TMB/MSI Placeholders | S | Should-have |
| **3 UX** | 3.1 Web UI (FastAPI + HTMX) | L | Should-have |
| | 3.2 Docker Packaging | M | Must-have |
| | 3.3 Report Polish | S | Should-have |
| | 3.4 Liftover Support | M | Should-have |
| **4 Production** | 4.1 Stress Testing | L | Must-have |
| | 4.2 gnomAD Coverage | M | Should-have |
| | 4.3 Test Suite Expansion | M | Must-have |
| | 4.4 Documentation | M | Should-have |
| | 4.5 Paperclip Cleanup | L | Nice-to-have |

**Must-haves:** 1.1, 1.2, 1.3, 2.1, 2.2, 3.2, 4.1, 4.3 (8 items)
**Should-haves:** 1.4, 2.3, 2.4, 3.1, 3.3, 3.4, 4.2, 4.4 (8 items)
**Nice-to-haves:** 4.5 (1 item)

---

## Risk Register

| Risk | Impact | Mitigation |
|------|--------|------------|
| NCBI API rate limits block stress testing | High | Obtain NCBI API key (10 req/sec); implement batch queries; add configurable delays |
| OMIM API requires academic license | Medium | Apply for OMIM API key early; fallback to ClinVar gene-phenotype data if OMIM unavailable |
| VEP REST API latency for large VCFs | Medium | Batch requests (200 variants/batch); fallback to local SnpEff; cache VEP results |
| WeasyPrint rendering differences across platforms | Low | Docker ensures consistent rendering; test on Linux + macOS |
| GeneReviews content extraction is fragile | Medium | Use structured NCBI Bookshelf API; fall back to AI-generated text with watermark |
| Korean PGx prevalence data sparse for new genes | Low | Use published meta-analyses; document confidence level per gene |

---

## New Dependencies (requirements.txt additions)

```
# Phase 1
pyyaml>=6.0            # config management
pyliftover>=0.4        # GRCh37->38 liftover (Phase 3.4, but install early)
ruff>=0.4.0            # linting (dev dependency)

# Phase 2
# (no new deps; uses existing requests for API calls)

# Phase 3
fastapi>=0.110.0       # web UI
uvicorn>=0.29.0        # ASGI server
python-multipart>=0.0.9 # file upload
aiosqlite>=0.20.0      # async SQLite for job tracking
```

---

## Files Affected (estimated)

| Phase | New Files | Modified Files |
|-------|-----------|----------------|
| Phase 1 | config.yaml, scripts/common/config.py, scripts/intake/annotate_vcf.py, .github/workflows/ci.yml, pyproject.toml | orchestrate.py, models.py, gene_knowledge.py, korean_pgx.py, query_krgdb.py, acmg_engine.py, generate_pdf.py, requirements.txt, templates/report.html (split into includes) |
| Phase 2 | scripts/clinical/hpo_matcher.py, scripts/clinical/query_omim.py, scripts/clinical/query_clingen.py, scripts/common/sources/genereviews.py, scripts/common/sources/cpic_api.py, scripts/common/sources/clinvar_knowledge.py, scripts/common/sources/pubmed.py, scripts/tools/refresh_knowledge.py, templates/rare-disease/report.html | orchestrate.py, korean_pgx.py, acmg_engine.py, gene_knowledge.json, korean_pgx_table.json, templates/shared/* |
| Phase 3 | scripts/web/app.py, scripts/web/models.py, templates/web/*.html, Dockerfile, docker-compose.yml, docker-entrypoint.sh, scripts/intake/liftover.py | orchestrate.py, requirements.txt, generate_pdf.py, templates/ |
| Phase 4 | data/test_vcf/*.vcf, scripts/tools/stress_test.py, docs/*.md, tests/test_acmg_rules_full.py, tests/test_pgx_expanded.py, tests/test_rare_disease.py, tests/test_web.py, CHANGELOG.md | tests/*, paperclip-config/*, README.md |

---

## Milestone Checkpoints

**v0.5 (Phase 1 complete):**
- Config file drives all settings
- VEP annotation produces real HGVS for each variant
- `--mode cancer` and `--mode rare-disease` flags accepted (rare-disease is placeholder)
- CI/CD runs on all PRs

**v0.7 (Phase 2 complete):**
- Rare disease report mode fully functional with HPO input, OMIM, ClinGen
- Gene knowledge backed by PMIDs; AI-generated text labeled
- 12+ PGx genes with Korean prevalence data
- Cancer mode has TMB/MSI placeholder fields

**v0.9 (Phase 3 complete):**
- Web UI serves VCF upload, progress, report download
- Docker image builds and runs both CLI and Web modes
- GRCh37 input auto-detected and lifted over
- PDF rendering polished with cover page

**v1.0 (Phase 4 complete):**
- 1,000-variant VCF processes without crash in < 15 minutes
- > 90% ACMG concordance with ClinVar expert panel
- 250+ tests, > 85% coverage
- Complete developer and clinical documentation
- Production Docker deployment guide
