# BIKO GenomeBoard

Korean Population-Aware Genomic Variant Interpretation Platform

> **Research Use Only** — Not for clinical diagnosis or medical decision-making.

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

### AI Clinical Board
- Local multi-agent diagnostic synthesis powered by **MedGemma 27B** (Ollama)
- 4 domain specialists: Variant Pathologist, Disease Geneticist, PGx Specialist, Literature Analyst
- Board Chair synthesis with consensus opinion and differential diagnoses
- Bilingual output (English / Korean)
- Deterministic classification is never altered by AI — AI provides interpretive synthesis only

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
├── tests/                          # 638 pytest tests
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

638 tests covering ACMG classification, in silico PP3/BP4 thresholds, ClinVar override, CIViC integration, OncoKB tiering, TMB calculation, CNV/SV parsing, HPO matching, Korean frequency comparison, PGx analysis, AI Clinical Board agents, batch deduplication, config validation, and report generation.

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
