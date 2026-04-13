# BIKO GenomeBoard

Korean Population-Aware Genomic Variant Interpretation Platform

> **Research Use Only** — Not for clinical diagnosis or medical decision-making.

---

## Features

- **Cancer Report (FoundationOne CDx style)**
  - OncoKB-based variant tiering (Tier 1–4)
  - CIViC curated treatment evidence with PMID citations
  - Cancer hotspot detection (CIViC-derived)
  - ClinVar direct classification override (expert panel / multi-submitter)
  - VUS filtering (`--hide-vus`) for focused clinical reports

- **Rare Disease Report**
  - HPO phenotype-driven candidate gene ranking
  - OMIM gene-disease associations
  - ClinGen gene validity scores
  - Inheritance pattern display (AD/AR/XL)

- **Data Sources**
  - ClinVar local SQLite DB (4.4M+ GRCh38 variants)
  - gnomAD exomes v4.1 via tabix (direct VCF query)
  - CIViC gene/variant knowledge (958 genes, 4812 evidence items)
  - KRGDB Korean population frequencies
  - 12 PGx genes with Korean vs Western prevalence comparison

- **Pipeline**
  - VEP/SnpEff pre-annotated VCF parsing (CSQ/ANN fields, rsID extraction)
  - ACMG/AMP 2015 classification engine (deterministic rule engine)
  - Batch processing with variant deduplication across samples
  - SQLite response cache (7-day TTL)
  - Docker packaging for on-premise deployment

---

## Quick Start

### Single Sample
```bash
pip install -r requirements.txt
python scripts/orchestrate.py sample.vcf -o report.html --json
```

### With Local Databases (Recommended)
```bash
# Build ClinVar DB (~4.4M variants)
python scripts/db/build_clinvar_db.py data/db/variant_summary.txt.gz

# Build CIViC DB (auto-downloads from civicdb.org)
python scripts/db/build_civic_db.py

# Run offline with local DBs and VUS filtering
python scripts/orchestrate.py sample.vcf --skip-api --hide-vus -o report.html
```

### Rare Disease Mode
```bash
python scripts/orchestrate.py patient.vcf --mode rare-disease \
  --hpo HP:0001250,HP:0001263 -o report.html
```

### Batch Processing
```bash
python scripts/orchestrate.py --batch vcf_dir/ --output-dir output/batch --workers 8 --hide-vus
```

### Docker
```bash
docker build -t genomeboard .
docker run -v ./data/db:/app/data/db -v ./input:/app/input -v ./output:/app/output \
  genomeboard /app/input/sample.vcf -o /app/output/report.html --skip-api --hide-vus
```

---

## Local Database Setup

### ClinVar
```bash
# Download from NCBI
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

# Build SQLite DB
python scripts/db/build_clinvar_db.py data/db/variant_summary.txt.gz
```

### gnomAD (tabix)
Download exomes v4.1 VCF + `.tbi` index files from https://gnomad.broadinstitute.org/downloads#v4.
Place in `data/db/gnomad_vcf/`. The pipeline queries them directly via tabix with no import step.

### CIViC
```bash
# Auto-downloads TSVs from civicdb.org and builds SQLite
python scripts/db/build_civic_db.py
```

---

## Configuration

All settings live in `config.yaml`. Key options:

| Key | Default | Description |
|-----|---------|-------------|
| `annotation.source` | `auto` | `local` / `api` / `auto` (local-first with API fallback) |
| `paths.clinvar_db` | `data/db/clinvar.sqlite3` | ClinVar SQLite path |
| `paths.gnomad_vcf_dir` | `data/db/gnomad_vcf` | gnomAD tabix VCF directory |
| `paths.civic_db` | `data/db/civic.sqlite3` | CIViC SQLite path |
| `paths.krgdb` | `data/krgdb_freq.tsv` | KRGDB frequency table |
| `thresholds.ba1` | `0.05` | Stand-alone benign allele frequency |
| `thresholds.bs1` | `0.01` | Strong benign allele frequency |
| `thresholds.pm2` | `0.001` | PM2_Supporting frequency cutoff |
| `pgx.genes` | 12 genes | PGx gene list (CYP2D6, HLA-B, etc.) |
| `cache.ttl_seconds` | `604800` | Variant response cache TTL (7 days) |

---

## Report Modes

### Cancer (default)
FoundationOne CDx-style tiered report:

| Tier | Criteria | Display |
|------|----------|---------|
| Tier 1 | Pathogenic/LP on OncoKB Level 1–2 gene; Drug Response; Risk Factor | Full detail + treatment evidence |
| Tier 2 | Pathogenic/LP on any cancer gene; Hotspot VUS | Full detail page |
| Tier 3 | VUS on cancer gene (non-hotspot) | Abbreviated table |
| Tier 4 | All other variants | Count only |

### Rare Disease
Germline variant report with phenotype matching:
- HPO phenotype scoring against associated gene lists
- Candidate gene ranking (Pathogenic first, then by HPO score)
- OMIM/ClinGen annotations per gene
- Inheritance patterns (AD/AR/XL)

---

## Project Structure

```
gb/
├── scripts/
│   ├── orchestrate.py          # Main CLI entry point
│   ├── classification/
│   │   └── acmg_engine.py      # Deterministic ACMG/AMP 2015 classifier
│   ├── clinical/
│   │   ├── oncokb.py           # OncoKB tiering + hotspot detection
│   │   ├── query_clinvar.py    # ClinVar API
│   │   ├── hpo_matcher.py      # HPO phenotype scoring
│   │   ├── query_omim.py       # OMIM gene-disease lookup
│   │   └── query_clingen.py    # ClinGen validity scores
│   ├── db/
│   │   ├── build_clinvar_db.py # ClinVar SQLite builder
│   │   ├── build_civic_db.py   # CIViC SQLite builder
│   │   ├── build_gnomad_db.py  # gnomAD SQLite builder
│   │   ├── query_local_clinvar.py
│   │   ├── query_tabix_gnomad.py
│   │   └── query_civic.py      # CIViC evidence + hotspot queries
│   ├── korean_pop/
│   │   ├── query_krgdb.py      # KRGDB local TSV lookup
│   │   ├── query_gnomad.py     # gnomAD GraphQL API
│   │   └── compare_freq.py     # 3-tier frequency comparison + ACMG codes
│   ├── pharma/
│   │   └── korean_pgx.py       # 12-gene PGx with Korean prevalence
│   ├── counselor/
│   │   └── generate_pdf.py     # Jinja2 HTML/PDF report generator
│   └── common/
│       ├── models.py            # Shared data models
│       ├── cache.py             # SQLite response cache
│       └── config.py            # config.yaml loader
├── templates/                   # Jinja2 report templates
├── data/
│   ├── krgdb_freq.tsv           # Korean population frequencies
│   ├── oncokb_cancer_genes.json # OncoKB gene list
│   ├── acmg_rules.json          # ACMG classification rules
│   └── db/                      # Local database files (build separately)
├── tests/                       # pytest suite (372 tests)
├── config.yaml
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```

---

## Testing

```bash
pip install -r requirements-dev.txt
python -m pytest tests/ -v
```

372 tests covering ACMG classification, ClinVar override logic, CIViC integration, batch deduplication, HPO matching, OncoKB tiering, tabix gnomAD queries, and report generation.

---

## Documentation

| Doc | Contents |
|-----|----------|
| [docs/SETUP.md](docs/SETUP.md) | Installation and database setup |
| [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) | System architecture and data flow |
| [docs/KOREAN_STRATEGY.md](docs/KOREAN_STRATEGY.md) | Korean population analysis strategy |
| [docs/API_KEYS.md](docs/API_KEYS.md) | API key configuration |

---

## License

MIT — see [LICENSE](LICENSE).
