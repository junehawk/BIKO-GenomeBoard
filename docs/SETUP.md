# GenomeBoard Setup Guide

## Prerequisites

| Tool | Version | Notes |
|------|---------|-------|
| Python | 3.10+ | Virtual environment recommended |
| bcftools | 1.17+ | Optional, for VCF pre-filtering |
| tabix | any | Required for gnomAD tabix queries |

---

## 1. Install Python dependencies

```bash
python -m venv .venv
source .venv/bin/activate    # Windows: .venv\Scripts\activate

pip install -r requirements.txt
```

Key dependencies:
- `cyvcf2` — VCF parsing
- `pysam` — tabix queries for gnomAD
- `WeasyPrint` + `Jinja2` — HTML/PDF report generation
- `requests` — ClinVar/gnomAD API fallback

**macOS build issues:**
```bash
# WeasyPrint
brew install pango cairo gdk-pixbuf libffi
pip install WeasyPrint

# cyvcf2
brew install htslib
pip install cyvcf2
```

---

## 2. Local Database Setup (Recommended)

Local databases eliminate API rate limits and enable fully offline operation. Build order: ClinVar → gnomAD → CIViC.

### ClinVar SQLite (4.4M variants)

```bash
# Download variant summary from NCBI
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

# Build SQLite (takes ~2–5 min, produces ~1.5 GB)
python scripts/db/build_clinvar_db.py data/db/variant_summary.txt.gz
# Output: data/db/clinvar.sqlite3
```

### gnomAD v4.1 Exomes (tabix)

GenomeBoard queries gnomAD VCFs directly via tabix — no import step required.

1. Download exomes v4.1 per-chromosome VCF + `.tbi` index files from:
   https://gnomad.broadinstitute.org/downloads#v4

2. Place files in `data/db/gnomad_vcf/`:
   ```
   data/db/gnomad_vcf/
   ├── gnomad.exomes.v4.1.sites.chr1.vcf.bgz
   ├── gnomad.exomes.v4.1.sites.chr1.vcf.bgz.tbi
   ├── gnomad.exomes.v4.1.sites.chr2.vcf.bgz
   ├── gnomad.exomes.v4.1.sites.chr2.vcf.bgz.tbi
   └── ...
   ```

3. Verify tabix is available: `tabix --version`

The pipeline auto-detects files in that directory via `paths.gnomad_vcf_dir` in `config.yaml`.

**Partial downloads:** Downloading all chromosomes (~700 GB) is optional. Any missing chromosome falls back to the gnomAD GraphQL API automatically when `annotation.source` is `auto`.

### CIViC (958 genes, 4812 evidence items)

```bash
# Auto-downloads TSV files from civicdb.org and builds SQLite
python scripts/db/build_civic_db.py
# Output: data/db/civic.sqlite3
```

### HPO Gene-Phenotype (Rare Disease mode)

Required for offline rare disease candidate ranking. Without this DB, HPO-gene associations depend on the JAX API.

```bash
# Download annotation file (~5MB)
curl -o data/db/genes_to_phenotype.txt \
  https://hpo.jax.org/data/annotations/genes_to_phenotype.txt

# Build SQLite
python scripts/db/build_hpo_db.py data/db/genes_to_phenotype.txt
# Output: data/db/hpo.sqlite3
```

### ClinGen Gene-Disease Validity (Rare Disease mode)

Replaces the minimal built-in static data (7 genes) with the full ClinGen gene-validity dataset.

1. Go to https://search.clinicalgenome.org/kb/gene-validity
2. Click "Download" → CSV format
3. Save to `data/db/clingen_gene_validity.csv`

```bash
python scripts/db/build_clingen_db.py data/db/clingen_gene_validity.csv
# Output: data/db/clingen.sqlite3
```

---

## 3. Configuration

All settings are in `config.yaml` in the project root.

```yaml
annotation:
  source: "auto"   # "local" | "api" | "auto" (local-first, API fallback)

paths:
  clinvar_db: "data/db/clinvar.sqlite3"
  gnomad_vcf_dir: "data/db/gnomad_vcf"
  civic_db: "data/db/civic.sqlite3"
  hpo_db: "data/db/hpo.sqlite3"
  clingen_db: "data/db/clingen.sqlite3"
  krgdb: "data/krgdb_freq.tsv"

thresholds:
  ba1: 0.05     # Stand-alone benign (>5%)
  bs1: 0.01     # Strong benign (>=1%)
  pm2: 0.001    # PM2_Supporting (<=0.1%)

cache:
  enabled: true
  ttl_seconds: 604800   # 7 days
```

For Docker/on-premise, override paths with container mount points (see section 5).

---

## 4. Running the Pipeline

### Single sample

```bash
# Cancer mode (default), with API
python scripts/orchestrate.py sample.vcf -o report.html

# Fully offline (local DBs only)
python scripts/orchestrate.py sample.vcf -o report.html --skip-api --hide-vus

# Rare disease mode with HPO phenotypes
python scripts/orchestrate.py patient.vcf --mode rare-disease \
  --hpo HP:0001250,HP:0001263 -o report.html

# Also write JSON output
python scripts/orchestrate.py sample.vcf -o report.html --json

# PDF output
python scripts/orchestrate.py sample.vcf -o report.pdf
```

### Batch processing

Batch mode parses all VCFs, deduplicates variants across samples (annotating each unique variant only once), then generates per-sample reports in parallel.

```bash
# Directory of VCFs
python scripts/orchestrate.py --batch vcf_dir/ --output-dir output/batch --workers 8

# Manifest CSV (columns: sample_id,vcf_path)
python scripts/orchestrate.py --batch manifest.csv --output-dir output/batch --hide-vus
```

Batch manifest format:
```csv
sample_id,vcf_path
SAMPLE001,/data/vcfs/s001.vcf
SAMPLE002,/data/vcfs/s002.vcf
```

### VCF pre-filtering (optional)

GenomeBoard processes PASS-filter variants by default. Pre-filter with bcftools if needed:

```bash
# PASS only
bcftools view -f PASS input.vcf > filtered.vcf

# Compress and index (required for tabix gnomAD queries on input VCFs)
bgzip filtered.vcf
tabix -p vcf filtered.vcf.gz
```

> GRCh38 is required. KRGDB and gnomAD v4.1 coordinates are GRCh38.

---

## 5. Docker

### Build

```bash
docker build -t genomeboard .
```

### Single sample

```bash
docker run \
  -v ./data/db:/app/data/db \
  -v ./input:/app/input \
  -v ./output:/app/output \
  genomeboard /app/input/sample.vcf -o /app/output/report.html --skip-api --hide-vus
```

### Batch mode

```bash
docker run \
  -v ./data/db:/app/data/db \
  -v ./input:/app/input \
  -v ./output:/app/output \
  genomeboard --batch /app/input/ --output-dir /app/output/batch --workers 4
```

### docker-compose

```bash
docker-compose up
```

See `docker-compose.yml` for volume and environment configuration.

---

## 6. Testing

```bash
pip install -r requirements-dev.txt
python -m pytest tests/ -v

# Single module
python -m pytest tests/test_oncokb.py -v
python -m pytest tests/test_batch.py -v
python -m pytest tests/test_rare_disease.py -v
```

372 tests pass. Coverage includes ACMG engine, ClinVar override, CIViC hotspot detection, batch deduplication, tabix gnomAD, HPO matching, and report generation.
