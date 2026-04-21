# BIKO GenomeBoard Setup Guide

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

> **Quick path:** `bash scripts/setup_databases.sh` automates ClinVar, gnomAD, CIViC, HPO, Orphanet, GeneReviews, KOVA v7 (Korean Variant Archive), and cancerhotspots. The per-step instructions below describe what the script does and which sources still require manual setup (OMIM genemap2, ClinGen, KOVA v7 bundle).

### ClinVar SQLite (4.4M variants)

```bash
# Download variant summary from NCBI
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

# Build SQLite (takes ~2–5 min, produces ~1.5 GB)
python scripts/db/build_clinvar_db.py data/db/variant_summary.txt.gz
# Output: data/db/clinvar.sqlite3
```

### gnomAD v4.1 Exomes (tabix)

BIKO GenomeBoard queries gnomAD VCFs directly via tabix — no import step required.

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

### Orphanet Disease Prevalence (Rare Disease mode)

Provides disease prevalence data per gene from Orphanet Product 9 (~15 MB XML).

```bash
# Download Orphanet Product 9 prevalence XML
curl -o data/db/en_product9_prev.xml \
  https://www.orphadata.com/data/xml/en_product9_prev.xml

# Build SQLite (~few seconds)
python scripts/db/build_orphanet_db.py data/db/en_product9_prev.xml
# Output: data/db/orphanet.sqlite3
```

### GeneReviews Gene→NBK Mapping (Rare Disease mode)

Maps gene symbols to GeneReviews book chapters (NBK IDs) and PMIDs for offline reference lookup.

```bash
# Download GeneReviews identifier files from NCBI FTP
curl -o data/db/GRshortname_NBKid_genesymbol_dzname.txt \
  ftp://ftp.ncbi.nlm.nih.gov/pub/GeneReviews/GRshortname_NBKid_genesymbol_dzname.txt
curl -o data/db/GRtitle_shortname_NBKid.txt \
  ftp://ftp.ncbi.nlm.nih.gov/pub/GeneReviews/GRtitle_shortname_NBKid.txt

# Build SQLite
python scripts/db/build_genreviews_db.py \
  --genes data/db/GRshortname_NBKid_genesymbol_dzname.txt \
  --titles data/db/GRtitle_shortname_NBKid.txt
# Output: data/db/genreviews.sqlite3
```

### KOVA v7 — Korean Variant Archive (Korean frequency cohort)

KOVA (Korean Variant Archive) v7 is the KOGO / gene2korea release used as
BIKO GenomeBoard's Korean population cohort. It supplies 43 million
variants with per-variant allele frequency **and homozygote counts**,
the latter being the feature that gives the strongest Korean-specific
signal for autosomal-recessive review (BS2 candidates).

```bash
# 1. Obtain the KOVA v7 bundle from KOGO / gene2korea
#    (license-gated; follow the publisher's distribution terms)
#    Place the VCF + tabix index locally, then convert to the BIKO TSV layout:
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC_Hom\t%INFO/AC_Het\t%INFO/AN\n' \
  kova_v7.vcf.gz > data/kova_freq.tsv

# 2. Register the TSV path in config.yaml (default):
#    paths:
#      kova: "data/kova_freq.tsv"

# 3. Verify the query path picks it up
python -m pytest tests/test_kova.py -v
```

KOVA v7 replaces the earlier multi-cohort Korean-frequency strategy
(see `docs/KOREAN_STRATEGY.md`). The historical multi-cohort design is
preserved in `docs/superpowers/specs/` under superseded notes for
audit continuity.

### OMIM Gene→MIM Mapping

Maps gene symbols to OMIM MIM numbers for cross-reference links.

```bash
# Download mim2gene.txt (free, no login required)
curl -o data/db/mim2gene.txt https://omim.org/static/omim/data/mim2gene.txt

# Build SQLite
python scripts/db/build_omim_mapping.py data/db/mim2gene.txt
# Output: data/db/omim_mapping.sqlite3
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
  kova: "data/kova_freq.tsv"       # KOVA v7 (KOGO / gene2korea Korean Variant Archive)

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

### Structural Variant / CNV Integration

BIKO GenomeBoard supports CNV/SV analysis from AnnotSV output:

```bash
# Run with SNV + SV
python scripts/orchestrate.py sample.vcf -o report.html \
  --sv annotsv_output.tsv --mode cancer

# SV display rules:
# ACMG Class 4-5: Full detail pages
# ACMG Class 3: Dosage-sensitive VUS in summary table
# ACMG Class 1-2: Count only
```

AnnotSV must be run separately before BIKO GenomeBoard. The `--sv` argument accepts the AnnotSV TSV output file (full annotation rows only; `SV_type=full` rows are automatically selected).

```bash
# Rare disease mode with HPO + SV
python scripts/orchestrate.py patient.vcf --mode rare-disease \
  --hpo HP:0001250,HP:0001263 \
  --sv annotsv_rare.tsv -o report.html
```

---

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

BIKO GenomeBoard processes PASS-filter variants by default. Pre-filter with bcftools if needed:

```bash
# PASS only
bcftools view -f PASS input.vcf > filtered.vcf

# Compress and index (required for tabix gnomAD queries on input VCFs)
bgzip filtered.vcf
tabix -p vcf filtered.vcf.gz
```

> GRCh38 is required. KOVA v7 and gnomAD v4.1 coordinates are GRCh38.

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
python -m pytest tests/ -q

# Single module
python -m pytest tests/test_oncokb.py -q
python -m pytest tests/test_batch.py -q
python -m pytest tests/test_rare_disease.py -q
```

901+ tests pass on Python 3.10 / 3.11 / 3.12. Coverage includes ACMG engine, in silico thresholds, ClinVar override, CIViC/OncoKB integration, TMB, AnnotSV CNV/SV, HPO matching, Korean frequency comparison (KOVA v7), PGx, AI Clinical Board (curate-then-narrate), variant selector, and report generation.

---

## 7. Gene Knowledge Base (Curated Sources)

Build the gene knowledge database from authoritative sources (NCBI Gene, CIViC, GeneReviews, CPIC):

```bash
# Build for specific genes
python -m scripts.tools.build_gene_knowledge --genes TP53,BRCA2,KRAS,BRAF,EGFR

# Build for OncoKB cancer genes (~300 genes, ~5 min)
python -m scripts.tools.build_gene_knowledge --source oncokb

# Build from VCF (extract gene list automatically)
python -m scripts.tools.build_gene_knowledge --vcf input/sample.vcf

# Build all (OncoKB + PGx genes)
python -m scripts.tools.build_gene_knowledge --source all

# Custom output path
python -m scripts.tools.build_gene_knowledge --genes TP53,BRCA2 --output data/gene_knowledge.json
```

Output is written to `data/gene_knowledge.json`. Source priority per gene:
1. **CPIC** (PGx genes: CYP2D6, CYP2C19, etc.) → `curated-cpic`
2. **CIViC local DB** (cancer genes with descriptions) → `curated-civic`
3. **NCBI Gene API** (universal fallback) → `curated-ncbi`
4. **Minimal entry** (all sources failed) → `auto-minimal`

Additional enrichment applied to all genes: GeneReviews PMID (local DB or API), ClinGen gene validity, CIViC treatment evidence, Orphanet disease prevalence, OMIM MIM cross-reference.
