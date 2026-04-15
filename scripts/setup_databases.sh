#!/usr/bin/env bash
# ============================================================================
# BIKO GenomeBoard — Database Setup Script
#
# Downloads external reference databases and builds SQLite indexes.
# Run from project root:  bash scripts/setup_databases.sh [--all | --skip-gnomad]
#
# Sources:
#   ClinVar     — NCBI (public, ~120 MB gz → ~1.5 GB SQLite)
#   CIViC       — civicdb.org (public, auto-download in build script)
#   HPO         — hpo.jax.org (public, ~20 MB)
#   OMIM mim2gene — omim.org (public, ~1 MB)
#   Orphanet    — orphadata.com (public, ~15 MB XML)
#   GeneReviews — NCBI FTP (public, ~200 KB)
#   ClinGen     — clinicalgenome.org (manual CSV export required)
#   OMIM genemap2 — omim.org (requires OMIM account, manual download)
#   gnomAD v4.1 — Broad Institute (optional, ~700 GB for full exome VCFs)
# ============================================================================

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA_DIR="$PROJECT_ROOT/data/db"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info()  { echo -e "${GREEN}[INFO]${NC} $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }
log_step()  { echo -e "\n${BLUE}━━━ $* ━━━${NC}"; }

mkdir -p "$DATA_DIR"
mkdir -p "$DATA_DIR/civic"

# ── Parse arguments ─────────────────────────────────────────────────────────

SKIP_GNOMAD=true
SKIP_MANUAL=true

for arg in "$@"; do
    case "$arg" in
        --all)        SKIP_GNOMAD=false; SKIP_MANUAL=false ;;
        --skip-gnomad) SKIP_GNOMAD=true ;;
        --include-gnomad) SKIP_GNOMAD=false ;;
        --help|-h)
            echo "Usage: bash scripts/setup_databases.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --all             Download everything including gnomAD (~700 GB)"
            echo "  --skip-gnomad     Skip gnomAD download (default)"
            echo "  --include-gnomad  Download gnomAD exome VCFs (requires ~700 GB disk)"
            echo "  -h, --help        Show this help"
            echo ""
            echo "Databases requiring manual download:"
            echo "  - OMIM genemap2.txt (requires OMIM account)"
            echo "  - ClinGen gene_validity CSV (manual export from website)"
            echo ""
            echo "Place manual files in: data/db/"
            exit 0
            ;;
        *)
            log_error "Unknown option: $arg"
            exit 1
            ;;
    esac
done

FAILED=()
SKIPPED=()
SUCCESS=()

# ── Helper: download with retry ─────────────────────────────────────────────

download_file() {
    local url="$1"
    local output="$2"
    local desc="$3"
    local max_retries=3

    if [ -f "$output" ]; then
        log_info "$desc — already exists, skipping download"
        return 0
    fi

    for attempt in $(seq 1 $max_retries); do
        log_info "$desc — downloading (attempt $attempt/$max_retries)..."
        if curl -fSL --progress-bar --connect-timeout 30 --max-time 600 \
            -o "$output.tmp" "$url" 2>&1; then
            mv "$output.tmp" "$output"
            log_info "$desc — downloaded ($(du -h "$output" | cut -f1))"
            return 0
        fi
        log_warn "Attempt $attempt failed, retrying..."
        rm -f "$output.tmp"
        sleep 2
    done

    log_error "$desc — download failed after $max_retries attempts"
    rm -f "$output.tmp"
    return 1
}

# ── Helper: build SQLite DB ─────────────────────────────────────────────────

build_db() {
    local script="$1"
    local desc="$2"
    shift 2

    log_info "Building $desc..."
    if python "$PROJECT_ROOT/$script" "$@" 2>&1; then
        log_info "$desc — build complete"
        return 0
    else
        log_error "$desc — build failed"
        return 1
    fi
}

# ============================================================================
# 1. ClinVar
# ============================================================================
log_step "1/8  ClinVar (NCBI)"

CLINVAR_GZ="$DATA_DIR/variant_summary.txt.gz"
CLINVAR_DB="$DATA_DIR/clinvar.sqlite3"

if [ -f "$CLINVAR_DB" ]; then
    log_info "ClinVar DB already exists, skipping"
    SUCCESS+=("ClinVar")
else
    if download_file \
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz" \
        "$CLINVAR_GZ" "ClinVar variant_summary.txt.gz"; then
        if build_db "scripts/db/build_clinvar_db.py" "ClinVar SQLite"; then
            SUCCESS+=("ClinVar")
        else
            FAILED+=("ClinVar (build)")
        fi
    else
        FAILED+=("ClinVar (download)")
    fi
fi

# ============================================================================
# 2. CIViC
# ============================================================================
log_step "2/8  CIViC (civicdb.org)"

CIVIC_DB="$DATA_DIR/civic.sqlite3"
CIVIC_DIR="$DATA_DIR/civic"
CIVIC_BASE="https://civicdb.org/downloads"

if [ -f "$CIVIC_DB" ]; then
    log_info "CIViC DB already exists, skipping"
    SUCCESS+=("CIViC")
else
    civic_ok=true
    for tsv in "01-Jan-2026/01-Jan-2026-GeneSummaries.tsv" \
               "01-Jan-2026/01-Jan-2026-VariantSummaries.tsv" \
               "01-Jan-2026/01-Jan-2026-ClinicalEvidenceSummaries.tsv"; do
        fname=$(basename "$tsv" | sed 's/[0-9]*-[A-Za-z]*-[0-9]*-//')
        if [ ! -f "$CIVIC_DIR/$fname" ]; then
            # CIViC download URLs change monthly; try the build script directly
            log_warn "CIViC TSV not found: $fname"
            civic_ok=false
        fi
    done

    if [ "$civic_ok" = true ]; then
        if build_db "scripts/db/build_civic_db.py" "CIViC SQLite"; then
            SUCCESS+=("CIViC")
        else
            FAILED+=("CIViC (build)")
        fi
    else
        log_warn "CIViC TSVs not found. Download manually from https://civicdb.org/downloads"
        log_warn "Place GeneSummaries.tsv, VariantSummaries.tsv, ClinicalEvidenceSummaries.tsv in $CIVIC_DIR/"
        SKIPPED+=("CIViC (manual download needed)")
    fi
fi

# ============================================================================
# 3. HPO (genes_to_phenotype)
# ============================================================================
log_step "3/8  HPO (hpo.jax.org)"

HPO_TXT="$DATA_DIR/genes_to_phenotype.txt"
HPO_DB="$DATA_DIR/hpo.sqlite3"

if [ -f "$HPO_DB" ]; then
    log_info "HPO DB already exists, skipping"
    SUCCESS+=("HPO")
else
    if download_file \
        "https://hpo.jax.org/data/annotations/genes_to_phenotype.txt" \
        "$HPO_TXT" "HPO genes_to_phenotype.txt"; then
        if build_db "scripts/db/build_hpo_db.py" "HPO SQLite"; then
            SUCCESS+=("HPO")
        else
            FAILED+=("HPO (build)")
        fi
    else
        FAILED+=("HPO (download)")
    fi
fi

# ============================================================================
# 4. OMIM mim2gene (public)
# ============================================================================
log_step "4/8  OMIM mim2gene (omim.org)"

MIM2GENE="$DATA_DIR/mim2gene.txt"
OMIM_MAP_DB="$DATA_DIR/omim_mapping.sqlite3"

if [ -f "$OMIM_MAP_DB" ]; then
    log_info "OMIM mapping DB already exists, skipping"
    SUCCESS+=("OMIM mim2gene")
else
    if download_file \
        "https://omim.org/static/omim/data/mim2gene.txt" \
        "$MIM2GENE" "OMIM mim2gene.txt"; then
        if build_db "scripts/db/build_omim_mapping.py" "OMIM mapping SQLite"; then
            SUCCESS+=("OMIM mim2gene")
        else
            FAILED+=("OMIM mim2gene (build)")
        fi
    else
        FAILED+=("OMIM mim2gene (download)")
    fi
fi

# ============================================================================
# 5. OMIM genemap2 (requires login — check if file exists)
# ============================================================================
log_step "5/8  OMIM genemap2 (requires OMIM account)"

GENEMAP2="$DATA_DIR/genemap2.txt"
OMIM_GM_DB="$DATA_DIR/omim_genemap.sqlite3"

if [ -f "$OMIM_GM_DB" ]; then
    log_info "OMIM genemap DB already exists, skipping"
    SUCCESS+=("OMIM genemap2")
elif [ -f "$GENEMAP2" ]; then
    if build_db "scripts/db/build_omim_genemap_db.py" "OMIM genemap SQLite"; then
        SUCCESS+=("OMIM genemap2")
    else
        FAILED+=("OMIM genemap2 (build)")
    fi
else
    log_warn "OMIM genemap2.txt not found — requires OMIM account"
    log_warn "  1. Register at https://omim.org/downloads"
    log_warn "  2. Download genemap2.txt"
    log_warn "  3. Place at: $GENEMAP2"
    log_warn "  4. Re-run this script"
    SKIPPED+=("OMIM genemap2 (manual download)")
fi

# ============================================================================
# 6. Orphanet (prevalence data)
# ============================================================================
log_step "6/8  Orphanet (orphadata.com)"

ORPHANET_XML="$DATA_DIR/en_product9_prev.xml"
ORPHANET_DB="$DATA_DIR/orphanet.sqlite3"

if [ -f "$ORPHANET_DB" ]; then
    log_info "Orphanet DB already exists, skipping"
    SUCCESS+=("Orphanet")
else
    if download_file \
        "https://www.orphadata.com/data/xml/en_product9_prev.xml" \
        "$ORPHANET_XML" "Orphanet en_product9_prev.xml"; then
        if build_db "scripts/db/build_orphanet_db.py" "Orphanet SQLite"; then
            SUCCESS+=("Orphanet")
        else
            FAILED+=("Orphanet (build)")
        fi
    else
        FAILED+=("Orphanet (download)")
    fi
fi

# ============================================================================
# 7. GeneReviews (NCBI FTP)
# ============================================================================
log_step "7/8  GeneReviews (NCBI FTP)"

GR_GENES="$DATA_DIR/GRshortname_NBKid_genesymbol_dzname.txt"
GR_TITLES="$DATA_DIR/GRtitle_shortname_NBKid.txt"
GR_DB="$DATA_DIR/genreviews.sqlite3"

if [ -f "$GR_DB" ]; then
    log_info "GeneReviews DB already exists, skipping"
    SUCCESS+=("GeneReviews")
else
    gr_ok=true
    download_file \
        "https://ftp.ncbi.nlm.nih.gov/pub/GeneReviews/GRshortname_NBKid_genesymbol_dzname.txt" \
        "$GR_GENES" "GeneReviews genes file" || gr_ok=false
    download_file \
        "https://ftp.ncbi.nlm.nih.gov/pub/GeneReviews/GRtitle_shortname_NBKid.txt" \
        "$GR_TITLES" "GeneReviews titles file" || gr_ok=false

    if [ "$gr_ok" = true ]; then
        if build_db "scripts/db/build_genreviews_db.py" "GeneReviews SQLite"; then
            SUCCESS+=("GeneReviews")
        else
            FAILED+=("GeneReviews (build)")
        fi
    else
        FAILED+=("GeneReviews (download)")
    fi
fi

# ============================================================================
# 8. ClinGen (manual CSV export)
# ============================================================================
log_step "8/8  ClinGen Gene Validity"

CLINGEN_CSV="$DATA_DIR/clingen_gene_validity.csv"
CLINGEN_DB="$DATA_DIR/clingen.sqlite3"

# An empty DB shell (file exists but has no gene_validity table) is the
# known catch-22 failure mode from an earlier run where the manual CSV
# export was missing. If we detect that shape, remove the shell so the
# CSV/build branch below can re-trigger with a fresh download.
if [ -f "$CLINGEN_DB" ]; then
    if python3 -c "
import sqlite3, sys
try:
    c = sqlite3.connect('$CLINGEN_DB')
    if c.execute(\"SELECT name FROM sqlite_master WHERE type='table' AND name='gene_validity'\").fetchone():
        sys.exit(0)
    sys.exit(1)
finally:
    c.close()
" 2>/dev/null; then
        log_info "ClinGen DB already exists, skipping"
        SUCCESS+=("ClinGen")
        # Skip to next section — need a subshell to preserve control flow.
        # We use a helper flag instead of an early return.
        CLINGEN_DONE=1
    else
        log_warn "ClinGen DB at $CLINGEN_DB is an empty shell (no gene_validity table)"
        log_warn "Removing it so the CSV build path can re-trigger on the next run."
        rm -f "$CLINGEN_DB"
    fi
fi

if [ -n "${CLINGEN_DONE:-}" ]; then
    :  # already handled above
elif [ -f "$CLINGEN_CSV" ]; then
    if build_db "scripts/db/build_clingen_db.py" "ClinGen SQLite"; then
        SUCCESS+=("ClinGen")
    else
        FAILED+=("ClinGen (build)")
    fi
else
    log_warn "ClinGen CSV not found — requires manual export"
    log_warn "  1. Visit https://search.clinicalgenome.org/kb/gene-validity"
    log_warn "  2. Click 'Download' → CSV format"
    log_warn "  3. Place at: $CLINGEN_CSV"
    log_warn "  4. Re-run this script"
    SKIPPED+=("ClinGen (manual download)")
fi

# ============================================================================
# Optional: gnomAD (very large)
# ============================================================================
if [ "$SKIP_GNOMAD" = false ]; then
    log_step "OPTIONAL  gnomAD v4.1 Exomes"
    log_warn "gnomAD download is very large (~700 GB total)"
    log_warn "BIKO GenomeBoard uses tabix (pysam) direct query on VCF files"
    log_warn "Download per-chromosome VCFs from: https://gnomad.broadinstitute.org/downloads#v4"

    GNOMAD_DIR="$DATA_DIR/gnomad_vcf"
    mkdir -p "$GNOMAD_DIR"

    GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes"
    for chr in $(seq 1 22) X Y; do
        vcf="gnomad.exomes.v4.1.sites.chr${chr}.vcf.bgz"
        tbi="${vcf}.tbi"
        if [ -f "$GNOMAD_DIR/$vcf" ]; then
            log_info "chr${chr} — already exists"
        else
            download_file "$GNOMAD_BASE/$vcf" "$GNOMAD_DIR/$vcf" "gnomAD chr${chr}" || true
            download_file "$GNOMAD_BASE/$tbi" "$GNOMAD_DIR/$tbi" "gnomAD chr${chr} index" || true
        fi
    done
    SUCCESS+=("gnomAD (partial)")
else
    log_step "OPTIONAL  gnomAD — skipped (use --include-gnomad to download)"
    SKIPPED+=("gnomAD (use --include-gnomad)")
fi

# ============================================================================
# Summary
# ============================================================================
echo ""
echo -e "${BLUE}════════════════════════════════════════════${NC}"
echo -e "${BLUE}  BIKO GenomeBoard — Database Setup Summary ${NC}"
echo -e "${BLUE}════════════════════════════════════════════${NC}"
echo ""

if [ ${#SUCCESS[@]} -gt 0 ]; then
    echo -e "${GREEN}Completed (${#SUCCESS[@]}):${NC}"
    for s in "${SUCCESS[@]}"; do echo -e "  ✓ $s"; done
fi

if [ ${#SKIPPED[@]} -gt 0 ]; then
    echo ""
    echo -e "${YELLOW}Skipped (${#SKIPPED[@]}):${NC}"
    for s in "${SKIPPED[@]}"; do echo -e "  ⊘ $s"; done
fi

if [ ${#FAILED[@]} -gt 0 ]; then
    echo ""
    echo -e "${RED}Failed (${#FAILED[@]}):${NC}"
    for s in "${FAILED[@]}"; do echo -e "  ✗ $s"; done
    echo ""
    exit 1
fi

echo ""
log_info "Done. Run 'python -m pytest tests/ -x -q' to verify."
