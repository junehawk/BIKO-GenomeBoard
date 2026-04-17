#!/usr/bin/env bash
# ============================================================================
# PharmCAT Setup — downloads PharmCAT JAR and validates dependencies
#
# Usage:  bash scripts/tools/setup_pharmcat.sh
#
# PharmCAT is a pharmacogenomics clinical annotation tool developed by
# Stanford/CPIC (https://pharmcat.org). It requires Java 17+ and
# optionally bcftools for VCF preprocessing.
#
# Installs to:  data/tools/pharmcat/
#   - pharmcat.jar              PharmCAT core JAR
#   - pharmcat_positions.vcf.bgz  PGx position reference VCF
# ============================================================================

set -euo pipefail

# ── Project paths ──────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
PHARMCAT_DIR="$PROJECT_ROOT/data/tools/pharmcat"

# ── Colors (consistent with setup_databases.sh) ───────────────────────────

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info()  { echo -e "${GREEN}[INFO]${NC} $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }
log_step()  { echo -e "\n${BLUE}--- $* ---${NC}"; }

# ── 1. Check Java 17+ ─────────────────────────────────────────────────────

log_step "1/5  Checking Java 17+"

JAVA_CMD=""

# Search order: JAVA_HOME > project-local JDK > system PATH
if [ -n "${JAVA_HOME:-}" ] && [ -x "$JAVA_HOME/bin/java" ]; then
    JAVA_CMD="$JAVA_HOME/bin/java"
else
    # Check project-local JDK (installed by this script or manually)
    LOCAL_JDK_DIR="$PROJECT_ROOT/data/tools/java"
    if [ -d "$LOCAL_JDK_DIR" ]; then
        for jdk_dir in "$LOCAL_JDK_DIR"/jdk-*/Contents/Home/bin/java \
                       "$LOCAL_JDK_DIR"/*/Contents/Home/bin/java; do
            if [ -x "$jdk_dir" ] 2>/dev/null; then
                JAVA_CMD="$jdk_dir"
                break
            fi
        done
    fi
    # Fallback to system PATH
    if [ -z "$JAVA_CMD" ] && command -v java >/dev/null 2>&1; then
        JAVA_CMD="java"
    fi
fi

if [ -z "$JAVA_CMD" ]; then
    log_warn "Java not found. Attempting to download Adoptium JDK 17..."
    LOCAL_JDK_DIR="$PROJECT_ROOT/data/tools/java"
    mkdir -p "$LOCAL_JDK_DIR"

    # Detect architecture
    ARCH=$(uname -m)
    case "$ARCH" in
        arm64|aarch64) ADOPT_ARCH="aarch64" ;;
        x86_64)        ADOPT_ARCH="x64" ;;
        *)             log_error "Unsupported architecture: $ARCH"; exit 1 ;;
    esac

    ADOPT_URL=$(curl -fsSL --connect-timeout 15 --max-time 30 \
        "https://api.adoptium.net/v3/assets/latest/17/hotspot?architecture=${ADOPT_ARCH}&image_type=jdk&os=mac&vendor=eclipse" 2>/dev/null \
        | python3 -c "
import json, sys
data = json.load(sys.stdin)
for item in data:
    pkg = item.get('binary', {}).get('package', {})
    if pkg.get('name', '').endswith('.tar.gz'):
        print(pkg['link'])
        break
" 2>/dev/null) || true

    if [ -n "$ADOPT_URL" ]; then
        log_info "Downloading JDK 17 from Adoptium..."
        JDK_TMP="/tmp/adoptium-jdk17.tar.gz"
        if curl -fSL --progress-bar --connect-timeout 30 --max-time 600 \
            -o "$JDK_TMP" "$ADOPT_URL" 2>&1; then
            cd "$LOCAL_JDK_DIR" && tar --no-same-owner --no-same-permissions -xzf "$JDK_TMP" 2>&1
            rm -f "$JDK_TMP"
            # Find the extracted java binary
            for jdk_dir in "$LOCAL_JDK_DIR"/jdk-*/Contents/Home/bin/java; do
                if [ -x "$jdk_dir" ] 2>/dev/null; then
                    JAVA_CMD="$jdk_dir"
                    break
                fi
            done
            if [ -n "$JAVA_CMD" ]; then
                log_info "JDK 17 installed to $LOCAL_JDK_DIR"
            fi
        fi
    fi

    if [ -z "$JAVA_CMD" ]; then
        log_error "Java not found and auto-download failed."
        log_error "PharmCAT requires Java 17+."
        log_error "Install options:"
        log_error "  - Adoptium (recommended): https://adoptium.net/"
        log_error "  - Homebrew:  brew install openjdk@17"
        log_error "  - SDKMAN:    sdk install java 17-tem"
        exit 1
    fi
fi

# Parse major version from 'java -version' output.
# Output varies by vendor:
#   openjdk version "17.0.10" 2024-01-16
#   java version "21.0.2" 2024-01-16
JAVA_VERSION_OUTPUT=$("$JAVA_CMD" -version 2>&1 | head -1)
JAVA_MAJOR=$(echo "$JAVA_VERSION_OUTPUT" | sed -E 's/.*"([0-9]+).*/\1/')

if [ -z "$JAVA_MAJOR" ] || [ "$JAVA_MAJOR" -lt 17 ] 2>/dev/null; then
    log_error "Java $JAVA_MAJOR detected, but PharmCAT requires Java 17+."
    log_error "Current: $JAVA_VERSION_OUTPUT"
    log_error "Install Java 17+: https://adoptium.net/"
    exit 1
fi

log_info "Java $JAVA_MAJOR OK ($JAVA_VERSION_OUTPUT)"

# ── 2. Check bcftools (optional) ──────────────────────────────────────────

log_step "2/5  Checking bcftools (optional)"

BCFTOOLS_VERSION=""
if command -v bcftools >/dev/null 2>&1; then
    BCFTOOLS_VERSION=$(bcftools --version 2>&1 | head -1)
    log_info "bcftools found: $BCFTOOLS_VERSION"
else
    log_warn "bcftools not found."
    log_warn "bcftools is optional but recommended for VCF preprocessing."
    log_warn "Install: brew install bcftools  (or conda install -c bioconda bcftools)"
fi

# ── 3. Download PharmCAT JAR ──────────────────────────────────────────────

log_step "3/5  Downloading PharmCAT JAR"

mkdir -p "$PHARMCAT_DIR"

PHARMCAT_JAR="$PHARMCAT_DIR/pharmcat.jar"

if [ -f "$PHARMCAT_JAR" ]; then
    log_info "PharmCAT JAR already exists at $PHARMCAT_JAR, skipping download"
else
    log_info "Querying GitHub API for latest PharmCAT release..."

    GITHUB_API="https://api.github.com/repos/PharmGKB/PharmCAT/releases/latest"

    # Fetch release metadata
    RELEASE_JSON=$(curl -fsSL --connect-timeout 15 --max-time 30 "$GITHUB_API" 2>&1) || {
        log_error "Failed to query GitHub API."
        log_error "Manual download:"
        log_error "  1. Visit https://github.com/PharmGKB/PharmCAT/releases/latest"
        log_error "  2. Download pharmcat-*.jar"
        log_error "  3. Place at: $PHARMCAT_JAR"
        exit 1
    }

    # Extract JAR download URL from assets
    JAR_URL=$(echo "$RELEASE_JSON" | python3 -c "
import json, sys
data = json.load(sys.stdin)
for asset in data.get('assets', []):
    name = asset.get('name', '')
    if name.endswith('.jar') and 'pharmcat' in name.lower():
        print(asset['browser_download_url'])
        break
" 2>/dev/null)

    if [ -z "$JAR_URL" ]; then
        log_error "Could not find PharmCAT JAR asset in latest release."
        log_error "Manual download:"
        log_error "  1. Visit https://github.com/PharmGKB/PharmCAT/releases/latest"
        log_error "  2. Download pharmcat-*.jar"
        log_error "  3. Place at: $PHARMCAT_JAR"
        exit 1
    fi

    RELEASE_TAG=$(echo "$RELEASE_JSON" | python3 -c "
import json, sys
data = json.load(sys.stdin)
print(data.get('tag_name', 'unknown'))
" 2>/dev/null)

    log_info "Latest release: $RELEASE_TAG"
    log_info "Downloading: $JAR_URL"

    curl -fSL --progress-bar --connect-timeout 30 --max-time 300 \
        -o "$PHARMCAT_JAR.tmp" "$JAR_URL" 2>&1 || {
        log_error "Failed to download PharmCAT JAR."
        log_error "Manual download:"
        log_error "  1. Visit https://github.com/PharmGKB/PharmCAT/releases/latest"
        log_error "  2. Download pharmcat-*.jar"
        log_error "  3. Place at: $PHARMCAT_JAR"
        rm -f "$PHARMCAT_JAR.tmp"
        exit 1
    }

    mv "$PHARMCAT_JAR.tmp" "$PHARMCAT_JAR"
    log_info "PharmCAT JAR saved to $PHARMCAT_JAR ($(du -h "$PHARMCAT_JAR" | cut -f1))"
fi

# ── 4. Download pharmcat_positions VCF ────────────────────────────────────

log_step "4/5  Downloading PharmCAT positions VCF"

POSITIONS_VCF="$PHARMCAT_DIR/pharmcat_positions.vcf.bgz"

if [ -f "$POSITIONS_VCF" ]; then
    log_info "Positions VCF already exists at $POSITIONS_VCF, skipping download"
else
    # Positions VCF is bundled inside the JAR in newer versions.
    # Try extracting from JAR first; if not found, download from release assets.
    log_info "Checking if positions VCF is bundled in JAR..."

    POSITIONS_IN_JAR=$("$JAVA_CMD" -jar "$PHARMCAT_JAR" -h 2>&1 | grep -i "position" || true)

    # Try to extract from release assets
    if [ -n "${RELEASE_JSON:-}" ]; then
        POSITIONS_URL=$(echo "$RELEASE_JSON" | python3 -c "
import json, sys
data = json.load(sys.stdin)
for asset in data.get('assets', []):
    name = asset.get('name', '')
    if 'position' in name.lower() and name.endswith(('.bgz', '.vcf.bgz', '.vcf.gz')):
        print(asset['browser_download_url'])
        break
" 2>/dev/null)
    else
        POSITIONS_URL=""
    fi

    if [ -n "$POSITIONS_URL" ]; then
        log_info "Downloading positions VCF from release assets..."
        curl -fSL --progress-bar --connect-timeout 30 --max-time 300 \
            -o "$POSITIONS_VCF.tmp" "$POSITIONS_URL" 2>&1 || {
            log_warn "Failed to download positions VCF from release assets."
            log_warn "PharmCAT may still work if positions are bundled in the JAR."
            rm -f "$POSITIONS_VCF.tmp"
        }
        if [ -f "$POSITIONS_VCF.tmp" ]; then
            mv "$POSITIONS_VCF.tmp" "$POSITIONS_VCF"
            log_info "Positions VCF saved to $POSITIONS_VCF ($(du -h "$POSITIONS_VCF" | cut -f1))"
        fi
    else
        log_info "Positions VCF not found in release assets."
        log_info "PharmCAT v2.x bundles positions inside the JAR — this is expected."
    fi
fi

# ── 5. Verify installation ───────────────────────────────────────────────

log_step "5/5  Verifying installation"

PHARMCAT_VERSION=$("$JAVA_CMD" -jar "$PHARMCAT_JAR" --version 2>&1 || true)

if [ -z "$PHARMCAT_VERSION" ]; then
    # Some PharmCAT versions use -version instead of --version
    PHARMCAT_VERSION=$("$JAVA_CMD" -jar "$PHARMCAT_JAR" -version 2>&1 || true)
fi

if echo "$PHARMCAT_VERSION" | grep -qi "pharmcat\|version\|[0-9]\+\.[0-9]\+"; then
    log_info "PharmCAT verification: OK"
else
    # PharmCAT may not have a --version flag; verify JAR is valid by listing main class
    if "$JAVA_CMD" -jar "$PHARMCAT_JAR" -h >/dev/null 2>&1 || "$JAVA_CMD" -jar "$PHARMCAT_JAR" --help >/dev/null 2>&1; then
        log_info "PharmCAT JAR is valid (help output OK)"
        PHARMCAT_VERSION="(version from help)"
    else
        log_warn "Could not verify PharmCAT version. JAR may still be functional."
        PHARMCAT_VERSION="unknown"
    fi
fi

# ── Summary ───────────────────────────────────────────────────────────────

echo ""
echo -e "${BLUE}============================================${NC}"
echo -e "${BLUE}  PharmCAT Setup Summary${NC}"
echo -e "${BLUE}============================================${NC}"
echo ""
echo -e "  Install path:    $PHARMCAT_DIR"
echo -e "  PharmCAT JAR:    $(du -h "$PHARMCAT_JAR" 2>/dev/null | cut -f1 || echo 'N/A')"
echo -e "  PharmCAT ver:    $PHARMCAT_VERSION"
echo -e "  Java version:    Java $JAVA_MAJOR ($JAVA_VERSION_OUTPUT)"

if [ -f "$POSITIONS_VCF" ]; then
    echo -e "  Positions VCF:   $(du -h "$POSITIONS_VCF" | cut -f1)"
else
    echo -e "  Positions VCF:   bundled in JAR (v2.x)"
fi

if [ -n "$BCFTOOLS_VERSION" ]; then
    echo -e "  bcftools:        $BCFTOOLS_VERSION"
else
    echo -e "  bcftools:        ${YELLOW}not installed (optional)${NC}"
fi

echo ""
log_info "PharmCAT setup complete."
