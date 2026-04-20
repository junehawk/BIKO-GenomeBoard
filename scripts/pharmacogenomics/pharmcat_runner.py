"""PharmCAT subprocess wrapper for pharmacogenomics analysis.

Runs PharmCAT (https://pharmcat.org) on a germline VCF to produce
star allele diplotypes, metabolizer phenotypes, and CPIC/DPWG drug
recommendations.  Falls back gracefully when PharmCAT (Java 17+ / JAR)
is not installed — the caller handles the builtin-PGx fallback path.
"""

import glob
import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from scripts.common.config import get

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass
class PharmCATResult:
    """Parsed output of a PharmCAT run."""

    diplotypes: Dict[str, str] = field(default_factory=dict)
    phenotypes: Dict[str, str] = field(default_factory=dict)
    drug_recommendations: List[Dict] = field(default_factory=list)
    source: str = "pharmcat"
    version: str = ""
    warnings: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Availability checks
# ---------------------------------------------------------------------------

_MIN_JAVA_VERSION = 17


def _find_java(config: dict | None = None) -> str | None:
    """Return path to a Java >= 17 binary, or None."""
    cfg_java = ""
    if config:
        cfg_java = config.get("java_path", "")
    if not cfg_java:
        cfg_java = get("pharmacogenomics.java_path", "")

    candidates: list[str] = []
    if cfg_java:
        candidates.append(cfg_java)
    # Auto-discover Adoptium JDK installed by setup_pharmcat.sh under
    # data/tools/java/jdk-*/. Checked before system PATH so the known-good
    # Java 17 takes precedence over an older system-wide installation.
    project_root = Path(__file__).resolve().parent.parent.parent
    local_java_dir = project_root / "data" / "tools" / "java"
    if local_java_dir.exists():
        for jdk_dir in sorted(local_java_dir.iterdir(), reverse=True):
            for bin_path in [
                jdk_dir / "Contents" / "Home" / "bin" / "java",  # macOS layout
                jdk_dir / "bin" / "java",  # Linux layout
            ]:
                if bin_path.exists():
                    candidates.append(str(bin_path))
                    break
    system_java = shutil.which("java")
    if system_java:
        candidates.append(system_java)

    for java in candidates:
        try:
            proc = subprocess.run(
                [java, "-version"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            # Java version string can be on stdout or stderr
            version_text = proc.stderr + proc.stdout
            match = re.search(r'"(\d+)[\._]', version_text)
            if match and int(match.group(1)) >= _MIN_JAVA_VERSION:
                return java
        except (OSError, subprocess.TimeoutExpired):
            continue
    return None


def _find_jar(config: dict | None = None) -> str | None:
    """Return path to the PharmCAT JAR, or None."""
    cfg_jar = ""
    if config:
        cfg_jar = config.get("pharmcat_jar", "")
    if not cfg_jar:
        cfg_jar = get("pharmacogenomics.pharmcat_jar", "")

    if cfg_jar and os.path.isfile(cfg_jar):
        return cfg_jar

    # Auto-discover under data/tools/pharmcat/
    # Match both "pharmcat.jar" (setup_pharmcat.sh default) and
    # "pharmcat-2.13.0.jar" (versioned naming from GitHub releases).
    project_root = Path(__file__).resolve().parent.parent.parent
    pattern = str(project_root / "data" / "tools" / "pharmcat" / "pharmcat*.jar")
    jars = sorted(glob.glob(pattern))
    if jars:
        return jars[-1]  # newest version

    # Environment variable fallback
    env_jar = os.environ.get("PHARMCAT_JAR", "")
    if env_jar and os.path.isfile(env_jar):
        return env_jar

    return None


def is_pharmcat_available(config: dict | None = None) -> bool:
    """Return True if both Java 17+ and the PharmCAT JAR are present."""
    return _find_java(config) is not None and _find_jar(config) is not None


def _pre_filter_pgx_positions(germline_vcf: str, jar_path: str, output_dir: str) -> Optional[str]:
    """Extract only pharmacogene positions from a large germline VCF.

    Reads ``pharmcat_positions.vcf.bgz`` (shipped alongside the JAR) to get
    the ~1000 positions PharmCAT actually inspects, then uses pysam tabix to
    pull matching records from the germline VCF into a tiny temp VCF that
    PharmCAT can process in seconds instead of minutes.

    Returns the path to the filtered VCF on success, or None if pre-filtering
    is not possible (pysam missing, positions file missing, germline not indexed).
    """
    try:
        import pysam
    except ImportError:
        return None

    # Locate pharmcat_positions.vcf.bgz next to the JAR
    jar_dir = Path(jar_path).parent
    positions_file = None
    for candidate in [
        jar_dir / "pharmcat_positions.vcf.bgz",
        jar_dir / "pharmcat_positions.vcf.gz",
    ]:
        if candidate.exists():
            positions_file = candidate
            break
    if not positions_file:
        return None

    # Check germline VCF has tabix index
    tbi = germline_vcf + ".tbi"
    csi = germline_vcf + ".csi"
    if not (Path(tbi).exists() or Path(csi).exists()):
        return None

    # Parse target positions from PharmCAT positions file
    targets: list[tuple[str, int]] = []
    try:
        import gzip

        opener = gzip.open if str(positions_file).endswith(".gz") or str(positions_file).endswith(".bgz") else open
        with opener(str(positions_file), "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t", 3)
                if len(parts) >= 2:
                    chrom = parts[0]
                    pos = int(parts[1])
                    targets.append((chrom, pos))
    except Exception:
        return None

    if not targets:
        return None

    # Extract matching records from germline VCF
    filtered_path = os.path.join(output_dir, "pgx_filtered.vcf")
    try:
        tbx = pysam.TabixFile(germline_vcf)
        header_lines = list(tbx.header)

        with open(filtered_path, "w", encoding="utf-8") as out:
            for h in header_lines:
                out.write(h + "\n")
            for chrom, pos in targets:
                try:
                    for record in tbx.fetch(chrom, pos - 1, pos):
                        out.write(record + "\n")
                except ValueError:
                    continue
        tbx.close()
    except Exception:
        return None

    return filtered_path


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


def run_pharmcat(
    germline_vcf: str,
    output_dir: str | None = None,
    config: dict | None = None,
) -> Optional[PharmCATResult]:
    """Run PharmCAT on *germline_vcf* and return parsed results.

    Returns ``None`` when PharmCAT is not installed or execution fails,
    so that the caller can fall back to the builtin PGx engine.
    """
    java = _find_java(config)
    jar = _find_jar(config)
    if not java or not jar:
        logger.info("PharmCAT not available (java=%s, jar=%s)", java, jar)
        return None

    if not os.path.isfile(germline_vcf):
        logger.warning("Germline VCF not found: %s", germline_vcf)
        return None

    use_tmpdir = output_dir is None
    if use_tmpdir:
        tmpdir = tempfile.mkdtemp(prefix="pharmcat_")
        output_dir = tmpdir
    else:
        os.makedirs(output_dir, exist_ok=True)
        tmpdir = None

    try:
        # ── 0. Pre-filter: extract only PGx positions from germline VCF ──
        # A full germline VCF (434 MB, ~5M variants) causes PharmCAT to
        # time out. PharmCAT only inspects ~1000 pharmacogene positions
        # listed in pharmcat_positions.vcf.bgz. We use pysam tabix to
        # extract just those positions into a tiny temp VCF (<100 KB),
        # then pass that to PharmCAT for near-instant processing.
        input_vcf = germline_vcf
        pre_filtered_vcf = None
        try:
            pre_filtered_vcf = _pre_filter_pgx_positions(germline_vcf, jar, output_dir)
            if pre_filtered_vcf:
                input_vcf = pre_filtered_vcf
                logger.info("Pre-filtered germline to PGx positions: %s", pre_filtered_vcf)
        except Exception as pf_err:
            logger.debug("PGx pre-filter failed (%s), using full VCF", pf_err)

        # ── 1. Run PharmCAT core ──────────────────────────────────────────
        cmd = [
            java,
            "-jar",
            jar,
            "-vcf",
            input_vcf,
            "-o",
            output_dir,
            "-reporterJson",
        ]
        logger.info("Running PharmCAT: %s", " ".join(cmd))

        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
        )

        if proc.returncode != 0:
            logger.warning(
                "PharmCAT exited with code %d: %s",
                proc.returncode,
                proc.stderr[:500],
            )
            return None

        # ── 2. Locate and parse the reporter JSON ─────────────────────────
        # Multi-sample VCFs produce one report per sample. Try to pick the
        # proband's report by matching the germline VCF filename stem against
        # the report filenames; fall back to the first report otherwise.
        report_jsons = sorted(glob.glob(os.path.join(output_dir, "*.report.json")))
        if not report_jsons:
            logger.warning("PharmCAT produced no report JSON in %s", output_dir)
            return None

        report_path = report_jsons[0]
        if len(report_jsons) > 1:
            vcf_stem = Path(germline_vcf).stem.replace(".vcf", "").replace(".gz", "").replace(".bgz", "")
            for rj in report_jsons:
                # PharmCAT names reports as <vcf_stem>.<sample_id>.report.json
                if vcf_stem.lower() in Path(rj).name.lower():
                    report_path = rj
                    break
            logger.info(
                "Multi-sample VCF: %d reports found, selected %s",
                len(report_jsons),
                Path(report_path).name,
            )

        return _parse_pharmcat_json(report_path)

    except subprocess.TimeoutExpired:
        logger.warning("PharmCAT timed out after 300 s")
        return None
    except Exception as e:
        logger.warning("PharmCAT execution failed: %s", e)
        return None
    finally:
        if tmpdir:
            shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# JSON parser
# ---------------------------------------------------------------------------


def _parse_pharmcat_json(report_path: str) -> Optional[PharmCATResult]:
    """Parse a PharmCAT ``*.report.json`` into a :class:`PharmCATResult`.

    The parser is deliberately lenient — it extracts only the fields BIKO
    needs and silently skips unexpected structures so that minor PharmCAT
    version differences do not crash the pipeline.
    """
    try:
        with open(report_path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("Failed to read PharmCAT report JSON: %s", exc)
        return None

    diplotypes: Dict[str, str] = {}
    phenotypes: Dict[str, str] = {}
    drug_recs: List[Dict] = []
    warnings: List[str] = []
    version = data.get("version", data.get("pharmcatVersion", ""))

    # PharmCAT v3 emits ``genes`` as a **dict** keyed by gene symbol,
    # not a list.  Each value is the gene-level result object.
    genes_raw = data.get("genes", {})
    genes_iter = genes_raw.values() if isinstance(genes_raw, dict) else genes_raw

    for gene_entry in genes_iter:
        if not isinstance(gene_entry, dict):
            continue
        gene = gene_entry.get("geneSymbol", gene_entry.get("gene", ""))
        if not gene:
            continue

        # Diplotype — from recommendationDiplotypes (preferred) or sourceDiplotypes
        rec_dips = gene_entry.get("recommendationDiplotypes", [])
        if rec_dips and isinstance(rec_dips[0], dict):
            label = rec_dips[0].get("label", "")
            if label:
                diplotypes[gene] = label
            # Phenotype lives inside recommendationDiplotypes[0].phenotypes
            pheno_list = rec_dips[0].get("phenotypes", [])
            if pheno_list:
                phenotypes[gene] = (
                    pheno_list[0] if isinstance(pheno_list[0], str) else pheno_list[0].get("phenotype", "")
                )

        if gene not in diplotypes:
            src_dips = gene_entry.get("sourceDiplotypes", [])
            if src_dips and isinstance(src_dips[0], dict):
                allele1 = src_dips[0].get("allele1", {}).get("name", "")
                allele2 = src_dips[0].get("allele2", {}).get("name", "")
                if allele1 and allele2:
                    diplotypes[gene] = f"{allele1}/{allele2}"

    # Drug recommendations — PharmCAT v3 has top-level ``drugs`` dict
    # keyed by guideline name, each containing drug-level entries.
    drugs_raw = data.get("drugs", {})
    drugs_iter = drugs_raw.values() if isinstance(drugs_raw, dict) else drugs_raw
    for guideline_group in drugs_iter:
        if not isinstance(guideline_group, dict):
            continue
        for drug_name, drug_data in guideline_group.items():
            if not isinstance(drug_data, dict):
                continue
            guideline = drug_data.get("guideline", drug_data.get("source", ""))
            for rec in drug_data.get("recommendations", []):
                if not isinstance(rec, dict):
                    continue
                classification = rec.get("classification", "")
                drug_recs.append(
                    {
                        "gene": ", ".join(rec.get("genes", [])) if rec.get("genes") else "",
                        "drug": drug_name,
                        "guideline": guideline,
                        "classification": classification,
                    }
                )

    return PharmCATResult(
        diplotypes=diplotypes,
        phenotypes=phenotypes,
        drug_recommendations=drug_recs,
        version=version,
        warnings=warnings,
    )
