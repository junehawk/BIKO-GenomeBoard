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
    project_root = Path(__file__).resolve().parent.parent.parent
    pattern = str(project_root / "data" / "tools" / "pharmcat" / "pharmcat-*.jar")
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
        # ── 1. Run PharmCAT core ──────────────────────────────────────────
        cmd = [
            java,
            "-jar",
            jar,
            "-vcf",
            germline_vcf,
            "-o",
            output_dir,
            "-reporterJson",
        ]
        logger.info("Running PharmCAT: %s", " ".join(cmd))

        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
        )

        if proc.returncode != 0:
            logger.warning(
                "PharmCAT exited with code %d: %s",
                proc.returncode,
                proc.stderr[:500],
            )
            return None

        # ── 2. Locate and parse the reporter JSON ─────────────────────────
        report_jsons = glob.glob(os.path.join(output_dir, "*.report.json"))
        if not report_jsons:
            logger.warning("PharmCAT produced no report JSON in %s", output_dir)
            return None

        report_path = report_jsons[0]
        return _parse_pharmcat_json(report_path)

    except subprocess.TimeoutExpired:
        logger.warning("PharmCAT timed out after 120 s")
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

    for gene_entry in data.get("genes", []):
        gene = gene_entry.get("gene", "")
        if not gene:
            continue

        # Diplotype
        src_dips = gene_entry.get("sourceDiplotypes", [])
        if src_dips:
            best = src_dips[0]
            diplotypes[gene] = best.get("name", "")
        rec_dip = gene_entry.get("recommendationDiplotype", "")
        if rec_dip and gene not in diplotypes:
            diplotypes[gene] = rec_dip

        # Phenotype
        pheno_list = gene_entry.get("phenotypes", [])
        if pheno_list:
            phenotypes[gene] = pheno_list[0].get("phenotype", "")

        # Drug recommendations
        for drug_entry in gene_entry.get("drugs", []):
            drug_name = drug_entry.get("name", drug_entry.get("drug", ""))
            guideline = drug_entry.get("guideline", "")
            classification = drug_entry.get("classification", "")
            drug_recs.append(
                {
                    "gene": gene,
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
