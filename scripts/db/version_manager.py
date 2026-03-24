"""Collect and manage database version metadata for reports."""
import os
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict
from scripts.common.config import get

logger = logging.getLogger(__name__)


def get_all_db_versions(skip_api: bool = False) -> Dict:
    """Collect version metadata from all data sources.
    Returns a structured dict for inclusion in report_data["db_versions"].
    """
    versions = {}

    # ClinVar local DB
    try:
        from scripts.db.query_local_clinvar import get_db_version as get_clinvar_version
        clinvar_meta = get_clinvar_version()
        if clinvar_meta.get("source") != "not available":
            versions["ClinVar"] = {
                "source": "local_db",
                "release": clinvar_meta.get("clinvar_release", "unknown"),
                "build_date": clinvar_meta.get("build_date", "unknown"),
                "variant_count": clinvar_meta.get("variant_count", "unknown"),
                "assembly": clinvar_meta.get("assembly", "GRCh38"),
            }
    except Exception:
        pass

    # If no local ClinVar and API was used
    if "ClinVar" not in versions and not skip_api:
        versions["ClinVar"] = {
            "source": "api",
            "release": f"E-utilities (queried {datetime.now(timezone.utc).strftime('%Y-%m-%d')})",
            "build_date": datetime.now(timezone.utc).strftime('%Y-%m-%d'),
            "assembly": "GRCh38",
        }
    elif "ClinVar" not in versions:
        versions["ClinVar"] = {"source": "not_available", "release": "N/A"}

    # gnomAD tabix (VCF.bgz) — checked first, takes priority over SQLite
    try:
        from scripts.db.query_tabix_gnomad import get_db_version as get_tabix_version
        tabix_meta = get_tabix_version()
        if tabix_meta.get("source") != "not_available":
            versions["gnomAD"] = tabix_meta
    except Exception:
        pass

    # gnomAD local DB (SQLite fallback)
    if "gnomAD" not in versions:
        try:
            from scripts.db.query_local_gnomad import get_db_version as get_gnomad_version
            gnomad_meta = get_gnomad_version()
            if gnomad_meta.get("source") != "not available":
                versions["gnomAD"] = {
                    "source": "local_db",
                    "version": gnomad_meta.get("gnomad_version", "unknown"),
                    "build_date": gnomad_meta.get("build_date", "unknown"),
                    "variant_count": gnomad_meta.get("variant_count", "unknown"),
                    "assembly": gnomad_meta.get("assembly", "GRCh38"),
                }
        except Exception:
            pass

    if "gnomAD" not in versions and not skip_api:
        versions["gnomAD"] = {
            "source": "api",
            "version": "v4.1 (r4) / v2.1 (r2_1 fallback)",
            "build_date": datetime.now(timezone.utc).strftime('%Y-%m-%d'),
            "assembly": "GRCh38/GRCh37",
        }
    elif "gnomAD" not in versions:
        versions["gnomAD"] = {"source": "not_available", "version": "N/A"}

    # KRGDB
    krgdb_path = get("paths.krgdb", "data/krgdb_freq.tsv")
    if Path(krgdb_path).exists():
        stat = Path(krgdb_path).stat()
        versions["KRGDB"] = {
            "source": "local_file",
            "path": str(krgdb_path),
            "modified": datetime.fromtimestamp(stat.st_mtime).strftime('%Y-%m-%d'),
            "size_bytes": stat.st_size,
        }
    else:
        versions["KRGDB"] = {"source": "not_available"}

    # ACMG Rules
    versions["ACMG"] = {
        "standard": "ACMG/AMP 2015 + ClinGen SVI updates",
        "source": get("paths.acmg_rules", "data/acmg_rules.json"),
    }

    # PGx
    versions["CPIC/PGx"] = {
        "genes_count": len(get("pgx.genes", [])),
        "source": get("paths.pgx_table", "data/korean_pgx_table.json"),
    }

    # Gene Knowledge
    gk_path = get("paths.gene_knowledge", "data/gene_knowledge.json")
    if Path(gk_path).exists():
        versions["Gene Knowledge"] = {
            "source": gk_path,
            "content_status": "ai-generated-with-references",
            "modified": datetime.fromtimestamp(Path(gk_path).stat().st_mtime).strftime('%Y-%m-%d'),
        }

    # Annotation source config
    versions["_annotation_source"] = get("annotation.source", "auto")

    return versions
