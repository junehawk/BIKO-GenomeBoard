"""Collect and manage database version metadata for reports."""

import json
import logging
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional

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
            "build_date": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
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
            "build_date": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
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
            "modified": datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d"),
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
            "modified": datetime.fromtimestamp(Path(gk_path).stat().st_mtime).strftime("%Y-%m-%d"),
        }

    # PM1 Hotspots — curated ACMG PM1 table (v2.2 A3)
    pm1_meta = _pm1_hotspots_version()
    if pm1_meta:
        versions["PM1_Hotspots"] = pm1_meta

    # DDG2P neurodevelopmental panel — full EBI Gene2Phenotype ingest (v2.3 T7)
    ddg2p_meta = _ddg2p_panel_version()
    if ddg2p_meta:
        versions["DDG2P"] = ddg2p_meta

    # CIViC — previously omitted despite having a metadata table (v2.2 bycatch C2-db-1)
    civic_meta = _civic_version()
    if civic_meta:
        versions["CIViC"] = civic_meta

    # cancerhotspots v2 single-residue TSV — informational; may be stubbed out
    ch_meta = _cancerhotspots_version()
    if ch_meta:
        versions["cancerhotspots_v2_single"] = ch_meta

    # gnomAD v4.1 gene constraint metrics — pLI / LOEUF / missense Z used by
    # the de novo carve-out admission OR-branch (v2.3-T8).
    try:
        from scripts.db.query_gnomad_constraint import get_db_version as get_gnomad_constraint_version

        gc_meta = get_gnomad_constraint_version()
        if gc_meta.get("source") != "not_available":
            versions["gnomAD_constraint"] = gc_meta
    except Exception:
        pass

    # Annotation source config
    versions["_annotation_source"] = get("annotation.source", "auto")

    return versions


# ---------------------------------------------------------------------------
# v2.2 helpers — individual source version lookups
# ---------------------------------------------------------------------------


def _pm1_hotspots_version() -> Optional[Dict]:
    """Read version metadata from data/pm1_hotspot_domains.json if present."""
    json_path = Path(get("paths.pm1_hotspots_json", "data/pm1_hotspot_domains.json"))
    if not json_path.exists():
        return None
    try:
        with json_path.open() as f:
            payload = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        logger.warning("PM1 hotspot JSON unreadable at %s: %s", json_path, e)
        return None

    return {
        "source": "local_json",
        "version": payload.get("version", "unknown"),
        "build_date": payload.get("build_date", "unknown"),
        "source_refs": payload.get("source_refs", []),
        "source_hash": payload.get("source_hash", ""),
        "record_count": payload.get("record_count", 0),
        "path": str(json_path),
    }


def _ddg2p_panel_version() -> Optional[Dict]:
    """Read version metadata from data/ddg2p_neurodev_genes.json if present.

    The file is built by ``scripts/tools/build_ddg2p_table.py`` from the EBI
    Gene2Phenotype FTP archive (CC0). Both the v1 hand-curated 30-gene starter
    set and the v2.3 full ~2200-gene ingest produce a JSON with the same
    top-level keys; this helper reads the build metadata so reports can
    surface DDG2P provenance alongside ClinVar/CIViC/etc.
    """
    json_path = Path(get("paths.ddg2p_panel_json", "data/ddg2p_neurodev_genes.json"))
    if not json_path.exists():
        return None
    try:
        with json_path.open() as f:
            payload = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        logger.warning("DDG2P panel JSON unreadable at %s: %s", json_path, e)
        return None

    genes = payload.get("genes") or {}
    return {
        "source": "local_json",
        "upstream": payload.get("source", "Gene2Phenotype DDG2P"),
        "license": payload.get("license", "CC0"),
        "build_date": payload.get("build_date", "unknown"),
        "source_url": payload.get("source_url", ""),
        "source_md5": payload.get("source_md5", ""),
        "gene_count": int(payload.get("record_count", len(genes)) or len(genes)),
        "admission_confidences": payload.get("admission_confidences", []),
        "path": str(json_path),
    }


def _civic_version() -> Optional[Dict]:
    """Read CIViC metadata table — previously missing from get_all_db_versions."""
    db_path = Path(get("paths.civic_db", "data/db/civic.sqlite3"))
    if not db_path.exists():
        return None
    try:
        conn = sqlite3.connect(str(db_path))
        try:
            rows = dict(conn.execute("SELECT key, value FROM metadata").fetchall())
        finally:
            conn.close()
    except sqlite3.Error as e:
        logger.warning("CIViC metadata read failed for %s: %s", db_path, e)
        return None

    return {
        "source": "local_db",
        "build_date": rows.get("build_date", "unknown"),
        "upstream": rows.get("source", "CIViC (civicdb.org)"),
        "gene_count": int(rows.get("gene_count", 0) or 0),
        "variant_count": int(rows.get("variant_count", 0) or 0),
        "evidence_count": int(rows.get("evidence_count", 0) or 0),
        "path": str(db_path),
    }


def _cancerhotspots_version() -> Optional[Dict]:
    """Surface cancerhotspots_v2_single TSV file status (may be a stub)."""
    tsv_path = Path(get("paths.cancerhotspots_tsv", "data/db/hotspots/hotspots_v2_single.tsv"))
    if not tsv_path.exists():
        return None
    stat = tsv_path.stat()
    size = stat.st_size
    # A non-stub TSV would be orders of magnitude larger than 1 KB.
    stubbed = size < 1024
    return {
        "source": "local_tsv" if not stubbed else "stub",
        "path": str(tsv_path),
        "size_bytes": size,
        "modified": datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d"),
        "note": "placeholder — upstream download not re-wired" if stubbed else "",
    }


def get_version(name: str, skip_api: bool = True) -> Optional[Dict]:
    """Return the version metadata dict for a single named source.

    Thin wrapper around `get_all_db_versions` for test ergonomics. Returns
    None when the source is unknown or unavailable.
    """
    versions = get_all_db_versions(skip_api=skip_api)
    return versions.get(name)
