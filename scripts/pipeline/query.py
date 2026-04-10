"""Database query orchestration for per-variant annotation."""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

from scripts.clinical.query_clinvar import query_clinvar
from scripts.db.query_local_clinvar import query_local_clinvar
from scripts.db.query_local_gnomad import query_local_gnomad
from scripts.db.query_tabix_gnomad import query_tabix_gnomad
from scripts.korean_pop.query_gnomad import query_gnomad
from scripts.korean_pop.query_krgdb import query_krgdb
from scripts.korean_pop.query_korea4k import query_korea4k
from scripts.korean_pop.query_nard2 import query_nard2
from scripts.pharma.korean_pgx import check_korean_pgx
from scripts.common.config import get

logger = logging.getLogger(__name__)


def query_variant_databases(variant, krgdb_path: str, skip_api: bool) -> dict:
    """Run ClinVar, gnomAD, KRGDB, Korea4K, NARD2, and PGx queries for a single variant.

    ClinVar and gnomAD are run in parallel (unless skip_api is True).
    Respects annotation.source config: "local", "api", or "auto" (local-first with API fallback).
    """
    clinvar_result = {"clinvar_significance": "Not Found", "acmg_codes": [], "api_available": False}
    gnomad_result = {"gnomad_all": None, "gnomad_eas": None, "api_available": False}
    krgdb_freq = None
    korea4k_freq = None
    nard2_freq = None
    pgx_result = None

    annotation_source = get("annotation.source", "auto")

    if skip_api:
        try:
            clinvar_result = query_local_clinvar(variant)
        except Exception as e:
            logger.warning(f"Local ClinVar lookup failed for {variant.variant_id}: {e}")
        try:
            gnomad_result = query_tabix_gnomad(variant)
            if gnomad_result["gnomad_all"] is None:
                gnomad_result = query_local_gnomad(variant)
        except Exception as e:
            logger.warning(f"Local gnomAD lookup failed for {variant.variant_id}: {e}")
        try:
            krgdb_freq = query_krgdb(variant, krgdb_path)
        except Exception as e:
            logger.warning(f"KRGDB lookup failed for {variant.variant_id}: {e}")
        try:
            korea4k_freq = query_korea4k(variant)
        except Exception as e:
            logger.warning(f"Korea4K lookup failed for {variant.variant_id}: {e}")
        try:
            nard2_freq = query_nard2(variant)
        except Exception as e:
            logger.warning(f"NARD2 lookup failed for {variant.variant_id}: {e}")
        try:
            pgx_result = check_korean_pgx(variant)
        except Exception as e:
            logger.warning(f"PGx check failed for {variant.variant_id}: {e}")
    else:
        def _run_clinvar():
            try:
                if annotation_source == "local":
                    return query_local_clinvar(variant)
                elif annotation_source == "api":
                    return query_clinvar(variant)
                else:  # auto: local first, API fallback
                    result = query_local_clinvar(variant)
                    if result["clinvar_significance"] == "Not Found":
                        logger.debug(f"Local ClinVar miss for {variant.variant_id}, falling back to API")
                        result = query_clinvar(variant)
                    return result
            except Exception as e:
                logger.warning(f"ClinVar query failed for {variant.variant_id}: {e}")
                return {"clinvar_significance": "Not Found", "acmg_codes": [], "api_available": False}

        def _run_gnomad():
            try:
                if annotation_source == "local":
                    result = query_tabix_gnomad(variant)
                    if result["gnomad_all"] is None:
                        result = query_local_gnomad(variant)
                    return result
                elif annotation_source == "api":
                    return query_gnomad(variant)
                else:  # auto: tabix first, then SQLite, then API fallback
                    result = query_tabix_gnomad(variant)
                    if result["gnomad_all"] is None:
                        result = query_local_gnomad(variant)
                    if result["gnomad_all"] is None:
                        logger.debug(f"Local gnomAD miss for {variant.variant_id}, falling back to API")
                        result = query_gnomad(variant)
                    return result
            except Exception as e:
                logger.warning(f"gnomAD query failed for {variant.variant_id}: {e}")
                return {"gnomad_all": None, "gnomad_eas": None, "api_available": False}

        def _run_krgdb():
            try:
                return query_krgdb(variant, krgdb_path)
            except Exception as e:
                logger.warning(f"KRGDB lookup failed for {variant.variant_id}: {e}")
                return None

        def _run_pgx():
            try:
                return check_korean_pgx(variant)
            except Exception as e:
                logger.warning(f"PGx check failed for {variant.variant_id}: {e}")
                return None

        def _run_korea4k():
            try:
                return query_korea4k(variant)
            except Exception as e:
                logger.warning(f"Korea4K lookup failed for {variant.variant_id}: {e}")
                return None

        def _run_nard2():
            try:
                return query_nard2(variant)
            except Exception as e:
                logger.warning(f"NARD2 lookup failed for {variant.variant_id}: {e}")
                return None

        with ThreadPoolExecutor(max_workers=6) as executor:
            futures = {
                executor.submit(_run_clinvar): "clinvar",
                executor.submit(_run_gnomad): "gnomad",
                executor.submit(_run_krgdb): "krgdb",
                executor.submit(_run_korea4k): "korea4k",
                executor.submit(_run_nard2): "nard2",
                executor.submit(_run_pgx): "pgx",
            }
            for future in as_completed(futures):
                key = futures[future]
                try:
                    result = future.result()
                    if key == "clinvar":
                        clinvar_result = result
                    elif key == "gnomad":
                        gnomad_result = result
                    elif key == "krgdb":
                        krgdb_freq = result
                    elif key == "korea4k":
                        korea4k_freq = result
                    elif key == "nard2":
                        nard2_freq = result
                    elif key == "pgx":
                        pgx_result = result
                except Exception as e:
                    logger.warning(f"{key} query raised unexpectedly for {variant.variant_id}: {e}")

    return {
        "clinvar": clinvar_result,
        "gnomad": gnomad_result,
        "krgdb_freq": krgdb_freq,
        "korea4k_freq": korea4k_freq,
        "nard2_freq": nard2_freq,
        "pgx": pgx_result,
    }
