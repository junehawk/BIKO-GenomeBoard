"""Database query orchestration for per-variant annotation."""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

from scripts.enrichment.query_clinvar import query_clinvar
from scripts.common.config import get
from scripts.storage.query_local_clinvar import query_local_clinvar
from scripts.storage.query_tabix_gnomad import query_tabix_gnomad
from scripts.population.query_gnomad import query_gnomad
from scripts.population.query_kova import query_kova
from scripts.pharmacogenomics.korean_pgx import check_korean_pgx

logger = logging.getLogger(__name__)


def query_variant_databases(variant, skip_api: bool) -> dict:
    """Run ClinVar, gnomAD, KOVA, and PGx queries for a single variant.

    ClinVar and gnomAD are run in parallel (unless skip_api is True).
    Respects annotation.source config: "local", "api", or "auto" (local-first with API fallback).

    KOVA v7 (Korean Variant Archive) is the sole Korean population source.
    """
    clinvar_result = {"clinvar_significance": "Not Found", "acmg_codes": [], "api_available": False}
    gnomad_result = {"gnomad_all": None, "gnomad_eas": None, "api_available": False}
    kova_freq = None
    kova_homozygote = None
    pgx_result = None

    annotation_source = get("annotation.source", "auto")

    def _unpack_kova(result):
        if result is None:
            return None, None
        return result.get("kova_af"), result.get("kova_homozygote")

    if skip_api:
        try:
            clinvar_result = query_local_clinvar(variant)
        except Exception as e:
            logger.warning(f"Local ClinVar lookup failed for {variant.variant_id}: {e}")
        try:
            gnomad_result = query_tabix_gnomad(variant)
        except Exception as e:
            logger.warning(f"Local gnomAD lookup failed for {variant.variant_id}: {e}")
        try:
            kova_freq, kova_homozygote = _unpack_kova(query_kova(variant))
        except Exception as e:
            logger.warning(f"KOVA lookup failed for {variant.variant_id}: {e}")
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
                    return query_tabix_gnomad(variant)
                elif annotation_source == "api":
                    return query_gnomad(variant)
                else:  # auto: tabix first, API fallback
                    result = query_tabix_gnomad(variant)
                    if result["gnomad_all"] is None:
                        logger.debug(f"Local gnomAD miss for {variant.variant_id}, falling back to API")
                        result = query_gnomad(variant)
                    return result
            except Exception as e:
                logger.warning(f"gnomAD query failed for {variant.variant_id}: {e}")
                return {"gnomad_all": None, "gnomad_eas": None, "api_available": False}

        def _run_kova():
            try:
                return query_kova(variant)
            except Exception as e:
                logger.warning(f"KOVA lookup failed for {variant.variant_id}: {e}")
                return None

        def _run_pgx():
            try:
                return check_korean_pgx(variant)
            except Exception as e:
                logger.warning(f"PGx check failed for {variant.variant_id}: {e}")
                return None

        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = {
                executor.submit(_run_clinvar): "clinvar",
                executor.submit(_run_gnomad): "gnomad",
                executor.submit(_run_kova): "kova",
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
                    elif key == "kova":
                        kova_freq, kova_homozygote = _unpack_kova(result)
                    elif key == "pgx":
                        pgx_result = result
                except Exception as e:
                    logger.warning(f"{key} query raised unexpectedly for {variant.variant_id}: {e}")

    return {
        "clinvar": clinvar_result,
        "gnomad": gnomad_result,
        "kova_freq": kova_freq,
        "kova_homozygote": kova_homozygote,
        "pgx": pgx_result,
    }
