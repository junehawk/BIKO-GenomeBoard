#!/usr/bin/env python3
"""GenomeBoard Standalone Pipeline — VCF → Report"""

import argparse
import csv
import json
import logging
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from datetime import date
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.intake.parse_vcf import parse_vcf
from scripts.clinical.query_clinvar import query_clinvar
from scripts.db.query_local_clinvar import query_local_clinvar
from scripts.db.query_local_gnomad import query_local_gnomad
from scripts.db.query_tabix_gnomad import query_tabix_gnomad
from scripts.db.version_manager import get_all_db_versions
from scripts.clinical.hpo_matcher import resolve_hpo_terms, calculate_hpo_score, get_matching_hpo_terms
from scripts.clinical.query_omim import query_omim
from scripts.clinical.query_clingen import get_gene_validity
from scripts.korean_pop.query_gnomad import query_gnomad
from scripts.korean_pop.query_krgdb import query_krgdb
from scripts.korean_pop.compare_freq import compare_frequencies
from scripts.pharma.korean_pgx import check_korean_pgx
from scripts.classification.acmg_engine import classify_variant, check_clinvar_conflict, apply_clinvar_override
from scripts.counselor.generate_pdf import generate_report_html, generate_pdf
from scripts.common.models import AcmgEvidence, FrequencyData
from scripts.common.config import get
from scripts.clinical.oncokb import assign_tier, get_cancer_gene_info, get_tier_label

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("genomeboard")


def _progress(msg: str) -> None:
    """Print a progress message to stderr."""
    print(msg, file=sys.stderr, flush=True)


def _query_variant_databases(variant, krgdb_path: str, skip_api: bool) -> dict:
    """Run ClinVar, gnomAD, KRGDB, and PGx queries for a single variant.
    ClinVar and gnomAD are run in parallel (unless skip_api is True).
    Respects annotation.source config: "local", "api", or "auto" (local-first with API fallback).
    """
    clinvar_result = {"clinvar_significance": "Not Found", "acmg_codes": [], "api_available": False}
    gnomad_result = {"gnomad_all": None, "gnomad_eas": None, "api_available": False}
    krgdb_freq = None
    pgx_result = None

    annotation_source = get("annotation.source", "auto")

    if skip_api:
        # Still run KRGDB (local) and PGx (local)
        # Also attempt local ClinVar and gnomAD DBs regardless of annotation.source
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
            pgx_result = check_korean_pgx(variant)
        except Exception as e:
            logger.warning(f"PGx check failed for {variant.variant_id}: {e}")
    else:
        # Determine ClinVar query function(s) based on annotation.source
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
                    # Try tabix first, then SQLite fallback
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

        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = {
                executor.submit(_run_clinvar): "clinvar",
                executor.submit(_run_gnomad): "gnomad",
                executor.submit(_run_krgdb): "krgdb",
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
                    elif key == "pgx":
                        pgx_result = result
                except Exception as e:
                    logger.warning(f"{key} query raised unexpectedly for {variant.variant_id}: {e}")

    return {
        "clinvar": clinvar_result,
        "gnomad": gnomad_result,
        "krgdb_freq": krgdb_freq,
        "pgx": pgx_result,
    }


def _build_variant_records(variants, db_results, freq_results, classification_results,
                           mode, hpo_results):
    """Build the variant_records list from per-variant results.

    This is the shared logic used by both run_pipeline() and _assemble_sample_report()
    to avoid duplication.
    """
    variant_records = []
    for variant in variants:
        db = db_results[variant.variant_id]
        freq = freq_results[variant.variant_id]
        classification = classification_results[variant.variant_id]
        clinvar = db["clinvar"]

        variant_records.append(
            {
                "variant": variant.variant_id,
                "gene": variant.gene,
                "chrom": variant.chrom,
                "pos": variant.pos,
                "ref": variant.ref,
                "alt": variant.alt,
                "classification": classification.classification,
                "acmg_codes": classification.evidence_codes,
                "conflict": classification.conflict,
                "clinvar_override": getattr(classification, "clinvar_override", False),
                "original_engine_classification": getattr(classification, "original_engine_classification", None),
                "clinvar_significance": clinvar.get("clinvar_significance", "Not Found"),
                "clinvar_id": clinvar.get("clinvar_id"),
                "review_status": clinvar.get("review_status"),
                "gnomad_all": db["gnomad"].get("gnomad_all"),
                "gnomad_eas": db["gnomad"].get("gnomad_eas"),
                "krgdb_freq": db["krgdb_freq"],
                "korean_flag": freq.get("korean_flag", ""),
                # Annotation fields (from pre-annotated VCF via VEP/SnpEff)
                "hgvsc": variant.hgvsc or "",
                "hgvsp": variant.hgvsp or "",
                "consequence": variant.consequence or "",
                "transcript": variant.transcript or "",
                "impact": variant.impact or "",
                "sift": variant.sift or "",
                "polyphen": variant.polyphen or "",
                "agents": {
                    "clinical": clinvar,
                    "korean_pop": {
                        "gnomad_all": db["gnomad"].get("gnomad_all"),
                        "gnomad_eas": db["gnomad"].get("gnomad_eas"),
                        "krgdb_freq": db["krgdb_freq"],
                        "korean_flag": freq.get("korean_flag", ""),
                        "api_available": db["gnomad"].get("api_available", False),
                    },
                },
            }
        )

    # Assign OncoKB tiers to all variants
    for v_result in variant_records:
        gene = v_result.get("gene", "")
        cls = v_result.get("classification", "VUS")
        clinvar_sig = v_result.get("clinvar_significance", "")
        tier = assign_tier(cls, gene, clinvar_sig)
        v_result["tier"] = tier
        v_result["tier_label"] = get_tier_label(tier)
        cancer_info = get_cancer_gene_info(gene)
        if cancer_info:
            v_result["cancer_gene_type"] = cancer_info.get("type", "")
            v_result["oncokb_level"] = cancer_info.get("level", "")
        else:
            v_result["cancer_gene_type"] = ""
            v_result["oncokb_level"] = ""

    # For rare disease mode, add HPO score and OMIM/ClinGen data to each variant
    if mode == "rare-disease":
        for v_result in variant_records:
            gene = v_result.get("gene", "")
            # HPO phenotype matching
            v_result["hpo_score"] = calculate_hpo_score(gene, hpo_results)
            v_result["matching_hpo"] = get_matching_hpo_terms(gene, hpo_results)
            # OMIM
            omim = query_omim(gene)
            if omim:
                v_result["omim_mim"] = omim.get("mim", "")
                v_result["omim_phenotypes"] = omim.get("phenotypes", [])
                v_result["inheritance"] = omim.get("inheritance", "")
            else:
                v_result["omim_mim"] = ""
                v_result["omim_phenotypes"] = []
                v_result["inheritance"] = ""
            # ClinGen
            v_result["clingen_validity"] = get_gene_validity(gene) or ""

        # Sort by: classification rank (Pathogenic first) then HPO score (descending)
        cls_rank = {
            "Pathogenic": 0,
            "Likely Pathogenic": 1,
            "VUS": 2,
            "Drug Response": 3,
            "Risk Factor": 4,
            "Likely Benign": 5,
            "Benign": 6,
        }
        variant_records.sort(
            key=lambda v: (cls_rank.get(v.get("classification", "VUS"), 2), -v.get("hpo_score", 0))
        )

    return variant_records


def _build_summary(variant_records):
    """Build summary counts dict from variant records."""
    total = len(variant_records)
    return {
        "total": total,
        "pathogenic": sum(1 for r in variant_records if r["classification"] == "Pathogenic"),
        "likely_pathogenic": sum(1 for r in variant_records if r["classification"] == "Likely Pathogenic"),
        "drug_response": sum(1 for r in variant_records if r["classification"] == "Drug Response"),
        "risk_factor": sum(1 for r in variant_records if r["classification"] == "Risk Factor"),
        "vus": sum(1 for r in variant_records if r["classification"] == "VUS"),
        "likely_benign": sum(1 for r in variant_records if r["classification"] == "Likely Benign"),
        "benign": sum(1 for r in variant_records if r["classification"] == "Benign"),
    }


def _classify_variants(variants, db_results, freq_results):
    """Run frequency comparison and ACMG classification for a list of variants.

    Returns (freq_results, classification_results) dicts keyed by variant_id.
    freq_results may be pre-populated; if provided, it is augmented in place
    for any missing keys (batch reuse case).
    """
    classification_results = {}
    for variant in variants:
        db = db_results[variant.variant_id]
        freq = freq_results[variant.variant_id]

        evidences = []
        for code in db["clinvar"].get("acmg_codes", []):
            evidences.append(AcmgEvidence(code=code, source="clinvar", description=""))
        for code in freq.get("acmg_codes", []):
            evidences.append(AcmgEvidence(code=code, source="freq_comparison", description=""))

        classification = classify_variant(evidences, gene=variant.gene)
        clinvar_sig = db["clinvar"].get("clinvar_significance", "Not Found")
        review_status = db["clinvar"].get("review_status", "")
        conflict = check_clinvar_conflict(classification.classification, clinvar_sig)
        classification.conflict = conflict

        # Apply ClinVar override when evidence is high-confidence
        original_classification = classification.classification
        final_classification = apply_clinvar_override(original_classification, clinvar_sig, review_status)
        if final_classification != original_classification:
            classification.classification = final_classification
            classification.clinvar_override = True
            classification.original_engine_classification = original_classification
        else:
            classification.clinvar_override = False
            classification.original_engine_classification = None

        classification_results[variant.variant_id] = classification

    return classification_results


def run_pipeline(
    vcf_path: str,
    output_path: str = "output/report.html",
    krgdb_path: str = "data/krgdb_freq.tsv",
    sample_id: str = None,
    json_output: str = None,
    skip_api: bool = False,
    mode: str = None,
    hpo_ids: list = None,
    hide_vus: bool = False,
) -> dict:
    """Run the full GenomeBoard analysis pipeline.

    Returns the assembled report data dict (or None on fatal error).
    """
    mode = mode or get("report.default_mode", "cancer")
    start_time = time.time()

    # HPO processing (rare disease mode) — always resolve, even in skip_api mode
    # HPO phenotype matching is essential for candidate ranking
    hpo_results = []
    if mode == "rare-disease" and hpo_ids:
        logger.info("[HPO] Resolving phenotype terms...")
        hpo_results = resolve_hpo_terms(hpo_ids)
        for h in hpo_results:
            logger.info(f"  → {h['id']}: {h['name']} ({len(h['genes'])} associated genes)")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Derive sample_id from filename if not provided
    if sample_id is None:
        sample_id = Path(vcf_path).stem.upper()

    # ── Step 1: Parse VCF ────────────────────────────────────────────────────
    _progress(f"[1/6] Parsing VCF: {vcf_path}")
    vcf_file = Path(vcf_path)
    if not vcf_file.exists():
        logger.error(f"VCF file not found: {vcf_path}")
        return None

    variants = parse_vcf(str(vcf_file))
    if not variants:
        logger.error(f"No variants parsed from {vcf_path}")
        return None
    _progress(f"  → Found {len(variants)} variants")

    # Check KRGDB availability
    krgdb_file = Path(krgdb_path)
    if not krgdb_file.exists():
        logger.warning(f"KRGDB file not found: {krgdb_path} — continuing without Korean frequency data")

    # ── Step 2: Query databases per variant ──────────────────────────────────
    _progress("[2/6] Querying databases (ClinVar, gnomAD, KRGDB)...")
    if skip_api:
        _progress("  [skip-api mode: ClinVar and gnomAD queries disabled]")

    db_results = {}
    for variant in variants:
        result = _query_variant_databases(variant, str(krgdb_file), skip_api)
        db_results[variant.variant_id] = result

        clinvar_sig = result["clinvar"]["clinvar_significance"] if not skip_api else "N/A (offline)"
        gnomad_all = result["gnomad"].get("gnomad_all")
        krgdb_val = result["krgdb_freq"]

        gnomad_str = f"{gnomad_all:.5f}" if gnomad_all is not None else "N/A"
        krgdb_str = f"{krgdb_val:.5f}" if krgdb_val is not None else "N/A"
        gene_label = variant.gene or variant.variant_id
        _progress(f"  → {gene_label}: ClinVar={clinvar_sig}, gnomAD={gnomad_str}, KRGDB={krgdb_str}")

    # ── Step 3: Frequency comparison ─────────────────────────────────────────
    _progress("[3/6] Running frequency comparison...")
    freq_results = {}
    for variant in variants:
        db = db_results[variant.variant_id]
        freq_data = FrequencyData(
            krgdb=db["krgdb_freq"],
            gnomad_eas=db["gnomad"].get("gnomad_eas"),
            gnomad_all=db["gnomad"].get("gnomad_all"),
        )
        freq_results[variant.variant_id] = compare_frequencies(freq_data)

    # ── Step 4: Pharmacogenomics ──────────────────────────────────────────────
    _progress("[4/6] Checking pharmacogenomics (PGx)...")
    pgx_hits = []
    for variant in variants:
        pgx = db_results[variant.variant_id]["pgx"]
        if pgx:
            pgx_hits.append(pgx)
            cpic_str = getattr(pgx, "cpic_level", "")
            _progress(f"  → {pgx.gene}: {pgx.phenotype}, CPIC Level {cpic_str}")

    # ── Step 5: ACMG classification ───────────────────────────────────────────
    _progress("[5/6] Running ACMG classification engine...")
    classification_results = _classify_variants(variants, db_results, freq_results)
    for variant in variants:
        classification = classification_results[variant.variant_id]
        codes_str = "+".join(classification.evidence_codes) if classification.evidence_codes else "no codes"
        gene_label = variant.gene or variant.variant_id
        _progress(f"  → {gene_label}: {classification.classification} ({codes_str})")

    # ── Step 6: Assemble report data ──────────────────────────────────────────
    variant_records = _build_variant_records(
        variants, db_results, freq_results, classification_results, mode, hpo_results
    )
    summary = _build_summary(variant_records)

    # Split variants for template rendering (VUS filtering)
    # Tier lists are always populated for cancer mode (drive summary-page display).
    # detailed_variants drives the per-variant detail pages and follows hide_vus logic.
    _vus_classes = ("VUS", "Benign", "Likely Benign")
    tier1_variants = [v for v in variant_records if v.get("tier") == 1]
    tier2_variants = [v for v in variant_records if v.get("tier") == 2]
    tier3_variants = [v for v in variant_records if v.get("tier") == 3]
    tier4_count = sum(1 for v in variant_records if v.get("tier") == 4)

    if hide_vus:
        detailed_variants = [v for v in variant_records if v["classification"] not in _vus_classes]
        omitted_variants = [v for v in variant_records if v["classification"] in _vus_classes]
    else:
        detailed_variants = variant_records
        omitted_variants = []

    report_data = {
        "sample_id": sample_id,
        "date": str(date.today()),
        "variants": variant_records,
        "detailed_variants": detailed_variants,
        "omitted_variants": omitted_variants,
        "hide_vus": hide_vus,
        "tier1_variants": tier1_variants,
        "tier2_variants": tier2_variants,
        "tier3_variants": tier3_variants,
        "tier4_count": tier4_count,
        "pgx_results": [
            {
                "gene": p.gene,
                "star_allele": p.star_allele,
                "phenotype": p.phenotype,
                "cpic_level": p.cpic_level,
                "korean_prevalence": p.korean_prevalence,
                "western_prevalence": p.western_prevalence,
                "clinical_impact": p.clinical_impact,
                "cpic_recommendation": p.cpic_recommendation,
                "korean_flag": p.korean_flag,
            }
            for p in pgx_hits
        ],
        "summary": summary,
        "db_versions": get_all_db_versions(skip_api=skip_api),
        "pipeline": {
            "skip_api": skip_api,
            "krgdb_path": str(krgdb_file),
        },
        "mode": mode,
        "hpo_results": hpo_results,
    }

    # ── Step 6: Generate report ────────────────────────────────────────────────
    _progress("[6/6] Generating report...")
    output_str = str(output_path)
    if output_str.endswith(".pdf"):
        actual_path = generate_pdf(report_data, output_str, mode=mode)
        final_path = Path(actual_path)
    else:
        html = generate_report_html(report_data, mode=mode)
        output_path.write_text(html, encoding="utf-8")
        final_path = output_path

    size_kb = final_path.stat().st_size // 1024
    suffix = final_path.suffix.upper().lstrip(".")
    _progress(f"  → {suffix}: {final_path} ({size_kb}KB)")

    # Optional JSON output
    if json_output:
        json_path = Path(json_output)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        # Serialize FrequencyData objects that may have slipped in
        json_path.write_text(
            json.dumps(report_data, indent=2, default=str, ensure_ascii=False),
            encoding="utf-8",
        )
        _progress(f"  → JSON: {json_path}")

    elapsed = time.time() - start_time
    _progress(f"\nDone! {summary['total']} variants analyzed in {elapsed:.1f}s")
    _progress(
        f"  Pathogenic: {summary['pathogenic']} | "
        f"Likely Pathogenic: {summary['likely_pathogenic']} | "
        f"Drug Response: {summary['drug_response']} | "
        f"VUS: {summary['vus']} | "
        f"Benign: {summary['benign'] + summary['likely_benign']}"
    )

    return report_data


# ── Batch pipeline helpers ─────────────────────────────────────────────────────

def _discover_samples(batch_path: str) -> list:
    """Discover VCF files from directory or manifest CSV.

    Returns: [{"sample_id": "S001", "vcf_path": "/path/to/s001.vcf"}, ...]
    """
    path = Path(batch_path)
    if path.is_dir():
        samples = []
        vcf_files = sorted(path.glob("*.vcf")) + sorted(path.glob("*.vcf.gz"))
        for vcf_file in vcf_files:
            # Strip .vcf or .vcf.gz suffix to get sample_id
            sample_id = vcf_file.name
            if sample_id.endswith(".vcf.gz"):
                sample_id = sample_id[:-7]
            elif sample_id.endswith(".vcf"):
                sample_id = sample_id[:-4]
            samples.append({"sample_id": sample_id, "vcf_path": str(vcf_file)})
        return samples
    elif path.suffix == ".csv":
        samples = []
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                samples.append({"sample_id": row["sample_id"], "vcf_path": row["vcf_path"]})
        return samples
    else:
        raise ValueError(f"Batch path must be a directory or CSV file: {batch_path}")


def _collect_unique_variants(samples: list) -> tuple:
    """Parse all VCFs, return (unique_variants_dict, sample_variant_map).

    unique_variants: {variant_key: Variant}
    sample_variant_map: {sample_id: [variant_key, ...]}
    """
    unique_variants = {}  # key -> Variant
    sample_map = {}       # sample_id -> [keys]

    for sample in samples:
        variants = parse_vcf(sample["vcf_path"])
        keys = []
        for v in variants:
            key = f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}"
            if key not in unique_variants:
                unique_variants[key] = v
            keys.append(key)
        sample_map[sample["sample_id"]] = keys

    return unique_variants, sample_map


def _bulk_annotate_variants(unique_variants: dict, krgdb_path: str, skip_api: bool,
                            max_workers: int = 8) -> dict:
    """Annotate all unique variants. Returns {variant_key: annotation_dict}."""
    annotations = {}
    rate_limiter = threading.Semaphore(10)  # Max 10 concurrent API calls

    def _annotate_one(key, variant):
        with rate_limiter:
            return key, _query_variant_databases(variant, krgdb_path, skip_api)

    _progress(f"[2/5] Annotating {len(unique_variants):,} unique variants ({max_workers} workers)...")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_annotate_one, k, v): k for k, v in unique_variants.items()}
        done = 0
        for future in as_completed(futures):
            key, result = future.result()
            annotations[key] = result
            done += 1
            if done % 100 == 0:
                _progress(f"  → {done:,}/{len(unique_variants):,} variants annotated")

    return annotations


def _assemble_sample_report(sample, variant_keys, annotations, unique_variants,
                            mode, hpo_ids, hpo_results, krgdb_path, hide_vus=False):
    """Build report_data dict for one sample using pre-computed annotations.

    Reuses the same ACMG classification and frequency-comparison logic as run_pipeline().
    """
    # Reconstruct the variant list for this sample (in order)
    variants = [unique_variants[k] for k in variant_keys]

    # Build db_results from cached annotations (keyed by variant_id)
    db_results = {}
    for key, variant in zip(variant_keys, variants):
        db_results[variant.variant_id] = annotations[key]

    # Frequency comparison
    freq_results = {}
    for variant in variants:
        db = db_results[variant.variant_id]
        freq_data = FrequencyData(
            krgdb=db["krgdb_freq"],
            gnomad_eas=db["gnomad"].get("gnomad_eas"),
            gnomad_all=db["gnomad"].get("gnomad_all"),
        )
        freq_results[variant.variant_id] = compare_frequencies(freq_data)

    # PGx hits
    pgx_hits = []
    for variant in variants:
        pgx = db_results[variant.variant_id]["pgx"]
        if pgx:
            pgx_hits.append(pgx)

    # ACMG classification (shared logic)
    classification_results = _classify_variants(variants, db_results, freq_results)

    # Build variant records (shared logic, including rare-disease HPO enrichment)
    variant_records = _build_variant_records(
        variants, db_results, freq_results, classification_results, mode, hpo_results
    )
    summary = _build_summary(variant_records)

    # Split variants for template rendering (VUS filtering)
    _vus_classes = ("VUS", "Benign", "Likely Benign")
    tier1_variants = [v for v in variant_records if v.get("tier") == 1]
    tier2_variants = [v for v in variant_records if v.get("tier") == 2]
    tier3_variants = [v for v in variant_records if v.get("tier") == 3]
    tier4_count = sum(1 for v in variant_records if v.get("tier") == 4)

    if hide_vus:
        detailed_variants = [v for v in variant_records if v["classification"] not in _vus_classes]
        omitted_variants = [v for v in variant_records if v["classification"] in _vus_classes]
    else:
        detailed_variants = variant_records
        omitted_variants = []

    return {
        "sample_id": sample["sample_id"],
        "date": str(date.today()),
        "variants": variant_records,
        "detailed_variants": detailed_variants,
        "omitted_variants": omitted_variants,
        "hide_vus": hide_vus,
        "tier1_variants": tier1_variants,
        "tier2_variants": tier2_variants,
        "tier3_variants": tier3_variants,
        "tier4_count": tier4_count,
        "pgx_results": [
            {
                "gene": p.gene,
                "star_allele": p.star_allele,
                "phenotype": p.phenotype,
                "cpic_level": p.cpic_level,
                "korean_prevalence": p.korean_prevalence,
                "western_prevalence": p.western_prevalence,
                "clinical_impact": p.clinical_impact,
                "cpic_recommendation": p.cpic_recommendation,
                "korean_flag": p.korean_flag,
            }
            for p in pgx_hits
        ],
        "summary": summary,
        "db_versions": get_all_db_versions(skip_api=False),
        "pipeline": {
            "skip_api": False,
            "krgdb_path": str(krgdb_path),
        },
        "mode": mode,
        "hpo_results": hpo_results,
    }


def _generate_single_report(report_data: dict, output_path: str, mode: str) -> str:
    """Generate a single HTML report (runs in subprocess for PDF parallelism)."""
    # Re-import needed because this may run in a subprocess via ProcessPoolExecutor
    import sys
    from pathlib import Path as _Path
    sys.path.insert(0, str(_Path(__file__).parent.parent))
    from scripts.counselor.generate_pdf import generate_report_html as _gen_html

    html = _gen_html(report_data, mode=mode)
    _Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    _Path(output_path).write_text(html, encoding="utf-8")
    return output_path


def _generate_reports_parallel(sample_reports: list, output_dir: str, mode: str,
                               workers: int) -> list:
    """Generate HTML reports using ProcessPoolExecutor."""
    _progress(f"[5/5] Generating {len(sample_reports):,} reports ({workers} workers)...")

    results = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {}
        for sample_data in sample_reports:
            output_path = Path(output_dir) / f"{sample_data['sample_id']}_report.html"
            futures[
                executor.submit(_generate_single_report, sample_data, str(output_path), mode)
            ] = sample_data["sample_id"]

        for future in as_completed(futures):
            sample_id = futures[future]
            try:
                path = future.result()
                results.append(path)
            except Exception as e:
                logger.error(f"Report generation failed for {sample_id}: {e}")

    return results


def run_batch_pipeline(
    batch_path: str,
    output_dir: str = "output/batch",
    mode: str = "cancer",
    workers: int = None,
    skip_api: bool = False,
    krgdb_path: str = None,
    hpo_ids: list = None,
    hide_vus: bool = False,
    _skip_reports: bool = False,
) -> dict:
    """Process multiple VCF files with variant deduplication.

    Returns: {
        "samples_processed": N,
        "total_variants": N,
        "unique_variants": N,
        "cache_hits": N,
        "reports_generated": ["path1.html", ...],
        "errors": [{"sample": "...", "error": "..."}],
        "elapsed_seconds": N
    }
    """
    import time as time_mod

    workers = workers or os.cpu_count() or 4
    start = time_mod.time()

    # Step 1: Discover samples
    _progress(f"[1/5] Discovering samples from {batch_path}...")
    samples = _discover_samples(batch_path)
    _progress(f"  → Found {len(samples)} samples")

    if not samples:
        return {
            "samples_processed": 0,
            "total_variants": 0,
            "unique_variants": 0,
            "cache_hits": 0,
            "reports_generated": [],
            "errors": [],
            "elapsed_seconds": time_mod.time() - start,
        }

    # Step 2: Collect unique variants across all samples
    unique_variants, sample_map = _collect_unique_variants(samples)
    total_variants = sum(len(keys) for keys in sample_map.values())
    cache_hits = total_variants - len(unique_variants)
    _progress(f"  → {total_variants:,} total variants, {len(unique_variants):,} unique "
              f"({cache_hits:,} deduplicated)")

    # Step 3: Bulk annotate unique variants
    krgdb_file = krgdb_path or get("paths.krgdb", "data/krgdb_freq.tsv")
    annotations = _bulk_annotate_variants(
        unique_variants, krgdb_file, skip_api,
        max_workers=min(workers, 10)
    )

    # Step 4: HPO (rare disease mode)
    hpo_results = []
    if mode == "rare-disease" and hpo_ids:
        hpo_results = resolve_hpo_terms(hpo_ids)

    # Step 5: Assemble per-sample reports
    _progress(f"[3/5] Assembling {len(samples)} sample reports...")
    sample_reports = []
    errors = []
    for sample in samples:
        try:
            report_data = _assemble_sample_report(
                sample,
                sample_map[sample["sample_id"]],
                annotations,
                unique_variants,
                mode,
                hpo_ids,
                hpo_results,
                krgdb_file,
                hide_vus=hide_vus,
            )
            sample_reports.append(report_data)
        except Exception as e:
            logger.error(f"Failed to assemble report for {sample['sample_id']}: {e}")
            errors.append({"sample": sample["sample_id"], "error": str(e)})

    _progress(f"[4/5] Assembled {len(sample_reports)} reports ({len(errors)} errors)")

    # Step 6: Generate reports in parallel (skipped in test-only mode)
    report_paths = []
    if not _skip_reports and output_dir is not None:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        report_paths = _generate_reports_parallel(sample_reports, output_dir, mode, workers)

    elapsed = time_mod.time() - start
    _progress(f"\nBatch complete! {len(report_paths)} reports in {elapsed:.1f}s")
    _progress(f"  Unique variants: {len(unique_variants):,} / Total: {total_variants:,}")
    if errors:
        _progress(f"  Errors: {len(errors)}")

    return {
        "samples_processed": len(samples),
        "total_variants": total_variants,
        "unique_variants": len(unique_variants),
        "cache_hits": cache_hits,
        "reports_generated": report_paths,
        "errors": errors,
        "elapsed_seconds": elapsed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="GenomeBoard — Korean Population-Aware Genomic Variant Interpretation Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
REPORT MODES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  cancer (default)    FoundationOne CDx-style somatic variant report.
                      Includes: ACMG classification, PGx (12 genes),
                      Korean population comparison, therapeutic context.

  rare-disease        Germline rare disease report with phenotype matching.
                      Includes: HPO phenotype scoring, OMIM gene-disease
                      associations, ClinGen validity, inheritance patterns,
                      candidate gene ranking. Use with --hpo flag.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
DATA SOURCES (config.yaml → annotation.source)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  auto (default)      Local SQLite DB first, API fallback if not found.
  local               Local DB only (on-premise, no internet required).
  api                 ClinVar/gnomAD API only (requires internet).

  To use local databases, build them first:
    python scripts/db/build_clinvar_db.py data/db/variant_summary.txt.gz
    python scripts/db/build_gnomad_db.py data/db/gnomad_vcfs/*.vcf.bgz

  Download sources:
    ClinVar:  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
    gnomAD:   https://gnomad.broadinstitute.org/downloads#v4 (exomes sites VCF)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
EXAMPLES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # Single sample (cancer mode, API)
  python scripts/orchestrate.py sample.vcf -o report.html --json

  # Single sample (offline, local DB only)
  python scripts/orchestrate.py sample.vcf -o report.html --skip-api

  # Rare disease mode with HPO phenotypes
  python scripts/orchestrate.py patient.vcf --mode rare-disease \\
    --hpo HP:0001250,HP:0001263 -o report.html

  # PDF output
  python scripts/orchestrate.py sample.vcf -o report.pdf

  # Batch processing (directory of VCFs)
  python scripts/orchestrate.py --batch vcf_dir/ --output-dir output/batch --workers 8

  # Batch processing (manifest CSV: sample_id,vcf_path)
  python scripts/orchestrate.py --batch manifest.csv --output-dir output/batch

  # Build local databases for on-premise deployment
  python scripts/db/build_clinvar_db.py data/db/variant_summary.txt.gz
  python scripts/db/build_gnomad_db.py chr1.vcf.bgz chr2.vcf.bgz --version 4.1

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
DOCKER
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  docker build -t genomeboard .
  docker run -v ./input:/app/input -v ./output:/app/output \\
    genomeboard /app/input/sample.vcf -o /app/output/report.html
        """,
    )
    parser.add_argument(
        "vcf_path",
        nargs="?",
        default=None,
        help="Path to input VCF file (omit when using --batch)",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="output/report.html",
        help="Output file path (default: output/report.html). Use .pdf extension for PDF output.",
    )
    parser.add_argument(
        "--krgdb",
        default="data/krgdb_freq.tsv",
        help="KRGDB frequency data file (default: data/krgdb_freq.tsv)",
    )
    parser.add_argument(
        "--sample-id",
        default=None,
        dest="sample_id",
        help="Sample ID for the report (default: derived from VCF filename)",
    )
    parser.add_argument(
        "--json",
        nargs="?",
        const=True,
        default=None,
        dest="json_flag",
        metavar="PATH",
        help="Also write raw JSON data. Optionally specify a path; defaults to <output>.json",
    )
    parser.add_argument(
        "--skip-api",
        action="store_true",
        dest="skip_api",
        help="Skip external API calls (ClinVar, gnomAD). Uses local SQLite DBs if available. For on-premise deployment.",
    )
    parser.add_argument(
        "--mode",
        choices=["cancer", "rare-disease"],
        default=None,  # will use config default
        help="Report mode (default: cancer). 'cancer' = somatic variant report (FoundationOne style). 'rare-disease' = germline report with HPO/OMIM/ClinGen. See REPORT MODES below.",
    )
    parser.add_argument(
        "--hpo",
        type=str,
        default=None,
        help="Comma-separated HPO IDs for rare disease mode. Enables phenotype-driven candidate gene ranking. (e.g., HP:0001250,HP:0001263)",
    )
    parser.add_argument(
        "--clear-cache",
        action="store_true",
        dest="clear_cache",
        help="Clear variant response cache before running",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose (DEBUG) logging",
    )
    # Batch mode arguments
    parser.add_argument(
        "--batch",
        type=str,
        default=None,
        help="Batch mode: path to directory of VCFs or manifest CSV (columns: sample_id,vcf_path)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel workers for batch report generation (default: CPU count)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output/batch",
        dest="output_dir",
        help="Output directory for batch reports (default: output/batch)",
    )
    parser.add_argument(
        "--hide-vus",
        action="store_true",
        dest="hide_vus",
        help="Hide VUS and Benign variants from detailed report pages (shown in summary count only)",
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)

    if args.clear_cache:
        from scripts.common.cache import clear_cache
        count = clear_cache()
        logger.info(f"Cache cleared: {count} entries removed")

    # Parse HPO IDs from comma-separated string
    hpo_ids = [h.strip() for h in args.hpo.split(",")] if args.hpo else None

    # Resolve mode (shared between single and batch)
    mode = args.mode or get("report.default_mode", "cancer")

    if args.batch:
        # ── Batch mode ────────────────────────────────────────────────────────
        result = run_batch_pipeline(
            batch_path=args.batch,
            output_dir=args.output_dir,
            mode=mode,
            workers=args.workers,
            skip_api=args.skip_api,
            krgdb_path=args.krgdb,
            hpo_ids=hpo_ids,
            hide_vus=args.hide_vus,
        )

        if args.json_flag is not None:
            if args.json_flag is True:
                json_out = str(Path(args.output_dir) / "batch_summary.json")
            else:
                json_out = args.json_flag
            Path(json_out).parent.mkdir(parents=True, exist_ok=True)
            Path(json_out).write_text(
                json.dumps(result, indent=2, default=str, ensure_ascii=False),
                encoding="utf-8",
            )
            _progress(f"  → Batch summary JSON: {json_out}")

        if result["errors"]:
            sys.exit(1)

    else:
        # ── Single-sample mode ────────────────────────────────────────────────
        if not args.vcf_path:
            parser.error("vcf_path is required when not using --batch")

        # Resolve JSON output path
        json_output = None
        if args.json_flag is not None:
            if args.json_flag is True:
                json_output = str(Path(args.output).with_suffix(".json"))
            else:
                json_output = args.json_flag

        result = run_pipeline(
            vcf_path=args.vcf_path,
            output_path=args.output,
            krgdb_path=args.krgdb,
            sample_id=args.sample_id,
            json_output=json_output,
            skip_api=args.skip_api,
            mode=mode,
            hpo_ids=hpo_ids,
            hide_vus=args.hide_vus,
        )

        if result is None:
            sys.exit(1)


if __name__ == "__main__":
    main()
