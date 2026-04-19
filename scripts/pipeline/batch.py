"""Batch processing pipeline — multi-sample VCF analysis with variant deduplication."""

import csv
import logging
import os
import threading
import time as time_mod
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path

from scripts.clinical.hpo_matcher import resolve_hpo_terms
from scripts.common.config import get
from scripts.common.models import FrequencyData
from scripts.db.version_manager import get_all_db_versions
from scripts.intake.parse_vcf import parse_vcf
from scripts.korean_pop.compare_freq import compare_frequencies
from scripts.pipeline.classify import (
    build_summary,
    build_variant_records,
    classify_variants,
    split_variants_for_display,
)
from scripts.pipeline.query import query_variant_databases

logger = logging.getLogger(__name__)


def _progress(msg: str) -> None:
    import sys

    print(msg, file=sys.stderr, flush=True)


def discover_samples(batch_path: str) -> list:
    """Discover VCF files from directory or manifest CSV."""
    path = Path(batch_path)
    if path.is_dir():
        samples = []
        vcf_files = sorted(path.glob("*.vcf")) + sorted(path.glob("*.vcf.gz"))
        for vcf_file in vcf_files:
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


def collect_unique_variants(samples: list) -> tuple:
    """Parse all VCFs, return (unique_variants_dict, sample_variant_map)."""
    unique_variants = {}
    sample_map = {}

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


def _bulk_annotate_variants(unique_variants: dict, krgdb_path: str, skip_api: bool, max_workers: int = 8) -> dict:
    """Annotate all unique variants. Returns {variant_key: annotation_dict}."""
    annotations = {}
    rate_limiter = threading.Semaphore(10)

    def _annotate_one(key, variant):
        with rate_limiter:
            return key, query_variant_databases(variant, krgdb_path, skip_api)

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


def _assemble_sample_report(
    sample, variant_keys, annotations, unique_variants, mode, hpo_ids, hpo_results, krgdb_path, hide_vus=False
):
    """Build report_data dict for one sample using pre-computed annotations."""
    variants = [unique_variants[k] for k in variant_keys]

    db_results = {}
    for key, variant in zip(variant_keys, variants):
        db_results[variant.variant_id] = annotations[key]

    freq_results = {}
    for variant in variants:
        db = db_results[variant.variant_id]
        freq_data = FrequencyData(
            krgdb=db["krgdb_freq"],
            gnomad_eas=db["gnomad"].get("gnomad_eas"),
            gnomad_all=db["gnomad"].get("gnomad_all"),
            korea4k=db.get("korea4k_freq"),
            nard2=db.get("nard2_freq"),
        )
        freq_results[variant.variant_id] = compare_frequencies(freq_data)

    pgx_hits = []
    for variant in variants:
        pgx = db_results[variant.variant_id]["pgx"]
        if pgx:
            pgx_hits.append(pgx)

    classification_results = classify_variants(variants, db_results, freq_results)
    variant_records = build_variant_records(
        variants, db_results, freq_results, classification_results, mode, hpo_results
    )
    summary = build_summary(variant_records)

    tier1, tier2, tier3, tier4_count, detailed_variants, omitted_variants = split_variants_for_display(
        variant_records, hide_vus
    )

    return {
        "sample_id": sample["sample_id"],
        "date": str(date.today()),
        "variants": variant_records,
        "detailed_variants": detailed_variants,
        "omitted_variants": omitted_variants,
        "hide_vus": hide_vus,
        "tier1_variants": tier1,
        "tier2_variants": tier2,
        "tier3_variants": tier3,
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
        "pipeline": {"skip_api": False, "krgdb_path": str(krgdb_path)},
        "mode": mode,
        "hpo_results": hpo_results,
    }


def _generate_single_report(report_data: dict, output_path: str, mode: str) -> str:
    """Generate a single HTML report (runs in subprocess for PDF parallelism)."""
    import sys
    from pathlib import Path as _Path

    sys.path.insert(0, str(_Path(__file__).parent.parent.parent))
    from scripts.counselor.generate_pdf import generate_report_html as _gen_html

    html = _gen_html(report_data, mode=mode)
    _Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    _Path(output_path).write_text(html, encoding="utf-8")
    return output_path


def _generate_reports_parallel(sample_reports: list, output_dir: str, mode: str, workers: int) -> list:
    """Generate HTML reports using ProcessPoolExecutor."""
    _progress(f"[5/5] Generating {len(sample_reports):,} reports ({workers} workers)...")

    results = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {}
        for sample_data in sample_reports:
            output_path = Path(output_dir) / f"{sample_data['sample_id']}_report.html"
            futures[executor.submit(_generate_single_report, sample_data, str(output_path), mode)] = sample_data[
                "sample_id"
            ]

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
    hide_vus: bool = True,
    _skip_reports: bool = False,
) -> dict:
    """Process multiple VCF files with variant deduplication."""
    workers = workers or os.cpu_count() or 4
    start = time_mod.time()

    _progress(f"[1/5] Discovering samples from {batch_path}...")
    samples = discover_samples(batch_path)
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

    unique_variants, sample_map = collect_unique_variants(samples)
    total_variants = sum(len(keys) for keys in sample_map.values())
    cache_hits = total_variants - len(unique_variants)
    _progress(f"  → {total_variants:,} total variants, {len(unique_variants):,} unique ({cache_hits:,} deduplicated)")

    krgdb_file = krgdb_path or get("paths.krgdb", "data/krgdb_freq.tsv")
    annotations = _bulk_annotate_variants(unique_variants, krgdb_file, skip_api, max_workers=min(workers, 10))

    hpo_results = []
    if mode == "rare-disease" and hpo_ids:
        hpo_results = resolve_hpo_terms(hpo_ids)

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
