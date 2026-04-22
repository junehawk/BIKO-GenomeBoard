"""Batch multi-sample pipeline — thin orchestration layer over canonical.

v2.5.4 rewrite. The previous batch implementation parsed every VCF into
one shared ``unique_variants`` dict keyed on ``chrom:pos:ref:alt`` and
re-used the **first sample's Variant object** for every later sample. The
Variant class carries per-sample state (canonical transcript, HGVSp,
trio genotype, selection_reason_list, source, ...) so that reuse leaked
annotation from one patient's report into another — a silent
correctness bug (v2.5.4 H1).

The new implementation calls ``canonical.build_sample_report`` once per
sample. Each sample re-parses its own VCF, so Variant identity does not
leak across samples. External API / local DB traffic is no longer
deduplicated at the Variant level; Phase 4 adds an optional process-
scoped query cache keyed on ``variant_id`` if batch throughput regresses.

Feature parity (v2.5.4 M6): per-sample germline / ped / hpo / intervar /
sv / clinical_note / panel_size_mb columns are supported via an extended
manifest CSV. The legacy two-column manifest (``sample_id,vcf_path``)
stays valid.

Provenance (M3): the pipeline block in each sample's report reflects the
``skip_api`` value the caller actually passed, not the hardcoded False
that pre-v2.5.4 batch wrote regardless of the CLI flag.
"""

from __future__ import annotations

import csv
import logging
import os
import time as time_mod
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from scripts.intake.parse_vcf import parse_vcf
from scripts.orchestration.canonical import build_sample_report, normalize_sample_id

logger = logging.getLogger(__name__)


def _progress(msg: str) -> None:
    import sys

    print(msg, file=sys.stderr, flush=True)


# ── BatchSample dataclass (M6) ───────────────────────────────────────────────


@dataclass
class BatchSample:
    """One row of a batch manifest.

    Only ``sample_id`` and ``vcf_path`` are required; every other field
    is optional and defaults to None / empty. When a field is supplied it
    is forwarded to ``canonical.build_sample_report`` for that one sample.
    """

    sample_id: str
    vcf_path: str
    germline_vcf: str | None = None
    ped_path: str | None = None
    hpo_ids: list[str] = field(default_factory=list)
    sv_path: str | None = None
    intervar_path: str | None = None
    clinical_note: str | None = None
    panel_size_mb: float | None = None


# Manifest CSV columns the parser recognises. Order here only controls
# the order new columns are documented; DictReader is tolerant of any
# header ordering. The legacy two-column manifest
# ("sample_id,vcf_path") continues to work.
_MANIFEST_COLUMNS = (
    "sample_id",
    "vcf_path",
    "germline_vcf",
    "ped_path",
    "hpo_ids",
    "sv_path",
    "intervar_path",
    "clinical_note",
    "panel_size_mb",
)


# ── Discovery (M5 + M6) ──────────────────────────────────────────────────────


def _row_to_batch_sample(row: dict[str, Any]) -> BatchSample:
    """Build a BatchSample from one CSV row."""
    raw_hpo = (row.get("hpo_ids") or "").strip()
    if raw_hpo:
        # Manifest writers may use either ; or | or , (comma would be
        # ambiguous inside an unquoted CSV cell but we allow it when the
        # column was quoted).
        separator = ";" if ";" in raw_hpo else ("|" if "|" in raw_hpo else ",")
        hpo_ids = [x.strip() for x in raw_hpo.split(separator) if x.strip()]
    else:
        hpo_ids = []

    panel_str = (row.get("panel_size_mb") or "").strip()
    panel_size_mb: float | None = None
    if panel_str:
        try:
            panel_size_mb = float(panel_str)
        except ValueError:
            logger.warning("Manifest panel_size_mb not numeric: %r (sample %s)", panel_str, row.get("sample_id"))

    return BatchSample(
        sample_id=row["sample_id"],
        vcf_path=row["vcf_path"],
        germline_vcf=(row.get("germline_vcf") or None) or None,
        ped_path=(row.get("ped_path") or None) or None,
        hpo_ids=hpo_ids,
        sv_path=(row.get("sv_path") or None) or None,
        intervar_path=(row.get("intervar_path") or None) or None,
        clinical_note=(row.get("clinical_note") or None) or None,
        panel_size_mb=panel_size_mb,
    )


def discover_samples(batch_path: str) -> list:
    """Discover VCF files from a directory or a manifest CSV.

    Directory mode: ``*.vcf``, ``*.vcf.gz``, ``*.vcf.bgz`` are picked up.
    sample_id is derived via ``normalize_sample_id`` so single-mode and
    batch-mode agree on the same id for the same file (v2.5.4 M4 fix).

    CSV mode: ``sample_id`` and ``vcf_path`` are required; other columns
    are optional (see ``_MANIFEST_COLUMNS``).

    The returned shape stays a list of dicts for backward compatibility
    with existing consumers (tests, external harnesses); each dict is a
    superset of the legacy {sample_id, vcf_path} schema.
    """
    path = Path(batch_path)
    if path.is_dir():
        samples: list[dict[str, Any]] = []
        vcf_files = sorted(path.glob("*.vcf")) + sorted(path.glob("*.vcf.gz")) + sorted(path.glob("*.vcf.bgz"))
        for vcf_file in vcf_files:
            samples.append(
                {
                    "sample_id": normalize_sample_id(str(vcf_file)),
                    "vcf_path": str(vcf_file),
                }
            )
        return samples
    elif path.suffix == ".csv":
        samples = []
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if not row.get("sample_id") or not row.get("vcf_path"):
                    logger.warning("Manifest row missing sample_id or vcf_path, skipping: %r", row)
                    continue
                bs = _row_to_batch_sample(row)
                # Preserve the legacy dict-shape list but include the
                # extra fields so downstream consumers can feature-detect.
                samples.append(
                    {
                        "sample_id": bs.sample_id,
                        "vcf_path": bs.vcf_path,
                        "germline_vcf": bs.germline_vcf,
                        "ped_path": bs.ped_path,
                        "hpo_ids": bs.hpo_ids,
                        "sv_path": bs.sv_path,
                        "intervar_path": bs.intervar_path,
                        "clinical_note": bs.clinical_note,
                        "panel_size_mb": bs.panel_size_mb,
                    }
                )
        return samples
    else:
        raise ValueError(f"Batch path must be a directory or CSV file: {batch_path}")


# ── Back-compat helper kept for tests (test_batch.py) ────────────────────────


def collect_unique_variants(samples: list) -> tuple:
    """Parse all VCFs and return (unique_variants_dict, sample_variant_map).

    Pre-v2.5.4 the batch pipeline used this to dedup Variant objects
    across samples — which caused the H1 annotation-leak bug. The new
    pipeline does **not** use this function for annotation sharing; it
    is retained only so existing tests that assert on dedup counts can
    still inspect the VCF overlap. New code should call
    ``canonical.build_sample_report`` per sample and not reuse Variant
    instances.
    """
    unique_variants: dict[str, Any] = {}
    sample_map: dict[str, list[str]] = {}

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


# ── Report generation (ProcessPoolExecutor over pre-built report_data) ───────


def _generate_single_report(report_data: dict, output_path: str, mode: str) -> str:
    """Generate a single HTML report (runs in subprocess for PDF parallelism)."""
    import sys as _sys
    from pathlib import Path as _Path

    _sys.path.insert(0, str(_Path(__file__).parent.parent.parent))
    from scripts.reporting.generate_pdf import generate_report_html as _gen_html

    html = _gen_html(report_data, mode=mode)
    _Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    _Path(output_path).write_text(html, encoding="utf-8")
    return output_path


def _generate_reports_parallel(sample_reports: list, output_dir: str, mode: str, workers: int) -> list:
    """Generate HTML reports using ProcessPoolExecutor."""
    _progress(f"[4/4] Generating {len(sample_reports):,} reports ({workers} workers)...")

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
                logger.error("Report generation failed for %s: %s", sample_id, e)

    return results


# ── Per-sample assembly (parallel over canonical) ────────────────────────────


def _assemble_one(
    sample: dict[str, Any],
    *,
    mode: str,
    skip_api: bool,
    hpo_ids: list | None,
    hide_vus: bool,
    clinical_board: bool,
    board_lang: str | None,
) -> dict | None:
    """Call ``canonical.build_sample_report`` for one manifest row.

    Per-sample overrides from the manifest (germline, ped, hpo, sv,
    intervar, clinical_note, panel_size_mb) take precedence over the
    batch-wide defaults; HPO IDs merge (manifest extends the batch list).
    """
    per_sample_hpo = list(sample.get("hpo_ids") or [])
    merged_hpo = list(hpo_ids or []) + per_sample_hpo if per_sample_hpo else hpo_ids

    return build_sample_report(
        vcf_path=sample["vcf_path"],
        sample_id=sample["sample_id"],
        mode=mode,
        skip_api=skip_api,
        hpo_ids=merged_hpo,
        hide_vus=hide_vus,
        sv_path=sample.get("sv_path"),
        panel_size_mb=sample.get("panel_size_mb"),
        intervar_path=sample.get("intervar_path"),
        clinical_board=clinical_board,
        board_lang=board_lang,
        clinical_note=sample.get("clinical_note"),
        germline_vcf=sample.get("germline_vcf"),
        ped_path=sample.get("ped_path"),
    )


def _assemble_all(
    samples: list,
    *,
    mode: str,
    skip_api: bool,
    hpo_ids: list | None,
    hide_vus: bool,
    clinical_board: bool,
    board_lang: str | None,
    workers: int,
) -> tuple[list[dict], list[dict]]:
    """Assemble every sample's report_data. Returns (reports, errors)."""
    sample_reports: list[dict] = []
    errors: list[dict] = []

    # Clinical Board internally drives an LLM serialisation layer; we do
    # not want N ThreadPools × one LLM each hammering a single Ollama
    # server. Force serial assembly when Board synthesis is enabled.
    effective_workers = 1 if clinical_board else max(1, workers)

    if effective_workers == 1:
        for sample in samples:
            try:
                rd = _assemble_one(
                    sample,
                    mode=mode,
                    skip_api=skip_api,
                    hpo_ids=hpo_ids,
                    hide_vus=hide_vus,
                    clinical_board=clinical_board,
                    board_lang=board_lang,
                )
                if rd is not None:
                    sample_reports.append(rd)
                else:
                    errors.append({"sample": sample["sample_id"], "error": "no variants parsed"})
            except Exception as e:
                logger.error("Failed to assemble report for %s: %s", sample["sample_id"], e)
                errors.append({"sample": sample["sample_id"], "error": str(e)})
        return sample_reports, errors

    with ThreadPoolExecutor(max_workers=effective_workers) as executor:
        futures = {
            executor.submit(
                _assemble_one,
                s,
                mode=mode,
                skip_api=skip_api,
                hpo_ids=hpo_ids,
                hide_vus=hide_vus,
                clinical_board=clinical_board,
                board_lang=board_lang,
            ): s
            for s in samples
        }
        for future in as_completed(futures):
            sample = futures[future]
            try:
                rd = future.result()
                if rd is not None:
                    sample_reports.append(rd)
                else:
                    errors.append({"sample": sample["sample_id"], "error": "no variants parsed"})
            except Exception as e:
                logger.error("Failed to assemble report for %s: %s", sample["sample_id"], e)
                errors.append({"sample": sample["sample_id"], "error": str(e)})
    return sample_reports, errors


# ── Top-level entry ──────────────────────────────────────────────────────────


def run_batch_pipeline(
    batch_path: str,
    output_dir: str = "output/batch",
    mode: str = "cancer",
    workers: int | None = None,
    skip_api: bool = False,
    hpo_ids: list | None = None,
    hide_vus: bool = True,
    clinical_board: bool = False,
    board_lang: str | None = None,
    _skip_reports: bool = False,
) -> dict:
    """Process multiple VCF files and emit one report per sample.

    Each sample goes through ``canonical.build_sample_report`` — the same
    function ``run_pipeline`` calls for single-sample mode. There is no
    variant-level deduplication across samples (that was the v2.5.4 H1
    correctness bug).
    """
    workers = workers or os.cpu_count() or 4
    start = time_mod.time()

    _progress(f"[1/4] Discovering samples from {batch_path}...")
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

    _progress(f"[2/4] Assembling {len(samples)} sample reports (workers={workers})...")
    sample_reports, errors = _assemble_all(
        samples,
        mode=mode,
        skip_api=skip_api,
        hpo_ids=hpo_ids,
        hide_vus=hide_vus,
        clinical_board=clinical_board,
        board_lang=board_lang,
        workers=workers,
    )

    _progress(f"[3/4] Assembled {len(sample_reports)} reports ({len(errors)} errors)")

    # Count variants from the actual report outputs; this number is
    # post-filter (the variants the reports actually render) which is
    # what downstream callers want.
    total_variants = sum(len(rd.get("variants", [])) for rd in sample_reports)

    report_paths: list[str] = []
    if not _skip_reports and output_dir is not None:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        report_paths = _generate_reports_parallel(sample_reports, output_dir, mode, workers)

    elapsed = time_mod.time() - start
    _progress(f"\nBatch complete! {len(report_paths)} reports in {elapsed:.1f}s")
    if errors:
        _progress(f"  Errors: {len(errors)}")

    # unique_variants / cache_hits are legacy keys maintained for result
    # schema backward compatibility. The values are now always 0 /
    # total_variants because the v2.5.4 pipeline runs every sample in
    # isolation — Phase 4 may reintroduce a correctness-safe variant_id
    # cache that would populate non-zero hit counts.
    return {
        "samples_processed": len(samples),
        "total_variants": total_variants,
        "unique_variants": total_variants,
        "cache_hits": 0,
        "reports_generated": report_paths,
        "errors": errors,
        "elapsed_seconds": elapsed,
    }
