#!/usr/bin/env python3
"""BIKO GenomeBoard Standalone Pipeline — VCF → Report

Entry point for single-sample and batch variant analysis.
Core logic is in scripts/orchestration/ modules:
  - pipeline.query: per-variant database queries
  - pipeline.classify: ACMG classification, variant record assembly
  - pipeline.batch: batch processing with deduplication
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.common.config import get
from scripts.reporting.generate_pdf import generate_pdf, generate_report_html
from scripts.orchestration.batch import (
    collect_unique_variants,
    discover_samples,
    run_batch_pipeline,
)
from scripts.orchestration.canonical import (
    build_sample_report,
    normalize_sample_id,
    report_data_for_json,
)
from scripts.orchestration.classify import (
    build_summary,
    build_variant_records,
    classify_variants,
    split_variants_for_display,  # re-exported for tests that monkey-patch this path
    sv_to_dict,
)

# Pipeline modules (extracted from this file)
from scripts.orchestration.query import query_variant_databases

__all__ = [
    "run_pipeline",
    "build_summary",
    "build_variant_records",
    "classify_variants",
    "split_variants_for_display",
    "sv_to_dict",
    "query_variant_databases",
    "normalize_sample_id",
    "collect_unique_variants",
    "discover_samples",
    "run_batch_pipeline",
]

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("genomeboard")

# Backward compatibility aliases — tests import _prefixed names from this module
_query_variant_databases = query_variant_databases
_classify_variants = classify_variants
_build_variant_records = build_variant_records
_build_summary = build_summary
_sv_to_dict = sv_to_dict
_discover_samples = discover_samples
_collect_unique_variants = collect_unique_variants


def _progress(msg: str) -> None:
    """Print a progress message to stderr."""
    print(msg, file=sys.stderr, flush=True)


def _format_board_summary(board_opinion) -> str:
    """Format a one-line Clinical Board summary, branching on opinion type."""
    from scripts.clinical_board.models import CancerBoardOpinion

    if isinstance(board_opinion, CancerBoardOpinion):
        return (
            f"  → Therapeutic implications: {board_opinion.therapeutic_implications} "
            f"({board_opinion.confidence} confidence)"
        )
    return f"  → Primary diagnosis: {board_opinion.primary_diagnosis} ({board_opinion.confidence} confidence)"


def run_pipeline(
    vcf_path: str,
    output_path: str = "output/report.html",
    sample_id: str = None,
    json_output: str = None,
    skip_api: bool = False,
    mode: str = None,
    hpo_ids: list = None,
    hide_vus: bool = True,
    sv_path: str = None,
    panel_size_mb: float = None,
    bed_path: str = None,
    intervar_path: str = None,
    clinical_board: bool = False,
    board_lang: str = None,
    clinical_note: str = None,
    germline_vcf: str = None,
    ped_path: str = None,
) -> dict:
    """Run the full BIKO GenomeBoard analysis pipeline.

    Thin CLI-facing wrapper over ``canonical.build_sample_report``. Handles
    progress logging, path setup, HTML / PDF emission, and JSON dump.
    The actual assembly logic lives in
    ``scripts.orchestration.canonical.build_sample_report``.

    Returns the assembled report_data dict (or ``None`` on fatal error).
    """
    mode = mode or get("report.default_mode", "cancer")
    start_time = time.time()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if sample_id is None:
        # v2.5.4 — previously used .stem.upper(); the canonical helper now
        # strips .vcf / .vcf.gz / .vcf.bgz and preserves case. Batch mode
        # uses the same helper so sample_id agrees between modes.
        sample_id = normalize_sample_id(vcf_path)

    _progress(f"[1/6] Parsing VCF: {vcf_path}")
    if ped_path:
        _progress(f"  → Trio resolution: PED file ({ped_path})")
    if skip_api:
        _progress("  [skip-api mode: external API calls disabled, using local DBs only]")

    # Warn when primary VCF is somatic/de novo and no germline provided.
    # This used to live inside the PGx step — now surfaced here before we
    # delegate so the warning order matches pre-refactor logs.
    if not germline_vcf and mode in ("cancer", "rare-disease"):
        logger.warning(
            "PGx results from somatic/de novo input may be incomplete. "
            "Provide --germline for accurate pharmacogenomics assessment."
        )

    report_data = build_sample_report(
        vcf_path=vcf_path,
        sample_id=sample_id,
        mode=mode,
        skip_api=skip_api,
        hpo_ids=hpo_ids,
        hide_vus=hide_vus,
        sv_path=sv_path,
        panel_size_mb=panel_size_mb,
        bed_path=bed_path,
        intervar_path=intervar_path,
        clinical_board=clinical_board,
        board_lang=board_lang,
        clinical_note=clinical_note,
        germline_vcf=germline_vcf,
        ped_path=ped_path,
    )

    if report_data is None:
        # build_sample_report already logged the reason (missing VCF or
        # zero parsed variants). Keep the CLI contract: return None so
        # main() can exit with a non-zero status.
        return None

    summary = report_data["summary"]
    _progress(f"  → Found {len(report_data['variants'])} variants")
    _progress(f"  [PGx engine: {report_data.get('pgx_source', '')}]")

    if clinical_board and report_data.get("clinical_board"):
        _progress(_format_board_summary(report_data["clinical_board"]))

    # ── Generate report ───────────────────────────────────────────────────────
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
        dump_data = report_data_for_json(report_data)
        json_path.write_text(
            json.dumps(dump_data, indent=2, default=str, ensure_ascii=False),
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


# ── CLI ──────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(
        description="BIKO GenomeBoard — Korean Population-Aware Genomic Variant Interpretation Pipeline",
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

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
EXAMPLES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  python scripts/orchestrate.py sample.vcf -o report.html --json
  python scripts/orchestrate.py sample.vcf -o report.html --skip-api
  python scripts/orchestrate.py patient.vcf --mode rare-disease --hpo HP:0001250 -o report.html
  python scripts/orchestrate.py --batch vcf_dir/ --output-dir output/batch --workers 8
        """,
    )
    parser.add_argument("vcf_path", nargs="?", default=None, help="Path to input VCF file (omit when using --batch)")
    parser.add_argument("--output", "-o", default="output/report.html", help="Output file path (.html or .pdf)")
    parser.add_argument("--sample-id", default=None, dest="sample_id", help="Sample ID for the report")
    parser.add_argument(
        "--json", nargs="?", const=True, default=None, dest="json_flag", metavar="PATH", help="Also write raw JSON data"
    )
    parser.add_argument(
        "--skip-api", action="store_true", dest="skip_api", help="Skip external API calls, use local DBs only"
    )
    parser.add_argument(
        "--mode", choices=["cancer", "rare-disease"], default=None, help="Report mode (default: cancer)"
    )
    parser.add_argument("--hpo", type=str, default=None, help="Comma-separated HPO IDs for rare disease mode")
    parser.add_argument("--clear-cache", action="store_true", dest="clear_cache", help="Clear variant response cache")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose (DEBUG) logging")
    parser.add_argument("--batch", type=str, default=None, help="Batch mode: directory of VCFs or manifest CSV")
    parser.add_argument("--workers", type=int, default=None, help="Parallel workers for batch mode")
    parser.add_argument(
        "--output-dir", type=str, default="output/batch", dest="output_dir", help="Batch output directory"
    )
    parser.add_argument(
        "--hide-vus", action="store_true", dest="hide_vus", help="(Default) Hide VUS/Benign from detail pages"
    )
    parser.add_argument(
        "--show-all-variants", action="store_true", dest="show_all_variants", help="Show ALL variants in detail pages"
    )
    parser.add_argument("--sv", dest="sv_path", help="AnnotSV TSV file for CNV/SV integration")
    parser.add_argument("--panel-size", type=float, dest="panel_size", help="Panel size in Mb for TMB calculation")
    parser.add_argument("--bed", dest="bed_path", help="BED file for panel size calculation (overrides --panel-size)")
    parser.add_argument("--intervar", dest="intervar_path", help="InterVar output TSV for ACMG evidence codes")
    parser.add_argument(
        "--clinical-board",
        action="store_true",
        dest="clinical_board",
        help="Enable Clinical Board diagnostic synthesis (requires Ollama)",
    )
    parser.add_argument(
        "--board-lang",
        default=None,
        dest="board_lang",
        choices=["en", "ko"],
        help="Clinical Board output language (default: en)",
    )
    parser.add_argument(
        "--germline",
        type=str,
        default=None,
        help="Germline VCF for pharmacogenomics analysis. PharmCAT will be used "
        "if available; falls back to built-in SNV matching otherwise. "
        "Accepts .vcf, .vcf.gz, .vcf.bgz.",
    )
    parser.add_argument(
        "--ped",
        type=str,
        default=None,
        dest="ped_path",
        help="Plink-style PED file for trio family structure. When provided, "
        "trio proband and parents are resolved from PED in STRICT MODE "
        "(fails if no trio can be identified against VCF samples). "
        "Without --ped, the parser falls back to the filename-based "
        "proband heuristic (3-sample VCFs only). PED mode supports "
        "quartet / N≥3 VCFs.",
    )
    note_group = parser.add_mutually_exclusive_group()
    note_group.add_argument(
        "--clinical-note",
        type=str,
        default=None,
        dest="clinical_note",
        help="Free-text clinical note (Korean or English) for the Clinical Board",
    )
    note_group.add_argument(
        "--clinical-note-file",
        type=str,
        default=None,
        dest="clinical_note_file",
        help="Path to file containing the clinical note",
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)

    if args.clear_cache:
        from scripts.common.cache import clear_cache

        count = clear_cache()
        logger.info(f"Cache cleared: {count} entries removed")

    hpo_ids = [h.strip() for h in args.hpo.split(",")] if args.hpo else None
    mode = args.mode or get("report.default_mode", "cancer")
    effective_hide_vus = not args.show_all_variants

    if args.batch:
        result = run_batch_pipeline(
            batch_path=args.batch,
            output_dir=args.output_dir,
            mode=mode,
            workers=args.workers,
            skip_api=args.skip_api,
            hpo_ids=hpo_ids,
            hide_vus=effective_hide_vus,
        )

        if args.json_flag is not None:
            json_out = (
                args.json_flag if args.json_flag is not True else str(Path(args.output_dir) / "batch_summary.json")
            )
            Path(json_out).parent.mkdir(parents=True, exist_ok=True)
            Path(json_out).write_text(json.dumps(result, indent=2, default=str, ensure_ascii=False), encoding="utf-8")
            _progress(f"  → Batch summary JSON: {json_out}")

        if result["errors"]:
            sys.exit(1)
    else:
        if not args.vcf_path:
            parser.error("vcf_path is required when not using --batch")

        json_output = None
        if args.json_flag is not None:
            json_output = args.json_flag if args.json_flag is not True else str(Path(args.output).with_suffix(".json"))

        clinical_note_text = getattr(args, "clinical_note", None)
        clinical_note_file = getattr(args, "clinical_note_file", None)
        if clinical_note_file:
            try:
                clinical_note_text = Path(clinical_note_file).read_text(encoding="utf-8")
            except OSError as e:
                logger.error(f"Failed to read clinical note file: {e}")
                sys.exit(1)

        result = run_pipeline(
            vcf_path=args.vcf_path,
            output_path=args.output,
            sample_id=args.sample_id,
            json_output=json_output,
            skip_api=args.skip_api,
            mode=mode,
            hpo_ids=hpo_ids,
            hide_vus=effective_hide_vus,
            sv_path=args.sv_path,
            panel_size_mb=args.panel_size,
            bed_path=args.bed_path,
            intervar_path=getattr(args, "intervar_path", None),
            clinical_board=getattr(args, "clinical_board", False),
            board_lang=getattr(args, "board_lang", None),
            clinical_note=clinical_note_text,
            germline_vcf=getattr(args, "germline", None),
            ped_path=getattr(args, "ped_path", None),
        )

        if result is None:
            sys.exit(1)


if __name__ == "__main__":
    main()
