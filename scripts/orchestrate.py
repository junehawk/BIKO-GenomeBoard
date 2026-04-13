#!/usr/bin/env python3
"""GenomeBoard Standalone Pipeline — VCF → Report

Entry point for single-sample and batch variant analysis.
Core logic is in scripts/pipeline/ modules:
  - pipeline.query: per-variant database queries
  - pipeline.classify: ACMG classification, variant record assembly
  - pipeline.batch: batch processing with deduplication
"""

import argparse
import json
import logging
import sys
import time
from datetime import date
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.intake.parse_vcf import parse_vcf
from scripts.korean_pop.compare_freq import compare_frequencies
from scripts.clinical.hpo_matcher import resolve_hpo_terms
from scripts.db.version_manager import get_all_db_versions
from scripts.counselor.generate_pdf import generate_report_html, generate_pdf
from scripts.common.models import FrequencyData
from scripts.common.config import get

# Pipeline modules (extracted from this file)
from scripts.pipeline.query import query_variant_databases
from scripts.pipeline.classify import (
    classify_variants, build_variant_records, build_summary,
    sv_to_dict, split_variants_for_display,
)
from scripts.pipeline.batch import (
    discover_samples, collect_unique_variants, run_batch_pipeline,
)

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


def run_pipeline(
    vcf_path: str,
    output_path: str = "output/report.html",
    krgdb_path: str = "data/krgdb_freq.tsv",
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
) -> dict:
    """Run the full GenomeBoard analysis pipeline.

    Returns the assembled report data dict (or None on fatal error).
    """
    mode = mode or get("report.default_mode", "cancer")
    start_time = time.time()

    # HPO processing (rare disease mode)
    hpo_results = []
    if mode == "rare-disease" and hpo_ids:
        logger.info("[HPO] Resolving phenotype terms...")
        hpo_results = resolve_hpo_terms(hpo_ids)
        for h in hpo_results:
            logger.info(f"  → {h['id']}: {h['name']} ({len(h['genes'])} associated genes)")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

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

    krgdb_file = Path(krgdb_path)
    if not krgdb_file.exists():
        logger.warning(f"KRGDB file not found: {krgdb_path} — continuing without Korean frequency data")

    # ── Step 2: Query databases per variant ──────────────────────────────────
    _progress("[2/6] Querying databases (ClinVar, gnomAD, KRGDB)...")
    if skip_api:
        _progress("  [skip-api mode: ClinVar and gnomAD queries disabled]")

    db_results = {}
    for variant in variants:
        result = query_variant_databases(variant, str(krgdb_file), skip_api)
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
            korea4k=db.get("korea4k_freq"),
            nard2=db.get("nard2_freq"),
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

    intervar_data = None
    if intervar_path:
        try:
            from scripts.intake.parse_intervar import parse_intervar
            intervar_data = parse_intervar(intervar_path)
            _progress(f"  [InterVar] Loaded evidence for {len(intervar_data)} variants")
        except Exception as e:
            logger.warning(f"InterVar parsing failed: {e}")

    classification_results = classify_variants(variants, db_results, freq_results, intervar_data=intervar_data)
    for variant in variants:
        classification = classification_results[variant.variant_id]
        codes_str = "+".join(classification.evidence_codes) if classification.evidence_codes else "no codes"
        gene_label = variant.gene or variant.variant_id
        _progress(f"  → {gene_label}: {classification.classification} ({codes_str})")

    # ── Step 6: Assemble report data ──────────────────────────────────────────
    variant_records = build_variant_records(
        variants, db_results, freq_results, classification_results, mode, hpo_results
    )
    summary = build_summary(variant_records)

    tier1, tier2, tier3, tier4_count, detailed_variants, omitted_variants = \
        split_variants_for_display(variant_records, hide_vus)

    report_data = {
        "sample_id": sample_id,
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
        "db_versions": get_all_db_versions(skip_api=skip_api),
        "pipeline": {
            "skip_api": skip_api,
            "krgdb_path": str(krgdb_file),
        },
        "mode": mode,
        "hpo_results": hpo_results,
    }

    # Parse structural variants if provided
    sv_variants = []
    if sv_path:
        from scripts.intake.parse_annotsv import parse_annotsv
        sv_variants = parse_annotsv(sv_path)
        logger.info(f"[SV] Parsed {len(sv_variants)} structural variants")

    sv_class45 = [sv for sv in sv_variants if sv.acmg_class in (4, 5)]
    sv_class3_all = [sv for sv in sv_variants if sv.acmg_class == 3]
    sv_class3_display = [sv for sv in sv_class3_all if sv.is_dosage_sensitive(mode)]
    sv_class3_hidden = len(sv_class3_all) - len(sv_class3_display)
    sv_benign_count = sum(1 for sv in sv_variants if sv.acmg_class in (1, 2))

    report_data["sv_variants"] = [sv_to_dict(sv) for sv in sv_variants]
    report_data["sv_class45"] = [sv_to_dict(sv) for sv in sv_class45]
    report_data["sv_class3_display"] = [sv_to_dict(sv) for sv in sv_class3_display]
    report_data["sv_class3_hidden"] = sv_class3_hidden
    report_data["sv_benign_count"] = sv_benign_count

    # Calculate TMB (cancer mode only)
    if mode == "cancer":
        from scripts.somatic.tmb import calculate_tmb, calculate_panel_size_from_bed
        tmb_panel = panel_size_mb
        if bed_path:
            tmb_panel = calculate_panel_size_from_bed(bed_path)
        tmb_result = calculate_tmb(variants, panel_size_mb=tmb_panel)
        report_data["tmb"] = {
            "score": tmb_result.score,
            "level": tmb_result.level,
            "variant_count": tmb_result.variant_count,
            "total_variants": tmb_result.total_variants,
            "panel_size_mb": tmb_result.panel_size_mb,
            "counted_consequences": tmb_result.counted_consequences,
        }
        logger.info(f"[TMB] {tmb_result.score} mut/Mb ({tmb_result.level}) — "
                     f"{tmb_result.variant_count}/{tmb_result.total_variants} variants, "
                     f"{tmb_result.panel_size_mb} Mb panel")
    else:
        report_data["tmb"] = None

    # ── Clinical Board (optional LLM diagnostic synthesis) ────────────────────
    if clinical_board:
        try:
            from scripts.clinical_board.runner import run_clinical_board
            from scripts.clinical_board.render import render_board_opinion_html
            _progress("[Board] Running Clinical Board diagnostic synthesis...")
            board_opinion = run_clinical_board(report_data, mode, language=board_lang)
            if board_opinion:
                report_data["clinical_board"] = board_opinion
                report_data["clinical_board_html"] = render_board_opinion_html(board_opinion, language=board_lang or get("clinical_board.language", "en"))
                _progress(f"  → Primary diagnosis: {board_opinion.primary_diagnosis} "
                         f"({board_opinion.confidence} confidence)")
            else:
                _progress("  → Clinical Board skipped (Ollama not available)")
        except Exception as e:
            logger.warning(f"Clinical Board failed: {e}")
            _progress(f"  → Clinical Board error: {e}")

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


# ── CLI ──────────────────────────────────────────────────────────────────────

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
    parser.add_argument("--krgdb", default="data/krgdb_freq.tsv", help="KRGDB frequency data file")
    parser.add_argument("--sample-id", default=None, dest="sample_id", help="Sample ID for the report")
    parser.add_argument("--json", nargs="?", const=True, default=None, dest="json_flag", metavar="PATH", help="Also write raw JSON data")
    parser.add_argument("--skip-api", action="store_true", dest="skip_api", help="Skip external API calls, use local DBs only")
    parser.add_argument("--mode", choices=["cancer", "rare-disease"], default=None, help="Report mode (default: cancer)")
    parser.add_argument("--hpo", type=str, default=None, help="Comma-separated HPO IDs for rare disease mode")
    parser.add_argument("--clear-cache", action="store_true", dest="clear_cache", help="Clear variant response cache")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose (DEBUG) logging")
    parser.add_argument("--batch", type=str, default=None, help="Batch mode: directory of VCFs or manifest CSV")
    parser.add_argument("--workers", type=int, default=None, help="Parallel workers for batch mode")
    parser.add_argument("--output-dir", type=str, default="output/batch", dest="output_dir", help="Batch output directory")
    parser.add_argument("--hide-vus", action="store_true", dest="hide_vus", help="(Default) Hide VUS/Benign from detail pages")
    parser.add_argument("--show-all-variants", action="store_true", dest="show_all_variants", help="Show ALL variants in detail pages")
    parser.add_argument("--sv", dest="sv_path", help="AnnotSV TSV file for CNV/SV integration")
    parser.add_argument("--panel-size", type=float, dest="panel_size", help="Panel size in Mb for TMB calculation")
    parser.add_argument("--bed", dest="bed_path", help="BED file for panel size calculation (overrides --panel-size)")
    parser.add_argument("--intervar", dest="intervar_path", help="InterVar output TSV for ACMG evidence codes")
    parser.add_argument("--clinical-board", action="store_true", dest="clinical_board", help="Enable Clinical Board diagnostic synthesis (requires Ollama)")
    parser.add_argument("--board-lang", default=None, dest="board_lang", choices=["en", "ko"], help="Clinical Board output language (default: en)")

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
            batch_path=args.batch, output_dir=args.output_dir, mode=mode,
            workers=args.workers, skip_api=args.skip_api, krgdb_path=args.krgdb,
            hpo_ids=hpo_ids, hide_vus=effective_hide_vus,
        )

        if args.json_flag is not None:
            json_out = args.json_flag if args.json_flag is not True else str(Path(args.output_dir) / "batch_summary.json")
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

        result = run_pipeline(
            vcf_path=args.vcf_path, output_path=args.output, krgdb_path=args.krgdb,
            sample_id=args.sample_id, json_output=json_output, skip_api=args.skip_api,
            mode=mode, hpo_ids=hpo_ids, hide_vus=effective_hide_vus,
            sv_path=args.sv_path, panel_size_mb=args.panel_size, bed_path=args.bed_path,
            intervar_path=getattr(args, "intervar_path", None),
            clinical_board=getattr(args, "clinical_board", False),
            board_lang=getattr(args, "board_lang", None),
        )

        if result is None:
            sys.exit(1)


if __name__ == "__main__":
    main()
