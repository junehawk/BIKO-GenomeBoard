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
from dataclasses import asdict, is_dataclass
from datetime import date
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.enrichment.hpo_matcher import resolve_hpo_terms
from scripts.common.config import get
from scripts.common.models import FrequencyData
from scripts.reporting.generate_pdf import generate_pdf, generate_report_html
from scripts.storage.version_manager import get_all_db_versions
from scripts.intake.parse_vcf import parse_vcf
from scripts.population.compare_freq import compare_frequencies
from scripts.orchestration.batch import (
    collect_unique_variants,
    discover_samples,
    run_batch_pipeline,
)
from scripts.orchestration.classify import (
    build_summary,
    build_variant_records,
    classify_variants,
    split_variants_for_display,
    sv_to_dict,
)

# Pipeline modules (extracted from this file)
from scripts.orchestration.query import query_variant_databases

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
    clinical_note: str = None,
    germline_vcf: str = None,
    ped_path: str = None,
) -> dict:
    """Run the full BIKO GenomeBoard analysis pipeline.

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

    variants = parse_vcf(str(vcf_file), ped_path=ped_path)
    if not variants:
        logger.error(f"No variants parsed from {vcf_path}")
        return None
    _progress(f"  → Found {len(variants)} variants")
    if ped_path:
        _progress(f"  → Trio resolution: PED file ({ped_path})")

    krgdb_file = Path(krgdb_path)
    if not krgdb_file.exists():
        logger.warning(f"KRGDB file not found: {krgdb_path} — continuing without Korean frequency data")

    # ── Step 2: Query databases per variant ──────────────────────────────────
    _progress("[2/6] Querying databases (ClinVar, gnomAD, KRGDB)...")
    if skip_api:
        _progress("  [skip-api mode: external API calls disabled, using local DBs only]")

    db_results = {}
    for variant in variants:
        result = query_variant_databases(variant, str(krgdb_file), skip_api)
        db_results[variant.variant_id] = result

        clinvar_sig = result["clinvar"].get("clinvar_significance", "Not Found")
        gnomad_all = result["gnomad"].get("gnomad_all")
        krgdb_val = result["krgdb_freq"]

        gnomad_str = f"{gnomad_all:.5f}" if gnomad_all is not None else "N/A"
        krgdb_str = f"{krgdb_val:.5f}" if krgdb_val is not None else "N/A"
        gene_label = variant.gene or variant.variant_id
        _progress(f"  → {gene_label}: ClinVar={clinvar_sig}, gnomAD={gnomad_str}, KRGDB={krgdb_str}")

    # ── Step 2b: Extract inherited variants from germline VCF (rare-disease) ─
    inherited_variants: list = []
    if germline_vcf and mode == "rare-disease":
        try:
            from scripts.orchestration.extract_germline import extract_inherited_variants

            primary_ids = {v.variant_id for v in variants}
            inherited_variants = extract_inherited_variants(
                germline_vcf,
                primary_variant_ids=primary_ids,
            )
            if inherited_variants:
                _progress(f"  [Germline] Extracted {len(inherited_variants)} inherited target variants")
                for iv in inherited_variants:
                    result = query_variant_databases(iv, str(krgdb_file), skip_api)
                    db_results[iv.variant_id] = result
                    gene_label = iv.gene or iv.variant_id
                    clinvar_sig = result["clinvar"].get("clinvar_significance", "Not Found")
                    _progress(f"  → {gene_label} (inherited): ClinVar={clinvar_sig}")
                # Merge into main variant list (primary takes precedence via dedup above)
                variants = list(variants) + inherited_variants
            else:
                _progress("  [Germline] No inherited target variants found")
        except Exception as e:
            logger.warning("Germline inherited variant extraction failed: %s", e)

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

    from scripts.pharmacogenomics.korean_pgx import get_pgx_results

    pgx_data = get_pgx_results(variants, germline_vcf=germline_vcf)
    pgx_source = pgx_data["pgx_source"]
    _progress(f"  [PGx engine: {pgx_source}]")

    if pgx_data["warnings"]:
        for w in pgx_data["warnings"]:
            _progress(f"  [PGx warning] {w}")

    # Warn when primary VCF is somatic/de novo and no germline provided
    if not germline_vcf and mode in ("cancer", "rare-disease"):
        logger.warning(
            "PGx results from somatic/de novo input may be incomplete. "
            "Provide --germline for accurate pharmacogenomics assessment."
        )

    # Backward-compat: collect PgxResult objects from per-variant db_results
    # for the legacy pgx_hits list (used when builtin engine runs via query)
    pgx_hits = []
    if pgx_source in ("builtin", "builtin_limited"):
        for variant in variants:
            pgx = db_results[variant.variant_id]["pgx"]
            if pgx:
                pgx_hits.append(pgx)
                cpic_str = getattr(pgx, "cpic_level", "")
                _progress(f"  -> {pgx.gene}: {pgx.phenotype}, CPIC Level {cpic_str}")
    else:
        # PharmCAT path: log from the unified results
        for hit in pgx_data["pgx_hits"]:
            _progress(f"  -> {hit['gene']}: {hit.get('phenotype', '')}")

    # ── Step 5: ACMG classification ───────────────────────────────────────────
    _progress("[5/6] Running ACMG classification engine...")

    intervar_data = None
    if intervar_path:
        try:
            from scripts.annotation.parse_intervar import parse_intervar

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

    tier1, tier2, tier3, tier4_count, detailed_variants, omitted_variants = split_variants_for_display(
        variant_records, hide_vus
    )

    # Build pgx_results list — prefer PharmCAT hits when available,
    # otherwise serialise from per-variant PgxResult objects.
    if pgx_source == "pharmcat":
        pgx_results_list = pgx_data["pgx_hits"]
    else:
        pgx_results_list = [
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
        ]

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
        "pgx_results": pgx_results_list,
        "pgx_source": pgx_source,
        "germline_provided": pgx_data["germline_provided"],
        "pharmcat_version": pgx_data["pharmcat_version"],
        "summary": summary,
        "db_versions": get_all_db_versions(skip_api=skip_api),
        "pipeline": {
            "skip_api": skip_api,
            "krgdb_path": str(krgdb_file),
            "ped_used": bool(ped_path),
            "ped_path": str(ped_path) if ped_path else None,
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
        from scripts.somatic.tmb import calculate_panel_size_from_bed, calculate_tmb

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
        logger.info(
            f"[TMB] {tmb_result.score} mut/Mb ({tmb_result.level}) — "
            f"{tmb_result.variant_count}/{tmb_result.total_variants} variants, "
            f"{tmb_result.panel_size_mb} Mb panel"
        )
    else:
        report_data["tmb"] = None

    # ── Clinical Board (optional LLM diagnostic synthesis) ────────────────────
    if clinical_note:
        report_data["clinical_note"] = clinical_note
    if clinical_board:
        try:
            from scripts.clinical_board.render import render_board_opinion_html
            from scripts.clinical_board.runner import run_clinical_board

            _progress("[Board] Running Clinical Board diagnostic synthesis...")
            board_opinion = run_clinical_board(report_data, mode, language=board_lang)
            if board_opinion:
                report_data["clinical_board"] = board_opinion
                report_data["clinical_board_html"] = render_board_opinion_html(
                    board_opinion, language=board_lang or get("clinical_board.language", "en")
                )
                _progress(_format_board_summary(board_opinion))
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
        # Unwrap the BoardOpinion/CancerBoardOpinion dataclass so that
        # rerender_report.py can round-trip it back to a dataclass instance.
        # json.dumps(default=str) would otherwise stringify it as a repr.
        dump_data = report_data
        cb = report_data.get("clinical_board")
        if cb is not None and is_dataclass(cb):
            dump_data = {**report_data, "clinical_board": asdict(cb)}
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
    parser.add_argument("--krgdb", default="data/krgdb_freq.tsv", help="KRGDB frequency data file")
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
            krgdb_path=args.krgdb,
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
            krgdb_path=args.krgdb,
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
