#!/usr/bin/env python3
"""GenomeBoard Standalone Pipeline — VCF → Report"""

import argparse
import json
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.intake.parse_vcf import parse_vcf
from scripts.clinical.query_clinvar import query_clinvar
from scripts.db.query_local_clinvar import query_local_clinvar, get_db_version as get_clinvar_db_version
from scripts.clinical.hpo_matcher import resolve_hpo_terms, calculate_hpo_score, get_matching_hpo_terms
from scripts.clinical.query_omim import query_omim
from scripts.clinical.query_clingen import get_gene_validity
from scripts.korean_pop.query_gnomad import query_gnomad
from scripts.korean_pop.query_krgdb import query_krgdb
from scripts.korean_pop.compare_freq import compare_frequencies
from scripts.pharma.korean_pgx import check_korean_pgx
from scripts.classification.acmg_engine import classify_variant, check_clinvar_conflict
from scripts.counselor.generate_pdf import generate_report_html, generate_pdf
from scripts.common.models import AcmgEvidence, FrequencyData
from scripts.common.config import get

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
        # Also attempt local ClinVar DB regardless of annotation.source
        try:
            clinvar_result = query_local_clinvar(variant)
        except Exception as e:
            logger.warning(f"Local ClinVar lookup failed for {variant.variant_id}: {e}")
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
                return query_gnomad(variant)
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


def run_pipeline(
    vcf_path: str,
    output_path: str = "output/report.html",
    krgdb_path: str = "data/krgdb_freq.tsv",
    sample_id: str = None,
    json_output: str = None,
    skip_api: bool = False,
    mode: str = None,
    hpo_ids: list = None,
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
        conflict = check_clinvar_conflict(classification.classification, clinvar_sig)
        classification.conflict = conflict

        classification_results[variant.variant_id] = classification
        codes_str = "+".join(classification.evidence_codes) if classification.evidence_codes else "no codes"
        gene_label = variant.gene or variant.variant_id
        _progress(f"  → {gene_label}: {classification.classification} ({codes_str})")

    # ── Step 6: Assemble report data ──────────────────────────────────────────
    variant_records = []
    for variant in variants:
        db = db_results[variant.variant_id]
        freq = freq_results[variant.variant_id]
        classification = classification_results[variant.variant_id]
        clinvar = db["clinvar"]

        freq_data_obj = FrequencyData(
            krgdb=db["krgdb_freq"],
            gnomad_eas=db["gnomad"].get("gnomad_eas"),
            gnomad_all=db["gnomad"].get("gnomad_all"),
        )

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

    # Build summary counts
    total = len(variant_records)
    pathogenic_count = sum(1 for r in variant_records if r["classification"] == "Pathogenic")
    likely_path_count = sum(1 for r in variant_records if r["classification"] == "Likely Pathogenic")
    drug_response_count = sum(1 for r in variant_records if r["classification"] == "Drug Response")
    risk_factor_count = sum(1 for r in variant_records if r["classification"] == "Risk Factor")
    vus_count = sum(1 for r in variant_records if r["classification"] == "VUS")
    likely_benign_count = sum(1 for r in variant_records if r["classification"] == "Likely Benign")
    benign_count = sum(1 for r in variant_records if r["classification"] == "Benign")

    report_data = {
        "sample_id": sample_id,
        "date": str(date.today()),
        "variants": variant_records,
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
        "summary": {
            "total": total,
            "pathogenic": pathogenic_count,
            "likely_pathogenic": likely_path_count,
            "drug_response": drug_response_count,
            "risk_factor": risk_factor_count,
            "vus": vus_count,
            "likely_benign": likely_benign_count,
            "benign": benign_count,
        },
        "db_versions": {
            "clinvar": str(date.today()),
            "gnomad": "4.0",
            "krgdb": "2026-03-01",
            "clinvar_local": get_clinvar_db_version(),
        },
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
    _progress(f"\nDone! {total} variants analyzed in {elapsed:.1f}s")
    _progress(
        f"  Pathogenic: {pathogenic_count} | "
        f"Likely Pathogenic: {likely_path_count} | "
        f"Drug Response: {drug_response_count} | "
        f"VUS: {vus_count} | "
        f"Benign: {benign_count + likely_benign_count}"
    )

    return report_data


def main() -> None:
    parser = argparse.ArgumentParser(
        description="GenomeBoard Standalone Pipeline — VCF to Report",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/orchestrate.py data/sample_vcf/demo_variants.vcf
  python scripts/orchestrate.py data/sample_vcf/demo_variants.vcf -o output/report.pdf
  python scripts/orchestrate.py data/sample_vcf/demo_variants.vcf --skip-api --json
        """,
    )
    parser.add_argument("vcf_path", help="Path to input VCF file")
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
        help="Skip external API calls (ClinVar, gnomAD) for offline mode",
    )
    parser.add_argument(
        "--mode",
        choices=["cancer", "rare-disease"],
        default=None,  # will use config default
        help="Report mode: cancer (somatic) or rare-disease (germline)",
    )
    parser.add_argument(
        "--hpo",
        type=str,
        default=None,
        help="Comma-separated HPO IDs for rare disease mode (e.g., HP:0001250,HP:0001263)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose (DEBUG) logging",
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)

    # Resolve JSON output path
    json_output = None
    if args.json_flag is not None:
        if args.json_flag is True:
            # --json with no path: replace extension with .json
            json_output = str(Path(args.output).with_suffix(".json"))
        else:
            json_output = args.json_flag

    # Parse HPO IDs from comma-separated string
    hpo_ids = [h.strip() for h in args.hpo.split(",")] if args.hpo else None

    result = run_pipeline(
        vcf_path=args.vcf_path,
        output_path=args.output,
        krgdb_path=args.krgdb,
        sample_id=args.sample_id,
        json_output=json_output,
        skip_api=args.skip_api,
        mode=args.mode,
        hpo_ids=hpo_ids,
    )

    if result is None:
        sys.exit(1)


if __name__ == "__main__":
    main()
