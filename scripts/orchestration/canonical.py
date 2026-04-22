"""Canonical single-sample pipeline.

This module is the single source of truth for building a report_data dict
from one VCF. ``scripts.orchestrate.run_pipeline`` is a thin wrapper that
handles CLI concerns (progress logging, PED parsing, output path resolution,
JSON dump), then delegates the assembly work to ``build_sample_report``.
``scripts.orchestration.batch.run_batch_pipeline`` calls the same function
once per sample.

Design contract:
    * Classification stays deterministic (no LLM on this path).
    * No mutation of caller state except a scoped ``report_data`` dict that
      is returned.
    * Variant objects are parsed per sample — batch mode must not reuse
      Variant instances across samples (v2.5.4 H1 bug).
    * CIViC + gene_knowledge enrichment is applied here (in assembly), not
      during render. Render stays pure from v2.5.4 onward (M2 / L7 fixes).
"""

from __future__ import annotations

import logging
import re
from dataclasses import asdict, is_dataclass
from datetime import date
from pathlib import Path
from typing import Any

from scripts.common.config import get
from scripts.common.gene_knowledge import get_gene_info
from scripts.common.hgvs_utils import hgvsp_to_civic_variant as _hgvsp_to_civic_variant
from scripts.common.models import FrequencyData
from scripts.enrichment.hpo_matcher import resolve_hpo_terms
from scripts.intake.parse_vcf import parse_vcf
from scripts.orchestration.classify import (
    build_summary,
    build_variant_records,
    classify_variants,
    split_variants_for_display,
    sv_to_dict,
)
from scripts.orchestration.query import query_variant_databases
from scripts.population.compare_freq import compare_frequencies
from scripts.storage.query_civic import get_gene_summary, get_treatment_summary, get_variant_evidence
from scripts.storage.version_manager import get_all_db_versions

logger = logging.getLogger(__name__)


# ── Narrative shaping helpers (moved from generate_pdf for render purity) ────


def _adjust_finding_summary(summary: str, classification: str) -> str:
    """Adjust gene knowledge finding_summary to match actual variant classification.

    Replacement order matters: replace longer phrases first to avoid
    partial matches (e.g. "A likely pathogenic variant" before "A pathogenic variant").
    """
    cls_lower = classification.lower()
    vus_label = "A variant of uncertain significance (VUS)"
    if "vus" in cls_lower or "uncertain" in cls_lower:
        summary = summary.replace("A likely pathogenic variant", vus_label)
        summary = summary.replace("A pathogenic variant", vus_label)
        summary = summary.replace("A pharmacogenomic variant", vus_label)
        if "further evidence is needed" not in summary.lower():
            summary = summary.replace(
                "was identified in this specimen",
                "was identified in this specimen. Further evidence is needed to determine clinical significance",
            )
    elif "benign" in cls_lower:
        summary = summary.replace("A likely pathogenic variant", "A benign variant")
        summary = summary.replace("A pathogenic variant", "A benign variant")
    elif "drug response" in cls_lower:
        summary = summary.replace("A likely pathogenic variant", "A pharmacogenomic variant")
        summary = summary.replace("A pathogenic variant", "A pharmacogenomic variant")
    elif "risk factor" in cls_lower:
        summary = summary.replace("A likely pathogenic variant", "A risk factor variant")
        summary = summary.replace("A pathogenic variant", "A risk factor variant")
    elif "likely pathogenic" in cls_lower:
        summary = summary.replace("A pathogenic variant", "A likely pathogenic variant")
    return summary


def _linkify_pmids(text: str) -> str:
    """Convert PMID references in text to PubMed links."""
    if not text:
        return text

    def _replace_pmid(match: re.Match) -> str:
        pmid = match.group(1)
        return (
            f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank" '
            f'rel="noopener" style="color:#1d4ed8;text-decoration:none;">PMID:{pmid}</a>'
        )

    return re.sub(r"PMID:\s*(\d+)", _replace_pmid, text)


def _build_frequency_text(v: dict) -> str:
    """Assemble a human-readable frequency summary from available data.

    Extracted from ``generate_report_html``; runs once per variant at
    assembly time rather than on every render.
    """
    freq_parts: list[str] = []
    if v.get("gnomad_all") is not None:
        freq_parts.append(f"gnomAD global AF: {v['gnomad_all']:.6f}")
    if v.get("gnomad_eas") is not None:
        freq_parts.append(f"gnomAD East Asian AF: {v['gnomad_eas']:.6f}")
    kf = ""
    if v.get("agents") and v["agents"].get("korean_pop"):
        kf = v["agents"]["korean_pop"].get("korean_flag", "")
    if kf and kf not in ("No notable findings", "No frequency data available"):
        freq_parts.append(f"Korean: {kf}")
    if v.get("kova_freq") is not None:
        freq_parts.append(f"KOVA Korean AF: {v['kova_freq']}")
    if v.get("kova_homozygote") is not None:
        freq_parts.append(f"KOVA homozygote count: {v['kova_homozygote']}")
    if freq_parts:
        return ". ".join(freq_parts) + "."
    if v.get("gnomad_all") is None:
        return "Not observed in gnomAD population database (absent or extremely rare)."
    return ""


# ── Shared helpers ───────────────────────────────────────────────────────────


def normalize_sample_id(vcf_path: str) -> str:
    """Canonical sample_id from a VCF filename.

    Strips ``.vcf``, ``.vcf.gz``, ``.vcf.bgz`` in that order of preference
    (longest first) and preserves the original case. This replaces the
    v2.5.3 ``Path(vcf_path).stem.upper()`` behaviour and the batch-mode
    ad-hoc extension stripping, which disagreed on the same file.

    Example:
        >>> normalize_sample_id("/data/patient_001.vcf.gz")
        'patient_001'
    """
    name = Path(vcf_path).name
    for ext in (".vcf.bgz", ".vcf.gz", ".vcf"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


# ── Pipeline stages ──────────────────────────────────────────────────────────


def _parse_variants(vcf_path: str, ped_path: str | None) -> list:
    """Parse a VCF into ``Variant`` objects. Returns empty list on failure."""
    vcf_file = Path(vcf_path)
    if not vcf_file.exists():
        logger.error("VCF file not found: %s", vcf_path)
        return []
    variants = parse_vcf(str(vcf_file), ped_path=ped_path)
    return list(variants) if variants else []


def _query_all(variants: list, skip_api: bool) -> dict:
    """Query ClinVar/gnomAD/KOVA/PGx for every variant. Returns {variant_id: dict}."""
    db_results: dict[str, Any] = {}
    for variant in variants:
        db_results[variant.variant_id] = query_variant_databases(variant, skip_api)
    return db_results


def _extract_germline_if_applicable(
    germline_vcf: str | None,
    mode: str,
    primary_variants: list,
    db_results: dict,
    skip_api: bool,
) -> list:
    """Extract inherited target variants from germline VCF (rare-disease only).

    Returns the list of inherited Variants (may be empty). Side effect:
    populates db_results with their annotations.
    """
    if not germline_vcf or mode != "rare-disease":
        return []
    try:
        from scripts.orchestration.extract_germline import extract_inherited_variants

        primary_ids = {v.variant_id for v in primary_variants}
        inherited = extract_inherited_variants(
            germline_vcf,
            primary_variant_ids=primary_ids,
        )
        if not inherited:
            return []
        for iv in inherited:
            db_results[iv.variant_id] = query_variant_databases(iv, skip_api)
        return list(inherited)
    except Exception as e:
        logger.warning("Germline inherited variant extraction failed: %s", e)
        return []


def _build_freq_results(variants: list, db_results: dict) -> dict:
    """Run frequency comparison for every variant."""
    freq_results: dict[str, Any] = {}
    for variant in variants:
        db = db_results[variant.variant_id]
        freq_data = FrequencyData(
            kova=db.get("kova_freq"),
            gnomad_eas=db["gnomad"].get("gnomad_eas"),
            gnomad_all=db["gnomad"].get("gnomad_all"),
            kova_homozygote=db.get("kova_homozygote"),
        )
        freq_results[variant.variant_id] = compare_frequencies(freq_data)
    return freq_results


def _run_pgx(variants: list, germline_vcf: str | None, db_results: dict) -> dict:
    """Run pharmacogenomics (PharmCAT if germline provided, builtin fallback)."""
    from scripts.pharmacogenomics.korean_pgx import get_pgx_results

    pgx_data = get_pgx_results(variants, germline_vcf=germline_vcf)
    pgx_source = pgx_data["pgx_source"]

    # Backward-compat: collect PgxResult objects from db_results when builtin
    pgx_hits = []
    if pgx_source in ("builtin", "builtin_limited"):
        for variant in variants:
            pgx = db_results[variant.variant_id].get("pgx")
            if pgx:
                pgx_hits.append(pgx)

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
    return {
        "pgx_results_list": pgx_results_list,
        "pgx_source": pgx_source,
        "germline_provided": pgx_data["germline_provided"],
        "pharmcat_version": pgx_data["pharmcat_version"],
        "warnings": pgx_data.get("warnings", []),
    }


def _classify_all(
    variants: list,
    db_results: dict,
    freq_results: dict,
    intervar_path: str | None,
) -> dict:
    """Classify every variant (ACMG + optional InterVar evidence codes)."""
    intervar_data = None
    if intervar_path:
        try:
            from scripts.annotation.parse_intervar import parse_intervar

            intervar_data = parse_intervar(intervar_path)
        except Exception as e:
            logger.warning("InterVar parsing failed: %s", e)

    return classify_variants(variants, db_results, freq_results, intervar_data=intervar_data)


# ── Post-classification enrichment (M2 + L7 fix) ─────────────────────────────


def _enrich_with_gene_knowledge_and_civic(variant_records: list[dict], mode: str) -> list[dict]:
    """Fill gene_knowledge + CIViC fields on every variant record.

    This was previously performed inside ``generate_report_html``. Moving it
    here serves two goals:

      * **Assembly vs render separation (L7)** — the renderer should not
        mutate its input. Equivalence tests (``test_render_does_not_mutate_input``)
        now assert byte-for-byte stability of the dict passed to
        ``generate_report_html``.
      * **CIViC priority (M2)** — the previous implementation chained
        ``setdefault()`` calls, which meant gene_knowledge values (set
        first) always won over CIViC values (set second). The spec and
        the v2.2 curate-then-narrate contract require the opposite.

    Cancer-mode behaviour:
      * ``finding_summary`` — CIViC description (truncated) if present,
        else gene_knowledge, else existing value.
      * ``treatment_strategies`` — CIViC treatment summary if present,
        else gene_knowledge, else existing value.
      * ``references`` + ``content_status`` — CIViC evidence rows set
        ``content_status="curated-civic"``; when absent we fall back to
        gene_knowledge (``ai-generated``).
      * Non-clinical fields (gene_full_name, associated_conditions,
        korean_specific_note, frequency_prognosis, hgvs) are populated
        from gene_knowledge only; CIViC doesn't carry them.

    Rare-disease mode: CIViC is skipped entirely (OMIM / ClinGen / HPO
    drive the enrichment via ``build_variant_records`` already). Only
    gene_knowledge fallback fields are filled here.

    The dicts are mutated in place (cheaper than deep-copying every
    record) and the same list is returned for a fluent API.
    """
    for v in variant_records:
        gene = v.get("gene")
        if not gene:
            continue

        gk = get_gene_info(gene) or {}

        civic_gene = None
        civic_treatment = None
        civic_evidence: list = []
        if mode == "cancer":
            try:
                civic_gene = get_gene_summary(gene)
            except Exception as e:
                logger.debug("CIViC gene summary lookup failed for %s: %s", gene, e)
            try:
                hgvsp = v.get("hgvsp", "")
                civic_variant_name = _hgvsp_to_civic_variant(hgvsp) if hgvsp else None
                civic_treatment = get_treatment_summary(gene, civic_variant_name)
                if civic_variant_name:
                    civic_evidence = get_variant_evidence(gene, civic_variant_name) or []
                if not civic_evidence:
                    civic_evidence = get_variant_evidence(gene) or []
            except Exception as e:
                logger.debug("CIViC treatment / evidence lookup failed for %s: %s", gene, e)

        # ── Clinical fields ──────────────────────────────────────────────────
        # Precedence (M2 fix): explicit caller-provided value > CIViC
        # > gene_knowledge > empty string. The pre-v2.5.4 setdefault chain
        # gave gene_knowledge the win over CIViC — this block flips the
        # gk-vs-civic tie in CIViC's favour without dropping explicit values
        # the Board / curator already wrote to the record.
        civic_summary = ""
        if civic_gene and civic_gene.get("description"):
            civic_summary = civic_gene["description"][:500]
        if not v.get("finding_summary"):
            v["finding_summary"] = civic_summary or gk.get("finding_summary", "")

        if not v.get("treatment_strategies"):
            v["treatment_strategies"] = civic_treatment or gk.get("treatment_strategies", "")

        # References: CIViC evidence rows replace (tagged with
        # ``content_status="curated-civic"``); otherwise fall back to
        # gene_knowledge. Caller-provided references are always kept.
        # ``content_status`` semantics (test_references coverage):
        #   - caller-supplied → preserved verbatim (plain "ai-generated"
        #     watermark survives even when references slot is empty)
        #   - no caller value + CIViC evidence attached → "curated-civic"
        #   - no caller value + gene_knowledge has references → inherit
        #     gk.content_status if present
        #   - no caller value + no CIViC + no gk data → leave unset so
        #     the watermark div stays out of the template (FAKEGENE case)
        caller_supplied_content_status = "content_status" in v and v.get("content_status")
        if not v.get("references"):
            if civic_evidence:
                refs = [
                    {
                        "pmid": e["pmid"],
                        "source": e["citation"],
                        "note": f"{e['evidence_type']} — {e['significance']}",
                    }
                    for e in civic_evidence[:5]
                    if e.get("pmid")
                ]
                if refs:
                    v["references"] = refs
                    if not caller_supplied_content_status:
                        v["content_status"] = "curated-civic"
            if not v.get("references"):
                gk_refs = gk.get("references", [])
                if gk_refs:
                    v["references"] = gk_refs
        if not v.get("content_status") and "content_status" not in v:
            gk_status = gk.get("content_status")
            if gk_status:
                v["content_status"] = gk_status

        # ── Non-clinical fallbacks (gene_knowledge only) ─────────────────────
        if not v.get("frequency_prognosis"):
            v["frequency_prognosis"] = gk.get("frequency_prognosis", "")
        if not v.get("gene_full_name"):
            v["gene_full_name"] = gk.get("full_name", "")
        if not v.get("associated_conditions"):
            v["associated_conditions"] = gk.get("associated_conditions", [])
        if not v.get("korean_specific_note"):
            v["korean_specific_note"] = gk.get("korean_specific_note")

        # ── HGVS assembly — prefer VCF annotation, gene_knowledge fallback ───
        if v.get("hgvsc") or v.get("hgvsp"):
            if not v.get("hgvs"):
                v["hgvs"] = {
                    "transcript": v.get("transcript", ""),
                    "cdna": v.get("hgvsc", ""),
                    "protein": v.get("hgvsp", ""),
                    "variant_effect": v.get("consequence", ""),
                }
            if not v.get("variant_effect"):
                v["variant_effect"] = v.get("consequence", "")
        else:
            if not v.get("hgvs"):
                v["hgvs"] = gk.get("hgvs", {})
            if not v.get("variant_effect"):
                v["variant_effect"] = gk.get("hgvs", {}).get("variant_effect", "")

        # ── Frequency narrative ──────────────────────────────────────────────
        if not v.get("frequency_prognosis"):
            v["frequency_prognosis"] = _build_frequency_text(v)

        # ── AMP tier defaults ────────────────────────────────────────────────
        v.setdefault("tier_evidence_source", "")
        v.setdefault("civic_match_level", "none")

        # ── Classification-aware finding_summary adjustment + PMID linkify ──
        classification = v.get("classification", "VUS")
        raw_summary = v.get("finding_summary", "")
        if raw_summary:
            v["finding_summary"] = _adjust_finding_summary(raw_summary, classification)

        if v.get("treatment_strategies") and "<a href=" not in v["treatment_strategies"]:
            v["treatment_strategies"] = _linkify_pmids(v["treatment_strategies"])
        if v.get("finding_summary") and "<a href=" not in v["finding_summary"]:
            v["finding_summary"] = _linkify_pmids(v["finding_summary"])

    return variant_records


def _enrich_pgx_with_gene_knowledge(pgx_results: list[dict]) -> list[dict]:
    """Fill HGVS / variant_effect on PGx rows from gene_knowledge (defaults).

    Moved from generate_report_html for render purity. Mutates and returns
    the same list.
    """
    for pgx in pgx_results:
        gene = pgx.get("gene")
        if not gene:
            continue
        info = get_gene_info(gene)
        if info:
            pgx.setdefault("hgvs", info.get("hgvs", {}))
            pgx.setdefault("variant_effect", info.get("hgvs", {}).get("variant_effect", ""))
    return pgx_results


# ── Top-level entry point ────────────────────────────────────────────────────


def build_sample_report(
    vcf_path: str,
    *,
    sample_id: str | None = None,
    mode: str = "cancer",
    skip_api: bool = False,
    hpo_ids: list | None = None,
    hide_vus: bool = True,
    sv_path: str | None = None,
    panel_size_mb: float | None = None,
    bed_path: str | None = None,
    intervar_path: str | None = None,
    clinical_board: bool = False,
    board_lang: str | None = None,
    clinical_note: str | None = None,
    germline_vcf: str | None = None,
    ped_path: str | None = None,
) -> dict | None:
    """Build the full ``report_data`` dict for a single VCF.

    This is the canonical single-sample pipeline. Both ``run_pipeline``
    (CLI wrapper) and ``run_batch_pipeline`` (per-sample loop) route
    through here. Returns ``None`` on fatal input failure (missing VCF,
    no variants parsed) — callers should treat that as a skip.
    """
    # ── Sample id + HPO setup ─────────────────────────────────────────────
    if sample_id is None:
        sample_id = normalize_sample_id(vcf_path)

    hpo_results: list = []
    if mode == "rare-disease" and hpo_ids:
        hpo_results = resolve_hpo_terms(hpo_ids)

    # ── Parse primary VCF ─────────────────────────────────────────────────
    variants = _parse_variants(vcf_path, ped_path)
    if not variants:
        logger.error("No variants parsed from %s", vcf_path)
        return None

    # ── Query databases ───────────────────────────────────────────────────
    db_results = _query_all(variants, skip_api)

    # ── Germline inherited extraction (rare-disease only) ─────────────────
    inherited = _extract_germline_if_applicable(germline_vcf, mode, variants, db_results, skip_api)
    if inherited:
        variants = variants + inherited

    # ── Frequency comparison ──────────────────────────────────────────────
    freq_results = _build_freq_results(variants, db_results)

    # ── PGx ───────────────────────────────────────────────────────────────
    pgx_bundle = _run_pgx(variants, germline_vcf, db_results)
    # Phase 2 (L7): enrich PGx rows here, not in render.
    _enrich_pgx_with_gene_knowledge(pgx_bundle["pgx_results_list"])

    # ── ACMG classification ───────────────────────────────────────────────
    classification_results = _classify_all(variants, db_results, freq_results, intervar_path)

    # ── Assemble variant records ──────────────────────────────────────────
    variant_records = build_variant_records(
        variants, db_results, freq_results, classification_results, mode, hpo_results
    )
    # v2.5.4 Phase 2 (M2 + L7 fix): enrichment happens here at assembly
    # time rather than inside generate_report_html. CIViC data now wins
    # over gene_knowledge for clinical fields (finding_summary,
    # treatment_strategies, references, content_status) — the previous
    # setdefault chain had the priorities inverted.
    _enrich_with_gene_knowledge_and_civic(variant_records, mode)

    summary = build_summary(variant_records)
    tier1, tier2, tier3, tier4_count, detailed_variants, omitted_variants = split_variants_for_display(
        variant_records, hide_vus
    )

    report_data: dict[str, Any] = {
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
        "pgx_results": pgx_bundle["pgx_results_list"],
        "pgx_source": pgx_bundle["pgx_source"],
        "germline_provided": pgx_bundle["germline_provided"],
        "pharmcat_version": pgx_bundle["pharmcat_version"],
        "summary": summary,
        "db_versions": get_all_db_versions(skip_api=skip_api),
        "pipeline": {
            "skip_api": skip_api,
            "ped_used": bool(ped_path),
            "ped_path": str(ped_path) if ped_path else None,
        },
        "mode": mode,
        "hpo_results": hpo_results,
    }

    # ── Structural variants (AnnotSV) ─────────────────────────────────────
    sv_variants: list = []
    if sv_path:
        try:
            from scripts.intake.parse_annotsv import parse_annotsv

            sv_variants = parse_annotsv(sv_path)
        except Exception as e:
            logger.warning("AnnotSV parse failed: %s", e)

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

    # ── TMB (cancer mode only) ────────────────────────────────────────────
    if mode == "cancer":
        try:
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
        except Exception as e:
            logger.warning("TMB calculation failed: %s", e)
            report_data["tmb"] = None
    else:
        report_data["tmb"] = None

    # ── Clinical Board (optional LLM diagnostic synthesis) ────────────────
    if clinical_note:
        report_data["clinical_note"] = clinical_note
    if clinical_board:
        try:
            from scripts.clinical_board.render import render_board_opinion_html
            from scripts.clinical_board.runner import run_clinical_board

            board_opinion = run_clinical_board(report_data, mode, language=board_lang)
            if board_opinion:
                report_data["clinical_board"] = board_opinion
                report_data["clinical_board_html"] = render_board_opinion_html(
                    board_opinion, language=board_lang or get("clinical_board.language", "en")
                )
        except Exception as e:
            logger.warning("Clinical Board failed: %s", e)

    return report_data


def report_data_for_json(report_data: dict) -> dict:
    """Return a version of ``report_data`` safe for ``json.dumps``.

    Unwraps the ``clinical_board`` dataclass so ``rerender_report`` can
    round-trip it back into a dataclass instance.
    """
    cb = report_data.get("clinical_board")
    if cb is not None and is_dataclass(cb):
        return {**report_data, "clinical_board": asdict(cb)}
    return report_data
