"""Build structured case briefings from pipeline results.

Converts run_pipeline() output into a well-formatted text document that
reads like a clinical case presentation.  LLM agents reason about this
briefing — it is NOT JSON but human-readable markdown-like text.

The briefing is capped at ~3 000 tokens (roughly 12 000 characters) to
stay within efficient LLM context windows.
"""
from __future__ import annotations

from typing import Optional

# Classification priority for sorting variants (lower = more significant)
_CLASSIFICATION_PRIORITY = {
    "Pathogenic": 0,
    "Likely Pathogenic": 1,
    "Drug Response": 2,
    "Risk Factor": 3,
    "VUS": 4,
    "Likely Benign": 5,
    "Benign": 6,
}

MAX_VARIANTS = 20
MAX_CHARS = 12_000  # ~3 000 tokens


def build_case_briefing(report_data: dict, mode: str) -> str:
    """Convert pipeline results to a structured case briefing for LLM agents.

    Parameters
    ----------
    report_data : dict
        The dict returned by ``run_pipeline()`` (see ``scripts/orchestrate.py``).
    mode : str
        Analysis mode — ``"cancer"`` or ``"rare-disease"``.

    Returns
    -------
    str
        A well-formatted text briefing suitable for LLM consumption.
    """
    if not report_data:
        return _minimal_briefing(mode)

    sections: list[str] = []

    # 1. Case Overview
    sections.append(_build_overview(report_data, mode))

    # 2. Classified Variants
    sections.append(_build_variants_section(report_data, mode))

    # 3. HPO Phenotype Matches (rare disease)
    if mode == "rare-disease":
        hpo_section = _build_hpo_section(report_data)
        if hpo_section:
            sections.append(hpo_section)

    # 4. OMIM Gene-Disease
    if mode == "rare-disease":
        omim_section = _build_omim_section(report_data)
        if omim_section:
            sections.append(omim_section)

    # 5. Pharmacogenomics
    pgx_section = _build_pgx_section(report_data)
    if pgx_section:
        sections.append(pgx_section)

    # 6. Structural Variants
    sv_section = _build_sv_section(report_data)
    if sv_section:
        sections.append(sv_section)

    # 7. TMB (cancer)
    if mode == "cancer":
        tmb_section = _build_tmb_section(report_data)
        if tmb_section:
            sections.append(tmb_section)

    # 8. Summary Statistics
    sections.append(_build_summary_section(report_data))

    briefing = "\n\n".join(s for s in sections if s)

    # Truncate if over budget
    if len(briefing) > MAX_CHARS:
        briefing = briefing[:MAX_CHARS] + "\n\n[...briefing truncated for length]"

    return briefing


# ── Section Builders ─────────────────────────────────────────────────────────

def _minimal_briefing(mode: str) -> str:
    """Return a minimal briefing when report_data is empty."""
    return (
        "= CASE BRIEFING =\n\n"
        f"Mode: {mode}\n"
        "No variant data available for analysis."
    )


def _build_overview(data: dict, mode: str) -> str:
    sample_id = data.get("sample_id", "Unknown")
    analysis_date = data.get("date", "Unknown")
    total = data.get("summary", {}).get("total", 0)
    mode_label = "Cancer Somatic" if mode == "cancer" else "Rare Disease Germline"

    lines = [
        "= CASE BRIEFING =",
        "",
        f"Sample ID: {sample_id}",
        f"Analysis Mode: {mode_label}",
        f"Date: {analysis_date}",
        f"Total Variants Analyzed: {total}",
        "",
        "IMPORTANT NOTES:",
        "- Zygosity (heterozygous/homozygous) is NOT available in this data.",
        "  For autosomal recessive (AR) conditions, do NOT assume affected status",
        "  from a single pathogenic variant — it may represent carrier status only.",
        "  Two pathogenic alleles are required for AR disease diagnosis.",
    ]
    return "\n".join(lines)


def _build_variants_section(data: dict, mode: str) -> str:
    variants = data.get("variants", [])
    if not variants:
        return "== CLASSIFIED VARIANTS ==\n\nNo variants found."

    # Sort by classification priority, then by tier if available
    sorted_variants = sorted(
        variants,
        key=lambda v: (
            _CLASSIFICATION_PRIORITY.get(v.get("classification", "VUS"), 99),
            v.get("tier", "IV"),
        ),
    )

    # Take top N most significant
    display_variants = sorted_variants[:MAX_VARIANTS]
    omitted = len(sorted_variants) - len(display_variants)

    lines = ["== CLASSIFIED VARIANTS ==", ""]

    for i, v in enumerate(display_variants, 1):
        gene = v.get("gene", "Unknown")
        variant_id = v.get("variant", "")
        classification = v.get("classification", "VUS")
        acmg_codes = v.get("acmg_codes", [])
        codes_str = ", ".join(acmg_codes) if acmg_codes else "none"

        lines.append(f"--- Variant {i}: {gene} ---")
        lines.append(f"  Position: {variant_id}")

        hgvsc = v.get("hgvsc", "")
        hgvsp = v.get("hgvsp", "")
        if hgvsc:
            lines.append(f"  HGVSc: {hgvsc}")
        if hgvsp:
            lines.append(f"  HGVSp: {hgvsp}")

        consequence = v.get("consequence", "")
        if consequence:
            lines.append(f"  Consequence: {consequence}")

        lines.append(f"  Classification: {classification}")
        lines.append(f"  Evidence Codes: {codes_str}")

        # Tier (cancer mode)
        if mode == "cancer":
            tier = v.get("tier", "")
            tier_label = v.get("tier_label", "")
            tier_source = v.get("tier_evidence_source", "")
            if tier:
                tier_line = f"  Tier: {tier}"
                if tier_label:
                    tier_line += f" ({tier_label})"
                if tier_source:
                    tier_line += f" — {tier_source}"
                lines.append(tier_line)

        # ClinVar
        clinvar_sig = v.get("clinvar_significance", "")
        if clinvar_sig and clinvar_sig != "Not Found":
            review = v.get("review_status", "")
            clinvar_line = f"  ClinVar: {clinvar_sig}"
            if review:
                clinvar_line += f" ({review})"
            lines.append(clinvar_line)

        # In silico scores
        in_silico = v.get("in_silico", {})
        if in_silico:
            scores = []
            for key in ("revel", "REVEL"):
                val = in_silico.get(key)
                if val is not None:
                    scores.append(f"REVEL={val}")
            for key in ("cadd_phred", "CADD_phred"):
                val = in_silico.get(key)
                if val is not None:
                    scores.append(f"CADD={val}")
            for key in ("am_pathogenicity", "AlphaMissense"):
                val = in_silico.get(key)
                if val is not None:
                    scores.append(f"AlphaMissense={val}")
            for key in ("spliceai_max", "SpliceAI"):
                val = in_silico.get(key)
                if val is not None:
                    scores.append(f"SpliceAI={val}")
            if scores:
                lines.append(f"  In Silico: {', '.join(scores)}")

        # SIFT / PolyPhen (from VEP annotation, separate from in_silico)
        sift = v.get("sift", "")
        polyphen = v.get("polyphen", "")
        if sift or polyphen:
            pred_parts = []
            if sift:
                pred_parts.append(f"SIFT={sift}")
            if polyphen:
                pred_parts.append(f"PolyPhen={polyphen}")
            lines.append(f"  Predictions: {', '.join(pred_parts)}")

        # Population frequencies
        freq_parts = []
        gnomad_all = v.get("gnomad_all")
        gnomad_eas = v.get("gnomad_eas")
        krgdb = v.get("krgdb_freq")
        korea4k = v.get("korea4k_freq")
        nard2 = v.get("nard2_freq")
        if gnomad_all is not None:
            freq_parts.append(f"gnomAD_all={gnomad_all:.5f}")
        if gnomad_eas is not None:
            freq_parts.append(f"gnomAD_EAS={gnomad_eas:.5f}")
        if krgdb is not None:
            freq_parts.append(f"KRGDB={krgdb:.5f}")
        if korea4k is not None:
            freq_parts.append(f"Korea4K={korea4k:.5f}")
        if nard2 is not None:
            freq_parts.append(f"NARD2={nard2:.5f}")
        if freq_parts:
            lines.append(f"  Frequencies: {', '.join(freq_parts)}")

        # HPO match (rare disease)
        if mode == "rare-disease":
            hpo_score = v.get("hpo_score")
            matching_hpo = v.get("matching_hpo", [])
            if hpo_score and hpo_score > 0:
                lines.append(f"  HPO Match Score: {hpo_score}")
            if matching_hpo:
                hpo_strs = []
                for h in matching_hpo[:3]:
                    if isinstance(h, dict):
                        hpo_strs.append(f"{h.get('id', '')}: {h.get('name', '')}")
                    else:
                        hpo_strs.append(str(h))
                lines.append(f"  Matching HPO: {'; '.join(hpo_strs)}")

        lines.append("")

    if omitted > 0:
        lines.append(f"[{omitted} additional lower-priority variants omitted]")
        lines.append("")

    return "\n".join(lines)


def _build_hpo_section(data: dict) -> Optional[str]:
    hpo_results = data.get("hpo_results", [])
    if not hpo_results:
        return None

    lines = ["== HPO PHENOTYPE MATCHES ==", ""]
    for hpo in hpo_results:
        hpo_id = hpo.get("id", "")
        name = hpo.get("name", "")
        genes = hpo.get("genes", [])
        gene_str = ", ".join(genes[:10]) if genes else "none"
        lines.append(f"  {hpo_id}: {name}")
        lines.append(f"    Associated genes: {gene_str}")
        if len(genes) > 10:
            lines.append(f"    (+{len(genes) - 10} more genes)")
        # Warn about non-specific broad HPO terms
        if len(genes) > 500:
            lines.append(f"    ⚠ WARNING: This HPO term is extremely broad ({len(genes)} genes).")
            lines.append(f"      HPO matches for this term have LOW specificity.")
            lines.append(f"      Do NOT use this match alone to support gene-disease correlation.")
    lines.append("")
    return "\n".join(lines)


def _build_omim_section(data: dict) -> Optional[str]:
    variants = data.get("variants", [])
    omim_entries = []
    for v in variants:
        omim_phenos = v.get("omim_phenotypes", [])
        inheritance = v.get("inheritance", "")
        if omim_phenos or inheritance:
            omim_entries.append({
                "gene": v.get("gene", "Unknown"),
                "phenotypes": omim_phenos,
                "inheritance": inheritance,
                "mim": v.get("omim_mim", ""),
            })

    if not omim_entries:
        return None

    lines = ["== OMIM GENE-DISEASE ASSOCIATIONS ==", ""]
    for entry in omim_entries:
        lines.append(f"  {entry['gene']} (MIM: {entry['mim']})")
        if entry["inheritance"]:
            lines.append(f"    Inheritance: {entry['inheritance']}")
        for pheno in entry["phenotypes"][:3]:
            if isinstance(pheno, dict):
                lines.append(f"    - {pheno.get('phenotype', pheno)}")
            else:
                lines.append(f"    - {pheno}")
    lines.append("")
    return "\n".join(lines)


def _build_pgx_section(data: dict) -> Optional[str]:
    pgx_results = data.get("pgx_results", [])
    if not pgx_results:
        return None

    lines = ["== PHARMACOGENOMICS ==", ""]
    for pgx in pgx_results:
        gene = pgx.get("gene", "Unknown")
        star = pgx.get("star_allele", "")
        pheno = pgx.get("phenotype", "")
        cpic = pgx.get("cpic_level", "")
        recommendation = pgx.get("cpic_recommendation", "")
        korean_flag = pgx.get("korean_flag", False)

        lines.append(f"  {gene} {star}: {pheno}")
        lines.append(f"    CPIC Level: {cpic}")
        if recommendation:
            lines.append(f"    Recommendation: {recommendation}")
        if korean_flag:
            lines.append("    * Korean population frequency significantly differs from Western populations")
    lines.append("")
    return "\n".join(lines)


def _build_sv_section(data: dict) -> Optional[str]:
    sv_class45 = data.get("sv_class45", [])
    sv_class3 = data.get("sv_class3_display", [])
    if not sv_class45 and not sv_class3:
        return None

    lines = ["== STRUCTURAL VARIANTS ==", ""]
    for sv in sv_class45:
        sv_type = sv.get("sv_type", "")
        gene = sv.get("gene_name", "Unknown")
        acmg_label = sv.get("acmg_label", "")
        size = sv.get("size_display", "")
        cytoband = sv.get("cytoband", "")
        lines.append(f"  {sv_type} — {gene} ({cytoband})")
        lines.append(f"    Classification: {acmg_label}, Size: {size}")
        pheno = sv.get("phenotypes", "")
        if pheno:
            lines.append(f"    Phenotype: {pheno}")

    for sv in sv_class3[:5]:
        sv_type = sv.get("sv_type", "")
        gene = sv.get("gene_name", "Unknown")
        size = sv.get("size_display", "")
        lines.append(f"  {sv_type} — {gene} (VUS, {size})")

    lines.append("")
    return "\n".join(lines)


def _build_tmb_section(data: dict) -> Optional[str]:
    tmb = data.get("tmb")
    if not tmb:
        return None

    score = tmb.get("score", 0)
    level = tmb.get("level", "Unknown")
    variant_count = tmb.get("variant_count", 0)
    total = tmb.get("total_variants", 0)
    panel_size = tmb.get("panel_size_mb", 0)

    lines = [
        "== TUMOR MUTATIONAL BURDEN ==",
        "",
        f"  TMB Score: {score:.1f} mutations/Mb",
        f"  Level: {level}",
        f"  Counted Variants: {variant_count}/{total} (nonsynonymous)",
        f"  Panel Size: {panel_size} Mb",
        "",
    ]
    return "\n".join(lines)


def _build_summary_section(data: dict) -> str:
    summary = data.get("summary", {})
    if not summary:
        return "== SUMMARY ==\n\nNo summary data available."

    lines = [
        "== SUMMARY STATISTICS ==",
        "",
        f"  Total variants: {summary.get('total', 0)}",
        f"  Pathogenic: {summary.get('pathogenic', 0)}",
        f"  Likely Pathogenic: {summary.get('likely_pathogenic', 0)}",
        f"  VUS: {summary.get('vus', 0)}",
        f"  Likely Benign: {summary.get('likely_benign', 0)}",
        f"  Benign: {summary.get('benign', 0)}",
        f"  Drug Response: {summary.get('drug_response', 0)}",
        f"  Risk Factor: {summary.get('risk_factor', 0)}",
    ]
    return "\n".join(lines)
