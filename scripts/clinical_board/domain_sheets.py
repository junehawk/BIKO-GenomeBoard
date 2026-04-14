"""Build per-agent, per-mode domain sheets from local DB data.

Each builder consumes the variants/report_data already collected by the
pipeline and formats a focused, human-readable text section for one agent.
No new DB queries are issued here — the pipeline has already collected the
relevant fields.

**Caller contract:** ``variants`` MUST be the board-selected variant list
produced by ``scripts.clinical_board.variant_selector.select_board_variants``,
not the raw ``report_data["variants"]`` from the pipeline. The domain sheets
are fed directly to per-domain LLM agents, and showing unfiltered passenger
variants to the agents erodes the clinical criteria applied in the selector.
``runner.py`` is responsible for calling the selector once and passing the
filtered list to every ``build_domain_sheet`` invocation.
"""

from __future__ import annotations

from typing import Callable

MAX_DOMAIN_CHARS = 16_000  # ~4K tokens at ~4 chars/token

_TRUNCATION_MARKER = "\n[...truncated for context limit]"


def build_domain_sheet(
    agent_domain: str,
    mode: str,
    variants: list,
    report_data: dict,
) -> str:
    """Build a domain-specific data sheet for one agent.

    Args:
        agent_domain: e.g. "variant_pathology", "disease_genetics".
        mode: "cancer" or "rare-disease".
        variants: list of variant dicts already enriched by the pipeline.
        report_data: the assembled report context (TMB, sample meta, etc).
    """
    if mode == "cancer":
        builders = _CANCER_BUILDERS
    else:
        builders = _RARE_DISEASE_BUILDERS

    builder = builders.get(agent_domain, _empty_sheet)
    sheet = builder(variants, report_data)

    if len(sheet) > MAX_DOMAIN_CHARS:
        keep = MAX_DOMAIN_CHARS - len(_TRUNCATION_MARKER)
        sheet = sheet[:keep] + _TRUNCATION_MARKER
    return sheet


def _empty_sheet(variants: list, report_data: dict) -> str:
    return ""


# ---------------------------------------------------------------------------
# Rare disease builders
# ---------------------------------------------------------------------------


def _rd_variant_pathologist(variants: list, report_data: dict) -> str:
    if not variants:
        return "## Variant Pathology Domain Data\n(No variants supplied.)\n"

    lines = ["## Variant Pathology Domain Data", ""]
    for v in variants:
        gene = v.get("gene", "?")
        classification = v.get("classification", "?")
        hgvsp = v.get("hgvsp") or v.get("hgvsc") or ""
        clinvar_sig = v.get("clinvar_significance", "n/a")
        review = v.get("review_status", "n/a")
        in_silico = v.get("in_silico", {}) or {}
        revel = in_silico.get("revel", "n/a")
        cadd = in_silico.get("cadd_phred", "n/a")
        spliceai = in_silico.get("spliceai_max", "n/a")

        lines.append(f"- {gene} {hgvsp}".rstrip())
        lines.append(f"    Classification: {classification}")
        lines.append(f"    ClinVar: {clinvar_sig} (review: {review})")
        lines.append(f"    In silico: REVEL={revel}, CADD={cadd}, SpliceAI={spliceai}")
        evidence = v.get("acmg_evidence") or v.get("evidence_codes") or []
        if evidence:
            lines.append(f"    ACMG codes: {', '.join(map(str, evidence))}")
        lines.append("")
    return "\n".join(lines)


def _rd_disease_geneticist(variants: list, report_data: dict) -> str:
    if not variants:
        return "## Disease Genetics Domain Data\n(No variants supplied.)\n"

    lines = ["## Disease Genetics Domain Data", ""]
    for v in variants:
        gene = v.get("gene", "?")
        omim_phenotypes = v.get("omim_phenotypes") or []
        inheritance = v.get("inheritance", "n/a")
        matching_hpo = v.get("matching_hpo") or []
        gene_reviews = v.get("gene_reviews_summary") or v.get("gene_reviews", "")
        clingen_haplo = v.get("clingen_haploinsufficiency") or v.get("clingen", "")

        lines.append(f"- {gene}")
        lines.append(f"    Inheritance: {inheritance}")
        if omim_phenotypes:
            phen_strs = [
                p if isinstance(p, str) else f"{p.get('phenotype', p.get('name', '?'))}" for p in omim_phenotypes
            ]
            lines.append(f"    OMIM phenotypes: {'; '.join(phen_strs)}")
        if matching_hpo:
            hpo_strs = [
                h if isinstance(h, str) else f"{h.get('id', '?')} {h.get('name', '')}".strip() for h in matching_hpo
            ]
            lines.append(f"    Matching HPO: {', '.join(hpo_strs)}")
        if gene_reviews:
            lines.append(f"    GeneReviews: {gene_reviews}")
        if clingen_haplo:
            lines.append(f"    ClinGen: {clingen_haplo}")
        lines.append("")
    return "\n".join(lines)


def _rd_pgx_specialist(variants: list, report_data: dict) -> str:
    pgx_results = report_data.get("pgx_results") or report_data.get("pgx") or []
    if not pgx_results:
        return "## Pharmacogenomics Domain Data\n(No PGx results.)\n"

    lines = ["## Pharmacogenomics Domain Data", ""]
    for r in pgx_results:
        gene = r.get("gene", "?")
        diplotype = r.get("diplotype") or r.get("genotype", "?")
        phenotype = r.get("phenotype", "?")
        cpic = r.get("cpic_recommendation") or r.get("recommendation", "")
        kor_freq = r.get("korean_frequency") or r.get("kor_freq", "")
        lines.append(f"- {gene} {diplotype} — {phenotype}")
        if cpic:
            lines.append(f"    CPIC: {cpic}")
        if kor_freq:
            lines.append(f"    Korean prevalence: {kor_freq}")
        lines.append("")
    return "\n".join(lines)


def _rd_literature_analyst(variants: list, report_data: dict) -> str:
    if not variants:
        return "## Literature Evidence Domain Data\n(No variants supplied.)\n"

    lines = ["## Literature Evidence Domain Data", ""]
    for v in variants:
        gene = v.get("gene", "?")
        hgvsp = v.get("hgvsp") or v.get("hgvsc") or ""
        civic = v.get("civic_evidence") or []
        clingen = v.get("clingen_validity") or ""
        gene_reviews = v.get("gene_reviews_summary") or ""
        pmids = v.get("pmids") or v.get("references") or []

        lines.append(f"- {gene} {hgvsp}".rstrip())
        if civic:
            for ev in civic[:5]:
                drug = ev.get("drug", "")
                level = ev.get("level", "")
                direction = ev.get("direction", "")
                lines.append(f"    CIViC: {drug} (level {level}, {direction})".rstrip())
        if clingen:
            lines.append(f"    ClinGen validity: {clingen}")
        if gene_reviews:
            lines.append(f"    GeneReviews: {gene_reviews}")
        if pmids:
            lines.append(f"    Refs: {', '.join(map(str, pmids[:8]))}")
        lines.append("")
    return "\n".join(lines)


_RARE_DISEASE_BUILDERS: dict[str, Callable[[list, dict], str]] = {
    "variant_pathology": _rd_variant_pathologist,
    "disease_genetics": _rd_disease_geneticist,
    "pharmacogenomics": _rd_pgx_specialist,
    "literature_evidence": _rd_literature_analyst,
}


# ---------------------------------------------------------------------------
# Cancer builders
# ---------------------------------------------------------------------------


def _cancer_therapeutic_target(variants: list, report_data: dict) -> str:
    if not variants:
        return "## Therapeutic Target Domain Data\n(No variants supplied.)\n"

    lines = ["## Therapeutic Target Domain Data", ""]
    for v in variants:
        gene = v.get("gene", "?")
        hgvsp = v.get("hgvsp") or v.get("hgvsc") or ""
        classification = v.get("classification", "?")
        clinvar_sig = v.get("clinvar_significance", "n/a")
        in_silico = v.get("in_silico", {}) or {}
        revel = in_silico.get("revel", "n/a")
        cadd = in_silico.get("cadd_phred", "n/a")
        civic = v.get("civic_evidence") or []
        oncokb = v.get("oncokb") or {}

        lines.append(f"- {gene} {hgvsp}".rstrip())
        lines.append(f"    Classification: {classification} (ClinVar: {clinvar_sig})")
        lines.append(f"    In silico: REVEL={revel}, CADD={cadd}")
        if oncokb:
            level = oncokb.get("level") or oncokb.get("therapeutic_level", "")
            lines.append(f"    OncoKB: {level}".rstrip())
        if civic:
            for ev in civic[:5]:
                drug = ev.get("drug", "")
                level = ev.get("level", "")
                direction = ev.get("direction", "")
                lines.append(f"    CIViC drug: {drug} (level {level}, {direction})".rstrip())
        resistance = v.get("resistance_notes") or v.get("resistance", "")
        if resistance:
            lines.append(f"    Resistance: {resistance}")
        lines.append("")
    return "\n".join(lines)


def _cancer_tumor_genomics(variants: list, report_data: dict) -> str:
    lines = ["## Tumor Genomics Domain Data", ""]
    tmb = report_data.get("tmb") or {}
    if tmb:
        score = tmb.get("score", "n/a")
        level = tmb.get("level", "n/a")
        lines.append(f"TMB: {score} mut/Mb ({level})")
    msi = report_data.get("msi") or {}
    if msi:
        lines.append(f"MSI: {msi.get('status', 'n/a')}")
    lines.append("")
    if not variants:
        lines.append("(No variants supplied.)")
        return "\n".join(lines)

    for v in variants:
        gene = v.get("gene", "?")
        hgvsp = v.get("hgvsp") or v.get("hgvsc") or ""
        vaf = v.get("vaf", "n/a")
        hotspot = v.get("hotspot") or v.get("is_hotspot")
        driver = v.get("driver_status") or v.get("driver", "")
        lines.append(f"- {gene} {hgvsp}".rstrip())
        lines.append(f"    VAF: {vaf}")
        if hotspot:
            lines.append(f"    Hotspot: {hotspot}")
        if driver:
            lines.append(f"    Driver: {driver}")
        co_occur = v.get("co_occurring") or []
        if co_occur:
            lines.append(f"    Co-occurring: {', '.join(map(str, co_occur))}")
        lines.append("")
    return "\n".join(lines)


def _cancer_pgx_specialist(variants: list, report_data: dict) -> str:
    # Same shape as rare disease PGx sheet — chemo drugs surface naturally
    return _rd_pgx_specialist(variants, report_data)


def _cancer_clinical_evidence(variants: list, report_data: dict) -> str:
    import os

    lines = ["## Clinical Evidence Domain Data", ""]
    if variants:
        for v in variants:
            gene = v.get("gene", "?")
            hgvsp = v.get("hgvsp") or v.get("hgvsc") or ""
            civic = v.get("civic_evidence") or []
            trials = v.get("clinical_trials") or []
            lines.append(f"- {gene} {hgvsp}".rstrip())
            for ev in civic[:5]:
                drug = ev.get("drug", "")
                level = ev.get("level", "")
                direction = ev.get("direction", "")
                disease = ev.get("disease", "")
                lines.append(f"    CIViC: {drug} — {disease} (level {level}, {direction})".rstrip())
            if trials:
                lines.append(f"    Trial markers: {', '.join(map(str, trials[:5]))}")
            lines.append("")

    kb_dir = report_data.get("_kb_treatments_dir")
    if kb_dir and os.path.isdir(kb_dir):
        lines.append("## Treatment Guidelines (from KB)")
        lines.append("")
        gene_set = {v.get("gene", "").upper() for v in variants if v.get("gene")}
        try:
            md_files = sorted(f for f in os.listdir(kb_dir) if f.endswith(".md"))
        except OSError:
            md_files = []
        for fname in md_files:
            fpath = os.path.join(kb_dir, fname)
            try:
                with open(fpath, encoding="utf-8") as f:
                    content = f.read()
            except OSError:
                continue
            # Include guidelines that mention any of the gene symbols, plus the
            # first available guideline as fallback context.
            content_upper = content.upper()
            include = not gene_set or any(g in content_upper for g in gene_set)
            if include:
                lines.append(f"### {fname}")
                lines.append(content.strip())
                lines.append("")
    return "\n".join(lines)


_CANCER_BUILDERS: dict[str, Callable[[list, dict], str]] = {
    "therapeutic_target": _cancer_therapeutic_target,
    "tumor_genomics": _cancer_tumor_genomics,
    "pharmacogenomics": _cancer_pgx_specialist,
    "clinical_evidence": _cancer_clinical_evidence,
}
