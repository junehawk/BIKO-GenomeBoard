"""Build per-agent, per-mode domain sheets from local DB data.

Each builder consumes the variants/report_data already collected by the
pipeline and formats a focused, human-readable text section for one agent.
No new DB queries are issued here — the pipeline has already collected the
relevant fields.
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
        lines.append(
            f"    In silico: REVEL={revel}, CADD={cadd}, SpliceAI={spliceai}"
        )
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
            lines.append(f"    OMIM phenotypes: {'; '.join(omim_phenotypes)}")
        if matching_hpo:
            hpo_strs = [
                f"{h.get('id', '?')} {h.get('name', '')}".strip()
                for h in matching_hpo
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
                lines.append(
                    f"    CIViC: {drug} (level {level}, {direction})".rstrip()
                )
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
# Cancer builders (populated in Task 4)
# ---------------------------------------------------------------------------

_CANCER_BUILDERS: dict[str, Callable[[list, dict], str]] = {}
