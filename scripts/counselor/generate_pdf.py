# scripts/counselor/generate_pdf.py
import os
import re
from pathlib import Path
from typing import Dict, Optional
from jinja2 import Environment, FileSystemLoader

from scripts.common.gene_knowledge import get_gene_info
from scripts.common.config import get
from scripts.db.query_civic import get_gene_summary, get_treatment_summary, get_variant_evidence
from scripts.common.hgvs_utils import hgvsp_to_civic_variant as _hgvsp_to_civic_variant


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


def _linkify_pmids(text):
    """Convert PMID references in text to PubMed links."""
    if not text:
        return text
    def _replace_pmid(match):
        pmid = match.group(1)
        return f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank" rel="noopener" style="color:#1d4ed8;text-decoration:none;">PMID:{pmid}</a>'
    return re.sub(r'PMID:\s*(\d+)', _replace_pmid, text)


def generate_report_html(report_data: Dict, mode: str = "cancer") -> str:
    # Determine template directory based on mode
    templates_base = get("paths.templates") or str(Path(__file__).parent.parent.parent / "templates")
    template_dir = os.path.join(templates_base, mode)

    # Fall back to cancer if mode-specific template doesn't exist
    if not os.path.exists(os.path.join(template_dir, "report.html")):
        template_dir = os.path.join(templates_base, "cancer")

    # Enrich variants with gene knowledge
    for v in report_data.get("variants", []):
        gene = v.get("gene")
        if gene:
            info = get_gene_info(gene)
            if info:
                v.setdefault("treatment_strategies", info.get("treatment_strategies", ""))
                v.setdefault("frequency_prognosis", info.get("frequency_prognosis", ""))
                v.setdefault("gene_full_name", info.get("full_name", ""))
                v.setdefault("associated_conditions", info.get("associated_conditions", []))
                v.setdefault("korean_specific_note", info.get("korean_specific_note"))
                v.setdefault("finding_summary", info.get("finding_summary", ""))
                v.setdefault("references", info.get("references", []))
                v.setdefault("content_status", info.get("content_status", "ai-generated"))
                # VCF annotation takes priority over static gene_knowledge for hgvs/variant_effect
                if v.get("hgvsc") or v.get("hgvsp"):
                    v.setdefault(
                        "hgvs",
                        {
                            "transcript": v.get("transcript", ""),
                            "cdna": v.get("hgvsc", ""),
                            "protein": v.get("hgvsp", ""),
                            "variant_effect": v.get("consequence", ""),
                        },
                    )
                    v.setdefault("variant_effect", v.get("consequence", ""))
                else:
                    v.setdefault("hgvs", info.get("hgvs", {}))
                    hgvs = info.get("hgvs", {})
                    v.setdefault("variant_effect", hgvs.get("variant_effect", ""))

            # CIViC enrichment (takes priority over gene_knowledge for clinical content)
            civic_gene = get_gene_summary(gene)
            if civic_gene and civic_gene.get("description"):
                v.setdefault("finding_summary", civic_gene["description"][:500])

            # Treatment from CIViC
            hgvsp = v.get("hgvsp", "")
            civic_variant_name = _hgvsp_to_civic_variant(hgvsp)
            civic_treatment = get_treatment_summary(gene, civic_variant_name)
            if civic_treatment:
                v.setdefault("treatment_strategies", civic_treatment)

            # Evidence references from CIViC
            civic_evidence = get_variant_evidence(gene, civic_variant_name) if civic_variant_name else []
            if not civic_evidence:
                civic_evidence = get_variant_evidence(gene)
            if civic_evidence:
                refs = [
                    {
                        "pmid": e["pmid"],
                        "source": e["citation"],
                        "note": f"{e['evidence_type']} — {e['significance']}",
                    }
                    for e in civic_evidence[:5]
                    if e["pmid"]
                ]
                if refs:
                    v.setdefault("references", refs)
                    v.setdefault("content_status", "curated-civic")

        # Build frequency text from available data if frequency_prognosis is empty
        if not v.get("frequency_prognosis"):
            freq_parts = []
            if v.get("gnomad_all") is not None:
                freq_parts.append(f"gnomAD global AF: {v['gnomad_all']:.6f}")
            if v.get("gnomad_eas") is not None:
                freq_parts.append(f"gnomAD East Asian AF: {v['gnomad_eas']:.6f}")
            kf = ""
            if v.get("agents") and v["agents"].get("korean_pop"):
                kf = v["agents"]["korean_pop"].get("korean_flag", "")
            if kf and kf not in ("No notable findings", "No frequency data available"):
                freq_parts.append(f"Korean: {kf}")
            if v.get("krgdb_freq") is not None:
                freq_parts.append(f"KRGDB Korean AF: {v['krgdb_freq']}")
            if freq_parts:
                v["frequency_prognosis"] = ". ".join(freq_parts) + "."
            elif v.get("gnomad_all") is None:
                v["frequency_prognosis"] = "Not observed in gnomAD population database (absent or extremely rare)."

        # Always apply classification-aware adjustment to finding_summary
        classification = v.get("classification", "VUS")
        raw_summary = v.get("finding_summary", "")
        if raw_summary:
            v["finding_summary"] = _adjust_finding_summary(raw_summary, classification)

        # Convert PMID references to clickable HTML links (safe — we control this text)
        if v.get("treatment_strategies"):
            v["treatment_strategies"] = _linkify_pmids(v["treatment_strategies"])
        if v.get("finding_summary"):
            v["finding_summary"] = _linkify_pmids(v["finding_summary"])

    for pgx in report_data.get("pgx_results", []):
        gene = pgx.get("gene")
        if gene:
            info = get_gene_info(gene)
            if info:
                pgx.setdefault("hgvs", info.get("hgvs", {}))
                pgx.setdefault("variant_effect", info.get("hgvs", {}).get("variant_effect", ""))

    env = Environment(loader=FileSystemLoader(template_dir), autoescape=True)
    template = env.get_template("report.html")
    return template.render(**report_data)


def generate_pdf(report_data: Dict, output_path: str, mode: str = "cancer") -> str:
    html = generate_report_html(report_data, mode=mode)
    try:
        from weasyprint import HTML

        HTML(string=html).write_pdf(output_path)
    except ImportError:
        import logging

        logging.getLogger(__name__).warning("WeasyPrint not available. Saving HTML report instead.")
        html_path = output_path.replace(".pdf", ".html")
        Path(html_path).write_text(html)
        return html_path
    return output_path
