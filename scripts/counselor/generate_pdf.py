# scripts/counselor/generate_pdf.py
from pathlib import Path
from typing import Dict
from jinja2 import Environment, FileSystemLoader

from scripts.common.gene_knowledge import get_gene_info
from scripts.common.config import get

TEMPLATE_DIR = Path(get("paths.templates") or str(Path(__file__).parent.parent.parent / "templates"))


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


def generate_report_html(report_data: Dict) -> str:
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
                # VCF annotation takes priority over static gene_knowledge for hgvs/variant_effect
                if v.get("hgvsc") or v.get("hgvsp"):
                    v.setdefault("hgvs", {
                        "transcript": v.get("transcript", ""),
                        "cdna": v.get("hgvsc", ""),
                        "protein": v.get("hgvsp", ""),
                        "variant_effect": v.get("consequence", ""),
                    })
                    v.setdefault("variant_effect", v.get("consequence", ""))
                else:
                    v.setdefault("hgvs", info.get("hgvs", {}))
                    hgvs = info.get("hgvs", {})
                    v.setdefault("variant_effect", hgvs.get("variant_effect", ""))

        # Always apply classification-aware adjustment to finding_summary
        classification = v.get("classification", "VUS")
        raw_summary = v.get("finding_summary", "")
        if raw_summary:
            v["finding_summary"] = _adjust_finding_summary(raw_summary, classification)

    for pgx in report_data.get("pgx_results", []):
        gene = pgx.get("gene")
        if gene:
            info = get_gene_info(gene)
            if info:
                pgx.setdefault("hgvs", info.get("hgvs", {}))
                pgx.setdefault("variant_effect", info.get("hgvs", {}).get("variant_effect", ""))

    env = Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)), autoescape=True)
    template = env.get_template("report.html")
    return template.render(**report_data)


def generate_pdf(report_data: Dict, output_path: str) -> str:
    html = generate_report_html(report_data)
    try:
        from weasyprint import HTML
        HTML(string=html).write_pdf(output_path)
    except ImportError:
        import logging
        logging.getLogger(__name__).warning("WeasyPrint not available. Saving HTML report instead.")
        html_path = output_path.replace('.pdf', '.html')
        Path(html_path).write_text(html)
        return html_path
    return output_path
