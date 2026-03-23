# scripts/counselor/generate_pdf.py
from pathlib import Path
from typing import Dict
from jinja2 import Environment, FileSystemLoader

from scripts.common.gene_knowledge import get_gene_info

TEMPLATE_DIR = Path(__file__).parent.parent.parent / "templates"


def generate_report_html(report_data: Dict) -> str:
    # Enrich variants with gene knowledge
    for v in report_data.get("variants", []):
        gene = v.get("gene")
        if gene:
            info = get_gene_info(gene)
            if info:
                v.setdefault("treatment_strategies", info.get("treatment_strategies", ""))
                v.setdefault("frequency_prognosis", info.get("frequency_prognosis", ""))
                v.setdefault("finding_summary", info.get("finding_summary", ""))
                v.setdefault("gene_full_name", info.get("full_name", ""))
                v.setdefault("associated_conditions", info.get("associated_conditions", []))
                v.setdefault("korean_specific_note", info.get("korean_specific_note"))

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
