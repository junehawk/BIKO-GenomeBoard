# scripts/counselor/generate_pdf.py
from pathlib import Path
from typing import Dict
from jinja2 import Environment, FileSystemLoader

TEMPLATE_DIR = Path(__file__).parent.parent.parent / "templates"


def generate_report_html(report_data: Dict) -> str:
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
