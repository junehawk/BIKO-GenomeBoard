"""HTML + PDF rendering.

v2.5.4 Phase 2: this module is a **pure renderer**. The variant and PGx
enrichment logic that previously lived here (gene_knowledge fallback,
CIViC override, PMID linkification, classification-aware finding_summary
rewrite, frequency text synthesis) has moved to
``scripts.orchestration.canonical``. ``generate_report_html`` now:

    * Does **not** mutate its input dict. A shallow copy of the
      ``variants`` + ``pgx_results`` lists (with deep-copied items) is
      used for any defensive enrichment.
    * Runs the canonical enrichment helpers only when the input looks
      un-enriched (legacy test fixtures, rerender_report callers that
      loaded a pre-Phase-2 JSON dump). When the input is already enriched
      (the default path from ``build_sample_report``) the helpers are
      idempotent — every field guard uses ``not v.get(...)`` rather than
      ``setdefault`` — so calling them twice is safe.
    * Never changes which data source wins for overlapping fields. The
      CIViC-over-gene_knowledge priority (M2 fix) is enforced inside the
      canonical enrichment helper.
"""

from __future__ import annotations

import copy
import os
from pathlib import Path
from typing import Any

from jinja2 import Environment, FileSystemLoader

from scripts.common.config import get
from scripts.orchestration.canonical import (
    _enrich_pgx_with_gene_knowledge,
    _enrich_with_gene_knowledge_and_civic,
)

# Re-export for tests that monkey-patch the old module paths.
from scripts.common.gene_knowledge import get_gene_info  # noqa: F401
from scripts.common.hgvs_utils import hgvsp_to_civic_variant as _hgvsp_to_civic_variant  # noqa: F401
from scripts.storage.query_civic import (  # noqa: F401
    get_gene_summary,
    get_treatment_summary,
    get_variant_evidence,
)


def _maybe_enrich_view_model(report_data: dict, mode: str) -> dict:
    """Build a render-ready view model without mutating ``report_data``.

    Deep-copies the ``variants`` and ``pgx_results`` lists so the caller's
    dict and its variant/pgx dicts stay untouched (L7 fix). Then runs the
    canonical enrichment helpers on the copies — idempotent if already
    enriched by ``build_sample_report``.
    """
    view: dict[str, Any] = dict(report_data)
    view["variants"] = copy.deepcopy(report_data.get("variants", []))
    view["pgx_results"] = copy.deepcopy(report_data.get("pgx_results", []))
    _enrich_with_gene_knowledge_and_civic(view["variants"], mode)
    _enrich_pgx_with_gene_knowledge(view["pgx_results"])
    return view


def generate_report_html(report_data: dict, mode: str = "cancer") -> str:
    """Render ``report_data`` to an HTML document.

    The input is not mutated. All variant / PGx enrichment is performed
    on a defensive copy returned to the Jinja environment.
    """
    # Determine template directory based on mode
    templates_base = get("paths.templates") or str(Path(__file__).parent.parent.parent / "templates")
    template_dir = os.path.join(templates_base, mode)
    if not os.path.exists(os.path.join(template_dir, "report.html")):
        template_dir = os.path.join(templates_base, "cancer")

    view = _maybe_enrich_view_model(report_data, mode)

    # SV / TMB default keys (render-time defaults, not mutation of caller).
    view.setdefault("sv_class45", [])
    view.setdefault("sv_class3_display", [])
    view.setdefault("sv_class3_hidden", 0)
    view.setdefault("sv_benign_count", 0)
    view.setdefault("sv_variants", [])
    view.setdefault("tmb", None)

    # Add shared templates directory to loader
    shared_dir = os.path.join(templates_base, "shared")
    loader = FileSystemLoader([template_dir, templates_base, shared_dir])

    env = Environment(loader=loader, autoescape=True)
    template = env.get_template("report.html")
    return template.render(**view)


def generate_pdf(report_data: dict, output_path: str, mode: str = "cancer") -> str:
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
