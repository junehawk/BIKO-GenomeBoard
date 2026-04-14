"""Re-render an orchestrate HTML report from its cached JSON output.

Usage: python scripts/tools/rerender_report.py <input.json> <output.html>

Loads the JSON dump produced by orchestrate.py and calls the current renderer
without re-running classification, Ollama agents, or the pipeline. Useful when
a render-layer change needs to be applied to existing reports without burning
GPU time.

Exit codes:
    0 — rendered successfully
    1 — I/O error (missing input file, bad JSON, etc.)
    2 — JSON is missing required fields (older dump format)
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import fields
from pathlib import Path
from typing import Any


def _filter_dataclass_kwargs(cls, data: dict) -> dict:
    """Keep only keys that correspond to fields on the dataclass.

    Newer dumps may contain extra keys; older dumps may miss some — both
    cases are tolerated so the renderer receives a valid instance.
    """
    allowed = {f.name for f in fields(cls)}
    return {k: v for k, v in data.items() if k in allowed}


def _reconstruct_board_opinion(raw: Any):
    """Rebuild a BoardOpinion / CancerBoardOpinion from a dict.

    Returns ``None`` if ``raw`` is not a dict (e.g. older dumps where the
    dataclass was serialised via ``default=str`` as a repr string).
    """
    if not isinstance(raw, dict):
        return None

    from scripts.clinical_board.models import (
        AgentOpinion,
        BoardOpinion,
        CancerBoardOpinion,
    )

    agent_opinions = []
    for ao in raw.get("agent_opinions", []) or []:
        if isinstance(ao, dict):
            agent_opinions.append(AgentOpinion(**_filter_dataclass_kwargs(AgentOpinion, ao)))

    # Cancer opinions have therapeutic_implications; rare-disease opinions have
    # primary_diagnosis. Fall back to the presence of either key to disambiguate.
    if "therapeutic_implications" in raw or "treatment_options" in raw:
        cls = CancerBoardOpinion
    else:
        cls = BoardOpinion

    kwargs = _filter_dataclass_kwargs(cls, raw)
    kwargs["agent_opinions"] = agent_opinions
    return cls(**kwargs)


def rerender(input_json_path: Path, output_html_path: Path) -> int:
    try:
        report_data = json.loads(input_json_path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        print(f"error: input JSON not found: {input_json_path}", file=sys.stderr)
        return 1
    except json.JSONDecodeError as e:
        print(f"error: invalid JSON in {input_json_path}: {e}", file=sys.stderr)
        return 1

    if not isinstance(report_data, dict):
        print("error: JSON root must be a report_data object", file=sys.stderr)
        return 2

    mode = report_data.get("mode") or report_data.get("report_mode")
    if not mode:
        print(
            "error: JSON is missing the 'mode' field — cannot determine template (older dump format?)",
            file=sys.stderr,
        )
        return 2

    # Re-render the Clinical Board HTML fragment if a structured opinion dict
    # is available. Older dumps that serialised the dataclass as a repr string
    # fall back to whatever clinical_board_html is already in the JSON.
    from scripts.clinical_board.render import render_board_opinion_html

    board_raw = report_data.get("clinical_board")
    rebuilt = _reconstruct_board_opinion(board_raw)
    if rebuilt is not None:
        language = report_data.get("clinical_board_language", "en")
        report_data["clinical_board"] = rebuilt
        report_data["clinical_board_html"] = render_board_opinion_html(rebuilt, language=language)

    from scripts.counselor.generate_pdf import generate_report_html

    html = generate_report_html(report_data, mode=mode)
    output_html_path.parent.mkdir(parents=True, exist_ok=True)
    output_html_path.write_text(html, encoding="utf-8")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Re-render an orchestrate HTML report from cached JSON.",
    )
    parser.add_argument("input_json", type=Path, help="Cached orchestrate JSON dump")
    parser.add_argument("output_html", type=Path, help="Destination HTML file")
    args = parser.parse_args(argv)
    return rerender(args.input_json, args.output_html)


if __name__ == "__main__":
    sys.exit(main())
