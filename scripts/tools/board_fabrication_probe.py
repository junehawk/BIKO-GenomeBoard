#!/usr/bin/env python3
"""Fabrication probe for Clinical Board reports (v2.8 Phase 3).

Turns the qualitative "did the board hallucinate?" question into a repeatable
NUMBER: the count of known gene symbols named in the board narrative that are NOT
in the case variant set. Run it on an ``orchestrate --json`` report — or several
at once for a before/after or deliberation-on/off A/B — to measure fabrication
objectively instead of eyeballing headlines.

Usage::

    python scripts/tools/board_fabrication_probe.py report1.json [report2.json ...]

The probe is Ollama-free (it analyses an already-produced report), so it is fast,
deterministic, and unit-testable.
"""

from __future__ import annotations

import json
import os
import sys
from typing import Iterable, Set

# Allow running as a bare script (python scripts/tools/board_fabrication_probe.py)
# by putting the repo root on sys.path before importing the package.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from scripts.clinical_board.grounding_scrubber import find_offbriefing_genes  # noqa: E402

_TEXT_KEYS = (
    "therapeutic_headline",
    "therapeutic_implications",
    "therapeutic_evidence",
    "immunotherapy_eligibility",
    "primary_diagnosis",
    "primary_diagnosis_evidence",
)
_LIST_KEYS = (
    "actionable_findings",
    "clinical_actions",
    "monitoring_plan",
    "key_findings",
    "recommendations",
    "follow_up",
    "dissenting_opinions",
)


def _board_texts(cb: dict) -> Iterable[str]:
    for k in _TEXT_KEYS:
        v = cb.get(k)
        if isinstance(v, str) and v:
            yield v
    for k in _LIST_KEYS:
        for item in cb.get(k) or []:
            if isinstance(item, str) and item:
                yield item
    for dx in cb.get("differential_diagnoses") or []:
        if isinstance(dx, dict):
            for kk in ("diagnosis", "evidence"):
                if isinstance(dx.get(kk), str) and dx[kk]:
                    yield dx[kk]


def _briefing_genes(report: dict) -> Set[str]:
    out: Set[str] = set()
    for key in ("_board_variants", "variants"):
        for v in report.get(key) or []:
            gene = (v.get("gene") if isinstance(v, dict) else None) or ""
            if gene:
                out.add(str(gene).upper())
    return out


def probe(report: dict) -> dict:
    """Return a fabrication metric for one report_data dict."""
    cb = report.get("clinical_board") or {}
    briefing = _briefing_genes(report)
    off = find_offbriefing_genes(_board_texts(cb), briefing)
    return {
        "sample_id": report.get("sample_id"),
        "mode": report.get("mode"),
        "case_genes": sorted(briefing),
        "offbriefing_genes": off,
        "fabrication_count": len(off),
        "headline": cb.get("therapeutic_headline") or cb.get("primary_diagnosis") or "",
        "grounding_flags": cb.get("grounding_flags") or [],
        "consensus": cb.get("agent_consensus"),
        "confidence": cb.get("confidence"),
    }


def main(argv: list) -> int:
    if len(argv) < 2:
        print("usage: board_fabrication_probe.py REPORT.json [REPORT2.json ...]")
        return 2
    for path in argv[1:]:
        try:
            with open(path) as f:
                report = json.load(f)
        except (OSError, json.JSONDecodeError) as exc:
            print(f"\n=== {path} ===\n  ERROR: {exc}")
            continue
        r = probe(report)
        print(f"\n=== {path} ===")
        print(f"  sample={r['sample_id']} mode={r['mode']} consensus={r['consensus']} confidence={r['confidence']}")
        print(f"  case genes: {', '.join(r['case_genes']) or '(none)'}")
        print(f"  headline: {r['headline'][:120]}")
        print(f"  FABRICATION COUNT (off-briefing genes): {r['fabrication_count']}  {r['offbriefing_genes']}")
        if r["grounding_flags"]:
            print(f"  grounding_flags: {r['grounding_flags']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
