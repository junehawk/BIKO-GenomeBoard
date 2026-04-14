"""Tests for scripts/tools/rerender_report.py — rebuild HTML from cached JSON."""

from __future__ import annotations

import json


def _minimal_cancer_report_json() -> dict:
    return {
        "mode": "cancer",
        "patient": {"name": "TEST-001", "age": 55, "sex": "F"},
        "sample": {"id": "S-001", "type": "Tumor"},
        "variants": [
            {
                "gene": "EGFR",
                "classification": "Pathogenic",
                "hgvsp": "p.Leu858Arg",
                "hgvsc": "c.2573T>G",
                "chrom": "chr7",
                "pos": 55259515,
                "ref": "T",
                "alt": "G",
                "tier": "I",
                "consequence": "missense_variant",
            }
        ],
        "pgx_results": [],
        "summary": {"total": 1, "pathogenic": 1},
        "clinical_board": {
            "therapeutic_implications": "EGFR L858R — EGFR TKI sensitive",
            "therapeutic_evidence": "CIViC Level A",
            "treatment_options": [{"drug": "Erlotinib", "evidence_level": "A", "resistance_notes": "T790M risk"}],
            "actionable_findings": ["Activating EGFR mutation"],
            "clinical_actions": ["Start erlotinib"],
            "immunotherapy_eligibility": "TMB-low",
            "agent_opinions": [
                {
                    "agent_name": "Variant Pathologist",
                    "domain": "variant_pathology",
                    "findings": [{"finding": "EGFR exon 21 L858R — activating"}],
                    "confidence": "high",
                }
            ],
            "agent_consensus": "unanimous",
            "confidence": "high",
            "selection_metadata": {
                "mode": "cancer",
                "total_input": 42,
                "selected": 7,
                "must_included": 3,
                "may_included": 4,
                "excluded": 35,
                "truncated": False,
                "n_dropped": 0,
                "hard_cap_applied": False,
                "empty": False,
                "empty_reason": "",
                "criteria_summary": "Tier I/II + OncoKB 1-2 + ClinVar P/LP",
                "by_selection_reason": {},
            },
        },
    }


def test_rerender_from_json_roundtrip(tmp_path):
    """Programmatic invocation: minimal JSON → rerender() → HTML file written."""
    from scripts.tools.rerender_report import rerender

    input_json = tmp_path / "report.json"
    output_html = tmp_path / "report.html"
    input_json.write_text(json.dumps(_minimal_cancer_report_json()), encoding="utf-8")

    rc = rerender(input_json, output_html)
    assert rc == 0, "rerender should succeed for a valid report JSON"
    assert output_html.exists()

    html = output_html.read_text(encoding="utf-8")
    # Report-level marker
    assert "EGFR" in html
    # Clinical board was re-rendered from the dict, so its fragment
    # should contain the treatment and the pre-analytic filter caption
    # added by Task 5 + Task 10.
    assert "Erlotinib" in html
    assert "Pre-analytic filtering" in html


def test_rerender_rejects_missing_mode(tmp_path):
    """JSON missing the 'mode' field should exit with code 2 (older format)."""
    from scripts.tools.rerender_report import rerender

    input_json = tmp_path / "bad.json"
    output_html = tmp_path / "bad.html"
    # Intentionally omit 'mode'
    input_json.write_text(json.dumps({"variants": []}), encoding="utf-8")

    rc = rerender(input_json, output_html)
    assert rc == 2
    assert not output_html.exists()


def test_rerender_reports_missing_input(tmp_path):
    """Missing input file is an I/O error (exit 1), not a format error."""
    from scripts.tools.rerender_report import rerender

    missing = tmp_path / "does-not-exist.json"
    output_html = tmp_path / "out.html"

    rc = rerender(missing, output_html)
    assert rc == 1
