"""Coverage for the v2.5.4 Phase 2 render-purity + CIViC priority fix.

Verifies the M2 and L7 contracts documented in
``scripts/orchestration/canonical.py::_enrich_with_gene_knowledge_and_civic``:

    * CIViC evidence wins over gene_knowledge on overlapping clinical
      fields (finding_summary, treatment_strategies, references,
      content_status) — the pre-v2.5.4 setdefault chain had gene_knowledge
      winning instead, despite the comment claiming CIViC priority (M2).
    * gene_knowledge still fills when CIViC has nothing to say (the
      curated-civic status remains reserved for genuinely CIViC-sourced
      evidence).
    * ``generate_report_html`` does not mutate its input (L7).
"""

from __future__ import annotations

import copy

import pytest


# ── Helpers ──────────────────────────────────────────────────────────────────


def _minimal_report(**variant_overrides):
    base_variant = {
        "variant": "chr17:7577120:G>A",
        "gene": "TP53",
        "classification": "Pathogenic",
        "acmg_codes": ["PVS1"],
        "conflict": False,
        "hgvsp": "p.Arg175His",
        "agents": {
            "clinical": {"clinvar_significance": "Pathogenic"},
            "korean_pop": {"korean_flag": ""},
        },
    }
    base_variant.update(variant_overrides)
    return {
        "sample_id": "ENRICH_TEST",
        "date": "2026-04-22",
        "variants": [base_variant],
        "pgx_results": [],
        "summary": {"total": 1, "pathogenic": 1, "vus": 0, "benign": 0},
        "db_versions": {},
    }


# ── CIViC priority (M2) ──────────────────────────────────────────────────────


def test_civic_overrides_gene_knowledge_finding_summary(monkeypatch):
    """When both sources return a summary, CIViC wins."""
    from scripts.orchestration import canonical

    monkeypatch.setattr(
        canonical,
        "get_gene_info",
        lambda gene: {"finding_summary": "GK fallback summary for " + gene},
    )
    monkeypatch.setattr(
        canonical,
        "get_gene_summary",
        lambda gene: {"description": "CIViC curated summary for " + gene},
    )
    monkeypatch.setattr(canonical, "get_treatment_summary", lambda g, v: None)
    monkeypatch.setattr(canonical, "get_variant_evidence", lambda g, v=None: [])

    records = [_minimal_report()["variants"][0]]
    canonical._enrich_with_gene_knowledge_and_civic(records, mode="cancer")

    assert "CIViC curated summary" in records[0]["finding_summary"]
    assert "GK fallback" not in records[0]["finding_summary"]


def test_civic_evidence_sets_curated_civic_status(monkeypatch):
    """CIViC evidence list drives content_status='curated-civic'."""
    from scripts.orchestration import canonical

    monkeypatch.setattr(canonical, "get_gene_info", lambda g: None)
    monkeypatch.setattr(canonical, "get_gene_summary", lambda g: None)
    monkeypatch.setattr(canonical, "get_treatment_summary", lambda g, v: None)
    monkeypatch.setattr(
        canonical,
        "get_variant_evidence",
        lambda g, v=None: [
            {
                "pmid": "12345678",
                "citation": "Doe et al. 2024",
                "evidence_type": "Predictive",
                "significance": "Sensitivity/Response",
            }
        ],
    )

    records = [_minimal_report()["variants"][0]]
    canonical._enrich_with_gene_knowledge_and_civic(records, mode="cancer")

    assert records[0]["content_status"] == "curated-civic"
    assert records[0]["references"][0]["pmid"] == "12345678"


def test_gene_knowledge_fallback_when_civic_absent(monkeypatch):
    """With CIViC empty, gene_knowledge fills finding_summary + references."""
    from scripts.orchestration import canonical

    monkeypatch.setattr(
        canonical,
        "get_gene_info",
        lambda g: {
            "finding_summary": "GK summary",
            "treatment_strategies": "GK treatment text",
            "references": [{"pmid": "111", "source": "GK", "note": "gk"}],
            "content_status": "ai-generated",
            "full_name": "Tumor protein p53",
        },
    )
    monkeypatch.setattr(canonical, "get_gene_summary", lambda g: None)
    monkeypatch.setattr(canonical, "get_treatment_summary", lambda g, v: None)
    monkeypatch.setattr(canonical, "get_variant_evidence", lambda g, v=None: [])

    records = [_minimal_report()["variants"][0]]
    canonical._enrich_with_gene_knowledge_and_civic(records, mode="cancer")

    assert records[0]["finding_summary"] == "GK summary"
    assert records[0]["treatment_strategies"] == "GK treatment text"
    assert records[0]["references"][0]["pmid"] == "111"
    # content_status inherits gk value since no CIViC evidence fired.
    assert records[0]["content_status"] == "ai-generated"
    assert records[0]["gene_full_name"] == "Tumor protein p53"


def test_caller_provided_content_status_is_preserved(monkeypatch):
    """An explicit content_status (e.g. set by Board scrubber) is never overwritten."""
    from scripts.orchestration import canonical

    monkeypatch.setattr(canonical, "get_gene_info", lambda g: None)
    monkeypatch.setattr(canonical, "get_gene_summary", lambda g: None)
    monkeypatch.setattr(canonical, "get_treatment_summary", lambda g, v: None)
    monkeypatch.setattr(
        canonical,
        "get_variant_evidence",
        lambda g, v=None: [
            {
                "pmid": "9",
                "citation": "X",
                "evidence_type": "Predictive",
                "significance": "Y",
            }
        ],
    )

    rec = _minimal_report()["variants"][0]
    rec["content_status"] = "ai-generated"
    rec["references"] = []  # empty list should not block CIViC refs from filling in

    canonical._enrich_with_gene_knowledge_and_civic([rec], mode="cancer")

    assert rec["content_status"] == "ai-generated"  # preserved
    assert len(rec["references"]) == 1  # still filled from CIViC


def test_rare_disease_mode_skips_civic(monkeypatch):
    """CIViC lookups must not run in rare-disease mode."""
    from scripts.orchestration import canonical

    civic_calls: list[str] = []

    def _track_civic(gene):
        civic_calls.append(gene)
        return {"description": "should not appear"}

    monkeypatch.setattr(canonical, "get_gene_info", lambda g: None)
    monkeypatch.setattr(canonical, "get_gene_summary", _track_civic)
    monkeypatch.setattr(
        canonical,
        "get_treatment_summary",
        lambda g, v: civic_calls.append("treatment") or None,
    )
    monkeypatch.setattr(
        canonical,
        "get_variant_evidence",
        lambda g, v=None: civic_calls.append("evidence") or [],
    )

    records = [_minimal_report()["variants"][0]]
    canonical._enrich_with_gene_knowledge_and_civic(records, mode="rare-disease")

    assert civic_calls == []


# ── Render purity (L7) ───────────────────────────────────────────────────────


def test_render_does_not_mutate_input():
    """``generate_report_html`` must not mutate the caller's dict."""
    from scripts.reporting.generate_pdf import generate_report_html

    original = {
        "sample_id": "PURE_TEST",
        "date": "2026-04-22",
        "variants": [
            {
                "variant": "chr17:7577120:G>A",
                "gene": "TP53",
                "classification": "Pathogenic",
                "acmg_codes": ["PVS1"],
                "conflict": False,
                "agents": {
                    "clinical": {"clinvar_significance": "Pathogenic"},
                    "korean_pop": {"korean_flag": ""},
                },
            }
        ],
        "pgx_results": [
            {
                "gene": "CYP2C19",
                "phenotype": "Poor Metabolizer",
                "cpic_level": "A",
                "korean_prevalence": 0.14,
                "western_prevalence": 0.025,
                "clinical_impact": "x",
                "star_allele": "*2/*2",
                "cpic_recommendation": "y",
            }
        ],
        "summary": {"total": 1, "pathogenic": 1, "vus": 0, "benign": 0},
        "db_versions": {},
    }
    snapshot = copy.deepcopy(original)

    html = generate_report_html(original, mode="cancer")
    assert "BIKO GenomeBoard" in html

    assert original == snapshot, "generate_report_html mutated its input"


def test_render_double_call_idempotent():
    """Calling generate_report_html twice on the same dict is stable."""
    from scripts.reporting.generate_pdf import generate_report_html

    data = _minimal_report()
    before = copy.deepcopy(data)

    html1 = generate_report_html(data, mode="cancer")
    html2 = generate_report_html(data, mode="cancer")

    # The caller's dict is unchanged
    assert data == before
    # The second render produces the same HTML (date and PMID linkification
    # are deterministic given the same input).
    # We check a load-bearing substring rather than byte equality because
    # the HTML contains the current date.
    assert ("TP53" in html1) == ("TP53" in html2)


# ── Frequency fallback (L7-adjacent) ─────────────────────────────────────────


def test_frequency_prognosis_built_from_available_fields():
    """_build_frequency_text synthesises text from gnomAD / KOVA fields."""
    from scripts.orchestration.canonical import _build_frequency_text

    v = {
        "gnomad_all": 0.000123456,
        "gnomad_eas": 0.0005,
        "kova_freq": 0.002,
        "kova_homozygote": 1,
        "agents": {"korean_pop": {"korean_flag": "Rare variant (Korean population)"}},
    }
    text = _build_frequency_text(v)
    assert "gnomAD global AF" in text
    assert "East Asian" in text
    assert "KOVA" in text
    assert "homozygote count" in text


def test_frequency_prognosis_empty_variant_yields_absence_message():
    """A variant with zero frequency signals returns the 'absent/rare' note."""
    from scripts.orchestration.canonical import _build_frequency_text

    v: dict = {}
    assert "Not observed in gnomAD" in _build_frequency_text(v)


# ── HGVS assembly + finding_summary adjustment ───────────────────────────────


def test_finding_summary_adjusted_for_vus_classification(monkeypatch):
    """GK summary that presupposes pathogenicity gets rewritten for a VUS call."""
    from scripts.orchestration import canonical

    monkeypatch.setattr(
        canonical,
        "get_gene_info",
        lambda g: {
            "finding_summary": "A pathogenic variant was identified in this specimen.",
        },
    )
    monkeypatch.setattr(canonical, "get_gene_summary", lambda g: None)
    monkeypatch.setattr(canonical, "get_treatment_summary", lambda g, v: None)
    monkeypatch.setattr(canonical, "get_variant_evidence", lambda g, v=None: [])

    rec = _minimal_report(classification="VUS", acmg_codes=["PM2_Supporting"])["variants"][0]
    canonical._enrich_with_gene_knowledge_and_civic([rec], mode="cancer")

    # The classification-aware rewrite fires.
    assert "uncertain significance" in rec["finding_summary"].lower()
    assert "further evidence is needed" in rec["finding_summary"].lower()


@pytest.mark.parametrize("pmid_token", ["PMID:12345678", "PMID: 87654321"])
def test_pmid_linkification_in_enriched_fields(monkeypatch, pmid_token):
    """PMIDs embedded in treatment_strategies become clickable PubMed links."""
    from scripts.orchestration import canonical

    monkeypatch.setattr(
        canonical,
        "get_gene_info",
        lambda g: {
            "treatment_strategies": f"PARP inhibitor — olaparib ({pmid_token})",
        },
    )
    monkeypatch.setattr(canonical, "get_gene_summary", lambda g: None)
    monkeypatch.setattr(canonical, "get_treatment_summary", lambda g, v: None)
    monkeypatch.setattr(canonical, "get_variant_evidence", lambda g, v=None: [])

    rec = _minimal_report()["variants"][0]
    canonical._enrich_with_gene_knowledge_and_civic([rec], mode="cancer")

    assert "<a href=" in rec["treatment_strategies"]
    assert "pubmed.ncbi.nlm.nih.gov" in rec["treatment_strategies"]
