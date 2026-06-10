"""Grounding scrubber (v2.8 Phase 1) — annotate off-briefing gene mentions.

Pins the deterministic enforcement layer that backs up the prompt grounding
clause: known gene symbols named in the board narrative but absent from the case
variant set are flagged (annotate-only) so the reviewer can verify. Motivated by
the empirical KRAS-G12D fabrication on a KRAS-absent case.
"""

from __future__ import annotations

from scripts.clinical_board.grounding_scrubber import annotate_grounding
from scripts.clinical_board.models import BoardOpinion, CancerBoardOpinion


def _report(*genes: str) -> dict:
    return {"variants": [{"gene": g} for g in genes]}


def test_flags_offbriefing_gene_in_cancer_headline():
    """KRAS named in the headline but not in the case → flagged."""
    op = CancerBoardOpinion(therapeutic_headline="KRAS G12D driver, no standard therapy")
    stats = annotate_grounding(op, _report("TP53", "BRCA2"))
    assert "KRAS" in stats["offbriefing_genes"]
    assert op.grounding_flags  # a note was appended
    assert "KRAS" in op.grounding_flags[0]


def test_does_not_flag_briefing_genes():
    """Genes actually in the case must NOT be flagged."""
    op = CancerBoardOpinion(
        therapeutic_headline="TP53/BRCA2 biallelic inactivation, HRD candidate",
        therapeutic_implications="BRCA2 frameshift supports PARP inhibition.",
    )
    stats = annotate_grounding(op, _report("TP53", "BRCA2"))
    assert stats["offbriefing_genes"] == []
    assert op.grounding_flags == []


def test_ignores_clinical_abbreviations():
    """TMB / HRD / MSI are neither genes nor tumour types and must not be flagged.

    (NSCLC is a tumour type — covered by the tumour-type tests below, not here.)
    """
    op = CancerBoardOpinion(
        therapeutic_headline="TMB-high, MSI-stable, HRD-positive",
        immunotherapy_eligibility="TMB-high → pembrolizumab candidate.",
    )
    stats = annotate_grounding(op, _report("TP53"))
    assert stats["offbriefing_genes"] == []
    assert op.grounding_flags == []


def test_flags_offbriefing_gene_in_rare_disease_opinion():
    op = BoardOpinion(
        primary_diagnosis="Noonan syndrome",
        primary_diagnosis_evidence="PTPN11 missense consistent with RASopathy.",
        key_findings=["SCN1A variant noted"],
    )
    # Case actually contains only SCN1A → PTPN11 is off-briefing.
    stats = annotate_grounding(op, _report("SCN1A"))
    assert "PTPN11" in stats["offbriefing_genes"]
    assert "SCN1A" not in stats["offbriefing_genes"]


def test_skips_treatment_option_drug_field():
    """treatment_options drug names are governed by narrative_scrubber, not here —
    a gene-like token there should not be scanned by the grounding scrubber."""
    op = CancerBoardOpinion(
        therapeutic_headline="TP53 LOF",
        treatment_options=[{"drug": "KRAS-pathway inhibitor", "curated_id": "x", "evidence_level": "B"}],
    )
    stats = annotate_grounding(op, _report("TP53"))
    assert stats["offbriefing_genes"] == []  # KRAS in a drug field is not scanned


def test_annotate_grounding_never_raises_without_grounding_flags_attr():
    """Defensive: an object lacking grounding_flags returns cleanly."""

    class _Bare:
        therapeutic_headline = "KRAS driver"

    stats = annotate_grounding(_Bare(), _report("TP53"))
    assert stats == {"offbriefing_genes": [], "offbriefing_tumor_terms": []}


# ── Phase 3.2: fabrication probe (quantitative metric) ────────────────────────


def test_fabrication_probe_counts_offbriefing_genes():
    from scripts.tools.board_fabrication_probe import probe

    report = {
        "sample_id": "PROBE-TEST",
        "mode": "cancer",
        "variants": [{"gene": "TP53"}, {"gene": "BRCA2"}],
        "clinical_board": {
            "therapeutic_headline": "KRAS G12D driver; TP53 LOF, BRCA2 frameshift",
            "agent_consensus": "majority",
            "confidence": "moderate",
        },
    }
    r = probe(report)
    assert r["fabrication_count"] == 1
    assert r["offbriefing_genes"] == ["KRAS"]  # TP53/BRCA2 are in the case
    assert set(r["case_genes"]) == {"TP53", "BRCA2"}


def test_fabrication_probe_zero_when_grounded():
    from scripts.tools.board_fabrication_probe import probe

    report = {
        "sample_id": "PROBE-CLEAN",
        "mode": "cancer",
        "variants": [{"gene": "TP53"}, {"gene": "BRCA2"}],
        "clinical_board": {"therapeutic_headline": "TP53/BRCA2 co-inactivation, HRD candidate"},
    }
    assert probe(report)["fabrication_count"] == 0


# ── Phase 4: tumour-type / stage fabrication (the v2.8 residual) ──────────────


def test_flags_tumor_type_when_no_clinical_note():
    """'Stage III NSCLC' with no clinical note → flagged (ungrounded tumour type)."""
    op = CancerBoardOpinion(therapeutic_headline="Stage III NSCLC — TP53/BRCA2 drivers, low TMB")
    stats = annotate_grounding(op, _report("TP53", "BRCA2"))  # no clinical_note key
    assert stats["offbriefing_tumor_terms"]  # non-empty
    joined = " ".join(stats["offbriefing_tumor_terms"]).lower()
    assert "nsclc" in joined and "stage" in joined
    assert any("tumour type" in f.lower() for f in op.grounding_flags)


def test_no_tumor_flag_when_clinical_note_present():
    """A clinical note grounds the tumour type → no flag."""
    op = CancerBoardOpinion(therapeutic_headline="Stage III NSCLC — TP53/BRCA2 drivers")
    report = {"variants": [{"gene": "TP53"}], "clinical_note": "Stage III NSCLC, prior platinum."}
    stats = annotate_grounding(op, report)
    assert stats["offbriefing_tumor_terms"] == []


def test_tumor_term_detector_matches_histology_and_stage():
    from scripts.clinical_board.grounding_scrubber import find_offbriefing_tumor_terms

    terms = find_offbriefing_tumor_terms(
        ["pancreatic adenocarcinoma, Stage IV", "glioblastoma", "no tumour here"],
        has_clinical_note=False,
    )
    low = " ".join(terms).lower()
    assert "adenocarcinoma" in low
    assert "stage iv" in low
    assert "glioblastoma" in low or "glioma" in low
    # When a note is present, nothing is flagged.
    assert find_offbriefing_tumor_terms(["Stage IV adenocarcinoma"], has_clinical_note=True) == []


def test_probe_counts_tumor_fabrication():
    from scripts.tools.board_fabrication_probe import probe

    report = {
        "sample_id": "T",
        "mode": "cancer",
        "variants": [{"gene": "TP53"}, {"gene": "BRCA2"}],
        "clinical_board": {"therapeutic_headline": "Stage III NSCLC — TP53/BRCA2 drivers"},
    }
    r = probe(report)
    assert r["tumor_fabrication_count"] >= 1
    assert r["gene_fabrication_count"] == 0  # TP53/BRCA2 are in the case
    assert r["fabrication_count"] == r["gene_fabrication_count"] + r["tumor_fabrication_count"]


# ── v2.8 (A): grounding flags rendered in the report ──────────────────────────


def test_report_renders_grounding_flags_banner_cancer():
    from scripts.clinical_board.render import render_board_opinion_html

    op = CancerBoardOpinion(
        therapeutic_headline="Stage III NSCLC — TP53/BRCA2",
        grounding_flags=["⚠ Grounding: gene KRAS not in case variant set"],
    )
    out = render_board_opinion_html(op, "en")
    assert "Grounding Warnings" in out
    assert "KRAS not in case variant set" in out


def test_report_renders_grounding_flags_banner_rare():
    from scripts.clinical_board.render import render_board_opinion_html

    op = BoardOpinion(primary_diagnosis="X", grounding_flags=["⚠ Grounding: PTPN11 not in case"])
    out = render_board_opinion_html(op, "en")
    assert "Grounding Warnings" in out
    assert "PTPN11" in out


def test_report_no_banner_without_flags():
    from scripts.clinical_board.render import render_board_opinion_html

    out = render_board_opinion_html(CancerBoardOpinion(therapeutic_headline="X", grounding_flags=[]), "en")
    assert "Grounding Warnings" not in out


def test_grounding_banner_escapes_html():
    from scripts.clinical_board.render import render_board_opinion_html

    op = CancerBoardOpinion(grounding_flags=["<script>alert(1)</script>"])
    out = render_board_opinion_html(op, "en")
    assert "<script>alert(1)</script>" not in out
    assert "&lt;script&gt;" in out


def test_flags_offbriefing_gene_in_agent_opinion():
    """A domain agent that names an off-case gene in its finding must be flagged too —
    agent opinions are rendered verbatim, so the grounding scrubber must scan them, not
    just the Chair narrative (BOARD-03 / T5-04)."""
    from scripts.clinical_board.models import AgentOpinion

    op = CancerBoardOpinion(
        therapeutic_headline="TP53-driven tumour",  # Chair narrative is clean
        agent_opinions=[
            AgentOpinion(
                agent_name="Variant Pathologist",
                domain="variant_pathology",
                findings=[{"finding": "This resembles an EGFR-mutant context", "evidence": ""}],
                concerns=["BRAF status unknown"],
            )
        ],
    )
    stats = annotate_grounding(op, _report("TP53"))
    assert "EGFR" in stats["offbriefing_genes"]
    assert "BRAF" in stats["offbriefing_genes"]
    assert op.grounding_flags
