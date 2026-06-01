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
    """TMB / HRD / MSI / NSCLC are not genes and must not be flagged."""
    op = CancerBoardOpinion(
        therapeutic_headline="NSCLC, TMB-high, MSI-stable, HRD-positive",
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
    assert stats == {"offbriefing_genes": []}
