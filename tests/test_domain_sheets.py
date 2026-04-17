"""Tests for domain sheet builders."""

from unittest.mock import patch

from scripts.clinical_board.domain_sheets import MAX_DOMAIN_CHARS, build_domain_sheet


def test_variant_pathologist_sheet_rare_disease():
    """Rare disease Variant Pathologist sheet includes ClinVar + in silico."""
    variants = [
        {
            "gene": "TP53",
            "classification": "Pathogenic",
            "clinvar_significance": "Pathogenic",
            "review_status": "expert panel",
            "in_silico": {"revel": 0.95, "cadd_phred": 35.0, "spliceai_max": 0.01},
        }
    ]
    sheet = build_domain_sheet("variant_pathology", "rare-disease", variants, {})
    assert "TP53" in sheet
    assert "REVEL" in sheet or "revel" in sheet.lower()
    assert "SpliceAI" in sheet or "spliceai" in sheet.lower()


def test_disease_geneticist_sheet_rare_disease():
    """Rare disease Disease Geneticist sheet includes OMIM + HPO."""
    variants = [
        {
            "gene": "CFTR",
            "omim_phenotypes": ["Cystic fibrosis"],
            "inheritance": "AR",
            "matching_hpo": [{"id": "HP:0002110", "name": "Bronchiectasis"}],
        }
    ]
    sheet = build_domain_sheet("disease_genetics", "rare-disease", variants, {})
    assert "CFTR" in sheet
    assert "Cystic fibrosis" in sheet
    assert "AR" in sheet


def test_disease_geneticist_sheet_accepts_str_hpo_and_phenotypes():
    """Real pipeline emits matching_hpo/omim_phenotypes as plain strings; builder must not crash."""
    variants = [
        {
            "gene": "TP53",
            "omim_phenotypes": ["Li-Fraumeni syndrome"],
            "inheritance": "AD",
            "matching_hpo": ["HP:0001263", "HP:0001250"],
        }
    ]
    sheet = build_domain_sheet("disease_genetics", "rare-disease", variants, {})
    assert "TP53" in sheet
    assert "HP:0001263" in sheet
    assert "Li-Fraumeni" in sheet


def test_domain_sheet_respects_char_limit():
    """Domain sheet is truncated to MAX_DOMAIN_CHARS."""
    variants = [{"gene": f"GENE{i}", "classification": "VUS"} for i in range(100)]
    sheet = build_domain_sheet("variant_pathology", "rare-disease", variants, {})
    assert len(sheet) <= MAX_DOMAIN_CHARS


def test_domain_sheet_empty_variants():
    """Empty variants produce minimal sheet."""
    sheet = build_domain_sheet("variant_pathology", "rare-disease", [], {})
    assert isinstance(sheet, str)
    assert len(sheet) < 200


def test_therapeutic_target_sheet_cancer():
    """Cancer Therapeutic Target sheet includes CIViC drug evidence."""
    variants = [
        {
            "gene": "EGFR",
            "classification": "Pathogenic",
            "civic_evidence": [{"drug": "Erlotinib", "level": "A", "direction": "Supports"}],
            "in_silico": {"revel": 0.9},
        }
    ]
    sheet = build_domain_sheet("therapeutic_target", "cancer", variants, {})
    assert "EGFR" in sheet
    assert "Erlotinib" in sheet


def test_tumor_genomics_sheet_cancer():
    """Cancer Tumor Genomics sheet includes TMB and VAF."""
    variants = [{"gene": "TP53", "vaf": 0.45}]
    report_data = {"tmb": {"score": 12.5, "level": "High"}}
    sheet = build_domain_sheet("tumor_genomics", "cancer", variants, report_data)
    assert "TMB" in sheet
    assert "12.5" in sheet


def test_clinical_evidence_sheet_loads_guidelines(tmp_path):
    """Clinical Evidence sheet loads treatment guidelines from KB."""
    guideline = tmp_path / "lung_cancer_nsclc.md"
    guideline.write_text("---\nsource: CIViC + public guidelines\n---\n# NSCLC\nEGFR: erlotinib first-line")
    variants = [{"gene": "EGFR"}]
    report_data = {"_kb_treatments_dir": str(tmp_path)}
    sheet = build_domain_sheet("clinical_evidence", "cancer", variants, report_data)
    assert "NSCLC" in sheet or "erlotinib" in sheet


# ---------------------------------------------------------------------------
# CIViC literature grounding tests (Phase 1)
# ---------------------------------------------------------------------------

_MOCK_CIVIC_EVIDENCE = [
    {
        "gene": "BRAF",
        "variant": "V600E",
        "disease": "Melanoma",
        "therapies": "Vemurafenib",
        "evidence_type": "Predictive",
        "evidence_level": "A",
        "significance": "Sensitivity/Response",
        "statement": "Vemurafenib showed significant clinical benefit in BRAF V600E melanoma patients.",
        "pmid": "22735384",
        "citation": "Chapman PB et al., N Engl J Med, 2011",
        "nct_ids": "",
    },
    {
        "gene": "BRAF",
        "variant": "V600E",
        "disease": "Colorectal Cancer",
        "therapies": "Encorafenib + Cetuximab",
        "evidence_type": "Predictive",
        "evidence_level": "B",
        "significance": "Sensitivity/Response",
        "statement": "BEACON CRC trial demonstrated improved OS with encorafenib plus cetuximab in BRAF V600E mCRC.",
        "pmid": "31566309",
        "citation": "Kopetz S et al., N Engl J Med, 2019",
        "nct_ids": "NCT02928224",
    },
]


def test_clinical_evidence_sheet_includes_civic_pmid():
    """Cancer clinical_evidence domain sheet includes CIViC PMID from DB query."""
    variants = [{"gene": "BRAF", "hgvsp": "p.Val600Glu"}]
    with patch(
        "scripts.db.query_civic.get_variant_evidence",
        return_value=_MOCK_CIVIC_EVIDENCE,
    ):
        sheet = build_domain_sheet("clinical_evidence", "cancer", variants, {})
    assert "PMID:22735384" in sheet
    assert "PMID:31566309" in sheet


def test_clinical_evidence_sheet_includes_evidence_statement():
    """Cancer clinical_evidence domain sheet includes evidence statement text."""
    variants = [{"gene": "BRAF"}]
    with patch(
        "scripts.db.query_civic.get_variant_evidence",
        return_value=_MOCK_CIVIC_EVIDENCE,
    ):
        sheet = build_domain_sheet("clinical_evidence", "cancer", variants, {})
    assert "Vemurafenib showed significant clinical benefit" in sheet
    assert "Level A" in sheet


def test_clinical_evidence_sheet_includes_citation():
    """Cancer clinical_evidence domain sheet includes citation information."""
    variants = [{"gene": "BRAF"}]
    with patch(
        "scripts.db.query_civic.get_variant_evidence",
        return_value=_MOCK_CIVIC_EVIDENCE,
    ):
        sheet = build_domain_sheet("clinical_evidence", "cancer", variants, {})
    assert "Chapman PB" in sheet
    assert "Kopetz S" in sheet


def test_literature_sheet_includes_civic_pmid():
    """Rare-disease literature_evidence domain sheet includes CIViC PMID."""
    variants = [{"gene": "BRAF", "hgvsp": "p.Val600Glu"}]
    with patch(
        "scripts.db.query_civic.get_variant_evidence",
        return_value=_MOCK_CIVIC_EVIDENCE,
    ):
        sheet = build_domain_sheet("literature_evidence", "rare-disease", variants, {})
    assert "PMID:22735384" in sheet
    assert "CIViC Literature Evidence" in sheet


def test_literature_sheet_includes_evidence_statement():
    """Rare-disease literature_evidence sheet includes evidence statement text."""
    variants = [{"gene": "BRAF"}]
    with patch(
        "scripts.db.query_civic.get_variant_evidence",
        return_value=_MOCK_CIVIC_EVIDENCE,
    ):
        sheet = build_domain_sheet("literature_evidence", "rare-disease", variants, {})
    assert "BEACON CRC trial" in sheet
    assert "Level B" in sheet


def test_domain_sheet_graceful_without_civic():
    """Domain sheet degrades gracefully when CIViC DB is unavailable."""
    variants = [{"gene": "TP53"}]
    with patch(
        "scripts.db.query_civic.get_variant_evidence",
        side_effect=Exception("DB not found"),
    ):
        sheet = build_domain_sheet("clinical_evidence", "cancer", variants, {})
    # Should not crash; should still produce a valid sheet
    assert "TP53" in sheet
    assert "CIViC Literature Evidence" not in sheet


def test_domain_sheet_graceful_import_error():
    """Domain sheet degrades gracefully when query_civic module is missing."""
    variants = [{"gene": "TP53"}]
    with patch(
        "scripts.clinical_board.domain_sheets._get_civic_literature_evidence",
        return_value=[],
    ):
        sheet = build_domain_sheet("literature_evidence", "rare-disease", variants, {})
    assert "TP53" in sheet
    assert isinstance(sheet, str)


def test_civic_evidence_dedup_across_same_gene():
    """CIViC evidence is queried once per gene even with multiple variants."""
    variants = [
        {"gene": "BRAF", "hgvsp": "p.Val600Glu"},
        {"gene": "BRAF", "hgvsp": "p.Val600Lys"},
    ]
    with patch(
        "scripts.db.query_civic.get_variant_evidence",
        return_value=_MOCK_CIVIC_EVIDENCE,
    ) as mock_query:
        sheet = build_domain_sheet("clinical_evidence", "cancer", variants, {})
    # Should call get_variant_evidence only once for BRAF
    assert mock_query.call_count == 1
    assert "CIViC Literature Evidence" in sheet


def test_civic_evidence_long_statement_truncated():
    """Long evidence statements are truncated to avoid bloating the domain sheet."""
    long_evidence = [
        {
            "gene": "EGFR",
            "variant": "L858R",
            "disease": "NSCLC",
            "therapies": "Erlotinib",
            "evidence_type": "Predictive",
            "evidence_level": "A",
            "significance": "Sensitivity/Response",
            "statement": "A" * 300,
            "pmid": "12345678",
            "citation": "Test Author, J Test, 2024",
            "nct_ids": "",
        }
    ]
    variants = [{"gene": "EGFR"}]
    with patch(
        "scripts.db.query_civic.get_variant_evidence",
        return_value=long_evidence,
    ):
        sheet = build_domain_sheet("clinical_evidence", "cancer", variants, {})
    # Statement should be truncated — no 300-char run of 'A'
    assert "A" * 300 not in sheet
    assert "A" * 197 + "..." in sheet
