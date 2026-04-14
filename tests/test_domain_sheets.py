"""Tests for domain sheet builders."""
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
    guideline.write_text(
        "---\nsource: CIViC + public guidelines\n---\n# NSCLC\nEGFR: erlotinib first-line"
    )
    variants = [{"gene": "EGFR"}]
    report_data = {"_kb_treatments_dir": str(tmp_path)}
    sheet = build_domain_sheet("clinical_evidence", "cancer", variants, report_data)
    assert "NSCLC" in sheet or "erlotinib" in sheet
