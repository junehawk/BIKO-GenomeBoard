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
