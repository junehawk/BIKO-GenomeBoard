# tests/test_generate_pdf.py
from scripts.counselor.generate_pdf import generate_report_html, generate_pdf
from scripts.common.models import Variant


# ── Shared fixture data ──────────────────────────────────────────────────────

MINIMAL_REPORT = {
    "sample_id": "TEST_001",
    "date": "2026-03-20",
    "variants": [
        {
            "variant": "chr17:7577120:G>A",
            "gene": "TP53",
            "classification": "Pathogenic",
            "acmg_codes": ["PVS1", "PS1"],
            "agents": {
                "clinical": {"clinvar_significance": "Pathogenic"},
                "korean_pop": {"korean_flag": "한국인 희귀 변이"},
            },
            "conflict": False,
        }
    ],
    "pgx_results": [],
    "summary": {"total": 1, "pathogenic": 1, "vus": 0, "benign": 0},
    "db_versions": {"clinvar": "2026-03-15", "gnomad": "4.0"},
}

FULL_REPORT = {
    "sample_id": "DEMO-KR-2026-001",
    "date": "2026-03-20",
    "summary": {"total": 3, "pathogenic": 1, "vus": 1, "benign": 1},
    "variants": [
        {
            "variant": "chr17:7577120:G>A",
            "gene": "TP53",
            "classification": "Pathogenic",
            "acmg_codes": ["PVS1", "PS1", "PM2_Supporting"],
            "conflict": False,
            "agents": {
                "clinical": {"clinvar_significance": "Pathogenic"},
                "korean_pop": {"korean_flag": "Korean freq 6.1x higher than global"},
            },
        },
        {
            "variant": "chr13:32337326:C>T",
            "gene": "BRCA2",
            "classification": "VUS",
            "acmg_codes": ["PM2_Supporting", "PP3"],
            "conflict": False,
            "agents": {
                "clinical": {"clinvar_significance": "Uncertain significance"},
                "korean_pop": {"korean_flag": ""},
            },
        },
        {
            "variant": "chr10:96541616:G>A",
            "gene": "CYP2C19",
            "classification": "Benign",
            "acmg_codes": ["BA1"],
            "conflict": False,
            "agents": {
                "clinical": {"clinvar_significance": "Drug Response"},
                "korean_pop": {"korean_flag": "Korean freq 5.6x higher (14pct vs 2.5pct)"},
            },
        },
    ],
    "pgx_results": [
        {
            "gene": "CYP2C19",
            "phenotype": "Poor Metabolizer (*2 carrier)",
            "cpic_level": "A",
            "korean_prevalence": 0.14,
            "western_prevalence": 0.025,
            "clinical_impact": "클로피도그렐 저반응 위험",
            "star_allele": "*2/*2",
            "cpic_recommendation": "Consider alternative antiplatelet therapy.",
        }
    ],
    "db_versions": {
        "ClinVar": "2026-03",
        "KRGDB": "v3.0 (2025-12)",
        "gnomAD": "v4.1 (2025)",
        "CPIC": "2022",
    },
}


# ── Core rendering tests ─────────────────────────────────────────────────────

def test_generate_html_minimal():
    """Template renders without error for minimal data."""
    html = generate_report_html(MINIMAL_REPORT)
    assert "GenomeBoard" in html
    assert "TP53" in html
    assert "Pathogenic" in html
    assert "Research Use Only" in html


def test_generate_html_contains_sample_id():
    html = generate_report_html(MINIMAL_REPORT)
    assert "TEST_001" in html


def test_generate_html_contains_date():
    html = generate_report_html(MINIMAL_REPORT)
    assert "2026-03-20" in html


# ── FoundationOne CDx structure tests ────────────────────────────────────────

def test_orange_accent_bar_present():
    """Header bar with orange accent color must appear on every page."""
    html = generate_report_html(MINIMAL_REPORT)
    assert "page-header-bar" in html
    assert "#E8712A" in html


def test_section_badges_present():
    """Section badge elements must appear for Genomic Findings."""
    html = generate_report_html(MINIMAL_REPORT)
    assert "section-badge" in html
    assert "Genomic Findings" in html


def test_variant_detail_page_gene_name():
    """Each variant gets a detail page with the gene name in large bold."""
    html = generate_report_html(MINIMAL_REPORT)
    assert "variant-gene-name" in html
    assert "TP53" in html


def test_three_column_detail_grid():
    """Detail grid with three columns must be present for variant pages."""
    html = generate_report_html(MINIMAL_REPORT)
    assert "detail-grid" in html
    assert "Clinical Significance" in html
    assert "Frequency" in html
    assert "Finding Summary" in html


def test_acmg_codes_rendered():
    """ACMG evidence codes appear as pills in detail section."""
    html = generate_report_html(MINIMAL_REPORT)
    assert "PVS1" in html
    assert "PS1" in html
    assert "acmg-code" in html


def test_clinvar_significance_rendered():
    html = generate_report_html(MINIMAL_REPORT)
    assert "ClinVar" in html


# ── Classification badge tests ────────────────────────────────────────────────

def test_pathogenic_badge():
    html = generate_report_html(MINIMAL_REPORT)
    assert "badge-pathogenic" in html


def test_vus_badge():
    data = {**MINIMAL_REPORT, "variants": [
        {**MINIMAL_REPORT["variants"][0], "classification": "VUS", "gene": "BRCA1"}
    ]}
    html = generate_report_html(data)
    assert "badge-vus" in html


def test_benign_badge():
    data = {**MINIMAL_REPORT, "variants": [
        {**MINIMAL_REPORT["variants"][0], "classification": "Benign", "gene": "APOE"}
    ]}
    html = generate_report_html(data)
    assert "badge-benign" in html


def test_likely_pathogenic_badge():
    data = {**MINIMAL_REPORT, "variants": [
        {**MINIMAL_REPORT["variants"][0], "classification": "Likely Pathogenic", "gene": "ATM"}
    ]}
    html = generate_report_html(data)
    assert "badge-likely-pathogenic" in html


# ── Korean population callout tests ──────────────────────────────────────────

def test_korean_callout_rendered_when_flag_present():
    """When korean_flag is populated, Korean callout box must appear."""
    html = generate_report_html(MINIMAL_REPORT)
    assert "korean-callout" in html
    assert "한국인 희귀 변이" in html


def test_korean_callout_absent_when_no_flag():
    """When korean_flag is empty, Korean callout section is suppressed."""
    data = {**MINIMAL_REPORT, "variants": [
        {**MINIMAL_REPORT["variants"][0],
         "agents": {"clinical": {"clinvar_significance": "Pathogenic"}, "korean_pop": {"korean_flag": ""}}}
    ]}
    html = generate_report_html(data)
    # The per-variant callout box should not appear (no korean population findings section on p1)
    assert "Korean-Specific Annotations" not in html


def test_korean_population_findings_section_on_summary():
    """Page 1 must show Korean Population Findings section when flags present."""
    html = generate_report_html(FULL_REPORT)
    assert "Korean Population Findings" in html
    assert "Korean-Specific Annotations" in html


# ── Clinically significant vs VUS/Benign split ───────────────────────────────

def test_clinically_significant_panel():
    html = generate_report_html(FULL_REPORT)
    assert "Clinically Significant" in html


def test_vus_benign_panel():
    html = generate_report_html(FULL_REPORT)
    assert "VUS" in html or "Benign" in html


def test_no_reportable_list_shows_vus_benign():
    html = generate_report_html(FULL_REPORT)
    assert "No Reportable Therapeutic Options" in html


# ── PGx section tests ────────────────────────────────────────────────────────

def test_pgx_section_absent_when_empty():
    html = generate_report_html(MINIMAL_REPORT)
    # PGx detail page should not appear when pgx_results is empty
    assert "Pharmacogenomics (PGx) Detailed Findings" not in html


def test_pgx_section_present_when_populated():
    html = generate_report_html(FULL_REPORT)
    assert "Pharmacogenomics" in html
    assert "CYP2C19" in html
    assert "CPIC Level A" in html


def test_pgx_korean_prevalence_rendered():
    html = generate_report_html(FULL_REPORT)
    # Korean prevalence 14% should appear formatted
    assert "14.0%" in html


def test_pgx_enrichment_flag():
    """When Korean prevalence >= 2x Western, enrichment flag must appear."""
    html = generate_report_html(FULL_REPORT)
    assert "higher" in html


# ── DB versions / methodology ─────────────────────────────────────────────────

def test_db_versions_rendered():
    html = generate_report_html(MINIMAL_REPORT)
    assert "clinvar" in html.lower() or "ClinVar" in html
    assert "gnomad" in html.lower() or "gnomAD" in html


def test_methodology_section_present():
    html = generate_report_html(MINIMAL_REPORT)
    assert "Methodology" in html
    assert "Pipeline" in html or "pipeline" in html


def test_limitations_section_present():
    html = generate_report_html(MINIMAL_REPORT)
    assert "Limitations" in html or "limitations" in html


# ── Optional new fields — graceful handling ───────────────────────────────────

def test_optional_treatment_strategies_absent_gracefully():
    """Template must not error when treatment_strategies is absent."""
    html = generate_report_html(MINIMAL_REPORT)
    assert html is not None
    assert len(html) > 100


def test_optional_treatment_strategies_rendered_when_present():
    data = {**MINIMAL_REPORT, "variants": [
        {**MINIMAL_REPORT["variants"][0], "treatment_strategies": "PARP inhibitor — olaparib"}
    ]}
    html = generate_report_html(data)
    assert "olaparib" in html


def test_optional_frequency_prognosis_rendered_when_present():
    data = {**MINIMAL_REPORT, "variants": [
        {**MINIMAL_REPORT["variants"][0], "frequency_prognosis": "Rare variant; poor prognosis associated"}
    ]}
    html = generate_report_html(data)
    assert "poor prognosis" in html


def test_optional_finding_summary_rendered_when_present():
    data = {**MINIMAL_REPORT, "variants": [
        {**MINIMAL_REPORT["variants"][0], "finding_summary": "TP53 is a well-established tumor suppressor."}
    ]}
    html = generate_report_html(data)
    assert "tumor suppressor" in html


# ── Missing agents data — graceful handling ──────────────────────────────────

def test_missing_agents_field_graceful():
    """Template must not error when agents field is entirely absent."""
    data = {**MINIMAL_REPORT, "variants": [
        {
            "variant": "chr1:100:A>G",
            "gene": "GENE1",
            "classification": "VUS",
            "acmg_codes": [],
            "conflict": False,
            # agents intentionally omitted
        }
    ]}
    html = generate_report_html(data)
    assert "GENE1" in html


def test_missing_korean_pop_agent_graceful():
    """Template must not error when korean_pop sub-key is absent."""
    data = {**MINIMAL_REPORT, "variants": [
        {
            "variant": "chr1:100:A>G",
            "gene": "GENE2",
            "classification": "Pathogenic",
            "acmg_codes": ["PVS1"],
            "conflict": False,
            "agents": {
                "clinical": {"clinvar_significance": "Pathogenic"},
                # korean_pop intentionally omitted
            },
        }
    ]}
    html = generate_report_html(data)
    assert "GENE2" in html


def test_empty_variants_list():
    """Template must render gracefully with zero variants."""
    data = {
        "sample_id": "EMPTY_001",
        "date": "2026-03-20",
        "variants": [],
        "pgx_results": [],
        "summary": {"total": 0, "pathogenic": 0, "vus": 0, "benign": 0},
        "db_versions": {},
    }
    html = generate_report_html(data)
    assert "EMPTY_001" in html
    assert "GenomeBoard" in html


def test_conflict_flag_rendered():
    """Conflict flag badge must appear when conflict=True."""
    data = {**MINIMAL_REPORT, "variants": [
        {**MINIMAL_REPORT["variants"][0], "conflict": True}
    ]}
    html = generate_report_html(data)
    assert "Conflict" in html


# ── Print / layout structural tests ──────────────────────────────────────────

def test_a4_page_size_in_css():
    html = generate_report_html(MINIMAL_REPORT)
    assert "A4" in html


def test_page_break_before_in_css():
    html = generate_report_html(MINIMAL_REPORT)
    assert "page-break-before" in html


def test_footer_on_every_rendered_page():
    """page-footer class must appear for each page (summary + variant + methodology)."""
    html = generate_report_html(MINIMAL_REPORT)
    # 1 summary page + 1 variant page + 1 methodology page = at least 3 footers
    assert html.count("page-footer") >= 3


def test_masthead_brand_present():
    html = generate_report_html(MINIMAL_REPORT)
    assert "masthead-brand" in html


def test_navy_color_in_css():
    html = generate_report_html(MINIMAL_REPORT)
    assert "#1c2b4a" in html
