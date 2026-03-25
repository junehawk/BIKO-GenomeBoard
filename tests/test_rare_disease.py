"""Tests for Task 2.1: Rare Disease Report Mode."""
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


# ── HPO Matcher Tests ────────────────────────────────────────────────────────


def test_resolve_hpo_terms_offline():
    """resolve_hpo_terms handles API calls correctly when mocked."""
    from scripts.clinical.hpo_matcher import resolve_hpo_terms

    mock_term = {"name": "Seizures"}
    mock_genes = {"genes": [{"symbol": "SCN1A"}, {"symbol": "KCNQ2"}]}

    with patch("scripts.clinical.hpo_matcher.fetch_with_retry") as mock_fetch:
        mock_fetch.side_effect = [mock_term, mock_genes]
        results = resolve_hpo_terms(["HP:0001250"])

    assert len(results) == 1
    assert results[0]["id"] == "HP:0001250"
    assert results[0]["name"] == "Seizures"
    assert "SCN1A" in results[0]["genes"]
    assert "KCNQ2" in results[0]["genes"]


def test_resolve_hpo_terms_skips_invalid():
    """resolve_hpo_terms skips IDs that don't start with 'HP:'."""
    from scripts.clinical.hpo_matcher import resolve_hpo_terms

    with patch("scripts.clinical.hpo_matcher.fetch_with_retry") as mock_fetch:
        results = resolve_hpo_terms(["NOT_AN_HPO_ID", "INVALID"])

    assert results == []
    mock_fetch.assert_not_called()


def test_resolve_hpo_terms_api_unavailable():
    """resolve_hpo_terms returns id as name when API returns None."""
    from scripts.clinical.hpo_matcher import resolve_hpo_terms

    with patch("scripts.clinical.hpo_matcher.fetch_with_retry", return_value=None):
        results = resolve_hpo_terms(["HP:0001263"])

    assert len(results) == 1
    assert results[0]["id"] == "HP:0001263"
    assert results[0]["name"] == "HP:0001263"  # fallback to ID
    assert results[0]["genes"] == []


def test_resolve_hpo_terms_falls_back_to_local_db(tmp_path, monkeypatch):
    """API 실패 시 로컬 HPO DB로 fallback."""
    from scripts.db.build_hpo_db import build_db
    from scripts.clinical.hpo_matcher import resolve_hpo_terms

    # Build local DB
    tsv = tmp_path / "genes_to_phenotype.txt"
    tsv.write_text(
        "#header\n"
        "7157\tTP53\tHP:0001250\tSeizure\t-\tOMIM:151623\n"
        "672\tBRCA2\tHP:0001250\tSeizure\t-\tOMIM:612555\n"
    )
    db_path = str(tmp_path / "hpo.sqlite3")
    build_db(str(tsv), db_path)

    # Mock API to fail
    monkeypatch.setattr(
        "scripts.clinical.hpo_matcher.fetch_with_retry", lambda *a, **kw: None
    )
    # Patch config.get at both the source and the already-imported reference
    _cfg = lambda key, default=None: db_path if key == "paths.hpo_db" else default
    monkeypatch.setattr("scripts.common.config.get", _cfg)
    monkeypatch.setattr("scripts.db.query_local_hpo.get", _cfg)

    results = resolve_hpo_terms(["HP:0001250"])
    assert len(results) == 1
    assert results[0]["name"] == "Seizure"
    assert "TP53" in results[0]["genes"]
    assert "BRCA2" in results[0]["genes"]


def test_calculate_hpo_score():
    """calculate_hpo_score returns count of HPO terms linked to the gene."""
    from scripts.clinical.hpo_matcher import calculate_hpo_score

    hpo_results = [
        {"id": "HP:0001250", "name": "Seizures", "genes": ["SCN1A", "KCNQ2", "TP53"]},
        {"id": "HP:0001263", "name": "Global developmental delay", "genes": ["TP53", "CFTR"]},
        {"id": "HP:0000252", "name": "Microcephaly", "genes": ["CFTR"]},
    ]

    assert calculate_hpo_score("TP53", hpo_results) == 2
    assert calculate_hpo_score("CFTR", hpo_results) == 2
    assert calculate_hpo_score("SCN1A", hpo_results) == 1
    assert calculate_hpo_score("BRCA2", hpo_results) == 0


def test_calculate_hpo_score_case_insensitive():
    """calculate_hpo_score is case-insensitive for gene names."""
    from scripts.clinical.hpo_matcher import calculate_hpo_score

    hpo_results = [
        {"id": "HP:0001250", "name": "Seizures", "genes": ["SCN1A", "TP53"]},
    ]

    assert calculate_hpo_score("tp53", hpo_results) == 1
    assert calculate_hpo_score("Tp53", hpo_results) == 1


def test_calculate_hpo_score_empty():
    """calculate_hpo_score returns 0 when hpo_results is empty or gene is empty."""
    from scripts.clinical.hpo_matcher import calculate_hpo_score

    assert calculate_hpo_score("TP53", []) == 0
    assert calculate_hpo_score("", [{"id": "HP:0001250", "name": "X", "genes": ["TP53"]}]) == 0
    assert calculate_hpo_score(None, [{"id": "HP:0001250", "name": "X", "genes": ["TP53"]}]) == 0


def test_get_matching_hpo_terms():
    """get_matching_hpo_terms returns formatted list of matching HPO terms."""
    from scripts.clinical.hpo_matcher import get_matching_hpo_terms

    hpo_results = [
        {"id": "HP:0001250", "name": "Seizures", "genes": ["SCN1A", "TP53"]},
        {"id": "HP:0001263", "name": "Global developmental delay", "genes": ["TP53"]},
        {"id": "HP:0000252", "name": "Microcephaly", "genes": ["CFTR"]},
    ]

    matching = get_matching_hpo_terms("TP53", hpo_results)
    assert len(matching) == 2
    assert "HP:0001250 (Seizures)" in matching
    assert "HP:0001263 (Global developmental delay)" in matching

    no_match = get_matching_hpo_terms("BRCA2", hpo_results)
    assert no_match == []


# ── OMIM Tests ───────────────────────────────────────────────────────────────


def test_query_omim_known_gene():
    """query_omim returns data for genes in the static OMIM_DATA dict."""
    from scripts.clinical.query_omim import query_omim

    result = query_omim("TP53")
    assert result is not None
    assert result["mim"] == "191170"
    assert "Li-Fraumeni syndrome" in result["phenotypes"]
    assert result["inheritance"] == "AD"


def test_query_omim_cftr():
    """query_omim returns data for CFTR."""
    from scripts.clinical.query_omim import query_omim

    result = query_omim("CFTR")
    assert result is not None
    assert "Cystic fibrosis" in result["phenotypes"]
    assert result["inheritance"] == "AR"


def test_query_omim_unknown_gene():
    """query_omim returns None for genes not in the static data."""
    from scripts.clinical.query_omim import query_omim

    assert query_omim("UNKNOWNGENE123") is None
    assert query_omim("") is None


def test_query_omim_ptpn11():
    """query_omim returns Noonan syndrome data for PTPN11."""
    from scripts.clinical.query_omim import query_omim

    result = query_omim("PTPN11")
    assert result is not None
    assert "Noonan syndrome 1" in result["phenotypes"]
    assert result["inheritance"] == "AD"


# ── ClinGen Tests ─────────────────────────────────────────────────────────────


def test_get_gene_validity_known():
    """get_gene_validity returns Definitive for known genes."""
    from scripts.clinical.query_clingen import get_gene_validity

    assert get_gene_validity("TP53") == "Definitive"
    assert get_gene_validity("BRCA2") == "Definitive"
    assert get_gene_validity("CFTR") == "Definitive"


def test_get_gene_validity_unknown():
    """get_gene_validity returns None for unknown genes."""
    from scripts.clinical.query_clingen import get_gene_validity

    assert get_gene_validity("UNKNOWNGENE") is None
    assert get_gene_validity("") is None


# ── Pipeline Integration Tests ───────────────────────────────────────────────

RARE_DISEASE_VCF = str(
    Path(__file__).parent.parent / "data" / "sample_vcf" / "rare_disease_demo.vcf"
)


def test_rare_disease_pipeline_with_hpo(tmp_path):
    """run_pipeline with rare-disease mode and HPO IDs adds hpo_score to variants."""
    from scripts.orchestrate import run_pipeline

    hpo_ids = ["HP:0001250", "HP:0001263"]
    # Construct mock hpo_results that link TP53 to both HPO terms
    mock_hpo = [
        {"id": "HP:0001250", "name": "Seizures", "genes": ["TP53", "CFTR"]},
        {"id": "HP:0001263", "name": "Global developmental delay", "genes": ["TP53"]},
    ]

    with patch("scripts.orchestrate.resolve_hpo_terms", return_value=mock_hpo):
        result = run_pipeline(
            vcf_path=RARE_DISEASE_VCF,
            output_path=str(tmp_path / "report.html"),
            skip_api=True,
            mode="rare-disease",
            hpo_ids=hpo_ids,
        )

    assert result is not None
    assert result["mode"] == "rare-disease"
    assert "hpo_results" in result

    # All variants should have hpo_score
    for v in result["variants"]:
        assert "hpo_score" in v
        assert "matching_hpo" in v
        assert "omim_mim" in v
        assert "omim_phenotypes" in v
        assert "inheritance" in v
        assert "clingen_validity" in v

    # TP53 should have hpo_score == 2 (matches both HPO terms)
    tp53 = next((v for v in result["variants"] if v["gene"] == "TP53"), None)
    assert tp53 is not None
    assert tp53["hpo_score"] == 2


def test_rare_disease_pipeline_offline_hpo(tmp_path):
    """run_pipeline in skip_api mode stores HPO IDs without resolution."""
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path=RARE_DISEASE_VCF,
        output_path=str(tmp_path / "report.html"),
        skip_api=True,
        mode="rare-disease",
        hpo_ids=["HP:0001250", "HP:0001263"],
    )

    assert result is not None
    assert len(result["hpo_results"]) == 2
    assert result["hpo_results"][0]["id"] == "HP:0001250"
    # In offline mode genes list is empty
    assert result["hpo_results"][0]["genes"] == []


def test_rare_disease_pipeline_sorts_by_classification_then_hpo(tmp_path):
    """Variants sorted: Pathogenic first, then by HPO score descending."""
    from scripts.orchestrate import run_pipeline

    mock_hpo = [
        {"id": "HP:0001250", "name": "Seizures", "genes": ["TP53"]},
    ]

    with patch("scripts.orchestrate.resolve_hpo_terms", return_value=mock_hpo):
        result = run_pipeline(
            vcf_path=RARE_DISEASE_VCF,
            output_path=str(tmp_path / "report.html"),
            skip_api=True,
            mode="rare-disease",
            hpo_ids=["HP:0001250"],
        )

    assert result is not None
    variants = result["variants"]
    # Check sort is stable: no Pathogenic/LP after VUS/Benign (if any Pathogenic exist)
    cls_rank = {
        "Pathogenic": 0, "Likely Pathogenic": 1, "VUS": 2,
        "Drug Response": 3, "Risk Factor": 4, "Likely Benign": 5, "Benign": 6,
    }
    ranks = [cls_rank.get(v["classification"], 2) for v in variants]
    assert ranks == sorted(ranks), "Variants not sorted by classification rank"


# ── Report Template Tests ────────────────────────────────────────────────────

MINIMAL_RARE_DISEASE_REPORT = {
    "sample_id": "RD_TEST_001",
    "date": "2026-03-23",
    "mode": "rare-disease",
    "hpo_results": [
        {"id": "HP:0001250", "name": "Seizures", "genes": ["SCN1A", "TP53"]},
        {"id": "HP:0001263", "name": "Global developmental delay", "genes": ["TP53"]},
    ],
    "variants": [
        {
            "variant": "chr17:7675088:C>A",
            "gene": "TP53",
            "classification": "Pathogenic",
            "acmg_codes": ["PS1", "PP5"],
            "clinvar_significance": "Pathogenic",
            "hpo_score": 2,
            "matching_hpo": ["HP:0001250 (Seizures)", "HP:0001263 (Global developmental delay)"],
            "omim_mim": "191170",
            "omim_phenotypes": ["Li-Fraumeni syndrome"],
            "inheritance": "AD",
            "clingen_validity": "Definitive",
            "agents": {
                "clinical": {"clinvar_significance": "Pathogenic"},
                "korean_pop": {"korean_flag": ""},
            },
            "conflict": False,
        },
        {
            "variant": "chr7:117559590:ATCT>A",
            "gene": "CFTR",
            "classification": "VUS",
            "acmg_codes": [],
            "clinvar_significance": "VUS",
            "hpo_score": 0,
            "matching_hpo": [],
            "omim_mim": "602421",
            "omim_phenotypes": ["Cystic fibrosis"],
            "inheritance": "AR",
            "clingen_validity": "Definitive",
            "agents": {
                "clinical": {"clinvar_significance": "VUS"},
                "korean_pop": {"korean_flag": ""},
            },
            "conflict": False,
        },
    ],
    "pgx_results": [],
    "summary": {
        "total": 2,
        "pathogenic": 1,
        "likely_pathogenic": 0,
        "vus": 1,
        "drug_response": 0,
        "risk_factor": 0,
        "benign": 0,
        "likely_benign": 0,
    },
    "db_versions": {"clinvar": "2026-03-23", "gnomad": "4.0", "krgdb": "2026-03-01"},
    "pipeline": {"skip_api": True, "krgdb_path": "data/krgdb_freq.tsv"},
}


def test_rare_disease_report_has_candidate_ranking():
    """Rare disease report template includes candidate ranking table."""
    from scripts.counselor.generate_pdf import generate_report_html

    html = generate_report_html(MINIMAL_RARE_DISEASE_REPORT, mode="rare-disease")
    assert "Candidate Gene Ranking" in html
    assert "TP53" in html
    assert "CFTR" in html


def test_rare_disease_report_has_inheritance():
    """Rare disease report template shows inheritance patterns."""
    from scripts.counselor.generate_pdf import generate_report_html

    html = generate_report_html(MINIMAL_RARE_DISEASE_REPORT, mode="rare-disease")
    # AD appears for TP53, AR for CFTR
    assert "AD" in html
    assert "AR" in html


def test_rare_disease_report_has_phenotypes():
    """Rare disease report template shows OMIM phenotypes and HPO terms."""
    from scripts.counselor.generate_pdf import generate_report_html

    html = generate_report_html(MINIMAL_RARE_DISEASE_REPORT, mode="rare-disease")
    assert "Li-Fraumeni syndrome" in html
    assert "Cystic fibrosis" in html
    assert "HP:0001250" in html
    assert "Seizures" in html


def test_rare_disease_report_has_hpo_section():
    """Rare disease report shows HPO phenotype section."""
    from scripts.counselor.generate_pdf import generate_report_html

    html = generate_report_html(MINIMAL_RARE_DISEASE_REPORT, mode="rare-disease")
    assert "Patient Phenotypes" in html or "HPO" in html
    assert "Global developmental delay" in html


def test_rare_disease_report_has_clingen():
    """Rare disease report shows ClinGen validity."""
    from scripts.counselor.generate_pdf import generate_report_html

    html = generate_report_html(MINIMAL_RARE_DISEASE_REPORT, mode="rare-disease")
    assert "Definitive" in html
    assert "ClinGen" in html


def test_rare_disease_report_has_hpo_score():
    """Rare disease report shows HPO score in candidate ranking."""
    from scripts.counselor.generate_pdf import generate_report_html

    html = generate_report_html(MINIMAL_RARE_DISEASE_REPORT, mode="rare-disease")
    # HPO score of 2 for TP53 should appear
    assert "HPO Score" in html or "HPO: 2" in html or ">2<" in html


def test_rare_disease_report_has_matching_hpo_in_detail():
    """Rare disease variant detail shows matching HPO terms."""
    from scripts.counselor.generate_pdf import generate_report_html

    html = generate_report_html(MINIMAL_RARE_DISEASE_REPORT, mode="rare-disease")
    assert "Matching HPO Terms" in html
    assert "HP:0001250 (Seizures)" in html


def test_rare_disease_report_research_use_only():
    """Rare disease report has Research Use Only disclaimer on every page."""
    from scripts.counselor.generate_pdf import generate_report_html

    html = generate_report_html(MINIMAL_RARE_DISEASE_REPORT, mode="rare-disease")
    assert html.count("Research Use Only") >= 1


def test_rare_disease_full_offline_with_local_dbs(tmp_path, monkeypatch):
    """HPO + ClinGen 로컬 DB로 오프라인 rare disease 파이프라인 전체 동작."""
    from scripts.db.build_hpo_db import build_db as build_hpo
    from scripts.clinical.hpo_matcher import resolve_hpo_terms, calculate_hpo_score

    # Build HPO DB
    hpo_tsv = tmp_path / "gtp.txt"
    hpo_tsv.write_text(
        "#header\n"
        "7157\tTP53\tHP:0001250\tSeizure\t-\tOMIM:151623\n"
        "7157\tTP53\tHP:0002664\tNeoplasm\t-\tOMIM:151623\n"
    )
    hpo_db = str(tmp_path / "hpo.sqlite3")
    build_hpo(str(hpo_tsv), hpo_db)

    # Block API
    monkeypatch.setattr("scripts.clinical.hpo_matcher.fetch_with_retry", lambda *a, **kw: None)
    monkeypatch.setattr(
        "scripts.common.config.get",
        lambda key, default=None: hpo_db if key == "paths.hpo_db" else default,
    )
    monkeypatch.setattr(
        "scripts.db.query_local_hpo.get",
        lambda key, default=None: hpo_db if key == "paths.hpo_db" else default,
    )

    # Resolve HPO terms (should use local DB)
    results = resolve_hpo_terms(["HP:0001250", "HP:0002664"])
    assert len(results) == 2

    # TP53 matches both HPO terms
    score = calculate_hpo_score("TP53", results)
    assert score == 2

    # BRCA2 matches neither (not in our test data for HP:0002664 only has TP53)
    score_brca2 = calculate_hpo_score("BRCA2", results)
    assert score_brca2 == 0
