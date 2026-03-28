# tests/test_references.py
"""Tests for PMID reference infrastructure and AI-generated content watermarking."""
import json
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

from scripts.counselor.generate_pdf import generate_report_html
from scripts.tools.fetch_references import (
    KNOWN_REFERENCES,
    fetch_gene_references,
    search_pubmed,
)


# ── Fixtures ─────────────────────────────────────────────────────────────────

MOCK_ESEARCH_RESPONSE = {
    "esearchresult": {
        "idlist": ["25741868", "31006110"],
    }
}

MOCK_ESUMMARY_RESPONSE = {
    "result": {
        "25741868": {
            "title": "Standards and guidelines for the interpretation of sequence variants",
            "fulljournalname": "Genetics in Medicine",
            "pubdate": "2015 May",
            "source": "Genet Med",
        },
        "31006110": {
            "title": "Li-Fraumeni Syndrome",
            "fulljournalname": "GeneReviews",
            "pubdate": "2019",
            "source": "GeneReviews",
        },
    }
}

VARIANT_WITH_REFS = {
    "variant": "chr17:7577120:G>A",
    "gene": "TP53",
    "classification": "Pathogenic",
    "acmg_codes": ["PVS1", "PS1"],
    "conflict": False,
    "agents": {
        "clinical": {"clinvar_significance": "Pathogenic"},
        "korean_pop": {"korean_flag": ""},
    },
    "references": [
        {"pmid": "25741868", "source": "Genet Med 2015", "note": "ACMG/AMP standards"},
        {"pmid": "31006110", "source": "GeneReviews", "note": "Li-Fraumeni Syndrome"},
    ],
    "content_status": "ai-generated-with-references",
}

MINIMAL_REPORT = {
    "sample_id": "REF_TEST_001",
    "date": "2026-03-24",
    "variants": [VARIANT_WITH_REFS],
    "pgx_results": [],
    "summary": {"total": 1, "pathogenic": 1, "vus": 0, "benign": 0},
    "db_versions": {"clinvar": "2026-03-15", "gnomad": "4.0"},
}


# ── test_search_pubmed (mock) ─────────────────────────────────────────────────


def test_search_pubmed_returns_articles():
    """search_pubmed parses esearch + esummary responses into article dicts."""
    with patch("scripts.tools.fetch_references.fetch_with_retry") as mock_fetch:
        mock_fetch.side_effect = [MOCK_ESEARCH_RESPONSE, MOCK_ESUMMARY_RESPONSE]
        results = search_pubmed("TP53 clinical significance review", max_results=5)

    assert len(results) == 2
    assert results[0]["pmid"] == "25741868"
    assert "Genetics in Medicine" in results[0]["source"]
    assert results[0]["year"] == "2015"


def test_search_pubmed_empty_on_no_ids():
    """search_pubmed returns [] when esearch returns no IDs."""
    with patch("scripts.tools.fetch_references.fetch_with_retry") as mock_fetch:
        mock_fetch.return_value = {"esearchresult": {"idlist": []}}
        results = search_pubmed("nonexistent gene query xyz")
    assert results == []


def test_search_pubmed_empty_on_api_failure():
    """search_pubmed returns [] when fetch_with_retry returns None."""
    with patch("scripts.tools.fetch_references.fetch_with_retry") as mock_fetch:
        mock_fetch.return_value = None
        results = search_pubmed("TP53")
    assert results == []


# ── test_fetch_gene_references (mock) ─────────────────────────────────────────


def test_fetch_gene_references_deduplicates():
    """fetch_gene_references deduplicates PMIDs that appear in multiple searches."""
    duplicate_article = [{"pmid": "25741868", "title": "dup", "source": "Genet Med", "year": "2015"}]
    with patch("scripts.tools.fetch_references.search_pubmed", return_value=duplicate_article):
        results = fetch_gene_references("TP53")

    pmids = [r["pmid"] for r in results]
    assert len(pmids) == len(set(pmids)), "Duplicate PMIDs found in results"


def test_fetch_gene_references_max_eight():
    """fetch_gene_references caps results at 8 entries."""
    many_articles = [
        {"pmid": str(i), "title": f"Article {i}", "source": "Journal", "year": "2020"}
        for i in range(1, 12)
    ]
    with patch("scripts.tools.fetch_references.search_pubmed", return_value=many_articles):
        results = fetch_gene_references("TP53")
    assert len(results) <= 8


def test_fetch_gene_references_empty_on_all_failures():
    """fetch_gene_references returns [] when all PubMed searches fail."""
    with patch("scripts.tools.fetch_references.search_pubmed", return_value=[]):
        results = fetch_gene_references("UNKNOWNGENE999")
    assert results == []


# ── test_gene_knowledge_has_references ────────────────────────────────────────


def test_gene_knowledge_has_references_field():
    """Every gene entry in gene_knowledge.json must have a 'references' list."""
    path = Path(__file__).parent.parent / "data" / "gene_knowledge.json"
    data = json.loads(path.read_text())
    for entry in data["genes"]:
        gene = entry["gene"]
        assert "references" in entry, f"{gene} missing 'references' field"
        assert isinstance(entry["references"], list), f"{gene} 'references' must be a list"


def test_gene_knowledge_has_content_status():
    """Every gene entry must have a 'content_status' field."""
    path = Path(__file__).parent.parent / "data" / "gene_knowledge.json"
    data = json.loads(path.read_text())
    valid_statuses = {
        "ai-generated", "ai-generated-with-references", "verified",
        "curated-civic", "curated-ncbi", "curated-cpic", "curated-genreviews", "auto-minimal",
    }
    for entry in data["genes"]:
        gene = entry["gene"]
        assert "content_status" in entry, f"{gene} missing 'content_status' field"
        assert entry["content_status"] in valid_statuses, (
            f"{gene} has invalid content_status: {entry['content_status']!r}"
        )


def test_gene_knowledge_known_genes_have_pmids():
    """Genes with KNOWN_REFERENCES must have at least one PMID in gene_knowledge.json."""
    path = Path(__file__).parent.parent / "data" / "gene_knowledge.json"
    data = json.loads(path.read_text())
    entries_by_gene = {e["gene"]: e for e in data["genes"]}

    for gene in KNOWN_REFERENCES:
        if gene in entries_by_gene:
            refs = entries_by_gene[gene].get("references", [])
            assert len(refs) > 0, f"{gene} has KNOWN_REFERENCES but empty references in JSON"
            pmids = [r.get("pmid") for r in refs if r.get("pmid")]
            assert len(pmids) > 0, f"{gene} has KNOWN_REFERENCES but no entry with a 'pmid'"


# ── test_report_shows_pmid ────────────────────────────────────────────────────


def test_report_shows_pmid_cancer_template():
    """Cancer report template renders PMID citations when references are provided."""
    html = generate_report_html(MINIMAL_REPORT, mode="cancer")
    assert "PMID:25741868" in html
    assert "PMID:31006110" in html


def test_report_shows_pmid_rare_disease_template():
    """Rare-disease report template renders PMID citations when references are provided."""
    html = generate_report_html(MINIMAL_REPORT, mode="rare-disease")
    assert "PMID:25741868" in html
    assert "PMID:31006110" in html


def test_report_shows_reference_note():
    """Reference note text appears alongside PMID in rendered report."""
    html = generate_report_html(MINIMAL_REPORT, mode="cancer")
    assert "ACMG/AMP standards" in html


def test_report_shows_reference_source():
    """Reference source (journal name) appears in rendered report."""
    html = generate_report_html(MINIMAL_REPORT, mode="cancer")
    assert "Genet Med 2015" in html


def test_report_no_references_section_when_empty():
    """References section is suppressed when variant has no references and gene not in knowledge base."""
    data = {
        **MINIMAL_REPORT,
        "variants": [
            {
                "variant": "chr99:1:A>T",
                "gene": "FAKEGENE",
                "classification": "VUS",
                "acmg_codes": [],
                "conflict": False,
                "references": [],
            }
        ],
    }
    html = generate_report_html(data, mode="cancer")
    # Variant-level reference block should not appear (appendix may still have PMIDs)
    assert "[PMID:" not in html.split("FAKEGENE")[1].split("</div>")[0] if "FAKEGENE" in html else True


# ── test_ai_generated_watermark ──────────────────────────────────────────────


def test_ai_generated_with_references_watermark():
    """'ai-generated-with-references' status produces the correct watermark text."""
    html = generate_report_html(MINIMAL_REPORT, mode="cancer")
    assert "Referenced PMIDs provided for verification" in html


def test_ai_generated_plain_watermark():
    """'ai-generated' status produces the fallback watermark text."""
    data = {
        **MINIMAL_REPORT,
        "variants": [
            {
                **VARIANT_WITH_REFS,
                "content_status": "ai-generated",
                "references": [],
            }
        ],
    }
    html = generate_report_html(data, mode="cancer")
    assert "Not peer-reviewed" in html


def test_no_watermark_when_content_status_absent():
    """No watermark div is rendered when gene is not in knowledge base and content_status not set."""
    data = {
        **MINIMAL_REPORT,
        "variants": [
            {
                "variant": "chr99:1:A>T",
                "gene": "FAKEGENE",
                "classification": "VUS",
                "acmg_codes": [],
                "conflict": False,
            }
        ],
    }
    html = generate_report_html(data, mode="cancer")
    assert "Referenced PMIDs provided for verification" not in html
    assert "Not peer-reviewed" not in html


def test_generate_pdf_enriches_references_from_gene_knowledge():
    """generate_report_html enriches variant with references from CIViC in cancer mode."""
    bare_variant = {
        "variant": "chr17:7577120:G>A",
        "gene": "TP53",
        "classification": "Pathogenic",
        "acmg_codes": ["PVS1"],
        "conflict": False,
        "agents": {
            "clinical": {"clinvar_significance": "Pathogenic"},
            "korean_pop": {"korean_flag": ""},
        },
        # No references / content_status — filled from CIViC in cancer mode
    }
    report = {**MINIMAL_REPORT, "variants": [bare_variant]}
    html = generate_report_html(report, mode="cancer")
    # TP53 has CIViC evidence — references populated in cancer mode
    assert "PMID:" in html
