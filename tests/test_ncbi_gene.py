import pytest


def test_fetch_gene_summary_known_gene(monkeypatch):
    """TP53 gene summary 조회."""
    from scripts.tools.sources.ncbi_gene import fetch_gene_summary

    # Mock esearch response
    mock_search = {"esearchresult": {"idlist": ["7157"]}}
    # Mock esummary response
    mock_summary = {"result": {"7157": {
        "uid": "7157",
        "name": "TP53",
        "description": "tumor protein p53",
        "summary": "This gene encodes a tumor protein that responds to diverse cellular stresses.",
        "otheraliases": "p53, LFS1",
    }}}

    call_count = [0]
    def mock_fetch(url, **kw):
        call_count[0] += 1
        if "esearch" in url:
            return mock_search
        if "esummary" in url:
            return mock_summary
        return None

    monkeypatch.setattr("scripts.tools.sources.ncbi_gene.fetch_with_retry", mock_fetch)

    result = fetch_gene_summary("TP53")
    assert result is not None
    assert result["gene"] == "TP53"
    assert result["entrez_id"] == "7157"
    assert "tumor protein" in result["summary"].lower() or "tumor protein" in result["full_name"]
    assert result["full_name"] == "tumor protein p53"


def test_fetch_gene_summary_unknown_gene(monkeypatch):
    """존재하지 않는 유전자는 None 반환."""
    from scripts.tools.sources.ncbi_gene import fetch_gene_summary

    mock_search = {"esearchresult": {"idlist": []}}
    monkeypatch.setattr(
        "scripts.tools.sources.ncbi_gene.fetch_with_retry",
        lambda url, **kw: mock_search,
    )

    result = fetch_gene_summary("FAKEGENE_XYZ")
    assert result is None


def test_fetch_gene_summary_api_failure(monkeypatch):
    """API 실패 시 None 반환 (graceful)."""
    from scripts.tools.sources.ncbi_gene import fetch_gene_summary

    monkeypatch.setattr(
        "scripts.tools.sources.ncbi_gene.fetch_with_retry",
        lambda url, **kw: None,
    )

    result = fetch_gene_summary("TP53")
    assert result is None
