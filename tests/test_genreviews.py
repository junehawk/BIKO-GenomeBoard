def test_fetch_genreviews_known_gene(monkeypatch):
    """TP53에 대한 GeneReviews PMID 조회."""
    from scripts.tools.sources.genreviews import fetch_genreviews_info

    mock_search = {"esearchresult": {"idlist": ["20301371"]}}
    mock_summary = {
        "result": {
            "20301371": {
                "uid": "20301371",
                "title": "Li-Fraumeni Syndrome",
                "source": "GeneReviews",
                "bookshelfaccession": "NBK1311",
            }
        }
    }

    call_count = [0]

    def mock_fetch(url, **kw):
        call_count[0] += 1
        if "esearch" in url:
            return mock_search
        if "esummary" in url:
            return mock_summary
        return None

    monkeypatch.setattr("scripts.tools.sources.genreviews.fetch_with_retry", mock_fetch)

    result = fetch_genreviews_info("TP53")
    assert result is not None
    assert result["pmid"] == "20301371"
    assert "Li-Fraumeni" in result["title"]
    assert result["nbk_id"] == "NBK1311"
    assert "ncbi.nlm.nih.gov/books/NBK1311" in result["url"]


def test_fetch_genreviews_no_entry(monkeypatch):
    """GeneReviews에 없는 유전자는 None."""
    from scripts.tools.sources.genreviews import fetch_genreviews_info

    mock_search = {"esearchresult": {"idlist": []}}
    monkeypatch.setattr(
        "scripts.tools.sources.genreviews.fetch_with_retry",
        lambda url, **kw: mock_search,
    )

    result = fetch_genreviews_info("OBSCURE_GENE_123")
    assert result is None


def test_fetch_genreviews_api_failure(monkeypatch):
    """API 실패 시 None."""
    from scripts.tools.sources.genreviews import fetch_genreviews_info

    monkeypatch.setattr(
        "scripts.tools.sources.genreviews.fetch_with_retry",
        lambda url, **kw: None,
    )
    assert fetch_genreviews_info("TP53") is None
