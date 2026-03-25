import json
import pytest

def test_fetch_cpic_gene_info(monkeypatch):
    """CPIC API에서 유전자 정보를 가져올 수 있어야 함."""
    from scripts.tools.build_gene_knowledge import fetch_cpic_gene

    mock_response = [{
        "symbol": "CYP2C19",
        "name": "cytochrome P450 family 2 subfamily C member 19",
        "cpicPgxGene": True,
        "url": "https://cpicpgx.org/genes-drugs/cyp2c19/",
        "functionExampleSubstratesDrugs": "clopidogrel, voriconazole",
    }]
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_with_retry",
        lambda url, **kw: mock_response,
    )
    # Also mock the guidelines call
    mock_guidelines = [{
        "drugname": "clopidogrel",
        "cpicLevel": "A",
        "pgkbLevel": "1A",
        "guideline": {"name": "CPIC CYP2C19 Clopidogrel 2022", "url": "https://..."},
    }]
    original_fetch = None
    call_count = [0]
    def mock_fetch(url, **kw):
        call_count[0] += 1
        if "gene?" in url:
            return mock_response
        if "pair?" in url:
            return mock_guidelines
        return None
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_with_retry",
        mock_fetch,
    )

    result = fetch_cpic_gene("CYP2C19")
    assert result is not None
    assert result["gene"] == "CYP2C19"
    assert result["content_status"] == "curated-cpic"
    assert any("CPIC" in r.get("source", "") for r in result["references"])

def test_build_gene_knowledge_merges(tmp_path, monkeypatch):
    """기존 gene_knowledge.json과 CPIC 데이터를 병합."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    existing = {
        "genes": [
            {"gene": "TP53", "finding_summary": "existing", "content_status": "ai-generated"},
            {"gene": "CYP2C19", "finding_summary": "old", "content_status": "ai-generated"},
        ]
    }
    existing_path = tmp_path / "gene_knowledge.json"
    existing_path.write_text(json.dumps(existing))

    # Mock CPIC to return data for CYP2C19 only
    def mock_cpic(gene):
        if gene == "CYP2C19":
            return {
                "gene": "CYP2C19",
                "full_name": "cytochrome P450 2C19",
                "finding_summary": "CPIC curated",
                "content_status": "curated-cpic",
                "references": [{"pmid": "34216116", "source": "CPIC 2022"}],
            }
        return None

    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_cpic_gene", mock_cpic)

    output_path = tmp_path / "output.json"
    build_knowledge(str(existing_path), str(output_path), cpic_genes=["CYP2C19"])

    result = json.loads(output_path.read_text())
    genes = {g["gene"]: g for g in result["genes"]}
    # CYP2C19 updated to curated
    assert genes["CYP2C19"]["content_status"] == "curated-cpic"
    assert genes["CYP2C19"]["finding_summary"] == "CPIC curated"
    # TP53 unchanged
    assert genes["TP53"]["content_status"] == "ai-generated"
