# tests/test_build_gene_knowledge_phase2b.py
import json
import pytest


def test_orphanet_prevalence_in_knowledge(tmp_path, monkeypatch):
    """Orphanet prevalence가 frequency_prognosis에 포함."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_summary",
        lambda g: {"gene": g, "description": "Test gene", "aliases": ""} if g == "CFTR" else None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_genreviews_info", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_validity_local", lambda g, **kw: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_treatment_summary", lambda g, v=None: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_variant_evidence", lambda g, v=None: [])
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_cpic_gene", lambda g: None)

    # Mock Orphanet
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_prevalence_text",
        lambda g: "Cystic fibrosis: 1-5 / 10 000 (Europe)" if g == "CFTR" else "")

    # Mock GeneReviews local DB
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_genreviews_for_gene_local",
        lambda g: None)

    # Mock OMIM
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_mim_for_gene",
        lambda g: None)

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["CFTR"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    cftr = {g["gene"]: g for g in data["genes"]}["CFTR"]
    assert "Cystic fibrosis" in cftr["frequency_prognosis"]
    assert "1-5 / 10 000" in cftr["frequency_prognosis"]


def test_genreviews_local_replaces_api(tmp_path, monkeypatch):
    """GeneReviews local DB가 API 대신 references에 추가됨."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_gene_summary",
        lambda g: {"gene": g, "entrez_id": "7157", "full_name": "tumor protein p53",
                    "summary": "Tumor suppressor.", "aliases": ""} if g == "TP53" else None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_genreviews_info", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_validity_local", lambda g, **kw: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_treatment_summary", lambda g, v=None: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_variant_evidence", lambda g, v=None: [])
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_cpic_gene", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_prevalence_text", lambda g: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_mim_for_gene", lambda g: None)

    # GeneReviews local DB returns data
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_genreviews_for_gene_local",
        lambda g: {"gene": g, "nbk_id": "NBK1311", "pmid": "20301371",
                    "title": "Li-Fraumeni Syndrome",
                    "url": "https://www.ncbi.nlm.nih.gov/books/NBK1311/"} if g == "TP53" else None)

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["TP53"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    tp53 = {g["gene"]: g for g in data["genes"]}["TP53"]
    pmids = [r.get("pmid") for r in tp53.get("references", [])]
    assert "20301371" in pmids


def test_omim_mim_in_references(tmp_path, monkeypatch):
    """OMIM MIM number가 references에 포함."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_gene_summary",
        lambda g: {"gene": g, "entrez_id": "7157", "full_name": "tumor protein p53",
                    "summary": "Tumor suppressor.", "aliases": ""} if g == "TP53" else None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_genreviews_info", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_validity_local", lambda g, **kw: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_treatment_summary", lambda g, v=None: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_variant_evidence", lambda g, v=None: [])
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_cpic_gene", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_prevalence_text", lambda g: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_genreviews_for_gene_local", lambda g: None)

    # OMIM returns MIM
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_mim_for_gene",
        lambda g: {"gene": g, "mim_number": "191170", "entry_type": "gene",
                    "entrez_id": "7157", "ensembl_id": "", "url": "https://omim.org/entry/191170"}
        if g == "TP53" else None)

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["TP53"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    tp53 = {g["gene"]: g for g in data["genes"]}["TP53"]
    sources = [r.get("source") for r in tp53.get("references", [])]
    assert any("OMIM" in s for s in sources)
