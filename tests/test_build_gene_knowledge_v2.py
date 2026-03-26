import json
import pytest


def test_build_from_civic_local(tmp_path, monkeypatch):
    """CIViC local DB에서 gene description을 가져와 knowledge 생성."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    # Mock CIViC gene summary
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_gene_summary",
        lambda gene: {"gene": gene, "description": f"{gene} is a key cancer gene.", "aliases": ""}
        if gene == "BRAF" else None,
    )
    # Mock NCBI Gene (fallback)
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_gene_summary",
        lambda gene: None,
    )
    # Mock GeneReviews
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_genreviews_info",
        lambda gene: None,
    )
    # Mock ClinGen
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_gene_validity_local",
        lambda gene, **kw: None,
    )
    # Mock CIViC treatment
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_treatment_summary",
        lambda gene, variant=None: "Level A: Vemurafenib — Sensitivity in Melanoma",
    )
    # Mock CIViC evidence for refs
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_variant_evidence",
        lambda gene, variant=None: [{"pmid": "20979469", "citation": "Chapman 2011",
         "evidence_type": "Predictive", "significance": "Sensitivity"}],
    )

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["BRAF"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    genes = {g["gene"]: g for g in data["genes"]}
    assert "BRAF" in genes
    assert genes["BRAF"]["content_status"] == "curated-civic"
    assert "cancer gene" in genes["BRAF"]["finding_summary"]
    assert "Vemurafenib" in genes["BRAF"]["treatment_strategies"]


def test_build_ncbi_fallback(tmp_path, monkeypatch):
    """CIViC에 없는 유전자는 NCBI Gene fallback."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_gene_summary",
        lambda gene: None,
    )
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_gene_summary",
        lambda gene: {"gene": gene, "entrez_id": "999", "full_name": "novel gene",
                      "summary": "This gene is involved in cell signaling.", "aliases": ""}
        if gene == "NOVELGENE" else None,
    )
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_genreviews_info",
        lambda gene: None,
    )
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_gene_validity_local",
        lambda gene, **kw: None,
    )
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_treatment_summary",
        lambda gene, variant=None: "",
    )
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_variant_evidence",
        lambda gene, variant=None: [],
    )

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["NOVELGENE"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    genes = {g["gene"]: g for g in data["genes"]}
    assert genes["NOVELGENE"]["content_status"] == "curated-ncbi"
    assert "cell signaling" in genes["NOVELGENE"]["finding_summary"]


def test_build_cpic_priority(tmp_path, monkeypatch):
    """PGx 유전자는 CPIC가 우선."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_cpic_gene",
        lambda gene: {"gene": gene, "full_name": "CYP2C19 enzyme",
                      "finding_summary": "CPIC curated", "content_status": "curated-cpic",
                      "references": [{"pmid": "34216116", "source": "CPIC"}],
                      "treatment_strategies": "CPIC guidelines", "frequency_prognosis": "",
                      "function_summary": "", "clinical_significance": "",
                      "associated_conditions": [], "korean_specific_note": None, "hgvs": {}}
        if gene == "CYP2C19" else None,
    )
    # Mock others to return data (should be ignored for PGx genes)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_genreviews_info", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_validity_local", lambda g, **kw: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_treatment_summary", lambda g, v=None: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_variant_evidence", lambda g, v=None: [])

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["CYP2C19"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    genes = {g["gene"]: g for g in data["genes"]}
    assert genes["CYP2C19"]["content_status"] == "curated-cpic"


def test_build_minimal_fallback(tmp_path, monkeypatch):
    """모든 소스 실패 시 auto-minimal entry 생성."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_genreviews_info", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_validity_local", lambda g, **kw: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_cpic_gene", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_treatment_summary", lambda g, v=None: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_variant_evidence", lambda g, v=None: [])

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["UNKNOWN_GENE"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    genes = {g["gene"]: g for g in data["genes"]}
    assert genes["UNKNOWN_GENE"]["content_status"] == "auto-minimal"


def test_build_genreviews_adds_reference(tmp_path, monkeypatch):
    """GeneReviews PMID가 references에 추가됨."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_summary", lambda g: None)
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_gene_summary",
        lambda g: {"gene": g, "entrez_id": "7157", "full_name": "tumor protein p53",
                    "summary": "Encodes a tumor suppressor.", "aliases": ""} if g == "TP53" else None,
    )
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_genreviews_info",
        lambda g: {"gene": g, "pmid": "20301371", "title": "Li-Fraumeni Syndrome",
                    "nbk_id": "NBK1311", "url": "https://www.ncbi.nlm.nih.gov/books/NBK1311/"}
        if g == "TP53" else None,
    )
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_validity_local", lambda g, **kw: "Definitive" if g == "TP53" else None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_cpic_gene", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_treatment_summary", lambda g, v=None: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_variant_evidence", lambda g, v=None: [])

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["TP53"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    tp53 = {g["gene"]: g for g in data["genes"]}["TP53"]
    pmids = [r.get("pmid") for r in tp53.get("references", [])]
    assert "20301371" in pmids
    assert "Definitive" in tp53["finding_summary"]
