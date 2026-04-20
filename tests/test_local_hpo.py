# tests/test_local_hpo.py
"""Tests for the local HPO gene-phenotype SQLite DB."""

import os
import sqlite3
import tempfile

import pytest


def _create_sample_hpo_tsv(path):
    """Create minimal genes_to_phenotype.txt for testing."""
    with open(path, "w") as f:
        f.write("#Format: ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n")
        f.write("7157\tTP53\tHP:0001250\tSeizure\t-\tOMIM:151623\n")
        f.write("7157\tTP53\tHP:0002664\tNeoplasm\t-\tOMIM:151623\n")
        f.write("672\tBRCA2\tHP:0002664\tNeoplasm\t-\tOMIM:612555\n")
        f.write("1080\tCFTR\tHP:0002110\tBronchiectasis\t-\tOMIM:219700\n")


@pytest.fixture
def tmp_hpo_db(tmp_path):
    """Build temp HPO DB for testing."""
    from scripts.storage.build_hpo_db import build_db

    tsv_path = tmp_path / "genes_to_phenotype.txt"
    with open(tsv_path, "w") as f:
        f.write("#Format: ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n")
        f.write("7157\tTP53\tHP:0001250\tSeizure\t-\tOMIM:151623\n")
        f.write("7157\tTP53\tHP:0002664\tNeoplasm\t-\tOMIM:151623\n")
        f.write("672\tBRCA2\tHP:0002664\tNeoplasm\t-\tOMIM:612555\n")
        f.write("1080\tCFTR\tHP:0002110\tBronchiectasis\t-\tOMIM:219700\n")
    db_path = str(tmp_path / "hpo.sqlite3")
    build_db(str(tsv_path), db_path)
    return db_path


def test_build_hpo_db():
    from scripts.storage.build_hpo_db import build_db

    with tempfile.TemporaryDirectory() as tmpdir:
        tsv_path = os.path.join(tmpdir, "genes_to_phenotype.txt")
        db_path = os.path.join(tmpdir, "hpo.sqlite3")
        _create_sample_hpo_tsv(tsv_path)
        result = build_db(tsv_path, db_path)
        assert result == db_path
        conn = sqlite3.connect(db_path)
        # TP53 has 2 HPO terms
        rows = conn.execute("SELECT hpo_id, hpo_name FROM gene_phenotype WHERE gene_symbol = 'TP53'").fetchall()
        assert len(rows) == 2
        hpo_ids = {r[0] for r in rows}
        assert "HP:0001250" in hpo_ids
        assert "HP:0002664" in hpo_ids
        # HP:0002664 maps to 2 genes
        rows = conn.execute("SELECT gene_symbol FROM gene_phenotype WHERE hpo_id = 'HP:0002664'").fetchall()
        genes = {r[0] for r in rows}
        assert genes == {"TP53", "BRCA2"}
        conn.close()


def test_get_genes_for_hpo(tmp_hpo_db):
    """HPO ID로 연관 유전자 목록 조회."""
    from scripts.storage.query_local_hpo import get_genes_for_hpo

    genes = get_genes_for_hpo("HP:0002664", tmp_hpo_db)
    assert set(genes) == {"TP53", "BRCA2"}


def test_get_hpo_for_gene(tmp_hpo_db):
    """유전자로 연관 HPO 목록 조회."""
    from scripts.storage.query_local_hpo import get_hpo_for_gene

    terms = get_hpo_for_gene("TP53", tmp_hpo_db)
    assert len(terms) == 2
    ids = {t["hpo_id"] for t in terms}
    assert ids == {"HP:0001250", "HP:0002664"}


def test_resolve_hpo_terms_local(tmp_hpo_db):
    """resolve_hpo_terms_local: HPO ID 리스트 → term name + genes (API 없이)."""
    from scripts.storage.query_local_hpo import resolve_hpo_terms_local

    results = resolve_hpo_terms_local(["HP:0001250", "HP:0002664"], tmp_hpo_db)
    assert len(results) == 2
    seizure = next(r for r in results if r["id"] == "HP:0001250")
    assert seizure["name"] == "Seizure"
    assert "TP53" in seizure["genes"]
