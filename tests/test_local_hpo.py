# tests/test_local_hpo.py
"""Tests for the local HPO gene-phenotype SQLite DB."""
import sqlite3
import tempfile
import os
import pytest
from pathlib import Path


def _create_sample_hpo_tsv(path):
    """Create minimal genes_to_phenotype.txt for testing."""
    with open(path, "w") as f:
        f.write("#Format: ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n")
        f.write("7157\tTP53\tHP:0001250\tSeizure\t-\tOMIM:151623\n")
        f.write("7157\tTP53\tHP:0002664\tNeoplasm\t-\tOMIM:151623\n")
        f.write("672\tBRCA2\tHP:0002664\tNeoplasm\t-\tOMIM:612555\n")
        f.write("1080\tCFTR\tHP:0002110\tBronchiectasis\t-\tOMIM:219700\n")


def test_build_hpo_db():
    from scripts.db.build_hpo_db import build_db
    with tempfile.TemporaryDirectory() as tmpdir:
        tsv_path = os.path.join(tmpdir, "genes_to_phenotype.txt")
        db_path = os.path.join(tmpdir, "hpo.sqlite3")
        _create_sample_hpo_tsv(tsv_path)
        result = build_db(tsv_path, db_path)
        assert result == db_path
        conn = sqlite3.connect(db_path)
        # TP53 has 2 HPO terms
        rows = conn.execute(
            "SELECT hpo_id, hpo_name FROM gene_phenotype WHERE gene_symbol = 'TP53'"
        ).fetchall()
        assert len(rows) == 2
        hpo_ids = {r[0] for r in rows}
        assert "HP:0001250" in hpo_ids
        assert "HP:0002664" in hpo_ids
        # HP:0002664 maps to 2 genes
        rows = conn.execute(
            "SELECT gene_symbol FROM gene_phenotype WHERE hpo_id = 'HP:0002664'"
        ).fetchall()
        genes = {r[0] for r in rows}
        assert genes == {"TP53", "BRCA2"}
        conn.close()
