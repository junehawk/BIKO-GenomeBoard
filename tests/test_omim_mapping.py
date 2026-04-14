import sqlite3
import pytest
from pathlib import Path


def _create_sample_mim2gene(path):
    """Minimal mim2gene.txt."""
    Path(path).write_text(
        "# Copyright OMIM\n"
        "# Generated: 2026-03-01\n"
        "# MIM Number\tMIM Entry Type\tEntrez Gene ID\tApproved Gene Symbol (HGNC)\tEnsembl Gene ID\n"
        "191170\tgene\t7157\tTP53\tENSG00000141510\n"
        "600185\tgene\t675\tBRCA2\tENSG00000139618\n"
        "602421\tgene\t1080\tCFTR\tENSG00000001626\n"
        "176876\tgene\t5781\tPTPN11\tENSG00000179295\n"
        "114480\tphenotype\t\t\t\n"
        "190685\tgene/phenotype\t7015\tTERT\tENSG00000164362\n"
    )


def test_build_omim_mapping(tmp_path):
    from scripts.db.build_omim_mapping import build_db

    txt_path = str(tmp_path / "mim2gene.txt")
    db_path = str(tmp_path / "omim_mapping.sqlite3")
    _create_sample_mim2gene(txt_path)
    result = build_db(txt_path, db_path)
    assert result == db_path
    conn = sqlite3.connect(db_path)
    rows = conn.execute("SELECT * FROM mim2gene WHERE gene_symbol = 'TP53'").fetchall()
    assert len(rows) >= 1
    conn.close()


def test_get_mim_for_gene(tmp_omim_db):
    from scripts.db.query_omim_mapping import get_mim_for_gene

    result = get_mim_for_gene("TP53", tmp_omim_db)
    assert result is not None
    assert result["mim_number"] == "191170"
    assert result["entrez_id"] == "7157"


def test_get_mim_unknown_gene(tmp_omim_db):
    from scripts.db.query_omim_mapping import get_mim_for_gene

    assert get_mim_for_gene("FAKEGENE", tmp_omim_db) is None


def test_phenotype_entries_excluded(tmp_omim_db):
    """phenotype-only entries (no gene symbol) are excluded."""
    conn = sqlite3.connect(tmp_omim_db)
    rows = conn.execute("SELECT * FROM mim2gene WHERE mim_number = '114480'").fetchall()
    conn.close()
    assert len(rows) == 0


@pytest.fixture
def tmp_omim_db(tmp_path):
    from scripts.db.build_omim_mapping import build_db

    txt_path = str(tmp_path / "mim2gene.txt")
    db_path = str(tmp_path / "omim_mapping.sqlite3")
    _create_sample_mim2gene(txt_path)
    build_db(txt_path, db_path)
    return db_path
