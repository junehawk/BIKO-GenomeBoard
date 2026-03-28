# tests/test_genreviews_db.py
import sqlite3
import pytest
from pathlib import Path


def _create_sample_genreviews_tsv(path):
    """Minimal GeneReviews identifier file (pipe-delimited, real NCBI format)."""
    # Real format: GR_shortname|NBK_id|gene_symbol|disease_name (no header)
    Path(path).write_text(
        "li-fraumeni|NBK1311|TP53|Li-Fraumeni Syndrome\n"
        "brca1-2|NBK1247|BRCA2|Hereditary Breast and Ovarian Cancer\n"
        "brca1-2|NBK1247|BRCA1|Hereditary Breast and Ovarian Cancer\n"
        "cftr|NBK1250|CFTR|CFTR-Related Disorders\n"
        "noonan|NBK1124|PTPN11|Noonan Syndrome\n"
    )


def _create_sample_titles_tsv(path):
    """Minimal GeneReviews titles TSV."""
    # Format: GR_shortname\tGR_Title\tNBK_id\tPMID
    Path(path).write_text(
        "#GR_shortname\tGR_Title\tNBK_id\tPMID\n"
        "li-fraumeni\tLi-Fraumeni Syndrome\tNBK1311\t20301371\n"
        "brca1-2\tBRCA1- and BRCA2-Associated Hereditary Breast and Ovarian Cancer\tNBK1247\t20301425\n"
        "cftr\tCFTR-Related Disorders\tNBK1250\t20301428\n"
        "noonan\tNoonan Syndrome\tNBK1124\t20301303\n"
    )


def test_build_genreviews_db(tmp_path):
    from scripts.db.build_genreviews_db import build_db
    genes_path = str(tmp_path / "genes.txt")
    titles_path = str(tmp_path / "titles.txt")
    db_path = str(tmp_path / "genreviews.sqlite3")
    _create_sample_genreviews_tsv(genes_path)
    _create_sample_titles_tsv(titles_path)
    result = build_db(genes_path, titles_path, db_path)
    assert result == db_path
    conn = sqlite3.connect(db_path)
    rows = conn.execute("SELECT * FROM genreviews WHERE gene_symbol = 'TP53'").fetchall()
    assert len(rows) >= 1
    conn.close()


def test_query_genreviews_by_gene(tmp_genreviews_db):
    from scripts.db.query_genreviews import get_genreviews_for_gene
    result = get_genreviews_for_gene("TP53", tmp_genreviews_db)
    assert result is not None
    assert result["nbk_id"] == "NBK1311"
    assert result["pmid"] == "20301371"
    assert "Li-Fraumeni" in result["title"]
    assert "ncbi.nlm.nih.gov/books/NBK1311" in result["url"]


def test_query_genreviews_unknown(tmp_genreviews_db):
    from scripts.db.query_genreviews import get_genreviews_for_gene
    assert get_genreviews_for_gene("FAKEGENE", tmp_genreviews_db) is None


def test_query_brca_shared_entry(tmp_genreviews_db):
    """BRCA1과 BRCA2는 같은 GeneReviews entry를 공유."""
    from scripts.db.query_genreviews import get_genreviews_for_gene
    r1 = get_genreviews_for_gene("BRCA1", tmp_genreviews_db)
    r2 = get_genreviews_for_gene("BRCA2", tmp_genreviews_db)
    assert r1["nbk_id"] == r2["nbk_id"] == "NBK1247"


@pytest.fixture
def tmp_genreviews_db(tmp_path):
    from scripts.db.build_genreviews_db import build_db
    genes_path = str(tmp_path / "genes.txt")
    titles_path = str(tmp_path / "titles.txt")
    db_path = str(tmp_path / "genreviews.sqlite3")
    _create_sample_genreviews_tsv(genes_path)
    _create_sample_titles_tsv(titles_path)
    build_db(genes_path, titles_path, db_path)
    return db_path
