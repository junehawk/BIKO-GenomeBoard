import sqlite3
import tempfile

import pytest


def _create_sample_clingen_csv(path):
    with open(path, "w") as f:
        f.write("CLINGEN GENE VALIDITY CURATIONS\n")
        f.write("FILE CREATED: 2026-03-01\n")
        f.write("WEBPAGE: https://search.clinicalgenome.org\n")
        f.write(
            '"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","SOP","CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"\n'
        )
        f.write(
            '"TP53","HGNC:11998","Li-Fraumeni syndrome","MONDO:0007903","SOP8","Definitive","https://...","2024-01-15","Cancer"\n'
        )
        f.write(
            '"BRCA2","HGNC:1101","Hereditary breast cancer","MONDO:0016419","SOP8","Definitive","https://...","2024-02-10","Cancer"\n'
        )
        f.write(
            '"SCN1A","HGNC:10585","Dravet syndrome","MONDO:0011073","SOP8","Definitive","https://...","2024-03-01","Epilepsy"\n'
        )
        f.write(
            '"FAKEGENE","HGNC:99999","Some disease","MONDO:0099999","SOP8","Disputed","https://...","2024-01-01","Other"\n'
        )


def test_build_clingen_db():
    from scripts.storage.build_clingen_db import build_db

    with tempfile.TemporaryDirectory() as tmpdir:
        import os

        csv_path = os.path.join(tmpdir, "clingen.csv")
        db_path = os.path.join(tmpdir, "clingen.sqlite3")
        _create_sample_clingen_csv(csv_path)
        result = build_db(csv_path, db_path)
        assert result == db_path
        conn = sqlite3.connect(db_path)
        rows = conn.execute("SELECT * FROM gene_validity WHERE gene_symbol = 'TP53'").fetchall()
        assert len(rows) >= 1
        conn.close()


def test_query_gene_validity_local(tmp_clingen_db):
    from scripts.storage.query_local_clingen import get_gene_validity_local

    assert get_gene_validity_local("TP53", tmp_clingen_db) == "Definitive"
    assert get_gene_validity_local("SCN1A", tmp_clingen_db) == "Definitive"
    assert get_gene_validity_local("FAKEGENE", tmp_clingen_db) == "Disputed"
    assert get_gene_validity_local("UNKNOWN", tmp_clingen_db) is None


def test_get_gene_disease_pairs(tmp_clingen_db):
    from scripts.storage.query_local_clingen import get_gene_disease_pairs

    pairs = get_gene_disease_pairs("TP53", tmp_clingen_db)
    assert len(pairs) >= 1
    assert pairs[0]["disease"] == "Li-Fraumeni syndrome"
    assert pairs[0]["classification"] == "Definitive"


@pytest.fixture
def tmp_clingen_db(tmp_path):
    from scripts.storage.build_clingen_db import build_db

    csv_path = tmp_path / "clingen.csv"
    _create_sample_clingen_csv(str(csv_path))
    db_path = str(tmp_path / "clingen.sqlite3")
    build_db(str(csv_path), db_path)
    return db_path


# ── T3 regression: unavailable-DB hardening (log-once, silent per variant) ──


def test_missing_clingen_db_returns_none_and_logs_once(tmp_path, caplog):
    """If `paths.clingen_db` points at a nonexistent file, every query
    must return None without raising, and the user must see exactly ONE
    warning per run (not one per variant)."""
    from scripts.storage.query_local_clingen import (
        get_gene_validity_local,
        reset_availability_cache,
    )

    reset_availability_cache()
    missing_path = str(tmp_path / "definitely_not_built.sqlite3")

    with caplog.at_level("WARNING"):
        assert get_gene_validity_local("TP53", missing_path) is None
        assert get_gene_validity_local("BRCA1", missing_path) is None
        assert get_gene_validity_local("KCNQ3", missing_path) is None
        assert get_gene_validity_local("ARID1B", missing_path) is None

    not_found_warnings = [r for r in caplog.records if "ClinGen local DB not found" in r.getMessage()]
    assert len(not_found_warnings) == 1


def test_empty_shell_clingen_db_returns_none_and_logs_once(tmp_path, caplog):
    """If the DB file exists but has no `gene_validity` table (the
    setup_databases.sh shell-without-CSV case), queries must still return
    None without raising, and emit exactly ONE warning per run."""
    import sqlite3 as _sq

    from scripts.storage.query_local_clingen import (
        get_gene_disease_pairs,
        get_gene_validity_local,
        reset_availability_cache,
    )

    reset_availability_cache()
    empty_path = str(tmp_path / "empty_shell.sqlite3")
    _sq.connect(empty_path).close()  # create zero-table shell

    with caplog.at_level("WARNING"):
        assert get_gene_validity_local("TP53", empty_path) is None
        assert get_gene_disease_pairs("BRCA1", empty_path) == []
        assert get_gene_validity_local("ARID1B", empty_path) is None

    shell_warnings = [r for r in caplog.records if "has no `gene_validity` table" in r.getMessage()]
    assert len(shell_warnings) == 1


def test_availability_cache_reset_refires_warning(tmp_path, caplog):
    """reset_availability_cache() must re-enable the one-shot warning so
    tests that rotate through multiple DB states can observe each
    transition."""
    from scripts.storage.query_local_clingen import (
        get_gene_validity_local,
        reset_availability_cache,
    )

    reset_availability_cache()
    missing = str(tmp_path / "never_built.sqlite3")

    with caplog.at_level("WARNING"):
        get_gene_validity_local("TP53", missing)
    first_count = sum(1 for r in caplog.records if "ClinGen local DB not found" in r.getMessage())
    assert first_count == 1

    reset_availability_cache()
    with caplog.at_level("WARNING"):
        get_gene_validity_local("TP53", missing)
    second_count = sum(1 for r in caplog.records if "ClinGen local DB not found" in r.getMessage())
    assert second_count == 2
