"""Tests for OMIM genemap2 build and query modules."""

from pathlib import Path

import pytest

from scripts.db.build_omim_genemap_db import _normalise_inheritance, _parse_phenotypes, build_db
from scripts.db.query_omim_genemap import get_gene_phenotypes, get_inheritance_patterns

SAMPLE_FILE = Path(__file__).resolve().parent.parent / "data" / "db" / "sample_genemap2.txt"


@pytest.fixture
def genemap_db(tmp_path):
    """Build a genemap DB from the sample file into a temp directory."""
    db_path = str(tmp_path / "omim_genemap_test.sqlite3")
    build_db(str(SAMPLE_FILE), db_path)
    return db_path


# ── Build tests ──────────────────────────────────────────────────────────


def test_build_omim_genemap_db(genemap_db):
    """DB should contain records for all genes in the sample file."""
    import sqlite3

    conn = sqlite3.connect(genemap_db)
    count = conn.execute("SELECT COUNT(*) FROM omim_genemap").fetchone()[0]
    gene_count = conn.execute("SELECT COUNT(DISTINCT gene) FROM omim_genemap").fetchone()[0]
    conn.close()
    # Sample has 51 gene lines, most with multiple phenotypes
    assert gene_count == 51
    assert count > 51  # more phenotype rows than gene lines


def test_build_metadata(genemap_db):
    """Metadata table should have build_date, source, record_count."""
    import sqlite3

    conn = sqlite3.connect(genemap_db)
    meta = dict(conn.execute("SELECT key, value FROM metadata").fetchall())
    conn.close()
    assert "build_date" in meta
    assert "record_count" in meta
    assert "source" in meta
    assert int(meta["record_count"]) > 0


# ── Query: gene phenotypes ───────────────────────────────────────────────


def test_query_gene_phenotypes_found(genemap_db):
    """BRCA1 should return phenotypes including breast-ovarian cancer."""
    results = get_gene_phenotypes("BRCA1", db_path=genemap_db)
    assert results is not None
    assert len(results) >= 2
    phenotype_names = [r["phenotype"] for r in results]
    assert any("Breast" in p or "breast" in p for p in phenotype_names)


def test_query_gene_phenotypes_not_found(genemap_db):
    """Nonexistent gene should return None."""
    result = get_gene_phenotypes("NONEXISTENT_GENE_XYZ", db_path=genemap_db)
    assert result is None


def test_query_no_db():
    """Querying a nonexistent DB path should return None gracefully."""
    result = get_gene_phenotypes("BRCA1", db_path="/nonexistent/path/db.sqlite3")
    assert result is None

    result2 = get_inheritance_patterns("BRCA1", db_path="/nonexistent/path/db.sqlite3")
    assert result2 is None


# ── Query: inheritance patterns ──────────────────────────────────────────


def test_query_inheritance_patterns(genemap_db):
    """CFTR should have AR inheritance."""
    patterns = get_inheritance_patterns("CFTR", db_path=genemap_db)
    assert patterns is not None
    assert "AR" in patterns


def test_query_inheritance_ad(genemap_db):
    """PTEN has only AD phenotypes."""
    patterns = get_inheritance_patterns("PTEN", db_path=genemap_db)
    assert patterns is not None
    assert patterns == ["AD"]


def test_query_inheritance_xlinked(genemap_db):
    """DMD should have XLR inheritance."""
    patterns = get_inheritance_patterns("DMD", db_path=genemap_db)
    assert patterns is not None
    assert "XLR" in patterns


def test_query_inheritance_mitochondrial(genemap_db):
    """MT-TF should have MT (mitochondrial) inheritance."""
    patterns = get_inheritance_patterns("MT-TF", db_path=genemap_db)
    assert patterns is not None
    assert "MT" in patterns


def test_query_inheritance_somatic(genemap_db):
    """ABL1 should have SOM (somatic mutation) inheritance."""
    patterns = get_inheritance_patterns("ABL1", db_path=genemap_db)
    assert patterns is not None
    assert "SOM" in patterns


def test_query_inheritance_digenic(genemap_db):
    """CDH23 should have DIG (digenic) inheritance."""
    patterns = get_inheritance_patterns("CDH23", db_path=genemap_db)
    assert patterns is not None
    assert "DIG" in patterns


# ── Multiple phenotypes ──────────────────────────────────────────────────


def test_multiple_phenotypes(genemap_db):
    """BRCA1 should have both breast-ovarian cancer and Fanconi anemia."""
    results = get_gene_phenotypes("BRCA1", db_path=genemap_db)
    assert results is not None
    phenotypes = [r["phenotype"] for r in results]
    assert any("Breast" in p for p in phenotypes)
    assert any("Fanconi" in p for p in phenotypes)


def test_multiple_phenotypes_gba(genemap_db):
    """GBA should have multiple Gaucher types plus Parkinson."""
    results = get_gene_phenotypes("GBA", db_path=genemap_db)
    assert results is not None
    assert len(results) >= 4
    phenotypes = [r["phenotype"] for r in results]
    assert any("Gaucher" in p for p in phenotypes)
    assert any("Parkinson" in p for p in phenotypes)


# ── Dual inheritance ─────────────────────────────────────────────────────


def test_dual_inheritance(genemap_db):
    """BRCA1 should have both AD and AR patterns (breast=AD, Fanconi=AR)."""
    patterns = get_inheritance_patterns("BRCA1", db_path=genemap_db)
    assert patterns is not None
    assert "AD" in patterns
    assert "AR" in patterns


def test_dual_inheritance_combined(genemap_db):
    """RYR1 has entries with 'Autosomal dominant/Autosomal recessive' → both AD and AR."""
    patterns = get_inheritance_patterns("RYR1", db_path=genemap_db)
    assert patterns is not None
    assert "AD" in patterns
    assert "AR" in patterns


def test_dual_inheritance_scn5a(genemap_db):
    """SCN5A has both pure AD entries and an AD/AR combined entry."""
    patterns = get_inheritance_patterns("SCN5A", db_path=genemap_db)
    assert patterns is not None
    assert "AD" in patterns
    assert "AR" in patterns


# ── Phenotype parsing unit tests ─────────────────────────────────────────


def test_parse_phenotypes_basic():
    """Parse a single phenotype entry."""
    result = _parse_phenotypes("{Cystic fibrosis}, 219700 (3), Autosomal recessive")
    assert len(result) == 1
    assert result[0]["phenotype"] == "Cystic fibrosis"
    assert result[0]["phenotype_mim"] == "219700"
    assert result[0]["inheritance"] == "AR"


def test_parse_phenotypes_multiple():
    """Parse semicolon-separated phenotype entries."""
    result = _parse_phenotypes(
        "{Breast-ovarian cancer, familial, 1}, 604370 (3), Autosomal dominant; "
        "{Fanconi anemia, complementation group S}, 617883 (3), Autosomal recessive"
    )
    assert len(result) == 2
    assert result[0]["inheritance"] == "AD"
    assert result[1]["inheritance"] == "AR"


def test_parse_phenotypes_empty():
    """Empty string should return empty list."""
    assert _parse_phenotypes("") == []
    assert _parse_phenotypes("   ") == []


def test_normalise_inheritance_combined():
    """AD/AR combined pattern."""
    assert _normalise_inheritance("Autosomal dominant/Autosomal recessive") == "AD/AR"


def test_normalise_inheritance_single():
    """Standard single patterns."""
    assert _normalise_inheritance("Autosomal dominant") == "AD"
    assert _normalise_inheritance("Autosomal recessive") == "AR"
    assert _normalise_inheritance("X-linked recessive") == "XLR"
    assert _normalise_inheritance("Mitochondrial") == "MT"
    assert _normalise_inheritance("Digenic") == "DIG"
    assert _normalise_inheritance("Somatic mutation") == "SOM"
