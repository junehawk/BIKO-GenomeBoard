# tests/test_local_gnomad.py
"""Tests for the local gnomAD SQLite DB (Task 3.2)."""
import sqlite3
import tempfile
from pathlib import Path

import pytest

from scripts.common.models import Variant


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_db(db_path: str):
    """Build the test gnomAD DB at db_path."""
    from scripts.db.build_test_gnomad_db import build_test_db
    build_test_db(db_path)
    return db_path


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def test_db_path(tmp_path_factory):
    """Build a fresh test gnomAD DB once per module in a temp directory."""
    d = tmp_path_factory.mktemp("gnomad_db")
    db_path = str(d / "gnomad_test.sqlite3")
    _make_db(db_path)
    return db_path


@pytest.fixture(autouse=True)
def reset_connection():
    """Reset the module-level DB connection before each test."""
    import scripts.db.query_local_gnomad as qmod
    qmod.close()
    qmod._conn = None
    yield
    qmod.close()
    qmod._conn = None


# ---------------------------------------------------------------------------
# 1. test_build_test_db
# ---------------------------------------------------------------------------

def test_build_test_db(tmp_path):
    db_path = str(tmp_path / "gnomad_test.sqlite3")
    _make_db(db_path)
    assert Path(db_path).exists()

    conn = sqlite3.connect(db_path)
    count = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    conn.close()
    assert count == 6


# ---------------------------------------------------------------------------
# 2. test_query_by_position_exact
# ---------------------------------------------------------------------------

def test_query_by_position_exact(test_db_path):
    import scripts.db.query_local_gnomad as qmod
    from scripts.common.config import reset
    reset()
    qmod._conn = sqlite3.connect(test_db_path)
    qmod._conn.row_factory = sqlite3.Row

    # TP53 variant — exact chrom+pos+ref+alt match
    variant = Variant(chrom="chr17", pos=7675088, ref="C", alt="A", gene="TP53", rsid=None)
    result = qmod.query_local_gnomad(variant)

    assert result["api_available"] is True
    assert result["gnomad_all"] == pytest.approx(0.000004, rel=1e-3)
    assert result["gnomad_eas"] == pytest.approx(0.000008, rel=1e-3)


# ---------------------------------------------------------------------------
# 3. test_query_by_rsid
# ---------------------------------------------------------------------------

def test_query_by_rsid(test_db_path):
    import scripts.db.query_local_gnomad as qmod
    from scripts.common.config import reset
    reset()
    qmod._conn = sqlite3.connect(test_db_path)
    qmod._conn.row_factory = sqlite3.Row

    # Use rsID fallback: provide wrong chrom to force rsid path
    variant = Variant(chrom="chr1", pos=999, ref="C", alt="A", gene="TP53", rsid="rs28934578")
    result = qmod.query_local_gnomad(variant)

    assert result["api_available"] is True
    assert result["gnomad_all"] == pytest.approx(0.000004, rel=1e-3)
    assert result["gnomad_eas"] == pytest.approx(0.000008, rel=1e-3)


# ---------------------------------------------------------------------------
# 4. test_query_not_found
# ---------------------------------------------------------------------------

def test_query_not_found(test_db_path):
    import scripts.db.query_local_gnomad as qmod
    from scripts.common.config import reset
    reset()
    qmod._conn = sqlite3.connect(test_db_path)
    qmod._conn.row_factory = sqlite3.Row

    variant = Variant(chrom="chr1", pos=999999999, ref="A", alt="T", rsid=None)
    result = qmod.query_local_gnomad(variant)

    assert result["gnomad_all"] is None
    assert result["gnomad_eas"] is None
    assert result["api_available"] is True  # DB is available, just no hit


# ---------------------------------------------------------------------------
# 5. test_query_returns_eas_af
# ---------------------------------------------------------------------------

def test_query_returns_eas_af(test_db_path):
    """EAS-enriched variant (NUDT15 rs116855232) has higher EAS AF than global."""
    import scripts.db.query_local_gnomad as qmod
    from scripts.common.config import reset
    reset()
    qmod._conn = sqlite3.connect(test_db_path)
    qmod._conn.row_factory = sqlite3.Row

    variant = Variant(chrom="chr13", pos=48045719, ref="C", alt="T", gene="NUDT15", rsid="rs116855232")
    result = qmod.query_local_gnomad(variant)

    assert result["api_available"] is True
    assert result["gnomad_all"] == pytest.approx(0.01257, rel=1e-3)
    # EAS AF should be notably higher than global for NUDT15
    assert result["gnomad_eas"] == pytest.approx(0.101, rel=1e-3)
    assert result["gnomad_eas"] > result["gnomad_all"]


# ---------------------------------------------------------------------------
# 6. test_get_db_version
# ---------------------------------------------------------------------------

def test_get_db_version(test_db_path):
    import scripts.db.query_local_gnomad as qmod
    from scripts.common.config import reset
    reset()
    qmod._conn = sqlite3.connect(test_db_path)
    qmod._conn.row_factory = sqlite3.Row

    meta = qmod.get_db_version()
    assert meta["gnomad_version"] == "4.1"
    assert meta["variant_count"] == "6"
    assert meta["assembly"] == "GRCh38"
    assert "build_date" in meta


# ---------------------------------------------------------------------------
# 7. test_pipeline_local_gnomad
# ---------------------------------------------------------------------------

def test_pipeline_local_gnomad(test_db_path, tmp_path, monkeypatch):
    """Run orchestrate.run_pipeline with annotation.source=local using test gnomAD DB."""
    import yaml
    from scripts.common.config import reset

    cfg = {
        "paths": {
            "gnomad_db": test_db_path,
            "clinvar_db": "data/db/clinvar.sqlite3",
            "krgdb": "data/krgdb_freq.tsv",
            "gene_knowledge": "data/gene_knowledge.json",
            "pgx_table": "data/korean_pgx_table.json",
            "acmg_rules": "data/acmg_rules.json",
            "templates": "templates",
        },
        "annotation": {"source": "local"},
        "thresholds": {"ba1": 0.05, "bs1": 0.01, "pm2": 0.001},
        "report": {"default_mode": "cancer", "default_genome_build": "GRCh38"},
        "pgx": {
            "genes": ["CYP2D6", "CYP2C19", "CYP2C9", "HLA-B", "HLA-A", "NUDT15", "TPMT", "DPYD"],
            "risk_factor_genes": ["APOE"],
        },
    }
    cfg_path = str(tmp_path / "test_config.yaml")
    with open(cfg_path, "w") as f:
        yaml.dump(cfg, f)

    monkeypatch.setenv("GB_CONFIG_PATH", cfg_path)
    reset()

    # Reset gnomAD connection so it picks up the new config
    import scripts.db.query_local_gnomad as qmod
    qmod.close()
    qmod._conn = None

    # Use a CYP2C19 variant that exists in the test gnomAD DB
    vcf_content = (
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr10\t94781859\trs4244285\tG\tA\t.\tPASS\tGENE=CYP2C19\n"
    )
    vcf_path = str(tmp_path / "test.vcf")
    Path(vcf_path).write_text(vcf_content)

    output_path = str(tmp_path / "report.html")

    from scripts.orchestrate import run_pipeline
    result = run_pipeline(
        vcf_path=vcf_path,
        output_path=output_path,
        krgdb_path="data/krgdb_freq.tsv",
        skip_api=True,
    )

    assert result is not None
    assert "variants" in result
    # gnomAD local version should appear in db_versions
    assert "gnomad_local" in result["db_versions"]

    reset()
    qmod.close()
    qmod._conn = None
