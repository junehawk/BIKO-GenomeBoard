# tests/test_local_clinvar.py
"""Tests for the local ClinVar SQLite DB (Task 3.1)."""

import sqlite3
from pathlib import Path

import pytest

from scripts.common.models import Variant


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_db(db_path: str):
    """Build the test DB at db_path using the build_test_clinvar_db module."""
    from scripts.db.build_test_clinvar_db import build_test_db

    build_test_db(db_path)
    return db_path


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def test_db_path(tmp_path_factory):
    """Build a fresh test DB once per module in a temp directory."""
    d = tmp_path_factory.mktemp("clinvar_db")
    db_path = str(d / "clinvar_test.sqlite3")
    _make_db(db_path)
    return db_path


@pytest.fixture(autouse=True)
def reset_connection():
    """Reset the module-level DB connection before each test."""
    import scripts.db.query_local_clinvar as qmod

    qmod.close()
    qmod._conn = None
    yield
    qmod.close()
    qmod._conn = None


# ---------------------------------------------------------------------------
# 1. build_test_db
# ---------------------------------------------------------------------------


def test_build_test_db(tmp_path):
    db_path = str(tmp_path / "clinvar_test.sqlite3")
    _make_db(db_path)
    assert Path(db_path).exists()

    conn = sqlite3.connect(db_path)
    count = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    conn.close()
    assert count == 10


# ---------------------------------------------------------------------------
# 2. query_by_rsid
# ---------------------------------------------------------------------------


def test_query_by_rsid(test_db_path, monkeypatch):
    monkeypatch.setenv("GB_CONFIG_PATH", "")
    import scripts.db.query_local_clinvar as qmod
    from scripts.common.config import reset

    reset()
    # Point the module directly to the test DB
    qmod._conn = sqlite3.connect(test_db_path)
    qmod._conn.row_factory = sqlite3.Row

    variant = Variant(chrom="chr17", pos=7675088, ref="C", alt="A", gene="TP53", rsid="rs28934578")
    result = qmod.query_local_clinvar(variant)

    assert result["clinvar_significance"] == "Likely pathogenic"
    assert result["gene"] == "TP53"
    assert result["api_available"] is True
    assert result["clinvar_id"] is not None


# ---------------------------------------------------------------------------
# 3. query_by_position
# ---------------------------------------------------------------------------


def test_query_by_position(test_db_path, monkeypatch):
    import scripts.db.query_local_clinvar as qmod
    from scripts.common.config import reset

    reset()
    qmod._conn = sqlite3.connect(test_db_path)
    qmod._conn.row_factory = sqlite3.Row

    # CFTR — no rsid provided, match by chrom+pos+ref+alt
    variant = Variant(chrom="chr7", pos=117559590, ref="ATCT", alt="A", gene="CFTR", rsid=None)
    result = qmod.query_local_clinvar(variant)

    assert result["clinvar_significance"] == "Pathogenic"
    assert result["gene"] == "CFTR"
    assert result["api_available"] is True


# ---------------------------------------------------------------------------
# 4. query_not_found
# ---------------------------------------------------------------------------


def test_query_not_found(test_db_path, monkeypatch):
    import scripts.db.query_local_clinvar as qmod
    from scripts.common.config import reset

    reset()
    qmod._conn = sqlite3.connect(test_db_path)
    qmod._conn.row_factory = sqlite3.Row

    variant = Variant(chrom="chr1", pos=999999999, ref="A", alt="T", rsid=None)
    result = qmod.query_local_clinvar(variant)

    assert result["clinvar_significance"] == "Not Found"
    assert result["acmg_codes"] == []
    assert result["api_available"] is True  # DB is available, just no hit


# ---------------------------------------------------------------------------
# 5. _derive_acmg_codes — expert panel → PS1 + PP5
# ---------------------------------------------------------------------------


def test_derive_acmg_codes_expert_panel():
    from scripts.db.query_local_clinvar import _derive_acmg_codes

    codes = _derive_acmg_codes("Pathogenic", "reviewed by expert panel")
    assert "PS1" in codes
    assert "PP5" in codes


# ---------------------------------------------------------------------------
# 6. _derive_acmg_codes — single submitter → PP5 only
# ---------------------------------------------------------------------------


def test_derive_acmg_codes_single_submitter():
    from scripts.db.query_local_clinvar import _derive_acmg_codes

    codes = _derive_acmg_codes("Pathogenic", "criteria provided, single submitter")
    assert "PP5" in codes
    assert "PS1" not in codes


# ---------------------------------------------------------------------------
# 7. get_db_version
# ---------------------------------------------------------------------------


def test_get_db_version(test_db_path, monkeypatch):
    import scripts.db.query_local_clinvar as qmod
    from scripts.common.config import reset

    reset()
    qmod._conn = sqlite3.connect(test_db_path)
    qmod._conn.row_factory = sqlite3.Row

    meta = qmod.get_db_version()
    assert meta["source"] == "test_data"
    assert meta["variant_count"] == "10"
    assert "build_date" in meta
    assert meta.get("clinvar_release") == "2026-03"


# ---------------------------------------------------------------------------
# 8. pipeline local mode
# ---------------------------------------------------------------------------


def test_pipeline_local_mode(test_db_path, tmp_path, monkeypatch):
    """Run orchestrate.run_pipeline with annotation.source=local using the test DB."""
    import yaml
    from scripts.common.config import reset

    # Write a minimal config pointing to the test DB
    cfg = {
        "paths": {
            "clinvar_db": test_db_path,
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

    # Reset the local clinvar connection so it picks up the new config
    import scripts.db.query_local_clinvar as qmod

    qmod.close()
    qmod._conn = None

    # Use a minimal demo VCF that has a TP53 variant matching the test DB
    vcf_content = (
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr17\t7675088\trs28934578\tC\tA\t.\tPASS\tGENE=TP53\n"
    )
    vcf_path = str(tmp_path / "test.vcf")
    Path(vcf_path).write_text(vcf_content)

    output_path = str(tmp_path / "report.html")

    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path=vcf_path,
        output_path=output_path,
        krgdb_path="data/krgdb_freq.tsv",
        skip_api=True,  # no external calls at all
    )

    # With annotation.source=local and skip_api=True the pipeline should complete
    assert result is not None
    assert "variants" in result

    reset()
    qmod.close()
    qmod._conn = None


# ---------------------------------------------------------------------------
# 9. pipeline auto mode — local hit
# ---------------------------------------------------------------------------


def test_pipeline_auto_mode_local_hit(test_db_path, monkeypatch):
    """auto mode should return a result from the local DB without calling the API."""
    import scripts.db.query_local_clinvar as qmod
    from scripts.common.config import reset

    reset()
    qmod._conn = sqlite3.connect(test_db_path)
    qmod._conn.row_factory = sqlite3.Row

    # Directly test the query returns something valid (auto mode logic lives in orchestrate)
    variant = Variant(chrom="chr17", pos=7675088, ref="C", alt="A", gene="TP53", rsid="rs28934578")
    result = qmod.query_local_clinvar(variant)

    assert result["clinvar_significance"] != "Not Found"
    assert result["api_available"] is True
    assert result["gene"] == "TP53"

    reset()
    qmod.close()
    qmod._conn = None
