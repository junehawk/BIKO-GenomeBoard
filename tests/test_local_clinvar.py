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
    from scripts.storage.build_test_clinvar_db import build_test_db

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
    import scripts.storage.query_local_clinvar as qmod

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
    import scripts.storage.query_local_clinvar as qmod
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
    import scripts.storage.query_local_clinvar as qmod
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
    import scripts.storage.query_local_clinvar as qmod
    from scripts.common.config import reset

    reset()
    qmod._conn = sqlite3.connect(test_db_path)
    qmod._conn.row_factory = sqlite3.Row

    variant = Variant(chrom="chr1", pos=999999999, ref="A", alt="T", rsid=None)
    result = qmod.query_local_clinvar(variant)

    assert result["clinvar_significance"] == "Not Found"
    assert result["acmg_codes"] == []


def _build_single_row_clinvar_db(db_path: str, *, chrom, pos, ref, alt, sig, review, gene="G"):
    """Build a 1-row ClinVar test DB (schema mirrors build_test_clinvar_db)."""
    conn = sqlite3.connect(db_path)
    conn.execute(
        """CREATE TABLE variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, rsid TEXT, gene TEXT,
            clinical_significance TEXT, review_status TEXT, phenotype_list TEXT,
            variation_id TEXT, allele_id TEXT, origin TEXT, assembly TEXT,
            last_evaluated TEXT, number_submitters INTEGER
        )"""
    )
    conn.execute(
        "INSERT INTO variants (chrom,pos,ref,alt,rsid,gene,clinical_significance,review_status,variation_id) "
        "VALUES (?,?,?,?,?,?,?,?,?)",
        (chrom, pos, ref, alt, None, gene, sig, review, "12345"),
    )
    conn.commit()
    conn.close()


def test_relaxed_position_match_does_not_attribute_foreign_allele_significance(tmp_path):
    """A Strategy-3 pos-only hit is a DIFFERENT allele at the same coordinate; its verdict
    must NOT be reported as this variant's classification-driving significance, so it can
    never feed apply_clinvar_override / conflict reconciliation."""
    import scripts.storage.query_local_clinvar as qmod
    from scripts.common.config import reset

    reset()
    db_path = str(tmp_path / "cv.sqlite3")
    _build_single_row_clinvar_db(
        db_path, chrom="chr7", pos=140753336, ref="A", alt="T", sig="Pathogenic", review="reviewed by expert panel"
    )
    qmod._conn = sqlite3.connect(db_path)
    qmod._conn.row_factory = sqlite3.Row

    # Query a DIFFERENT allele (A>C) at the same position — exact ref/alt misses, pos-only hits.
    variant = Variant(chrom="chr7", pos=140753336, ref="A", alt="C", gene="BRAF", rsid=None)
    result = qmod.query_local_clinvar(variant)

    assert result["clinvar_significance"] == "Not Found"
    assert result["acmg_codes"] == []
    assert result["relaxed_match"] is True
    assert result["relaxed_significance"] == "Pathogenic"


def test_exact_match_is_not_flagged_relaxed(tmp_path):
    """An exact ref/alt match (Strategy 2) is a true hit — relaxed_match must be falsy and
    the significance is attributed normally."""
    import scripts.storage.query_local_clinvar as qmod
    from scripts.common.config import reset

    reset()
    db_path = str(tmp_path / "cv.sqlite3")
    _build_single_row_clinvar_db(
        db_path, chrom="chr7", pos=140753336, ref="A", alt="T", sig="Pathogenic", review="reviewed by expert panel"
    )
    qmod._conn = sqlite3.connect(db_path)
    qmod._conn.row_factory = sqlite3.Row

    variant = Variant(chrom="chr7", pos=140753336, ref="A", alt="T", gene="BRAF", rsid=None)
    result = qmod.query_local_clinvar(variant)

    assert result["clinvar_significance"] == "Pathogenic"
    assert not result.get("relaxed_match")
    assert result["api_available"] is True  # DB is available, just no hit


# ---------------------------------------------------------------------------
# 5. _derive_acmg_codes — pathogenic hit → PP5 only (v2.7 CLAS-04)
# ---------------------------------------------------------------------------


def test_derive_acmg_codes_expert_panel():
    """Expert-panel Pathogenic → PP5 annotation only (no generic PS1).

    PS1 specifically means "same amino-acid change as a known pathogenic
    variant", which this gene+position lookup does not verify, and it collides
    with the engine's self-computed PM5. The classification weight of the
    expert-panel verdict is applied independently by apply_clinvar_override.
    """
    from scripts.storage.query_local_clinvar import _derive_acmg_codes

    codes = _derive_acmg_codes("Pathogenic", "reviewed by expert panel")
    assert codes == ["PP5"]
    assert "PS1" not in codes


# ---------------------------------------------------------------------------
# 6. _derive_acmg_codes — single submitter → PP5 only
# ---------------------------------------------------------------------------


def test_derive_acmg_codes_single_submitter():
    from scripts.storage.query_local_clinvar import _derive_acmg_codes

    codes = _derive_acmg_codes("Pathogenic", "criteria provided, single submitter")
    assert "PP5" in codes
    assert "PS1" not in codes


def test_derive_acmg_codes_multiple_submitters_no_ps1():
    """Multiple-submitter Pathogenic → PP5 only (no PS1).

    v2.7 CLAS-04: this tier previously emitted PS1; the multiple-submitter
    classification weight is now applied via apply_clinvar_override instead.
    """
    from scripts.storage.query_local_clinvar import _derive_acmg_codes

    codes = _derive_acmg_codes("Pathogenic", "criteria provided, multiple submitters, no conflicts")
    assert codes == ["PP5"]


def test_derive_acmg_codes_conflicting_emits_nothing():
    """A conflicting ClinVar entry derives no ACMG code (unchanged)."""
    from scripts.storage.query_local_clinvar import _derive_acmg_codes

    codes = _derive_acmg_codes("Conflicting classifications of pathogenicity", "multiple submitters")
    assert codes == []


# ---------------------------------------------------------------------------
# 7. get_db_version
# ---------------------------------------------------------------------------


def test_get_db_version(test_db_path, monkeypatch):
    import scripts.storage.query_local_clinvar as qmod
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
            "kova_db": "data/db/kova.sqlite3",
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
    import scripts.storage.query_local_clinvar as qmod

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
    import scripts.storage.query_local_clinvar as qmod
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
