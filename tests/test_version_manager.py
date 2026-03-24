# tests/test_version_manager.py
"""Tests for the centralized DB version manager (Task 3.5)."""
import sqlite3
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_clinvar_db(db_path: str):
    from scripts.db.build_test_clinvar_db import build_test_db
    build_test_db(db_path)
    return db_path


def _make_gnomad_db(db_path: str):
    from scripts.db.build_test_gnomad_db import build_test_db
    build_test_db(db_path)
    return db_path


# ---------------------------------------------------------------------------
# 1. test_get_all_db_versions_with_local_dbs
# ---------------------------------------------------------------------------

def test_get_all_db_versions_with_local_dbs(tmp_path, monkeypatch):
    """With both local DBs present, versions report local_db source for ClinVar & gnomAD."""
    import yaml
    from scripts.common.config import reset

    clinvar_db = str(tmp_path / "clinvar_test.sqlite3")
    gnomad_db = str(tmp_path / "gnomad_test.sqlite3")
    _make_clinvar_db(clinvar_db)
    _make_gnomad_db(gnomad_db)

    cfg = {
        "paths": {
            "clinvar_db": clinvar_db,
            "gnomad_db": gnomad_db,
            "gnomad_vcf_dir": str(tmp_path / "empty_gnomad_vcf"),
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
            "genes": ["CYP2D6", "CYP2C19"],
            "risk_factor_genes": ["APOE"],
        },
    }
    cfg_path = str(tmp_path / "test_config.yaml")
    with open(cfg_path, "w") as f:
        yaml.dump(cfg, f)

    monkeypatch.setenv("GB_CONFIG_PATH", cfg_path)
    reset()

    # Reset local DB connections so they pick up new config
    import scripts.db.query_local_clinvar as cv_mod
    import scripts.db.query_local_gnomad as gn_mod
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None

    from scripts.db.version_manager import get_all_db_versions
    versions = get_all_db_versions(skip_api=True)

    assert "ClinVar" in versions
    assert versions["ClinVar"]["source"] == "local_db"
    assert versions["ClinVar"]["release"] == "2026-03"
    assert versions["ClinVar"]["assembly"] == "GRCh38"

    assert "gnomAD" in versions
    assert versions["gnomAD"]["source"] == "local_db"
    assert versions["gnomAD"]["version"] == "4.1"
    assert versions["gnomAD"]["assembly"] == "GRCh38"

    reset()
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None


# ---------------------------------------------------------------------------
# 2. test_get_all_db_versions_api_mode
# ---------------------------------------------------------------------------

def test_get_all_db_versions_api_mode(tmp_path, monkeypatch):
    """When no local DBs are configured and skip_api=False, versions show api source."""
    import yaml
    from scripts.common.config import reset

    cfg = {
        "paths": {
            "clinvar_db": str(tmp_path / "nonexistent_clinvar.sqlite3"),
            "gnomad_db": str(tmp_path / "nonexistent_gnomad.sqlite3"),
            "gnomad_vcf_dir": str(tmp_path / "nonexistent_gnomad_vcf"),
            "krgdb": str(tmp_path / "nonexistent_krgdb.tsv"),
            "gene_knowledge": "data/gene_knowledge.json",
            "pgx_table": "data/korean_pgx_table.json",
            "acmg_rules": "data/acmg_rules.json",
            "templates": "templates",
        },
        "annotation": {"source": "api"},
        "thresholds": {"ba1": 0.05, "bs1": 0.01, "pm2": 0.001},
        "report": {"default_mode": "cancer", "default_genome_build": "GRCh38"},
        "pgx": {"genes": [], "risk_factor_genes": []},
    }
    cfg_path = str(tmp_path / "test_config.yaml")
    with open(cfg_path, "w") as f:
        yaml.dump(cfg, f)

    monkeypatch.setenv("GB_CONFIG_PATH", cfg_path)
    reset()

    import scripts.db.query_local_clinvar as cv_mod
    import scripts.db.query_local_gnomad as gn_mod
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None

    from scripts.db.version_manager import get_all_db_versions
    versions = get_all_db_versions(skip_api=False)

    assert versions["ClinVar"]["source"] == "api"
    assert "E-utilities" in versions["ClinVar"]["release"]
    assert versions["gnomAD"]["source"] == "api"
    assert "v4.1" in versions["gnomAD"]["version"]

    reset()
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None


# ---------------------------------------------------------------------------
# 3. test_get_all_db_versions_skip_api_no_local
# ---------------------------------------------------------------------------

def test_get_all_db_versions_skip_api_no_local(tmp_path, monkeypatch):
    """skip_api=True with no local DBs → not_available source."""
    import yaml
    from scripts.common.config import reset

    cfg = {
        "paths": {
            "clinvar_db": str(tmp_path / "nonexistent_clinvar.sqlite3"),
            "gnomad_db": str(tmp_path / "nonexistent_gnomad.sqlite3"),
            "gnomad_vcf_dir": str(tmp_path / "nonexistent_gnomad_vcf"),
            "krgdb": str(tmp_path / "nonexistent_krgdb.tsv"),
            "gene_knowledge": "data/gene_knowledge.json",
            "pgx_table": "data/korean_pgx_table.json",
            "acmg_rules": "data/acmg_rules.json",
            "templates": "templates",
        },
        "annotation": {"source": "local"},
        "thresholds": {"ba1": 0.05, "bs1": 0.01, "pm2": 0.001},
        "report": {"default_mode": "cancer", "default_genome_build": "GRCh38"},
        "pgx": {"genes": [], "risk_factor_genes": []},
    }
    cfg_path = str(tmp_path / "test_config.yaml")
    with open(cfg_path, "w") as f:
        yaml.dump(cfg, f)

    monkeypatch.setenv("GB_CONFIG_PATH", cfg_path)
    reset()

    import scripts.db.query_local_clinvar as cv_mod
    import scripts.db.query_local_gnomad as gn_mod
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None

    from scripts.db.version_manager import get_all_db_versions
    versions = get_all_db_versions(skip_api=True)

    assert versions["ClinVar"]["source"] == "not_available"
    assert versions["gnomAD"]["source"] == "not_available"

    reset()
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None


# ---------------------------------------------------------------------------
# 4. test_get_all_db_versions_includes_krgdb
# ---------------------------------------------------------------------------

def test_get_all_db_versions_includes_krgdb(tmp_path, monkeypatch):
    """KRGDB file present → versions includes KRGDB with local_file source."""
    import yaml
    from scripts.common.config import reset

    krgdb_file = tmp_path / "krgdb_freq.tsv"
    krgdb_file.write_text("chrom\tpos\tref\talt\tfreq\n")

    cfg = {
        "paths": {
            "clinvar_db": str(tmp_path / "no_clinvar.sqlite3"),
            "gnomad_db": str(tmp_path / "no_gnomad.sqlite3"),
            "krgdb": str(krgdb_file),
            "gene_knowledge": "data/gene_knowledge.json",
            "pgx_table": "data/korean_pgx_table.json",
            "acmg_rules": "data/acmg_rules.json",
            "templates": "templates",
        },
        "annotation": {"source": "local"},
        "thresholds": {"ba1": 0.05, "bs1": 0.01, "pm2": 0.001},
        "report": {"default_mode": "cancer", "default_genome_build": "GRCh38"},
        "pgx": {"genes": [], "risk_factor_genes": []},
    }
    cfg_path = str(tmp_path / "test_config.yaml")
    with open(cfg_path, "w") as f:
        yaml.dump(cfg, f)

    monkeypatch.setenv("GB_CONFIG_PATH", cfg_path)
    reset()

    import scripts.db.query_local_clinvar as cv_mod
    import scripts.db.query_local_gnomad as gn_mod
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None

    from scripts.db.version_manager import get_all_db_versions
    versions = get_all_db_versions(skip_api=True)

    assert "KRGDB" in versions
    assert versions["KRGDB"]["source"] == "local_file"
    assert "modified" in versions["KRGDB"]
    assert "size_bytes" in versions["KRGDB"]

    reset()
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None


# ---------------------------------------------------------------------------
# 5. test_get_all_db_versions_includes_acmg
# ---------------------------------------------------------------------------

def test_get_all_db_versions_includes_acmg(tmp_path, monkeypatch):
    """ACMG entry is always present in versions dict."""
    import yaml
    from scripts.common.config import reset

    cfg = {
        "paths": {
            "clinvar_db": str(tmp_path / "no_clinvar.sqlite3"),
            "gnomad_db": str(tmp_path / "no_gnomad.sqlite3"),
            "krgdb": str(tmp_path / "no_krgdb.tsv"),
            "gene_knowledge": "data/gene_knowledge.json",
            "pgx_table": "data/korean_pgx_table.json",
            "acmg_rules": "data/acmg_rules.json",
            "templates": "templates",
        },
        "annotation": {"source": "local"},
        "thresholds": {"ba1": 0.05, "bs1": 0.01, "pm2": 0.001},
        "report": {"default_mode": "cancer", "default_genome_build": "GRCh38"},
        "pgx": {"genes": ["CYP2D6", "CYP2C19"], "risk_factor_genes": ["APOE"]},
    }
    cfg_path = str(tmp_path / "test_config.yaml")
    with open(cfg_path, "w") as f:
        yaml.dump(cfg, f)

    monkeypatch.setenv("GB_CONFIG_PATH", cfg_path)
    reset()

    import scripts.db.query_local_clinvar as cv_mod
    import scripts.db.query_local_gnomad as gn_mod
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None

    from scripts.db.version_manager import get_all_db_versions
    versions = get_all_db_versions(skip_api=True)

    assert "ACMG" in versions
    assert "ACMG/AMP 2015" in versions["ACMG"]["standard"]
    assert "CPIC/PGx" in versions
    assert versions["CPIC/PGx"]["genes_count"] == 2

    reset()
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None


# ---------------------------------------------------------------------------
# 6. test_report_shows_db_source_type
# ---------------------------------------------------------------------------

def test_report_shows_db_source_type(tmp_path, monkeypatch):
    """Rendered report HTML contains 'Local DB' text when local DBs are used."""
    import yaml
    from scripts.common.config import reset

    clinvar_db = str(tmp_path / "clinvar_test.sqlite3")
    gnomad_db = str(tmp_path / "gnomad_test.sqlite3")
    _make_clinvar_db(clinvar_db)
    _make_gnomad_db(gnomad_db)

    cfg = {
        "paths": {
            "clinvar_db": clinvar_db,
            "gnomad_db": gnomad_db,
            "gnomad_vcf_dir": str(tmp_path / "empty_gnomad_vcf"),
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
            "genes": ["CYP2D6", "CYP2C19"],
            "risk_factor_genes": ["APOE"],
        },
    }
    cfg_path = str(tmp_path / "test_config.yaml")
    with open(cfg_path, "w") as f:
        yaml.dump(cfg, f)

    monkeypatch.setenv("GB_CONFIG_PATH", cfg_path)
    reset()

    import scripts.db.query_local_clinvar as cv_mod
    import scripts.db.query_local_gnomad as gn_mod
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None

    from scripts.db.version_manager import get_all_db_versions
    from scripts.counselor.generate_pdf import generate_report_html

    db_versions = get_all_db_versions(skip_api=True)
    report_data = {
        "sample_id": "TEST001",
        "date": "2026-03-24",
        "variants": [],
        "pgx_results": [],
        "summary": {
            "total": 0, "pathogenic": 0, "likely_pathogenic": 0,
            "drug_response": 0, "risk_factor": 0, "vus": 0,
            "likely_benign": 0, "benign": 0,
        },
        "db_versions": db_versions,
        "pipeline": {"skip_api": True, "krgdb_path": ""},
        "mode": "cancer",
        "hpo_results": [],
    }

    html = generate_report_html(report_data, mode="cancer")
    assert "Local DB" in html

    reset()
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None


# ---------------------------------------------------------------------------
# 7. test_report_shows_build_date
# ---------------------------------------------------------------------------

def test_report_shows_build_date(tmp_path, monkeypatch):
    """Rendered report HTML contains the build date from the local ClinVar DB."""
    import yaml
    from scripts.common.config import reset

    clinvar_db = str(tmp_path / "clinvar_test.sqlite3")
    gnomad_db = str(tmp_path / "gnomad_test.sqlite3")
    _make_clinvar_db(clinvar_db)
    _make_gnomad_db(gnomad_db)

    cfg = {
        "paths": {
            "clinvar_db": clinvar_db,
            "gnomad_db": gnomad_db,
            "gnomad_vcf_dir": str(tmp_path / "empty_gnomad_vcf"),
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
            "genes": ["CYP2D6", "CYP2C19"],
            "risk_factor_genes": ["APOE"],
        },
    }
    cfg_path = str(tmp_path / "test_config.yaml")
    with open(cfg_path, "w") as f:
        yaml.dump(cfg, f)

    monkeypatch.setenv("GB_CONFIG_PATH", cfg_path)
    reset()

    import scripts.db.query_local_clinvar as cv_mod
    import scripts.db.query_local_gnomad as gn_mod
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None

    from scripts.db.version_manager import get_all_db_versions
    from scripts.counselor.generate_pdf import generate_report_html

    db_versions = get_all_db_versions(skip_api=True)

    # The test ClinVar DB has build_date = '2026-03-24T00:00:00'
    assert db_versions["ClinVar"]["build_date"] == "2026-03-24T00:00:00"

    report_data = {
        "sample_id": "TEST002",
        "date": "2026-03-24",
        "variants": [],
        "pgx_results": [],
        "summary": {
            "total": 0, "pathogenic": 0, "likely_pathogenic": 0,
            "drug_response": 0, "risk_factor": 0, "vus": 0,
            "likely_benign": 0, "benign": 0,
        },
        "db_versions": db_versions,
        "pipeline": {"skip_api": True, "krgdb_path": ""},
        "mode": "cancer",
        "hpo_results": [],
    }

    html = generate_report_html(report_data, mode="cancer")
    assert "2026-03-24" in html

    reset()
    cv_mod.close()
    cv_mod._conn = None
    gn_mod.close()
    gn_mod._conn = None
