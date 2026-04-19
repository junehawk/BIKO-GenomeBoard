"""Tests for centralized configuration management."""

import os

import pytest


@pytest.fixture(autouse=True)
def reset_config():
    """Reset config cache before and after each test."""
    from scripts.common import config as cfg_module

    cfg_module.reset()
    # Also clear any GB_CONFIG_PATH env var set by tests
    old_env = os.environ.pop("GB_CONFIG_PATH", None)
    old_ncbi = os.environ.pop("NCBI_API_KEY", None)
    yield
    cfg_module.reset()
    if old_env is not None:
        os.environ["GB_CONFIG_PATH"] = old_env
    if old_ncbi is not None:
        os.environ["NCBI_API_KEY"] = old_ncbi


def test_load_config_defaults_without_file(tmp_path, monkeypatch):
    """When no config file exists, all defaults must be returned."""
    from scripts.common.config import load_config

    # Point to a nonexistent config path
    nonexistent = str(tmp_path / "no_config.yaml")
    config = load_config(config_path=nonexistent)
    assert config["thresholds"]["ba1"] == 0.05
    assert config["thresholds"]["bs1"] == 0.01
    assert config["thresholds"]["pm2"] == 0.001
    assert config["api"]["timeout"] == 30
    assert "CYP2D6" in config["pgx"]["genes"]
    assert "APOE" in config["pgx"]["risk_factor_genes"]


def test_load_config_from_yaml(tmp_path):
    """Values from a YAML file should override defaults."""
    from scripts.common.config import load_config

    config_file = tmp_path / "config.yaml"
    config_file.write_text("thresholds:\n  ba1: 0.10\n  bs1: 0.02\n  pm2: 0.0005\napi:\n  timeout: 60\n")
    config = load_config(config_path=str(config_file))
    assert config["thresholds"]["ba1"] == 0.10
    assert config["thresholds"]["bs1"] == 0.02
    assert config["thresholds"]["pm2"] == 0.0005
    assert config["api"]["timeout"] == 60
    # Defaults for unspecified keys must still be present
    assert config["api"]["max_retries"] == 3
    assert "CYP2D6" in config["pgx"]["genes"]


def test_load_config_gb_config_path_env(tmp_path, monkeypatch):
    """GB_CONFIG_PATH env var should point to custom config."""
    from scripts.common import config as cfg_module

    config_file = tmp_path / "custom.yaml"
    config_file.write_text("thresholds:\n  ba1: 0.07\n")
    monkeypatch.setenv("GB_CONFIG_PATH", str(config_file))
    cfg_module.reset()
    config = cfg_module.load_config()
    assert config["thresholds"]["ba1"] == 0.07
    # Other defaults still present
    assert config["thresholds"]["bs1"] == 0.01


def test_ncbi_api_key_env_override(tmp_path, monkeypatch):
    """NCBI_API_KEY env var should override ncbi_api_key in config."""
    from scripts.common.config import load_config

    nonexistent = str(tmp_path / "no_config.yaml")
    monkeypatch.setenv("NCBI_API_KEY", "test_key_123")
    config = load_config(config_path=nonexistent)
    assert config["api"]["ncbi_api_key"] == "test_key_123"


def test_get_dot_notation_simple():
    """get() with dot notation should return scalar config values."""
    from scripts.common.config import get

    assert get("thresholds.ba1") == 0.05
    assert get("thresholds.bs1") == 0.01
    assert get("api.timeout") == 30


def test_get_dot_notation_missing_returns_default():
    """get() with unknown key should return the default value."""
    from scripts.common.config import get

    assert get("nonexistent.key", 42) == 42
    assert get("thresholds.nonexistent", 99) == 99


def test_get_dot_notation_missing_returns_none():
    """get() with unknown key and no default should return None."""
    from scripts.common.config import get

    assert get("totally.made.up") is None


def test_get_nested_list():
    """get() should return list values correctly."""
    from scripts.common.config import get

    datasets = get("api.gnomad_datasets")
    assert isinstance(datasets, list)
    assert "gnomad_r4" in datasets
    assert "gnomad_r2_1" in datasets


def test_get_pgx_genes():
    """get() should return the PGx gene list."""
    from scripts.common.config import get

    genes = get("pgx.genes")
    assert isinstance(genes, list)
    for gene in ["CYP2D6", "CYP2C19", "CYP2C9", "HLA-B", "HLA-A", "NUDT15", "TPMT", "DPYD"]:
        assert gene in genes


def test_config_is_cached():
    """Second call to load_config() should return the same object."""
    from scripts.common.config import load_config

    cfg1 = load_config()
    cfg2 = load_config()
    assert cfg1 is cfg2


def test_reset_clears_cache():
    """reset() should force a fresh load on next call."""
    from scripts.common import config as cfg_module

    cfg1 = cfg_module.load_config()
    cfg_module.reset()
    cfg2 = cfg_module.load_config()
    # After reset, a new dict is created (not the same object)
    assert cfg1 is not cfg2


def test_paths_resolved_to_absolute(tmp_path):
    """Relative paths in config should be resolved to absolute paths."""
    from scripts.common.config import load_config

    config_file = tmp_path / "config.yaml"
    config_file.write_text("paths:\n  krgdb: 'data/krgdb_freq.tsv'\n")
    config = load_config(config_path=str(config_file))
    assert os.path.isabs(config["paths"]["krgdb"])


def test_load_config_with_custom_gb_config_path(tmp_path, monkeypatch):
    """Test full round-trip: custom YAML via GB_CONFIG_PATH env var."""
    from scripts.common import config as cfg_module

    config_file = tmp_path / "gb_custom.yaml"
    config_file.write_text(
        "api:\n"
        "  timeout: 45\n"
        "  pharmgkb_rate_limit: 1.0\n"
        "pgx:\n"
        "  genes: ['CYP2D6', 'CYP2C19']\n"
        "  risk_factor_genes: ['APOE', 'BRCA1']\n"
    )
    monkeypatch.setenv("GB_CONFIG_PATH", str(config_file))
    cfg_module.reset()
    config = cfg_module.load_config()
    assert config["api"]["timeout"] == 45
    assert config["api"]["pharmgkb_rate_limit"] == 1.0
    assert config["pgx"]["genes"] == ["CYP2D6", "CYP2C19"]
    assert "BRCA1" in config["pgx"]["risk_factor_genes"]
    # Defaults still filled in
    assert config["thresholds"]["ba1"] == 0.05
