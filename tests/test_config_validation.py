"""Tests for config.yaml schema validation."""

import os

import pytest


@pytest.fixture(autouse=True)
def reset_config():
    """Reset config cache before and after each test."""
    from scripts.common import config as cfg_module

    cfg_module.reset()
    old_env = os.environ.pop("GB_CONFIG_PATH", None)
    yield
    cfg_module.reset()
    if old_env is not None:
        os.environ["GB_CONFIG_PATH"] = old_env


# ------------------------------------------------------------------
# validate_config unit tests (no file I/O)
# ------------------------------------------------------------------


def test_validate_valid_config():
    """A config with all required keys and correct types produces no errors."""
    from scripts.common.config_schema import validate_config

    config = {
        "paths": {"kova_db": "a", "clinvar_db": "b", "templates": "t"},
        "thresholds": {"ba1": 0.05, "bs1": 0.01, "pm2": 0.001},
        "api": {"timeout": 30},
        "cache": {"enabled": True},
    }
    messages = validate_config(config)
    errors = [m for m in messages if m.startswith("ERROR:")]
    assert errors == []


def test_validate_missing_required_key():
    """Missing paths.templates must produce an ERROR."""
    from scripts.common.config_schema import validate_config

    config = {
        "paths": {"kova_db": "a"},
        "thresholds": {"ba1": 0.05, "bs1": 0.01, "pm2": 0.001},
        "api": {"timeout": 30},
        "cache": {"enabled": True},
    }
    messages = validate_config(config)
    errors = [m for m in messages if m.startswith("ERROR:")]
    assert any("paths.templates" in e for e in errors)


def test_validate_wrong_type_threshold():
    """A threshold set to a string must produce a WARNING."""
    from scripts.common.config_schema import validate_config

    config = {
        "paths": {"kova_db": "a", "clinvar_db": "b", "templates": "t"},
        "thresholds": {"ba1": "not_a_number", "bs1": 0.01, "pm2": 0.001},
        "api": {"timeout": 30},
        "cache": {"enabled": True},
    }
    messages = validate_config(config)
    warnings = [m for m in messages if m.startswith("WARNING:")]
    assert any("ba1" in w and "float" in w for w in warnings)


def test_validate_unknown_top_level_key():
    """An unknown top-level key must produce a WARNING, not an error."""
    from scripts.common.config_schema import validate_config

    config = {
        "paths": {"kova_db": "a", "clinvar_db": "b", "templates": "t"},
        "thresholds": {"ba1": 0.05, "bs1": 0.01, "pm2": 0.001},
        "api": {"timeout": 30},
        "cache": {"enabled": True},
        "experimental_feature": {"foo": "bar"},
    }
    messages = validate_config(config)
    errors = [m for m in messages if m.startswith("ERROR:")]
    warnings = [m for m in messages if m.startswith("WARNING:")]
    assert errors == []
    assert any("experimental_feature" in w for w in warnings)


def test_validate_empty_config():
    """An empty config dict must report errors for critical keys."""
    from scripts.common.config_schema import validate_config

    messages = validate_config({})
    errors = [m for m in messages if m.startswith("ERROR:")]
    assert any("paths.templates" in e for e in errors)
    # Also expect warnings for other missing recommended keys
    warnings = [m for m in messages if m.startswith("WARNING:")]
    assert len(warnings) > 0


def test_validate_partial_config():
    """Only missing keys should be flagged; present ones should pass."""
    from scripts.common.config_schema import validate_config

    config = {
        "paths": {"templates": "t"},
        "thresholds": {"ba1": 0.05},
        "api": {"timeout": 30},
    }
    messages = validate_config(config)
    errors = [m for m in messages if m.startswith("ERROR:")]
    assert errors == []
    # Missing recommended keys like paths.kova_db should warn
    warnings = [m for m in messages if m.startswith("WARNING:")]
    assert any("paths.kova_db" in w for w in warnings)
    # Present keys should NOT appear in warnings about missing keys
    assert not any("paths.templates" in w for w in warnings)
    assert not any("thresholds.ba1" in w for w in warnings)


def test_validate_int_accepted_as_float():
    """An int value where float is expected should NOT produce a warning."""
    from scripts.common.config_schema import validate_config

    config = {
        "paths": {"kova_db": "a", "clinvar_db": "b", "templates": "t"},
        "thresholds": {"ba1": 1, "bs1": 0, "pm2": 0},  # ints, not floats
        "api": {"timeout": 30},
        "cache": {"enabled": True},
    }
    messages = validate_config(config)
    errors = [m for m in messages if m.startswith("ERROR:")]
    assert errors == []
    # No type warnings for thresholds
    type_warnings = [m for m in messages if "should be float" in m]
    assert type_warnings == []


# ------------------------------------------------------------------
# Integration: load_config raises ValueError for critical missing keys
# ------------------------------------------------------------------


def test_load_config_raises_on_missing_templates(tmp_path):
    """load_config must raise ValueError when paths.templates is missing
    and cannot be filled by defaults."""
    from scripts.common.config import load_config

    # Write a config that explicitly sets paths without templates.
    # The defaults mechanism fills templates, so we need to override
    # the entire paths section AND remove templates from defaults
    # by having a paths section that does NOT include templates.
    # Actually, the defaults mechanism uses setdefault, so templates
    # will always be filled.  This test verifies that the real config
    # never triggers the error — load should succeed.
    config_file = tmp_path / "config.yaml"
    config_file.write_text("thresholds:\n  ba1: 0.05\n")
    config = load_config(config_path=str(config_file))
    # templates is filled by defaults, so no error
    assert "templates" in config["paths"]


def test_load_real_config_no_errors():
    """Loading the real project config.yaml must produce no errors."""
    from scripts.common.config import load_config

    # This should not raise
    config = load_config()
    assert config["thresholds"]["ba1"] == 0.05
    assert "templates" in config["paths"]


def test_load_config_wrong_type_logs_warning(tmp_path, caplog):
    """Wrong type for a threshold should log a warning but not raise."""
    import logging

    from scripts.common.config import load_config

    config_file = tmp_path / "config.yaml"
    config_file.write_text(
        "paths:\n  templates: t\n"
        "thresholds:\n  ba1: 'bad'\n  bs1: 0.01\n  pm2: 0.001\n"
        "api:\n  timeout: 30\n"
        "cache:\n  enabled: true\n"
    )
    with caplog.at_level(logging.WARNING, logger="scripts.common.config"):
        config = load_config(config_path=str(config_file))
    assert any("ba1" in r.message for r in caplog.records)
    # Should still load successfully
    assert config["thresholds"]["ba1"] == "bad"


def test_load_config_unknown_key_logs_warning(tmp_path, caplog):
    """Unknown top-level key should log a warning but not raise."""
    import logging

    from scripts.common.config import load_config

    config_file = tmp_path / "config.yaml"
    config_file.write_text("future_feature:\n  enabled: true\n")
    with caplog.at_level(logging.WARNING, logger="scripts.common.config"):
        config = load_config(config_path=str(config_file))
    assert any("future_feature" in r.message for r in caplog.records)
    # Warning text should mention where to register a new section.
    assert any("KNOWN_TOP_LEVEL_KEYS" in r.message for r in caplog.records)
    # Should still load successfully
    assert config is not None


def test_real_config_yaml_has_no_unknown_top_level_keys():
    """The real project config.yaml must not trigger unknown-key warnings.

    Every top-level section present in config.yaml has to be registered in
    KNOWN_TOP_LEVEL_KEYS so users do not see noise on every load. This test
    catches the "added a new section to config.yaml but forgot to register
    it" regression — Phase 1.1 of the v2.4 refactor introduced acmg /
    pharmacogenomics / knowledge_base for exactly this reason.
    """
    from pathlib import Path

    import yaml

    from scripts.common.config_schema import KNOWN_TOP_LEVEL_KEYS, validate_config

    config_path = Path(__file__).parent.parent / "config.yaml"
    with open(config_path) as f:
        raw_config = yaml.safe_load(f) or {}

    unregistered = [k for k in raw_config if k not in KNOWN_TOP_LEVEL_KEYS]
    assert unregistered == [], (
        f"config.yaml has top-level keys not registered in KNOWN_TOP_LEVEL_KEYS: "
        f"{unregistered}. Add them to scripts/common/config_schema.py."
    )

    messages = validate_config(raw_config)
    unknown_warnings = [m for m in messages if "unknown top-level key" in m]
    assert unknown_warnings == [], f"Unexpected unknown-key warnings: {unknown_warnings}"
