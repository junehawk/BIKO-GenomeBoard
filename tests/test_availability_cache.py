"""Unit tests for :mod:`scripts.common.availability_cache`.

These cover the log-once + short-circuit contract consolidated from the
v2.2 T3 ClinGen hardening pattern. They do not exercise any specific
DB — the behaviour is verified with ``tmp_path`` and fake probes so the
tests stay fast and don't depend on optional local DBs (CI guard).
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest


# ── AvailabilityCache class ──────────────────────────────────────────


def test_probe_runs_once_per_path():
    from scripts.common.availability_cache import AvailabilityCache

    calls: list[str] = []

    def probe(p: str) -> bool:
        calls.append(p)
        return True

    cache = AvailabilityCache()
    assert cache.check("/tmp/a", probe) is True
    assert cache.check("/tmp/a", probe) is True
    assert cache.check("/tmp/a", probe) is True
    assert calls == ["/tmp/a"], "probe must run at most once per path"


def test_probe_false_is_cached_like_true():
    from scripts.common.availability_cache import AvailabilityCache

    calls = {"n": 0}

    def probe(_p: str) -> bool:
        calls["n"] += 1
        return False

    cache = AvailabilityCache()
    for _ in range(5):
        assert cache.check("/tmp/missing", probe) is False
    assert calls["n"] == 1


def test_reset_refires_probe():
    from scripts.common.availability_cache import AvailabilityCache

    calls = {"n": 0}

    def probe(_p: str) -> bool:
        calls["n"] += 1
        return False

    cache = AvailabilityCache()
    cache.check("/tmp/x", probe)
    cache.check("/tmp/x", probe)
    assert calls["n"] == 1

    cache.reset()
    cache.check("/tmp/x", probe)
    assert calls["n"] == 2


def test_distinct_paths_probed_independently():
    from scripts.common.availability_cache import AvailabilityCache

    seen: list[str] = []

    def probe(p: str) -> bool:
        seen.append(p)
        return True

    cache = AvailabilityCache()
    cache.check("/tmp/a", probe)
    cache.check("/tmp/b", probe)
    cache.check("/tmp/a", probe)
    assert sorted(seen) == ["/tmp/a", "/tmp/b"]


# ── check_availability() module-level helper ─────────────────────────


def test_check_availability_missing_file_logs_once(tmp_path, caplog):
    from scripts.common.availability_cache import check_availability, reset_shared_cache

    reset_shared_cache()
    missing = str(tmp_path / "never_built.sqlite3")

    with caplog.at_level("WARNING"):
        assert check_availability("Fakesource", missing) is False
        assert check_availability("Fakesource", missing) is False
        assert check_availability("Fakesource", missing) is False

    warnings = [r for r in caplog.records if "Fakesource local DB not found" in r.getMessage()]
    assert len(warnings) == 1


def test_check_availability_existing_file_returns_true(tmp_path, caplog):
    from scripts.common.availability_cache import check_availability, reset_shared_cache

    reset_shared_cache()
    p = tmp_path / "present.txt"
    p.write_text("data")

    with caplog.at_level("WARNING"):
        assert check_availability("Fakelabel", str(p)) is True

    assert not any("not found" in r.getMessage() for r in caplog.records)


def test_check_availability_empty_sqlite_fails_default_probe(tmp_path, caplog):
    from scripts.common.availability_cache import check_availability, reset_shared_cache

    reset_shared_cache()
    empty = tmp_path / "empty.sqlite3"
    sqlite3.connect(str(empty)).close()  # no tables

    with caplog.at_level("WARNING"):
        assert check_availability("SQLprobe", str(empty)) is False

    assert any("SQLprobe local DB not found" in r.getMessage() for r in caplog.records)


def test_check_availability_custom_probe_and_message(tmp_path, caplog):
    from scripts.common.availability_cache import check_availability, reset_shared_cache

    reset_shared_cache()
    target = str(tmp_path / "custom.blob")

    calls = {"n": 0}

    def probe(p: str) -> bool:
        calls["n"] += 1
        return False

    custom_msg = "CUSTOM_MESSAGE_SENTINEL for unit test"
    with caplog.at_level("WARNING"):
        assert check_availability("X", target, probe=probe, message=custom_msg) is False
        assert check_availability("X", target, probe=probe, message=custom_msg) is False

    assert calls["n"] == 1
    custom = [r for r in caplog.records if "CUSTOM_MESSAGE_SENTINEL" in r.getMessage()]
    assert len(custom) == 1


def test_reset_shared_cache_refires_warning(tmp_path, caplog):
    from scripts.common.availability_cache import check_availability, reset_shared_cache

    reset_shared_cache()
    missing = str(tmp_path / "still_missing.sqlite3")

    with caplog.at_level("WARNING"):
        check_availability("Resetprobe", missing)
    first = sum(1 for r in caplog.records if "Resetprobe local DB not found" in r.getMessage())
    assert first == 1

    reset_shared_cache()
    with caplog.at_level("WARNING"):
        check_availability("Resetprobe", missing)
    second = sum(1 for r in caplog.records if "Resetprobe local DB not found" in r.getMessage())
    assert second == 2


# ── Integration with existing call-sites ─────────────────────────────


def test_clingen_still_logs_once_after_refactor(tmp_path, caplog):
    """Regression: the ClinGen module's public warning-log contract
    must survive the internal refactor onto AvailabilityCache."""
    from scripts.storage.query_local_clingen import (
        get_gene_validity_local,
        reset_availability_cache,
    )

    reset_availability_cache()
    missing = str(tmp_path / "gone.sqlite3")

    with caplog.at_level("WARNING"):
        for _ in range(5):
            assert get_gene_validity_local("TP53", missing) is None

    not_found = [r for r in caplog.records if "ClinGen local DB not found" in r.getMessage()]
    assert len(not_found) == 1


def test_gnomad_constraint_still_logs_once_after_refactor(tmp_path, caplog):
    """Regression: the gnomAD-constraint module's public warning-log
    contract must survive the internal refactor onto AvailabilityCache."""
    from scripts.storage.query_gnomad_constraint import (
        get_constraint,
        reset_availability_cache,
    )

    reset_availability_cache()
    missing = str(tmp_path / "no_constraint.sqlite3")

    with caplog.at_level("WARNING"):
        for _ in range(5):
            assert get_constraint("SCN1A", missing) is None

    not_found = [r for r in caplog.records if "gnomAD constraint DB not found" in r.getMessage()]
    assert len(not_found) == 1


def test_caches_are_independent_between_modules(tmp_path, caplog):
    """Resetting one call-site's cache must not disturb others — each
    module owns its own AvailabilityCache instance."""
    from scripts.storage.query_gnomad_constraint import (
        reset_availability_cache as reset_constraint,
    )
    from scripts.storage.query_local_clingen import (
        get_gene_validity_local,
        reset_availability_cache as reset_clingen,
    )

    reset_clingen()
    reset_constraint()
    missing = str(tmp_path / "shared_missing.sqlite3")

    with caplog.at_level("WARNING"):
        assert get_gene_validity_local("TP53", missing) is None

    # Reset only the constraint cache; ClinGen's memoised False should
    # remain, so a second ClinGen call must NOT re-emit the warning.
    caplog.clear()
    reset_constraint()
    with caplog.at_level("WARNING"):
        assert get_gene_validity_local("TP53", missing) is None

    clingen_warnings = [r for r in caplog.records if "ClinGen local DB not found" in r.getMessage()]
    assert clingen_warnings == [], "ClinGen cache must not be cleared by reset_constraint"


def test_default_probe_rejects_missing_then_present(tmp_path):
    """``_default_probe`` is the fallback probe used by
    :func:`check_availability` when no custom probe is supplied."""
    from scripts.common.availability_cache import _default_probe

    missing = tmp_path / "nope"
    assert _default_probe(missing) is False

    existing = tmp_path / "data.txt"
    existing.write_text("hello")
    assert _default_probe(existing) is True


def test_default_probe_sqlite_needs_tables(tmp_path):
    from scripts.common.availability_cache import _default_probe

    empty = tmp_path / "empty.sqlite3"
    sqlite3.connect(str(empty)).close()
    assert _default_probe(empty) is False

    populated = tmp_path / "populated.sqlite3"
    conn = sqlite3.connect(str(populated))
    conn.execute("CREATE TABLE t (x INTEGER)")
    conn.commit()
    conn.close()
    assert _default_probe(populated) is True


def test_default_probe_accepts_path_and_str(tmp_path):
    from scripts.common.availability_cache import _default_probe

    f = tmp_path / "a.txt"
    f.write_text("x")
    assert _default_probe(str(f)) is True
    assert _default_probe(Path(f)) is True


@pytest.mark.parametrize("suffix", [".sqlite3", ".db"])
def test_default_probe_sqlite_suffix_triggers_table_check(tmp_path, suffix):
    from scripts.common.availability_cache import _default_probe

    empty = tmp_path / f"shell{suffix}"
    sqlite3.connect(str(empty)).close()
    assert _default_probe(empty) is False
