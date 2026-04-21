"""Unit + smoke coverage for the KOVA v7 build script and query module.

The real KOVA v7 TSV is 3.3 GB, so these tests operate against a tiny
synthetic TSV generated in ``tmp_path``. They exercise:

* Build-script round trip (schema creation, row filtering on missing
  ``kova_AF``, metadata rows, PRIMARY KEY uniqueness).
* Query-module return shape against a populated DB.
* Missing-DB log-once hardening (mirrors the regression suite already
  in place for :mod:`scripts.storage.query_local_clingen` and
  :mod:`scripts.storage.query_gnomad_constraint`).
"""

from __future__ import annotations

import gzip
import sqlite3
from pathlib import Path

import pytest

from scripts.common.models import Variant

_KOVA_HEADER = "chrom\tpos\tref_allele\talt_allele\tkova_AC\tkova_AF\tkova_AN\tKOVA_homozygote_count\n"

_SAMPLE_ROWS = [
    # (chrom, pos, ref, alt, kova_AC, kova_AF, kova_AN, homozygote)
    ("chr1", "100", "A", "G", "5", "6.1350e-03", "326", "0"),
    # kova_AF=NA → must be skipped by the builder.
    ("chr1", "200", "C", "T", "", "NA", "100", "0"),
    # kova_AF=. → also skipped.
    ("chr2", "55", "G", "A", "1", ".", "200", "0"),
    ("chr17", "7675088", "C", "A", "12", "3.7000e-03", "3240", "1"),
    # Malformed pos — must be skipped but not crash.
    ("chr3", "notapos", "T", "G", "3", "0.001", "500", "0"),
]


def _write_synthetic_tsv(tsv_path: Path, *, gz: bool = True) -> Path:
    lines = [_KOVA_HEADER]
    for row in _SAMPLE_ROWS:
        lines.append("\t".join(row) + "\n")
    payload = "".join(lines)
    if gz:
        with gzip.open(tsv_path, "wt", encoding="utf-8") as f:
            f.write(payload)
    else:
        tsv_path.write_text(payload, encoding="utf-8")
    return tsv_path


@pytest.fixture
def built_kova_db(tmp_path: Path) -> Path:
    """Build a small KOVA SQLite DB in tmp_path and return its filesystem path."""
    from scripts.storage.build_kova_db import build_db

    tsv = _write_synthetic_tsv(tmp_path / "kova.tsv.gz")
    db_path = tmp_path / "kova.sqlite3"
    build_db(tsv, db_path, source_url="test://kova", kova_version="v7-test")
    return db_path


# ── build_kova_db ───────────────────────────────────────────────────────────


def test_build_kova_db_roundtrip(built_kova_db: Path) -> None:
    conn = sqlite3.connect(built_kova_db)
    try:
        rows = conn.execute(
            "SELECT chrom, pos, ref, alt, kova_af, kova_ac, kova_an, kova_homozygote FROM kova_af ORDER BY chrom, pos"
        ).fetchall()
        meta = dict(conn.execute("SELECT key, value FROM metadata").fetchall())
    finally:
        conn.close()

    # Expect exactly two rows — the NA, ".", and bad-pos rows were skipped.
    assert rows == [
        ("chr1", 100, "A", "G", 0.006135, 5, 326, 0),
        ("chr17", 7675088, "C", "A", 0.0037, 12, 3240, 1),
    ]
    assert meta["record_count"] == "2"
    assert meta["kova_version"] == "v7-test"
    assert meta["source"] == "test://kova"
    assert meta["assembly"] == "GRCh38"
    assert "build_date" in meta


def test_build_kova_db_idempotent(tmp_path: Path) -> None:
    """Re-running the build on the same input yields the same row count —
    DROP TABLE IF EXISTS + INSERT OR REPLACE keep it safe."""
    from scripts.storage.build_kova_db import build_db

    tsv = _write_synthetic_tsv(tmp_path / "k.tsv.gz")
    db_path = tmp_path / "k.sqlite3"

    build_db(tsv, db_path)
    build_db(tsv, db_path)  # second run must not raise

    conn = sqlite3.connect(db_path)
    try:
        (count,) = conn.execute("SELECT COUNT(*) FROM kova_af").fetchone()
    finally:
        conn.close()
    assert count == 2


def test_build_kova_db_missing_column_raises(tmp_path: Path) -> None:
    """If the KOVA TSV is missing a required column, the builder surfaces a
    clear ValueError instead of silently producing an empty DB."""
    from scripts.storage.build_kova_db import build_db

    bad = tmp_path / "bad.tsv.gz"
    with gzip.open(bad, "wt", encoding="utf-8") as f:
        f.write("chrom\tpos\tref_allele\talt_allele\n")  # missing kova_AF etc.
        f.write("chr1\t1\tA\tG\n")
    db_path = tmp_path / "bad.sqlite3"

    with pytest.raises(ValueError, match="missing required columns"):
        build_db(bad, db_path)


# ── query_kova ──────────────────────────────────────────────────────────────


def test_query_kova_hit(built_kova_db: Path) -> None:
    from scripts.population.query_kova import query_kova, reset_availability_cache

    reset_availability_cache()
    v = Variant(chrom="chr17", pos=7675088, ref="C", alt="A")
    result = query_kova(v, str(built_kova_db))
    assert result is not None
    assert result["kova_af"] == pytest.approx(0.0037)
    assert result["kova_ac"] == 12
    assert result["kova_an"] == 3240
    assert result["kova_homozygote"] == 1
    assert result["api_available"] is True


def test_query_kova_miss_returns_none(built_kova_db: Path) -> None:
    from scripts.population.query_kova import query_kova, reset_availability_cache

    reset_availability_cache()
    v = Variant(chrom="chrX", pos=99999999, ref="A", alt="T")
    assert query_kova(v, str(built_kova_db)) is None


def test_query_kova_accepts_bare_chrom(built_kova_db: Path) -> None:
    """The pipeline stores chroms with the chr prefix but some callers pass
    the bare form. The query module must normalise either shape."""
    from scripts.population.query_kova import query_kova, reset_availability_cache

    reset_availability_cache()
    v = Variant(chrom="17", pos=7675088, ref="C", alt="A")
    result = query_kova(v, str(built_kova_db))
    assert result is not None
    assert result["kova_ac"] == 12


def test_missing_kova_db_returns_none_and_logs_once(tmp_path, caplog) -> None:
    """Missing paths.kova_db: every query returns None without raising,
    and the user sees exactly ONE warning for the run (not one per variant)."""
    from scripts.population.query_kova import query_kova, reset_availability_cache

    reset_availability_cache()
    missing_path = str(tmp_path / "definitely_not_built.sqlite3")

    variants = [
        Variant(chrom="chr1", pos=100, ref="A", alt="G"),
        Variant(chrom="chr2", pos=200, ref="C", alt="T"),
        Variant(chrom="chr17", pos=7675088, ref="C", alt="A"),
        Variant(chrom="chrX", pos=1, ref="A", alt="T"),
    ]
    with caplog.at_level("WARNING"):
        for v in variants:
            assert query_kova(v, missing_path) is None

    not_found = [r for r in caplog.records if "KOVA local DB not found" in r.getMessage()]
    assert len(not_found) == 1


def test_empty_shell_kova_db_logs_once(tmp_path, caplog) -> None:
    """DB file exists but has no kova_af table — log once, return None."""
    from scripts.population.query_kova import query_kova, reset_availability_cache

    reset_availability_cache()
    shell = tmp_path / "empty.sqlite3"
    sqlite3.connect(shell).close()  # zero-table shell

    v = Variant(chrom="chr1", pos=100, ref="A", alt="G")
    with caplog.at_level("WARNING"):
        assert query_kova(v, str(shell)) is None
        assert query_kova(v, str(shell)) is None
        assert query_kova(v, str(shell)) is None

    shell_warnings = [r for r in caplog.records if "has no `kova_af` table" in r.getMessage()]
    assert len(shell_warnings) == 1


def test_get_db_version_returns_metadata(built_kova_db: Path) -> None:
    from scripts.population.query_kova import get_db_version, reset_availability_cache

    reset_availability_cache()
    meta = get_db_version(str(built_kova_db))
    assert meta["source"] == "local_db"
    assert meta["version"] == "v7-test"
    assert meta["record_count"] == "2"
    assert meta["assembly"] == "GRCh38"


def test_get_db_version_missing_returns_not_available(tmp_path) -> None:
    from scripts.population.query_kova import get_db_version, reset_availability_cache

    reset_availability_cache()
    meta = get_db_version(str(tmp_path / "never_built.sqlite3"))
    assert meta == {"source": "not_available"}


def test_availability_cache_reset_refires(tmp_path, caplog) -> None:
    from scripts.population.query_kova import query_kova, reset_availability_cache

    reset_availability_cache()
    missing = str(tmp_path / "never_built.sqlite3")
    v = Variant(chrom="chr1", pos=1, ref="A", alt="G")

    with caplog.at_level("WARNING"):
        query_kova(v, missing)
    first = sum(1 for r in caplog.records if "KOVA local DB not found" in r.getMessage())
    assert first == 1

    reset_availability_cache()
    with caplog.at_level("WARNING"):
        query_kova(v, missing)
    second = sum(1 for r in caplog.records if "KOVA local DB not found" in r.getMessage())
    assert second == 2
