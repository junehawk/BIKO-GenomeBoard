"""Tests for the gnomAD v4.1 gene constraint DB integration (v2.3-T8).

Covers:

1. The ``build_gnomad_constraint_db`` builder against a 10-gene fixture TSV.
2. ``query_gnomad_constraint.is_constrained`` thresholds — TP53 (low
   missense_z, NOT constrained), SCN1A (canonical NDD, constrained),
   borderline genes.
3. The "log once then short-circuit" availability cache contract
   (matches the v2.2 T3 ClinGen pattern).
4. End-to-end: variant_selector's de novo carve-out admits a
   protein-impacting de novo VUS in a constrained-but-non-DDG2P gene
   purely via the gnomAD constraint OR-branch.

Fixture TSV uses the exact column shape from
``gnomad.v4.1.constraint_metrics.tsv``: tab-separated, dotted column
names (``lof.pLI``, ``mis.z_score``, ...). DictReader handles them
transparently.
"""

from __future__ import annotations

import sqlite3

import pytest

# ---------------------------------------------------------------------------
# Fixture TSV — minimal columns, real gnomAD field names
# ---------------------------------------------------------------------------

_FIXTURE_HEADERS = [
    "gene",
    "transcript",
    "mane_select",
    "lof.pLI",
    "lof.oe",
    "lof.oe_ci.upper",
    "mis.z_score",
    "syn.z_score",
    "mis.oe",
]

# (gene, transcript, mane_select, pli, oe_lof, loeuf, mis_z, syn_z, oe_mis)
_FIXTURE_ROWS = [
    # Constrained NDD gene — passes both thresholds
    ("SCN1A", "ENST00000674923.1", "true", "1.0", "0.05", "0.10", "7.62", "1.77", "0.56"),
    # Constrained NDD gene — passes both thresholds
    ("STXBP1", "ENST00000373302.7", "true", "1.0", "0.06", "0.12", "5.50", "0.20", "0.45"),
    # TP53 — high pLI but very low missense Z (tumour suppressor; tolerates
    # missense). MUST return is_constrained == False.
    ("TP53", "ENST00000269305.9", "true", "0.998", "0.28", "0.45", "1.15", "0.96", "0.86"),
    # Borderline pLI gene — pli below 0.9, missense_z high → False
    ("BORDER1", "ENST00000000001.1", "true", "0.50", "0.40", "0.55", "5.00", "0.10", "0.40"),
    # Borderline mis_z gene — pli high, missense_z below 3.09 → False
    ("BORDER2", "ENST00000000002.1", "true", "0.95", "0.20", "0.35", "2.50", "0.10", "0.50"),
    # Non-MANE row — must be skipped entirely
    ("SKIPME", "ENST00000000003.1", "false", "1.0", "0.01", "0.05", "9.99", "0.10", "0.10"),
    # Generic constrained gene with NaN syn_z (still admit pLI/mis_z)
    ("GOODGENE", "ENST00000000004.1", "true", "0.99", "0.10", "0.20", "4.50", "NA", "0.40"),
    # Empty pLI value — should round-trip to None and is_constrained False
    ("EMPTYPLI", "ENST00000000005.1", "true", "", "0.30", "0.45", "5.00", "0.20", "0.50"),
]


def _write_fixture_tsv(path) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(_FIXTURE_HEADERS) + "\n")
        for row in _FIXTURE_ROWS:
            f.write("\t".join(row) + "\n")


@pytest.fixture
def tmp_constraint_db(tmp_path):
    from scripts.db.build_gnomad_constraint_db import build_db
    from scripts.db.query_gnomad_constraint import reset_availability_cache

    tsv_path = tmp_path / "gnomad_constraint.tsv"
    db_path = tmp_path / "gnomad_constraint.sqlite3"
    _write_fixture_tsv(tsv_path)
    build_db(str(tsv_path), str(db_path))
    # Probes are memoised per-path, but tests rotate through fresh tmp_paths;
    # clear the cache so each test starts from a clean availability state.
    reset_availability_cache()
    return str(db_path)


# ---------------------------------------------------------------------------
# 1. Builder
# ---------------------------------------------------------------------------


def test_build_gnomad_constraint_db(tmp_path):
    from scripts.db.build_gnomad_constraint_db import build_db

    tsv_path = tmp_path / "gnomad_constraint.tsv"
    db_path = tmp_path / "out.sqlite3"
    _write_fixture_tsv(tsv_path)

    result = build_db(str(tsv_path), str(db_path))
    assert db_path.exists()

    # SKIPME (mane_select=false) is dropped → 7 rows from 8 fixture rows.
    assert result["records"] == 7

    conn = sqlite3.connect(str(db_path))
    try:
        # Schema sanity — columns we promised the variant_selector.
        cols = [r[1] for r in conn.execute("PRAGMA table_info(constraint_metrics)").fetchall()]
        for needed in ("gene", "pli", "loeuf", "mis_z", "syn_z", "oe_lof", "oe_mis"):
            assert needed in cols

        # Index is in place
        idx = conn.execute("SELECT name FROM sqlite_master WHERE type='index' AND name='idx_gnomad_gene'").fetchone()
        assert idx is not None

        # Gene-level lookups
        scn1a = conn.execute("SELECT pli, mis_z FROM constraint_metrics WHERE gene='SCN1A'").fetchone()
        assert scn1a == (1.0, 7.62)

        tp53 = conn.execute("SELECT pli, mis_z FROM constraint_metrics WHERE gene='TP53'").fetchone()
        assert tp53 == (0.998, 1.15)

        # Non-MANE row absent
        skipme = conn.execute("SELECT 1 FROM constraint_metrics WHERE gene='SKIPME'").fetchone()
        assert skipme is None

        # NaN/empty handling
        emptypli = conn.execute("SELECT pli, mis_z FROM constraint_metrics WHERE gene='EMPTYPLI'").fetchone()
        assert emptypli == (None, 5.0)

        goodgene = conn.execute("SELECT syn_z FROM constraint_metrics WHERE gene='GOODGENE'").fetchone()
        assert goodgene == (None,)  # 'NA' → None

        # Metadata table populated
        meta = dict(conn.execute("SELECT key, value FROM metadata").fetchall())
        assert meta["gnomad_version"] == "v4.1"
        assert meta["record_count"] == "7"
        assert "constraint_metrics.tsv" in meta["source"]
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# 2. is_constrained / get_constraint thresholds
# ---------------------------------------------------------------------------


def test_get_constraint_returns_metrics_dict(tmp_constraint_db):
    from scripts.db.query_gnomad_constraint import get_constraint

    metrics = get_constraint("SCN1A", tmp_constraint_db)
    assert metrics is not None
    assert metrics["pli"] == 1.0
    assert metrics["missense_z"] == 7.62
    assert metrics["loeuf"] == 0.10
    assert metrics["oe_lof"] == 0.05


def test_is_constrained_scn1a_true(tmp_constraint_db):
    from scripts.db.query_gnomad_constraint import is_constrained

    # SCN1A: pLI=1.0 (>=0.9), mis_z=7.62 (>=3.09) → True
    assert is_constrained("SCN1A", tmp_constraint_db) is True


def test_is_constrained_tp53_false_due_to_low_missense_z(tmp_constraint_db):
    """TP53 has high pLI but tolerates missense (low mis_z); admission
    requires BOTH thresholds, so the de novo carve-out must reject it."""
    from scripts.db.query_gnomad_constraint import is_constrained

    assert is_constrained("TP53", tmp_constraint_db) is False


def test_is_constrained_borderline_pli_false(tmp_constraint_db):
    """BORDER1 has mis_z above 3.09 but pLI below 0.9 — must be False."""
    from scripts.db.query_gnomad_constraint import is_constrained

    assert is_constrained("BORDER1", tmp_constraint_db) is False


def test_is_constrained_borderline_missense_z_false(tmp_constraint_db):
    """BORDER2 has pLI above 0.9 but mis_z below 3.09 — must be False."""
    from scripts.db.query_gnomad_constraint import is_constrained

    assert is_constrained("BORDER2", tmp_constraint_db) is False


def test_is_constrained_unknown_gene_false(tmp_constraint_db):
    from scripts.db.query_gnomad_constraint import get_constraint, is_constrained

    assert get_constraint("DOES_NOT_EXIST_XYZ", tmp_constraint_db) is None
    assert is_constrained("DOES_NOT_EXIST_XYZ", tmp_constraint_db) is False


def test_is_constrained_empty_pli_false(tmp_constraint_db):
    """A row with NULL pLI must not crash and must return False."""
    from scripts.db.query_gnomad_constraint import is_constrained

    assert is_constrained("EMPTYPLI", tmp_constraint_db) is False


def test_is_constrained_none_or_empty_gene(tmp_constraint_db):
    from scripts.db.query_gnomad_constraint import get_constraint, is_constrained

    assert is_constrained("", tmp_constraint_db) is False
    assert is_constrained(None, tmp_constraint_db) is False  # type: ignore[arg-type]
    assert get_constraint("", tmp_constraint_db) is None


# ---------------------------------------------------------------------------
# 3. Log-once / availability-cache contract (v2.2 T3 pattern)
# ---------------------------------------------------------------------------


def test_missing_db_returns_none_and_logs_once(tmp_path, caplog):
    from scripts.db.query_gnomad_constraint import (
        get_constraint,
        is_constrained,
        reset_availability_cache,
    )

    reset_availability_cache()
    missing = str(tmp_path / "never_built.sqlite3")

    with caplog.at_level("WARNING"):
        assert get_constraint("SCN1A", missing) is None
        assert is_constrained("STXBP1", missing) is False
        assert get_constraint("TP53", missing) is None
        assert is_constrained("FAKE", missing) is False

    not_found = [r for r in caplog.records if "gnomAD constraint DB not found" in r.getMessage()]
    assert len(not_found) == 1


def test_empty_shell_db_returns_none_and_logs_once(tmp_path, caplog):
    """A SQLite file with no constraint_metrics table is the
    setup_databases.sh half-failure mode. Must log once and short-circuit."""
    from scripts.db.query_gnomad_constraint import (
        get_constraint,
        is_constrained,
        reset_availability_cache,
    )

    reset_availability_cache()
    empty = str(tmp_path / "empty_shell.sqlite3")
    sqlite3.connect(empty).close()  # zero-table shell

    with caplog.at_level("WARNING"):
        assert get_constraint("SCN1A", empty) is None
        assert is_constrained("STXBP1", empty) is False
        assert get_constraint("TP53", empty) is None

    shell_warnings = [r for r in caplog.records if "no `constraint_metrics` table" in r.getMessage()]
    assert len(shell_warnings) == 1


def test_reset_cache_refires_warning(tmp_path, caplog):
    from scripts.db.query_gnomad_constraint import (
        get_constraint,
        reset_availability_cache,
    )

    reset_availability_cache()
    missing = str(tmp_path / "still_missing.sqlite3")

    with caplog.at_level("WARNING"):
        get_constraint("SCN1A", missing)
    first = sum(1 for r in caplog.records if "gnomAD constraint DB not found" in r.getMessage())
    assert first == 1

    reset_availability_cache()
    with caplog.at_level("WARNING"):
        get_constraint("SCN1A", missing)
    second = sum(1 for r in caplog.records if "gnomAD constraint DB not found" in r.getMessage())
    assert second == 2


# ---------------------------------------------------------------------------
# 4. End-to-end: variant_selector's OR-branch is now LIVE
# ---------------------------------------------------------------------------


def test_variant_selector_carve_out_uses_constraint(monkeypatch, tmp_constraint_db):
    """De novo missense VUS in a constrained-but-non-DDG2P gene must reach
    the rare-disease MAY arm via the constraint branch alone.

    Picks ``GOODGENE`` (a fictional symbol guaranteed to be absent from the
    DDG2P panel). Pre-T8 the OR-branch was dead, so this variant would have
    been rejected; post-T8 it admits with reason ``VUS_denovo_neurodev``.
    """
    from scripts.clinical_board import variant_selector
    from scripts.db import query_gnomad_constraint
    from scripts.db.query_gnomad_constraint import reset_availability_cache

    # Wire the variant_selector's lazy import to the fixture DB.
    monkeypatch.setattr(
        query_gnomad_constraint,
        "_get_db_path",
        lambda db_path=None: tmp_constraint_db,
    )
    reset_availability_cache()

    # Sanity: GOODGENE is not in DDG2P → admission must come from constraint.
    from scripts.common.ddg2p_panel import is_admitted_neurodev_gene

    assert is_admitted_neurodev_gene("GOODGENE") is False
    assert query_gnomad_constraint.is_constrained("GOODGENE") is True

    variant = {
        "gene": "GOODGENE",
        "classification": "VUS",
        "consequence": "Missense",
        "hgvsp": "p.Arg101Leu",
        "variant_inheritance": "de_novo",
        "confirmed_denovo": False,
        "hpo_score": 0,
        "gnomad_af": 0.0,
        "in_silico": {},
        "tier": "",
        "cancer_gene_type": "",
        "oncokb_level": "",
    }
    selected, meta = variant_selector.select_board_variants([variant], mode="rare-disease")
    assert len(selected) == 1
    assert selected[0]["selection_reason"] == "VUS_denovo_neurodev"
    assert meta["selected"] == 1


def test_variant_selector_rejects_unconstrained_non_ddg2p(monkeypatch, tmp_constraint_db):
    """Companion sanity test: a de novo missense in an unlisted, low-pLI gene
    must still be rejected — the OR-branch should not become a free pass."""
    from scripts.clinical_board import variant_selector
    from scripts.db import query_gnomad_constraint
    from scripts.db.query_gnomad_constraint import reset_availability_cache

    monkeypatch.setattr(
        query_gnomad_constraint,
        "_get_db_path",
        lambda db_path=None: tmp_constraint_db,
    )
    reset_availability_cache()

    # BORDER1 has low pLI → constraint admission denied; not in DDG2P either.
    from scripts.common.ddg2p_panel import is_admitted_neurodev_gene

    assert is_admitted_neurodev_gene("BORDER1") is False
    assert query_gnomad_constraint.is_constrained("BORDER1") is False

    variant = {
        "gene": "BORDER1",
        "classification": "VUS",
        "consequence": "Missense",
        "hgvsp": "p.Arg101Leu",
        "variant_inheritance": "de_novo",
        "confirmed_denovo": False,
        "hpo_score": 0,
        "gnomad_af": 0.0,
        "in_silico": {},
        "tier": "",
        "cancer_gene_type": "",
        "oncokb_level": "",
    }
    selected, _ = variant_selector.select_board_variants([variant], mode="rare-disease")
    assert len(selected) == 0
