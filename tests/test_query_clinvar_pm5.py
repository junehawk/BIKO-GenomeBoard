"""Regression tests for v2.3-T6 self-computed PM5 path.

Covers ``scripts/storage/query_local_clinvar.get_clinvar_pathogenic_positions``
and the wiring inside ``scripts/pipeline/classify.classify_variants`` that
feeds it into ``evidence_collector.collect_additional_evidence``.

These tests build a tiny throw-away ClinVar SQLite DB with the v2.3-T6
schema (including the ``hgvsp`` column) and patch ``paths.clinvar_db`` to
point at it, so the test suite never depends on the real 4 GB build.
"""

from __future__ import annotations

import sqlite3

import pytest

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _write_test_db(db_path: str) -> None:
    """Build a minimal ClinVar SQLite that mirrors the v2.3-T6 schema."""
    conn = sqlite3.connect(db_path)
    conn.execute(
        """
        CREATE TABLE variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            rsid TEXT,
            gene TEXT,
            clinical_significance TEXT,
            review_status TEXT,
            phenotype_list TEXT,
            variation_id TEXT,
            allele_id TEXT,
            origin TEXT,
            assembly TEXT DEFAULT 'GRCh38',
            last_evaluated TEXT,
            number_submitters INTEGER,
            hgvsp TEXT
        )
        """
    )
    conn.execute("CREATE INDEX idx_clinvar_gene_hgvsp ON variants(gene, hgvsp)")
    conn.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT)")

    rows = [
        # TP53 known hotspots — Pathogenic, with HGVSp.
        ("chr17", 7676154, "C", "T", "TP53", "Pathogenic", "criteria provided", "p.Arg175His"),
        ("chr17", 7674220, "G", "A", "TP53", "Likely pathogenic", "criteria provided", "p.Arg248Gln"),
        ("chr17", 7674894, "G", "A", "TP53", "Pathogenic", "criteria provided", "p.Arg273His"),
        # TP53 conflicting — must NOT contribute to the pathogenic position set.
        (
            "chr17",
            7674230,
            "G",
            "A",
            "TP53",
            "Conflicting interpretations of pathogenicity",
            "criteria provided",
            "p.Arg999Lys",
        ),
        # TP53 entry without HGVSp (e.g. CNV / intronic) — must be skipped.
        ("chr17", 7676000, "A", "G", "TP53", "Pathogenic", "criteria provided", None),
        # BRCA1 unrelated — distinct gene, must stay isolated.
        ("chr17", 43093000, "C", "T", "BRCA1", "Pathogenic", "criteria provided", "p.Cys61Gly"),
    ]
    conn.executemany(
        """
        INSERT INTO variants
            (chrom, pos, ref, alt, gene, clinical_significance, review_status, hgvsp)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )
    conn.execute("INSERT INTO metadata VALUES ('schema_version', '2')")
    conn.execute("INSERT INTO metadata VALUES ('hgvsp_count', '5')")
    conn.commit()
    conn.close()


@pytest.fixture
def patched_clinvar_db(tmp_path, monkeypatch):
    """Write a minimal v2.3-T6 ClinVar DB and point config at it."""
    db_path = tmp_path / "clinvar_t6.sqlite3"
    _write_test_db(str(db_path))

    import scripts.storage.query_local_clinvar as qmod
    from scripts.common import config as cfgmod

    # Force the lazy config getter to return our path. We patch via monkeypatch
    # on the module-level ``get`` so we do not have to mutate the global yaml.
    real_get = cfgmod.get

    def _patched_get(key, default=None):
        if key == "paths.clinvar_db":
            return str(db_path)
        return real_get(key, default)

    monkeypatch.setattr("scripts.storage.query_local_clinvar.get", _patched_get)

    qmod.close()
    qmod.reset_cache_for_tests()

    yield db_path

    qmod.close()
    qmod.reset_cache_for_tests()


# ---------------------------------------------------------------------------
# get_clinvar_pathogenic_positions
# ---------------------------------------------------------------------------


def test_get_pathogenic_positions_tp53(patched_clinvar_db):
    from scripts.storage.query_local_clinvar import get_clinvar_pathogenic_positions

    positions = get_clinvar_pathogenic_positions("TP53")
    # All three P/LP entries with parseable HGVSp must appear.
    assert 175 in positions
    assert 248 in positions
    assert 273 in positions


def test_get_pathogenic_positions_ignores_conflict(patched_clinvar_db):
    from scripts.storage.query_local_clinvar import get_clinvar_pathogenic_positions

    positions = get_clinvar_pathogenic_positions("TP53")
    # 999 came from a "Conflicting interpretations of pathogenicity" row and
    # must NOT contribute to the set, even though the HGVSp parses cleanly.
    assert 999 not in positions


def test_get_pathogenic_positions_unknown_gene(patched_clinvar_db):
    from scripts.storage.query_local_clinvar import get_clinvar_pathogenic_positions

    assert get_clinvar_pathogenic_positions("BOGUS_GENE_XYZ") == set()


def test_get_pathogenic_positions_no_hgvsp_rows_skipped(patched_clinvar_db):
    from scripts.storage.query_local_clinvar import get_clinvar_pathogenic_positions

    # The TP53 row at 7676000 has hgvsp=NULL — it must be silently skipped
    # rather than crashing extract_protein_position with a None argument.
    positions = get_clinvar_pathogenic_positions("TP53")
    # Exactly the three valid rows survive.
    assert positions == {175, 248, 273}


def test_get_pathogenic_positions_isolates_genes(patched_clinvar_db):
    from scripts.storage.query_local_clinvar import get_clinvar_pathogenic_positions

    brca = get_clinvar_pathogenic_positions("BRCA1")
    assert brca == {61}
    # A second-gene query must not pollute the first-gene cache.
    tp53 = get_clinvar_pathogenic_positions("TP53")
    assert 61 not in tp53


def test_cache_reused(patched_clinvar_db, monkeypatch):
    """A second call for the same gene must NOT re-issue the SELECT."""
    import scripts.storage.query_local_clinvar as qmod
    from scripts.storage.query_local_clinvar import get_clinvar_pathogenic_positions

    # Prime the cache.
    first = get_clinvar_pathogenic_positions("TP53")
    assert first == {175, 248, 273}

    # Replace the active connection with one that raises on any execute().
    class _BoomConn:
        def execute(self, *a, **kw):
            raise AssertionError("cache miss — get_clinvar_pathogenic_positions should not re-query")

        def close(self):
            pass

    qmod._conn = _BoomConn()

    second = get_clinvar_pathogenic_positions("TP53")
    assert second == first


def test_availability_cache_reset(patched_clinvar_db):
    """reset_cache_for_tests() must let the next call re-query the DB."""
    import scripts.storage.query_local_clinvar as qmod
    from scripts.storage.query_local_clinvar import get_clinvar_pathogenic_positions

    first = get_clinvar_pathogenic_positions("TP53")
    qmod.reset_cache_for_tests()
    # After reset both the per-gene memo and the hgvsp-availability probe
    # must be cleared, so the next call must hit SQLite again and produce
    # the same result.
    assert qmod._HGVSP_AVAILABLE is None
    assert qmod._PATHOGENIC_POS_CACHE == {}
    second = get_clinvar_pathogenic_positions("TP53")
    assert second == first


def test_db_missing_returns_empty(monkeypatch, tmp_path):
    """When the DB file does not exist the function must return an empty set."""
    import scripts.storage.query_local_clinvar as qmod

    qmod.close()
    qmod.reset_cache_for_tests()

    missing = tmp_path / "no_such_db.sqlite3"

    def _patched_get(key, default=None):
        if key == "paths.clinvar_db":
            return str(missing)
        return default

    monkeypatch.setattr("scripts.storage.query_local_clinvar.get", _patched_get)

    from scripts.storage.query_local_clinvar import get_clinvar_pathogenic_positions

    assert get_clinvar_pathogenic_positions("TP53") == set()


def test_legacy_db_without_hgvsp_column(tmp_path, monkeypatch):
    """A legacy DB (pre v2.3-T6) without the hgvsp column must degrade gracefully."""
    import scripts.storage.query_local_clinvar as qmod

    db_path = tmp_path / "legacy_clinvar.sqlite3"
    conn = sqlite3.connect(str(db_path))
    conn.execute(
        """
        CREATE TABLE variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
            gene TEXT, clinical_significance TEXT
        )
        """
    )
    conn.execute(
        "INSERT INTO variants (chrom, pos, ref, alt, gene, clinical_significance) "
        "VALUES ('chr17', 1, 'A', 'T', 'TP53', 'Pathogenic')"
    )
    conn.commit()
    conn.close()

    def _patched_get(key, default=None):
        if key == "paths.clinvar_db":
            return str(db_path)
        return default

    monkeypatch.setattr("scripts.storage.query_local_clinvar.get", _patched_get)
    qmod.close()
    qmod.reset_cache_for_tests()

    from scripts.storage.query_local_clinvar import get_clinvar_pathogenic_positions

    # No crash, no Exception — just an empty set.
    assert get_clinvar_pathogenic_positions("TP53") == set()
    qmod.close()
    qmod.reset_cache_for_tests()


# ---------------------------------------------------------------------------
# End-to-end: PM5 fires through classify_variants
# ---------------------------------------------------------------------------


def test_pm5_fires_via_self_computed_path(patched_clinvar_db, monkeypatch):
    """A novel TP53 missense at residue 175 must collect PM5 through the
    full classify_variants -> collect_additional_evidence wiring."""
    from scripts.common.models import Variant as VariantModel
    from scripts.storage.query_local_clinvar import get_clinvar_pathogenic_positions
    from scripts.pipeline.classify import classify_variants

    # Sanity — fixture DB really does have residue 175 in the set.
    assert 175 in get_clinvar_pathogenic_positions("TP53")

    # Novel missense at a hotspot residue: AA change differs from the
    # ClinVar entry (His → Cys instead of His). evidence_collector's PM5
    # check is purely position-based, so this must fire PM5.
    variant = VariantModel(
        chrom="chr17",
        pos=7676155,
        ref="C",
        alt="A",
        gene="TP53",
        hgvsp="p.Arg175Cys",
        consequence="missense_variant",
    )

    db_results = {
        variant.variant_id: {
            "clinvar": {
                "clinvar_significance": "Not Found",
                "clinvar_id": None,
                "acmg_codes": [],
                "review_status": "",
            },
            "gnomad": {"gnomad_all": None, "gnomad_eas": None, "api_available": False},
            "krgdb_freq": None,
        }
    }
    freq_results = {variant.variant_id: {"acmg_codes": [], "korean_flag": ""}}

    results = classify_variants([variant], db_results, freq_results, intervar_data=None)
    assert variant.variant_id in results
    codes = {c.upper() for c in results[variant.variant_id].evidence_codes}
    assert "PM5" in codes, f"PM5 missing from self-computed path. codes={codes}"
