# tests/test_cache.py
"""Tests for SQLite variant cache (scripts/common/cache.py)."""
import time
import pytest
from unittest.mock import patch
from scripts.common.models import Variant


@pytest.fixture(autouse=True)
def isolated_cache(tmp_path):
    """Each test gets a fresh in-memory/temp SQLite cache."""
    import scripts.common.cache as cache_mod

    # Close and reset the global connection before each test
    cache_mod.close()
    cache_mod._conn = None

    # Point cache at a temp file so tests don't pollute each other
    tmp_db = str(tmp_path / "test_cache.sqlite3")
    with patch("scripts.common.cache.get", side_effect=lambda key, default=None: {
        "cache.path": tmp_db,
        "cache.ttl_seconds": 604800,
    }.get(key, default)):
        yield cache_mod

    # Cleanup after test
    cache_mod.close()
    cache_mod._conn = None


def test_set_and_get_cached(isolated_cache):
    cache = isolated_cache
    data = {"clinvar_significance": "Pathogenic", "acmg_codes": ["PS1"]}
    cache.set_cached("chr17", 7577120, "G", "A", "clinvar", data)
    result = cache.get_cached("chr17", 7577120, "G", "A", "clinvar")
    assert result is not None
    assert result["clinvar_significance"] == "Pathogenic"
    assert result["acmg_codes"] == ["PS1"]


def test_cache_miss_returns_none(isolated_cache):
    cache = isolated_cache
    result = cache.get_cached("chr1", 99999, "A", "T", "clinvar")
    assert result is None


def test_cache_expired_returns_none(tmp_path):
    """Entry stored with TTL=0 should be treated as immediately expired."""
    import scripts.common.cache as cache_mod
    cache_mod.close()
    cache_mod._conn = None

    tmp_db = str(tmp_path / "expired_cache.sqlite3")

    # Store with normal TTL
    with patch("scripts.common.cache.get", side_effect=lambda key, default=None: {
        "cache.path": tmp_db,
        "cache.ttl_seconds": 604800,
    }.get(key, default)):
        cache_mod._conn = None
        cache_mod.set_cached("chr1", 12345, "A", "T", "clinvar", {"clinvar_significance": "Benign"})

    cache_mod.close()
    cache_mod._conn = None

    # Retrieve with TTL=0 (everything is expired)
    with patch("scripts.common.cache.get", side_effect=lambda key, default=None: {
        "cache.path": tmp_db,
        "cache.ttl_seconds": 0,
    }.get(key, default)):
        cache_mod._conn = None
        result = cache_mod.get_cached("chr1", 12345, "A", "T", "clinvar")

    assert result is None
    cache_mod.close()
    cache_mod._conn = None


def test_clear_cache(isolated_cache):
    cache = isolated_cache
    cache.set_cached("chr1", 111, "A", "T", "clinvar", {"x": 1})
    cache.set_cached("chr2", 222, "G", "C", "gnomad", {"y": 2})
    count = cache.clear_cache()
    assert count == 2
    assert cache.get_cached("chr1", 111, "A", "T", "clinvar") is None
    assert cache.get_cached("chr2", 222, "G", "C", "gnomad") is None


def test_purge_expired(tmp_path):
    """Entries older than TTL should be purged."""
    import scripts.common.cache as cache_mod
    cache_mod.close()
    cache_mod._conn = None

    tmp_db = str(tmp_path / "purge_cache.sqlite3")

    with patch("scripts.common.cache.get", side_effect=lambda key, default=None: {
        "cache.path": tmp_db,
        "cache.ttl_seconds": 604800,
    }.get(key, default)):
        cache_mod._conn = None
        # Insert one fresh entry
        cache_mod.set_cached("chr1", 100, "A", "T", "clinvar", {"val": "fresh"})
        # Manually insert an old entry directly into the DB
        conn = cache_mod._get_connection()
        old_ts = time.time() - 604800 - 10  # 7 days + 10 seconds ago
        conn.execute(
            "INSERT OR REPLACE INTO cache (key, value, source, created_at) VALUES (?, ?, ?, ?)",
            ("clinvar:chr2:200:G:C", '{"val":"old"}', "clinvar", old_ts),
        )
        conn.commit()

        purged = cache_mod.purge_expired()
        assert purged == 1
        # Fresh entry still exists
        assert cache_mod.get_cached("chr1", 100, "A", "T", "clinvar") is not None

    cache_mod.close()
    cache_mod._conn = None


def test_cache_stats(isolated_cache):
    cache = isolated_cache
    cache.set_cached("chr1", 1, "A", "T", "clinvar", {"v": 1})
    cache.set_cached("chr2", 2, "G", "C", "gnomad", {"v": 2})
    stats = cache.cache_stats()
    assert stats["total"] == 2
    assert stats["valid"] == 2
    assert stats["expired"] == 0
    assert stats["clinvar"] == 1
    assert stats["gnomad"] == 1


def test_cache_key_format(isolated_cache):
    cache = isolated_cache
    key = cache._make_key("chr17", 7577120, "G", "A", "clinvar")
    assert key == "clinvar:chr17:7577120:G:A"

    key2 = cache._make_key("chrX", 1000, "C", "T", "gnomad")
    assert key2 == "gnomad:chrX:1000:C:T"


def test_clinvar_uses_cache(tmp_path, mocker):
    """Second call to query_clinvar should hit cache and not call the API."""
    import scripts.common.cache as cache_mod
    cache_mod.close()
    cache_mod._conn = None

    tmp_db = str(tmp_path / "clinvar_cache.sqlite3")

    with patch("scripts.common.cache.get", side_effect=lambda key, default=None: {
        "cache.path": tmp_db,
        "cache.ttl_seconds": 604800,
    }.get(key, default)):
        cache_mod._conn = None

        mock_search = mocker.patch(
            "scripts.clinical.query_clinvar._search_clinvar_variant",
            return_value={
                "clinical_significance": {"description": "Pathogenic"},
                "gene": {"symbol": "TP53"},
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "variation_id": "12375",
            },
        )
        # Patch cache functions used by query_clinvar to use our isolated cache
        mocker.patch("scripts.clinical.query_clinvar.get_cached", side_effect=cache_mod.get_cached)
        mocker.patch("scripts.clinical.query_clinvar.set_cached", side_effect=cache_mod.set_cached)

        from scripts.clinical.query_clinvar import query_clinvar
        variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")

        # First call — API hit
        result1 = query_clinvar(variant)
        assert result1["clinvar_significance"] == "Pathogenic"
        assert mock_search.call_count == 1

        # Second call — should come from cache
        result2 = query_clinvar(variant)
        assert result2["clinvar_significance"] == "Pathogenic"
        assert mock_search.call_count == 1  # still 1, not 2

    cache_mod.close()
    cache_mod._conn = None


def test_gnomad_uses_cache(tmp_path, mocker):
    """Second call to query_gnomad should hit cache and not call the API."""
    import scripts.common.cache as cache_mod
    cache_mod.close()
    cache_mod._conn = None

    tmp_db = str(tmp_path / "gnomad_cache.sqlite3")

    with patch("scripts.common.cache.get", side_effect=lambda key, default=None: {
        "cache.path": tmp_db,
        "cache.ttl_seconds": 604800,
    }.get(key, default)):
        cache_mod._conn = None

        mock_graphql = mocker.patch(
            "scripts.korean_pop.query_gnomad._graphql_query",
            return_value={
                "data": {
                    "variant": {
                        "variant_id": "17-7577120-G-A",
                        "genome": {
                            "ac": 200,
                            "an": 1000000,
                            "populations": [{"id": "eas", "ac": 3, "an": 10000}],
                        },
                        "exome": None,
                    }
                }
            },
        )
        mocker.patch("scripts.korean_pop.query_gnomad.get_cached", side_effect=cache_mod.get_cached)
        mocker.patch("scripts.korean_pop.query_gnomad.set_cached", side_effect=cache_mod.set_cached)

        from scripts.korean_pop.query_gnomad import query_gnomad
        variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")

        # First call — API hit
        result1 = query_gnomad(variant)
        assert result1["gnomad_all"] is not None
        assert mock_graphql.call_count >= 1
        first_count = mock_graphql.call_count

        # Second call — should come from cache
        result2 = query_gnomad(variant)
        assert result2["gnomad_all"] == result1["gnomad_all"]
        assert mock_graphql.call_count == first_count  # no new API calls

    cache_mod.close()
    cache_mod._conn = None
