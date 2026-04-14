"""OncoKB free-tier client — timeout / 429 / connection-refused / offline mode.

The client is a thin wrapper around ``/annotate/mutations/byProteinChange``.
v2.2 A1 requires graceful degradation: any transient network failure must
raise ``OncoKBUnavailable`` so the curator can fall back to CIViC-only.
"""

from __future__ import annotations

import time
from unittest.mock import MagicMock

import pytest


@pytest.fixture(autouse=True)
def clear_ns_cache():
    """Wipe the oncokb cache namespace between tests so stubs don't leak."""
    from scripts.common import cache

    conn = cache._get_connection()
    conn.execute("DELETE FROM cache WHERE namespace = ?", ("oncokb",))
    conn.commit()
    yield
    conn.execute("DELETE FROM cache WHERE namespace = ?", ("oncokb",))
    conn.commit()


def test_offline_mode_shortcircuits_without_http(monkeypatch):
    """When offline_mode=True, the client returns [] without touching requests."""
    from scripts.clinical import oncokb_client

    called = []

    def fake_get(*a, **kw):
        called.append(a)
        raise RuntimeError("network should not be touched")

    monkeypatch.setattr(oncokb_client.requests, "get", fake_get)
    result = oncokb_client.annotate_protein_change("KRAS", "G12D", offline_mode=True)
    assert result == []
    assert called == []


def test_http_429_raises_oncokb_unavailable(monkeypatch):
    from scripts.clinical import oncokb_client

    response = MagicMock()
    response.status_code = 429
    response.json.return_value = {"error": "rate limit"}
    monkeypatch.setattr(oncokb_client.requests, "get", lambda *a, **kw: response)

    with pytest.raises(oncokb_client.OncoKBUnavailable):
        oncokb_client.annotate_protein_change("KRAS", "G12D", offline_mode=False)


def test_connection_refused_raises_oncokb_unavailable(monkeypatch):
    from scripts.clinical import oncokb_client

    def raise_conn_err(*a, **kw):
        raise oncokb_client.requests.exceptions.ConnectionError("refused")

    monkeypatch.setattr(oncokb_client.requests, "get", raise_conn_err)
    with pytest.raises(oncokb_client.OncoKBUnavailable):
        oncokb_client.annotate_protein_change("KRAS", "G12D", offline_mode=False)


def test_timeout_raises_oncokb_unavailable(monkeypatch):
    from scripts.clinical import oncokb_client

    def raise_timeout(*a, **kw):
        raise oncokb_client.requests.exceptions.Timeout("slow")

    monkeypatch.setattr(oncokb_client.requests, "get", raise_timeout)
    with pytest.raises(oncokb_client.OncoKBUnavailable):
        oncokb_client.annotate_protein_change("KRAS", "G12D", offline_mode=False)


def test_happy_path_returns_hits_and_caches(monkeypatch):
    from scripts.clinical import oncokb_client

    payload = {
        "query": {"hugoSymbol": "KRAS", "alteration": "G12D"},
        "treatments": [
            {
                "drugs": [{"drugName": "Sotorasib"}],
                "level": "LEVEL_1",
                "pmids": ["32955176"],
                "indication": "Colorectal adenocarcinoma",
            },
        ],
        "mutationEffect": {"knownEffect": "Gain-of-function"},
    }
    response = MagicMock()
    response.status_code = 200
    response.json.return_value = payload
    calls = {"count": 0}

    def fake_get(*a, **kw):
        calls["count"] += 1
        return response

    monkeypatch.setattr(oncokb_client.requests, "get", fake_get)
    result = oncokb_client.annotate_protein_change("KRAS", "G12D", offline_mode=False)
    assert any("sotorasib" in (d.get("drug") or "").lower() for d in result)
    assert all(d.get("pmids") for d in result)

    # Second call served from cache — requests.get not hit again
    result2 = oncokb_client.annotate_protein_change("KRAS", "G12D", offline_mode=False)
    assert result2 == result
    assert calls["count"] == 1


def test_5xx_raises_oncokb_unavailable(monkeypatch):
    from scripts.clinical import oncokb_client

    response = MagicMock()
    response.status_code = 503
    response.json.return_value = {"error": "unavailable"}
    monkeypatch.setattr(oncokb_client.requests, "get", lambda *a, **kw: response)
    with pytest.raises(oncokb_client.OncoKBUnavailable):
        oncokb_client.annotate_protein_change("KRAS", "G12D", offline_mode=False)


def test_walltime_happy_path_cached_is_under_1s(monkeypatch):
    """Cached lookup must be near-instant (well under 1 s)."""
    from scripts.clinical import oncokb_client
    from scripts.common import cache

    cache.set_cached_ns("oncokb", "KRAS::G12D", [{"drug": "Sotorasib", "pmids": ["32955176"], "level": "LEVEL_1"}])

    def fail(*a, **kw):
        raise AssertionError("should hit cache")

    monkeypatch.setattr(oncokb_client.requests, "get", fail)
    t0 = time.monotonic()
    for _ in range(30):
        oncokb_client.annotate_protein_change("KRAS", "G12D", offline_mode=False)
    elapsed = time.monotonic() - t0
    assert elapsed < 1.0, f"cached hot-path too slow: {elapsed:.2f}s"
