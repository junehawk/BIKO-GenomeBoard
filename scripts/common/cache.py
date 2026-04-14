"""Variant-level response cache (SQLite) to avoid redundant API calls.

Backward-compatible namespace support — the legacy `get_cached`/`set_cached`
API continues to serve variant-coordinate callers (clinvar/gnomad), while new
`get_cached_ns`/`set_cached_ns` functions share the same underlying database
via a `namespace` column. OncoKB/CIViC and other per-source caches live here
rather than spinning up parallel cache trees.
"""

import json
import sqlite3
import threading
import time
import logging
from pathlib import Path
from typing import Any, Optional, Dict
from scripts.common.config import get

logger = logging.getLogger(__name__)

_conn = None
_lock = threading.Lock()

DEFAULT_TTL_SECONDS = 7 * 24 * 3600  # 7 days
DEFAULT_NAMESPACE = "variant"


def _column_names(conn: sqlite3.Connection, table: str) -> set[str]:
    return {row[1] for row in conn.execute(f"PRAGMA table_info({table})").fetchall()}


def _migrate_namespace_column(conn: sqlite3.Connection) -> None:
    """Add the `namespace` column to pre-namespace cache databases.

    Idempotent: safe to call on fresh and already-migrated DBs.
    Existing rows inherit `namespace = source` so clinvar/gnomad rows
    remain addressable under their source name.
    """
    cols = _column_names(conn, "cache")
    if "namespace" in cols:
        return
    conn.execute(f"ALTER TABLE cache ADD COLUMN namespace TEXT NOT NULL DEFAULT '{DEFAULT_NAMESPACE}'")
    # Carry existing source labels forward into namespace so they stay queryable.
    conn.execute("UPDATE cache SET namespace = source WHERE namespace = ?", (DEFAULT_NAMESPACE,))
    conn.execute("CREATE INDEX IF NOT EXISTS idx_cache_namespace ON cache(namespace)")
    conn.commit()


def _get_connection() -> sqlite3.Connection:
    global _conn
    with _lock:
        if _conn is not None:
            return _conn

        cache_path = get("cache.path", "data/cache/variant_cache.sqlite3")
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)

        _conn = sqlite3.connect(cache_path, check_same_thread=False)
        _conn.execute("PRAGMA journal_mode=WAL")
        _conn.execute(
            f"""
            CREATE TABLE IF NOT EXISTS cache (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL,
                source TEXT NOT NULL,
                created_at REAL NOT NULL,
                namespace TEXT NOT NULL DEFAULT '{DEFAULT_NAMESPACE}'
            )
            """
        )
        _migrate_namespace_column(_conn)
        _conn.commit()
        return _conn


def _make_key(chrom: str, pos: int, ref: str, alt: str, source: str) -> str:
    """Create cache key from variant coordinates + source (clinvar/gnomad)."""
    return f"{source}:{chrom}:{pos}:{ref}:{alt}"


def _make_ns_key(namespace: str, key: str) -> str:
    """Build a namespaced cache key. Uses `::` separator to avoid collisions
    with legacy single-colon variant keys (clinvar:chr17:7577120:G:A)."""
    return f"{namespace}::{key}"


def get_cached(chrom: str, pos: int, ref: str, alt: str, source: str) -> Optional[Dict]:
    """Get cached response. Returns None if not cached or expired."""
    conn = _get_connection()
    ttl = get("cache.ttl_seconds", DEFAULT_TTL_SECONDS)
    key = _make_key(chrom, pos, ref, alt, source)

    cursor = conn.execute("SELECT value, created_at FROM cache WHERE key = ?", (key,))
    row = cursor.fetchone()
    if row is None:
        return None

    value, created_at = row
    if time.time() - created_at > ttl:
        # Expired
        conn.execute("DELETE FROM cache WHERE key = ?", (key,))
        conn.commit()
        return None

    return json.loads(value)


def set_cached(chrom: str, pos: int, ref: str, alt: str, source: str, data: Dict) -> None:
    """Store response in cache."""
    conn = _get_connection()
    key = _make_key(chrom, pos, ref, alt, source)

    conn.execute(
        "INSERT OR REPLACE INTO cache (key, value, source, created_at, namespace) VALUES (?, ?, ?, ?, ?)",
        (key, json.dumps(data, default=str), source, time.time(), source),
    )
    conn.commit()


def get_cached_ns(namespace: str, key: str) -> Optional[Any]:
    """Namespaced lookup for non-variant-coordinate caches (e.g. OncoKB).

    Returns the deserialised value, or None when missing/expired.
    """
    if not namespace or not key:
        return None
    conn = _get_connection()
    ttl = get("cache.ttl_seconds", DEFAULT_TTL_SECONDS)
    composite = _make_ns_key(namespace, key)

    cursor = conn.execute(
        "SELECT value, created_at FROM cache WHERE key = ? AND namespace = ?",
        (composite, namespace),
    )
    row = cursor.fetchone()
    if row is None:
        return None

    value, created_at = row
    if time.time() - created_at > ttl:
        conn.execute(
            "DELETE FROM cache WHERE key = ? AND namespace = ?",
            (composite, namespace),
        )
        conn.commit()
        return None

    return json.loads(value)


def set_cached_ns(namespace: str, key: str, value: Any) -> None:
    """Store arbitrary JSON-serialisable payload under (namespace, key)."""
    if not namespace or not key:
        raise ValueError("namespace and key must be non-empty")
    conn = _get_connection()
    composite = _make_ns_key(namespace, key)
    conn.execute(
        "INSERT OR REPLACE INTO cache (key, value, source, created_at, namespace) VALUES (?, ?, ?, ?, ?)",
        (composite, json.dumps(value, default=str), namespace, time.time(), namespace),
    )
    conn.commit()


def clear_cache() -> int:
    """Clear all cached entries. Returns count of deleted entries."""
    conn = _get_connection()
    cursor = conn.execute("SELECT COUNT(*) FROM cache")
    count = cursor.fetchone()[0]
    conn.execute("DELETE FROM cache")
    conn.commit()
    return count


def purge_expired() -> int:
    """Remove expired entries. Returns count of purged entries."""
    conn = _get_connection()
    ttl = get("cache.ttl_seconds", DEFAULT_TTL_SECONDS)
    cutoff = time.time() - ttl
    cursor = conn.execute("SELECT COUNT(*) FROM cache WHERE created_at < ?", (cutoff,))
    count = cursor.fetchone()[0]
    conn.execute("DELETE FROM cache WHERE created_at < ?", (cutoff,))
    conn.commit()
    return count


def cache_stats() -> Dict:
    """Get cache statistics."""
    conn = _get_connection()
    total = conn.execute("SELECT COUNT(*) FROM cache").fetchone()[0]
    ttl = get("cache.ttl_seconds", DEFAULT_TTL_SECONDS)
    cutoff = time.time() - ttl
    valid = conn.execute("SELECT COUNT(*) FROM cache WHERE created_at >= ?", (cutoff,)).fetchone()[0]
    expired = total - valid

    clinvar = conn.execute("SELECT COUNT(*) FROM cache WHERE source = 'clinvar'").fetchone()[0]
    gnomad = conn.execute("SELECT COUNT(*) FROM cache WHERE source = 'gnomad'").fetchone()[0]

    by_namespace: Dict[str, int] = {}
    for ns, n in conn.execute("SELECT namespace, COUNT(*) FROM cache GROUP BY namespace"):
        by_namespace[ns] = n

    return {
        "total": total,
        "valid": valid,
        "expired": expired,
        "clinvar": clinvar,
        "gnomad": gnomad,
        "by_namespace": by_namespace,
    }


def close():
    global _conn
    with _lock:
        if _conn:
            _conn.close()
            _conn = None
