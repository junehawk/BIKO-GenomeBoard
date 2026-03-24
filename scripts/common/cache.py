"""Variant-level response cache (SQLite) to avoid redundant API calls."""
import json
import sqlite3
import threading
import time
import logging
from pathlib import Path
from typing import Optional, Dict
from scripts.common.config import get

logger = logging.getLogger(__name__)

_conn = None
_lock = threading.Lock()

DEFAULT_TTL_SECONDS = 7 * 24 * 3600  # 7 days


def _get_connection() -> sqlite3.Connection:
    global _conn
    with _lock:
        if _conn is not None:
            return _conn

        cache_path = get("cache.path", "data/cache/variant_cache.sqlite3")
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)

        _conn = sqlite3.connect(cache_path, check_same_thread=False)
        _conn.execute("PRAGMA journal_mode=WAL")
        _conn.execute("""
            CREATE TABLE IF NOT EXISTS cache (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL,
                source TEXT NOT NULL,
                created_at REAL NOT NULL
            )
        """)
        _conn.commit()
        return _conn


def _make_key(chrom: str, pos: int, ref: str, alt: str, source: str) -> str:
    """Create cache key from variant coordinates + source (clinvar/gnomad)."""
    return f"{source}:{chrom}:{pos}:{ref}:{alt}"


def get_cached(chrom: str, pos: int, ref: str, alt: str, source: str) -> Optional[Dict]:
    """Get cached response. Returns None if not cached or expired."""
    conn = _get_connection()
    ttl = get("cache.ttl_seconds", DEFAULT_TTL_SECONDS)
    key = _make_key(chrom, pos, ref, alt, source)

    cursor = conn.execute(
        "SELECT value, created_at FROM cache WHERE key = ?", (key,)
    )
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
        "INSERT OR REPLACE INTO cache (key, value, source, created_at) VALUES (?, ?, ?, ?)",
        (key, json.dumps(data, default=str), source, time.time())
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

    return {"total": total, "valid": valid, "expired": expired, "clinvar": clinvar, "gnomad": gnomad}


def close():
    global _conn
    with _lock:
        if _conn:
            _conn.close()
            _conn = None
