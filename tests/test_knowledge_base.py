"""Tests for Knowledge Base SQLite and Wiki management."""
import sqlite3
from pathlib import Path

from scripts.db.build_kb_db import build_kb_db


def test_build_kb_db_creates_tables(tmp_path):
    """KB database has board_decisions table and variant_stats view."""
    db_path = tmp_path / "kb.sqlite3"
    build_kb_db(str(db_path))
    conn = sqlite3.connect(str(db_path))
    tables = [r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'"
    ).fetchall()]
    assert "board_decisions" in tables
    views = [r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='view'"
    ).fetchall()]
    assert "variant_stats" in views
    conn.close()


def test_build_kb_db_idempotent(tmp_path):
    """Building twice doesn't error."""
    db_path = tmp_path / "kb.sqlite3"
    build_kb_db(str(db_path))
    build_kb_db(str(db_path))


def test_kb_indexes_exist(tmp_path):
    """Required indexes are created."""
    db_path = tmp_path / "kb.sqlite3"
    build_kb_db(str(db_path))
    conn = sqlite3.connect(str(db_path))
    indexes = [r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='index'"
    ).fetchall()]
    assert "idx_bd_variant" in indexes
    assert "idx_bd_hgvsp" in indexes
    conn.close()
