"""Tests for Knowledge Base SQLite and Wiki management."""
import sqlite3
from pathlib import Path

import pytest

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


# ---------------------------------------------------------------------------
# Task 11 — KnowledgeBase storage
# ---------------------------------------------------------------------------


def test_knowledge_base_auto_initialises_empty_db(tmp_path):
    """Regression: a 0-byte or missing kb.sqlite3 must auto-initialise on first use.

    Previously, pointing KnowledgeBase at an empty file caused save_decision to
    fail with `no such table: board_decisions`. KnowledgeBase now ensures the
    schema exists during __init__.
    """
    from scripts.clinical_board.knowledge_base import KnowledgeBase

    db_path = tmp_path / "kb.sqlite3"
    db_path.write_bytes(b"")  # 0-byte file, like the gitignored placeholder
    assert db_path.exists() and db_path.stat().st_size == 0

    kb = KnowledgeBase(str(db_path), str(tmp_path))
    kb.save_decision(
        sample_id="S001",
        mode="cancer",
        gene="EGFR",
        variant="chr7:55259515:T:G",
        board_diagnosis="TKI sensitive",
        board_confidence="high",
        agent_consensus="unanimous",
    )

    conn = sqlite3.connect(str(db_path))
    rows = conn.execute("SELECT sample_id FROM board_decisions").fetchall()
    conn.close()
    assert len(rows) == 1
    assert rows[0][0] == "S001"


def test_knowledge_base_auto_initialises_missing_db(tmp_path):
    """KnowledgeBase works against a db_path that doesn't exist yet."""
    from scripts.clinical_board.knowledge_base import KnowledgeBase

    db_path = tmp_path / "nested" / "kb.sqlite3"
    assert not db_path.exists()

    kb = KnowledgeBase(str(db_path), str(tmp_path))
    kb.save_decision(
        sample_id="S042",
        mode="rare-disease",
        gene="BRCA2",
        variant="chr13:32936732:A:G",
        board_diagnosis="Likely pathogenic",
        board_confidence="high",
        agent_consensus="unanimous",
    )
    assert db_path.exists()


def test_save_board_decision(tmp_path):
    """Saving a decision inserts into SQLite."""
    from scripts.clinical_board.knowledge_base import KnowledgeBase

    db_path = tmp_path / "kb.sqlite3"
    build_kb_db(str(db_path))
    kb = KnowledgeBase(str(db_path), str(tmp_path))
    kb.save_decision(
        sample_id="S001",
        mode="cancer",
        gene="EGFR",
        variant="chr7:55259515:T:G",
        hgvsp="p.L858R",
        classification="Pathogenic",
        board_diagnosis="TKI sensitive",
        board_confidence="high",
        clinical_context_summary="age:60s|sex:M|dx:nsclc",
        agent_consensus="unanimous",
        raw_opinion_json="{}",
    )
    conn = sqlite3.connect(str(db_path))
    rows = conn.execute("SELECT * FROM board_decisions").fetchall()
    assert len(rows) == 1
    conn.close()


def test_variant_stats_aggregation(tmp_path):
    """variant_stats view aggregates correctly."""
    from scripts.clinical_board.knowledge_base import KnowledgeBase

    db_path = tmp_path / "kb.sqlite3"
    build_kb_db(str(db_path))
    kb = KnowledgeBase(str(db_path), str(tmp_path))
    for i in range(3):
        kb.save_decision(
            sample_id=f"S00{i}",
            mode="cancer",
            gene="EGFR",
            variant="chr7:55259515:T:G",
            board_diagnosis="TKI sensitive",
            board_confidence="high",
            agent_consensus="unanimous",
        )
    conn = sqlite3.connect(str(db_path))
    stats = conn.execute(
        "SELECT gene, variant, mode, total_cases, high_confidence_count "
        "FROM variant_stats WHERE gene='EGFR'"
    ).fetchone()
    assert stats is not None
    assert stats[3] == 3  # total_cases
    assert stats[4] == 3  # high_confidence_count
    conn.close()


def test_cross_mode_isolation(tmp_path):
    """Cancer decisions don't appear in rare-disease queries."""
    from scripts.clinical_board.knowledge_base import KnowledgeBase

    db_path = tmp_path / "kb.sqlite3"
    build_kb_db(str(db_path))
    kb = KnowledgeBase(str(db_path), str(tmp_path))
    kb.save_decision(
        sample_id="S001",
        mode="cancer",
        gene="TP53",
        variant="chr17:7675088:C:A",
        board_diagnosis="Somatic driver",
        board_confidence="high",
        agent_consensus="unanimous",
    )
    conn = sqlite3.connect(str(db_path))
    rare = conn.execute(
        "SELECT * FROM variant_stats "
        "WHERE gene='TP53' AND mode='rare-disease'"
    ).fetchall()
    assert len(rare) == 0
    conn.close()


def test_generate_gene_wiki_writes_markdown(tmp_path):
    """generate_gene_wiki writes a per-mode gene markdown file."""
    from scripts.clinical_board.knowledge_base import KnowledgeBase

    db_path = tmp_path / "kb.sqlite3"
    build_kb_db(str(db_path))
    kb = KnowledgeBase(str(db_path), str(tmp_path))
    kb.save_decision(
        sample_id="S001",
        mode="cancer",
        gene="EGFR",
        variant="chr7:55259515:T:G",
        hgvsp="p.L858R",
        board_diagnosis="TKI sensitive",
        board_confidence="high",
        agent_consensus="unanimous",
    )
    kb.generate_gene_wiki("EGFR", "cancer")
    wiki = tmp_path / "cancer" / "genes" / "EGFR.md"
    assert wiki.exists()
    text = wiki.read_text(encoding="utf-8")
    assert "EGFR" in text
    assert "TKI sensitive" in text


def test_update_log_appends(tmp_path):
    """update_log appends a line to log.md."""
    from scripts.clinical_board.knowledge_base import KnowledgeBase

    db_path = tmp_path / "kb.sqlite3"
    build_kb_db(str(db_path))
    kb = KnowledgeBase(str(db_path), str(tmp_path))
    kb.update_log("S001", "cancer", date="2026-04-14")
    kb.update_log("S002", "rare-disease", date="2026-04-14")
    log = tmp_path / "log.md"
    assert log.exists()
    text = log.read_text(encoding="utf-8")
    assert "S001" in text
    assert "S002" in text


def test_update_index_lists_gene_wikis(tmp_path):
    """update_index regenerates index.md referencing gene wiki files."""
    from scripts.clinical_board.knowledge_base import KnowledgeBase

    db_path = tmp_path / "kb.sqlite3"
    build_kb_db(str(db_path))
    kb = KnowledgeBase(str(db_path), str(tmp_path))
    kb.save_decision(
        sample_id="S001",
        mode="cancer",
        gene="EGFR",
        variant="chr7:55259515:T:G",
        board_diagnosis="TKI sensitive",
        board_confidence="high",
        agent_consensus="unanimous",
    )
    kb.generate_gene_wiki("EGFR", "cancer")
    kb.update_index()
    index = tmp_path / "index.md"
    assert index.exists()
    text = index.read_text(encoding="utf-8")
    assert "EGFR" in text


# ---------------------------------------------------------------------------
# Task 12 — Prior knowledge query
# ---------------------------------------------------------------------------


def test_query_prior_knowledge_found(tmp_path):
    """Prior knowledge returns formatted text for known variants."""
    from scripts.clinical_board.knowledge_base import KnowledgeBase
    from scripts.clinical_board.kb_query import query_prior_knowledge

    db_path = tmp_path / "kb.sqlite3"
    build_kb_db(str(db_path))
    kb = KnowledgeBase(str(db_path), str(tmp_path))
    kb.save_decision(
        sample_id="S001",
        mode="cancer",
        gene="EGFR",
        variant="chr7:55259515:T:G",
        board_diagnosis="TKI sensitive",
        board_confidence="high",
        agent_consensus="unanimous",
    )
    result = query_prior_knowledge(
        str(db_path), ["chr7:55259515:T:G"], "cancer"
    )
    assert "PRIOR BOARD KNOWLEDGE" in result
    assert "TKI sensitive" in result
    assert "참고 자료" in result or "reference only" in result.lower()


def test_query_prior_knowledge_empty_kb(tmp_path):
    """Empty KB returns empty string."""
    from scripts.clinical_board.kb_query import query_prior_knowledge

    db_path = tmp_path / "kb.sqlite3"
    build_kb_db(str(db_path))
    result = query_prior_knowledge(
        str(db_path), ["chr7:55259515:T:G"], "cancer"
    )
    assert result == ""


def test_query_prior_knowledge_no_db():
    """Missing DB file returns empty string (no crash)."""
    from scripts.clinical_board.kb_query import query_prior_knowledge

    result = query_prior_knowledge(
        "/nonexistent/kb.sqlite3", ["chr7:55259515:T:G"], "cancer"
    )
    assert result == ""


def test_query_prior_knowledge_mode_isolated(tmp_path):
    """Cancer entries do not leak into rare-disease queries."""
    from scripts.clinical_board.knowledge_base import KnowledgeBase
    from scripts.clinical_board.kb_query import query_prior_knowledge

    db_path = tmp_path / "kb.sqlite3"
    build_kb_db(str(db_path))
    kb = KnowledgeBase(str(db_path), str(tmp_path))
    kb.save_decision(
        sample_id="S001",
        mode="cancer",
        gene="EGFR",
        variant="chr7:55259515:T:G",
        board_diagnosis="TKI sensitive",
        board_confidence="high",
        agent_consensus="unanimous",
    )
    result = query_prior_knowledge(
        str(db_path), ["chr7:55259515:T:G"], "rare-disease"
    )
    assert result == ""


def test_query_prior_knowledge_truncates(tmp_path):
    """max_chars truncation is respected."""
    from scripts.clinical_board.knowledge_base import KnowledgeBase
    from scripts.clinical_board.kb_query import query_prior_knowledge

    db_path = tmp_path / "kb.sqlite3"
    build_kb_db(str(db_path))
    kb = KnowledgeBase(str(db_path), str(tmp_path))
    variants = [f"chr1:{100 + i}:A:T" for i in range(50)]
    for idx, v in enumerate(variants):
        kb.save_decision(
            sample_id=f"S{idx:03d}",
            mode="cancer",
            gene="EGFR",
            variant=v,
            board_diagnosis="TKI sensitive with a fairly long description",
            board_confidence="high",
            agent_consensus="unanimous",
        )
    result = query_prior_knowledge(
        str(db_path), variants, "cancer", max_chars=500
    )
    assert len(result) <= 500
