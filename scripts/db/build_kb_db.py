"""Initialize Knowledge Base SQLite schema.

Creates the board_decisions table and variant_stats view used by the
AI Clinical Board v2 to persist and aggregate per-variant board outcomes.

Usage:
    python scripts/db/build_kb_db.py [db_path]
"""
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_DB_PATH = "data/knowledge_base/kb.sqlite3"


def build_kb_db(db_path: str = DEFAULT_DB_PATH) -> dict:
    """Initialize (or upgrade) the Knowledge Base SQLite schema.

    Idempotent: safe to call repeatedly on an existing DB.

    Returns dict with build statistics.
    """
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS board_decisions (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id TEXT NOT NULL,
            mode TEXT NOT NULL,
            date TEXT NOT NULL,
            gene TEXT NOT NULL,
            variant TEXT NOT NULL,
            hgvsc TEXT,
            hgvsp TEXT,
            classification TEXT,
            board_diagnosis TEXT,
            board_confidence TEXT,
            clinical_context_summary TEXT,
            agent_consensus TEXT,
            raw_opinion_json TEXT
        )
    """)

    cursor.execute(
        "CREATE INDEX IF NOT EXISTS idx_bd_variant "
        "ON board_decisions(gene, variant, mode)"
    )
    cursor.execute(
        "CREATE INDEX IF NOT EXISTS idx_bd_gene "
        "ON board_decisions(gene, mode)"
    )
    cursor.execute(
        "CREATE INDEX IF NOT EXISTS idx_bd_date "
        "ON board_decisions(date)"
    )
    cursor.execute(
        "CREATE INDEX IF NOT EXISTS idx_bd_hgvsp "
        "ON board_decisions(gene, hgvsp, mode)"
    )

    cursor.execute("DROP VIEW IF EXISTS variant_stats")
    cursor.execute("""
        CREATE VIEW variant_stats AS
        SELECT
            gene,
            variant,
            mode,
            COUNT(*) AS total_cases,
            SUM(CASE WHEN board_confidence = 'high' THEN 1 ELSE 0 END)
                AS high_confidence_count,
            GROUP_CONCAT(DISTINCT board_diagnosis) AS diagnoses_seen,
            MAX(date) AS last_seen
        FROM board_decisions
        GROUP BY gene, variant, mode
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT
        )
    """)
    cursor.execute(
        "INSERT OR REPLACE INTO metadata VALUES ('schema', 'board_decisions_v1')"
    )
    cursor.execute(
        "INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)",
        (datetime.now(timezone.utc).isoformat(),),
    )

    conn.commit()
    count = cursor.execute(
        "SELECT COUNT(*) FROM board_decisions"
    ).fetchone()[0]
    conn.close()

    return {"db_path": db_path, "records": count}


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_DB_PATH
    stats = build_kb_db(target)
    print(f"Built {stats['db_path']}: {stats['records']} existing records")
