"""Prior knowledge retrieval for AI Clinical Board v2.

Queries the Knowledge Base SQLite `variant_stats` view for previously
seen variants and renders a compact block of prior board outcomes to
pass into agents as reference-only context. A bilingual anchoring-bias
warning is always included so the agents don't over-weight history.
"""
from __future__ import annotations

import os
import sqlite3
from typing import Iterable


_HEADER = "== PRIOR BOARD KNOWLEDGE =="
_ANCHOR_WARNING = (
    "⚠️  참고 자료 (reference only) — 과거 보드 판단은 앵커링 편향을 피하기 위해 "
    "절대적 기준이 아니라 참고용으로만 사용하십시오. "
    "Historical board calls are reference only; do not anchor on them."
)


def query_prior_knowledge(
    db_path: str,
    variant_ids: Iterable[str],
    mode: str,
    max_chars: int = 4000,
) -> str:
    """Return a formatted prior-knowledge block or "" if nothing to show.

    - Returns "" if the DB file doesn't exist (graceful degradation).
    - Returns "" if the `variant_stats` view holds no rows for any of the
      requested variants in the given mode.
    - Output is truncated to `max_chars` characters.
    """
    if not db_path or not os.path.exists(db_path):
        return ""

    variant_list = [v for v in variant_ids if v]
    if not variant_list:
        return ""

    try:
        conn = sqlite3.connect(db_path)
    except sqlite3.Error:
        return ""
    conn.row_factory = sqlite3.Row

    try:
        placeholders = ", ".join("?" for _ in variant_list)
        rows = conn.execute(
            "SELECT gene, variant, mode, total_cases, high_confidence_count, "
            "diagnoses_seen, last_seen "
            f"FROM variant_stats WHERE mode = ? AND variant IN ({placeholders}) "
            "ORDER BY total_cases DESC, gene",
            (mode, *variant_list),
        ).fetchall()
    except sqlite3.Error:
        conn.close()
        return ""
    finally:
        try:
            conn.close()
        except sqlite3.Error:
            pass

    if not rows:
        return ""

    lines = [_HEADER, ""]
    for row in rows:
        diagnoses = (row["diagnoses_seen"] or "").replace("|", "/")
        lines.append(
            f"- {row['gene']} {row['variant']} "
            f"(mode={row['mode']}): "
            f"total_cases={row['total_cases']}, "
            f"high_confidence={row['high_confidence_count']}, "
            f"diagnoses=[{diagnoses}], "
            f"last_seen={row['last_seen'] or '-'}"
        )
    lines.extend(["", _ANCHOR_WARNING, ""])
    text = "\n".join(lines)

    if len(text) > max_chars:
        text = text[: max(0, max_chars)]
    return text
