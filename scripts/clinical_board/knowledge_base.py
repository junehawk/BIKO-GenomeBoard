"""Knowledge Base storage for AI Clinical Board v2.

Persists board decisions into a SQLite database and mirrors per-gene
summaries into a markdown wiki tree (`{wiki_dir}/{mode}/genes/{gene}.md`)
so clinicians can browse historical board outcomes.
"""
from __future__ import annotations

import json
import os
import sqlite3
from datetime import date as _date
from pathlib import Path
from typing import Any


_ALLOWED_FIELDS = (
    "sample_id",
    "mode",
    "date",
    "gene",
    "variant",
    "hgvsc",
    "hgvsp",
    "classification",
    "board_diagnosis",
    "board_confidence",
    "clinical_context_summary",
    "agent_consensus",
    "raw_opinion_json",
)

_REQUIRED_FIELDS = ("sample_id", "mode", "gene", "variant")


class KnowledgeBase:
    """SQLite-backed store of board decisions with a markdown wiki mirror."""

    def __init__(self, db_path: str, wiki_dir: str):
        self.db_path = db_path
        self.wiki_dir = Path(wiki_dir)
        self.wiki_dir.mkdir(parents=True, exist_ok=True)

    def _connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def save_decision(self, **fields: Any) -> int:
        """Insert a single board decision. Returns the new row id."""
        for required in _REQUIRED_FIELDS:
            if not fields.get(required):
                raise ValueError(f"save_decision missing required field: {required}")

        row = {name: fields.get(name) for name in _ALLOWED_FIELDS}
        if not row.get("date"):
            row["date"] = _date.today().isoformat()

        if isinstance(row.get("raw_opinion_json"), (dict, list)):
            row["raw_opinion_json"] = json.dumps(row["raw_opinion_json"], ensure_ascii=False)

        cols = ", ".join(_ALLOWED_FIELDS)
        placeholders = ", ".join("?" for _ in _ALLOWED_FIELDS)
        values = tuple(row[name] for name in _ALLOWED_FIELDS)

        conn = self._connect()
        try:
            cursor = conn.execute(
                f"INSERT INTO board_decisions ({cols}) VALUES ({placeholders})",
                values,
            )
            conn.commit()
            return cursor.lastrowid or 0
        finally:
            conn.close()

    def generate_gene_wiki(self, gene: str, mode: str) -> Path:
        """Regenerate the markdown wiki page for a (gene, mode) pair."""
        conn = self._connect()
        try:
            stats_rows = conn.execute(
                "SELECT variant, total_cases, high_confidence_count, "
                "diagnoses_seen, last_seen "
                "FROM variant_stats WHERE gene = ? AND mode = ? "
                "ORDER BY total_cases DESC, variant",
                (gene, mode),
            ).fetchall()
            decisions = conn.execute(
                "SELECT date, sample_id, variant, hgvsp, classification, "
                "board_diagnosis, board_confidence, agent_consensus "
                "FROM board_decisions WHERE gene = ? AND mode = ? "
                "ORDER BY date DESC, id DESC",
                (gene, mode),
            ).fetchall()
        finally:
            conn.close()

        out_dir = self.wiki_dir / mode / "genes"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{gene}.md"

        lines: list[str] = [
            f"# {gene} — {mode}",
            "",
            f"_Total recorded decisions: {len(decisions)}_",
            "",
            "## Variant summary",
            "",
        ]
        if stats_rows:
            lines.append("| Variant | Cases | High-confidence | Diagnoses | Last seen |")
            lines.append("|---|---|---|---|---|")
            for row in stats_rows:
                diagnoses = (row["diagnoses_seen"] or "").replace("|", "/")
                lines.append(
                    f"| {row['variant']} | {row['total_cases']} | "
                    f"{row['high_confidence_count']} | {diagnoses} | "
                    f"{row['last_seen'] or ''} |"
                )
        else:
            lines.append("_No decisions recorded yet._")

        lines.extend(["", "## Decisions", ""])
        if decisions:
            for row in decisions:
                lines.append(
                    f"- **{row['date']}** · {row['sample_id']} · "
                    f"{row['variant']} "
                    f"({row['hgvsp'] or '-'}) · "
                    f"{row['classification'] or '-'} → "
                    f"**{row['board_diagnosis'] or '-'}** "
                    f"(confidence: {row['board_confidence'] or '-'}, "
                    f"consensus: {row['agent_consensus'] or '-'})"
                )
        else:
            lines.append("_No decisions recorded yet._")
        lines.append("")

        out_path.write_text("\n".join(lines), encoding="utf-8")
        return out_path

    def update_log(self, sample_id: str, mode: str, date: str | None = None) -> Path:
        """Append a one-line entry to `{wiki_dir}/log.md`."""
        entry_date = date or _date.today().isoformat()
        log_path = self.wiki_dir / "log.md"
        header_needed = not log_path.exists()
        with log_path.open("a", encoding="utf-8") as fh:
            if header_needed:
                fh.write("# Board decision log\n\n")
            fh.write(f"- {entry_date} · {mode} · {sample_id}\n")
        return log_path

    def update_index(self) -> Path:
        """Regenerate `{wiki_dir}/index.md` listing all gene wiki files."""
        entries: list[tuple[str, str, Path]] = []
        for mode_dir in sorted(self.wiki_dir.iterdir()):
            if not mode_dir.is_dir():
                continue
            genes_dir = mode_dir / "genes"
            if not genes_dir.is_dir():
                continue
            for md in sorted(genes_dir.glob("*.md")):
                entries.append((mode_dir.name, md.stem, md))

        lines = ["# Knowledge Base Index", ""]
        if not entries:
            lines.append("_No gene wiki pages yet._")
        else:
            current_mode: str | None = None
            for mode_name, gene_name, md_path in entries:
                if mode_name != current_mode:
                    lines.extend(["", f"## {mode_name}", ""])
                    current_mode = mode_name
                rel = os.path.relpath(md_path, self.wiki_dir)
                lines.append(f"- [{gene_name}]({rel})")
        lines.append("")

        index_path = self.wiki_dir / "index.md"
        index_path.write_text("\n".join(lines), encoding="utf-8")
        return index_path
