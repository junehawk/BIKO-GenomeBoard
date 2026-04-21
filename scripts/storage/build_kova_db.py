#!/usr/bin/env python3
"""Build the KOVA v7 allele-frequency SQLite database.

KOVA v7 ships as a single gzipped TSV (``1_KOVA.v7.tsv.gz`` — ~3.3 GB,
~43.3 M rows) whose schema is described in the KOVA consortium release
notes. This builder streams the TSV once, extracts the Korean cohort
frequency columns (``chrom``, ``pos``, ``ref_allele``, ``alt_allele``,
``kova_AF``, ``kova_AC``, ``kova_AN``, ``KOVA_homozygote_count``), and
writes them to a compact SQLite table at ``data/db/kova.sqlite3``. No
other KOVA columns are retained — callers that need CADD / gnomAD / VEP
annotations read those from their upstream sources.

Design notes:

* **Streaming** — the TSV is consumed with ``gzip.open`` + ``csv.DictReader``
  so memory stays constant regardless of input size.
* **Idempotent** — ``DROP TABLE IF EXISTS`` before create; the metadata
  table is upserted via ``INSERT OR REPLACE``. Re-running on the same
  input produces a byte-comparable DB (ignoring the ``build_date``
  metadata row).
* **Batch insert** — rows are staged into 50K-element lists and flushed
  via ``executemany`` inside a single transaction. ``journal_mode=OFF``
  and ``synchronous=OFF`` are used during the build for ~10× throughput;
  both are restored (WAL) and the file is ``VACUUM``-ed before exit.
* **Row filtering** — any row whose ``kova_AF`` is missing (``""`` /
  ``NA`` / ``.``) is skipped. ``kova_AC`` / ``kova_AN`` / homozygote
  are allowed to be NULL independently.
* **Version metadata** — a ``metadata`` key/value table records the
  source URL, build date, record count, KOVA release version and
  assembly so ``version_manager.get_all_db_versions`` can surface
  provenance in the report.

Usage::

    python scripts/storage/build_kova_db.py \\
        --input /Users/JL/Downloads/1_KOVA.v7.tsv.gz \\
        --output data/db/kova.sqlite3

Expected wall time: 30–60 min on a 43.3 M-row input (SSD, Apple Silicon).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import sqlite3
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

logger = logging.getLogger(__name__)

DEFAULT_DB_PATH = "data/db/kova.sqlite3"
DEFAULT_SOURCE_URL = "https://www.kobic.re.kr/kova/"
DEFAULT_KOVA_VERSION = "v7"
DEFAULT_ASSEMBLY = "GRCh38"

_BATCH_SIZE = 50_000
_PROGRESS_EVERY = 1_000_000

_REQUIRED_COLUMNS = (
    "chrom",
    "pos",
    "ref_allele",
    "alt_allele",
    "kova_AF",
    "kova_AC",
    "kova_AN",
    "KOVA_homozygote_count",
)

# Row-level sentinels that indicate a missing value in KOVA TSV cells.
_NULL_TOKENS = frozenset(("", ".", "NA", "nan", "NaN"))


def _coerce_float(value: str) -> Optional[float]:
    if value is None:
        return None
    if value in _NULL_TOKENS:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _coerce_int(value: str) -> Optional[int]:
    if value is None:
        return None
    if value in _NULL_TOKENS:
        return None
    try:
        return int(value)
    except ValueError:
        try:
            return int(float(value))
        except ValueError:
            return None


def _normalize_chrom(chrom: str) -> str:
    if not chrom:
        return chrom
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def _iter_rows(
    input_path: Path,
) -> Iterable[Tuple[str, int, str, str, Optional[float], Optional[int], Optional[int], Optional[int]]]:
    """Yield one ``(chrom, pos, ref, alt, af, ac, an, homo)`` tuple per
    KOVA v7 TSV row. Rows with missing ``kova_AF`` are skipped."""
    open_fn = gzip.open if str(input_path).endswith(".gz") else open
    with open_fn(input_path, mode="rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in _REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(
                f"KOVA TSV {input_path} is missing required columns: {missing}. "
                f"Found columns: {reader.fieldnames[:12] if reader.fieldnames else None}..."
            )
        for raw in reader:
            af = _coerce_float(raw.get("kova_AF", ""))
            if af is None:
                continue
            pos_raw = raw.get("pos", "")
            try:
                pos = int(pos_raw)
            except (TypeError, ValueError):
                continue
            chrom = _normalize_chrom(raw.get("chrom", ""))
            ref = raw.get("ref_allele", "")
            alt = raw.get("alt_allele", "")
            if not chrom or not ref or not alt:
                continue
            yield (
                chrom,
                pos,
                ref,
                alt,
                af,
                _coerce_int(raw.get("kova_AC", "")),
                _coerce_int(raw.get("kova_AN", "")),
                _coerce_int(raw.get("KOVA_homozygote_count", "")),
            )


def _create_schema(conn: sqlite3.Connection) -> None:
    conn.execute("DROP TABLE IF EXISTS kova_af")
    conn.execute(
        """
        CREATE TABLE kova_af (
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            kova_af REAL,
            kova_ac INTEGER,
            kova_an INTEGER,
            kova_homozygote INTEGER,
            PRIMARY KEY(chrom, pos, ref, alt)
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT
        )
        """
    )


def _write_metadata(
    conn: sqlite3.Connection,
    *,
    source_url: str,
    kova_version: str,
    assembly: str,
    record_count: int,
    input_path: Path,
) -> None:
    now = datetime.now(timezone.utc).isoformat()
    rows = [
        ("source", source_url),
        ("kova_version", kova_version),
        ("assembly", assembly),
        ("build_date", now),
        ("record_count", str(record_count)),
        ("input_file", input_path.name),
    ]
    conn.executemany("INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)", rows)


def build_db(
    input_path: str | Path,
    db_path: str | Path = DEFAULT_DB_PATH,
    *,
    source_url: str = DEFAULT_SOURCE_URL,
    kova_version: str = DEFAULT_KOVA_VERSION,
    assembly: str = DEFAULT_ASSEMBLY,
) -> str:
    """Build a KOVA SQLite DB at ``db_path`` from the KOVA v7 TSV.

    Parameters
    ----------
    input_path:
        Path to ``1_KOVA.v7.tsv.gz`` (or the uncompressed TSV).
    db_path:
        Target SQLite file. Parent directories are created if missing.
    source_url / kova_version / assembly:
        Recorded in the ``metadata`` table for ``version_manager``.

    Returns
    -------
    str
        The target ``db_path`` (for pipeline-style chaining).
    """
    input_p = Path(input_path)
    db_p = Path(db_path)
    if not input_p.exists():
        raise FileNotFoundError(f"KOVA input TSV not found: {input_p}")
    db_p.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    conn = sqlite3.connect(str(db_p))
    try:
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        conn.execute("PRAGMA temp_store=MEMORY")
        conn.execute("PRAGMA cache_size=-200000")  # ~200 MB page cache
        _create_schema(conn)

        insert_sql = (
            "INSERT OR REPLACE INTO kova_af "
            "(chrom, pos, ref, alt, kova_af, kova_ac, kova_an, kova_homozygote) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
        )

        batch: List[Tuple] = []
        total = 0
        conn.execute("BEGIN")
        for row in _iter_rows(input_p):
            batch.append(row)
            if len(batch) >= _BATCH_SIZE:
                conn.executemany(insert_sql, batch)
                total += len(batch)
                batch.clear()
                if total % _PROGRESS_EVERY == 0:
                    elapsed = time.time() - t0
                    rate = total / elapsed if elapsed else 0
                    logger.info("  %s rows inserted (%.0f rows/s, %.1f min elapsed)", f"{total:,}", rate, elapsed / 60)
        if batch:
            conn.executemany(insert_sql, batch)
            total += len(batch)
            batch.clear()
        conn.execute("COMMIT")

        logger.info("Creating indexes ...")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_kova_pos ON kova_af(chrom, pos)")

        _write_metadata(
            conn,
            source_url=source_url,
            kova_version=kova_version,
            assembly=assembly,
            record_count=total,
            input_path=input_p,
        )
        conn.commit()
    finally:
        conn.close()

    # Post-load compaction — switch back to WAL and shrink. VACUUM must run
    # outside of any open transaction, hence the fresh connection.
    # ``sqlite3.connect`` as a context manager commits/rollbacks on exit but
    # does NOT close the connection (CPython quirk), so close explicitly to
    # release the file lock for subsequent callers (notably the idempotent
    # rebuild path).
    logger.info("VACUUMing %s ...", db_p)
    vac = sqlite3.connect(str(db_p))
    try:
        vac.execute("PRAGMA journal_mode=WAL")
        vac.execute("VACUUM")
    finally:
        vac.close()

    elapsed = time.time() - t0
    size_mb = db_p.stat().st_size / (1024 * 1024)
    logger.info(
        "KOVA DB built: %s rows → %s (%.1f MB) in %.1f min",
        f"{total:,}",
        db_p,
        size_mb,
        elapsed / 60,
    )
    return str(db_p)


def _parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build KOVA v7 allele-frequency SQLite DB from 1_KOVA.v7.tsv.gz",
    )
    parser.add_argument("--input", "-i", required=True, help="Path to 1_KOVA.v7.tsv.gz (or uncompressed TSV)")
    parser.add_argument(
        "--output", "-o", default=DEFAULT_DB_PATH, help=f"Target SQLite file (default {DEFAULT_DB_PATH})"
    )
    parser.add_argument("--source-url", default=DEFAULT_SOURCE_URL, help="URL recorded in metadata.source")
    parser.add_argument("--kova-version", default=DEFAULT_KOVA_VERSION, help="KOVA release version (default v7)")
    parser.add_argument("--assembly", default=DEFAULT_ASSEMBLY, help="Reference assembly (default GRCh38)")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args = _parse_args(argv)
    try:
        build_db(
            args.input,
            args.output,
            source_url=args.source_url,
            kova_version=args.kova_version,
            assembly=args.assembly,
        )
    except FileNotFoundError as e:
        logger.error("%s", e)
        return 2
    except ValueError as e:
        logger.error("Input schema error: %s", e)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
