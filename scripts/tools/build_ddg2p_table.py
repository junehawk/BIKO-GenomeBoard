"""Build the DDG2P neurodevelopmental gene panel JSON from EBI Gene2Phenotype.

v2.3-T7 — replaces the v1 hand-curated 30-gene starter set
(``data/ddg2p_neurodev_genes.json``) with the full DDG2P panel from EBI's
public FTP archive (CC0 license).

The build script:

1. Resolves the latest DDG2P snapshot under
   ``https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/``
   by listing the directory index and picking the most recent ``YYYY_MM_DD``
   subfolder that contains a ``DDG2P_YYYY-MM-DD.csv.gz`` file. Several legacy
   URL aliases under ``www.ebi.ac.uk/gene2phenotype/downloads/`` are also tried
   first, in case a future build is pinned to a stable alias.
2. Streams the gzipped CSV, normalises gene records, collapses multiple disease
   rows per gene into a single record (taking the highest confidence and
   joining disease names), and filters to the admission set
   (``definitive`` / ``strong`` / ``moderate``).
3. Writes the result to ``data/ddg2p_neurodev_genes.json`` using exactly the
   schema the existing ``scripts/common/ddg2p_panel.py`` loader expects, plus
   the metadata fields the ``version_manager`` reads.
4. If the network or the EBI archive is unreachable, the existing
   ``data/ddg2p_neurodev_genes.json`` is left untouched and the script exits
   normally after recording the failure to
   ``_workspace/v23-engineering/t7_blockers.md`` and printing a summary.

Usage::

    python scripts/tools/build_ddg2p_table.py            # full build
    python scripts/tools/build_ddg2p_table.py --dry-run  # probe URL only
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import logging
import os
import re
import sys
import urllib.error
import urllib.request
from datetime import date, datetime, timezone
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PANEL_JSON_PATH = os.path.join(REPO_ROOT, "data", "ddg2p_neurodev_genes.json")
BLOCKER_PATH = os.path.join(REPO_ROOT, "_workspace", "v23-engineering", "t7_blockers.md")

# Admission set (spec Q2 — keep aligned with v1 JSON and the panel loader).
ADMITTED_CONFIDENCES: Tuple[str, ...] = ("definitive", "strong", "moderate")

# Confidence ranking (highest first) — used to collapse duplicate gene rows.
CONFIDENCE_RANK: Dict[str, int] = {
    "definitive": 5,
    "strong": 4,
    "moderate": 3,
    "limited": 2,
    "disputed": 1,
    "refuted": 0,
}

# Stable EBI FTP archive root (used for both directory listing and download).
EBI_FTP_ROOT = "https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/"

# Optional alias URLs that may exist on www.ebi.ac.uk in the future. These are
# tried first; on 404 / non-CSV response we fall back to the FTP archive scan.
ALIAS_URL_CANDIDATES: Tuple[str, ...] = (
    "https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz",
    "https://www.ebi.ac.uk/gene2phenotype/downloads/G2P_DD_panel.csv.gz",
)

USER_AGENT = "BIKO-GenomeBoard/2.3 (build_ddg2p_table.py; +https://www.ebi.ac.uk/gene2phenotype/)"

REQUEST_TIMEOUT = 60  # seconds


# ---------------------------------------------------------------------------
# URL discovery
# ---------------------------------------------------------------------------


def _http_get(url: str, timeout: int = REQUEST_TIMEOUT) -> bytes:
    """GET ``url`` and return the response body as bytes.

    Raises :class:`urllib.error.URLError` / :class:`urllib.error.HTTPError`
    on any failure. The caller is expected to treat those as "EBI unreachable".
    """
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _looks_like_csv_gz(payload: bytes) -> bool:
    """Return True if ``payload`` is plausibly a gzip-compressed CSV.

    EBI's www.ebi.ac.uk SPA returns a 1.1 KB HTML stub for unknown paths even
    when the URL ends in ``.csv.gz`` — we use the gzip magic bytes plus a
    minimum size threshold to detect that case.
    """
    if len(payload) < 1024:
        return False
    return payload[:2] == b"\x1f\x8b"


def _list_ftp_snapshots(index_html: str) -> List[str]:
    """Return ``YYYY_MM_DD`` snapshot directory names sorted newest-first."""
    pattern = re.compile(r'href="(\d{4}_\d{2}_\d{2})/"')
    snaps = pattern.findall(index_html)
    return sorted(set(snaps), reverse=True)


def _list_snapshot_files(index_html: str) -> List[str]:
    pattern = re.compile(r'href="([^"]+\.csv\.gz)"')
    return pattern.findall(index_html)


def discover_ddg2p_url() -> Tuple[str, bytes]:
    """Locate the latest DDG2P CSV.gz and return ``(url, payload)``.

    Strategy:

    1. Try each :data:`ALIAS_URL_CANDIDATES` in order; if any returns a valid
       gzip payload, use it.
    2. Otherwise scan :data:`EBI_FTP_ROOT` for the newest snapshot directory
       containing a ``DDG2P_*.csv.gz`` and download that file.

    Raises :class:`RuntimeError` if no candidate succeeds.
    """
    last_error: Optional[Exception] = None

    for url in ALIAS_URL_CANDIDATES:
        try:
            payload = _http_get(url)
        except (urllib.error.URLError, urllib.error.HTTPError) as exc:
            last_error = exc
            logger.info("alias %s unreachable: %s", url, exc)
            continue
        if _looks_like_csv_gz(payload):
            logger.info("alias %s returned a valid gzip payload (%d bytes)", url, len(payload))
            return url, payload
        logger.info("alias %s returned non-CSV payload (%d bytes) — skipping", url, len(payload))

    # Fallback: scan FTP archive
    try:
        index = _http_get(EBI_FTP_ROOT).decode("utf-8", errors="replace")
    except (urllib.error.URLError, urllib.error.HTTPError) as exc:
        raise RuntimeError(f"could not reach EBI FTP archive at {EBI_FTP_ROOT}: {exc}") from exc

    snapshots = _list_ftp_snapshots(index)
    if not snapshots:
        raise RuntimeError(
            f"no YYYY_MM_DD snapshot directories found at {EBI_FTP_ROOT} (last alias error: {last_error})"
        )

    for snap in snapshots:
        snap_url = f"{EBI_FTP_ROOT}{snap}/"
        try:
            snap_index = _http_get(snap_url).decode("utf-8", errors="replace")
        except (urllib.error.URLError, urllib.error.HTTPError) as exc:
            logger.info("snapshot %s unreachable: %s", snap_url, exc)
            continue
        files = _list_snapshot_files(snap_index)
        ddg2p_files = [f for f in files if f.startswith("DDG2P_") and f.endswith(".csv.gz")]
        if not ddg2p_files:
            continue
        ddg2p_files.sort(reverse=True)
        target_url = f"{snap_url}{ddg2p_files[0]}"
        try:
            payload = _http_get(target_url)
        except (urllib.error.URLError, urllib.error.HTTPError) as exc:
            logger.info("download %s failed: %s", target_url, exc)
            continue
        if _looks_like_csv_gz(payload):
            return target_url, payload
        logger.info("file %s is not a valid gzip — skipping", target_url)

    raise RuntimeError(f"could not find a usable DDG2P_*.csv.gz under {EBI_FTP_ROOT} (last alias error: {last_error})")


# ---------------------------------------------------------------------------
# CSV parsing & normalisation
# ---------------------------------------------------------------------------


def _better_confidence(current: Optional[str], candidate: Optional[str]) -> str:
    """Return the higher-ranked of two confidence strings."""
    cur = (current or "").lower()
    cand = (candidate or "").lower()
    if not cur:
        return cand
    if not cand:
        return cur
    return cand if CONFIDENCE_RANK.get(cand, -1) > CONFIDENCE_RANK.get(cur, -1) else cur


def parse_ddg2p_csv(payload: bytes) -> Dict[str, dict]:
    """Parse a DDG2P CSV.gz payload into ``{gene: record}``.

    Multiple disease rows for the same gene are collapsed:

    - ``confidence`` keeps the highest rank seen.
    - ``disease`` becomes a ``"; "``-joined string of unique disease names.
    - ``allelic_requirement`` is the value associated with the highest-confidence row.
    - ``mutation_consequence`` likewise.
    - ``pmids`` is the union of all publications across rows.
    """
    with gzip.open(io.BytesIO(payload), mode="rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise RuntimeError("DDG2P CSV has no header row")
        # Normalise expected columns — the spec lowercases & space-separates them.
        expected = {"gene symbol", "disease name", "allelic requirement", "confidence"}
        missing = expected - set(reader.fieldnames)
        if missing:
            raise RuntimeError(f"DDG2P CSV missing required columns: {sorted(missing)} (found: {reader.fieldnames})")

        records: Dict[str, dict] = {}
        for row in reader:
            gene = (row.get("gene symbol") or "").strip()
            if not gene:
                continue
            confidence = (row.get("confidence") or "").strip().lower()
            if not confidence:
                continue

            disease = (row.get("disease name") or "").strip()
            allelic = (row.get("allelic requirement") or "").strip()
            mech = (row.get("variant consequence") or row.get("molecular mechanism") or "").strip()
            pubs = (row.get("publications") or "").strip()
            pmids = [p.strip() for p in re.split(r"[;,\s]+", pubs) if p.strip().isdigit()]

            existing = records.get(gene)
            if existing is None:
                records[gene] = {
                    "confidence": confidence,
                    "disease": [disease] if disease else [],
                    "allelic_requirement": allelic,
                    "mutation_consequence": mech,
                    "pmids": list(pmids),
                }
                continue

            new_conf = _better_confidence(existing["confidence"], confidence)
            promoted = new_conf == confidence and confidence != existing["confidence"]
            existing["confidence"] = new_conf
            if disease and disease not in existing["disease"]:
                existing["disease"].append(disease)
            if promoted:
                # Adopt the allelic / mechanism of the now top-ranked row.
                if allelic:
                    existing["allelic_requirement"] = allelic
                if mech:
                    existing["mutation_consequence"] = mech
            else:
                if not existing["allelic_requirement"] and allelic:
                    existing["allelic_requirement"] = allelic
                if not existing["mutation_consequence"] and mech:
                    existing["mutation_consequence"] = mech
            for pmid in pmids:
                if pmid not in existing["pmids"]:
                    existing["pmids"].append(pmid)

        # Finalise: join disease lists into the "; "-joined form the v1 schema uses.
        for rec in records.values():
            rec["disease"] = "; ".join(rec["disease"]) if rec["disease"] else ""

        return records


def filter_admitted(records: Dict[str, dict]) -> Dict[str, dict]:
    return {g: r for g, r in records.items() if (r.get("confidence") or "").lower() in ADMITTED_CONFIDENCES}


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------


def build_panel_json(
    *,
    source_url: str,
    payload: bytes,
    admitted_records: Dict[str, dict],
) -> dict:
    today = date.today().isoformat()
    return {
        "source": "Gene2Phenotype DDG2P",
        "license": "CC0",
        "build_date": today,
        "source_url": source_url,
        "source_download_size": f"{len(payload)} bytes",
        "source_md5": hashlib.md5(payload).hexdigest(),
        "admission_confidences": list(ADMITTED_CONFIDENCES),
        "build_note": (
            "Full DDG2P ingest (v2.3-T7) from EBI Gene2Phenotype FTP archive. "
            "Replaces the v1 hand-curated 30-gene starter set. License CC0. "
            "Admission set restricted to definitive/strong/moderate confidence "
            "per spec Q2; limited/disputed/refuted records are excluded."
        ),
        "record_count": len(admitted_records),
        "genes": admitted_records,
    }


def write_panel_json(panel: dict, path: str = PANEL_JSON_PATH) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(panel, fh, indent=2, ensure_ascii=False, sort_keys=False)
        fh.write("\n")
    os.replace(tmp, path)


def write_blocker(reason: str) -> None:
    """Append a fallback note to t7_blockers.md so the lead can see why."""
    os.makedirs(os.path.dirname(BLOCKER_PATH), exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    line = f"- {stamp} build_ddg2p_table fallback: {reason}\n"
    with open(BLOCKER_PATH, "a", encoding="utf-8") as fh:
        fh.write(line)


# ---------------------------------------------------------------------------
# CLI entry-point
# ---------------------------------------------------------------------------


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build the DDG2P neurodevelopmental gene panel JSON from EBI Gene2Phenotype."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Probe the EBI source URL and exit without writing the JSON file.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable INFO-level logging during URL discovery.",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(message)s",
    )

    started = datetime.now()
    print(f"[ddg2p] starting build (dry_run={args.dry_run})")

    try:
        source_url, payload = discover_ddg2p_url()
    except (RuntimeError, urllib.error.URLError) as exc:
        msg = f"DDG2P download failed: {exc}"
        print(f"[ddg2p] {msg}")
        write_blocker(msg)
        # v1 JSON is left in place — script exits normally.
        if os.path.exists(PANEL_JSON_PATH):
            print(f"[ddg2p] keeping existing panel at {PANEL_JSON_PATH}")
        else:
            print(f"[ddg2p] WARNING: no existing panel at {PANEL_JSON_PATH}")
        return 0

    print(f"[ddg2p] source: {source_url}")
    print(f"[ddg2p] download size: {len(payload)} bytes")
    print(f"[ddg2p] md5: {hashlib.md5(payload).hexdigest()}")

    if args.dry_run:
        print("[ddg2p] --dry-run: not writing JSON")
        return 0

    try:
        all_records = parse_ddg2p_csv(payload)
    except RuntimeError as exc:
        msg = f"DDG2P parse failed: {exc}"
        print(f"[ddg2p] {msg}")
        write_blocker(msg)
        return 0

    admitted = filter_admitted(all_records)
    panel = build_panel_json(
        source_url=source_url,
        payload=payload,
        admitted_records=admitted,
    )
    write_panel_json(panel)

    elapsed = (datetime.now() - started).total_seconds()
    print(
        f"[ddg2p] wrote {PANEL_JSON_PATH}: "
        f"{len(admitted)} admitted genes "
        f"(out of {len(all_records)} unique genes in DDG2P) "
        f"in {elapsed:.1f}s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
