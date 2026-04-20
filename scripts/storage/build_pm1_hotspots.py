#!/usr/bin/env python3
"""Build the ACMG PM1 hotspot table — `data/pm1_hotspot_domains.json`.

This file is the **single ACMG-authoritative source** for `PM1` evidence in the
classification engine. It supersedes `civic.sqlite3.hotspots` and
`data/db/hotspots/hotspots_v2_single.tsv`, which remain in place for Tier-III
selector / downstream signalling but **must not** be consulted by
`evidence_collector._pm1_hotspot_match`.

Merge layers (lowest → highest precedence):

1. **cancerhotspots v2 single-residue TSV** — `data/db/hotspots/hotspots_v2_single.tsv`.
   Parsed opportunistically: if the file is missing or shorter than the header
   row (v2.2 repo ships a 404 stub until the upstream download is re-wired),
   this layer contributes zero rows and the build proceeds.
2. **ClinGen VCEP / AMP 2017 baseline** — hard-coded constant in this script,
   mirroring the v2.2 clinical-review Table §1.1. This is the minimum set every
   BIKO deployment ships with.
3. **clinical-advisor overrides** — `data/pm1_hotspot_overrides.yaml`, authored
   by clinical-advisor and copied from `_workspace/v22-phaseA/artifacts/`. This
   layer is the clinical lock — `locked: true` entries cannot be downgraded by
   a future cancerhotspots re-merge.

Output (`data/pm1_hotspot_domains.json`):

```json
{
  "$schema": "https://biko/schema/pm1_hotspot_domains.v1.json",
  "version": "1.0",
  "build_date": "2026-04-14",
  "source_refs": ["cancerhotspots_v2", "ClinGen_SVI_TP53_v2.0",
                  "Giacomelli_2018_PMID_30224644", "AMP_ASCO_CAP_2017_PMID_27993330"],
  "source_hash": "sha256:<hex>",
  "record_count": <int>,
  "genes": {
    "TP53": [
      {"range": [175, 175], "domain": "DBD β-sandwich (R175 structural hotspot)",
       "strength": "moderate", "source": "PMID 30224644",
       "rationale": "...", "locked": true}
    ],
    ...
  }
}
```

Invariants (enforced by this script and by `tests/test_build_pm1_hotspots.py`):

- Every entry has a non-empty `source` matching `r"PMID \\d+"`.
- `strength` ∈ {"moderate", "supporting"}.
- `range` is a 2-element integer list with `start <= end`.
- `record_count` equals the total number of entries across all genes.
- `source_hash` is deterministic over the `genes` payload.
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import re
from datetime import date
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger(__name__)


REPO_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_OVERRIDES_PATH = REPO_ROOT / "data" / "pm1_hotspot_overrides.yaml"
DEFAULT_CANCERHOTSPOTS_TSV = REPO_ROOT / "data" / "db" / "hotspots" / "hotspots_v2_single.tsv"
DEFAULT_OUTPUT_PATH = REPO_ROOT / "data" / "pm1_hotspot_domains.json"

PMID_PATTERN = re.compile(r"PMID\s+\d+")
ALLOWED_STRENGTHS = {"moderate", "supporting"}

# v2.2 source reference block — human-readable anchor list for the JSON
# top-level `source_refs`. Individual entries carry their own PMID in
# the `source` field.
SOURCE_REFS_BASE = [
    "pm1_hotspot_overrides.yaml_v1.0_2026-04-14",
    "ClinGen_TP53_VCEP_v2.0",
    "Giacomelli_2018_PMID_30224644",
    "AMP_ASCO_CAP_2017_PMID_27993330",
    "Garrett_2021_PMID_33280026",
]
SOURCE_REF_CANCERHOTSPOTS = "cancerhotspots_v2_single_PMID_29247016"


# Required coverage for v2.2 — every (gene, residue) pair below MUST be
# satisfied by a range entry in the merged output. Clinical-advisor decision
# 2026-04-14: this set protects against a future edit silently dropping a
# required residue. A missing pair fails the build.
REQUIRED_COVERAGE: dict[str, list[int]] = {
    "TP53": [175, 245, 246, 247, 248, 249, 273, 282],
    "KRAS": [12, 13, 61, 117, 146],
    "NRAS": [12, 13, 61],
    "HRAS": [12, 13, 61],
    "BRAF": [600, 601],
    "PIK3CA": [542, 545, 1047],
    "EGFR": [719, 746, 747, 748, 749, 750, 790, 858, 861],
    "IDH1": [132],
    "IDH2": [140, 172],
}


# Hard-coded ClinGen VCEP / AMP 2017 baseline. Every entry is mirrored in the
# clinical-advisor overrides YAML for v2.2; this table is the fallback used
# when overrides are missing or empty, guaranteeing a known-good minimum set.
CLINGEN_BASELINE: dict[str, list[dict[str, Any]]] = {
    "TP53": [
        {
            "range": [175, 175],
            "domain": "DBD β-sandwich (R175 structural hotspot)",
            "strength": "moderate",
            "source": "PMID 30224644",
        },
        {"range": [245, 249], "domain": "DBD L3 loop (aa 245-254)", "strength": "moderate", "source": "PMID 30224644"},
        {
            "range": [273, 273],
            "domain": "DBD DNA-contact surface (R273)",
            "strength": "moderate",
            "source": "PMID 30224644",
        },
        {
            "range": [282, 282],
            "domain": "DBD L3 loop-adjacent (R282)",
            "strength": "moderate",
            "source": "PMID 30224644",
        },
    ],
    "KRAS": [
        {"range": [12, 13], "domain": "G-domain P-loop / switch I", "strength": "moderate", "source": "PMID 27993330"},
        {"range": [61, 61], "domain": "G-domain switch II", "strength": "moderate", "source": "PMID 27993330"},
        {
            "range": [117, 117],
            "domain": "G-domain nucleotide binding",
            "strength": "supporting",
            "source": "PMID 27993330",
        },
        {
            "range": [146, 146],
            "domain": "G-domain nucleotide binding",
            "strength": "supporting",
            "source": "PMID 27993330",
        },
    ],
    "NRAS": [
        {"range": [12, 13], "domain": "G-domain P-loop / switch I", "strength": "moderate", "source": "PMID 27993330"},
        {"range": [61, 61], "domain": "G-domain switch II", "strength": "moderate", "source": "PMID 27993330"},
    ],
    "HRAS": [
        {"range": [12, 13], "domain": "G-domain P-loop / switch I", "strength": "moderate", "source": "PMID 27993330"},
        {"range": [61, 61], "domain": "G-domain switch II", "strength": "moderate", "source": "PMID 27993330"},
    ],
    "BRAF": [
        {
            "range": [600, 600],
            "domain": "Kinase activation loop (V600)",
            "strength": "moderate",
            "source": "PMID 27993330",
        },
        {
            "range": [601, 601],
            "domain": "Kinase activation loop (adjacent to V600)",
            "strength": "supporting",
            "source": "PMID 27993330",
        },
    ],
    "PIK3CA": [
        {"range": [542, 542], "domain": "Helical domain (E542)", "strength": "moderate", "source": "PMID 27993330"},
        {"range": [545, 545], "domain": "Helical domain (E545)", "strength": "moderate", "source": "PMID 27993330"},
        {"range": [1047, 1047], "domain": "Kinase domain (H1047)", "strength": "moderate", "source": "PMID 27993330"},
    ],
    "EGFR": [
        {"range": [719, 719], "domain": "ATP-binding P-loop (G719)", "strength": "moderate", "source": "PMID 27993330"},
        {"range": [746, 750], "domain": "Exon 19 deletion region", "strength": "moderate", "source": "PMID 27993330"},
        {
            "range": [790, 790],
            "domain": "Gatekeeper T790M (resistance)",
            "strength": "moderate",
            "source": "PMID 27993330",
        },
        {"range": [858, 858], "domain": "Activation loop (L858R)", "strength": "moderate", "source": "PMID 27993330"},
        {"range": [861, 861], "domain": "Activation loop (L861Q)", "strength": "supporting", "source": "PMID 27993330"},
    ],
    "IDH1": [
        {"range": [132, 132], "domain": "Active site (R132)", "strength": "moderate", "source": "PMID 27993330"},
    ],
    "IDH2": [
        {"range": [140, 140], "domain": "Active site (R140)", "strength": "moderate", "source": "PMID 27993330"},
        {"range": [172, 172], "domain": "Active site (R172)", "strength": "supporting", "source": "PMID 27993330"},
    ],
}


# ---------------------------------------------------------------------------
# Layer loaders
# ---------------------------------------------------------------------------


def _load_cancerhotspots(tsv_path: Path) -> dict[str, list[dict[str, Any]]]:
    """Best-effort parse of the cancerhotspots v2 single-residue TSV.

    Returns an empty dict if the file is missing, empty, or does not look
    like a tab-separated table (v2.2 repo ships a 404 stub placeholder).
    """
    entries: dict[str, list[dict[str, Any]]] = {}
    if not tsv_path.exists():
        return entries

    try:
        text = tsv_path.read_text()
    except OSError as e:
        logger.warning("cancerhotspots TSV unreadable (%s): %s — skipping layer", tsv_path, e)
        return entries

    stripped = text.strip()
    if not stripped or "\t" not in stripped.splitlines()[0]:
        logger.warning(
            "cancerhotspots_v2 TSV unavailable (404 stub detected, %d bytes); "
            "falling back to YAML-only PM1 source layer "
            "(clinical-advisor pre-approved for v2.2)",
            len(stripped),
        )
        return entries

    try:
        reader = csv.DictReader(stripped.splitlines(), delimiter="\t")
        for row in reader:
            gene = (row.get("Hugo_Symbol") or row.get("gene") or "").strip()
            pos_raw = (row.get("Amino_Acid_Position") or row.get("residue") or "").strip()
            if not gene or not pos_raw.isdigit():
                continue
            pos = int(pos_raw)
            entries.setdefault(gene, []).append(
                {
                    "range": [pos, pos],
                    "domain": f"cancerhotspots v2 single ({gene} {pos})",
                    "strength": "supporting",
                    "source": "PMID 29247016",  # cancerhotspots.org v2 reference
                }
            )
    except csv.Error as e:
        logger.warning("cancerhotspots TSV parse failed: %s — skipping layer", e)
        return {}
    return entries


def _load_overrides(yaml_path: Path) -> dict[str, Any]:
    """Load the clinical-advisor overrides YAML.

    Missing or empty file → returns a blank overrides dict. That is permitted
    during bootstrap but discouraged; the acceptance tests require the YAML
    to be present for v2.2 so the clinical lock layer is always engaged.
    """
    if not yaml_path.exists():
        logger.warning("Overrides YAML missing at %s — building without clinical lock", yaml_path)
        return {"version": "1.0", "additions": {}, "exclusions": {}}
    with yaml_path.open() as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Overrides YAML at {yaml_path} did not parse as a mapping")
    data.setdefault("additions", {})
    data.setdefault("exclusions", {})
    return data


# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------


def _entry_key(entry: dict[str, Any]) -> tuple[int, int]:
    rng = entry.get("range") or [0, 0]
    return (int(rng[0]), int(rng[1]))


def _merge_layers(
    cancerhotspots: dict[str, list[dict[str, Any]]],
    baseline: dict[str, list[dict[str, Any]]],
    overrides: dict[str, Any],
) -> dict[str, list[dict[str, Any]]]:
    """Merge the three layers into a single gene → entries dict.

    Precedence (higher wins on identical `range` key):
        overrides  >  ClinGen baseline  >  cancerhotspots

    Exclusions from the overrides YAML remove any overlapping entries
    introduced by the lower layers.
    """
    merged: dict[str, dict[tuple[int, int], dict[str, Any]]] = {}

    def _install(gene: str, entry: dict[str, Any]) -> None:
        gene_map = merged.setdefault(gene, {})
        gene_map[_entry_key(entry)] = {
            "range": list(entry["range"]),
            "domain": entry.get("domain", ""),
            "strength": entry.get("strength", "supporting"),
            "source": entry.get("source", ""),
            "rationale": entry.get("rationale", ""),
            "locked": bool(entry.get("locked", False)),
        }

    for gene, entries in cancerhotspots.items():
        for e in entries:
            _install(gene, e)

    for gene, entries in baseline.items():
        for e in entries:
            _install(gene, e)

    additions = overrides.get("additions") or {}
    for gene, entries in additions.items():
        for e in entries:
            _install(gene, e)

    exclusions = overrides.get("exclusions") or {}
    for gene, ex_entries in exclusions.items():
        for ex in ex_entries or []:
            key = _entry_key({"range": ex["range"]})
            if gene in merged and key in merged[gene]:
                del merged[gene][key]

    # Flatten and sort for a deterministic output (critical for source_hash).
    flat: dict[str, list[dict[str, Any]]] = {}
    for gene in sorted(merged.keys()):
        gene_entries = sorted(merged[gene].values(), key=lambda e: (e["range"][0], e["range"][1]))
        flat[gene] = gene_entries
    return flat


# ---------------------------------------------------------------------------
# Validation + hashing
# ---------------------------------------------------------------------------


def _validate(genes: dict[str, list[dict[str, Any]]]) -> None:
    if not genes:
        raise ValueError("PM1 hotspot merge produced an empty genes payload")
    for gene, entries in genes.items():
        if not entries:
            raise ValueError(f"Gene {gene} has no entries after merge")
        for entry in entries:
            rng = entry.get("range")
            if not isinstance(rng, list) or len(rng) != 2 or not all(isinstance(x, int) for x in rng):
                raise ValueError(f"{gene}: range must be [start, end] integers, got {rng!r}")
            if rng[0] > rng[1]:
                raise ValueError(f"{gene}: range {rng} has start > end")
            strength = entry.get("strength")
            if strength not in ALLOWED_STRENGTHS:
                raise ValueError(f"{gene} range {rng}: strength must be one of {ALLOWED_STRENGTHS}, got {strength!r}")
            source = entry.get("source", "")
            if not PMID_PATTERN.search(source):
                raise ValueError(f"{gene} range {rng}: source {source!r} does not contain a PMID reference")


def _assert_coverage(genes: dict[str, list[dict[str, Any]]]) -> None:
    """Fail the build if any required (gene, residue) pair is missing.

    Clinical-advisor 2026-04-14: this assertion protects against a future edit
    accidentally dropping a gene or residue and silently shipping an incomplete
    PM1 table.
    """
    missing: list[str] = []
    for gene, required_residues in REQUIRED_COVERAGE.items():
        entries = genes.get(gene, [])
        for residue in required_residues:
            covered = any(e["range"][0] <= residue <= e["range"][1] for e in entries)
            if not covered:
                missing.append(f"{gene}:{residue}")
    if missing:
        raise ValueError(
            "PM1 hotspot table is missing required coverage for "
            f"{len(missing)} residue(s): {', '.join(missing)}. "
            "Edit data/pm1_hotspot_overrides.yaml or CLINGEN_BASELINE to restore."
        )


def _compute_source_hash(genes: dict[str, list[dict[str, Any]]]) -> str:
    """Deterministic sha256 over the genes payload. Sort keys so semantically
    identical payloads produce identical hashes regardless of dict order."""
    canonical = json.dumps(genes, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    return f"sha256:{digest}"


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def build_pm1_hotspots(
    overrides_path: Path | str = DEFAULT_OVERRIDES_PATH,
    cancerhotspots_tsv: Path | str = DEFAULT_CANCERHOTSPOTS_TSV,
    output_path: Path | str = DEFAULT_OUTPUT_PATH,
) -> dict[str, Any]:
    """Build the PM1 hotspot JSON. Idempotent; safe to re-run."""
    overrides_path = Path(overrides_path)
    cancerhotspots_tsv = Path(cancerhotspots_tsv)
    output_path = Path(output_path)

    cancerhotspots = _load_cancerhotspots(cancerhotspots_tsv)
    overrides = _load_overrides(overrides_path)

    genes = _merge_layers(cancerhotspots, CLINGEN_BASELINE, overrides)
    _validate(genes)
    _assert_coverage(genes)

    record_count = sum(len(v) for v in genes.values())
    source_hash = _compute_source_hash(genes)

    source_refs = list(SOURCE_REFS_BASE)
    if cancerhotspots:
        source_refs.append(SOURCE_REF_CANCERHOTSPOTS)

    payload: dict[str, Any] = {
        "$schema": "https://biko/schema/pm1_hotspot_domains.v1.json",
        "version": overrides.get("version", "1.0"),
        "build_date": date.today().isoformat(),
        "source_refs": source_refs,
        "source_hash": source_hash,
        "record_count": record_count,
        "genes": genes,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
        f.write("\n")

    logger.info(
        "PM1 hotspot JSON built: %d entries across %d genes → %s (hash=%s)",
        record_count,
        len(genes),
        output_path,
        source_hash[:18],
    )
    return payload


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build_pm1_hotspots()
