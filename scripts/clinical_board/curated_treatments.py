"""Curate-then-narrate treatment resolver (AI Clinical Board v2.2 · A1).

Deterministic curator that merges OncoKB free-tier hits with local CIViC
evidence into PMID-bearing treatment rows keyed by ``{chrom}:{pos}:{ref}:{alt}``.
Every row fed into the downstream narrator carries an audit trail so the
LLM cannot invent drug/target pairings (the futibatinib failure mode).

Patient-safety guarantees enforced here:

* Variants without genomic coordinates raise ``ValueError`` immediately.
* OncoKB network blips degrade to CIViC-only (single warning per session).
* Merges prefer stable ``therapy_ids`` set intersection, falling back to
  normalised drug-name matching with a small alias table.
* Every ``CuratedTreatment`` gets a stable 12-char ``curated_id`` the
  narrator cites back.
"""
from __future__ import annotations

import hashlib
import logging
import os
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from scripts.clinical import oncokb_client
from scripts.common.config import get
from scripts.common.hgvs_utils import normalize_hgvsp_for_civic

logger = logging.getLogger(__name__)

# Drug brand/generic alias pairs used by the drug-name fallback matcher.
# Kept intentionally short — extend in _workspace/v22-phaseA if new misses
# show up in the unmerged_near_matches.log during the pilot.
_DRUG_ALIASES: Dict[str, str] = {
    "zelboraf": "vemurafenib",
    "vemurafenib": "vemurafenib",
    "5-fu": "fluorouracil",
    "5-fluorouracil": "fluorouracil",
    "fluorouracil": "fluorouracil",
    "herceptin": "trastuzumab",
    "trastuzumab": "trastuzumab",
    "gleevec": "imatinib",
    "imatinib": "imatinib",
    "lumakras": "sotorasib",
    "sotorasib": "sotorasib",
    "tagrisso": "osimertinib",
    "osimertinib": "osimertinib",
}

_UNMERGED_LOG_PATH = "_workspace/v22-phaseA/unmerged_near_matches.log"

_LEVEL_RANK: Dict[str, int] = {"A": 0, "B": 1, "C": 2, "D": 3}

_DEGRADATION_WARNED = False


@dataclass
class CuratedTreatment:
    """A single deterministic treatment row the narrator may cite.

    ``curated_id`` is a stable short hash so the Board Chair's JSON output
    can reference the row by ID without leaking internal shapes.
    """

    curated_id: str
    variant_key: str
    drug: str
    target: str
    evidence_level: str            # Normalised AMP: A/B/C/D
    source: str                    # "oncokb" | "civic" | "both"
    pmids: List[str] = field(default_factory=list)
    disease_context: str = ""
    significance: str = "sensitivity"  # sensitivity | resistance | adverse
    therapy_ids: str = ""
    raw_row: Dict[str, Any] = field(default_factory=dict)


# ── Helpers ──────────────────────────────────────────────────────────────


def _variant_key(v: Dict[str, Any]) -> str:
    missing = [k for k in ("chrom", "pos", "ref", "alt") if k not in v or v[k] in (None, "")]
    if missing:
        raise ValueError(
            f"curate_treatments requires genomic coordinates; missing {missing} on variant {v.get('gene')} {v.get('hgvsp')}"
        )
    return f"{v['chrom']}:{v['pos']}:{v['ref']}:{v['alt']}"


def _canon_drug(name: Optional[str]) -> str:
    if not name:
        return ""
    key = str(name).strip().lower()
    # Strip common suffixes
    for suf in (" hcl", " hydrochloride", " sulfate", " mesylate"):
        if key.endswith(suf):
            key = key[: -len(suf)]
    return _DRUG_ALIASES.get(key, key)


def _curated_id(variant_key: str, drug: str, source: str) -> str:
    blob = f"{variant_key}|{_canon_drug(drug)}|{source}".encode()
    return hashlib.sha1(blob).hexdigest()[:12]


def _level_sort_key(row: "CuratedTreatment") -> Tuple[int, str]:
    return (_LEVEL_RANK.get(row.evidence_level, 9), row.drug.lower())


def _warn_degraded_once(gene: str, reason: str) -> None:
    global _DEGRADATION_WARNED
    if not _DEGRADATION_WARNED:
        logger.warning(
            "[curated_treatments] OncoKB unavailable (%s); degrading to CIViC-only for this session",
            reason,
        )
        _DEGRADATION_WARNED = True
    logger.debug("[curated_treatments] OncoKB miss for %s: %s", gene, reason)


def _log_unmerged(variant_key: str, oncokb_drug: str, civic_drug: str) -> None:
    try:
        path = Path(_UNMERGED_LOG_PATH)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", encoding="utf-8") as fh:
            fh.write(f"{variant_key}\t{oncokb_drug}\t{civic_drug}\n")
    except OSError:
        logger.debug("[curated_treatments] unmerged near-match log not writable")


# ── CIViC query (wrapped so tests can monkeypatch a single symbol) ────────


def _query_civic_for_variant(gene: str, hgvsp: str, *, offline_mode: bool = False) -> List[Dict[str, Any]]:
    """Return CIViC predictive-evidence rows for ``(gene, hgvsp)``.

    Kept as a single module-level function (not a method) so test suites
    can monkeypatch it with a lambda without rewiring the curator's internal
    control flow. Falls through to gene-level when the HGVSp → CIViC-name
    normaliser can't resolve a residue.
    """
    db_path = get("paths.civic_db", "data/db/civic.sqlite3")
    if not os.path.exists(db_path):
        return []

    variant_name = normalize_hgvsp_for_civic(hgvsp)
    try:
        conn = sqlite3.connect(db_path, check_same_thread=False)
        conn.row_factory = sqlite3.Row
        has_therapy_ids = any(
            row[1] == "therapy_ids"
            for row in conn.execute("PRAGMA table_info(evidence)").fetchall()
        )
        cols_common = (
            "gene, variant, disease, therapies, evidence_type, evidence_level, "
            "significance, evidence_statement, citation_id"
        )
        cols = cols_common + (", therapy_ids" if has_therapy_ids else "")
        rows: List[sqlite3.Row] = []
        if variant_name:
            rows = list(conn.execute(
                f"SELECT {cols} FROM evidence WHERE gene = ? AND variant = ? "
                "AND evidence_type = 'Predictive' ORDER BY evidence_level",
                (gene, variant_name),
            ))
        if not rows:
            rows = list(conn.execute(
                f"SELECT {cols} FROM evidence WHERE gene = ? AND evidence_type = 'Predictive' "
                "ORDER BY evidence_level",
                (gene,),
            ))
        conn.close()
    except sqlite3.Error as e:
        logger.debug("[curated_treatments] CIViC query error for %s: %s", gene, e)
        return []

    out: List[Dict[str, Any]] = []
    for r in rows:
        therapy_ids = r["therapy_ids"] if has_therapy_ids else ""
        pmids_str = r["citation_id"] or ""
        pmids = [p.strip() for p in pmids_str.split(",") if p.strip()]
        out.append({
            "drug": r["therapies"] or "",
            "level": (r["evidence_level"] or "D").upper(),
            "pmids": pmids,
            "disease": r["disease"] or "",
            "significance": (r["significance"] or "sensitivity").lower(),
            "therapy_ids": therapy_ids or "",
            "raw_row": dict(r),
        })
    return out


# ── Merge ────────────────────────────────────────────────────────────────


def _therapy_id_set(raw: str) -> set:
    if not raw:
        return set()
    return {t.strip() for t in str(raw).split(",") if t.strip()}


def _merge(
    variant_key: str,
    oncokb_rows: List[Dict[str, Any]],
    civic_rows: List[Dict[str, Any]],
) -> List[CuratedTreatment]:
    merged: List[CuratedTreatment] = []
    civic_consumed: set = set()

    for ok in oncokb_rows:
        ok_drug = ok.get("drug") or ""
        ok_canon = _canon_drug(ok_drug)
        ok_therapy_ids = _therapy_id_set(ok.get("therapy_ids", ""))

        match_idx: Optional[int] = None
        for idx, civ in enumerate(civic_rows):
            if idx in civic_consumed:
                continue
            civ_therapy_ids = _therapy_id_set(civ.get("therapy_ids", ""))
            if ok_therapy_ids and civ_therapy_ids and (ok_therapy_ids & civ_therapy_ids):
                match_idx = idx
                break
            if _canon_drug(civ.get("drug")) == ok_canon and ok_canon:
                match_idx = idx
                break

        if match_idx is not None:
            civ = civic_rows[match_idx]
            civic_consumed.add(match_idx)
            pmids = sorted(set(ok.get("pmids", []) + civ.get("pmids", [])))
            drug_label = ok_drug or civ.get("drug") or ""
            level = min(
                (ok.get("level") or "D", civ.get("level") or "D"),
                key=lambda lv: _LEVEL_RANK.get(lv, 9),
            )
            merged.append(CuratedTreatment(
                curated_id=_curated_id(variant_key, drug_label, "both"),
                variant_key=variant_key,
                drug=drug_label,
                target="",
                evidence_level=level,
                source="both",
                pmids=pmids,
                disease_context=civ.get("disease") or ok.get("disease") or "",
                significance=ok.get("significance") or civ.get("significance") or "sensitivity",
                therapy_ids=",".join(sorted(ok_therapy_ids | _therapy_id_set(civ.get("therapy_ids", "")))),
                raw_row={"oncokb": ok.get("raw_row", {}), "civic": civ.get("raw_row", {})},
            ))
        else:
            merged.append(CuratedTreatment(
                curated_id=_curated_id(variant_key, ok_drug, "oncokb"),
                variant_key=variant_key,
                drug=ok_drug,
                target="",
                evidence_level=ok.get("level") or "D",
                source="oncokb",
                pmids=list(ok.get("pmids", [])),
                disease_context=ok.get("disease", ""),
                significance=ok.get("significance") or "sensitivity",
                therapy_ids=",".join(sorted(ok_therapy_ids)),
                raw_row=ok.get("raw_row", {}),
            ))
            # Near-match logging for fuzzy drug-name misses
            for civ in civic_rows:
                civ_canon = _canon_drug(civ.get("drug"))
                if civ_canon and ok_canon and civ_canon != ok_canon and (
                    civ_canon.startswith(ok_canon[:4]) or ok_canon.startswith(civ_canon[:4])
                ):
                    _log_unmerged(variant_key, ok_drug, civ.get("drug") or "")
                    break

    for idx, civ in enumerate(civic_rows):
        if idx in civic_consumed:
            continue
        civ_drug = civ.get("drug") or ""
        merged.append(CuratedTreatment(
            curated_id=_curated_id(variant_key, civ_drug, "civic"),
            variant_key=variant_key,
            drug=civ_drug,
            target="",
            evidence_level=(civ.get("level") or "D").upper()[:1] if civ.get("level") else "D",
            source="civic",
            pmids=list(civ.get("pmids", [])),
            disease_context=civ.get("disease", ""),
            significance=civ.get("significance") or "sensitivity",
            therapy_ids=civ.get("therapy_ids", "") or "",
            raw_row=civ.get("raw_row", {}),
        ))

    # Normalise civic-only level letters (CIViC uses A..E; clamp unknown → D)
    for row in merged:
        if row.evidence_level not in _LEVEL_RANK:
            first = (row.evidence_level or "D")[:1].upper()
            row.evidence_level = first if first in _LEVEL_RANK else "D"

    merged.sort(key=_level_sort_key)
    return merged


# ── Public entrypoint ────────────────────────────────────────────────────


def curate_treatments(
    variants: List[Dict[str, Any]],
    offline_mode: Optional[bool] = None,
) -> Dict[str, List[CuratedTreatment]]:
    """Curate treatment rows for a pre-selected variant list.

    Args:
        variants: dicts with ``gene``, ``hgvsp``, ``chrom``, ``pos``,
            ``ref``, ``alt``. Genomic coords are MANDATORY (db-dev REQ A1-db-5).
        offline_mode: force-skip OncoKB. When ``None`` the config key
            ``clinical_board.curated_treatments.offline_mode`` is consulted.

    Returns:
        Mapping ``variant_key -> list[CuratedTreatment]``. Empty list when
        neither OncoKB nor CIViC have curated evidence for the variant.

    Raises:
        ValueError: when any variant lacks genomic coordinates.
    """
    global _DEGRADATION_WARNED
    _DEGRADATION_WARNED = False  # reset per call so tests see a fresh warning

    if offline_mode is None:
        offline_mode = bool(get("clinical_board.curated_treatments.offline_mode", False))

    out: Dict[str, List[CuratedTreatment]] = {}
    for v in variants:
        key = _variant_key(v)
        gene = v.get("gene") or ""
        hgvsp = v.get("hgvsp") or ""

        # OncoKB path (guarded by offline_mode)
        oncokb_rows: List[Dict[str, Any]] = []
        if offline_mode:
            oncokb_rows = oncokb_client.annotate_protein_change(gene, hgvsp, offline_mode=True)
        else:
            civic_name = normalize_hgvsp_for_civic(hgvsp) or hgvsp
            try:
                oncokb_rows = oncokb_client.annotate_protein_change(gene, civic_name, offline_mode=False)
            except oncokb_client.OncoKBUnavailable as e:
                _warn_degraded_once(gene, str(e))
                oncokb_rows = []

        # CIViC path — always local; safe in offline mode
        civic_rows = _query_civic_for_variant(gene, hgvsp, offline_mode=offline_mode)

        out[key] = _merge(key, oncokb_rows, civic_rows)
    return out


__all__ = ["CuratedTreatment", "curate_treatments"]
