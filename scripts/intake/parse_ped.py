# scripts/intake/parse_ped.py
"""Plink-style PED file parser for trio proband identification.

This module exists separately from ``parse_vcf.py`` (single-responsibility):
PED-file parsing and trio resolution logic is decoupled from VCF tokenisation
so each can be unit-tested independently.

Trio handling in BIKO prior to v2.4 relied on a filename heuristic
(``_detect_trio_proband`` in ``parse_vcf.py``) that (a) required exactly 3
samples in the VCF and (b) matched the VCF filename against sample IDs.
v2.4 Quick Win C introduces explicit PED-based resolution so cohort /
quartet / N≥3 VCFs can declare the trio via a standard PLINK PED file,
with **strict** failure (no silent fallback) when a PED is provided but
fails to resolve — getting the wrong proband silently propagates into
de novo classification and misleads the reviewing researcher.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Optional

logger = logging.getLogger(__name__)


@dataclass
class PedEntry:
    """One row of a Plink-style PED file.

    Attributes
    ----------
    family_id:
        Pedigree / family identifier (col 1).
    individual_id:
        Sample identifier (col 2), must match a VCF ``#CHROM`` sample column.
    father_id:
        Paternal sample ID or ``"0"`` when missing/unknown (col 3).
    mother_id:
        Maternal sample ID or ``"0"`` when missing/unknown (col 4).
    sex:
        1=male, 2=female, 0=unknown (col 5).
    affected:
        1=unaffected, 2=affected, 0 or -9=unknown (col 6).
    """

    family_id: str
    individual_id: str
    father_id: str
    mother_id: str
    sex: int
    affected: int


def parse_ped(ped_path: str) -> dict[str, PedEntry]:
    """Parse a Plink-style PED file into a ``{individual_id: PedEntry}`` dict.

    Format (whitespace-delimited, 6 columns):
        FID  IID  PID  MID  Sex  Phenotype

    Conventions handled:
    - ``"0"`` for missing parental IDs (standard PLINK convention)
    - ``"-9"`` or ``"0"`` for unknown phenotype
    - Lines starting with ``#`` are treated as comments and skipped
    - Empty / blank lines are skipped
    - UTF-8 BOM on the first line is stripped transparently
    - Windows CRLF line endings are handled by Python's text-mode reader

    A malformed row (< 6 whitespace-delimited tokens) raises
    ``ValueError`` with the 1-based line number — swallowing a malformed
    PED would let a silent misparse propagate into proband selection.

    Parameters
    ----------
    ped_path:
        Path to the PED file. An empty file produces an empty dict; a
        missing file raises ``FileNotFoundError`` from ``open()``.

    Returns
    -------
    dict[str, PedEntry]
        Keyed by ``individual_id``. When multiple rows share the same
        IID the *last* one wins and a warning is logged — BIKO's PED
        support does not currently cover multi-sample families with
        duplicate IIDs.
    """
    entries: dict[str, PedEntry] = {}
    if not os.path.exists(ped_path):
        raise FileNotFoundError(f"PED file not found: {ped_path}")

    with open(ped_path, encoding="utf-8-sig") as f:
        # "utf-8-sig" strips a UTF-8 BOM from the first line automatically;
        # CRLF is collapsed to LF by text-mode reading.
        for lineno, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            # PED is canonically tab-delimited but many tools emit spaces.
            # Use a permissive whitespace split to match PLINK behaviour.
            fields = line.split()
            if len(fields) < 6:
                raise ValueError(
                    f"Malformed PED line {lineno} in {ped_path}: "
                    f"expected ≥6 whitespace-delimited columns, got {len(fields)} "
                    f"({line!r})"
                )
            try:
                sex = int(fields[4])
            except ValueError:
                raise ValueError(
                    f"Malformed PED line {lineno} in {ped_path}: Sex column must be an integer, got {fields[4]!r}"
                ) from None
            try:
                affected = int(fields[5])
            except ValueError:
                raise ValueError(
                    f"Malformed PED line {lineno} in {ped_path}: Phenotype column must be an integer, got {fields[5]!r}"
                ) from None

            entry = PedEntry(
                family_id=fields[0],
                individual_id=fields[1],
                father_id=fields[2],
                mother_id=fields[3],
                sex=sex,
                affected=affected,
            )
            if entry.individual_id in entries:
                logger.warning(
                    "Duplicate individual_id %s on line %d of %s — last-wins",
                    entry.individual_id,
                    lineno,
                    ped_path,
                )
            entries[entry.individual_id] = entry

    return entries


def resolve_trio(
    ped: dict[str, PedEntry],
    sample_ids: list[str],
) -> tuple[Optional[int], list[int]]:
    """Identify proband + parents from a PED given the VCF sample order.

    Resolution rules (applied in order):

    1. Candidates = PED individuals with ``affected == 2`` (affected) that
       also appear in ``sample_ids``.
    2. Among candidates, prefer those whose ``father_id`` **and**
       ``mother_id`` are:
       (a) not ``"0"`` (not missing), AND
       (b) present in ``sample_ids``.
       A candidate missing either parent in the VCF is disqualified
       because de novo analysis needs the parental GT columns.
    3. Exactly one qualifying candidate → return
       ``(sample_ids.index(candidate_iid),
          [sample_ids.index(father), sample_ids.index(mother)])``.
    4. Zero qualifying candidates → return ``(None, [])``. The caller
       decides whether to strict-fail or fall back to a heuristic; this
       function never does the fallback itself.
    5. Multiple qualifying candidates → pick the first one alphabetically
       by ``individual_id`` and emit a WARNING listing all candidates.
       This keeps the choice deterministic for reproducibility while
       making the ambiguity visible in the log stream.

    Supports N ≥ 3 sample VCFs (quartet, cohort): the returned indices
    point into the provided ``sample_ids`` list, so a quartet VCF whose
    PED declares a 3-person trio subset resolves correctly.

    Parameters
    ----------
    ped:
        Output of :func:`parse_ped`, keyed by individual_id.
    sample_ids:
        The VCF's ``#CHROM`` sample columns, in order.

    Returns
    -------
    tuple[Optional[int], list[int]]
        ``(proband_idx, parent_idxs)`` or ``(None, [])`` on no-resolve.
        ``parent_idxs`` is always length 2 when ``proband_idx`` is set.
    """
    if not ped or not sample_ids:
        return None, []

    sample_set = set(sample_ids)

    # Step 1: affected individuals that also appear in the VCF.
    affected_in_vcf = [entry for entry in ped.values() if entry.affected == 2 and entry.individual_id in sample_set]
    if not affected_in_vcf:
        return None, []

    # Step 2: restrict to candidates with both parents present (not "0")
    # AND both parental sample IDs present in the VCF sample list.
    def _parents_usable(entry: PedEntry) -> bool:
        if entry.father_id == "0" or entry.mother_id == "0":
            return False
        if entry.father_id not in sample_set:
            return False
        if entry.mother_id not in sample_set:
            return False
        return True

    qualifying = [e for e in affected_in_vcf if _parents_usable(e)]
    if not qualifying:
        return None, []

    # Step 5 (applied before 3/4): deterministic alpha-sort so repeated
    # runs pick the same proband on multi-affected pedigrees.
    qualifying.sort(key=lambda e: e.individual_id)

    if len(qualifying) > 1:
        logger.warning(
            "PED has %d affected individuals with usable parents: %s. "
            "Choosing %s (first alphabetical) as proband. To disambiguate, "
            "provide a PED containing only the intended trio.",
            len(qualifying),
            [e.individual_id for e in qualifying],
            qualifying[0].individual_id,
        )

    proband = qualifying[0]
    proband_idx = sample_ids.index(proband.individual_id)
    father_idx = sample_ids.index(proband.father_id)
    mother_idx = sample_ids.index(proband.mother_id)
    return proband_idx, [father_idx, mother_idx]
