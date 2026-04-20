"""Shared availability cache — log a DB's unavailability once, skip thereafter.

BIKO loads several optional local resources (ClinGen SQLite, gnomAD v4.1
constraint SQLite, DDG2P JSON, gnomAD population-frequency SQLite). When
one of those resources is missing or empty on disk, the first query
should emit a single WARNING and every subsequent query should return
``None`` / ``False`` silently — a 100-variant rare-disease run must not
produce 100 copies of the same "DB not found" warning.

This module consolidates the pattern first shipped in
``scripts/storage/query_local_clingen.py`` (v2.2 T3 hardening) and
``scripts/storage/query_gnomad_constraint.py``. New call-sites should use
either the :class:`AvailabilityCache` class (for DB probes with
multiple failure modes — missing vs empty-shell vs probe error — each
with a different warning message) or the :func:`check_availability`
module-level helper (for the simple missing-file-or-empty case).

## 3-tier graceful-degradation policy

BIKO's fallback chain splits every input into three tiers, each with
a different failure posture. The tier determines whether a missing /
malformed input warns, skips, or aborts the run:

1. **External API** (ClinVar API, OncoKB, gnomAD API, PharmGKB).
   Timeout + retry + fallback to the local DB. Handled in
   ``scripts/common/api_utils.py``; this module is not involved.

2. **Local DB** (ClinGen, ClinVar, gnomAD SQLite, DDG2P JSON,
   gnomAD-constraint SQLite). Log-once WARNING via this module +
   return empty / ``None`` / ``False``. The pipeline continues; the
   corresponding evidence code (PM1 hotspot, PP2 constraint, de novo
   neurodev carve-out, …) silently disables itself for the run.

3. **Required input** (primary VCF, config file). Fail-loud via
   :class:`scripts.common.exceptions.MalformedVCF` or
   :class:`scripts.common.exceptions.InvalidInput` — these inherit
   from ``ValueError`` so existing ``pytest.raises(ValueError)`` tests
   keep passing.

**Optional input** (``--ped``, ``--germline``) sits between tiers 2
and 3: warn + skip when the flag is absent; fail-loud when the user
explicitly provides the flag but it cannot be resolved (strict mode,
see ``scripts/intake/parse_vcf.py``). This asymmetry is deliberate —
a silently-misresolved trio would propagate into de novo
classification and mislead the reviewing researcher.
"""

from __future__ import annotations

import logging
import sqlite3
import threading
from pathlib import Path
from typing import Callable, Optional

logger = logging.getLogger(__name__)


class AvailabilityCache:
    """Per-namespace log-once cache for DB availability decisions.

    Each instance memoises a boolean availability decision per-path.
    Thread-safe. The probe callable is responsible for emitting the
    one-shot warning on failure — this class only guarantees the
    probe runs at most once per path per :meth:`reset` cycle.

    Typical usage in a DB-query module::

        _AVAILABILITY = AvailabilityCache()

        def reset_availability_cache() -> None:
            _AVAILABILITY.reset()

        def _probe(path: str) -> bool:
            if not Path(path).exists():
                logger.warning("ClinGen local DB not found at %s", path)
                return False
            return True

        def get_something(gene: str, path: str) -> Optional[str]:
            if not _AVAILABILITY.check(path, _probe):
                return None
            # ... real query ...
    """

    def __init__(self) -> None:
        self._cache: dict[str, bool] = {}
        self._lock = threading.Lock()

    def check(self, path: str, probe: Callable[[str], bool]) -> bool:
        """Return the cached availability for ``path``, running ``probe``
        at most once per path per reset cycle.

        The probe receives the path string and must return ``True`` iff
        the resource is usable. Any one-shot WARNING logging is the
        probe's responsibility; this class only guarantees the probe
        is invoked at most once.
        """
        cached = self._cache.get(path)
        if cached is not None:
            return cached
        with self._lock:
            cached = self._cache.get(path)
            if cached is not None:
                return cached
            ok = probe(path)
            self._cache[path] = ok
            return ok

    def reset(self) -> None:
        """Clear the memoised decisions — exposed so tests can rotate
        through multiple DB states in a single process."""
        with self._lock:
            self._cache.clear()


# Module-level shared cache for :func:`check_availability`. Call-sites
# that need their own isolated reset surface should instantiate
# :class:`AvailabilityCache` directly instead.
_SHARED_CACHE = AvailabilityCache()


def check_availability(
    label: str,
    path: str | Path,
    probe: Optional[Callable[[str], bool]] = None,
    message: Optional[str] = None,
) -> bool:
    """Log a WARNING once if the resource at ``path`` is missing / empty,
    then short-circuit every subsequent call for the same ``(label, path)``.

    ``label``
        Human-readable identifier (e.g. ``"ClinGen"``) used both to
        namespace the cache key and to build the default warning
        message.
    ``path``
        Filesystem path to probe. Stringified into the cache key so
        the same label+path across runs collapses to a single entry.
    ``probe``
        Optional custom probe ``(str) -> bool``. When ``None`` the
        default probe is :func:`_default_probe`: file exists, and for
        SQLite files has at least one table.
    ``message``
        Optional WARNING text. When ``None`` the default is
        ``f"{label} local DB not found at {path}"``.

    Returns ``True`` iff the resource is available. On first
    unavailability the single WARNING is emitted; subsequent calls
    silently return ``False``.
    """
    key = f"{label}:{path}"
    path_str = str(path)

    def _probe(_cache_key: str) -> bool:
        # ``_cache_key`` is the namespaced cache key used by
        # AvailabilityCache; the user-facing probe must receive the raw
        # filesystem path, so we close over ``path_str`` explicitly.
        ok = _default_probe(path) if probe is None else probe(path_str)
        if not ok:
            logger.warning(message or f"{label} local DB not found at {path}")
        return ok

    return _SHARED_CACHE.check(key, _probe)


def reset_shared_cache() -> None:
    """Clear the module-level cache backing :func:`check_availability`.

    Call-sites with their own :class:`AvailabilityCache` instance
    should call ``.reset()`` on that instance directly instead.
    """
    _SHARED_CACHE.reset()


def _default_probe(path: str | Path) -> bool:
    """Default availability probe: file exists, and for SQLite files
    has at least one table in ``sqlite_master``."""
    p = Path(path)
    if not p.exists():
        return False
    if p.suffix in (".sqlite3", ".db"):
        try:
            conn = sqlite3.connect(f"file:{p}?mode=ro", uri=True)
            try:
                row = conn.execute("SELECT name FROM sqlite_master LIMIT 1").fetchone()
            finally:
                conn.close()
            return row is not None
        except sqlite3.Error:
            return False
    return p.stat().st_size > 0
