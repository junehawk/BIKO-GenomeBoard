"""Domain-specific exception hierarchy for BIKO GenomeBoard.

Rationale
---------

Before this module, input-validation failures in the intake / classification /
curator layers raised bare ``ValueError`` / ``RuntimeError`` / ``Exception``.
That made two things hard:

1. Distinguishing a BIKO domain error (malformed VCF, missing required
   config key, external API exhausted) from an unrelated standard-library
   ``ValueError`` in a broad ``except`` block.
2. Mapping errors to the 3-tier graceful-degradation policy documented
   in :mod:`scripts.common` — tier 3 (required input) should fail-loud
   with a typed error the orchestrator can surface as a human-readable
   report.

This module introduces :class:`BikoError` as the common root so callers
can write::

    try:
        run_pipeline(...)
    except BikoError as exc:
        # BIKO-specific handling — user error, not a bug
        ...

Backwards compatibility
-----------------------

:class:`InvalidInput` deliberately multi-inherits from
``(BikoError, ValueError)`` so every existing ``pytest.raises(ValueError)``
site (in ``tests/test_parse_ped.py``, ``tests/test_cache.py``,
``tests/test_curated_treatments.py``, ...) keeps passing unchanged. Only
the internal raise sites in ``scripts/intake/parse_vcf.py`` and
``scripts/intake/parse_ped.py`` are promoted to the typed class — the
catch side still sees a ``ValueError``.

Do **not** promote every ``raise ValueError`` in the codebase
indiscriminately. Promote only:

- Intake-layer failures where the user supplied malformed input
  (VCF, PED, config).
- External API exhaustion after retries (``ExternalAPIError``).
- Required local DB verifiably missing when it is *not* allowed to
  degrade silently (``DatabaseUnavailable``).

Leave bare ``ValueError`` in defensive coercion paths (argparse
converters, cache-key validation, dict-key checks) where the domain
meaning is ``"bad argument"`` rather than ``"bad user input"``.
"""


class BikoError(Exception):
    """Common root for all BIKO GenomeBoard domain errors.

    Catch this to recover from any BIKO-specific failure while still
    letting standard-library or third-party exceptions propagate.
    """


class InvalidInput(BikoError, ValueError):
    """User-supplied input is malformed or semantically invalid.

    Multi-inherits from ``ValueError`` so legacy
    ``except ValueError`` / ``pytest.raises(ValueError)`` call-sites
    continue to match — the BIKO-specific class is strictly narrower.
    """


class MalformedVCF(InvalidInput):
    """The primary (or germline) VCF file cannot be parsed.

    Used when the file exists but violates the VCF spec in a way that
    prevents intake from producing a per-variant record (missing
    ``#CHROM`` header, bad FORMAT field count, non-integer genotype, …).
    """


class DatabaseUnavailable(BikoError):
    """A local SQLite/JSON DB required for the current run is missing.

    Reserve this for the rare case where graceful degradation is
    *not* appropriate (e.g. the user explicitly requested a feature
    that has no fallback). The tier-2 default — log-once WARNING via
    :mod:`scripts.common.availability_cache` — stays the norm.
    """


class ExternalAPIError(BikoError):
    """An outbound HTTP call failed after the configured retry budget.

    Raised by :mod:`scripts.common.api_utils` only when the call-site
    explicitly opted out of the local-DB fallback chain. The default
    posture is still "return None / empty and let the local DB take
    over"; this exception signals that the caller asked for a fatal
    outcome on API failure.
    """
