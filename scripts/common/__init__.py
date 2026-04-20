"""Shared infrastructure for the BIKO GenomeBoard pipeline.

Modules in this package provide cross-cutting services consumed by
``scripts/intake``, ``scripts/classification``, ``scripts/clinical_board``,
``scripts/storage``, and ``scripts/orchestrate.py``. They are intentionally
free of any ACMG / AMP / PGx domain logic — anything here is expected
to be used from at least two pipeline stages.

## 3-tier graceful-degradation policy

BIKO classifies every runtime input by failure posture:

1. **External API** (ClinVar REST, OncoKB, gnomAD JSON, PharmGKB).
   Policy: timeout + bounded retry + fallback to the local DB.
   Owner: :mod:`scripts.common.api_utils`. Never fatal — the
   pipeline always has a local-DB route for the same evidence.

2. **Local DB** (ClinGen SQLite, ClinVar SQLite, gnomAD population
   SQLite, gnomAD v4.1 constraint SQLite, DDG2P JSON panel).
   Policy: log-once WARNING on first unavailability, return empty /
   ``None`` / ``False`` silently thereafter. The corresponding
   evidence code (PM1 hotspot, PP2 constraint, de novo neurodev
   carve-out, …) disables itself for the run.
   Owner: :mod:`scripts.common.availability_cache`
   (``AvailabilityCache`` class + ``check_availability`` helper).

3. **Required input** (primary VCF, ``config.yaml``). Policy:
   fail-loud via :mod:`scripts.common.exceptions`
   (:class:`MalformedVCF`, :class:`InvalidInput`). Domain exceptions
   inherit from both :class:`BikoError` and :class:`ValueError`, so
   existing ``pytest.raises(ValueError)`` call-sites stay green
   after the exception-hierarchy refactor.

**Optional input** (``--ped``, ``--germline``) sits between tiers 2
and 3. The posture depends on whether the user supplied the flag:

- Flag absent → warn + skip (tier 2 behaviour).
- Flag supplied but unresolvable → fail-loud (tier 3 behaviour).

This asymmetry is deliberate — a silently-misresolved trio would
propagate into de novo classification and mislead the reviewing
researcher. See ``scripts/intake/parse_vcf.py`` strict-mode handling.

New call-sites that query an optional local DB **must** route their
availability check through :mod:`scripts.common.availability_cache`
(either the class or the module-level helper) rather than reinventing
the log-once boolean. Do not extend the fail-loud set (tier 3)
without an explicit ticket — widening fail-loud raises the bar for
the pipeline to even start.
"""
