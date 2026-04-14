"""Lightweight config.yaml schema validation for BIKO GenomeBoard.

No external dependencies — validates types and required keys at load time.
"""

import logging
from typing import Any

logger = logging.getLogger(__name__)

# Schema: {dotted_key: (expected_type, required)}
# required=True means load_config will raise ValueError if missing.
SCHEMA: dict[str, tuple[type, bool]] = {
    # paths — templates is critical (used by report generation)
    "paths.krgdb": (str, False),
    "paths.clinvar_db": (str, False),
    "paths.templates": (str, True),  # critical — report generation fails without it
    # thresholds — must be numeric
    "thresholds.ba1": (float, False),
    "thresholds.bs1": (float, False),
    "thresholds.pm2": (float, False),
    # api
    "api.timeout": (int, False),
    # cache
    "cache.enabled": (bool, False),
}

# Known top-level sections (for unknown-key warnings)
KNOWN_TOP_LEVEL_KEYS = {
    "paths",
    "thresholds",
    "api",
    "annotation",
    "report",
    "pgx",
    "somatic",
    "in_silico",
    "cache",
    "clinical_board",
    "_project_root",
}


def _resolve(config: dict, dotted_key: str) -> tuple[bool, Any]:
    """Walk *config* using a dotted key. Returns (found, value)."""
    parts = dotted_key.split(".")
    current = config
    for part in parts:
        if isinstance(current, dict) and part in current:
            current = current[part]
        else:
            return False, None
    return True, current


def validate_config(config: dict) -> list[str]:
    """Validate *config* against the schema.

    Returns a list of human-readable diagnostic strings.  Each string is
    prefixed with ``ERROR:`` (missing critical key) or ``WARNING:``.

    The caller decides how to handle them — ``load_config`` logs warnings
    and raises ``ValueError`` for errors.
    """
    messages: list[str] = []

    # 1. Check required keys and types
    for dotted_key, (expected_type, required) in SCHEMA.items():
        found, value = _resolve(config, dotted_key)

        if not found:
            if required:
                messages.append(f"ERROR: missing required key '{dotted_key}'")
            else:
                messages.append(f"WARNING: missing recommended key '{dotted_key}'")
            continue

        # Type check — int is acceptable where float is expected
        if expected_type is float:
            if not isinstance(value, (int, float)):
                messages.append(
                    f"WARNING: '{dotted_key}' should be {expected_type.__name__}, got {type(value).__name__}"
                )
        elif not isinstance(value, expected_type):
            messages.append(f"WARNING: '{dotted_key}' should be {expected_type.__name__}, got {type(value).__name__}")

    # 2. Warn on unknown top-level keys
    for key in config:
        if key not in KNOWN_TOP_LEVEL_KEYS:
            messages.append(f"WARNING: unknown top-level key '{key}'")

    return messages
