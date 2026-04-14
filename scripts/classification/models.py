"""Public re-export of classification dataclasses.

``ClassificationResult`` lives in ``scripts.classification.acmg_engine`` for
historical reasons. This module re-exports it under the more natural
``scripts.classification.models`` path so downstream smoke tests and future
consumers have a stable import point.
"""

from scripts.classification.acmg_engine import ClassificationResult

__all__ = ["ClassificationResult"]
