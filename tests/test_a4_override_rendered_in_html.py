"""A4 render-layer regression: `clinvar_override_reason` surfaces in HTML.

End-to-end assertion that the A4 ClinVar-conflict override reason populated
by `apply_hotspot_conflict_reconciliation` flows all the way through
`build_variant_records` → Jinja render → the cancer `templates/cancer/report.html`
variant-detail card as a `.override-notice` block.

Covers qa-engineer REC A4 audit-log propagation and the report-dev HIGH A4
file-map correction (override rendering lives in the Jinja template, NOT in
`scripts/clinical_board/render.py`). The render seam uses a dedicated
`templates/macros/override.html` macro imported from the cancer template.
"""

from __future__ import annotations

import re

from scripts.classification.acmg_engine import ClassificationResult
from scripts.common.models import Variant
from scripts.reporting.generate_pdf import generate_report_html
from scripts.orchestration.classify import build_variant_records

# ─── Helpers ────────────────────────────────────────────────────────────────

# The canonical override-reason string from
# _workspace/v22-phaseA/artifacts/00_clinical_review.md §4 (165 chars, under cap).
TP53_R249M_REASON = (
    "engine LP override: PM1 (TP53 DBD L3 loop, PMID 30224644) + "
    "PM5 (R249S ClinVar Pathogenic); ClinVar R249M shows conflicting submitters"
)


def _tp53_r249m_variant() -> Variant:
    """Minimal TP53 R249M Variant dataclass for a smoke end-to-end render."""
    return Variant(
        chrom="chr17",
        pos=7577534,
        ref="C",
        alt="A",
        gene="TP53",
        hgvsc="c.746G>T",
        hgvsp="p.Arg249Met",
        consequence="missense_variant",
        transcript="NM_000546.6",
        impact="MODERATE",
    )


def _tp53_r249m_classification(reason: str = TP53_R249M_REASON) -> ClassificationResult:
    """ClassificationResult mimicking a post-A4-override TP53 R249M."""
    return ClassificationResult(
        classification="Likely Pathogenic",
        evidence_codes=["PM1", "PM5", "PM2_Supporting", "PP3"],
        conflict=False,
        clinvar_override_reason=reason,
    )


def _db_stub() -> dict:
    return {
        "clinvar": {
            "clinvar_significance": "Conflicting classifications of pathogenicity",
            "review_status": "criteria provided, conflicting interpretations",
            "clinvar_id": "RCV000123456",
        },
        "gnomad": {"gnomad_all": None, "gnomad_eas": None, "api_available": False},
        "kova_freq": None,
        "kova_homozygote": None,
    }


def _freq_stub() -> dict:
    return {"korean_flag": ""}


def _build_records(reason: str = TP53_R249M_REASON) -> list[dict]:
    variant = _tp53_r249m_variant()
    classification = _tp53_r249m_classification(reason=reason)
    vid = variant.variant_id
    records = build_variant_records(
        variants=[variant],
        db_results={vid: _db_stub()},
        freq_results={vid: _freq_stub()},
        classification_results={vid: classification},
        mode="cancer",
        hpo_results={},
    )
    return records


def _render_cancer(records: list[dict]) -> str:
    report_data = {
        "sample_id": "A4-SMOKE-001",
        "date": "2026-04-14",
        "variants": records,
        "detailed_variants": records,
        "pgx_results": [],
        "summary": {"total": len(records), "pathogenic": 0, "vus": 0, "benign": 0},
        "db_versions": {},
    }
    return generate_report_html(report_data, mode="cancer")


def _override_notice_blocks(html: str) -> list[str]:
    """Extract every `.override-notice` element's raw inner HTML."""
    return re.findall(
        r'<div[^>]*class="[^"]*override-notice[^"]*"[^>]*>(.*?)</div>',
        html,
        flags=re.DOTALL,
    )


# ─── Tests: data plumbing ───────────────────────────────────────────────────


def test_build_variant_records_surfaces_clinvar_override_reason():
    """build_variant_records must copy ClassificationResult.clinvar_override_reason
    onto the variant dict so the Jinja template can render it from `v.clinvar_override_reason`."""
    records = _build_records()
    assert len(records) == 1
    assert records[0].get("clinvar_override_reason") == TP53_R249M_REASON


def test_build_variant_records_override_reason_empty_when_absent():
    """When the engine did not fire A4, the field must still exist as an empty
    string so the template `{% if v.clinvar_override_reason %}` guard is safe."""
    records = _build_records(reason="")
    assert records[0].get("clinvar_override_reason", None) == ""


# ─── Tests: rendered HTML ──────────────────────────────────────────────────


def test_override_notice_rendered_in_variant_card_when_reason_present():
    """The cancer report HTML must contain a `.override-notice` block with
    'engine LP override' text when `v.clinvar_override_reason` is populated."""
    html = _render_cancer(_build_records())

    blocks = _override_notice_blocks(html)
    assert blocks, (
        "Expected at least one .override-notice block in rendered cancer HTML. "
        "The A4 render seam lives in templates/cancer/report.html + templates/macros/override.html "
        "(per report-dev HIGH A4 file-map correction)."
    )
    assert any("engine LP override" in b for b in blocks), (
        f"None of the {len(blocks)} .override-notice blocks contained the "
        f"'engine LP override' substring. Blocks: {blocks!r}"
    )


def test_override_notice_absent_when_reason_empty():
    """Negative case: no `.override-notice` block when reason is blank."""
    html = _render_cancer(_build_records(reason=""))
    assert not _override_notice_blocks(html), "override-notice must not render when clinvar_override_reason is empty"


def test_override_notice_has_page_break_inside_avoid():
    """The `.override-notice` block must carry `page-break-inside: avoid`
    (directly on the element or via a CSS class rule in the document) so the
    reason text cannot straddle an A4 page break."""
    html = _render_cancer(_build_records())

    # Inline style on at least one override-notice element …
    inline_match = re.search(
        r'<div[^>]*class="[^"]*override-notice[^"]*"[^>]*style="[^"]*page-break-inside\s*:\s*avoid',
        html,
    )
    # … OR a CSS rule keyed on `.override-notice { ... page-break-inside: avoid; ... }`.
    css_match = re.search(
        r"\.override-notice\s*\{[^}]*page-break-inside\s*:\s*avoid",
        html,
    )
    assert inline_match or css_match, (
        "override-notice must set page-break-inside: avoid either via inline "
        "style or via a CSS rule in the cancer template."
    )


def test_override_notice_bilingual_ko_label():
    """When language='ko' the override notice must use the Korean label
    '엔진 재조정' (per plan §A4 bilingual spec)."""
    records = _build_records()
    report_data = {
        "sample_id": "A4-SMOKE-KO",
        "date": "2026-04-14",
        "variants": records,
        "detailed_variants": records,
        "pgx_results": [],
        "summary": {"total": 1, "pathogenic": 0, "vus": 0, "benign": 0},
        "db_versions": {},
        "language": "ko",
    }
    html = generate_report_html(report_data, mode="cancer")
    blocks = _override_notice_blocks(html)
    assert blocks, "Expected an override-notice block under language=ko"
    assert any("엔진 재조정" in b for b in blocks), (
        f"Korean label '엔진 재조정' missing from override-notice blocks: {blocks!r}"
    )


def test_override_notice_truncation_tail_when_over_200_chars():
    """If the reason somehow exceeds 200 chars at the render layer, the macro
    must display the first 200 chars + a 'full reason in methodology' tail
    (defence-in-depth; acmg_engine already caps at 200)."""
    long_reason = (
        "engine LP override: PM1 (TP53 DBD L3 loop, PMID 30224644) + "
        "PM5 (R249S ClinVar Pathogenic); ClinVar R249M shows conflicting submitters "
        "and further supporting narrative that intentionally pushes this string "
        "well beyond the two-hundred-character visual cap to exercise the tail."
    )
    assert len(long_reason) > 200
    records = _build_records(reason=long_reason)
    html = _render_cancer(records)
    blocks = _override_notice_blocks(html)
    assert blocks, "Expected an override-notice block even when reason is long"
    combined = " ".join(blocks)
    assert "full reason in methodology" in combined, "Long reason must render the English truncation tail"
    # The excess tail text must NOT appear verbatim in the notice (it was cut).
    assert "exercise the tail" not in combined, "Text beyond the 200-char cap must be truncated from the notice body"
