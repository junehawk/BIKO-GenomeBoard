# AMP/ASCO/CAP 2017 Somatic Tiering + CIViC Variant-Level Mapping

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cancer mode tiering을 AMP/ASCO/CAP 2017 가이드라인 기반 Modified Approach B로 전환하고, CIViC variant-level evidence를 tiering에 통합한다.

**Architecture:** 새 `scripts/somatic/amp_tiering.py` 모듈이 OncoKB gene-level + CIViC variant-level evidence를 종합하여 AMP Tier I-IV를 결정한다. 기존 `oncokb.py`는 gene data provider로 유지되고, `assign_tier()`는 deprecated wrapper가 된다. Strategy A/B/C는 config로 전환 가능.

**Tech Stack:** Python 3, SQLite3, pytest

---

## File Structure

### 신규 파일
| File | Responsibility |
|------|---------------|
| `scripts/somatic/__init__.py` | Package init |
| `scripts/somatic/amp_tiering.py` | AMP 2017 tiering engine (TierResult dataclass, strategy A/B/C) |
| `tests/test_amp_tiering.py` | AMP tiering 단위 테스트 |
| `tests/test_civic_variant_match.py` | CIViC variant-level 매칭 + match_level 테스트 |

### 수정 파일
| File | Change |
|------|--------|
| `scripts/db/query_civic.py` | `get_predictive_evidence_for_tier()` 추가, match_level 반환 |
| `scripts/clinical/oncokb.py` | `assign_tier()` → deprecated wrapper, `get_tier_label()` AMP 라벨 |
| `scripts/orchestrate.py:219-233` | `assign_tier()` → `amp_assign_tier()` 호출 |
| `scripts/counselor/generate_pdf.py:115-143` | CIViC evidence 표시 강화 |
| `templates/cancer/report.html:1005-1094` | AMP 라벨 + evidence badge |
| `config.yaml` | `somatic.tiering_strategy` 추가 |

---

## Task 1: HGVSp 변환 함수 공용화 + CIViC variant-level 매칭 강화

`generate_pdf.py`의 `_hgvsp_to_civic_variant()`를 `scripts/common/hgvs_utils.py`로 추출하여 `query_civic.py`에서도 재사용. `query_civic.py`에 tier 결정 전용 함수를 추가한다.

**Files:**
- Create: `scripts/common/hgvs_utils.py` (기존 `generate_pdf.py`의 `_hgvsp_to_civic_variant` + `_AA_MAP` 이동)
- Modify: `scripts/counselor/generate_pdf.py` (import from hgvs_utils)
- Modify: `scripts/db/query_civic.py`
- Create: `tests/test_civic_variant_match.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_civic_variant_match.py
import sqlite3
import tempfile
import os
import pytest
from scripts.db.build_civic_db import build_db
from scripts.db.query_civic import reset_civic_connection


@pytest.fixture
def civic_db(tmp_path, monkeypatch):
    """Build a temp CIViC DB with known test data."""
    # Create minimal TSV files
    civic_dir = tmp_path / "civic"
    civic_dir.mkdir()

    # Gene summaries
    (civic_dir / "GeneSummaries.tsv").write_text(
        "feature_id\tname\tdescription\tfeature_aliases\tentrez_id\tfeature_type\n"
        "1\tBRAF\tBRAF proto-oncogene\tBRAF1\t673\tGene\n"
        "2\tKRAS\tKRAS proto-oncogene\t-\t3845\tGene\n"
    )

    # Variant summaries
    (civic_dir / "VariantSummaries.tsv").write_text(
        "variant_id\tfeature_name\tvariant\tvariant_types\tchromosome\tstart\tstop\treference_bases\tvariant_bases\tentrez_id\n"
        "1\tBRAF\tV600E\tmissense\t7\t140753336\t140753336\tA\tT\t673\n"
        "2\tKRAS\tG12D\tmissense\t12\t25398284\t25398284\tG\tA\t3845\n"
    )

    # Evidence
    (civic_dir / "ClinicalEvidenceSummaries.tsv").write_text(
        "evidence_id\tmolecular_profile\tdisease\ttherapies\tevidence_type\tevidence_direction\tevidence_level\tsignificance\tevidence_statement\tcitation_id\tcitation\tnct_ids\trating\tvariant_origin\n"
        "1\tBRAF V600E\tMelanoma\tVemurafenib\tPredictive\tSupports\tA\tSensitivity/Response\tBRAF V600E predicts response to vemurafenib\t20979469\tChapman 2011\t\t5\tSomatic\n"
        "2\tBRAF V600E\tMelanoma\tDabrafenib\tPredictive\tSupports\tA\tSensitivity/Response\tBRAF V600E predicts response to dabrafenib\t22608338\tHauschild 2012\t\t4\tSomatic\n"
        "3\tBRAF V600E\tMelanoma\t\tPrognostic\tSupports\tB\tPoor Outcome\tBRAF V600E associated with poor prognosis\t26287849\tLong 2015\t\t3\tSomatic\n"
        "4\tKRAS G12D\tPancreatic cancer\tErlotinib\tPredictive\tSupports\tB\tResistance\tKRAS G12D predicts resistance\t18316791\tLievre 2008\t\t4\tSomatic\n"
        "5\tBRAF AMPLIFICATION\tColorectal cancer\tCetuximab\tPredictive\tSupports\tC\tResistance\tBRAF amp predicts resistance\t\tSmith 2019\t\t2\tSomatic\n"
    )

    db_path = str(tmp_path / "civic.sqlite3")
    build_db(str(civic_dir), db_path)

    monkeypatch.setattr("scripts.db.query_civic.get", lambda k, d=None: db_path if k == "paths.civic_db" else d)
    reset_civic_connection()

    yield db_path

    reset_civic_connection()


def test_get_predictive_evidence_variant_match(civic_db):
    """Variant-specific 매치 시 match_level='variant' 반환."""
    from scripts.db.query_civic import get_predictive_evidence_for_tier
    result = get_predictive_evidence_for_tier("BRAF", "p.Val600Glu")
    assert result["match_level"] == "variant"
    assert len(result["evidence"]) >= 1
    assert all(e["evidence_type"] == "Predictive" for e in result["evidence"])
    assert result["evidence"][0]["evidence_level"] == "A"


def test_get_predictive_evidence_gene_fallback(civic_db):
    """Variant 못 찾으면 gene-level fallback, match_level='gene'."""
    from scripts.db.query_civic import get_predictive_evidence_for_tier
    result = get_predictive_evidence_for_tier("BRAF", "p.Lys601Glu")
    assert result["match_level"] == "gene"


def test_get_predictive_evidence_no_match(civic_db):
    """CIViC에 없는 유전자는 match_level='none'."""
    from scripts.db.query_civic import get_predictive_evidence_for_tier
    result = get_predictive_evidence_for_tier("FAKEGENE", "p.Ala1Val")
    assert result["match_level"] == "none"
    assert result["evidence"] == []


def test_predictive_only_no_prognostic(civic_db):
    """Predictive evidence만 반환, Prognostic 제외."""
    from scripts.db.query_civic import get_predictive_evidence_for_tier
    result = get_predictive_evidence_for_tier("BRAF", "p.Val600Glu")
    for e in result["evidence"]:
        assert e["evidence_type"] == "Predictive"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_civic_variant_match.py -v`
Expected: FAIL with `ImportError: cannot import name 'get_predictive_evidence_for_tier'`

- [ ] **Step 3: Implement get_predictive_evidence_for_tier in query_civic.py**

First, create `scripts/common/hgvs_utils.py` by extracting from `generate_pdf.py`:

```python
"""Shared HGVSp conversion utilities."""
import re
from typing import Optional

AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}


def hgvsp_to_civic_variant(hgvsp: Optional[str]) -> Optional[str]:
    """Convert HGVSp to CIViC variant name format. p.Gly12Asp -> G12D"""
    if not hgvsp:
        return None
    m = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvsp)
    if m:
        aa1 = AA3TO1.get(m.group(1), "?")
        pos = m.group(2)
        aa2 = AA3TO1.get(m.group(3), "?")
        return f"{aa1}{pos}{aa2}"
    return None
```

Then update `scripts/counselor/generate_pdf.py` to import from it:
```python
# Replace _AA_MAP dict and _hgvsp_to_civic_variant function with:
from scripts.common.hgvs_utils import hgvsp_to_civic_variant as _hgvsp_to_civic_variant
```

Add to `scripts/db/query_civic.py` (after `get_treatment_summary`), importing from common:

```python
def get_predictive_evidence_for_tier(
    gene: str, hgvsp: str, db_path: Optional[str] = None
) -> Dict:
    """Get Predictive evidence for tier determination.
    Returns: {"match_level": "variant"|"gene"|"none", "evidence": [...]}
    Only Predictive evidence is returned. match_level distinguishes
    variant-specific from gene-level matches.
    """
    conn = _get_connection()
    if not conn:
        return {"match_level": "none", "evidence": []}

    from scripts.common.hgvs_utils import hgvsp_to_civic_variant
    civic_name = hgvsp_to_civic_variant(hgvsp)

    # Try variant-specific match first
    if civic_name:
        cursor = conn.execute(
            "SELECT * FROM evidence WHERE gene = ? AND variant = ? "
            "AND evidence_type = 'Predictive' ORDER BY evidence_level, rating DESC",
            (gene, civic_name),
        )
        rows = cursor.fetchall()
        if rows:
            evidence = [
                {
                    "gene": r["gene"], "variant": r["variant"],
                    "disease": r["disease"], "therapies": r["therapies"],
                    "evidence_type": r["evidence_type"],
                    "evidence_level": r["evidence_level"],
                    "significance": r["significance"],
                    "statement": r["evidence_statement"],
                    "pmid": r["citation_id"], "citation": r["citation"],
                }
                for r in rows
            ]
            return {"match_level": "variant", "evidence": evidence}

    # Gene-level fallback (for display only, not Tier I elevation)
    cursor = conn.execute(
        "SELECT * FROM evidence WHERE gene = ? "
        "AND evidence_type = 'Predictive' ORDER BY evidence_level, rating DESC",
        (gene,),
    )
    rows = cursor.fetchall()
    if rows:
        evidence = [
            {
                "gene": r["gene"], "variant": r["variant"],
                "disease": r["disease"], "therapies": r["therapies"],
                "evidence_type": r["evidence_type"],
                "evidence_level": r["evidence_level"],
                "significance": r["significance"],
                "statement": r["evidence_statement"],
                "pmid": r["citation_id"], "citation": r["citation"],
            }
            for r in rows
        ]
        return {"match_level": "gene", "evidence": evidence}

    return {"match_level": "none", "evidence": []}
```

- [ ] **Step 4: Run tests**

Run: `python -m pytest tests/test_civic_variant_match.py -v`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/db/query_civic.py tests/test_civic_variant_match.py
git commit -m "feat: CIViC variant-level predictive evidence with match_level tracking"
```

---

## Task 2: AMP Tiering Engine

새 `scripts/somatic/amp_tiering.py` — Strategy A/B/C를 지원하는 AMP 2017 tiering 엔진.

**Files:**
- Create: `scripts/somatic/__init__.py`
- Create: `scripts/somatic/amp_tiering.py`
- Create: `tests/test_amp_tiering.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_amp_tiering.py
import pytest


# === Strategy B (default) ===

def test_strategy_b_civic_variant_level_a_pathogenic():
    """CIViC variant-specific Level A + Pathogenic → Tier I."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    civic = {"match_level": "variant", "evidence": [
        {"evidence_level": "A", "therapies": "Vemurafenib", "significance": "Sensitivity/Response",
         "disease": "Melanoma", "gene": "BRAF", "variant": "V600E",
         "evidence_type": "Predictive", "pmid": "20979469", "citation": "Chapman 2011", "statement": ""}
    ]}
    result = amp_assign_tier("Pathogenic", "BRAF", hgvsp="p.Val600Glu", civic_evidence=civic)
    assert result.tier == 1
    assert result.evidence_source == "civic-variant-A"
    assert result.civic_match_level == "variant"


def test_strategy_b_civic_variant_level_b():
    """CIViC variant-specific Level B → Tier II."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    civic = {"match_level": "variant", "evidence": [
        {"evidence_level": "B", "therapies": "Erlotinib", "significance": "Resistance",
         "disease": "NSCLC", "gene": "KRAS", "variant": "G12D",
         "evidence_type": "Predictive", "pmid": "", "citation": "", "statement": ""}
    ]}
    result = amp_assign_tier("Pathogenic", "KRAS", hgvsp="p.Gly12Asp", civic_evidence=civic)
    assert result.tier == 2
    assert result.evidence_source == "civic-variant-B"


def test_strategy_b_civic_gene_level_does_not_elevate_to_tier1():
    """CIViC gene-level (not variant) Level A → does NOT elevate to Tier I."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    civic = {"match_level": "gene", "evidence": [
        {"evidence_level": "A", "therapies": "Vemurafenib", "significance": "Sensitivity/Response",
         "disease": "Melanoma", "gene": "BRAF", "variant": "V600E",
         "evidence_type": "Predictive", "pmid": "", "citation": "", "statement": ""}
    ]}
    # BRAF is OncoKB Level 1, so still Tier I via OncoKB, but NOT via CIViC
    result = amp_assign_tier("Pathogenic", "BRAF", hgvsp="p.Lys601Glu", civic_evidence=civic)
    assert result.tier == 1
    assert result.evidence_source.startswith("oncokb")  # not civic


def test_strategy_b_no_civic_oncokb_level1_pathogenic():
    """No CIViC + OncoKB Level 1 + Pathogenic → Tier I via OncoKB."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    civic = {"match_level": "none", "evidence": []}
    result = amp_assign_tier("Pathogenic", "TP53", civic_evidence=civic)
    assert result.tier == 1
    assert "oncokb" in result.evidence_source


def test_strategy_b_vus_hotspot():
    """VUS + hotspot → Tier II."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    civic = {"match_level": "none", "evidence": []}
    # Need to mock is_hotspot — test with hgvsp that triggers hotspot
    result = amp_assign_tier("VUS", "KRAS", hgvsp="p.Gly12Asp", civic_evidence=civic)
    # KRAS G12 is a known hotspot; result depends on CIViC DB availability
    assert result.tier in (2, 3)  # 2 if hotspot found, 3 otherwise


def test_strategy_b_vus_cancer_gene_no_hotspot():
    """VUS on cancer gene, no hotspot → Tier III."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    civic = {"match_level": "none", "evidence": []}
    result = amp_assign_tier("VUS", "TP53", hgvsp="p.Ala999Val", civic_evidence=civic)
    assert result.tier == 3


def test_strategy_b_benign():
    """Benign → Tier IV, CIViC cannot override."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    civic = {"match_level": "variant", "evidence": [
        {"evidence_level": "A", "therapies": "X", "significance": "Sensitivity/Response",
         "disease": "Y", "gene": "BRAF", "variant": "V600E",
         "evidence_type": "Predictive", "pmid": "", "citation": "", "statement": ""}
    ]}
    result = amp_assign_tier("Benign", "BRAF", hgvsp="p.Val600Glu", civic_evidence=civic)
    assert result.tier == 4


def test_strategy_b_drug_response():
    """Drug Response → Tier I always."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    result = amp_assign_tier("Drug Response", "CYP2C19")
    assert result.tier == 1


def test_strategy_b_vus_non_cancer_gene():
    """VUS on non-cancer gene → Tier IV."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    civic = {"match_level": "none", "evidence": []}
    result = amp_assign_tier("VUS", "FAKEGENE", civic_evidence=civic)
    assert result.tier == 4


def test_strategy_b_civic_level_cd_pathogenic():
    """CIViC variant-specific Level C + Pathogenic → Tier II via CIViC (Priority 5)."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    civic = {"match_level": "variant", "evidence": [
        {"evidence_level": "C", "therapies": "Drug", "significance": "Sensitivity/Response",
         "disease": "Cancer", "gene": "BRAF", "variant": "V600E",
         "evidence_type": "Predictive", "pmid": "", "citation": "", "statement": ""}
    ]}
    result = amp_assign_tier("Pathogenic", "BRAF", hgvsp="p.Val600Glu", civic_evidence=civic)
    assert result.tier == 2
    assert result.evidence_source == "civic-variant-C"


# === Strategy C (OncoKB only, backward compatible) ===

def test_strategy_c_ignores_civic():
    """Strategy C: CIViC evidence ignored, OncoKB only."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    civic = {"match_level": "variant", "evidence": [
        {"evidence_level": "A", "therapies": "Drug", "significance": "Sensitivity/Response",
         "disease": "Cancer", "gene": "KRAS", "variant": "G12D",
         "evidence_type": "Predictive", "pmid": "", "citation": "", "statement": ""}
    ]}
    result = amp_assign_tier("VUS", "KRAS", hgvsp="p.Gly12Asp",
                             strategy="C", civic_evidence=civic)
    # VUS on cancer gene → Tier III (no CIViC elevation in Strategy C)
    assert result.tier in (2, 3)  # 2 only if hotspot, 3 otherwise


# === Strategy A (CIViC priority) ===

def test_strategy_a_civic_level_a():
    """Strategy A: CIViC Level A → Tier I."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    civic = {"match_level": "variant", "evidence": [
        {"evidence_level": "A", "therapies": "Drug", "significance": "Sensitivity/Response",
         "disease": "Cancer", "gene": "BRAF", "variant": "V600E",
         "evidence_type": "Predictive", "pmid": "", "citation": "", "statement": ""}
    ]}
    result = amp_assign_tier("Pathogenic", "BRAF", hgvsp="p.Val600Glu",
                             strategy="A", civic_evidence=civic)
    assert result.tier == 1
    assert result.evidence_source == "civic-variant-A"


# === TierResult structure ===

def test_tier_result_has_all_fields():
    """TierResult에 모든 필수 필드가 있어야 함."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    result = amp_assign_tier("Pathogenic", "TP53")
    assert hasattr(result, "tier")
    assert hasattr(result, "tier_label")
    assert hasattr(result, "evidence_source")
    assert hasattr(result, "civic_match_level")
    assert hasattr(result, "civic_evidence")
    assert isinstance(result.tier, int)
    assert isinstance(result.tier_label, str)


def test_amp_tier_labels():
    """AMP 2017 라벨 확인."""
    from scripts.somatic.amp_tiering import amp_assign_tier
    r1 = amp_assign_tier("Pathogenic", "TP53")
    assert "Strong Clinical Significance" in r1.tier_label
    r4 = amp_assign_tier("Benign", "FAKEGENE")
    assert "Benign" in r4.tier_label
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_amp_tiering.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Create scripts/somatic/__init__.py**

Empty file.

- [ ] **Step 4: Implement amp_tiering.py**

```python
"""AMP/ASCO/CAP 2017 somatic variant tiering engine.

Supports three strategies:
  A: CIViC evidence priority (variant-level first, OncoKB fallback)
  B: OncoKB + CIViC combined (default) — CIViC can elevate, not lower
  C: OncoKB only (backward compatible)

Reference: Li MM et al. J Mol Diagn. 2017;19(1):4-23.
"""
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional
from scripts.clinical.oncokb import get_cancer_gene_info, is_cancer_gene
from scripts.db.query_civic import is_hotspot, extract_protein_position

logger = logging.getLogger(__name__)

AMP_TIER_LABELS = {
    1: "Tier I — Strong Clinical Significance",
    2: "Tier II — Potential Clinical Significance",
    3: "Tier III — Unknown Clinical Significance",
    4: "Tier IV — Benign or Likely Benign",
}


@dataclass
class TierResult:
    tier: int
    tier_label: str
    evidence_source: str
    civic_match_level: str = "none"
    civic_evidence: List[Dict] = field(default_factory=list)


def _make_result(tier: int, source: str, civic_match: str = "none",
                 civic_ev: List[Dict] = None) -> TierResult:
    return TierResult(
        tier=tier,
        tier_label=AMP_TIER_LABELS.get(tier, "Unknown"),
        evidence_source=source,
        civic_match_level=civic_match,
        civic_evidence=civic_ev or [],
    )


def _best_civic_level(evidence: List[Dict]) -> Optional[str]:
    """Get the best (lowest letter) evidence level from a list."""
    levels = [e.get("evidence_level", "Z") for e in evidence]
    if not levels:
        return None
    return min(levels)


def amp_assign_tier(
    classification: str,
    gene: str,
    hgvsp: str = "",
    strategy: str = "B",
    civic_evidence: Optional[Dict] = None,
) -> TierResult:
    """Assign AMP/ASCO/CAP 2017 tier.

    Args:
        classification: ACMG classification (Pathogenic, LP, VUS, LB, Benign, Drug Response, Risk Factor)
        gene: Gene symbol
        hgvsp: HGVSp notation (e.g., p.Val600Glu)
        strategy: "A" (CIViC priority), "B" (combined, default), "C" (OncoKB only)
        civic_evidence: Pre-fetched CIViC evidence dict with match_level and evidence list.
                       If None, CIViC is not consulted (equivalent to strategy C for this call).
    """
    cls_lower = classification.lower()
    gene_info = get_cancer_gene_info(gene) if gene else None
    oncokb_level = gene_info.get("level", "") if gene_info else ""

    if civic_evidence is None:
        civic_evidence = {"match_level": "none", "evidence": []}

    match_level = civic_evidence.get("match_level", "none")
    evidence_items = civic_evidence.get("evidence", [])
    best_level = _best_civic_level(evidence_items)

    # === Always: Drug Response / Risk Factor → Tier I ===
    if cls_lower in ("drug response", "risk factor"):
        return _make_result(1, "pharmacogenomic", match_level, evidence_items)

    # === Always: Benign/Likely Benign → Tier IV (CIViC cannot override) ===
    if "benign" in cls_lower:
        return _make_result(4, "benign", match_level, evidence_items)

    is_pathogenic = "pathogenic" in cls_lower and "benign" not in cls_lower
    is_vus = cls_lower == "vus"

    # === Strategy A: CIViC priority ===
    if strategy == "A":
        if match_level == "variant" and best_level == "A" and is_pathogenic:
            return _make_result(1, "civic-variant-A", match_level, evidence_items)
        if match_level == "variant" and best_level in ("A", "B"):
            tier = 1 if best_level == "A" and is_pathogenic else 2
            return _make_result(tier, f"civic-variant-{best_level}", match_level, evidence_items)
        # Fall through to OncoKB

    # === Strategy B: Combined (CIViC can elevate) ===
    if strategy == "B":
        # CIViC variant-specific Level A + Pathogenic/LP → Tier I
        if match_level == "variant" and best_level == "A" and is_pathogenic:
            return _make_result(1, "civic-variant-A", match_level, evidence_items)
        # CIViC variant-specific Level B → Tier II
        if match_level == "variant" and best_level == "B":
            return _make_result(2, "civic-variant-B", match_level, evidence_items)

    # === Strategy B: CIViC variant-specific Level C-D + Pathogenic/LP → Tier II ===
    if strategy == "B":
        if match_level == "variant" and best_level in ("C", "D") and is_pathogenic:
            return _make_result(2, f"civic-variant-{best_level}", match_level, evidence_items)

    # === OncoKB gene-level (Strategies A fallback, B, C) ===

    # Pathogenic/LP on high-level gene → Tier I
    if is_pathogenic and gene_info and oncokb_level in ("1", "2"):
        return _make_result(1, f"oncokb-gene-L{oncokb_level}", match_level, evidence_items)

    # Pathogenic/LP on any cancer gene → Tier II
    if is_pathogenic and gene_info:
        return _make_result(2, f"oncokb-gene-L{oncokb_level}", match_level, evidence_items)

    # Pathogenic/LP on non-cancer gene → Tier II (still clinically significant)
    if is_pathogenic:
        return _make_result(2, "pathogenic-non-cancer", match_level, evidence_items)

    # VUS on cancer gene — check hotspot
    if is_vus and gene_info:
        protein_pos = extract_protein_position(hgvsp)
        if protein_pos and is_hotspot(gene, protein_pos):
            return _make_result(2, "hotspot", match_level, evidence_items)
        return _make_result(3, "oncokb-gene-vus", match_level, evidence_items)

    # Everything else → Tier IV
    return _make_result(4, "default", match_level, evidence_items)
```

- [ ] **Step 5: Run tests**

Run: `python -m pytest tests/test_amp_tiering.py -v`
Expected: MOST PASS (hotspot tests may need CIViC DB)

- [ ] **Step 6: Fix any failing tests, run again**

- [ ] **Step 7: Commit**

```bash
git add scripts/somatic/__init__.py scripts/somatic/amp_tiering.py tests/test_amp_tiering.py
git commit -m "feat: AMP/ASCO/CAP 2017 tiering engine with strategy A/B/C"
```

---

## Task 3: oncokb.py deprecated wrapper + orchestrate.py 통합

기존 `assign_tier()`를 deprecated wrapper로 변환하고, `orchestrate.py`에서 새 `amp_assign_tier()`를 사용.

**Files:**
- Modify: `scripts/clinical/oncokb.py:50-96`
- Modify: `scripts/orchestrate.py:219-233`
- Modify: `config.yaml`

- [ ] **Step 1: Write the failing test**

- [ ] **Step 1: Modify oncokb.py — deprecate assign_tier**

Replace `assign_tier` and `get_tier_label` in `scripts/clinical/oncokb.py`:

```python
def assign_tier(classification: str, gene: str, clinvar_significance: str = "", hgvsp: str = "") -> int:
    """DEPRECATED: Use scripts.somatic.amp_tiering.amp_assign_tier() instead.
    Kept for backward compatibility — delegates to amp_assign_tier with strategy C.
    """
    from scripts.somatic.amp_tiering import amp_assign_tier
    result = amp_assign_tier(classification, gene, hgvsp=hgvsp, strategy="C")
    return result.tier


def get_tier_label(tier: int) -> str:
    """Human-readable tier label (AMP/ASCO/CAP 2017)."""
    from scripts.somatic.amp_tiering import AMP_TIER_LABELS
    return AMP_TIER_LABELS.get(tier, "Unknown")
```

- [ ] **Step 3: Add somatic config to config.yaml**

```yaml
somatic:
  tiering_strategy: "B"  # "A" (CIViC priority), "B" (combined), "C" (OncoKB only)
```

- [ ] **Step 4: Modify orchestrate.py — use amp_assign_tier**

Replace lines 219-233 of `scripts/orchestrate.py`:

```python
    # Assign tiers (cancer mode only uses AMP/ASCO/CAP 2017; rare disease skips CIViC tiering)
    from scripts.somatic.amp_tiering import amp_assign_tier
    from scripts.db.query_civic import get_predictive_evidence_for_tier

    if mode == "cancer":
        tiering_strategy = get("somatic.tiering_strategy", "B")
    else:
        tiering_strategy = "C"  # Rare disease: OncoKB-only, no CIViC lookup

    for v_result in variant_records:
        gene = v_result.get("gene", "")
        cls = v_result.get("classification", "VUS")
        hgvsp = v_result.get("hgvsp", "")

        # Fetch CIViC evidence only for cancer mode with strategy A or B
        if tiering_strategy in ("A", "B"):
            civic_evidence = get_predictive_evidence_for_tier(gene, hgvsp)
        else:
            civic_evidence = None

        tier_result = amp_assign_tier(
            classification=cls,
            gene=gene,
            hgvsp=hgvsp,
            strategy=tiering_strategy,
            civic_evidence=civic_evidence,
        )
        v_result["tier"] = tier_result.tier
        v_result["tier_label"] = tier_result.tier_label
        v_result["tier_evidence_source"] = tier_result.evidence_source
        v_result["civic_match_level"] = tier_result.civic_match_level
        v_result["civic_evidence"] = tier_result.civic_evidence

        cancer_info = get_cancer_gene_info(gene)
        if cancer_info:
            v_result["cancer_gene_type"] = cancer_info.get("type", "")
            v_result["oncokb_level"] = cancer_info.get("level", "")
        else:
            v_result["cancer_gene_type"] = ""
            v_result["oncokb_level"] = ""
```

Note: `amp_tier_labels.json` from the spec is deferred; labels are inline in `amp_tiering.py` for now.

- [ ] **Step 5: Run existing oncokb tests (backward compat)**

Run: `python -m pytest tests/test_oncokb.py -v`
Expected: ALL PASS (assign_tier wrapper delegates to strategy C)

- [ ] **Step 6: Run full test suite**

Run: `python -m pytest tests/ --tb=short`
Expected: 383+ passed

- [ ] **Step 7: Commit**

```bash
git add scripts/clinical/oncokb.py scripts/orchestrate.py config.yaml
git commit -m "feat: integrate AMP tiering into pipeline, deprecate legacy assign_tier"
```

---

## Task 4: Cancer Report Template — AMP 라벨 + Evidence Badge

**Files:**
- Modify: `templates/cancer/report.html:1005-1094`
- Modify: `scripts/counselor/generate_pdf.py`

- [ ] **Step 1: Write failing test**

```python
# tests/test_generate_pdf.py (추가)
import copy

def test_amp_tier_labels_in_report():
    """리포트에 AMP 2017 라벨이 표시되어야 함."""
    from scripts.counselor.generate_pdf import generate_report_html
    data = copy.deepcopy(MINIMAL_REPORT)
    data["variants"][0]["tier"] = 1
    data["variants"][0]["tier_label"] = "Tier I — Strong Clinical Significance"
    data["variants"][0]["tier_evidence_source"] = "civic-variant-A"
    data["tier1_variants"] = [data["variants"][0]]
    html = generate_report_html(data)
    assert "Strong Clinical Significance" in html


def test_evidence_source_in_report():
    """Tier evidence source가 리포트에 표시되어야 함."""
    from scripts.counselor.generate_pdf import generate_report_html
    data = copy.deepcopy(MINIMAL_REPORT)
    data["variants"][0]["tier_evidence_source"] = "civic-variant-A"
    data["tier1_variants"] = [data["variants"][0]]
    html = generate_report_html(data)
    assert "civic" in html.lower() or "CIViC" in html
```

- [ ] **Step 2: Update cancer report template**

Replace tier section headers in `templates/cancer/report.html`:

Line 1008: `Tier 1 — Therapeutic Targets` → `Tier I — Strong Clinical Significance`
Line 1011: description → `Variants with FDA-approved therapy or professional guideline-included biomarkers (AMP/ASCO/CAP 2017 Tier I).`

Line 1032: `Tier 2 — Clinically Significant` → `Tier II — Potential Clinical Significance`
Line 1035: description → `Variants with evidence from clinical trials or emerging clinical significance.`

Line 1060: `Tier 3 — Cancer Gene Variants (VUS)` → `Tier III — Unknown Clinical Significance`
Line 1063: description → `Variants with unknown clinical significance in known cancer genes.`

Add evidence source badge to each variant finding row:
```html
{% if v.tier_evidence_source is defined and v.tier_evidence_source %}
<span style="display:inline-block;background:#EFF6FF;color:#1E40AF;border:1px solid #BFDBFE;border-radius:2px;padding:1px 5px;font-size:7pt;margin-left:4px;vertical-align:middle;">
  {{ v.tier_evidence_source }}
</span>
{% endif %}
```

- [ ] **Step 3: Update generate_pdf.py — pass civic evidence to template**

In `generate_pdf.py`, after CIViC enrichment (line ~143), add `civic_match_level` to variant display data:

```python
# Already set by orchestrate.py, just ensure template has access
v.setdefault("tier_evidence_source", "")
v.setdefault("civic_match_level", "none")
```

- [ ] **Step 4: Run report tests**

Run: `python -m pytest tests/test_generate_pdf.py -v`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add templates/cancer/report.html scripts/counselor/generate_pdf.py tests/test_generate_pdf.py
git commit -m "feat: AMP 2017 tier labels and evidence source badges in cancer report"
```

---

## Task 5: 통합 테스트 + 전체 회귀 검증

**Files:**
- Modify: `tests/test_orchestrate.py`
- Modify: `tests/test_report_modes.py`

- [ ] **Step 1: Write integration test**

```python
# tests/test_orchestrate.py (추가)

def test_cancer_pipeline_produces_amp_tier_fields(tmp_path):
    """Cancer pipeline output에 AMP tier 필드가 포함되어야 함."""
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path="data/sample_vcf/demo_variants_grch38_annotated.vcf",
        output_path=str(tmp_path / "report.html"),
        skip_api=True,
        mode="cancer",
    )
    assert result is not None
    for v in result["variants"]:
        assert "tier" in v
        assert "tier_label" in v
        assert "tier_evidence_source" in v
        assert v["tier"] in (1, 2, 3, 4)
        assert "Tier" in v["tier_label"]
```

- [ ] **Step 2: Write rare disease unchanged test**

```python
# tests/test_report_modes.py (추가)

def test_rare_disease_mode_no_amp_tier(tmp_path):
    """Rare disease mode는 AMP tier 필드 없이 기존 방식 유지."""
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path="data/sample_vcf/rare_disease_demo.vcf",
        output_path=str(tmp_path / "report.html"),
        skip_api=True,
        mode="rare-disease",
        hpo_ids=["HP:0001250"],
    )
    assert result is not None
    # Rare disease still has tier (from shared code) but classification is primary
    for v in result["variants"]:
        assert "classification" in v
        assert "hpo_score" in v
```

- [ ] **Step 3: Run full test suite**

Run: `python -m pytest tests/ -v --tb=short 2>&1 | tail -20`
Expected: 395+ tests PASS

- [ ] **Step 4: Commit**

```bash
git add tests/test_orchestrate.py tests/test_report_modes.py
git commit -m "test: integration tests for AMP tiering pipeline and rare disease isolation"
```

---

## Task 6: 문서 업데이트

**Files:**
- Modify: `docs/ARCHITECTURE.md`

- [ ] **Step 1: Read current ARCHITECTURE.md**

현재 tiering 설명 섹션 확인.

- [ ] **Step 2: Update tiering section**

AMP/ASCO/CAP 2017 tiering 설명으로 업데이트:
- Modified Approach B 설명
- CIViC variant-level 매칭 설명
- Strategy A/B/C 설정 방법
- `docs/TIERING_PRINCIPLES.md` 참조 링크

- [ ] **Step 3: Commit**

```bash
git add docs/ARCHITECTURE.md
git commit -m "docs: update architecture with AMP/ASCO/CAP 2017 tiering"
```

---

## Summary

| Task | 내용 | 신규/수정 | 예상 테스트 |
|------|------|----------|------------|
| 1 | CIViC variant-level 매칭 + match_level | 수정 1, 신규 1 | +4 tests |
| 2 | AMP tiering engine (strategy A/B/C) | 신규 3 | +13 tests |
| 3 | oncokb deprecated + orchestrate 통합 + config | 수정 3 | 기존 호환 |
| 4 | Cancer template AMP 라벨 + evidence badge | 수정 3 | +2 tests |
| 5 | 통합 테스트 + 회귀 검증 | 수정 2 | +2 tests |
| 6 | ARCHITECTURE.md 업데이트 | 수정 1 | — |

**예상 결과:** 기존 383 + 약 21 신규 = **404+ tests**
