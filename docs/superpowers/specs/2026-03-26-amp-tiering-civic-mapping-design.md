# AMP/ASCO/CAP 2017 Somatic Tiering + CIViC Variant-Level Mapping

## Goal

Cancer mode의 variant tiering을 OncoKB gene-level 단독에서 AMP/ASCO/CAP 2017 가이드라인 기반으로 전환하고, CIViC variant-level treatment evidence를 tiering에 통합한다. Rare disease mode는 영향 없이 기존 ACMG/AMP 2015 분류를 유지한다.

## Scope

**In scope:**
- AMP/ASCO/CAP 2017 somatic tiering engine (strategy A/B/C 전환 가능)
- CIViC variant-level evidence 매칭 강화 (variant vs gene-level 구분)
- `tier_evidence_source` 필드 추가 (근거 추적성)
- Cancer 리포트 템플릿: AMP 2017 용어 + evidence level 표시
- config.yaml에 `tiering_strategy` 설정

**Out of scope:**
- Rare disease mode 변경 (기존 유지)
- TMB/MSI (향후)
- OncoKB API 실시간 조회 (로컬 JSON 유지)

---

## Architecture

### 현재 데이터 흐름

```
VCF → ACMG Engine → assign_tier(oncokb.py) → Tier 1-4
                          ↑
                     OncoKB gene-level only
```

### 새 데이터 흐름

```
VCF → ACMG Engine → amp_assign_tier(amp_tiering.py) → Tier I-IV
                          ↑              ↑
                   OncoKB gene-level   CIViC variant-level
                                       (match_level 구분)
```

### 신규/수정 파일

| File | Action | Responsibility |
|------|--------|---------------|
| `scripts/somatic/__init__.py` | Create | Package init |
| `scripts/somatic/amp_tiering.py` | Create | AMP 2017 tiering engine, strategy A/B/C |
| `scripts/db/query_civic.py` | Modify | variant-level 매칭 강화, match_level 반환 |
| `scripts/clinical/oncokb.py` | Modify | assign_tier() deprecated → amp_tiering 호출 |
| `scripts/orchestrate.py` | Modify | CIViC evidence를 tiering에 전달 |
| `scripts/counselor/generate_pdf.py` | Modify | tier_evidence_source, civic evidence 표시 |
| `templates/cancer/report.html` | Modify | AMP 2017 라벨, evidence badges |
| `config.yaml` | Modify | tiering_strategy 설정 추가 |
| `data/amp_tier_labels.json` | Create | Tier I-IV 라벨/설명 (i18n 대비) |
| `tests/test_amp_tiering.py` | Create | AMP tiering 단위 테스트 |
| `tests/test_civic_variant_match.py` | Create | Variant-level 매칭 테스트 |

---

## Component Design

### 1. `scripts/somatic/amp_tiering.py`

핵심 모듈. Strategy pattern으로 A/B/C 전환 가능.

```python
@dataclass
class TierResult:
    tier: int                    # 1-4 (I-IV)
    tier_label: str              # "Tier I — Strong Clinical Significance"
    evidence_source: str         # "civic-variant-A", "oncokb-gene-L1", "hotspot", etc.
    civic_match_level: str       # "variant", "gene", "none"
    civic_evidence: List[Dict]   # matched evidence items

def amp_assign_tier(
    classification: str,         # ACMG classification
    gene: str,
    hgvsp: str = "",
    strategy: str = "B",         # "A", "B", "C"
    civic_db_path: str = None,
) -> TierResult:
```

**Strategy B logic (default):**
1. CIViC variant-specific lookup (HGVSp → CIViC variant name)
2. If Level A Predictive found → Tier I, source="civic-variant-A"
3. If Level B Predictive found → Tier II, source="civic-variant-B"
4. OncoKB gene-level check (existing logic)
5. If OncoKB Level 1-2 + Path/LP → Tier I, source="oncokb-gene-L{level}"
6. If OncoKB Level 1-2 + VUS + hotspot → Tier II, source="hotspot"
7. If Level C-D + Path/LP → Tier II, source="civic-variant-{level}"
8. If Path/LP on cancer gene → Tier II, source="oncokb-gene-L{level}"
9. If VUS on cancer gene → Tier III, source="oncokb-gene-vus"
10. Else → Tier IV, source="default"

**Strategy A:** CIViC steps 1-3 only, then OncoKB fallback only if no CIViC match.
**Strategy C:** Skip CIViC steps entirely, OncoKB only (current behavior).

### 2. `scripts/db/query_civic.py` 수정

`get_variant_evidence()` 반환값에 `match_level` 추가:

```python
def get_variant_evidence(gene: str, variant_name: str = None) -> Dict:
    """Returns: {
        "match_level": "variant" | "gene" | "none",
        "evidence": [list of evidence items]
    }"""
```

새 함수 추가:

```python
def get_predictive_evidence_for_tier(
    gene: str, hgvsp: str, db_path: str = None
) -> Dict:
    """Tier 결정 전용. Predictive evidence만 반환, match_level 구분."""
```

### 3. `scripts/orchestrate.py` 수정

`_build_variant_records()` (line ~224)에서:

```python
# Before (current):
tier = assign_tier(classification, gene, clinvar_sig, hgvsp)

# After:
from scripts.somatic.amp_tiering import amp_assign_tier
tier_result = amp_assign_tier(
    classification=classification,
    gene=gene,
    hgvsp=hgvsp,
    strategy=config.get("somatic.tiering_strategy", "B"),
)
v_result["tier"] = tier_result.tier
v_result["tier_label"] = tier_result.tier_label
v_result["tier_evidence_source"] = tier_result.evidence_source
v_result["civic_match_level"] = tier_result.civic_match_level
v_result["civic_evidence"] = tier_result.civic_evidence
```

### 4. Cancer Report Template 변경

**Tier 섹션 라벨 변경:**
- "Tier 1 — Therapeutic Target" → "Tier I — Strong Clinical Significance"
- "Tier 2 — Clinically Significant" → "Tier II — Potential Clinical Significance"
- "Tier 3 — VUS in Cancer Gene" → "Tier III — Unknown Clinical Significance"
- Tier IV는 기존처럼 count만 표시

**Evidence 표시 추가:**
- CIViC evidence level badge (A/B/C/D/E) — 각 variant detail에
- Tier evidence source 표시 (e.g., "Based on: CIViC Level A — Vemurafenib sensitivity")
- Diagnostic/Prognostic evidence 별도 섹션

### 5. config.yaml 추가

```yaml
somatic:
  tiering_strategy: "B"  # "A" (CIViC priority), "B" (combined), "C" (OncoKB only)
  civic_tier_elevation: true  # CIViC evidence로 tier 상승 허용
```

---

## Rare Disease Mode — 영향 없음

Rare disease mode는 이 변경에 영향받지 않는다:
- `mode == "rare-disease"`일 때 `amp_assign_tier()` 호출하지 않음
- 기존 ACMG/AMP 2015 classification + HPO ranking 유지
- ClinGen validity, OMIM, HPO 로컬 DB는 그대로 사용

단, `orchestrate.py`에서 모드 분기를 명확히 하여 향후 rare disease에도 tier 시스템을 도입할 수 있는 구조를 유지한다:

```python
if mode == "cancer":
    tier_result = amp_assign_tier(...)
elif mode == "rare-disease":
    # ACMG classification + HPO ranking (existing)
    pass
```

---

## Data Flow Example

**Input:** BRAF p.Val600Glu, ACMG classification = Pathogenic

**Strategy B processing:**
1. HGVSp → CIViC name: `p.Val600Glu` → `V600E`
2. CIViC query: BRAF V600E → Level A Predictive evidence found
   - "Vemurafenib — Sensitivity in Melanoma (PMID: 20979469)"
   - match_level = "variant"
3. **Result: Tier I**, source = "civic-variant-A"
4. OncoKB check skipped (CIViC already provided Tier I)

**Input:** Novel TP53 variant (no CIViC entry), ACMG = Pathogenic

**Strategy B processing:**
1. HGVSp → CIViC name: no match
2. CIViC variant query: none found, match_level = "none"
3. OncoKB: TP53 is Level 1 gene + Pathogenic
4. **Result: Tier I**, source = "oncokb-gene-L1"

**Input:** KRAS G12D, ACMG = VUS

**Strategy B processing:**
1. CIViC: KRAS G12D → Level B Predictive evidence
2. But classification is VUS, Level B cannot elevate VUS to Tier I
3. Hotspot check: KRAS position 12 = hotspot
4. **Result: Tier II**, source = "hotspot"

---

## Testing Strategy

### Unit Tests (`tests/test_amp_tiering.py`)
- Strategy B: CIViC Level A + Path → Tier I
- Strategy B: CIViC Level B + Path → Tier II
- Strategy B: No CIViC + OncoKB L1 + Path → Tier I
- Strategy B: CIViC Level A + VUS → Tier II (not I — VUS cannot be Tier I via CIViC alone)
- Strategy B: VUS + hotspot → Tier II
- Strategy B: VUS + non-cancer gene → Tier IV
- Strategy B: Benign → Tier IV (CIViC cannot override)
- Strategy A: CIViC only
- Strategy C: OncoKB only (backward compatible with current)
- gene-level CIViC match does NOT elevate to Tier I

### Integration Tests
- Full pipeline with real CIViC DB: BRAF V600E → Tier I via CIViC
- Full pipeline with TP53 novel variant → Tier I via OncoKB
- Full pipeline rare disease mode → no tier assignment (unchanged)
- Strategy switching via config

### Regression
- All existing 383 tests pass
- Existing cancer report output validated

---

## Migration

1. `oncokb.py`의 `assign_tier()`는 deprecated wrapper로 유지 (내부에서 `amp_assign_tier(strategy="C")` 호출)
2. 기존 테스트는 Strategy C로 통과 (backward compatible)
3. 새 기본값은 Strategy B
