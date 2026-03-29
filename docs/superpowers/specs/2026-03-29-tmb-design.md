# TMB (Tumor Mutational Burden) Implementation

## Goal

Cancer mode 리포트에 TMB를 자동 계산하여 표시한다. Nonsynonymous coding variant count / coding region Mb로 계산하며, 계산 방법을 리포트에 명시한다.

## Scope

**In scope:**
- TMB 계산 엔진 (`scripts/somatic/tmb.py`)
- Cancer mode 자동 계산, rare disease 스킵
- 리포트 TMB 뱃지 + Methodology 섹션 기록
- config.yaml로 threshold/consequence/panel size 설정 가능
- `--panel-size` / `--bed` CLI 옵션

**Out of scope:**
- MSI (향후 별도)
- TMB 기반 면역치료 자동 추천
- Tumor-normal pair differential TMB

---

## Architecture

```
VCF → parse_vcf() → variants
                       ↓
                  calculate_tmb(variants, panel_size)
                       ↓
                  TmbResult(score, level, count, panel_size, consequences)
                       ↓
                  report_data["tmb"] → 리포트 표시
```

---

## Data Model

### TmbResult

```python
@dataclass
class TmbResult:
    score: float              # mutations per Mb (e.g., 12.3)
    level: str                # "High", "Intermediate", "Low"
    variant_count: int        # number of counted variants
    total_variants: int       # total variants before filtering
    panel_size_mb: float      # coding region size used (e.g., 33.0)
    counted_consequences: List[str]  # which consequence types were counted
```

---

## Calculation Logic

### Variant Counting

VCF에서 파싱된 variant 리스트를 받아:
1. `consequence` 필드가 config의 `tmb_counted_consequences` 목록에 포함되는 것만 카운트
2. 기본 카운트 대상: missense_variant, stop_gained, frameshift_variant, splice_donor_variant, splice_acceptor_variant, inframe_deletion, inframe_insertion

### TMB Score

```
TMB = counted_variants / panel_size_mb
```

### Level Classification

| Level | Threshold | Color |
|-------|-----------|-------|
| **High** | ≥ 10 mut/Mb | Red (#991B1B) |
| **Intermediate** | ≥ 6 mut/Mb | Orange (#92400E) |
| **Low** | < 6 mut/Mb | Green (#166534) |

Threshold는 `config.yaml`에서 조정 가능.

---

## Configuration

```yaml
somatic:
  tiering_strategy: "B"
  tmb_high_threshold: 10
  tmb_intermediate_threshold: 6
  tmb_default_panel_size_mb: 33.0
  tmb_counted_consequences:
    - missense_variant
    - stop_gained
    - frameshift_variant
    - splice_donor_variant
    - splice_acceptor_variant
    - inframe_deletion
    - inframe_insertion
```

---

## CLI Options

```bash
# WGS (default: 33 Mb coding region)
python scripts/orchestrate.py sample.vcf -o report.html --mode cancer

# Custom panel size
python scripts/orchestrate.py sample.vcf -o report.html --panel-size 1.1

# BED file (auto-calculate panel size)
python scripts/orchestrate.py sample.vcf -o report.html --bed target_regions.bed
```

---

## Report Display

### Summary Page — TMB Badge

Cancer 리포트 summary page에 TMB 결과를 표시 (Genomic Findings 위 또는 아래):

```
┌──────────────────────────────────────┐
│ Tumor Mutational Burden              │
│ ████████████████░░░░  12.3 mut/Mb    │
│ TMB-High (≥10)                       │
│ 406 nonsynonymous variants / 33.0 Mb │
└──────────────────────────────────────┘
```

색상:
- TMB-High: 빨강 배경
- TMB-Intermediate: 주황 배경
- TMB-Low: 초록 배경

### Methodology Page — 계산 방법 기록

Methodology 섹션에 추가:

```
Tumor Mutational Burden (TMB)
TMB was calculated as the number of nonsynonymous coding variants
(missense, nonsense, frameshift, splice site, in-frame indel) per
megabase of coding region. Coding region size: 33.0 Mb (RefSeq exome).
TMB-High threshold: ≥10 mutations/Mb, consistent with the FDA-approved
FoundationOne CDx companion diagnostic methodology.
Counted consequence types: missense_variant, stop_gained, frameshift_variant,
splice_donor_variant, splice_acceptor_variant, inframe_deletion, inframe_insertion.
```

### JSON Output

```json
{
  "tmb": {
    "score": 12.3,
    "level": "High",
    "variant_count": 406,
    "total_variants": 777,
    "panel_size_mb": 33.0,
    "counted_consequences": ["missense_variant", "stop_gained", ...]
  }
}
```

---

## File Structure

### 신규 파일
| File | Responsibility |
|------|---------------|
| `scripts/somatic/tmb.py` | TMB 계산 엔진 |
| `tests/test_tmb.py` | TMB 단위 테스트 |

### 수정 파일
| File | Change |
|------|--------|
| `scripts/orchestrate.py` | `--panel-size`, `--bed` CLI, TMB 계산 호출 |
| `scripts/counselor/generate_pdf.py` | TMB 뱃지 렌더링 기본값 |
| `templates/cancer/report.html` | TMB 뱃지 섹션 + methodology 기록 |
| `config.yaml` | TMB threshold/consequences/panel size 설정 |

---

## Testing Strategy

### Unit Tests (`tests/test_tmb.py`)
- TMB-High (≥10): 330+ nonsynonymous / 33 Mb
- TMB-Intermediate (6-9.9): 200-329 / 33 Mb
- TMB-Low (<6): <200 / 33 Mb
- Empty variant list → score 0.0, level Low
- Custom panel size
- Consequence filtering (synonymous 제외 확인)
- BED file → panel size 계산

### Integration Tests
- Cancer pipeline TMB 자동 계산 확인
- Rare disease pipeline TMB 미계산 확인
- HTML 리포트에 TMB 뱃지 표시 확인
- JSON에 tmb 필드 포함 확인
- Methodology 섹션 기록 확인

### Regression
- 기존 449 tests 유지
