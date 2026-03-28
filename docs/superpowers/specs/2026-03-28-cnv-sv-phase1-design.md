# CNV/SV Phase 1: AnnotSV Parser + Report Integration

## Goal

AnnotSV TSV 출력을 파싱하여 CNV/SV 결과를 기존 SNV 리포트에 통합한다. AnnotSV가 제공하는 ACMG CNV classification (Class 1-5)을 그대로 사용하며, 자체 분류 엔진은 구축하지 않는다 (Phase 2).

## Scope

**In scope:**
- `StructuralVariant` dataclass
- AnnotSV TSV 파서 (`parse_annotsv.py`)
- CLI `--sv` 옵션 추가
- 리포트에 CNV/SV 섹션 추가 (cancer + rare disease 양쪽)
- Class별 3단계 표시 전략
- Dosage-sensitive VUS 필터링 + 플래그

**Out of scope (Phase 2+):**
- 독자적 ACMG CNV 2020 분류 엔진
- Fusion 전용 처리 (breakend → gene fusion detection)
- ClinGen dosage sensitivity API 자체 조회
- Family trio CNV analysis

---

## Architecture

```
기존 SNV 파이프라인 (변경 없음):
  VCF → parse_vcf() → Variant → ACMG/AMP → Tier I-IV → Report SNV 섹션

신규 SV 파이프라인 (별도):
  AnnotSV TSV → parse_annotsv() → StructuralVariant → Class 분류 (AnnotSV 제공)
                                                         ↓
                                                    Report SV 섹션

orchestrate.py:
  --sv annotsv.tsv 옵션 → parse_annotsv() 호출
  SNV report_data + SV report_data → 하나의 리포트 생성
```

---

## Data Model

### `StructuralVariant` dataclass

```python
@dataclass
class StructuralVariant:
    annotsv_id: str
    chrom: str
    start: int
    end: int
    length: int
    sv_type: str              # DEL, DUP, INV, BND, INS
    sample_id: str
    acmg_class: int           # 1-5
    ranking_score: float      # -0.99 to +0.99
    cytoband: str
    gene_name: str            # "ERBB2" or "TBX1;COMT;HIRA"
    gene_count: int
    # Gene detail (from split rows)
    gene_details: List[Dict]  # [{gene, transcript, cds_percent, frameshift, location, hi, ts, pli}]
    # Pathogenic evidence
    p_gain_phen: str
    p_gain_hpo: str
    p_gain_source: str
    p_loss_phen: str
    p_loss_hpo: str
    p_loss_source: str
    # Benign evidence
    b_gain_af_max: Optional[float]
    b_loss_af_max: Optional[float]
    # OMIM
    omim_morbid: bool
```

---

## AnnotSV TSV Parser

### `scripts/intake/parse_annotsv.py`

- TSV 파싱 (tab-delimited 또는 pipe-delimited — 실 데이터에서 구분자 자동 감지)
- **full rows** → `StructuralVariant` 1개 생성
- **split rows** → 해당 SV의 `gene_details` 리스트에 추가
- 파서는 full row를 기준으로 SV를 생성하고, split rows에서 gene-level detail을 수집

```python
def parse_annotsv(tsv_path: str) -> List[StructuralVariant]:
    """Parse AnnotSV TSV output into StructuralVariant objects.
    Returns one StructuralVariant per SV (from full rows),
    enriched with gene details (from split rows).
    """
```

---

## Class별 표시 전략

### Class 4-5 (Pathogenic / Likely Pathogenic)

**표시:** 상세 블록 (SNV detail page와 유사)
- SV 정보: type, coordinates, size, cytoband
- Gene overlap: gene별 transcript, CDS%, frameshift, HI/TS, pLI
- Evidence: pathogenic sources, HPO terms, phenotypes
- Benign overlap AF (있으면)

### Class 3 (VUS) — Dosage-Sensitive 필터링

**표시:** 필터링된 요약 테이블

**테이블에 포함 조건 (하나라도 해당):**
- DEL이고 `HI >= 2` (ClinGen haploinsufficiency)
- DUP이고 `TS >= 2` (ClinGen triplosensitivity)
- `pLI >= 0.9` (gnomAD loss-of-function intolerance)
- `length > 1,000,000` (>1Mb) AND `gene_count >= 3` AND `omim_morbid = yes`
- Cancer mode: `gene_name`이 OncoKB cancer gene list에 포함

**Rare disease mode:** threshold를 낮춤 (HI >= 1, pLI >= 0.8)

**테이블 컬럼:**
```
Cytoband | Gene(s) | Type | Size | CDS% | HI/TS | pLI | OMIM | AF | Flag
```

**Dosage-sensitive flag:** 주황색 뱃지 `⚠ DS` (Dosage Sensitive)

**나머지 Class 3:** count만 표시 ("X additional VUS structural variants not shown")

**테이블 하단 disclaimer:**
"Classified as VUS per ACMG/ClinGen CNV criteria. Not actionable at this time."

### Class 1-2 (Benign / Likely Benign)

**표시:** count 한 줄만
"X benign/likely benign structural variant(s) identified but not reported."

---

## Report Template Layout

### Summary Page

```
[기존 SNV Genomic Findings — Tier I-IV]

── Structural Variants / Copy Number Alterations ──

  Pathogenic/Likely Pathogenic (Class 4-5):
  ┌──────────────────────────────────────────────┐
  │ ERBB2  DUP  17q12  40.5kb  ACMG Class 5     │
  │ CDKN2A DEL  9p21.3 26.7kb  ACMG Class 5     │
  └──────────────────────────────────────────────┘

  VUS — Dosage Sensitive (Class 3):
  ┌────────────────────────────────────────────────┐
  │ Cytoband│Gene │Type│Size │CDS%│HI/TS│pLI│Flag │
  │ 3q26.32 │PIK3CA│DUP│41kb│100%│0/0 │0.02│    │
  └────────────────────────────────────────────────┘
  2 additional VUS structural variants not shown.

  3 benign structural variants identified but not reported.
```

### Detail Pages (Class 4-5 only)

각 Pathogenic/LP SV가 별도 상세 블록을 가짐:

```
┌─ ERBB2 Amplification ──────────────────────────┐
│ Type: DUP | chr17:37,844,393-37,884,925        │
│ Size: 40,532 bp | Cytoband: 17q12              │
│ ACMG Class: 5 (Pathogenic) | Score: 0.99       │
│                                                 │
│ Gene Overlap:                                   │
│  ERBB2 (NM_004448.4) — 100% CDS, in-frame     │
│  HI: 2  TS: 3  pLI: 0.15                      │
│                                                 │
│ Associated Phenotypes:                          │
│  Breast cancer; Gastric cancer                  │
│  HPO: HP:0003002; HP:0012126                   │
│                                                 │
│ Evidence: ClinVar; ClinGen                      │
└─────────────────────────────────────────────────┘
```

### SV Type 색상 코드

| Type | 색상 | 의미 |
|------|------|------|
| DEL | 빨강 (#991B1B) | 결실 — copy number loss |
| DUP | 파랑 (#1E40AF) | 증폭 — copy number gain |
| INV | 회색 (#475569) | 역위 |
| BND | 주황 (#92400E) | Breakend / 전좌 |
| INS | 보라 (#6B21A8) | 삽입 |

---

## CLI Interface

```bash
# SNV + SV 통합 리포트
python scripts/orchestrate.py sample.vcf -o report.html --sv annotsv_output.tsv

# SV만 (SNV 없이)
python scripts/orchestrate.py --sv-only annotsv_output.tsv -o sv_report.html

# 옵션
--sv <path>           AnnotSV TSV 파일 경로
--sv-only <path>      SV만 리포트 (VCF 없이)
--sv-show-all-vus     Class 3 VUS 전부 표시 (dosage filter 무시)
```

---

## File Structure

### 신규 파일

| File | Responsibility |
|------|---------------|
| `scripts/intake/parse_annotsv.py` | AnnotSV TSV → List[StructuralVariant] |
| `scripts/common/models.py` 수정 | StructuralVariant dataclass 추가 |
| `templates/shared/sv_section.html` | CNV/SV 리포트 섹션 (양쪽 모드 공용) |
| `tests/test_parse_annotsv.py` | 파서 테스트 |
| `tests/test_sv_report.py` | SV 리포트 생성 테스트 |

### 수정 파일

| File | Change |
|------|--------|
| `scripts/orchestrate.py` | `--sv` CLI 옵션, SV 파싱 + report_data 통합 |
| `scripts/counselor/generate_pdf.py` | SV 섹션 렌더링, dosage filter 로직 |
| `templates/cancer/report.html` | SV 섹션 include |
| `templates/rare-disease/report.html` | SV 섹션 include |

---

## Testing Strategy

### Unit Tests
- `test_parse_annotsv.py`: cancer/rare disease TSV 파싱, full/split 매핑, 빈 파일 처리
- StructuralVariant 모델 필드 검증
- Dosage-sensitive filter 로직 (HI/TS/pLI threshold)

### Integration Tests
- 전체 파이프라인: VCF + AnnotSV TSV → 통합 리포트 HTML
- SV-only 모드: AnnotSV TSV만으로 리포트 생성
- 기존 SNV-only 파이프라인 회귀 없음 (--sv 옵션 없이 기존 동작 동일)

### Test Data
- `data/sample_sv/cancer_somatic_annotsv.tsv` (8 SVs: Class 5×4, 4×2, 3×1, 1×1)
- `data/sample_sv/rare_disease_annotsv.tsv` (7 SVs: Class 5×4, 3×2, 1×1)
