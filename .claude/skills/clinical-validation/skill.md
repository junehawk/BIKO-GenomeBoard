---
name: clinical-validation
description: "GenomeBoard의 변이 분류, tiering, 임상적 의사결정에 대한 가이드라인 기반 검증 절차. ACMG/AMP 2015, AMP/ASCO/CAP 2017, ACMG CNV 2020 가이드라인 요약과 검증 체크리스트. '임상 검토', '가이드라인', 'ACMG', 'AMP', 'tiering', '분류 로직', 'clinical review', '임상적 의미', '임상 자문' 관련 작업 시 사용."
---

# Clinical Validation — 임상 검증 절차

GenomeBoard의 모든 변이 분류와 tiering이 국제 가이드라인에 부합하는지 검증하는 절차.

## 가이드라인 요약

### 1. ACMG/AMP 2015 (Germline)

**적용 대상:** rare-disease 모드의 SNV/Indel 분류

**5단계 분류:** Pathogenic → Likely Pathogenic → VUS → Likely Benign → Benign

**근거 코드 체계:**
- Pathogenic: PVS1, PS1-4, PM1-6, PP1-5 (Very Strong → Supporting)
- Benign: BA1, BS1-4, BP1-7 (Stand-alone → Supporting)

**한국인 특화 빈도 기준:**
- BA1: AF > 5% (gnomAD ALL 또는 KRGDB)
- BS1: AF > 1%
- PM2_Supporting: AF < 0.1% (gnomAD EAS 또는 KRGDB 기준)

**ClinVar override 규칙:**
- Expert panel/practice guideline → 직접 분류 적용 (ACMG 엔진 우회)
- Multiple submitters (no conflict) + ≥2 stars → 분류 적용
- Conflict 있으면 override 하지 않음

### 2. AMP/ASCO/CAP 2017 (Somatic)

**적용 대상:** cancer 모드의 SNV/Indel tiering

**4단계 분류:**
| Tier | 의미 | 임상 활용 |
|------|------|----------|
| I | Strong clinical significance | FDA 승인 치료제/동반진단 |
| II | Potential clinical significance | 임상시험/가이드라인 기반 치료 |
| III | Unknown significance | 임상적 의미 불확실 |
| IV | Benign/Likely benign | 양성 변이 |

**GenomeBoard 전략 (현재 B):**

| 우선순위 | 조건 | Tier |
|---------|------|------|
| 1 | ClinVar Pathogenic/LP (Expert Panel) | I |
| 2 | CIViC variant-level Evidence A | I |
| 3 | CIViC variant-level Evidence B | II |
| 4 | OncoKB Level 1-2 | I |
| 5 | CIViC C-D + ClinVar Path/LP | II |
| 6 | OncoKB Level 3A-4 | II |
| 7 | Cancer hotspot | II |
| 8 | Cancer gene VUS | III |
| 9 | 기타 | IV |

**전략 전환 조건:**
- A→B: CIViC DB 규모가 충분히 클 때 (현재 상태)
- B→C: OncoKB 라이선스 비용 문제 발생 시
- B→A: 특정 임상 파트너가 OncoKB 단독 사용 요구 시

### 3. ACMG CNV 2020 (구조변이)

**적용 대상:** CNV/SV (AnnotSV 결과)

**5단계 분류:** Class 5 (Pathogenic) → Class 4 (Likely Pathogenic) → Class 3 (VUS) → Class 2 (Likely Benign) → Class 1 (Benign)

**리포트 표시 전략:**
- Class 4-5: 주요 소견(Findings) 패널에 표시
- Class 3 + dosage-sensitive: 별도 테이블에 ⚠ DS 플래그와 함께 표시
- Class 3 (non-dosage): 통계만 (건수)
- Class 1-2: 통계만

**Dosage sensitivity 기준:**
- ClinGen HI/TS score ≥ 2
- pLI ≥ 0.9 (haploinsufficiency)
- Loss-of-function observed/expected upper bound ≤ 0.35

## 검증 체크리스트

### 분류 로직 검증

- [ ] 모든 ACMG evidence code가 올바르게 수집되는가?
- [ ] 한국인 빈도 기준(BA1/BS1/PM2)이 config.yaml의 thresholds와 일치하는가?
- [ ] ClinVar expert panel override가 conflict 없는 경우에만 적용되는가?
- [ ] PGx 유전자(CYP2C19, CYP2D6 등)와 APOE가 ACMG 분류에서 제외되는가?

### Tiering 로직 검증

- [ ] 현재 전략(B)의 우선순위 체인이 올바르게 구현되었는가?
- [ ] CIViC variant-level evidence와 gene-level evidence가 구분되는가?
- [ ] match_level (variant/gene/none)이 정확히 추적되는가?
- [ ] cancer 모드에서만 AMP tiering이 적용되는가?
- [ ] rare-disease 모드에서 CIViC enrichment가 호출되지 않는가?

### 리포트 검증

- [ ] VUS가 기본적으로 숨겨지는가? (hide_vus=True)
- [ ] Tier 순서(I→II→III→IV)로 정렬되는가?
- [ ] CNV Class 4-5가 주요 소견에 표시되는가?
- [ ] TMB가 cancer 모드에서만 계산/표시되는가?
- [ ] Methodology 섹션에 계산 방법이 기록되는가?

### 새 기능 임상 평가 기준

새 데이터 소스나 분석 기능 추가 시:
1. **임상적 타당성** — 이 정보가 환자 관리에 실질적 영향을 주는가?
2. **근거 수준** — 어떤 수준의 근거(Level A-D)로 뒷받침되는가?
3. **가이드라인 정합성** — 관련 학회 가이드라인과 일치하는가?
4. **기존 체계 영향** — 현재 분류/tiering 체계를 방해하지 않는가?
5. **오류 가능성** — false positive/negative의 임상적 위험은?
