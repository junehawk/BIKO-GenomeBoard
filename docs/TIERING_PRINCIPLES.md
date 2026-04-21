# BIKO GenomeBoard Variant Tiering Principles

> **PM1 hotspot + ClinVar reconciliation note** — Cancer 모드 ACMG engine에 PMID 기반 PM1 hotspot table(`data/pm1_hotspot_domains.json`, TP53/KRAS/NRAS/BRAF/PIK3CA/EGFR/IDH1-2 등)이 포함되어 VEP `DOMAINS`가 비어 있어도 PM1이 발화합니다. 좁은 4-조건 ClinVar override(`engine ≥ LP` + ClinVar Conflicting + PM1 + PM5)가 post-classify seam에 배치되어 canonical hotspot이 Conflicting submitters 때문에 VUS로 떨어지는 것을 방지합니다. Variant Selector는 protein-impacting consequence gate(P/LP bypass + SpliceAI ≥ 0.2 rescue)와 MMR/Lynch carve-out을 포함합니다. 설계 근거는 [`CLAUDE.md`](../CLAUDE.md) 참조.

## 개요

BIKO GenomeBoard는 두 가지 분석 모드에 따라 서로 다른 분류 체계를 사용한다.

| 모드 | 가이드라인 | 분류 체계 | 목적 |
|------|-----------|----------|------|
| **Cancer (Somatic)** | AMP/ASCO/CAP 2017 | Tier I–IV | 치료 결정 지원 |
| **Rare Disease (Germline)** | ACMG/AMP 2015 | Pathogenic → Benign (5단계) + HPO ranking | 진단 지원 |

---

## Cancer Mode: AMP/ASCO/CAP 2017 Somatic Tiering

### Tier 정의

| Tier | 의미 | 근거 수준 | 임상 행동 |
|------|------|----------|----------|
| **Tier I** | Strong Clinical Significance | FDA 승인 치료제, 전문 가이드라인 포함 | 치료 결정에 직접 활용 |
| **Tier II** | Potential Clinical Significance | 임상시험, 소규모 연구, 전문가 합의 | 임상시험 등록 고려, 전문가 논의 |
| **Tier III** | Unknown Clinical Significance | 임상적 의의 미확인 | 추가 연구 필요, 재분석 대상 |
| **Tier IV** | Benign or Likely Benign | 양성 근거 확보 | 임상적 행동 불필요 |

### Tiering 전략: Modified Approach B

두 가지 데이터 소스를 종합하여 tier를 결정한다:

1. **OncoKB (gene-level)** — 유전자 수준의 임상 근거 등급 (Level 1–4)
2. **CIViC (variant-level)** — 변이 수준의 치료 근거 (Evidence Level A–E)

#### 핵심 원칙

```
원칙 1: CIViC는 tier를 올릴 수만 있고, 내릴 수 없다.
원칙 2: Tier I 상승은 CIViC variant-specific 매치 필수 (gene-level fallback 불가).
원칙 3: Predictive evidence만 tier 결정에 반영한다.
원칙 4: 모든 tier 결정에는 근거 출처(evidence source)를 기록한다.
```

#### Tier 결정 우선순위

```
우선순위 1: CIViC variant-specific Level A Predictive     → Tier I
우선순위 2: CIViC variant-specific Level B Predictive     → Tier II
우선순위 3: OncoKB Level 1-2 + Pathogenic/LP              → Tier I
우선순위 4: OncoKB Level 1-2 + VUS + hotspot              → Tier II
우선순위 5: CIViC variant-specific Level C-D + Path/LP    → Tier II
우선순위 6: Pathogenic/LP on any cancer gene              → Tier II
우선순위 7: VUS on cancer gene                            → Tier III
우선순위 8: 나머지 (Benign/LB, 비암유전자 VUS)              → Tier IV
```

#### 전략 전환 조건

시스템은 config.yaml의 `tiering_strategy` 설정으로 A/B/C 전략을 전환할 수 있다:

| 전략 | 설명 | 전환 조건 |
|------|------|----------|
| **A** | CIViC evidence 우선 | CIViC 커버리지가 타겟 패널의 90%+ 확인 시 |
| **B** | OncoKB + CIViC 종합 (기본값) | 현재 권장 전략 |
| **C** | OncoKB 유지 + CIViC 표시만 | MFDS 등 규제기관이 고정 lookup table 요구 시 |

### Evidence 매칭 규칙

CIViC에서 변이를 찾을 때:

1. **Variant-specific match** — HGVSp를 CIViC variant name으로 변환하여 정확히 매칭
   - 예: `p.Val600Glu` → `V600E` → CIViC에서 BRAF V600E 검색
   - Tier I 상승에 사용 가능

2. **Gene-level fallback** — variant 매칭 실패 시 유전자 수준 evidence 조회
   - Tier 상승에 사용 불가
   - 리포트에 참고 정보로 표시

### Evidence 유형별 용도

| Evidence Type | Tier 결정 | 리포트 표시 | 설명 |
|--------------|----------|-----------|------|
| **Predictive** | O | O | 약물 반응/내성 예측 |
| **Diagnostic** | X | O | 진단 보조 |
| **Prognostic** | X | O | 예후 예측 |
| **Oncogenic** | X | O | 종양 유발 근거 |

---

## Rare Disease Mode: ACMG/AMP 2015 Germline Classification

### 분류 체계

| Classification | 의미 | 임상 행동 |
|---------------|------|----------|
| **Pathogenic** | 병원성 | 진단 확정, 유전 상담 |
| **Likely Pathogenic** | 병원성 가능 | 임상 관리 시작, 추가 검증 |
| **VUS** | 의미 불확실 | 재분석 대상, 임상 결정에 사용 불가 |
| **Likely Benign** | 양성 가능 | 임상적 행동 불필요 |
| **Benign** | 양성 | 임상적 행동 불필요 |

### Candidate Ranking

ACMG 분류 후, 다음 기준으로 후보 유전자를 순위화한다:

```
1차 정렬: Classification (Pathogenic > LP > VUS > LB > Benign)
2차 정렬: HPO Score (환자 표현형과 유전자 연관도, 내림차순)
3차 정렬: ClinGen Gene-Disease Validity (Definitive > Strong > Moderate > ...)
```

### 보조 데이터

| 데이터 | 출처 | 용도 |
|--------|------|------|
| HPO gene-phenotype | 로컬 SQLite (329K associations) | 후보 순위화 |
| ClinGen validity | 로컬 SQLite | 유전자-질환 연관 신뢰도 |
| OMIM | 정적 데이터 + API | 질환 표현형, 유전 패턴 |
| Korean frequency (KOVA v7) | 로컬 AF + homozygote count (43M variants) | 한국인 특이 변이 플래그, AR BS2 후보 |

---

## 공통 원칙

### 1. 로컬 데이터 우선

온프레미스 환경에서 API 의존 없이 동작해야 한다. 모든 핵심 데이터(ClinVar, gnomAD, CIViC, HPO, ClinGen)는 로컬 DB로 운영한다.

### 2. 근거 추적성 (Evidence Traceability)

모든 분류/tier 결정에는 근거 출처를 기록한다:
- Cancer: `tier_evidence_source` 필드 (예: `"civic-variant-A"`, `"oncokb-gene-L1"`)
- Rare Disease: `acmg_codes` 리스트 (예: `["PS1", "PM2_Supporting"]`)

### 3. 한국인 특이성

- **KOVA v7** (Korean Variant Archive, KOGO / gene2korea) 빈도로 한국인 특이 변이 플래그 — 43M variants, homozygote count 포함
- gnomAD EAS vs ALL 비교로 동아시아 빈도 이상 탐지
- 한국인 PGx 12개 유전자 특이 유병률 비교
- KOVA v7 homozygote count는 AR 질환 판정 시 BS2 후보 플래그로 활용

### 4. 보수적 분류

불확실할 때는 낮은 tier / 약한 분류를 부여한다. Over-tiering(과분류)은 불필요한 치료로 이어질 수 있고, under-tiering(저분류)은 재분석에서 발견될 수 있다.

---

## 참고 문헌

- Li MM, et al. Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer. *J Mol Diagn.* 2017;19(1):4-23. (AMP/ASCO/CAP 2017)
- Richards S, et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med.* 2015;17(5):405-424. (ACMG/AMP 2015)
- Chakravarty D, et al. OncoKB: A Precision Oncology Knowledge Base. *JCO Precis Oncol.* 2017.
- Griffith M, et al. CIViC is a community knowledgebase for expert crowdsourcing the clinical interpretation of variants in cancer. *Nat Genet.* 2017;49(2):170-174.
