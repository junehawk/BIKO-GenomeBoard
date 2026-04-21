# Rare Disease & General Population Analysis Enhancements

> **2026-04-21 update**: The Korean-frequency expansion described in
> Feature 1 below (Korea4K AF + NARD2) is **superseded by KOVA v7**
> (KOGO / gene2korea Korean Variant Archive, 43 M variants with
> homozygote counts). The earlier multi-cohort sources are not used in
> the shipped pipeline. This spec is retained for historical reference
> only — see `docs/KOREAN_STRATEGY.md` for the active design.

**Date:** 2026-04-09
**Status:** Superseded (Korean-frequency path replaced by KOVA v7 on 2026-04-21)
**Scope:** 4 features for rare disease / germline analysis improvement

---

## Overview

BIKO GenomeBoard의 희귀질환 및 일반인 유전체 분석 역량을 강화하는 4가지 기능 구현.
현재 상태: SNV+CNV/SV+TMB 파이프라인 완성, 463 tests, main branch.

---

## Feature 1: Korea4K AF + NARD2 — 한국인 빈도 확장

### 목적
한국인 인구집단 대립유전자 빈도(AF) 소스를 KRGDB 1개 → 3개(+Korea4K, NARD2)로 확장하여 ACMG 빈도 판정(BA1/BS1/PM2) 정확도를 높인다.

### 설계

**데이터 형식:** KRGDB와 동일한 TSV 패턴 (chrom, pos, ref, alt, frequency)

**새 파일:**
- `scripts/korean_pop/query_korea4k.py` — TSV → in-memory dict, thread-safe 싱글톤
- `scripts/korean_pop/query_nard2.py` — 동일 패턴
- `data/korea4k_freq.tsv` — 샘플 데이터 (테스트용)
- `data/nard2_freq.tsv` — 샘플 데이터 (테스트용)
- `tests/test_korea4k.py`, `tests/test_nard2.py`

**수정 파일:**
- `scripts/korean_pop/compare_freq.py` — 3-tier → 5-tier 비교 확장
  - 현재: KRGDB / gnomAD EAS / gnomAD ALL
  - 확장: Korea4K / NARD2 / KRGDB / gnomAD EAS / gnomAD ALL
  - BA1/BS1/PM2 판정에 Korea4K, NARD2 빈도도 참조
  - FrequencyData 모델에 korea4k_af, nard2_af 필드 추가
- `scripts/orchestrate.py` — _query_variant_databases()에 Korea4K/NARD2 쿼리 추가
- `config.yaml` — paths.korea4k_freq, paths.nard2_freq 추가

### 빈도 판정 로직
```
korean_max = max(krgdb_af, korea4k_af, nard2_af)  # 한국인 소스 중 최대값
BA1: korean_max > 0.05 OR gnomad_all > 0.05
BS1: korean_max > 0.01 OR gnomad_all > 0.01
PM2_Supporting: korean_max < 0.001 AND gnomad_eas < 0.001 AND gnomad_all < 0.001
```

---

## Feature 2: OMIM genemap2.txt — 유전패턴 확장

### 목적
OMIM 유전자-질환-유전패턴 매핑을 11개 하드코딩 → genemap2.txt 전체(~4,500 유전자)로 확장. 희귀질환 유전 상담의 핵심 정보.

### 설계

**genemap2.txt 포맷** (tab-separated, `#` comment lines):
```
# Chromosome  Genomic Position Start  ...  Gene Symbols  Gene Name  MIM Number  ...  Phenotypes
chr1  1000000  2000000  ...  GENE1  Gene 1 Name  100000  ...  {Phenotype 1}, 100001 (3), Autosomal dominant; {Phenotype 2}, 100002 (3), Autosomal recessive
```

Phenotypes 컬럼에서 유전패턴 파싱:
- `Autosomal dominant` → AD
- `Autosomal recessive` → AR
- `X-linked` → XL
- `X-linked dominant` → XLD
- `X-linked recessive` → XLR
- `Mitochondrial` → MT
- `Digenic` → DIG
- `Somatic` → SOM

**새 파일:**
- `scripts/db/build_omim_genemap_db.py` — genemap2.txt → SQLite (omim_genemap.sqlite3)
  - Table: omim_genemap (gene, mim_number, phenotype, inheritance, phenotype_mim)
  - Table: metadata (source, build_date, record_count)
- `scripts/db/query_omim_genemap.py` — get_gene_phenotypes(gene) → list[dict]
- `data/db/sample_genemap2.txt` — 테스트용 샘플 (50개 유전자)
- `tests/test_omim_genemap.py`

**수정 파일:**
- `scripts/clinical/query_omim.py` — static OMIM_DATA dict → omim_genemap DB 우선 조회, static fallback
- `config.yaml` — paths.omim_genemap_db 추가
- `templates/rare-disease/report.html` — inheritance pattern 표시 강화 (AD/AR/XL 배지)

---

## Feature 3: In Silico Predictions — VEP CSQ 확장

### 목적
VEP CSQ에서 REVEL, CADD, AlphaMissense, SpliceAI 점수를 파싱하여 PP3/BP4 evidence code 생성 및 리포트 표시. VUS 재분류에 결정적 역할.

### 설계

**VEP CSQ 필드 매핑:**
| CSQ 필드 | 파싱 대상 | 타입 |
|----------|----------|------|
| REVEL_score | REVEL 점수 | float (0-1) |
| CADD_PHRED | CADD phred-scaled | float |
| am_class | AlphaMissense 분류 | str (likely_pathogenic/ambiguous/likely_benign) |
| am_pathogenicity | AlphaMissense 점수 | float (0-1) |
| SpliceAI_pred_DS_AG | SpliceAI acceptor gain | float (0-1) |
| SpliceAI_pred_DS_AL | SpliceAI acceptor loss | float (0-1) |
| SpliceAI_pred_DS_DG | SpliceAI donor gain | float (0-1) |
| SpliceAI_pred_DS_DL | SpliceAI donor loss | float (0-1) |

**새 파일:**
- `scripts/classification/in_silico.py` — InSilicoScores dataclass + PP3/BP4 evidence 생성 함수
- `tests/test_in_silico.py`

**수정 파일:**
- `scripts/intake/parse_annotation.py` — CSQ 파싱에 위 필드 추가
- `scripts/common/models.py` — Variant에 in_silico: dict 필드 추가
- `scripts/orchestrate.py` — _classify_variants()에서 in silico → PP3/BP4 evidence 수집
- `config.yaml` — in_silico thresholds 섹션 추가
- `templates/cancer/report.html` — detail page에 scores 표시
- `templates/rare-disease/report.html` — detail page에 scores 표시

**PP3/BP4 Threshold (ClinGen SVI 2022 권고 기반):**
```yaml
in_silico:
  pp3:
    revel_moderate: 0.644      # PP3_Moderate
    revel_strong: 0.932        # PP3_Strong
    spliceai_moderate: 0.2     # PP3_Moderate (splice)
    spliceai_strong: 0.5       # PP3_Strong (splice)
  bp4:
    revel_moderate: 0.183      # BP4_Moderate
    revel_strong: 0.016        # BP4_Strong
    spliceai_benign: 0.1       # BP4 (splice, max delta < 0.1)
```

**Evidence 생성 로직:**
1. Missense variant: REVEL 기반 (primary), CADD/AlphaMissense (보조 표시)
2. Splice region variant: SpliceAI 기반 (max of 4 delta scores)
3. PP3와 BP4는 상호 배타적 (둘 다 해당되면 더 강한 쪽 선택)

---

## Feature 4: InterVar ACMG — 분류 강화

### 목적
ACMG evidence 수집을 강화하여 분류 정확도 향상. 이중 전략: InterVar TSV 파싱(선택적) + 자체 규칙 확장.

### 설계

**A. InterVar TSV 파싱 (선택적 입력)**

InterVar 출력 포맷 (tab-separated):
```
#Chr  Start  End  Ref  Alt  ...  InterVar: Pathogenic  PVS1  PS1  PM1  PM2  ...  PP3  BP4  ...
```

- `scripts/intake/parse_intervar.py` — InterVar TSV 파서
  - `parse_intervar(path) → dict[variant_key, list[str]]` (variant → evidence codes)
- CLI: `--intervar <path>` 옵션 추가
- orchestrate.py: InterVar evidence를 _classify_variants()에서 AcmgEvidence에 추가

**B. 자체 Evidence 수집 강화**

InterVar 없이도 추가 수집 가능한 규칙:

| Evidence | 조건 | 소스 |
|----------|------|------|
| PVS1 | stop_gained, frameshift_variant in LOF-intolerant gene (pLI ≥ 0.9) | VEP consequence + gnomAD pLI |
| PM1 | variant in established functional domain (InterPro) | VEP DOMAINS 필드 |
| PM4 | in-frame insertion/deletion in non-repeat region | VEP consequence |
| PM5 | novel missense at same position as known pathogenic | ClinVar local DB |
| PP2 | missense in gene with low benign missense rate | gnomAD constraint (missense Z ≥ 3.09) |
| PP3 | computational evidence (pathogenic) | Feature 3에서 구현 |
| BP4 | computational evidence (benign) | Feature 3에서 구현 |
| BP7 | synonymous with no splice impact | VEP consequence + SpliceAI |

**새 파일:**
- `scripts/intake/parse_intervar.py` — InterVar TSV 파서
- `scripts/classification/evidence_collector.py` — 자체 evidence 수집 로직 (PVS1, PM1, PM4, PM5, PP2, BP7)
- `tests/test_intervar.py`
- `tests/test_evidence_collector.py`

**수정 파일:**
- `scripts/orchestrate.py` — CLI --intervar 추가, evidence 수집 파이프라인 확장
- `scripts/classification/acmg_engine.py` — evidence strength suffix 처리 확장 (_Moderate, _Strong, _VeryStrong)

**우선순위:** InterVar 파싱은 선택적 기능. 자체 evidence 수집이 핵심이며, PP3/BP4(Feature 3)와 PVS1이 가장 높은 임팩트.

---

## 공통 사항

### 테스트 전략
- 각 기능별 단위 테스트 (build, query, evidence 생성)
- 통합 테스트 (전체 파이프라인 run_pipeline with 새 기능)
- 기존 463 테스트 회귀 방지
- 샘플 VCF에 in silico CSQ 필드를 추가한 테스트 VCF 생성

### 실행 순서
- **Phase 1** (병렬): Korea4K + NARD2 + OMIM genemap2
- **Phase 2** (병렬): In silico predictions + InterVar ACMG
- Phase 2는 Phase 1과 독립적이므로, 팀 재구성 후 병렬 진행 가능

### 에이전트 배정
- Phase 1: db-dev (빌드/쿼리) + pipeline-dev (통합) + qa-engineer (테스트)
- Phase 2: pipeline-dev (파싱/evidence) + clinical-advisor (threshold 검증) + report-dev (리포트) + qa-engineer (테스트)
