# 한국인 특화 분석 전략

BIKO GenomeBoard의 핵심 차별점은 한국인 집단 특이 데이터를 ACMG 분류에 통합하는 것이다. 이 문서는 그 전략과 구현 방식을 설명한다.

---

## 왜 한국인 특화인가?

### 대립유전자 빈도 차이

한국인 집단은 유럽인 중심으로 구축된 gnomAD 전체(ALL) 데이터와 유의미한 빈도 차이를 보이는 변이들이 존재한다. 유럽인 집단에서 희귀한 변이가 한국인에서는 흔하거나, 반대의 경우도 있다. gnomAD EAS(동아시아) 데이터도 있지만 한국인 단일 집단 데이터인 **KOVA v7** (Korean Variant Archive, KOGO / gene2korea 발표)가 더 정확한 한국인 빈도를 제공한다.

### Homozygote count의 임상적 의미

KOVA v7은 단순 대립유전자 빈도(AF) 뿐 아니라 각 변이의 **homozygote count**를 제공한다. 이 정보는 상염색체 열성(AR) 질환 해석에 결정적이다:

- gnomAD도 homozygote count를 제공하지만 한국인 특이적 AR 변이가 과소대표될 수 있다
- KOVA에서 homozygote가 다수 관찰되면 해당 변이는 고빈도 AR 보인자 변이일 가능성이 높다 (BS2 계열 증거)
- 반대로 homozygote가 0이면서 heterozygote만 관찰되는 경우 AR에서 병원성 가능성이 높아진다

### 약물유전체(PGx) 대사 차이

주요 PGx 유전자의 기능 변이 분포가 동아시아인과 유럽인 사이에 다르다. 예를 들어:
- `CYP2C19*2`, `*3`는 동아시아인에서 훨씬 높은 빈도로 나타난다
- `NUDT15*3`는 한국인에서 티오퓨린 독성과 강하게 연관된다
- `HLA-B*58:01`은 한국인에서 알로퓨리놀 과민반응의 주요 위험 인자다

---

## 3-tier 빈도 비교 전략 (KOVA v7 기반)

`scripts/population/compare_freq.py`가 구현하는 3-tier 한국인-인지 빈도 비교:

```
Tier 1: KOVA v7       →  한국인 Korean Variant Archive (43M variants, homozygote counts 포함)
Tier 2: gnomAD EAS    →  동아시아 집단 빈도
Tier 3: gnomAD ALL    →  전세계 집단 빈도
```

세 출처를 동시에 비교하고, 가용한 최대 빈도값을 ACMG 분류 기준(BA1/BS1/PM2)에 적용한다. 한국인 단독 코호트인 **KOVA v7** 빈도가 있으면 gnomAD 보다 우선시한다. Korean enrichment ratio(한국인 빈도 ÷ gnomAD ALL 빈도)는 각 변이마다 자동으로 계산되어 리포트에 반영된다.

> **과거 멀티-코호트 전략에서의 이행**: 프로젝트 초창기 설계에서는 한국인 단독 코호트 여러 종을 병렬 참조하도록 설계했었다. 공개 접근이 안정적으로 보장되지 않거나 배포 조건이 제한적이어서, 2026-04-21부터 KOGO / gene2korea가 발표한 **Korean Variant Archive v7 (KOVA v7)** 단일 소스로 통합되었다. KOVA v7은 43M variants와 homozygote count를 함께 제공하며, 이전 한국인 단독 코호트들이 커버하던 판정 시나리오를 모두 대체한다. 과거 설계에 대한 감사 추적은 `docs/superpowers/specs/` 아래 설계 스펙 파일에서 "superseded by KOVA v7" 주석과 함께 보존된다.

---

## ACMG 빈도 임계값

| ACMG 코드 | 방향 | 빈도 기준 | 의미 |
|-----------|------|-----------|------|
| `BA1` | Benign (Stand-alone) | > 5% | 매우 흔한 변이 — 단독으로 양성 판정 |
| `BS1` | Benign (Strong) | ≥ 1% | 흔한 변이 |
| `BS2` | Benign (Strong) | 관찰된 homozygote 다수 | AR / XL 질환에서 건강한 homozygote가 관찰 |
| `PM2_Supporting` | Pathogenic (Supporting) | ≤ 0.1% | 희귀 변이 |

> **한국인 특화 플래그**: KOVA v7 빈도가 gnomAD ALL 대비 5배 이상 높으면 "한국인 빈도 글로벌 대비 5배 이상 높음" 플래그가 추가된다. 이는 유럽인 기준 데이터베이스에서 과소평가될 수 있는 변이를 식별한다.

> **Homozygote 플래그**: KOVA v7에서 homozygote가 다수 관찰되는 변이는 AR 질환 해석 시 BS2 증거 후보로 표시된다. 단독 판정은 아니며, 표현형/유전 패턴 확인이 병행되어야 한다.

---

## 한국인 PGx 5대 유전자

`scripts/pharmacogenomics/korean_pgx.py`가 관리하는 한국인 특화 PGx 유전자:

### 1. CYP2D6 — 마약성 진통제, 항우울제

한국인에서 poor metabolizer(*4, *5)와 ultra-rapid metabolizer 빈도가 유럽인과 다르다. 코데인, 트라마돌 용량 조절에 영향.

### 2. CYP2C19 — 클로피도그렐, 프로톤펌프억제제

`*2`(rs4244285)와 `*3`(rs4986893)가 동아시아인에서 매우 흔하다. 클로피도그렐 내성, 항혈소판 치료 실패와 연관.

### 3. CYP2C9 — 와파린, NSAIDs

`*3`(rs1057910) 빈도가 한국인에서 유럽인 대비 낮지만, 와파린 용량 계산 시 반드시 고려해야 한다.

### 4. HLA-B — 약물 과민반응

`HLA-B*58:01`: 알로퓨리놀 유발 Stevens-Johnson 증후군/중독성 표피괴사용해(SJS/TEN). 한국인에서 빈도 약 6–7%.

`HLA-B*15:02`: 카르바마제핀 유발 SJS/TEN. 한국인보다 동남아시아인에서 더 흔하지만 한국인에서도 모니터링 필요.

### 5. NUDT15 — 티오퓨린 (아자티오프린, 6-MP)

`NUDT15*3`(c.415C>T, p.Arg139Cys)는 동아시아인 특이 변이로, 유럽인 데이터베이스에 과소대표되어 있다. 한국인 IBD, ALL 환자 티오퓨린 투약 시 반드시 스크리닝해야 한다.

---

## KOVA v7 데이터 관리

### 데이터 형식

KOVA v7 빈도 테이블은 `config.yaml`의 `paths.kova`로 등록된다. 레코드는 대립유전자 빈도(AF)와 homozygote count를 함께 제공한다.

```
chrom   pos      ref   alt   af_korean    hom_count    het_count
chr1    925952   G     A     0.0023       0            12
chr17   43057051 A     G     0.00041      0            2
...
```

### 조회 방식

`scripts/population/query_kova.py`가 KOVA v7 로컬 소스를 직접 읽어 `KovaRecord(af, hom_count, het_count, total_an)` 형태로 반환한다. 외부 API 의존 없이 빠른 응답이 가능하다.

### 업데이트 절차

1. KOGO / gene2korea 공식 KOVA 배포본 (v7 이상) 입수
2. 배포 VCF → 내부 TSV로 변환:
   ```bash
   bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC_Hom\t%INFO/AC_Het\t%INFO/AN\n' \
     kova_v7.vcf.gz > data/kova_freq.tsv
   ```
3. 기존 파일 교체 후 테스트 실행:
   ```bash
   pytest tests/test_kova.py -v
   ```

---

## 향후 개발 계획

| 항목 | 설명 |
|------|------|
| **KOGES 통합** | 한국인유전체역학조사사업(KoGES) 데이터 연동으로 표현형-유전형 상관관계 분석 강화 |
| **Korean HGMD 큐레이션** | 한국인 병원에서 보고된 변이 데이터베이스 구축 및 통합 |
| **한국인 암 변이** | KCDC 암 게놈 데이터와의 연동 |
| **SpliceAI 한국인 보정** | 스플라이싱 예측 모델의 한국인 특이 보정 |
| **KOVA 상위 버전 추적** | KOGO / gene2korea의 후속 릴리스(v8 이상) 추적 및 자동 업데이트 경로 |
