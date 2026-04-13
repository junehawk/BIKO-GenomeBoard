# API 키 레퍼런스

BIKO GenomeBoard가 사용하는 외부 서비스 및 API 키 정보를 정리한다.

---

## 환경 변수 목록

`.env.example`을 복사해 `.env`를 생성하고 아래 값을 채운다.

```bash
cp .env.example .env
```

---

## NCBI_API_KEY

- **필수 여부**: 선택 (권장)
- **환경 변수**: `NCBI_API_KEY`

### 속도 제한

| 상태 | 요청 한도 |
|------|-----------|
| 키 없음 | 3 req/sec |
| 키 있음 | 10 req/sec |

대량 변이 분석 시 키 없이 사용하면 `HTTP 429` 오류가 빈번하게 발생한다. 키 발급을 권장한다.

### 키 발급 방법

1. https://www.ncbi.nlm.nih.gov/account/settings/ 접속
2. NCBI 계정 로그인 (없으면 무료 가입)
3. "API Key Management" 섹션에서 키 생성
4. `.env`에 입력:

```
NCBI_API_KEY=your_key_here
```

### 사용 위치

`scripts/clinical/query_clinvar.py` — ClinVar E-utilities 호출 시 사용.

---

## gnomAD

- **필수 여부**: 불필요 (API 키 없음)
- **환경 변수**: 없음

gnomAD는 공개 GraphQL API를 제공하며 인증이 필요하지 않다.

### 주의사항

- **배치 쿼리 성능**: gnomAD API는 배치 쿼리 시 응답이 느릴 수 있다 (변이당 1–3초).
- **대규모 분석**: 100개 이상의 변이를 처리할 경우 [gnomAD 로컬 데이터 다운로드](https://gnomad.broadinstitute.org/downloads)를 권장한다.

```bash
# gnomAD v4 exomes VCF 다운로드 예시
gsutil cp gs://gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr17.vcf.bgz .
```

### 사용 위치

`scripts/korean_pop/query_gnomad.py` — EAS 및 ALL 집단 빈도 조회.

---

## PharmGKB

- **필수 여부**: 불필요 (API 키 없음)
- **환경 변수**: 없음

PharmGKB REST API는 공개 접근을 허용한다.

### 속도 제한

PharmGKB는 공식 Rate Limit을 명시하지 않지만, 코드 내에서 요청 간 **0.5초 간격**을 강제한다 (`scripts/pharma/query_pharmgkb.py`).

```python
# 코드 내 적용 예시
time.sleep(0.5)  # PharmGKB rate limit 준수
```

대규모 PGx 분석 시 로컬 PharmGKB 데이터베이스 다운로드를 고려할 수 있다: https://www.pharmgkb.org/downloads

### 사용 위치

`scripts/pharma/query_pharmgkb.py` — 약물-유전자 상호작용, CPIC 가이드라인 조회.

---

## ANTHROPIC_API_KEY

- **필수 여부**: **필수**
- **환경 변수**: `ANTHROPIC_API_KEY`

Claude 기반 에이전트 구동에 필요하다.

### 해당 에이전트

| 에이전트 | 역할 |
|---------|------|
| CEO | 사용자 소통, 분석 요청 접수 |
| CTO | 분석 오케스트레이션, ACMG 판정 |
| Counselor | 최종 리포트 작성 |

### 키 발급

1. https://console.anthropic.com/ 접속
2. "API Keys" 메뉴에서 키 생성
3. `.env`에 입력:

```
ANTHROPIC_API_KEY=sk-ant-...
```

---

## OPENAI_API_KEY

- **필수 여부**: **필수** (Pharma 에이전트 사용 시)
- **환경 변수**: `OPENAI_API_KEY`

Codex 기반 에이전트 구동에 필요하다.

### 해당 에이전트

| 에이전트 | 역할 |
|---------|------|
| Pharma (Pharmacogenomicist) | PharmGKB/CPIC 조회, PGx 패턴 매칭 |

### 키 발급

1. https://platform.openai.com/api-keys 접속
2. 새 키 생성
3. `.env`에 입력:

```
OPENAI_API_KEY=sk-...
```

---

## 요약 표

| 키 | 필수 | 서비스 | 비고 |
|----|------|--------|------|
| `NCBI_API_KEY` | 선택 | ClinVar E-utilities | 없으면 3 req/sec 제한 |
| gnomAD | 불필요 | gnomAD API | 키 없음, 배치 시 느릴 수 있음 |
| PharmGKB | 불필요 | PharmGKB API | 키 없음, 0.5s 간격 강제 |
| `ANTHROPIC_API_KEY` | **필수** | CEO, CTO, Counselor | Claude Sonnet |
| `OPENAI_API_KEY` | **필수** | Pharma | Codex |
