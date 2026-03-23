# GenomeBoard 개선 로드맵 (Improvement Roadmap)

**작성일:** 2026-03-20
**버전:** v0.1 → v1.0 전환 계획
**작성자:** 리뷰팀 종합 (Scientist × Code Quality × Architecture)
**대상 독자:** 개발팀 전원

---

## 목차

1. [검토 요약](#검토-요약)
2. [우선순위별 이슈 목록](#우선순위별-이슈-목록)
   - [Critical — 즉시 수정 필요](#critical--즉시-수정-필요)
   - [Important — 빠른 시일 내 수정](#important--빠른-시일-내-수정)
   - [Nice-to-have — 최적화 단계](#nice-to-have--최적화-단계)
3. [단계별 개선 계획](#단계별-개선-계획)
   - [v0.2 — 크래시 방지 & 임상 정확도](#v02--크래시-방지--임상-정확도)
   - [v0.3 — 품질 & 보안](#v03--품질--보안)
   - [v1.0 — 아키텍처 최적화](#v10--아키텍처-최적화)
4. [구체적 코드 변경 예시](#구체적-코드-변경-예시)

---

## 검토 요약

세 팀(데이터 정확도 · 코드 품질 · 아키텍처)이 GenomeBoard v0.1을 독립적으로 검토한 결과, **임상 안전성 위협 2건**, **런타임 크래시 3건**, **보안 취약점 1건**, **아키텍처 비효율 다수**를 발견했다. 현재 상태로는 임상 환경 배포가 불가하다.

| 카테고리 | Critical | Important | Nice-to-have |
|---|---|---|---|
| 데이터 정확도 | 2 | 3 | 3 |
| 코드 품질 | 3 | 12 | 4 |
| 아키텍처 | 0 | 4 | 5 |
| **합계** | **5** | **19** | **12** |

---

## 우선순위별 이슈 목록

### Critical — 즉시 수정 필요

> 임상적 위험 또는 런타임 크래시를 유발하는 이슈. v0.2 배포 전 반드시 해결.

---

#### C-1. PGx 변이에 ACMG BA1 규칙 오적용 (임상 위험)
- **출처:** Scientist 리뷰
- **노력(Effort):** Medium
- **위험도:** 임상 위험 (Clinical Safety Risk)

CYP2C19\*2 (rs4244285)와 HLA-B\*5701 (rs2395029) 두 변이가 population frequency BA1(>5%) 기준에 해당한다는 이유로 **"Benign"으로 분류**되고 있다. 그러나 ACMG/AMP 가이드라인은 pharmacogenomic (PGx) 변이에는 적용하지 않도록 명시하고 있다. HLA-B\*5701은 아바카비르(abacavir) 과민반응의 FDA Black Box Warning 바이오마커이며 "Benign" 출력은 **임상적으로 매우 위험하다.**

**수정 방향:**
```python
# acmg_engine.py 또는 classify_variant() 내
PGX_GENES = {"CYP2C19", "CYP2D6", "CYP3A5", "TPMT", "DPYD", "HLA-B", "HLA-A"}

def classify_variant(gene: str, ...):
    if gene in PGX_GENES:
        return classify_pgx_variant(gene, ...)  # ACMG 우회, CPIC 기준 적용
    # 기존 ACMG 로직
```

---

#### C-2. CYP2C19\*2 표현형 오분류 — Poor Metabolizer 출력
- **출처:** Scientist 리뷰
- **노력:** Small
- **위험도:** 임상 오류 (Clinical Error)

CYP2C19\*1/\*2 이형접합(heterozygous) 보인자를 **"Poor Metabolizer (PM)"으로 분류** 중이다. CPIC 가이드라인에 따르면 \*1/\*2는 **"Intermediate Metabolizer (IM)"**이다. PM은 \*2/\*2 동형접합체에만 해당한다. 클로피도그렐(clopidogrel) 용량 조정 지침이 달라지므로 임상 오류다.

**수정 방향:**
```python
# korean_pgx.py — diplotype_to_phenotype()
CYPC19_TABLE = {
    "*1/*1": "Normal Metabolizer",
    "*1/*2": "Intermediate Metabolizer",  # 현재 잘못됨: "Poor Metabolizer"
    "*2/*2": "Poor Metabolizer",
    # ...
}
```

---

#### C-3. FileNotFoundError — query_krgdb.py 크래시
- **출처:** Code Quality 리뷰
- **노력:** Small
- **위험도:** 런타임 크래시

KRGDB TSV 파일 경로가 하드코딩된 절대경로이며, 파일이 없을 때 예외 처리 없이 크래시된다.

**수정 방향:**
```python
# query_krgdb.py
from pathlib import Path

KRGDB_PATH = Path(__file__).parent.parent / "data" / "krgdb_freq.tsv"

def load_krgdb():
    if not KRGDB_PATH.exists():
        raise FileNotFoundError(
            f"KRGDB 데이터 파일을 찾을 수 없습니다: {KRGDB_PATH}\n"
            "data/krgdb_freq.tsv 파일을 프로젝트 루트에 배치하세요."
        )
    return pd.read_csv(KRGDB_PATH, sep="\t")
```

---

#### C-4. FileNotFoundError — korean_pgx.py 크래시
- **출처:** Code Quality 리뷰
- **노력:** Small
- **위험도:** 런타임 크래시

`korean_pgx_table.json` 로드 시 파일 미존재 예외 처리 없음.

**수정 방향:** C-3와 동일한 패턴 — `Path(__file__).parent.parent` 기준 상대 경로 + 명확한 오류 메시지.

---

#### C-5. WeasyPrint ImportError — generate_pdf.py 크래시
- **출처:** Code Quality 리뷰
- **노력:** Small
- **위험도:** 런타임 크래시

`import weasyprint` 시 시스템 라이브러리(cairo, pango) 미설치 환경에서 ImportError 발생. 아무 fallback 없이 전체 파이프라인이 종료된다.

**수정 방향:**
```python
# generate_pdf.py
try:
    import weasyprint
    WEASYPRINT_AVAILABLE = True
except ImportError:
    WEASYPRINT_AVAILABLE = False
    import logging
    logging.warning("WeasyPrint 미설치. PDF 출력이 비활성화됩니다. HTML 리포트만 생성됩니다.")

def generate_pdf(html_path, output_path):
    if not WEASYPRINT_AVAILABLE:
        raise RuntimeError("PDF 생성 불가: WeasyPrint가 설치되지 않았습니다.")
    # ...
```

---

### Important — 빠른 시일 내 수정

> 보안, 데이터 품질, 안정성에 영향을 미치는 이슈. v0.3 목표.

---

#### I-1. Thread-unsafe 캐시 (보안 · 안정성)
- **출처:** Code Quality
- **노력:** Small
- **위험도:** Race condition, 잘못된 데이터 반환

글로벌 dict를 캐시로 사용하는 코드가 멀티스레드 환경에서 race condition을 유발한다.

**수정 방향:** `functools.lru_cache` 또는 `threading.Lock` 적용.
```python
import threading
_cache_lock = threading.Lock()
_acmg_rules_cache = {}

def get_acmg_rules():
    with _cache_lock:
        if "rules" not in _acmg_rules_cache:
            _acmg_rules_cache["rules"] = _load_rules_from_disk()
        return _acmg_rules_cache["rules"]
```

---

#### I-2. ACMG rules 매 호출마다 디스크에서 로드
- **출처:** Code Quality
- **노력:** Small

`acmg_engine.py`가 변이 분류 호출마다 JSON 파일을 읽는다. 10개 변이 × N회 = 과도한 I/O.

**수정 방향:** 모듈 수준 singleton 또는 `@lru_cache(maxsize=1)` 적용.

---

#### I-3. Jinja2 XSS 취약점
- **출처:** Code Quality
- **노력:** Small
- **위험도:** 보안 (Security)

리포트 생성 시 사용자 입력 데이터(유전자명, 변이명 등)가 Jinja2 템플릿에 `|safe` 필터 없이 또는 auto-escape 비활성화 상태로 삽입될 위험.

**수정 방향:**
```python
# generate_report.py
from jinja2 import Environment, FileSystemLoader, select_autoescape

env = Environment(
    loader=FileSystemLoader("templates"),
    autoescape=select_autoescape(["html", "xml"])  # XSS 방지 필수
)
```

---

#### I-4. PTPN11 ClinVar 충돌 플래그 누락
- **출처:** Scientist 리뷰
- **노력:** Small

분류 엔진은 PTPN11을 "Likely Pathogenic"으로 출력하지만, ClinVar는 "Pathogenic"이다. `conflict` 플래그가 `false`로 표시되어 임상의가 불일치를 인지하지 못한다.

**수정 방향:** ClinVar 조회 결과와 내부 분류 비교 로직에서 등급 차이(LP vs P) 감지 시 `conflict: true` 설정.

---

#### I-5. APOE "Benign" vs ClinVar "Risk Factor" 불일치
- **출처:** Scientist 리뷰
- **노력:** Small

APOE ε4 (rs429358)에 BA1이 적용되어 "Benign" 출력. ClinVar는 알츠하이머 위험인자("Risk Factor")로 분류. PGx와 동일하게 risk-factor 변이는 별도 처리 로직 필요.

**수정 방향:** `RISK_FACTOR_GENES` 목록 도입, BA1 적용 제외.

---

#### I-5b. compare_freq.py: 빈도 임계값 gap — ACMG 코드 미할당 구간
- **출처:** Scientist 리뷰 (Task #1)
- **노력:** Small
- **위험도:** 분류 누락 (Silent miss)

`compare_freq.py`의 임계값 구조에 빈틈이 있다.

```
BA1:            freq > 0.05    → Benign (stand-alone)
BS1:            freq >= 0.01   → Benign (strong)
PM2_Supporting: freq <= 0.001  → Pathogenic supporting

⚠️ 0.001 < freq < 0.01 구간: ACMG 코드 미할당
```

이 범위(0.1%~1%)에 해당하는 변이는 어떤 빈도 코드도 부여받지 못한다.
ACMG/AMP 2015 기준에서 이 구간은 `PM2` (moderate)에 해당한다.

**수정 방향:**
```python
# compare_freq.py — threshold 사이 gap 처리 추가
elif max_freq <= PM2_THRESHOLD:
    acmg_codes.append("PM2_Supporting")
elif max_freq < BS1_THRESHOLD:   # 새 구간: 0.001 < freq < 0.01
    acmg_codes.append("PM2")     # moderate 수준 희귀 변이
```

---

#### I-5c. API 불가 시 리포트 미표기 — 데이터 출처 불명확
- **출처:** Scientist 리뷰 (Task #1)
- **노력:** Small
- **위험도:** 데이터 신뢰성 (Data Reliability)

현 실행에서 ClinVar/gnomAD API가 불가한 상태에서 ACMG 코드(PP5, PM2_Supporting)가 수동 생성되었으나, 리포트의 DB 버전 표기(`ClinVar: 2026-03`, `gnomAD: v4.1`)는 정상 조회된 것처럼 표시된다. 임상의는 이 값들이 실제 API 조회 결과가 아님을 알 수 없다.

**수정 방향:**
- 파이프라인에서 API 실패 시 `api_status` 플래그를 결과에 기록
- 리포트 DB 버전 섹션에 조회 실패 여부 표시
```json
// consilium_data.json에 추가
"api_status": {
  "ClinVar": "unavailable — codes manually assigned",
  "gnomAD": "unavailable — PM2 thresholds unverified"
}
```

---

#### I-6. Indiscriminate retry 로직
- **출처:** Code Quality
- **노력:** Small

HTTP 오류 전체에 retry를 적용하여 4xx(잘못된 요청) 오류도 재시도. 불필요한 API 호출 유발.

**수정 방향:**
```python
# 5xx만 재시도, 4xx는 즉시 실패
retryable_codes = {500, 502, 503, 504}
if response.status_code not in retryable_codes:
    raise HTTPError(response)
```

---

#### I-7. parse_vcf 오류 처리 없음
- **출처:** Code Quality
- **노력:** Small

VCF 파싱 중 형식 오류, 누락 컬럼, 인코딩 오류에 대한 처리 없음. 잘못된 VCF 파일 입력 시 비명시적 오류 발생.

---

#### I-8. GraphQL 에러 응답 미처리
- **출처:** Code Quality
- **노력:** Small

ClinVar GraphQL 쿼리 응답의 `errors` 필드를 검사하지 않아 부분적 실패가 무시된다.

---

#### I-9. Paperclip 에이전트 침묵 실패
- **출처:** Code Quality · Architecture
- **노력:** Medium

Paperclip 스킬 호출 실패 시 아무런 오류 메시지 없이 빈 결과 반환. 임상의는 데이터 누락을 인지하지 못한다.

**수정 방향:** 에이전트 응답 검증 + 실패 시 명시적 오류 전파.

---

#### I-10. Bash(*) 퍼미션 과도하게 넓음
- **출처:** Architecture
- **노력:** Small

`.claude/settings.json`의 `allowedTools`에 `Bash(*)` 패턴이 포함되어 있어 모든 셸 명령 실행이 허용된다.

**수정 방향:**
```json
// .claude/settings.json
"allowedTools": [
  "Bash(python src/**)",
  "Bash(python scripts/**)",
  "Bash(pip install*)",
  "Read", "Write", "Edit", "Glob", "Grep"
]
```

---

#### I-11. 에이전트 중복 등록 (Duplicate agent import)
- **출처:** Architecture
- **노력:** Small

동일 에이전트를 여러 번 import 시 Paperclip에 중복 등록된다.

**수정 방향:** 등록 전 기존 에이전트 존재 여부 확인하는 idempotent 전략 도입.

---

#### I-12. Hardcoded CWD — 이식성 파괴
- **출처:** Architecture
- **노력:** Small

`cwd: "/Users/JL/Research/gb"` 형태의 절대경로 하드코딩. 다른 개발자 환경에서 즉시 실패.

**수정 방향:**
```typescript
// agent 설정
cwd: process.env.GB_PROJECT_ROOT ?? path.resolve(__dirname, "../..")
```

---

#### I-13. PM2_Supporting 기여 점수 오집계
- **출처:** Scientist 리뷰
- **노력:** Small

`acmg_engine.py`에서 `PM2_Supporting`을 `PM` (Pathogenic Moderate) 수준으로 집계하는 잠재적 버그. 현재 데모 케이스에서는 분류 결과에 영향 없지만, 경계 케이스에서 오분류를 유발할 수 있다.

---

#### I-14. ARCHITECTURE.md가 삭제된 파일 참조
- **출처:** Architecture
- **노력:** Small

`ARCHITECTURE.md`가 `company.json`, `org-chart.json` 등 삭제된 파일을 참조. 신규 개발자 온보딩 혼란 유발.

---

### Nice-to-have — 최적화 단계

> 성능, 비용, 유지보수성 개선. v1.0 목표.

---

#### N-1. CEO→CTO→Counselor 최악 경로 10분 Heartbeat 지연
- **출처:** Architecture
- **노력:** Medium

다중 에이전트 체인에서 응답 최악 경우 10분. 비동기 팬아웃(fan-out) 패턴으로 개선 가능.

---

#### N-2. CTO 에이전트 Opus 모델 — 스크립트 실행에 과도한 비용
- **출처:** Architecture
- **노력:** Small

Python 스크립트 실행만 담당하는 CTO 에이전트에 `claude-opus-4-6` 사용. `claude-sonnet-4-6`으로 충분하며 비용이 대폭 절감된다.

**수정 방향:**
```json
// agents/cto.json
{ "model": "claude-sonnet-4-6" }
```

---

#### N-3. LLM 불필요한 에이전트 4개 (7개 중)
- **출처:** Architecture
- **노력:** Large

현재 7개 에이전트 중 4개(KRGDB 조회, VCF 파싱, PDF 생성, ClinVar 쿼리)는 결정론적 로직으로 LLM 없이 구현 가능. **Option B: 3개 LLM 에이전트로 통합** 권장.

| 현재 에이전트 | 권장 처리 방식 |
|---|---|
| krgdb-agent | Python 함수 (Direct call) |
| vcf-parser-agent | Python 함수 (Direct call) |
| pdf-generator-agent | Python 함수 (Direct call) |
| clinvar-agent | Python 함수 (Direct call) |
| counselor-agent | LLM 유지 (자연어 해석) |
| cto-agent | LLM 유지 (오케스트레이션) |
| ceo-agent | LLM 유지 (임상 판단) |

---

#### N-4. CPIC 버전 2022 — 최신화 필요
- **출처:** Scientist 리뷰
- **노력:** Small

`korean_pgx_table.json`의 CPIC 데이터가 2022년판. 2024년 업데이트 반영 필요 (특히 CYP2C19 + SSRIs 지침).

---

#### N-5. HLA-B 한국인 빈도 6배 강화 — 리포트 미표시
- **출처:** Scientist 리뷰
- **노력:** Small

HLA-B\*5701의 한국인 빈도가 글로벌 대비 6배 높지만 리포트에 표시되지 않는다. 한국인 특이 enrichment/depletion 정보는 리포트의 핵심 가치다.

---

#### N-6. Agent 프롬프트의 "CRITICAL" 강조 의존 — 취약한 패턴
- **출처:** Architecture
- **노력:** Medium

에이전트 프롬프트에서 `/paperclip` 스킬을 호출하게 만들기 위해 "CRITICAL: YOU MUST USE /paperclip" 형태의 대문자 강조에 의존. LLM 출력 비결정성으로 인해 취약하다.

**수정 방향:** `wrapWithPreamble()` Worker Preamble Protocol 적용으로 스킬 호출을 구조화.

---

#### N-7. PharmGKB 응답 스키마 불일치
- **출처:** Code Quality
- **노력:** Small

PharmGKB API 응답 파싱 시 스키마 버전 변화에 취약. 런타임 KeyError 위험.

---

#### N-8. 테스트 커버리지 부족 & 느슨한 어설션
- **출처:** Code Quality
- **노력:** Medium

핵심 분류 로직(acmg_engine, korean_pgx) 단위 테스트 부재. 기존 테스트도 반환값 타입만 확인하는 느슨한 어설션 사용.

---

#### N-9. CLI 입력 오류 처리 일관성 없음
- **출처:** Code Quality
- **노력:** Small

일부 스크립트는 잘못된 인자에서 Python 트레이스백 전체를 출력, 다른 스크립트는 조용히 종료. `argparse` + 일관된 오류 메시지 표준화 필요.

---

#### N-10. models.py 내 Dead Conditional
- **출처:** Code Quality
- **노력:** Small

도달 불가능한 조건문이 코드에 존재. 로직 이해를 어렵게 하고 테스트 커버리지를 낮춘다.

---

## 단계별 개선 계획

### v0.2 — 크래시 방지 & 임상 정확도
**목표:** 임상 환경에서 안전하게 실행 가능한 최소 수준
**기간 예상:** 1주
**담당:** 전체 팀

| 이슈 | ID | Effort |
|---|---|---|
| PGx 변이 ACMG BA1 오적용 제거 | C-1 | Medium |
| CYP2C19 표현형 IM/PM 수정 | C-2 | Small |
| query_krgdb.py FileNotFoundError 처리 | C-3 | Small |
| korean_pgx.py FileNotFoundError 처리 | C-4 | Small |
| generate_pdf.py WeasyPrint fallback | C-5 | Small |
| PTPN11 conflict 플래그 수정 | I-4 | Small |
| APOE Risk Factor 별도 처리 | I-5 | Small |
| compare_freq.py threshold gap 수정 | I-5b | Small |
| API 불가 시 리포트 경고 표시 | I-5c | Small |

**완료 기준:**
- `python src/acmg_engine.py` — PGx 변이에 BA1 미적용 확인
- `python src/korean_pgx.py` — CYP2C19\*1/\*2 → "Intermediate Metabolizer" 출력 확인
- 데이터 파일 없는 환경에서 명확한 오류 메시지 출력 확인
- Consilium 리포트 재생성 후 HLA-B "Benign" 출력 없음 확인

---

### v0.3 — 품질 & 보안
**목표:** 프로덕션 코드 품질 기준 달성
**기간 예상:** 2주
**담당:** 백엔드 개발자

| 이슈 | ID | Effort |
|---|---|---|
| Thread-safe 캐시 도입 | I-1 | Small |
| ACMG rules singleton 캐싱 | I-2 | Small |
| Jinja2 auto-escape 활성화 | I-3 | Small |
| Retry 로직 4xx 제외 | I-6 | Small |
| parse_vcf 오류 처리 추가 | I-7 | Small |
| GraphQL errors 필드 검사 | I-8 | Small |
| Paperclip 실패 명시적 전파 | I-9 | Medium |
| Bash(*) 퍼미션 범위 축소 | I-10 | Small |
| 에이전트 중복 등록 방지 | I-11 | Small |
| Hardcoded CWD 제거 | I-12 | Small |
| PM2_Supporting 점수 버그 수정 | I-13 | Small |
| ARCHITECTURE.md 참조 정정 | I-14 | Small |
| CLI 오류 처리 표준화 | N-9 | Small |
| Dead conditional 제거 | N-10 | Small |
| PharmGKB 스키마 방어적 파싱 | N-7 | Small |

**완료 기준:**
- `pytest` 전체 통과
- Jinja2 autoescape 설정 확인
- 멀티스레드 환경에서 캐시 race condition 없음

---

### v1.0 — 아키텍처 최적화
**목표:** 비용 효율적이고 확장 가능한 시스템
**기간 예상:** 3-4주
**담당:** 아키텍처 팀

| 이슈 | ID | Effort |
|---|---|---|
| LLM 불필요 에이전트 4개 → Python 함수 전환 | N-3 | Large |
| Heartbeat 비동기 팬아웃 최적화 | N-1 | Medium |
| CTO 에이전트 모델 Opus→Sonnet 다운그레이드 | N-2 | Small |
| CPIC 2024 업데이트 | N-4 | Small |
| HLA-B 한국인 enrichment 리포트 표시 | N-5 | Small |
| Worker Preamble Protocol 적용 | N-6 | Medium |
| 핵심 분류 로직 단위 테스트 작성 | N-8 | Medium |

**완료 기준:**
- 에이전트 수 7개 → 3개로 감소
- 전체 파이프라인 응답 시간 < 3분 (최악 경로)
- CTO 모델 변경 후 Opus 사용 비용 제거 확인
- CPIC 2024 기반 표현형 분류 정확도 검증

---

## 구체적 코드 변경 예시

### 예시 1: PGx 변이 ACMG 우회 (C-1)

```python
# src/acmg_engine.py

# 기존 (잘못됨)
def classify_variant(gene, variant_id, freq, ...):
    if freq > 0.05:
        return Classification(code="BA1", pathogenicity="Benign")

# 수정
PGX_GENES = frozenset({
    "CYP2C19", "CYP2D6", "CYP3A5", "CYP1A2",
    "TPMT", "DPYD", "UGT1A1", "SLCO1B1",
    "HLA-B", "HLA-A", "NUDT15", "G6PD"
})

RISK_FACTOR_GENES = frozenset({"APOE", "MTHFR"})

def classify_variant(gene, variant_id, freq, clinvar_sig, ...):
    if gene in PGX_GENES:
        # ACMG BA1/BS1 적용 금지 — CPIC 기준으로 분리 처리
        return classify_pgx_variant(gene, variant_id, ...)
    if gene in RISK_FACTOR_GENES and clinvar_sig == "risk factor":
        return Classification(
            code="RF",
            pathogenicity="Risk Factor",
            note="ACMG BA1 미적용 — 위험인자 변이"
        )
    # 기존 ACMG 로직
    if freq > 0.05:
        return Classification(code="BA1", pathogenicity="Benign")
    ...
```

---

### 예시 2: CYP2C19 표현형 수정 (C-2)

```python
# src/korean_pgx.py

# 기존 (잘못됨)
CYP2C19_PHENOTYPE = {
    "*1/*1": "Normal Metabolizer",
    "*1/*2": "Poor Metabolizer",   # 오류
    "*2/*2": "Poor Metabolizer",
}

# 수정 (CPIC 2024 기준)
CYP2C19_PHENOTYPE = {
    "*1/*1": "Normal Metabolizer",
    "*1/*2": "Intermediate Metabolizer",   # 수정
    "*1/*3": "Intermediate Metabolizer",
    "*2/*2": "Poor Metabolizer",
    "*2/*3": "Poor Metabolizer",
    "*17/*17": "Ultrarapid Metabolizer",
    "*1/*17": "Rapid Metabolizer",
}
```

---

### 예시 3: 경로 하드코딩 제거 (I-12)

```typescript
// src/agents/base-agent.ts

// 기존 (잘못됨)
const CWD = "/Users/JL/Research/gb";

// 수정
import path from "path";
const CWD = process.env.GB_PROJECT_ROOT
  ?? path.resolve(__dirname, "../../");

// 사용
const agent = new PaperclipAgent({
  cwd: CWD,
  // ...
});
```

---

### 예시 4: Bash 퍼미션 범위 축소 (I-10)

```json
// .claude/settings.json (수정 전)
{
  "allowedTools": ["Bash(*)"]
}

// .claude/settings.json (수정 후)
{
  "allowedTools": [
    "Bash(python src/**)",
    "Bash(python scripts/**)",
    "Bash(python -m pytest**)",
    "Bash(pip install *)",
    "Bash(cat data/**)",
    "Read", "Write", "Edit", "Glob", "Grep",
    "WebFetch"
  ]
}
```

---

### 예시 5: WeasyPrint Graceful Degradation (C-5)

```python
# src/generate_pdf.py

try:
    from weasyprint import HTML as WeasyHTML
    _PDF_AVAILABLE = True
except ImportError:
    _PDF_AVAILABLE = False
    import logging
    logging.warning(
        "WeasyPrint를 가져올 수 없습니다. "
        "PDF 생성이 비활성화됩니다. "
        "설치 방법: pip install weasyprint (+ 시스템 라이브러리 필요)"
    )

def generate_pdf(html_path: str, output_path: str) -> str:
    """HTML을 PDF로 변환합니다. WeasyPrint 미설치 시 HTML 경로를 반환합니다."""
    if not _PDF_AVAILABLE:
        logging.info(f"PDF 건너뜀. HTML 리포트 사용: {html_path}")
        return html_path
    WeasyHTML(filename=html_path).write_pdf(output_path)
    return output_path
```

---

## 최종 요약

| 단계 | 핵심 목표 | 이슈 수 | 예상 기간 |
|---|---|---|---|
| **v0.2** | 임상 안전성 + 크래시 방지 | 9개 (Critical 5 + 중요 4) | 1주 |
| **v0.3** | 코드 품질 + 보안 | 15개 (Important) | 2주 |
| **v1.0** | 아키텍처 최적화 + 비용 절감 | 8개 (Nice-to-have) | 3-4주 |

**v0.2는 임상 환경 배포의 최소 조건이다.** HLA-B "Benign" 오분류와 CYP2C19 IM/PM 오분류는 실제 환자 치료 결정에 영향을 줄 수 있으므로, 이 두 이슈는 다른 모든 작업보다 우선시되어야 한다.

---

#### N-9. Gene Knowledge Base — LLM 생성 텍스트를 공인 출처 기반으로 교체 (Critical for Production)
- **출처:** 개발 과정 리뷰
- **노력:** Large
- **위험도:** 데이터 신뢰성 (Data Reliability)

현재 `data/gene_knowledge.json`의 모든 임상 텍스트(치료 전략, 빈도/예후, 소견 요약)는 **LLM(Claude Opus)이 학습 데이터를 기반으로 작성**한 것이며, 실시간 API 조회나 peer-reviewed 출처 인용이 없다. FoundationOne CDx 보고서는 모든 문장에 PubMed 인용번호(PMID)를 포함한다.

**수정 방향:**
1. `gene_knowledge.json`에 `references` 필드 추가: `["PMID:12345", "ClinVar:VCV000012375"]`
2. GeneReviews (NCBI) API/크롤링으로 유전자별 공인 임상 요약 자동 수집
3. CPIC API → PGx 가이드라인 원문 + PMID 인용 자동 연동
4. ClinVar Variation API → 변이별 submitter evidence 수집
5. 리포트 템플릿에 각 문장별 인용 번호 표시 (FoundationOne 스타일)
6. 모든 AI 생성 텍스트에 "AI-generated, not peer-reviewed" 워터마크 표시 (공인 출처 확보 전까지)

**단계별 접근:**
- Phase 1: `references` 필드 추가 + AI-generated 워터마크 (Small)
- Phase 2: CPIC/PharmGKB API 자동 연동 (Medium)
- Phase 3: GeneReviews + ClinVar + PubMed 통합 파이프라인 (Large)
