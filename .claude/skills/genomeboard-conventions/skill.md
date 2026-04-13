---
name: genomeboard-conventions
description: "BIKO GenomeBoard 프로젝트의 코드 컨벤션, 모듈 구조, 파일 명명 규칙, 테스트 패턴, config 구조, DB 통합 패턴을 정의하는 공유 참조 스킬. 에이전트가 BIKO GenomeBoard 코드를 작성하거나 수정하기 전에 반드시 이 스킬을 참조. 코딩, 구현, 개발, 테스트 작성 시 자동 로드."
---

# BIKO GenomeBoard Conventions

BIKO GenomeBoard 프로젝트의 코드 작성 표준. 모든 개발 에이전트는 코드 작성 전 이 문서를 참조한다.

## 프로젝트 구조

```
scripts/
├── orchestrate.py           — 메인 파이프라인 (run_pipeline, run_batch_pipeline)
├── common/                  — 공유 모델, 설정, 캐시, API 유틸
│   ├── models.py            — Variant, StructuralVariant, AcmgEvidence 등 dataclass
│   ├── config.py            — thread-safe 싱글톤 설정 로더
│   ├── cache.py             — SQLite 응답 캐시 (7일 TTL)
│   ├── api_utils.py         — HTTP retry/backoff
│   ├── gene_knowledge.py    — gene_knowledge.json 로더
│   └── hgvs_utils.py        — HGVSp → CIViC variant 변환
├── intake/                  — 입력 파싱 (VCF, AnnotSV, annotation)
├── clinical/                — 임상 DB 쿼리 (ClinVar API, OncoKB, HPO, OMIM, ClinGen)
├── db/                      — 로컬 DB 빌드/쿼리 (10 build + 9 query 모듈)
├── classification/          — ACMG/AMP 2015 분류 엔진
├── somatic/                 — AMP tiering + TMB
├── korean_pop/              — 한국인 빈도 (KRGDB, gnomAD, 비교)
├── pharma/                  — PGx (Korean PGx, PharmGKB)
├── counselor/               — 리포트 생성 (Jinja2 → HTML → PDF)
└── tools/                   — 빌드 유틸 (gene knowledge, sources)

templates/
├── cancer/report.html       — 암 리포트 (AMP Tier I-IV)
├── rare-disease/report.html — 희귀질환 리포트 (ACMG P/LP/VUS)
└── shared/sv_section.html   — CNV/SV 공유 섹션

tests/
└── test_{module}.py         — 44개 파일, 463개 테스트

config.yaml                  — 전체 설정 (paths, thresholds, api, somatic, cache)
```

## 코드 컨벤션

### Python 스타일
- **Python ≥ 3.10**, type hints 사용
- **Ruff** 린팅 (120 char line length, E/F/W rules)
- dataclass 기반 모델 (Pydantic 아님)
- f-string 선호, .format() 지양

### 명명 규칙
- 모듈: `snake_case.py`
- 클래스: `PascalCase` (Variant, TierResult, ClassificationResult)
- 함수: `snake_case` (query_local_clinvar, amp_assign_tier)
- DB 빌드: `build_{source}_db.py`
- DB 쿼리: `query_{source}.py` 또는 `query_local_{source}.py`
- 테스트: `test_{module}.py`

### 핵심 원칙
- **분류 경로에 LLM 금지** — 모든 분류/tiering은 결정적 Python 로직
- **fallback chain** — local DB → API → skip (graceful degradation)
- **thread-safe** — config, cache, DB 연결에 lock 사용
- **config 기본값** — 새 설정은 항상 기본값 포함, 기존 사용자 영향 없도록

## Config 패턴

```yaml
# config.yaml에 새 항목 추가 시
paths:
  new_source_db: data/db/new_source.sqlite3    # 상대경로 → load_config()에서 절대경로로 변환

# Python에서 사용
from scripts.common.config import get
db_path = get("paths.new_source_db")  # dot-notation 접근
```

## 테스트 패턴

```python
# tests/test_{module}.py
"""Module description."""
import pytest

def test_happy_path():
    from scripts.module import function
    result = function(valid_input)
    assert result.expected_field == expected_value

def test_edge_case():
    ...

def test_error_handling():
    ...

# 통합 테스트 (CI에서 제외)
def test_integration_pipeline(tmp_path):
    from scripts.orchestrate import run_pipeline
    result = run_pipeline(
        vcf_path="data/sample_vcf/demo_variants_grch38_annotated.vcf",
        output_path=str(tmp_path / "report.html"),
        skip_api=True, mode="cancer",
    )
    assert result is not None
```

**테스트 실행:** `cd /Users/JL/Research/gb && python -m pytest tests/ -x -q`

**기존 fixture** (conftest.py): sample_tp53_variant, sample_brca2_variant, sample_cyp2c19_variant, mock_clinvar_pathogenic 등

## DB 통합 패턴

새 데이터 소스를 통합할 때 따르는 표준 패턴. 상세 예시는 `references/db-integration-pattern.md` 참조.

### 빌드 스크립트 (`scripts/db/build_{source}.py`)

```python
"""Build {Source} SQLite database from {format} file."""
import sqlite3

def build_db(source_path: str, db_path: str) -> dict:
    conn = sqlite3.connect(db_path)
    conn.execute("DROP TABLE IF EXISTS {table}")
    conn.execute("CREATE TABLE {table} (...)")
    # batch insert with executemany
    conn.execute("""
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY, value TEXT
        )
    """)
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('source', ?)", (source_url,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (datetime.now().isoformat(),))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('record_count', ?)", (str(count),))
    conn.commit()
    conn.close()
    return {"records": count}

if __name__ == "__main__":
    build_db(source_path, db_path)
```

### 쿼리 모듈 (`scripts/db/query_{source}.py`)

```python
"""Query {Source} local database."""
import sqlite3
from scripts.common.config import get

def query_{source}(gene: str) -> dict | None:
    db_path = get("paths.{source}_db")
    if not db_path or not os.path.exists(db_path):
        return None
    conn = sqlite3.connect(db_path)
    # query logic
    conn.close()
    return result
```

### 등록 체크리스트

1. `config.yaml` — `paths.{source}_db` 추가
2. `scripts/db/build_{source}.py` — 빌드 스크립트
3. `scripts/db/query_{source}.py` — 쿼리 모듈
4. `scripts/orchestrate.py` — 파이프라인에서 쿼리 호출
5. `tests/test_{source}.py` — 단위 + 통합 테스트
6. `scripts/db/version_manager.py` — get_all_db_versions()에 추가 (선택)

## 리포트 템플릿 패턴

- Jinja2 autoescape=True
- `{% set ns_pages = namespace(current=1) %}` — 페이지 번호 관리
- `{% include "shared/section.html" ignore missing %}` — 공유 섹션
- `{{ var|default('N/A') }}` — 안전한 변수 출력
- `@media print { page-break-inside: avoid; }` — 인쇄 최적화

## 모드 분기

```python
# cancer vs rare-disease 모드 분기 패턴
if mode == "cancer":
    # AMP tiering, CIViC enrichment, TMB
    tier = amp_assign_tier(classification, gene, hgvsp, strategy=config_strategy)
elif mode == "rare-disease":
    # HPO scoring, OMIM, ClinGen, inheritance pattern
    hpo_score = calculate_hpo_score(gene, hpo_ids)
```

cancer 전용 기능(CIViC, TMB, AMP tiering)은 반드시 `if mode == "cancer"` 가드를 적용한다.
