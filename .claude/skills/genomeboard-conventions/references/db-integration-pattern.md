# DB Integration Pattern — 상세 예시

BIKO GenomeBoard에 새 데이터 소스를 통합할 때 따르는 표준 패턴의 상세 예시.

## 목차

1. [빌드 스크립트 상세](#빌드-스크립트-상세)
2. [쿼리 모듈 상세](#쿼리-모듈-상세)
3. [config.yaml 등록](#configyaml-등록)
4. [파이프라인 통합](#파이프라인-통합)
5. [테스트 패턴](#테스트-패턴)
6. [실제 예시: CIViC DB](#실제-예시-civic-db)
7. [소스 형식별 파싱 팁](#소스-형식별-파싱-팁)

---

## 빌드 스크립트 상세

```python
"""Build {SourceName} SQLite database from {format} file.

Usage:
    python scripts/db/build_{source}_db.py [source_path] [db_path]

Source: {URL where data can be downloaded}
Format: {TSV/XML/VCF/GFF, encoding, delimiter}
"""
import sqlite3
import os
import sys
from datetime import datetime


SOURCE_URL = "https://example.com/data.tsv"


def build_db(source_path: str, db_path: str) -> dict:
    """Build SQLite database from source file.

    Returns dict with build statistics.
    """
    if not os.path.exists(source_path):
        raise FileNotFoundError(f"Source file not found: {source_path}")

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # 테이블 생성 (멱등)
    cursor.execute("DROP TABLE IF EXISTS main_table")
    cursor.execute("""
        CREATE TABLE main_table (
            gene TEXT NOT NULL,
            field1 TEXT,
            field2 REAL,
            PRIMARY KEY (gene, field1)
        )
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_gene ON main_table(gene)")

    # 데이터 임포트 (batch insert)
    count = 0
    batch = []
    with open(source_path, encoding="utf-8") as f:
        header = f.readline()  # skip header
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue
            batch.append((fields[0], fields[1], float(fields[2])))
            if len(batch) >= 10000:
                cursor.executemany(
                    "INSERT OR REPLACE INTO main_table VALUES (?, ?, ?)", batch
                )
                count += len(batch)
                batch = []
    if batch:
        cursor.executemany(
            "INSERT OR REPLACE INTO main_table VALUES (?, ?, ?)", batch
        )
        count += len(batch)

    # 메타데이터 기록 (필수)
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY, value TEXT
        )
    """)
    cursor.execute("INSERT OR REPLACE INTO metadata VALUES ('source_url', ?)", (SOURCE_URL,))
    cursor.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (datetime.utcnow().isoformat(),))
    cursor.execute("INSERT OR REPLACE INTO metadata VALUES ('record_count', ?)", (str(count),))

    conn.commit()
    conn.close()
    return {"records": count, "db_path": db_path}


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "data/db/source_data.tsv"
    db = sys.argv[2] if len(sys.argv) > 2 else "data/db/source.sqlite3"
    stats = build_db(src, db)
    print(f"Built {db}: {stats['records']} records")
```

### 빌드 스크립트 필수 요소

1. **멱등성** — DROP TABLE IF EXISTS + CREATE TABLE
2. **batch INSERT** — 10,000건 단위 executemany (대용량 성능)
3. **인덱스** — 주요 쿼리 필드에 CREATE INDEX
4. **metadata 테이블** — source_url, build_date, record_count 필수
5. **CLI 인터페이스** — `if __name__ == "__main__"` 블록

---

## 쿼리 모듈 상세

```python
"""Query {SourceName} local database."""
import os
import sqlite3
from scripts.common.config import get


def get_{source}_for_gene(gene: str) -> dict | None:
    """Get {source} data for a gene.

    Returns None if DB not available or gene not found.
    """
    db_path = get("paths.{source}_db")
    if not db_path or not os.path.exists(db_path):
        return None

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    cursor.execute("SELECT * FROM main_table WHERE gene = ?", (gene,))
    rows = cursor.fetchall()
    conn.close()

    if not rows:
        return None

    return {
        "gene": gene,
        "entries": [dict(row) for row in rows],
    }
```

### 쿼리 모듈 필수 요소

1. **DB 부재 시 None 반환** — FileNotFoundError/ImportError 금지
2. **config에서 경로 로드** — `get("paths.{source}_db")`
3. **sqlite3.Row** — dict-like 접근 가능
4. **결과 없으면 None** — 빈 list가 아닌 None
5. **conn.close()** — 연결 반드시 닫기

---

## config.yaml 등록

```yaml
paths:
  # 기존 경로들...
  new_source_db: data/db/new_source.sqlite3
```

config.py의 `load_config()`가 상대경로를 자동으로 절대경로로 변환한다.

---

## 파이프라인 통합

orchestrate.py의 `_query_variant_databases()` 또는 `_build_variant_records()`에 추가:

```python
# _query_variant_databases() 내부 — 변이별 쿼리
from scripts.db.query_new_source import get_new_source_for_gene
new_source_data = get_new_source_for_gene(variant.gene)
if new_source_data:
    result["new_source"] = new_source_data
```

---

## 테스트 패턴

```python
"""Tests for {source} database integration."""
import pytest
import os


def test_build_db(tmp_path):
    """빌드 스크립트가 올바른 DB를 생성하는지."""
    from scripts.db.build_new_source_db import build_db

    # 테스트용 소스 파일 생성
    source = tmp_path / "test_data.tsv"
    source.write_text("gene\tfield1\tfield2\nBRCA1\tval1\t0.5\nTP53\tval2\t0.8\n")

    db_path = str(tmp_path / "test.sqlite3")
    stats = build_db(str(source), db_path)

    assert stats["records"] == 2
    assert os.path.exists(db_path)


def test_query_found(tmp_path, monkeypatch):
    """존재하는 유전자 쿼리."""
    # DB 빌드 후 쿼리
    ...


def test_query_not_found(tmp_path, monkeypatch):
    """존재하지 않는 유전자 쿼리 → None."""
    ...


def test_query_no_db(monkeypatch):
    """DB 파일 없을 때 → None (에러 아님)."""
    monkeypatch.setattr("scripts.common.config.get", lambda k: "/nonexistent/path.sqlite3")
    from scripts.db.query_new_source import get_new_source_for_gene
    assert get_new_source_for_gene("BRCA1") is None
```

---

## 실제 예시: CIViC DB

가장 복잡한 기존 DB 통합 사례. 자동 다운로드 + 다중 테이블 + 다양한 쿼리 함수.

**빌드:** `build_civic_db.py` — TSV 자동 다운로드, genes/variants/evidence/hotspots 4개 테이블
**쿼리:** `query_civic.py` — get_gene_summary(), get_variant_evidence(), is_hotspot(), get_predictive_evidence_for_tier()

새 DB 통합 시 CIViC 수준의 복잡도가 필요한 경우는 드물다. 대부분 단일 테이블 + 단일 쿼리 함수 패턴(KOVA v7 Korean-frequency 모듈이 대표 예시)으로 충분하다.

---

## 소스 형식별 파싱 팁

| 형식 | 인코딩 | 구분자 | 주의사항 |
|------|--------|--------|---------|
| TSV (일반) | utf-8 | `\t` | 헤더 라인 skip |
| TSV (NCBI) | utf-8 | `\t` | `#` 접두사 헤더 가능 |
| TSV (GeneReviews) | **latin-1** | **`\|` (pipe)** | 인코딩+구분자 모두 비표준 |
| XML (Orphanet) | utf-8 | — | ElementTree 파싱, 깊은 중첩 |
| VCF.bgz | — | `\t` | pysam.TabixFile 사용, SQLite 변환 불필요 |
| GFF3 | utf-8 | `\t` | attributes 필드는 `;`+`=` 파싱 |
| CSV | utf-8 | `,` | csv.DictReader 사용, None 값 주의 `.lower()` 등 |
