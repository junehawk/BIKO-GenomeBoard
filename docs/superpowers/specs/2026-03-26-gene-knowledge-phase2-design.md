# Gene Knowledge Phase 2: Dynamic Builder with Curated Sources

## Goal

기존 18개 LLM-generated gene knowledge를 폐기하고, VCF에서 만나는 모든 유전자를 공인 소스(NCBI Gene, CIViC, GeneReviews, ClinGen, CPIC)에서 동적으로 빌드하는 시스템을 구축한다. Cancer/Rare Disease 모드별로 소스를 분리한다.

## Scope

**In scope (Phase 2a):**
- `build_gene_knowledge.py` 리라이트 — 모드별 소스 분리, 동적 유전자 리스트
- NCBI Gene E-utilities 조회 모듈 (`ncbi_gene.py`)
- GeneReviews PMID 자동 조회 모듈 (`genreviews.py`)
- `generate_pdf.py` mode gate — rare disease에서 CIViC enrichment 비활성화
- 기존 `data/gene_knowledge.json` LLM 데이터 교체
- content_status 체계 확장

**Out of scope (Phase 2b, 향후):**
- Orphanet prevalence XML → SQLite
- GeneReviews XML dump → Management 텍스트 파싱
- OMIM API 통합 (키 확보 시)

---

## Architecture

### 소스 우선순위 (모드별)

#### Cancer Mode

| 필드 | 1순위 (로컬) | 2순위 (API) | content_status |
|------|-------------|------------|----------------|
| finding_summary | CIViC gene description | NCBI Gene summary | curated-civic / curated-ncbi |
| treatment_strategies | CIViC Predictive evidence | — (없으면 빈칸) | curated-civic |
| frequency_prognosis | — (runtime: gnomAD/KRGDB) | — | — |
| references | CIViC evidence PMIDs | PubMed targeted search | — |

#### Rare Disease Mode

| 필드 | 1순위 | 2순위 | content_status |
|------|-------|-------|----------------|
| finding_summary | NCBI Gene summary + ClinGen validity | UniProt function | curated-ncbi |
| treatment_strategies | GeneReviews reference (PMID+URL) | — (없으면 빈칸) | curated-genreviews |
| frequency_prognosis | — (runtime: gnomAD/KRGDB) | — | — |
| references | GeneReviews PMID + ClinGen | PubMed review search | — |

#### PGx Genes (양쪽 공통)

| 필드 | 소스 | content_status |
|------|------|----------------|
| 모든 필드 | CPIC API (기존 유지) | curated-cpic |

### 데이터 흐름

```
[빌드 타임] python scripts/tools/build_gene_knowledge.py --genes <list|vcf|all>

  Input: gene list (CLI arg, VCF 추출, 또는 OncoKB+ClinGen 전체)
    ↓
  For each gene:
    1. CPIC check → PGx gene이면 CPIC 소스 사용 (기존 로직)
    2. CIViC local DB check → description 있으면 cancer 소스로 사용
    3. NCBI Gene API → summary 조회 (fallback 또는 rare disease primary)
    4. GeneReviews PMID → PubMed 자동 검색
    5. ClinGen local DB → gene-disease validity
    6. Compose fields → gene_knowledge entry 생성
    ↓
  Output: data/gene_knowledge.json

[리포트 타임] generate_pdf.py
  Cancer: gene_knowledge + CIViC runtime enrichment
  Rare Disease: gene_knowledge only (CIViC enrichment 비활성화)
```

---

## File Structure

### 신규 파일

| File | Responsibility |
|------|---------------|
| `scripts/tools/sources/__init__.py` | Sources package |
| `scripts/tools/sources/ncbi_gene.py` | NCBI Gene E-utilities 조회 (gene summary, Entrez ID) |
| `scripts/tools/sources/genreviews.py` | GeneReviews PMID 자동 조회 (PubMed E-utilities) |
| `tests/test_ncbi_gene.py` | NCBI Gene 모듈 테스트 |
| `tests/test_genreviews.py` | GeneReviews 모듈 테스트 |
| `tests/test_build_gene_knowledge_v2.py` | 리라이트된 빌더 통합 테스트 |

### 수정 파일

| File | Change |
|------|--------|
| `scripts/tools/build_gene_knowledge.py` | 리라이트 — 모드별 소스 분리, 동적 유전자 리스트 |
| `scripts/counselor/generate_pdf.py` | CIViC enrichment mode gate 추가 |
| `data/gene_knowledge.json` | LLM 데이터 → 빌드 스크립트 생성 데이터로 교체 |

---

## Component Design

### 1. `scripts/tools/sources/ncbi_gene.py`

```python
def fetch_gene_summary(gene_symbol: str, api_key: str = None) -> Optional[Dict]:
    """NCBI Gene E-utilities로 gene summary 조회.
    Returns: {
        "gene": "TP53",
        "entrez_id": "7157",
        "summary": "This gene encodes a tumor protein...",
        "full_name": "tumor protein p53",
        "aliases": ["p53", "LFS1"],
    }
    """
```

- `esearch.fcgi?db=gene&term={SYMBOL}[Gene Name]+AND+Homo+sapiens[Organism]` → Gene ID
- `esummary.fcgi?db=gene&id={ID}&retmode=json` → summary, full name
- Rate limit: 3 req/sec (무키), 10 req/sec (NCBI_API_KEY)
- 에러 시 None 반환 (graceful degradation)

### 2. `scripts/tools/sources/genreviews.py`

```python
def fetch_genreviews_pmid(gene_symbol: str) -> Optional[Dict]:
    """GeneReviews PMID 자동 조회.
    Returns: {
        "gene": "TP53",
        "pmid": "20301371",
        "title": "Li-Fraumeni Syndrome",
        "nbk_id": "NBK1311",
        "url": "https://www.ncbi.nlm.nih.gov/books/NBK1311/",
    }
    """
```

- `esearch.fcgi?db=pubmed&term={GENE}[Title]+AND+GeneReviews[Book]` → PMID
- `esummary.fcgi?db=pubmed&id={PMID}&retmode=json` → title, BookshelfID
- ~800 유전자에 GeneReviews 있음

### 3. `scripts/tools/build_gene_knowledge.py` (리라이트)

```python
def build_knowledge(
    gene_list: List[str],
    output_path: str = "data/gene_knowledge.json",
    civic_db_path: str = None,
    clingen_db_path: str = None,
) -> str:
    """Build gene_knowledge.json from curated sources.

    Source priority per gene:
    1. CPIC (if PGx gene) → content_status: curated-cpic
    2. CIViC local DB (if description exists) → curated-civic
    3. NCBI Gene API → curated-ncbi
    4. Minimal entry → auto-minimal
    """
```

각 유전자에 대해:
1. CPIC PGx 유전자 → 기존 `fetch_cpic_gene()` 사용
2. CIViC local DB에 description 있으면 → finding_summary로 사용
3. NCBI Gene API로 summary 조회 → finding_summary fallback
4. GeneReviews PMID 조회 → references에 추가
5. ClinGen local DB validity → finding_summary에 포함
6. 결과 JSON에 `content_status` 기록

CLI 인터페이스:
```bash
# OncoKB cancer gene list에서 빌드
python scripts/tools/build_gene_knowledge.py --source oncokb

# VCF에서 유전자 추출하여 빌드
python scripts/tools/build_gene_knowledge.py --vcf data/sample_vcf/demo.vcf

# 특정 유전자만 빌드
python scripts/tools/build_gene_knowledge.py --genes TP53,BRCA2,KRAS

# 전체 (OncoKB + ClinGen Definitive/Strong)
python scripts/tools/build_gene_knowledge.py --source all
```

### 4. `generate_pdf.py` Mode Gate

현재 (문제):
```python
# Line 96-98: CIViC가 rare disease에서도 finding_summary 덮어씀
civic_gene = get_gene_summary(gene)
if civic_gene and civic_gene.get("description"):
    v.setdefault("finding_summary", civic_gene["description"][:500])
```

수정:
```python
# CIViC enrichment은 cancer mode에서만
if mode == "cancer":
    civic_gene = get_gene_summary(gene)
    if civic_gene and civic_gene.get("description"):
        v.setdefault("finding_summary", civic_gene["description"][:500])

    # CIViC treatment + references도 cancer only
    civic_variant_name = _hgvsp_to_civic_variant(hgvsp)
    civic_treatment = get_treatment_summary(gene, civic_variant_name)
    if civic_treatment:
        v.setdefault("treatment_strategies", civic_treatment)
    # ... (existing CIViC evidence refs logic)
```

### 5. content_status 체계

| Status | 의미 | 워터마크 |
|--------|------|---------|
| `curated-civic` | CIViC expert description | 없음 |
| `curated-ncbi` | NCBI Gene summary | 없음 |
| `curated-cpic` | CPIC guideline | 없음 |
| `curated-genreviews` | GeneReviews referenced | 없음 |
| `auto-minimal` | 최소 정보 (이름/aliases만) | "Limited data available" |
| `ai-generated-with-references` | LLM 생성 (폐기 대상) | AI watermark |

---

## Testing Strategy

### Unit Tests
- `test_ncbi_gene.py` — NCBI Gene API 조회 (mock), Entrez ID 추출, summary 파싱
- `test_genreviews.py` — GeneReviews PMID 조회 (mock), NBK ID 추출
- `test_build_gene_knowledge_v2.py` — 빌더 통합 (CPIC/CIViC/NCBI 소스 우선순위)

### Integration Tests
- 실제 NCBI API 호출 테스트 (TP53, BRCA2 — 잘 알려진 유전자)
- CIViC local DB → gene_knowledge 생성 확인
- generate_pdf.py mode gate 확인 (rare disease에서 CIViC 비활성화)
- 기존 405 tests 회귀 없음

### Regression
- 기존 리포트 출력과 비교 (cancer, rare disease)
- content_status 기반 워터마크 정상 동작

---

## Migration

1. 기존 `data/gene_knowledge.json` 백업 → `data/gene_knowledge.legacy.json`
2. 빌드 스크립트로 새 `gene_knowledge.json` 생성
3. 기존 `tests/test_build_gene_knowledge.py` → 리라이트된 빌더에 맞게 업데이트
4. `generate_pdf.py`의 `get_gene_info()` fallback chain 유지 (knowledge 없으면 CIViC/빈칸)

---

## Rate Limit 고려

| API | Rate Limit | 300 유전자 예상 시간 |
|-----|-----------|-------------------|
| NCBI Gene | 3 req/sec (무키) | ~2분 (search + summary = 2 calls/gene) |
| PubMed (GeneReviews) | 3 req/sec | ~2분 |
| CPIC | 제한 없음 | ~10초 |
| CIViC | 로컬 DB | 즉시 |
| ClinGen | 로컬 DB | 즉시 |

NCBI_API_KEY 사용 시 10 req/sec → 300 유전자 ~40초.
전체 빌드 예상 시간: **4-5분 (무키) / 1-2분 (키 사용)**
