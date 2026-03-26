# Gene Knowledge Phase 2: Dynamic Builder Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** LLM-generated gene knowledge를 폐기하고, NCBI Gene/CIViC/GeneReviews/ClinGen/CPIC에서 동적으로 빌드하는 시스템을 구축한다.

**Architecture:** `scripts/tools/sources/` 패키지에 소스별 조회 모듈을 분리하고, `build_gene_knowledge.py`를 리라이트하여 소스 우선순위 체인으로 유전자별 knowledge를 조합한다. `generate_pdf.py`는 mode gate로 cancer/rare disease 소스 분리.

**Tech Stack:** Python 3, NCBI E-utilities, SQLite3 (CIViC/ClinGen local), pytest

---

## File Structure

### 신규 파일
| File | Responsibility |
|------|---------------|
| `scripts/tools/sources/__init__.py` | Sources package init |
| `scripts/tools/sources/ncbi_gene.py` | NCBI Gene E-utilities 조회 (summary, Entrez ID) |
| `scripts/tools/sources/genreviews.py` | GeneReviews PMID 자동 조회 (PubMed E-utilities) |
| `tests/test_ncbi_gene.py` | NCBI Gene 모듈 테스트 |
| `tests/test_genreviews.py` | GeneReviews 모듈 테스트 |
| `tests/test_build_gene_knowledge_v2.py` | 빌더 리라이트 통합 테스트 |

### 수정 파일
| File | Change |
|------|--------|
| `scripts/tools/build_gene_knowledge.py` | 리라이트 — 소스 우선순위 체인 |
| `scripts/counselor/generate_pdf.py:95-123` | CIViC enrichment mode gate |
| `data/gene_knowledge.json` | LLM 데이터 → 빈 JSON `{"genes":[]}` 으로 초기화 |
| `tests/test_build_gene_knowledge.py` | 기존 테스트 업데이트 |

---

## Task 1: NCBI Gene E-utilities 조회 모듈

**Files:**
- Create: `scripts/tools/sources/__init__.py`
- Create: `scripts/tools/sources/ncbi_gene.py`
- Create: `tests/test_ncbi_gene.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_ncbi_gene.py
import pytest


def test_fetch_gene_summary_known_gene(monkeypatch):
    """TP53 gene summary 조회."""
    from scripts.tools.sources.ncbi_gene import fetch_gene_summary

    # Mock esearch response
    mock_search = {"esearchresult": {"idlist": ["7157"]}}
    # Mock esummary response
    mock_summary = {"result": {"7157": {
        "uid": "7157",
        "name": "TP53",
        "description": "tumor protein p53",
        "summary": "This gene encodes a tumor protein that responds to diverse cellular stresses.",
        "otheraliases": "p53, LFS1",
    }}}

    call_count = [0]
    def mock_fetch(url, **kw):
        call_count[0] += 1
        if "esearch" in url:
            return mock_search
        if "esummary" in url:
            return mock_summary
        return None

    monkeypatch.setattr("scripts.tools.sources.ncbi_gene.fetch_with_retry", mock_fetch)

    result = fetch_gene_summary("TP53")
    assert result is not None
    assert result["gene"] == "TP53"
    assert result["entrez_id"] == "7157"
    assert "tumor protein" in result["summary"].lower() or "tumor protein" in result["full_name"]
    assert result["full_name"] == "tumor protein p53"


def test_fetch_gene_summary_unknown_gene(monkeypatch):
    """존재하지 않는 유전자는 None 반환."""
    from scripts.tools.sources.ncbi_gene import fetch_gene_summary

    mock_search = {"esearchresult": {"idlist": []}}
    monkeypatch.setattr(
        "scripts.tools.sources.ncbi_gene.fetch_with_retry",
        lambda url, **kw: mock_search,
    )

    result = fetch_gene_summary("FAKEGENE_XYZ")
    assert result is None


def test_fetch_gene_summary_api_failure(monkeypatch):
    """API 실패 시 None 반환 (graceful)."""
    from scripts.tools.sources.ncbi_gene import fetch_gene_summary

    monkeypatch.setattr(
        "scripts.tools.sources.ncbi_gene.fetch_with_retry",
        lambda url, **kw: None,
    )

    result = fetch_gene_summary("TP53")
    assert result is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_ncbi_gene.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Create sources package and implement ncbi_gene.py**

```python
# scripts/tools/sources/__init__.py
# (empty)
```

```python
# scripts/tools/sources/ncbi_gene.py
"""NCBI Gene E-utilities — fetch gene summaries from NCBI Gene database."""

import logging
import time
from typing import Dict, Optional
from scripts.common.api_utils import fetch_with_retry
from scripts.common.config import get

logger = logging.getLogger(__name__)

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def _get_api_key() -> str:
    return get("api.ncbi_api_key", "") or ""


def fetch_gene_summary(gene_symbol: str) -> Optional[Dict]:
    """Fetch gene summary from NCBI Gene database.

    Returns: {
        "gene": "TP53",
        "entrez_id": "7157",
        "full_name": "tumor protein p53",
        "summary": "This gene encodes...",
        "aliases": "p53, LFS1",
    } or None on failure.
    """
    api_key = _get_api_key()
    key_param = f"&api_key={api_key}" if api_key else ""

    # Step 1: Search for gene ID
    search_url = (
        f"{ESEARCH}?db=gene&term={gene_symbol}[Gene Name]"
        f"+AND+Homo+sapiens[Organism]&retmode=json{key_param}"
    )
    search_data = fetch_with_retry(search_url)
    if not search_data:
        return None

    id_list = search_data.get("esearchresult", {}).get("idlist", [])
    if not id_list:
        return None

    gene_id = id_list[0]

    # Rate limit: 3/sec without key, 10/sec with key
    time.sleep(0.35 if not api_key else 0.1)

    # Step 2: Fetch summary
    summary_url = f"{ESUMMARY}?db=gene&id={gene_id}&retmode=json{key_param}"
    summary_data = fetch_with_retry(summary_url)
    if not summary_data:
        return None

    result = summary_data.get("result", {}).get(gene_id)
    if not result:
        return None

    return {
        "gene": gene_symbol,
        "entrez_id": gene_id,
        "full_name": result.get("description", ""),
        "summary": result.get("summary", ""),
        "aliases": result.get("otheraliases", ""),
    }
```

- [ ] **Step 4: Run tests**

Run: `python -m pytest tests/test_ncbi_gene.py -v`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/tools/sources/__init__.py scripts/tools/sources/ncbi_gene.py tests/test_ncbi_gene.py
git commit -m "feat: NCBI Gene E-utilities module for gene summary lookup"
```

---

## Task 2: GeneReviews PMID 자동 조회 모듈

**Files:**
- Create: `scripts/tools/sources/genreviews.py`
- Create: `tests/test_genreviews.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_genreviews.py
import pytest


def test_fetch_genreviews_known_gene(monkeypatch):
    """TP53에 대한 GeneReviews PMID 조회."""
    from scripts.tools.sources.genreviews import fetch_genreviews_info

    mock_search = {"esearchresult": {"idlist": ["20301371"]}}
    mock_summary = {"result": {"20301371": {
        "uid": "20301371",
        "title": "Li-Fraumeni Syndrome",
        "source": "GeneReviews",
        "bookshelfaccession": "NBK1311",
    }}}

    call_count = [0]
    def mock_fetch(url, **kw):
        call_count[0] += 1
        if "esearch" in url:
            return mock_search
        if "esummary" in url:
            return mock_summary
        return None

    monkeypatch.setattr("scripts.tools.sources.genreviews.fetch_with_retry", mock_fetch)

    result = fetch_genreviews_info("TP53")
    assert result is not None
    assert result["pmid"] == "20301371"
    assert "Li-Fraumeni" in result["title"]
    assert result["nbk_id"] == "NBK1311"
    assert "ncbi.nlm.nih.gov/books/NBK1311" in result["url"]


def test_fetch_genreviews_no_entry(monkeypatch):
    """GeneReviews에 없는 유전자는 None."""
    from scripts.tools.sources.genreviews import fetch_genreviews_info

    mock_search = {"esearchresult": {"idlist": []}}
    monkeypatch.setattr(
        "scripts.tools.sources.genreviews.fetch_with_retry",
        lambda url, **kw: mock_search,
    )

    result = fetch_genreviews_info("OBSCURE_GENE_123")
    assert result is None


def test_fetch_genreviews_api_failure(monkeypatch):
    """API 실패 시 None."""
    from scripts.tools.sources.genreviews import fetch_genreviews_info

    monkeypatch.setattr(
        "scripts.tools.sources.genreviews.fetch_with_retry",
        lambda url, **kw: None,
    )
    assert fetch_genreviews_info("TP53") is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_genreviews.py -v`
Expected: FAIL

- [ ] **Step 3: Implement genreviews.py**

```python
# scripts/tools/sources/genreviews.py
"""GeneReviews PMID lookup via PubMed E-utilities."""

import logging
import time
from typing import Dict, Optional
from scripts.common.api_utils import fetch_with_retry
from scripts.common.config import get

logger = logging.getLogger(__name__)

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def fetch_genreviews_info(gene_symbol: str) -> Optional[Dict]:
    """Search PubMed for GeneReviews entry for a gene.

    Returns: {
        "gene": "TP53",
        "pmid": "20301371",
        "title": "Li-Fraumeni Syndrome",
        "nbk_id": "NBK1311",
        "url": "https://www.ncbi.nlm.nih.gov/books/NBK1311/",
    } or None.
    """
    api_key = get("api.ncbi_api_key", "") or ""
    key_param = f"&api_key={api_key}" if api_key else ""

    # Search PubMed for GeneReviews entries mentioning this gene
    search_url = (
        f"{ESEARCH}?db=pubmed"
        f"&term={gene_symbol}[Title]+AND+GeneReviews[Book]"
        f"&retmode=json&retmax=1{key_param}"
    )
    search_data = fetch_with_retry(search_url)
    if not search_data:
        return None

    id_list = search_data.get("esearchresult", {}).get("idlist", [])
    if not id_list:
        return None

    pmid = id_list[0]

    time.sleep(0.35 if not api_key else 0.1)

    # Fetch article summary
    summary_url = f"{ESUMMARY}?db=pubmed&id={pmid}&retmode=json{key_param}"
    summary_data = fetch_with_retry(summary_url)
    if not summary_data:
        return {"gene": gene_symbol, "pmid": pmid, "title": "", "nbk_id": "", "url": ""}

    article = summary_data.get("result", {}).get(pmid, {})
    title = article.get("title", "")
    nbk_id = article.get("bookshelfaccession", "")

    # Try to extract NBK from elocationid if bookshelfaccession not present
    if not nbk_id:
        eloc = article.get("elocationid", "")
        if "NBK" in eloc:
            nbk_id = eloc.split("NBK")[-1].split(".")[0]
            nbk_id = f"NBK{nbk_id}"

    url = f"https://www.ncbi.nlm.nih.gov/books/{nbk_id}/" if nbk_id else ""

    return {
        "gene": gene_symbol,
        "pmid": pmid,
        "title": title,
        "nbk_id": nbk_id,
        "url": url,
    }
```

- [ ] **Step 4: Run tests**

Run: `python -m pytest tests/test_genreviews.py -v`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/tools/sources/genreviews.py tests/test_genreviews.py
git commit -m "feat: GeneReviews PMID auto-lookup via PubMed E-utilities"
```

---

## Task 3: build_gene_knowledge.py 리라이트

기존 CPIC-only 빌더를 소스 우선순위 체인으로 리라이트.

**Files:**
- Modify: `scripts/tools/build_gene_knowledge.py` (전체 리라이트)
- Create: `tests/test_build_gene_knowledge_v2.py`
- Modify: `tests/test_build_gene_knowledge.py` (기존 테스트 업데이트)

- [ ] **Step 1: Write failing tests**

```python
# tests/test_build_gene_knowledge_v2.py
import json
import pytest


def test_build_from_civic_local(tmp_path, monkeypatch):
    """CIViC local DB에서 gene description을 가져와 knowledge 생성."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    # Mock CIViC gene summary
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_gene_summary",
        lambda gene: {"gene": gene, "description": f"{gene} is a key cancer gene.", "aliases": ""}
        if gene == "BRAF" else None,
    )
    # Mock NCBI Gene (fallback)
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_gene_summary",
        lambda gene: None,
    )
    # Mock GeneReviews
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_genreviews_info",
        lambda gene: None,
    )
    # Mock ClinGen
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_gene_validity_local",
        lambda gene, **kw: None,
    )
    # Mock CIViC treatment
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_treatment_summary",
        lambda gene, variant=None: "Level A: Vemurafenib — Sensitivity in Melanoma",
    )
    # Mock CIViC evidence for refs
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_variant_evidence",
        lambda gene, variant=None: [{"pmid": "20979469", "citation": "Chapman 2011",
         "evidence_type": "Predictive", "significance": "Sensitivity"}],
    )

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["BRAF"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    genes = {g["gene"]: g for g in data["genes"]}
    assert "BRAF" in genes
    assert genes["BRAF"]["content_status"] == "curated-civic"
    assert "cancer gene" in genes["BRAF"]["finding_summary"]
    assert "Vemurafenib" in genes["BRAF"]["treatment_strategies"]


def test_build_ncbi_fallback(tmp_path, monkeypatch):
    """CIViC에 없는 유전자는 NCBI Gene fallback."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_gene_summary",
        lambda gene: None,
    )
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_gene_summary",
        lambda gene: {"gene": gene, "entrez_id": "999", "full_name": "novel gene",
                      "summary": "This gene is involved in cell signaling.", "aliases": ""}
        if gene == "NOVELGENE" else None,
    )
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_genreviews_info",
        lambda gene: None,
    )
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_gene_validity_local",
        lambda gene, **kw: None,
    )
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_treatment_summary",
        lambda gene, variant=None: "",
    )
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.get_variant_evidence",
        lambda gene, variant=None: [],
    )

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["NOVELGENE"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    genes = {g["gene"]: g for g in data["genes"]}
    assert genes["NOVELGENE"]["content_status"] == "curated-ncbi"
    assert "cell signaling" in genes["NOVELGENE"]["finding_summary"]


def test_build_cpic_priority(tmp_path, monkeypatch):
    """PGx 유전자는 CPIC가 우선."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_cpic_gene",
        lambda gene: {"gene": gene, "full_name": "CYP2C19 enzyme",
                      "finding_summary": "CPIC curated", "content_status": "curated-cpic",
                      "references": [{"pmid": "34216116", "source": "CPIC"}],
                      "treatment_strategies": "CPIC guidelines", "frequency_prognosis": "",
                      "function_summary": "", "clinical_significance": "",
                      "associated_conditions": [], "korean_specific_note": None, "hgvs": {}}
        if gene == "CYP2C19" else None,
    )
    # Mock others to return data (should be ignored for PGx genes)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_genreviews_info", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_validity_local", lambda g, **kw: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_treatment_summary", lambda g, v=None: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_variant_evidence", lambda g, v=None: [])

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["CYP2C19"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    genes = {g["gene"]: g for g in data["genes"]}
    assert genes["CYP2C19"]["content_status"] == "curated-cpic"


def test_build_minimal_fallback(tmp_path, monkeypatch):
    """모든 소스 실패 시 auto-minimal entry 생성."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_genreviews_info", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_validity_local", lambda g, **kw: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_cpic_gene", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_treatment_summary", lambda g, v=None: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_variant_evidence", lambda g, v=None: [])

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["UNKNOWN_GENE"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    genes = {g["gene"]: g for g in data["genes"]}
    assert genes["UNKNOWN_GENE"]["content_status"] == "auto-minimal"


def test_build_genreviews_adds_reference(tmp_path, monkeypatch):
    """GeneReviews PMID가 references에 추가됨."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_summary", lambda g: None)
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_gene_summary",
        lambda g: {"gene": g, "entrez_id": "7157", "full_name": "tumor protein p53",
                    "summary": "Encodes a tumor suppressor.", "aliases": ""} if g == "TP53" else None,
    )
    monkeypatch.setattr(
        "scripts.tools.build_gene_knowledge.fetch_genreviews_info",
        lambda g: {"gene": g, "pmid": "20301371", "title": "Li-Fraumeni Syndrome",
                    "nbk_id": "NBK1311", "url": "https://www.ncbi.nlm.nih.gov/books/NBK1311/"}
        if g == "TP53" else None,
    )
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_validity_local", lambda g, **kw: "Definitive" if g == "TP53" else None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_cpic_gene", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_treatment_summary", lambda g, v=None: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_variant_evidence", lambda g, v=None: [])

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["TP53"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    tp53 = {g["gene"]: g for g in data["genes"]}["TP53"]
    pmids = [r.get("pmid") for r in tp53.get("references", [])]
    assert "20301371" in pmids
    assert "Definitive" in tp53["finding_summary"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_build_gene_knowledge_v2.py -v`
Expected: FAIL

- [ ] **Step 3: Rewrite build_gene_knowledge.py**

```python
#!/usr/bin/env python3
"""Build gene_knowledge.json from curated sources.

Source priority per gene:
1. CPIC (if PGx gene) → content_status: curated-cpic
2. CIViC local DB (if description exists) → curated-civic
3. NCBI Gene API (summary) → curated-ncbi
4. Minimal entry → auto-minimal

Additional enrichment (all genes):
- GeneReviews PMID → references
- ClinGen validity → finding_summary
- CIViC treatment evidence → treatment_strategies
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional
from scripts.common.api_utils import fetch_with_retry
from scripts.common.config import get
from scripts.db.query_civic import (
    get_gene_summary, get_treatment_summary, get_variant_evidence,
)
from scripts.tools.sources.ncbi_gene import fetch_gene_summary
from scripts.tools.sources.genreviews import fetch_genreviews_info

logger = logging.getLogger(__name__)

CPIC_API = "https://api.cpicpgx.org/v1"
CPIC_PGX_GENES = [
    "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5", "DPYD",
    "NUDT15", "TPMT", "UGT1A1", "SLCO1B1", "VKORC1",
    "HLA-B", "HLA-A",
]


def fetch_cpic_gene(gene: str) -> Optional[Dict]:
    """Fetch gene info from CPIC API (unchanged from existing)."""
    data = fetch_with_retry(f"{CPIC_API}/gene?symbol=eq.{gene}&select=*")
    if not data:
        return None
    record = data[0] if isinstance(data, list) and data else data
    if not record or not record.get("symbol"):
        return None

    guidelines = fetch_with_retry(
        f"{CPIC_API}/pair?genesymbol=eq.{gene}&select=drugname,cpicLevel,pgkbLevel,guideline(name,url)"
    )
    drugs, refs = [], []
    if guidelines and isinstance(guidelines, list):
        for g in guidelines:
            drug = g.get("drugname", "")
            if drug:
                drugs.append(drug)
            gl = g.get("guideline")
            if gl and isinstance(gl, dict) and gl.get("name"):
                refs.append({"source": f"CPIC Guideline: {gl['name']}",
                             "note": f"CPIC Level {g.get('cpicLevel', 'N/A')}", "pmid": ""})

    treatment = f"CPIC guidelines available for: {', '.join(drugs)}." if drugs else ""
    return {
        "gene": gene, "full_name": record.get("name", ""),
        "function_summary": record.get("functionExampleSubstratesDrugs", ""),
        "clinical_significance": "Pharmacogenomically relevant gene with CPIC guidelines.",
        "associated_conditions": [f"{d} response" for d in drugs[:5]],
        "treatment_strategies": treatment, "frequency_prognosis": "",
        "finding_summary": f"{gene} is a pharmacogene with CPIC-level evidence for drug dosing.",
        "korean_specific_note": None, "hgvs": record.get("hgvs", {}),
        "references": refs or [{"source": "CPIC (cpicpgx.org)", "pmid": "", "note": "Auto-sourced"}],
        "content_status": "curated-cpic",
    }


def _try_clingen_validity(gene: str) -> Optional[str]:
    """Try local ClinGen DB for gene validity."""
    try:
        from scripts.db.query_local_clingen import get_gene_validity_local
        return get_gene_validity_local(gene)
    except Exception:
        return None

# Alias for monkeypatching in tests
get_gene_validity_local = _try_clingen_validity


def _build_gene_entry(gene: str) -> Dict:
    """Build a single gene knowledge entry using source priority chain."""

    # 1. CPIC (PGx genes)
    if gene in CPIC_PGX_GENES:
        cpic = fetch_cpic_gene(gene)
        if cpic:
            logger.info(f"  {gene}: CPIC (PGx)")
            return cpic

    # 2. CIViC local DB (cancer genes with descriptions)
    civic = get_gene_summary(gene)
    civic_description = civic.get("description", "") if civic else ""

    # 3. NCBI Gene API (universal fallback)
    ncbi = fetch_gene_summary(gene)

    # 4. GeneReviews PMID
    genreviews = fetch_genreviews_info(gene)

    # 5. ClinGen validity
    clingen = get_gene_validity_local(gene)

    # 6. CIViC treatment evidence
    treatment = get_treatment_summary(gene)

    # 7. CIViC evidence references
    civic_evidence = get_variant_evidence(gene)
    refs = []
    if civic_evidence:
        for e in civic_evidence[:5]:
            if e.get("pmid"):
                refs.append({"pmid": e["pmid"], "source": e.get("citation", ""),
                             "note": f"{e.get('evidence_type', '')} — {e.get('significance', '')}"})
    if genreviews:
        refs.append({"pmid": genreviews["pmid"], "source": "GeneReviews",
                      "note": genreviews.get("title", "")})

    # Compose finding_summary
    if civic_description:
        finding_summary = civic_description[:500]
        content_status = "curated-civic"
        source_label = "CIViC"
    elif ncbi and ncbi.get("summary"):
        finding_summary = ncbi["summary"][:500]
        content_status = "curated-ncbi"
        source_label = "NCBI Gene"
    else:
        finding_summary = f"{gene} — limited data available."
        content_status = "auto-minimal"
        source_label = "minimal"

    # Append ClinGen validity to finding_summary
    if clingen:
        finding_summary += f" ClinGen gene-disease validity: {clingen}."

    full_name = ""
    if ncbi:
        full_name = ncbi.get("full_name", "")
    elif civic:
        full_name = civic.get("aliases", "")

    logger.info(f"  {gene}: {source_label}" +
                (f" + GeneReviews" if genreviews else "") +
                (f" + ClinGen:{clingen}" if clingen else ""))

    return {
        "gene": gene,
        "full_name": full_name,
        "function_summary": ncbi.get("summary", "")[:300] if ncbi else "",
        "clinical_significance": "",
        "associated_conditions": [],
        "treatment_strategies": treatment or "",
        "frequency_prognosis": "",
        "finding_summary": finding_summary,
        "korean_specific_note": None,
        "hgvs": {},
        "references": refs,
        "content_status": content_status,
    }


def build_knowledge(
    gene_list: List[str],
    output_path: str = "data/gene_knowledge.json",
) -> str:
    """Build gene_knowledge.json from curated sources."""
    logger.info(f"Building knowledge for {len(gene_list)} genes...")

    genes = []
    for gene in gene_list:
        entry = _build_gene_entry(gene)
        genes.append(entry)

    result = {"genes": genes}
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    curated = sum(1 for g in genes if "curated" in g.get("content_status", ""))
    logger.info(f"Gene knowledge built: {len(genes)} genes ({curated} curated) → {output_path}")
    return output_path


def _load_gene_list_from_vcf(vcf_path: str) -> List[str]:
    """Extract unique gene symbols from a VCF file."""
    from scripts.intake.parse_vcf import parse_vcf
    variants = parse_vcf(vcf_path)
    return sorted(set(v.gene for v in variants if v.gene))


def _load_oncokb_genes() -> List[str]:
    """Load gene list from OncoKB cancer genes JSON."""
    path = get("paths.oncokb_genes", "data/oncokb_cancer_genes.json")
    try:
        with open(path) as f:
            data = json.load(f)
        return sorted(data.get("genes", {}).keys())
    except (FileNotFoundError, json.JSONDecodeError):
        return []


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description="Build gene_knowledge.json from curated sources")
    parser.add_argument("--genes", help="Comma-separated gene list")
    parser.add_argument("--vcf", help="Extract genes from VCF file")
    parser.add_argument("--source", choices=["oncokb", "all"], help="Use predefined gene list")
    parser.add_argument("--output", default="data/gene_knowledge.json")
    args = parser.parse_args()

    if args.genes:
        gene_list = [g.strip() for g in args.genes.split(",")]
    elif args.vcf:
        gene_list = _load_gene_list_from_vcf(args.vcf)
    elif args.source == "oncokb":
        gene_list = _load_oncokb_genes()
    elif args.source == "all":
        gene_list = _load_oncokb_genes() + CPIC_PGX_GENES
        gene_list = sorted(set(gene_list))
    else:
        parser.error("Specify --genes, --vcf, or --source")

    build_knowledge(gene_list, args.output)
```

- [ ] **Step 4: Update existing test file**

`tests/test_build_gene_knowledge.py`의 기존 2개 테스트가 새 함수 시그니처에 맞도록 업데이트:
- `test_fetch_cpic_gene_info` — 변경 없음 (fetch_cpic_gene은 동일)
- `test_build_gene_knowledge_merges` — `build_knowledge(existing_path, output_path, cpic_genes)` → `build_knowledge(gene_list, output_path)` 시그니처 변경에 맞게 수정

- [ ] **Step 5: Run all tests**

Run: `python -m pytest tests/test_build_gene_knowledge_v2.py tests/test_build_gene_knowledge.py -v`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add scripts/tools/build_gene_knowledge.py tests/test_build_gene_knowledge_v2.py tests/test_build_gene_knowledge.py
git commit -m "feat: rewrite gene knowledge builder — dynamic sourcing from CIViC/NCBI/GeneReviews/ClinGen"
```

---

## Task 4: generate_pdf.py Mode Gate

CIViC enrichment을 cancer mode에서만 적용.

**Files:**
- Modify: `scripts/counselor/generate_pdf.py:95-123`
- Modify: `tests/test_generate_pdf.py`

- [ ] **Step 1: Write failing test**

```python
# tests/test_generate_pdf.py (추가)

def test_rare_disease_no_civic_enrichment():
    """Rare disease mode에서 CIViC finding_summary가 적용되지 않아야 함."""
    import copy
    from scripts.counselor.generate_pdf import generate_report_html

    data = copy.deepcopy(MINIMAL_REPORT)
    data["variants"][0]["finding_summary"] = "Original gene knowledge text"
    data["variants"][0]["gene"] = "TP53"
    html = generate_report_html(data, mode="rare-disease")
    # In rare disease mode, CIViC should NOT override finding_summary
    assert "Original gene knowledge text" in html
```

- [ ] **Step 2: Modify generate_pdf.py — add mode parameter and gate**

현재 `generate_report_html(report_data, mode="cancer")` 시그니처 확인 후, CIViC enrichment 블록을 mode gate로 감싸기:

```python
        # CIViC enrichment — cancer mode only
        if mode == "cancer":
            civic_gene = get_gene_summary(gene)
            if civic_gene and civic_gene.get("description"):
                v.setdefault("finding_summary", civic_gene["description"][:500])

            # Treatment from CIViC
            hgvsp = v.get("hgvsp", "")
            civic_variant_name = _hgvsp_to_civic_variant(hgvsp)
            civic_treatment = get_treatment_summary(gene, civic_variant_name)
            if civic_treatment:
                v.setdefault("treatment_strategies", civic_treatment)

            # Evidence references from CIViC
            civic_evidence = get_variant_evidence(gene, civic_variant_name) if civic_variant_name else []
            if not civic_evidence:
                civic_evidence = get_variant_evidence(gene)
            if civic_evidence:
                refs = [
                    {"pmid": e["pmid"], "source": e["citation"],
                     "note": f"{e['evidence_type']} — {e['significance']}"}
                    for e in civic_evidence[:5] if e["pmid"]
                ]
                if refs:
                    v.setdefault("references", refs)
                    v.setdefault("content_status", "curated-civic")
```

- [ ] **Step 3: Run tests**

Run: `python -m pytest tests/test_generate_pdf.py -v`
Expected: ALL PASS

- [ ] **Step 4: Commit**

```bash
git add scripts/counselor/generate_pdf.py tests/test_generate_pdf.py
git commit -m "feat: CIViC enrichment gated to cancer mode only — rare disease uses gene_knowledge directly"
```

---

## Task 5: LLM gene_knowledge.json 교체 + 통합 테스트

**Files:**
- Modify: `data/gene_knowledge.json` (LLM → 빈 JSON)
- Modify: `tests/test_orchestrate.py` (통합 테스트)

- [ ] **Step 1: Backup and replace gene_knowledge.json**

```bash
cp data/gene_knowledge.json data/gene_knowledge.legacy.json
echo '{"genes":[]}' > data/gene_knowledge.json
```

- [ ] **Step 2: Write integration test**

```python
# tests/test_orchestrate.py (추가)

def test_pipeline_works_with_empty_gene_knowledge(tmp_path):
    """gene_knowledge.json이 비어있어도 파이프라인 정상 동작."""
    from scripts.orchestrate import run_pipeline

    result = run_pipeline(
        vcf_path="data/sample_vcf/demo_variants_grch38_annotated.vcf",
        output_path=str(tmp_path / "report.html"),
        skip_api=True,
        mode="cancer",
    )
    assert result is not None
    assert len(result["variants"]) > 0
    # CIViC enrichment still works at runtime for cancer
    for v in result["variants"]:
        assert "classification" in v
        assert "tier" in v
```

- [ ] **Step 3: Run full test suite**

Run: `python -m pytest tests/ --tb=short`
Expected: 405+ passed (gene_knowledge 비어있어도 CIViC runtime fallback으로 동작)

- [ ] **Step 4: Commit**

```bash
git add data/gene_knowledge.json data/gene_knowledge.legacy.json tests/test_orchestrate.py
git commit -m "feat: replace LLM gene_knowledge with empty JSON — curated build replaces static data"
```

---

## Task 6: 실 데이터 빌드 테스트 + 문서

OncoKB 유전자 리스트로 실제 knowledge를 빌드하고 검증.

**Files:**
- Run: `scripts/tools/build_gene_knowledge.py --source oncokb`
- Modify: `docs/SETUP.md` (빌드 가이드)

- [ ] **Step 1: Run actual build with a few test genes**

```bash
python scripts/tools/build_gene_knowledge.py --genes TP53,BRCA2,KRAS,BRAF,EGFR --output data/gene_knowledge.json
```

- [ ] **Step 2: Verify built JSON**

```bash
python -c "
import json
with open('data/gene_knowledge.json') as f:
    data = json.load(f)
for g in data['genes']:
    print(f\"{g['gene']:10s} status={g['content_status']:20s} refs={len(g.get('references',[]))} summary={g['finding_summary'][:60]}...\")
"
```

- [ ] **Step 3: Run reports with built knowledge**

```bash
python scripts/orchestrate.py data/sample_vcf/demo_variants_grch38_annotated.vcf \
  -o output/test_curated_knowledge.html --json --mode cancer --skip-api
```

- [ ] **Step 4: Add build instructions to SETUP.md**

```markdown
### Gene Knowledge Base (Curated Sources)

Build the gene knowledge database from authoritative sources (NCBI Gene, CIViC, GeneReviews, CPIC):

\`\`\`bash
# Build for OncoKB cancer genes (~300 genes, ~5 min)
python scripts/tools/build_gene_knowledge.py --source oncokb

# Build for specific genes
python scripts/tools/build_gene_knowledge.py --genes TP53,BRCA2,KRAS

# Build from VCF (extract gene list automatically)
python scripts/tools/build_gene_knowledge.py --vcf input/sample.vcf

# Build all (OncoKB + PGx genes)
python scripts/tools/build_gene_knowledge.py --source all
\`\`\`
```

- [ ] **Step 5: Commit**

```bash
git add data/gene_knowledge.json docs/SETUP.md
git commit -m "feat: curated gene knowledge for test genes + build instructions"
```

---

## Summary

| Task | 내용 | 신규/수정 | 예상 테스트 |
|------|------|----------|------------|
| 1 | NCBI Gene E-utilities 모듈 | 신규 3 | +3 tests |
| 2 | GeneReviews PMID 조회 모듈 | 신규 2 | +3 tests |
| 3 | build_gene_knowledge.py 리라이트 | 수정 1, 신규 1, 수정 1 | +5 tests |
| 4 | generate_pdf.py mode gate | 수정 2 | +1 test |
| 5 | LLM 데이터 교체 + 통합 테스트 | 수정 2, 신규 1 | +1 test |
| 6 | 실 빌드 + 문서 | 수정 2 | — |

**예상 결과:** 기존 405 + 약 13 신규 = **418+ tests**
