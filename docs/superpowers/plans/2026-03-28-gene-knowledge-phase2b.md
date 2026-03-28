# Gene Knowledge Phase 2b: Orphanet + GeneReviews Files + OMIM Mapping

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Orphanet 질환 유병률, GeneReviews gene→NBK 매핑, OMIM gene→MIM 매핑을 로컬 DB로 구축하고 gene knowledge builder에 통합한다.

**Architecture:** 3개 독립 데이터 소스를 각각 다운로드→파싱→SQLite 패턴으로 구축. `build_gene_knowledge.py`에 새 소스를 추가하여 `frequency_prognosis`(Orphanet 유병률)와 `references`(GeneReviews NBK)를 강화.

**Tech Stack:** Python 3, SQLite3, XML (ElementTree), requests, pytest

---

## File Structure

### 신규 파일
| File | Responsibility |
|------|---------------|
| `scripts/db/build_orphanet_db.py` | Orphanet prevalence XML → SQLite |
| `scripts/db/query_orphanet.py` | Orphanet SQLite 조회 (질환 유병률) |
| `scripts/db/build_genreviews_db.py` | GeneReviews TSV 파일 → SQLite (gene→NBK→PMID) |
| `scripts/db/query_genreviews.py` | GeneReviews SQLite 조회 |
| `scripts/db/build_omim_mapping.py` | OMIM mim2gene.txt → SQLite (gene→MIM) |
| `scripts/db/query_omim_mapping.py` | OMIM mapping SQLite 조회 |
| `tests/test_orphanet.py` | Orphanet 빌드+쿼리 테스트 |
| `tests/test_genreviews_db.py` | GeneReviews 빌드+쿼리 테스트 |
| `tests/test_omim_mapping.py` | OMIM mapping 빌드+쿼리 테스트 |

### 수정 파일
| File | Change |
|------|--------|
| `scripts/tools/build_gene_knowledge.py` | Orphanet prevalence + GeneReviews local DB + OMIM MIM 통합 |
| `config.yaml` | orphanet_db, genreviews_db, omim_mapping_db 경로 추가 |
| `docs/SETUP.md` | 다운로드+빌드 가이드 추가 |

---

## Task 1: Orphanet Prevalence XML → SQLite

Orphanet Product 9 (en_product9_prev.xml, ~15MB)를 다운로드하여 질환별 유병률 SQLite로 변환.

**Files:**
- Create: `scripts/db/build_orphanet_db.py`
- Create: `scripts/db/query_orphanet.py`
- Create: `tests/test_orphanet.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_orphanet.py
import sqlite3
import pytest
from pathlib import Path


def _create_sample_orphanet_xml(path):
    """Minimal Orphanet Product 9 XML for testing."""
    Path(path).write_text('''<?xml version="1.0" encoding="UTF-8"?>
<JDBOR date="2025-12-01" version="1.0">
  <DisorderList count="2">
    <Disorder id="1">
      <OrphaCode>166024</OrphaCode>
      <Name lang="en">Cystic fibrosis</Name>
      <DisorderType><Name lang="en">Disease</Name></DisorderType>
      <DisorderGeneAssociationList count="1">
        <DisorderGeneAssociation>
          <Gene><Symbol>CFTR</Symbol><Name lang="en">CF transmembrane conductance regulator</Name></Gene>
        </DisorderGeneAssociation>
      </DisorderGeneAssociationList>
      <PrevalenceList count="1">
        <Prevalence>
          <PrevalenceType><Name lang="en">Point prevalence</Name></PrevalenceType>
          <PrevalenceQualification><Name lang="en">Value and class</Name></PrevalenceQualification>
          <PrevalenceClass><Name lang="en">1-5 / 10 000</Name></PrevalenceClass>
          <ValMoy>3.0</ValMoy>
          <PrevalenceGeographic><Name lang="en">Europe</Name></PrevalenceGeographic>
        </Prevalence>
      </PrevalenceList>
    </Disorder>
    <Disorder id="2">
      <OrphaCode>558</OrphaCode>
      <Name lang="en">Marfan syndrome</Name>
      <DisorderType><Name lang="en">Disease</Name></DisorderType>
      <DisorderGeneAssociationList count="1">
        <DisorderGeneAssociation>
          <Gene><Symbol>FBN1</Symbol><Name lang="en">fibrillin 1</Name></Gene>
        </DisorderGeneAssociation>
      </DisorderGeneAssociationList>
      <PrevalenceList count="1">
        <Prevalence>
          <PrevalenceType><Name lang="en">Point prevalence</Name></PrevalenceType>
          <PrevalenceQualification><Name lang="en">Value and class</Name></PrevalenceQualification>
          <PrevalenceClass><Name lang="en">1-9 / 100 000</Name></PrevalenceClass>
          <ValMoy>1.5</ValMoy>
          <PrevalenceGeographic><Name lang="en">Worldwide</Name></PrevalenceGeographic>
        </Prevalence>
      </PrevalenceList>
    </Disorder>
  </DisorderList>
</JDBOR>''')


def test_build_orphanet_db(tmp_path):
    from scripts.db.build_orphanet_db import build_db
    xml_path = str(tmp_path / "orphanet.xml")
    db_path = str(tmp_path / "orphanet.sqlite3")
    _create_sample_orphanet_xml(xml_path)
    result = build_db(xml_path, db_path)
    assert result == db_path
    conn = sqlite3.connect(db_path)
    rows = conn.execute("SELECT * FROM prevalence WHERE gene_symbol = 'CFTR'").fetchall()
    assert len(rows) >= 1
    conn.close()


def test_query_prevalence_by_gene(tmp_orphanet_db):
    from scripts.db.query_orphanet import get_prevalence_by_gene
    result = get_prevalence_by_gene("CFTR", tmp_orphanet_db)
    assert len(result) >= 1
    assert result[0]["disease_name"] == "Cystic fibrosis"
    assert "1-5 / 10 000" in result[0]["prevalence_class"]


def test_get_prevalence_text(tmp_orphanet_db):
    from scripts.db.query_orphanet import get_prevalence_text
    text = get_prevalence_text("CFTR", tmp_orphanet_db)
    assert "Cystic fibrosis" in text
    assert "1-5 / 10 000" in text


def test_query_unknown_gene(tmp_orphanet_db):
    from scripts.db.query_orphanet import get_prevalence_by_gene
    result = get_prevalence_by_gene("FAKEGENE", tmp_orphanet_db)
    assert result == []


@pytest.fixture
def tmp_orphanet_db(tmp_path):
    from scripts.db.build_orphanet_db import build_db
    xml_path = str(tmp_path / "orphanet.xml")
    db_path = str(tmp_path / "orphanet.sqlite3")
    _create_sample_orphanet_xml(xml_path)
    build_db(xml_path, db_path)
    return db_path
```

- [ ] **Step 2: Implement build_orphanet_db.py**

```python
#!/usr/bin/env python3
"""Build Orphanet prevalence SQLite from Product 9 XML."""

import sqlite3
import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime

logger = logging.getLogger(__name__)
DEFAULT_DB_PATH = "data/db/orphanet.sqlite3"


def build_db(xml_path: str, db_path: str = DEFAULT_DB_PATH) -> str:
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)

    conn.execute("""CREATE TABLE IF NOT EXISTS prevalence (
        orpha_code TEXT,
        disease_name TEXT NOT NULL,
        gene_symbol TEXT,
        gene_name TEXT,
        prevalence_type TEXT,
        prevalence_class TEXT,
        val_moy REAL,
        geographic TEXT
    )""")
    conn.execute("DELETE FROM prevalence")

    tree = ET.parse(xml_path)
    root = tree.getroot()

    for disorder in root.iter("Disorder"):
        orpha = disorder.findtext("OrphaCode", "")
        name = disorder.findtext("Name", "")

        # Extract gene associations
        genes = []
        for ga in disorder.iter("DisorderGeneAssociation"):
            gene_el = ga.find("Gene")
            if gene_el is not None:
                genes.append({
                    "symbol": gene_el.findtext("Symbol", ""),
                    "name": gene_el.findtext("Name", ""),
                })

        # Extract prevalence data
        for prev in disorder.iter("Prevalence"):
            prev_type = ""
            pt = prev.find("PrevalenceType")
            if pt is not None:
                prev_type = pt.findtext("Name", "")
            prev_class = ""
            pc = prev.find("PrevalenceClass")
            if pc is not None:
                prev_class = pc.findtext("Name", "")
            val_moy = prev.findtext("ValMoy", "")
            geo = ""
            pg = prev.find("PrevalenceGeographic")
            if pg is not None:
                geo = pg.findtext("Name", "")

            for gene in (genes or [{"symbol": "", "name": ""}]):
                conn.execute(
                    "INSERT INTO prevalence VALUES (?,?,?,?,?,?,?,?)",
                    (orpha, name, gene["symbol"], gene["name"],
                     prev_type, prev_class,
                     float(val_moy) if val_moy else None, geo),
                )

    conn.execute("CREATE INDEX IF NOT EXISTS idx_prev_gene ON prevalence(gene_symbol)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_prev_orpha ON prevalence(orpha_code)")

    now = datetime.utcnow().isoformat()
    conn.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (now,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('source', 'Orphanet (orphadata.com)')")
    count = conn.execute("SELECT COUNT(*) FROM prevalence").fetchone()[0]
    gene_count = conn.execute("SELECT COUNT(DISTINCT gene_symbol) FROM prevalence WHERE gene_symbol != ''").fetchone()[0]
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('row_count', ?)", (str(count),))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('gene_count', ?)", (str(gene_count),))

    conn.commit()
    conn.close()
    logger.info(f"Orphanet DB built: {gene_count} genes, {count} prevalence entries → {db_path}")
    return db_path


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument("xml_path", nargs="?", default="data/db/en_product9_prev.xml")
    parser.add_argument("--db-path", default=DEFAULT_DB_PATH)
    args = parser.parse_args()
    build_db(args.xml_path, args.db_path)
```

- [ ] **Step 3: Implement query_orphanet.py**

```python
"""Query Orphanet prevalence SQLite database."""

import sqlite3
import logging
from typing import Dict, List, Optional
from scripts.common.config import get

logger = logging.getLogger(__name__)


def _get_db_path(db_path: Optional[str] = None) -> str:
    return db_path or get("paths.orphanet_db") or "data/db/orphanet.sqlite3"


def get_prevalence_by_gene(gene: str, db_path: Optional[str] = None) -> List[Dict]:
    path = _get_db_path(db_path)
    try:
        conn = sqlite3.connect(path)
        rows = conn.execute(
            "SELECT disease_name, prevalence_type, prevalence_class, val_moy, geographic "
            "FROM prevalence WHERE gene_symbol = ? ORDER BY disease_name",
            (gene,),
        ).fetchall()
        conn.close()
        return [
            {"disease_name": r[0], "prevalence_type": r[1],
             "prevalence_class": r[2], "val_moy": r[3], "geographic": r[4]}
            for r in rows
        ]
    except Exception as e:
        logger.debug(f"Orphanet query failed for {gene}: {e}")
        return []


def get_prevalence_text(gene: str, db_path: Optional[str] = None) -> str:
    """Build a human-readable prevalence summary for a gene."""
    entries = get_prevalence_by_gene(gene, db_path)
    if not entries:
        return ""
    parts = []
    seen = set()
    for e in entries:
        key = (e["disease_name"], e["prevalence_class"])
        if key in seen:
            continue
        seen.add(key)
        geo = f" ({e['geographic']})" if e["geographic"] else ""
        parts.append(f"{e['disease_name']}: {e['prevalence_class']}{geo}")
    return "; ".join(parts[:3])
```

- [ ] **Step 4: Run tests, commit**

Run: `python -m pytest tests/test_orphanet.py -v`
Commit: `feat: Orphanet prevalence XML → SQLite build and query`

---

## Task 2: GeneReviews Identifier Files → SQLite

GeneReviews FTP의 TSV 파일을 다운로드하여 gene→NBK→PMID 매핑 SQLite로 변환.

**Data source:** `ftp://ftp.ncbi.nlm.nih.gov/pub/GeneReviews/GRshortname_NBKid_genesymbol_dzname.txt`

**Files:**
- Create: `scripts/db/build_genreviews_db.py`
- Create: `scripts/db/query_genreviews.py`
- Create: `tests/test_genreviews_db.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_genreviews_db.py
import sqlite3
import pytest
from pathlib import Path


def _create_sample_genreviews_tsv(path):
    """Minimal GeneReviews identifier TSV."""
    # Format: GR_shortname\tNBK_id\tgene_symbol\tdisease_name
    Path(path).write_text(
        "#GR_shortname\tNBK_id\tgene_symbol\tdisease_name\n"
        "li-fraumeni\tNBK1311\tTP53\tLi-Fraumeni Syndrome\n"
        "brca1-2\tNBK1247\tBRCA2\tHereditary Breast and Ovarian Cancer\n"
        "brca1-2\tNBK1247\tBRCA1\tHereditary Breast and Ovarian Cancer\n"
        "cftr\tNBK1250\tCFTR\tCFTR-Related Disorders\n"
        "noonan\tNBK1124\tPTPN11\tNoonan Syndrome\n"
    )


def _create_sample_titles_tsv(path):
    """Minimal GeneReviews titles TSV."""
    # Format: GR_shortname\tGR_Title\tNBK_id\tPMID
    Path(path).write_text(
        "#GR_shortname\tGR_Title\tNBK_id\tPMID\n"
        "li-fraumeni\tLi-Fraumeni Syndrome\tNBK1311\t20301371\n"
        "brca1-2\tBRCA1- and BRCA2-Associated Hereditary Breast and Ovarian Cancer\tNBK1247\t20301425\n"
        "cftr\tCFTR-Related Disorders\tNBK1250\t20301428\n"
        "noonan\tNoonan Syndrome\tNBK1124\t20301303\n"
    )


def test_build_genreviews_db(tmp_path):
    from scripts.db.build_genreviews_db import build_db
    genes_path = str(tmp_path / "genes.txt")
    titles_path = str(tmp_path / "titles.txt")
    db_path = str(tmp_path / "genreviews.sqlite3")
    _create_sample_genreviews_tsv(genes_path)
    _create_sample_titles_tsv(titles_path)
    result = build_db(genes_path, titles_path, db_path)
    assert result == db_path
    conn = sqlite3.connect(db_path)
    rows = conn.execute("SELECT * FROM genreviews WHERE gene_symbol = 'TP53'").fetchall()
    assert len(rows) >= 1
    conn.close()


def test_query_genreviews_by_gene(tmp_genreviews_db):
    from scripts.db.query_genreviews import get_genreviews_for_gene
    result = get_genreviews_for_gene("TP53", tmp_genreviews_db)
    assert result is not None
    assert result["nbk_id"] == "NBK1311"
    assert result["pmid"] == "20301371"
    assert "Li-Fraumeni" in result["title"]
    assert "ncbi.nlm.nih.gov/books/NBK1311" in result["url"]


def test_query_genreviews_unknown(tmp_genreviews_db):
    from scripts.db.query_genreviews import get_genreviews_for_gene
    assert get_genreviews_for_gene("FAKEGENE", tmp_genreviews_db) is None


def test_query_brca_shared_entry(tmp_genreviews_db):
    """BRCA1과 BRCA2는 같은 GeneReviews entry를 공유."""
    from scripts.db.query_genreviews import get_genreviews_for_gene
    r1 = get_genreviews_for_gene("BRCA1", tmp_genreviews_db)
    r2 = get_genreviews_for_gene("BRCA2", tmp_genreviews_db)
    assert r1["nbk_id"] == r2["nbk_id"] == "NBK1247"


@pytest.fixture
def tmp_genreviews_db(tmp_path):
    from scripts.db.build_genreviews_db import build_db
    genes_path = str(tmp_path / "genes.txt")
    titles_path = str(tmp_path / "titles.txt")
    db_path = str(tmp_path / "genreviews.sqlite3")
    _create_sample_genreviews_tsv(genes_path)
    _create_sample_titles_tsv(titles_path)
    build_db(genes_path, titles_path, db_path)
    return db_path
```

- [ ] **Step 2: Implement build_genreviews_db.py**

```python
#!/usr/bin/env python3
"""Build GeneReviews SQLite from NCBI FTP identifier files."""

import sqlite3
import logging
from pathlib import Path
from datetime import datetime

logger = logging.getLogger(__name__)
DEFAULT_DB_PATH = "data/db/genreviews.sqlite3"


def build_db(genes_path: str, titles_path: str, db_path: str = DEFAULT_DB_PATH) -> str:
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)

    conn.execute("""CREATE TABLE IF NOT EXISTS genreviews (
        gene_symbol TEXT NOT NULL,
        nbk_id TEXT NOT NULL,
        gr_shortname TEXT,
        disease_name TEXT,
        title TEXT,
        pmid TEXT
    )""")
    conn.execute("DELETE FROM genreviews")

    # Load titles (shortname → title + PMID)
    titles = {}
    with open(titles_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                titles[parts[0]] = {"title": parts[1], "nbk_id": parts[2], "pmid": parts[3]}

    # Load gene-disease mappings
    with open(genes_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            shortname, nbk_id, gene, disease = parts[0], parts[1], parts[2], parts[3]
            title_info = titles.get(shortname, {})
            conn.execute(
                "INSERT INTO genreviews VALUES (?,?,?,?,?,?)",
                (gene, nbk_id, shortname, disease,
                 title_info.get("title", disease),
                 title_info.get("pmid", "")),
            )

    conn.execute("CREATE INDEX IF NOT EXISTS idx_gr_gene ON genreviews(gene_symbol)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_gr_nbk ON genreviews(nbk_id)")

    now = datetime.utcnow().isoformat()
    conn.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (now,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('source', 'GeneReviews (NCBI)')")
    gene_count = conn.execute("SELECT COUNT(DISTINCT gene_symbol) FROM genreviews").fetchone()[0]
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('gene_count', ?)", (str(gene_count),))

    conn.commit()
    conn.close()
    logger.info(f"GeneReviews DB built: {gene_count} genes → {db_path}")
    return db_path


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", default="data/db/GRshortname_NBKid_genesymbol_dzname.txt")
    parser.add_argument("--titles", default="data/db/GRtitle_shortname_NBKid.txt")
    parser.add_argument("--db-path", default=DEFAULT_DB_PATH)
    args = parser.parse_args()
    build_db(args.genes, args.titles, args.db_path)
```

- [ ] **Step 3: Implement query_genreviews.py**

```python
"""Query GeneReviews local SQLite database."""

import sqlite3
import logging
from typing import Dict, Optional
from scripts.common.config import get

logger = logging.getLogger(__name__)


def _get_db_path(db_path: Optional[str] = None) -> str:
    return db_path or get("paths.genreviews_db") or "data/db/genreviews.sqlite3"


def get_genreviews_for_gene(gene: str, db_path: Optional[str] = None) -> Optional[Dict]:
    """Get GeneReviews entry for a gene (NBK ID, PMID, title, URL)."""
    path = _get_db_path(db_path)
    try:
        conn = sqlite3.connect(path)
        row = conn.execute(
            "SELECT nbk_id, pmid, title, disease_name FROM genreviews WHERE gene_symbol = ? LIMIT 1",
            (gene,),
        ).fetchone()
        conn.close()
        if not row:
            return None
        return {
            "gene": gene,
            "nbk_id": row[0],
            "pmid": row[1],
            "title": row[2],
            "disease_name": row[3],
            "url": f"https://www.ncbi.nlm.nih.gov/books/{row[0]}/",
        }
    except Exception as e:
        logger.debug(f"GeneReviews query failed for {gene}: {e}")
        return None
```

- [ ] **Step 4: Run tests, commit**

Run: `python -m pytest tests/test_genreviews_db.py -v`
Commit: `feat: GeneReviews identifier files → SQLite build and query`

---

## Task 3: OMIM mim2gene.txt → SQLite

무료로 제공되는 `mim2gene.txt`를 파싱하여 gene→MIM 매핑 SQLite 구축.

**Data source:** `https://omim.org/static/omim/data/mim2gene.txt`

**Files:**
- Create: `scripts/db/build_omim_mapping.py`
- Create: `scripts/db/query_omim_mapping.py`
- Create: `tests/test_omim_mapping.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_omim_mapping.py
import sqlite3
import pytest
from pathlib import Path


def _create_sample_mim2gene(path):
    """Minimal mim2gene.txt."""
    Path(path).write_text(
        "# Copyright OMIM\n"
        "# Generated: 2026-03-01\n"
        "# MIM Number\tMIM Entry Type\tEntrez Gene ID\tApproved Gene Symbol (HGNC)\tEnsembl Gene ID\n"
        "191170\tgene\t7157\tTP53\tENSG00000141510\n"
        "600185\tgene\t675\tBRCA2\tENSG00000139618\n"
        "602421\tgene\t1080\tCFTR\tENSG00000001626\n"
        "176876\tgene\t5781\tPTPN11\tENSG00000179295\n"
        "114480\tphenotype\t\t\t\n"
        "190685\tgene/phenotype\t7015\tTERT\tENSG00000164362\n"
    )


def test_build_omim_mapping(tmp_path):
    from scripts.db.build_omim_mapping import build_db
    txt_path = str(tmp_path / "mim2gene.txt")
    db_path = str(tmp_path / "omim_mapping.sqlite3")
    _create_sample_mim2gene(txt_path)
    result = build_db(txt_path, db_path)
    assert result == db_path
    conn = sqlite3.connect(db_path)
    rows = conn.execute("SELECT * FROM mim2gene WHERE gene_symbol = 'TP53'").fetchall()
    assert len(rows) >= 1
    conn.close()


def test_get_mim_for_gene(tmp_omim_db):
    from scripts.db.query_omim_mapping import get_mim_for_gene
    result = get_mim_for_gene("TP53", tmp_omim_db)
    assert result is not None
    assert result["mim_number"] == "191170"
    assert result["entrez_id"] == "7157"


def test_get_mim_unknown_gene(tmp_omim_db):
    from scripts.db.query_omim_mapping import get_mim_for_gene
    assert get_mim_for_gene("FAKEGENE", tmp_omim_db) is None


def test_phenotype_entries_excluded(tmp_omim_db):
    """phenotype-only entries (no gene symbol) are excluded."""
    from scripts.db.query_omim_mapping import get_mim_for_gene
    # 114480 is phenotype-only, no gene
    import sqlite3
    conn = sqlite3.connect(tmp_omim_db)
    rows = conn.execute("SELECT * FROM mim2gene WHERE mim_number = '114480'").fetchall()
    conn.close()
    assert len(rows) == 0


@pytest.fixture
def tmp_omim_db(tmp_path):
    from scripts.db.build_omim_mapping import build_db
    txt_path = str(tmp_path / "mim2gene.txt")
    db_path = str(tmp_path / "omim_mapping.sqlite3")
    _create_sample_mim2gene(txt_path)
    build_db(txt_path, db_path)
    return db_path
```

- [ ] **Step 2: Implement build_omim_mapping.py + query_omim_mapping.py**

`build_omim_mapping.py`:
```python
#!/usr/bin/env python3
"""Build OMIM gene mapping SQLite from mim2gene.txt."""

import sqlite3
import logging
from pathlib import Path
from datetime import datetime

logger = logging.getLogger(__name__)
DEFAULT_DB_PATH = "data/db/omim_mapping.sqlite3"


def build_db(txt_path: str, db_path: str = DEFAULT_DB_PATH) -> str:
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)

    conn.execute("""CREATE TABLE IF NOT EXISTS mim2gene (
        mim_number TEXT NOT NULL,
        entry_type TEXT,
        entrez_id TEXT,
        gene_symbol TEXT NOT NULL,
        ensembl_id TEXT
    )""")
    conn.execute("DELETE FROM mim2gene")

    with open(txt_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            gene_symbol = parts[3].strip()
            if not gene_symbol:
                continue
            conn.execute(
                "INSERT INTO mim2gene VALUES (?,?,?,?,?)",
                (parts[0], parts[1], parts[2], gene_symbol,
                 parts[4] if len(parts) > 4 else ""),
            )

    conn.execute("CREATE INDEX IF NOT EXISTS idx_omim_gene ON mim2gene(gene_symbol)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_omim_mim ON mim2gene(mim_number)")

    now = datetime.utcnow().isoformat()
    conn.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (now,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('source', 'OMIM (omim.org)')")
    count = conn.execute("SELECT COUNT(*) FROM mim2gene").fetchone()[0]
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('gene_count', ?)", (str(count),))

    conn.commit()
    conn.close()
    logger.info(f"OMIM mapping DB built: {count} gene entries → {db_path}")
    return db_path


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument("txt_path", nargs="?", default="data/db/mim2gene.txt")
    parser.add_argument("--db-path", default=DEFAULT_DB_PATH)
    args = parser.parse_args()
    build_db(args.txt_path, args.db_path)
```

`query_omim_mapping.py`:
```python
"""Query OMIM gene mapping SQLite database."""

import sqlite3
import logging
from typing import Dict, Optional
from scripts.common.config import get

logger = logging.getLogger(__name__)


def _get_db_path(db_path: Optional[str] = None) -> str:
    return db_path or get("paths.omim_mapping_db") or "data/db/omim_mapping.sqlite3"


def get_mim_for_gene(gene: str, db_path: Optional[str] = None) -> Optional[Dict]:
    path = _get_db_path(db_path)
    try:
        conn = sqlite3.connect(path)
        row = conn.execute(
            "SELECT mim_number, entry_type, entrez_id, ensembl_id "
            "FROM mim2gene WHERE gene_symbol = ? LIMIT 1",
            (gene,),
        ).fetchone()
        conn.close()
        if not row:
            return None
        return {
            "gene": gene,
            "mim_number": row[0],
            "entry_type": row[1],
            "entrez_id": row[2],
            "ensembl_id": row[3],
            "url": f"https://omim.org/entry/{row[0]}",
        }
    except Exception as e:
        logger.debug(f"OMIM mapping query failed for {gene}: {e}")
        return None
```

- [ ] **Step 3: Run tests, commit**

Run: `python -m pytest tests/test_omim_mapping.py -v`
Commit: `feat: OMIM mim2gene.txt → SQLite mapping build and query`

---

## Task 4: build_gene_knowledge.py에 새 소스 통합

Orphanet prevalence → `frequency_prognosis`, GeneReviews local DB → `references` (API 대체), OMIM MIM → `references` URL.

**Files:**
- Modify: `scripts/tools/build_gene_knowledge.py`
- Modify: `config.yaml`
- Create: `tests/test_build_gene_knowledge_phase2b.py`

- [ ] **Step 1: Write failing test**

```python
# tests/test_build_gene_knowledge_phase2b.py
import json
import pytest


def test_orphanet_prevalence_in_knowledge(tmp_path, monkeypatch):
    """Orphanet prevalence가 frequency_prognosis에 포함."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_summary",
        lambda g: {"gene": g, "description": "Test gene", "aliases": ""} if g == "CFTR" else None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_genreviews_info", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_validity_local", lambda g, **kw: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_treatment_summary", lambda g, v=None: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_variant_evidence", lambda g, v=None: [])
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_cpic_gene", lambda g: None)

    # Mock Orphanet
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_prevalence_text",
        lambda g: "Cystic fibrosis: 1-5 / 10 000 (Europe)" if g == "CFTR" else "")

    # Mock GeneReviews local DB
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_genreviews_for_gene_local",
        lambda g: None)

    # Mock OMIM
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_mim_for_gene",
        lambda g: None)

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["CFTR"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    cftr = {g["gene"]: g for g in data["genes"]}["CFTR"]
    assert "Cystic fibrosis" in cftr["frequency_prognosis"]
    assert "1-5 / 10 000" in cftr["frequency_prognosis"]


def test_genreviews_local_replaces_api(tmp_path, monkeypatch):
    """GeneReviews local DB가 API 대신 references에 추가됨."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_gene_summary",
        lambda g: {"gene": g, "entrez_id": "7157", "full_name": "tumor protein p53",
                    "summary": "Tumor suppressor.", "aliases": ""} if g == "TP53" else None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_genreviews_info", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_validity_local", lambda g, **kw: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_treatment_summary", lambda g, v=None: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_variant_evidence", lambda g, v=None: [])
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_cpic_gene", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_prevalence_text", lambda g: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_mim_for_gene", lambda g: None)

    # GeneReviews local DB returns data
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_genreviews_for_gene_local",
        lambda g: {"gene": g, "nbk_id": "NBK1311", "pmid": "20301371",
                    "title": "Li-Fraumeni Syndrome",
                    "url": "https://www.ncbi.nlm.nih.gov/books/NBK1311/"} if g == "TP53" else None)

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["TP53"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    tp53 = {g["gene"]: g for g in data["genes"]}["TP53"]
    pmids = [r.get("pmid") for r in tp53.get("references", [])]
    assert "20301371" in pmids


def test_omim_mim_in_references(tmp_path, monkeypatch):
    """OMIM MIM number가 references에 포함."""
    from scripts.tools.build_gene_knowledge import build_knowledge

    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_summary", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_gene_summary",
        lambda g: {"gene": g, "entrez_id": "7157", "full_name": "tumor protein p53",
                    "summary": "Tumor suppressor.", "aliases": ""} if g == "TP53" else None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_genreviews_info", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_gene_validity_local", lambda g, **kw: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_treatment_summary", lambda g, v=None: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_variant_evidence", lambda g, v=None: [])
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.fetch_cpic_gene", lambda g: None)
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_prevalence_text", lambda g: "")
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_genreviews_for_gene_local", lambda g: None)

    # OMIM returns MIM
    monkeypatch.setattr("scripts.tools.build_gene_knowledge.get_mim_for_gene",
        lambda g: {"gene": g, "mim_number": "191170", "entry_type": "gene",
                    "entrez_id": "7157", "ensembl_id": "", "url": "https://omim.org/entry/191170"}
        if g == "TP53" else None)

    output = str(tmp_path / "knowledge.json")
    build_knowledge(["TP53"], output)

    data = json.loads((tmp_path / "knowledge.json").read_text())
    tp53 = {g["gene"]: g for g in data["genes"]}["TP53"]
    sources = [r.get("source") for r in tp53.get("references", [])]
    assert any("OMIM" in s for s in sources)
```

- [ ] **Step 2: Update build_gene_knowledge.py**

Add imports at top:
```python
from scripts.db.query_orphanet import get_prevalence_text
from scripts.db.query_genreviews import get_genreviews_for_gene as get_genreviews_for_gene_local
from scripts.db.query_omim_mapping import get_mim_for_gene
```

In `_build_gene_entry()`, after existing logic, add:

```python
    # 8. Orphanet prevalence → frequency_prognosis
    orphanet_text = get_prevalence_text(gene)

    # 9. GeneReviews local DB (replaces API call when available)
    genreviews_local = get_genreviews_for_gene_local(gene)
    if genreviews_local and not genreviews:
        # Use local DB instead of API
        genreviews = genreviews_local

    # 10. OMIM MIM mapping
    omim = get_mim_for_gene(gene)
    if omim:
        refs.append({"pmid": "", "source": f"OMIM #{omim['mim_number']}",
                      "note": omim.get("url", "")})
```

And set `frequency_prognosis`:
```python
    "frequency_prognosis": orphanet_text or "",
```

- [ ] **Step 3: Add paths to config.yaml**

```yaml
  orphanet_db: "data/db/orphanet.sqlite3"
  genreviews_db: "data/db/genreviews.sqlite3"
  omim_mapping_db: "data/db/omim_mapping.sqlite3"
```

- [ ] **Step 4: Run tests, commit**

Run: `python -m pytest tests/test_build_gene_knowledge_phase2b.py -v`
Run: `python -m pytest tests/ --tb=short` (full suite)
Commit: `feat: integrate Orphanet/GeneReviews/OMIM into gene knowledge builder`

---

## Task 5: 실 데이터 다운로드 + 빌드 + SETUP.md

**Files:**
- Modify: `docs/SETUP.md`

- [ ] **Step 1: Download Orphanet Product 9**

```bash
curl -o data/db/en_product9_prev.xml https://www.orphadata.com/data/xml/en_product9_prev.xml
python scripts/db/build_orphanet_db.py data/db/en_product9_prev.xml
```

- [ ] **Step 2: Download GeneReviews identifier files**

```bash
curl -o data/db/GRshortname_NBKid_genesymbol_dzname.txt \
  ftp://ftp.ncbi.nlm.nih.gov/pub/GeneReviews/GRshortname_NBKid_genesymbol_dzname.txt
curl -o data/db/GRtitle_shortname_NBKid.txt \
  ftp://ftp.ncbi.nlm.nih.gov/pub/GeneReviews/GRtitle_shortname_NBKid.txt
python scripts/db/build_genreviews_db.py
```

- [ ] **Step 3: Download OMIM mim2gene.txt**

```bash
curl -o data/db/mim2gene.txt https://omim.org/static/omim/data/mim2gene.txt
python scripts/db/build_omim_mapping.py data/db/mim2gene.txt
```

- [ ] **Step 4: Rebuild gene knowledge with new sources**

```bash
python scripts/tools/build_gene_knowledge.py --genes TP53,BRCA2,KRAS,BRAF,EGFR,CFTR,PTPN11 --output data/gene_knowledge.json
```

- [ ] **Step 5: Test reports**

```bash
python scripts/orchestrate.py data/sample_vcf/rare_disease_demo.vcf \
  -o output/test_phase2b_rare.html --json --mode rare-disease \
  --hpo HP:0001250,HP:0001263 --skip-api
python scripts/orchestrate.py data/sample_vcf/codegen-Tumor_WB.mutect.passed.vep.vcf \
  -o output/test_phase2b_mutect.html --json --mode cancer --skip-api
```

- [ ] **Step 6: Update SETUP.md with download instructions**

- [ ] **Step 7: Run full suite, commit**

Run: `python -m pytest tests/ --tb=short`
Commit: `feat: Orphanet/GeneReviews/OMIM data download + build instructions`

---

## Summary

| Task | 내용 | 파일 | 테스트 |
|------|------|------|--------|
| 1 | Orphanet prevalence XML → SQLite | 3 신규 | +4 |
| 2 | GeneReviews identifiers → SQLite | 3 신규 | +4 |
| 3 | OMIM mim2gene → SQLite | 3 신규 | +4 |
| 4 | build_gene_knowledge 통합 | 수정 2, 신규 1 | +3 |
| 5 | 실 다운로드 + 빌드 + 문서 | 수정 1 | — |

**예상:** 419 + ~15 = **434+ tests**
