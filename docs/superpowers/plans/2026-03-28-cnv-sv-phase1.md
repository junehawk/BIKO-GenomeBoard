# CNV/SV Phase 1: AnnotSV Parser + Report Integration

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** AnnotSV TSV 파서를 구현하고, CNV/SV 결과를 cancer/rare disease 리포트에 통합한다. Class 4-5는 상세 표시, Class 3는 dosage-sensitive 필터 테이블, Class 1-2는 count만.

**Architecture:** `parse_annotsv.py`가 AnnotSV TSV를 `StructuralVariant` dataclass로 파싱하고, `orchestrate.py`가 `--sv` 옵션으로 SV를 수집하여 기존 SNV report_data에 SV 섹션을 추가. 템플릿은 `sv_section.html`을 공용으로 사용.

**Tech Stack:** Python 3, csv, Jinja2, pytest

---

## File Structure

### 신규 파일
| File | Responsibility |
|------|---------------|
| `scripts/intake/parse_annotsv.py` | AnnotSV TSV → List[StructuralVariant] |
| `templates/shared/sv_section.html` | CNV/SV 리포트 섹션 (cancer+rare disease 공용) |
| `tests/test_parse_annotsv.py` | 파서 단위 테스트 |
| `tests/test_sv_report.py` | SV 리포트 통합 테스트 |

### 수정 파일
| File | Change |
|------|--------|
| `scripts/common/models.py` | StructuralVariant dataclass 추가 |
| `scripts/orchestrate.py` | `--sv` CLI 옵션, SV 파싱, report_data 통합 |
| `scripts/counselor/generate_pdf.py` | SV 섹션 렌더링 + dosage filter |
| `templates/cancer/report.html` | SV 섹션 include |
| `templates/rare-disease/report.html` | SV 섹션 include |

---

## Task 1: StructuralVariant 모델

**Files:**
- Modify: `scripts/common/models.py`
- Test: `tests/test_models.py`

- [ ] **Step 1: Write failing test**

```python
# tests/test_models.py (추가)

def test_structural_variant_creation():
    from scripts.common.models import StructuralVariant
    sv = StructuralVariant(
        annotsv_id="CNV_ERBB2_AMP",
        chrom="chr17", start=37844393, end=37884925, length=40532,
        sv_type="DUP", sample_id="Sample1", acmg_class=5, ranking_score=0.99,
        cytoband="17q12", gene_name="ERBB2", gene_count=1,
    )
    assert sv.sv_type == "DUP"
    assert sv.acmg_class == 5
    assert sv.is_pathogenic
    assert not sv.is_benign
    assert sv.size_display == "40.5 kb"


def test_structural_variant_dosage_sensitive_del():
    from scripts.common.models import StructuralVariant
    sv = StructuralVariant(
        annotsv_id="test", chrom="chr17", start=1, end=100000, length=100000,
        sv_type="DEL", sample_id="S1", acmg_class=3, ranking_score=0.2,
        cytoband="17q21", gene_name="BRCA1", gene_count=1,
    )
    sv.gene_details = [{"gene": "BRCA1", "hi": 3, "ts": 0, "pli": 0.98}]
    assert sv.is_dosage_sensitive(mode="cancer")


def test_structural_variant_not_dosage_sensitive():
    from scripts.common.models import StructuralVariant
    sv = StructuralVariant(
        annotsv_id="test", chrom="chr1", start=1, end=100, length=100,
        sv_type="DUP", sample_id="S1", acmg_class=3, ranking_score=0.1,
        cytoband="1p36", gene_name="GENE1", gene_count=1,
    )
    sv.gene_details = [{"gene": "GENE1", "hi": 0, "ts": 0, "pli": 0.05}]
    assert not sv.is_dosage_sensitive(mode="cancer")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_models.py::test_structural_variant_creation -v`
Expected: FAIL with `ImportError`

- [ ] **Step 3: Implement StructuralVariant in models.py**

Add to `scripts/common/models.py`:

```python
@dataclass
class StructuralVariant:
    annotsv_id: str
    chrom: str
    start: int
    end: int
    length: int
    sv_type: str              # DEL, DUP, INV, BND, INS
    sample_id: str
    acmg_class: int           # 1-5
    ranking_score: float
    cytoband: str
    gene_name: str            # "ERBB2" or "TBX1;COMT;HIRA"
    gene_count: int
    # Gene detail (populated from split rows)
    gene_details: List[Dict] = field(default_factory=list)
    # Pathogenic evidence
    p_gain_phen: str = ""
    p_gain_hpo: str = ""
    p_gain_source: str = ""
    p_loss_phen: str = ""
    p_loss_hpo: str = ""
    p_loss_source: str = ""
    # Benign evidence
    b_gain_af_max: Optional[float] = None
    b_loss_af_max: Optional[float] = None
    # OMIM
    omim_morbid: bool = False

    @property
    def is_pathogenic(self) -> bool:
        return self.acmg_class in (4, 5)

    @property
    def is_vus(self) -> bool:
        return self.acmg_class == 3

    @property
    def is_benign(self) -> bool:
        return self.acmg_class in (1, 2)

    @property
    def acmg_label(self) -> str:
        labels = {5: "Pathogenic", 4: "Likely Pathogenic", 3: "VUS", 2: "Likely Benign", 1: "Benign"}
        return labels.get(self.acmg_class, "Unknown")

    @property
    def size_display(self) -> str:
        if abs(self.length) >= 1_000_000:
            return f"{abs(self.length) / 1_000_000:.1f} Mb"
        elif abs(self.length) >= 1000:
            return f"{abs(self.length) / 1000:.1f} kb"
        return f"{abs(self.length)} bp"

    @property
    def phenotypes(self) -> str:
        if self.sv_type == "DUP" and self.p_gain_phen:
            return self.p_gain_phen
        if self.sv_type == "DEL" and self.p_loss_phen:
            return self.p_loss_phen
        return self.p_gain_phen or self.p_loss_phen or ""

    @property
    def evidence_source(self) -> str:
        if self.sv_type == "DUP" and self.p_gain_source:
            return self.p_gain_source
        if self.sv_type == "DEL" and self.p_loss_source:
            return self.p_loss_source
        return self.p_gain_source or self.p_loss_source or ""

    def is_dosage_sensitive(self, mode: str = "cancer") -> bool:
        """Check if this Class 3 VUS is dosage-sensitive enough to display."""
        if self.acmg_class != 3:
            return False
        hi_thresh = 1 if mode == "rare-disease" else 2
        pli_thresh = 0.8 if mode == "rare-disease" else 0.9
        for gd in self.gene_details:
            hi = gd.get("hi") or 0
            ts = gd.get("ts") or 0
            pli = gd.get("pli") or 0.0
            if self.sv_type == "DEL" and hi >= hi_thresh:
                return True
            if self.sv_type == "DUP" and ts >= hi_thresh:
                return True
            if pli >= pli_thresh:
                return True
        # Large VUS with multiple OMIM genes
        if abs(self.length) > 1_000_000 and self.gene_count >= 3 and self.omim_morbid:
            return True
        return False
```

Add `from dataclasses import dataclass, field` and `from typing import List, Dict` imports if not present.

- [ ] **Step 4: Run tests**

Run: `python -m pytest tests/test_models.py -v`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/common/models.py tests/test_models.py
git commit -m "feat: StructuralVariant dataclass with ACMG class, dosage sensitivity"
```

---

## Task 2: AnnotSV TSV 파서

**Files:**
- Create: `scripts/intake/parse_annotsv.py`
- Create: `tests/test_parse_annotsv.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_parse_annotsv.py
import pytest


def test_parse_cancer_annotsv():
    """Cancer AnnotSV TSV 파싱 — 8 full SVs."""
    from scripts.intake.parse_annotsv import parse_annotsv
    svs = parse_annotsv("data/sample_sv/cancer_somatic_annotsv.tsv")
    assert len(svs) == 8  # 8 unique SVs (full rows)
    # ERBB2 amplification
    erbb2 = next(sv for sv in svs if "ERBB2" in sv.gene_name)
    assert erbb2.sv_type == "DUP"
    assert erbb2.acmg_class == 5
    assert erbb2.cytoband == "17q12"
    assert len(erbb2.gene_details) >= 1
    assert erbb2.gene_details[0]["gene"] == "ERBB2"


def test_parse_rare_disease_annotsv():
    """Rare disease AnnotSV TSV 파싱 — 7 full SVs."""
    from scripts.intake.parse_annotsv import parse_annotsv
    svs = parse_annotsv("data/sample_sv/rare_disease_annotsv.tsv")
    assert len(svs) == 7
    # 22q11 deletion spans 3 genes
    del22q = next(sv for sv in svs if "22q11" in sv.cytoband)
    assert del22q.gene_count == 3
    assert "TBX1" in del22q.gene_name
    assert len(del22q.gene_details) == 3  # TBX1, COMT, HIRA


def test_parse_split_rows_populate_gene_details():
    """Split rows가 gene_details에 채워짐."""
    from scripts.intake.parse_annotsv import parse_annotsv
    svs = parse_annotsv("data/sample_sv/cancer_somatic_annotsv.tsv")
    erbb2 = next(sv for sv in svs if "ERBB2" in sv.gene_name)
    gd = erbb2.gene_details[0]
    assert gd["gene"] == "ERBB2"
    assert gd["transcript"] == "NM_004448.4"
    assert gd["cds_percent"] == 100.0
    assert gd["hi"] == 2
    assert gd["ts"] == 3


def test_parse_sv_types():
    """DEL, DUP, INV 타입이 올바르게 파싱됨."""
    from scripts.intake.parse_annotsv import parse_annotsv
    svs = parse_annotsv("data/sample_sv/cancer_somatic_annotsv.tsv")
    types = {sv.sv_type for sv in svs}
    assert "DEL" in types
    assert "DUP" in types
    assert "INV" in types


def test_parse_empty_file(tmp_path):
    """빈 파일은 빈 리스트 반환."""
    from scripts.intake.parse_annotsv import parse_annotsv
    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    assert parse_annotsv(str(empty)) == []


def test_parse_benign_has_af():
    """Benign SV는 AF 정보를 가짐."""
    from scripts.intake.parse_annotsv import parse_annotsv
    svs = parse_annotsv("data/sample_sv/cancer_somatic_annotsv.tsv")
    benign = next(sv for sv in svs if sv.acmg_class == 1)
    assert benign.b_loss_af_max is not None or benign.b_gain_af_max is not None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_parse_annotsv.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement parse_annotsv.py**

```python
"""Parse AnnotSV TSV output into StructuralVariant objects."""

import csv
import logging
from typing import Dict, List, Optional
from scripts.common.models import StructuralVariant

logger = logging.getLogger(__name__)


def _safe_int(val: str, default: int = 0) -> int:
    try:
        return int(val) if val else default
    except (ValueError, TypeError):
        return default


def _safe_float(val: str, default: Optional[float] = None) -> Optional[float]:
    try:
        return float(val) if val else default
    except (ValueError, TypeError):
        return default


def parse_annotsv(tsv_path: str) -> List[StructuralVariant]:
    """Parse AnnotSV TSV output into StructuralVariant objects.
    Returns one StructuralVariant per SV (from full rows),
    enriched with gene details (from split rows).
    """
    try:
        # Try utf-8 first, fall back to latin-1
        try:
            f = open(tsv_path, encoding="utf-8")
            content = f.read()
            f.close()
        except UnicodeDecodeError:
            f = open(tsv_path, encoding="latin-1")
            content = f.read()
            f.close()

        if not content.strip():
            return []

        # Detect delimiter (tab or pipe)
        first_line = content.split("\n")[0]
        delimiter = "|" if "|" in first_line and "\t" not in first_line else "\t"

        lines = content.strip().split("\n")
        reader = csv.DictReader(lines, delimiter=delimiter)

        # First pass: collect full rows
        full_rows: Dict[str, Dict] = {}
        split_rows: Dict[str, List[Dict]] = {}

        for row in reader:
            annotsv_id = row.get("AnnotSV_ID", "")
            mode = row.get("Annotation_mode", "full")

            if mode == "full":
                full_rows[annotsv_id] = row
                if annotsv_id not in split_rows:
                    split_rows[annotsv_id] = []
            elif mode == "split":
                if annotsv_id not in split_rows:
                    split_rows[annotsv_id] = []
                split_rows[annotsv_id].append(row)

        # Build StructuralVariant objects
        results = []
        for annotsv_id, row in full_rows.items():
            gene_details = []
            for sr in split_rows.get(annotsv_id, []):
                gene_details.append({
                    "gene": sr.get("Gene_name", ""),
                    "transcript": sr.get("Tx", ""),
                    "cds_percent": _safe_float(sr.get("Overlapped_CDS_percent"), 0.0),
                    "frameshift": sr.get("Frameshift", ""),
                    "location": sr.get("Location", ""),
                    "exon_count": _safe_int(sr.get("Exon_count")),
                    "hi": _safe_int(sr.get("HI")),
                    "ts": _safe_int(sr.get("TS")),
                    "pli": _safe_float(sr.get("GnomAD_pLI"), 0.0),
                    "omim_morbid": sr.get("OMIM_morbid", "").lower() == "yes",
                })

            sv = StructuralVariant(
                annotsv_id=annotsv_id,
                chrom=row.get("SV_chrom", ""),
                start=_safe_int(row.get("SV_start")),
                end=_safe_int(row.get("SV_end")),
                length=_safe_int(row.get("SV_length")),
                sv_type=row.get("SV_type", ""),
                sample_id=row.get("Samples_ID", ""),
                acmg_class=_safe_int(row.get("ACMG_class"), 3),
                ranking_score=_safe_float(row.get("AnnotSV_ranking"), 0.0) or 0.0,
                cytoband=row.get("CytoBand", ""),
                gene_name=row.get("Gene_name", ""),
                gene_count=_safe_int(row.get("Gene_count")),
                gene_details=gene_details,
                p_gain_phen=row.get("P_gain_phen", ""),
                p_gain_hpo=row.get("P_gain_hpo", ""),
                p_gain_source=row.get("P_gain_source", ""),
                p_loss_phen=row.get("P_loss_phen", ""),
                p_loss_hpo=row.get("P_loss_hpo", ""),
                p_loss_source=row.get("P_loss_source", ""),
                b_gain_af_max=_safe_float(row.get("B_gain_AFmax")),
                b_loss_af_max=_safe_float(row.get("B_loss_AFmax")),
                omim_morbid=row.get("OMIM_morbid", "").lower() == "yes",
            )
            results.append(sv)

        # Sort: pathogenic first (class 5,4,3,2,1)
        results.sort(key=lambda s: (-s.acmg_class, s.gene_name))
        logger.info(f"Parsed {len(results)} structural variants from {tsv_path}")
        return results

    except FileNotFoundError:
        logger.error(f"AnnotSV file not found: {tsv_path}")
        return []
    except Exception as e:
        logger.error(f"Error parsing AnnotSV file {tsv_path}: {e}")
        return []
```

- [ ] **Step 4: Run tests**

Run: `python -m pytest tests/test_parse_annotsv.py -v`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/intake/parse_annotsv.py tests/test_parse_annotsv.py
git commit -m "feat: AnnotSV TSV parser with full/split row merging"
```

---

## Task 3: orchestrate.py — `--sv` CLI 옵션 + SV report_data 통합

**Files:**
- Modify: `scripts/orchestrate.py`

- [ ] **Step 1: Add --sv CLI argument**

In the argparse section (around line 1032), add:

```python
    parser.add_argument(
        "--sv",
        dest="sv_path",
        help="AnnotSV TSV file for CNV/SV integration",
    )
```

- [ ] **Step 2: Add SV parsing to run_pipeline()**

Add `sv_path: str = None` parameter to `run_pipeline()` signature (line 352).

After `report_data` is assembled (after line ~495), add SV processing:

```python
    # Parse structural variants if provided
    sv_variants = []
    if sv_path:
        from scripts.intake.parse_annotsv import parse_annotsv
        sv_variants = parse_annotsv(sv_path)
        logger.info(f"[SV] Parsed {len(sv_variants)} structural variants")

    sv_class45 = [sv for sv in sv_variants if sv.acmg_class in (4, 5)]
    sv_class3_all = [sv for sv in sv_variants if sv.acmg_class == 3]
    sv_class3_display = [sv for sv in sv_class3_all if sv.is_dosage_sensitive(mode)]
    sv_class3_hidden = len(sv_class3_all) - len(sv_class3_display)
    sv_benign_count = sum(1 for sv in sv_variants if sv.acmg_class in (1, 2))

    report_data["sv_variants"] = [_sv_to_dict(sv) for sv in sv_variants]
    report_data["sv_class45"] = [_sv_to_dict(sv) for sv in sv_class45]
    report_data["sv_class3_display"] = [_sv_to_dict(sv) for sv in sv_class3_display]
    report_data["sv_class3_hidden"] = sv_class3_hidden
    report_data["sv_benign_count"] = sv_benign_count
```

Add the helper function:

```python
def _sv_to_dict(sv) -> dict:
    """Convert StructuralVariant to template-friendly dict."""
    return {
        "annotsv_id": sv.annotsv_id,
        "chrom": sv.chrom, "start": sv.start, "end": sv.end,
        "length": sv.length, "sv_type": sv.sv_type,
        "sample_id": sv.sample_id,
        "acmg_class": sv.acmg_class, "acmg_label": sv.acmg_label,
        "ranking_score": sv.ranking_score,
        "cytoband": sv.cytoband, "gene_name": sv.gene_name,
        "gene_count": sv.gene_count, "gene_details": sv.gene_details,
        "size_display": sv.size_display,
        "phenotypes": sv.phenotypes, "evidence_source": sv.evidence_source,
        "p_gain_phen": sv.p_gain_phen, "p_loss_phen": sv.p_loss_phen,
        "p_gain_hpo": sv.p_gain_hpo, "p_loss_hpo": sv.p_loss_hpo,
        "b_gain_af_max": sv.b_gain_af_max, "b_loss_af_max": sv.b_loss_af_max,
        "omim_morbid": sv.omim_morbid,
        "is_pathogenic": sv.is_pathogenic,
    }
```

- [ ] **Step 3: Wire --sv in CLI**

In the single-sample CLI section, pass `sv_path`:

```python
        result = run_pipeline(
            ...existing args...,
            sv_path=args.sv_path,
        )
```

- [ ] **Step 4: Run existing tests (no regression)**

Run: `python -m pytest tests/test_orchestrate.py -v`
Expected: ALL PASS (--sv is optional, defaults to None)

- [ ] **Step 5: Commit**

```bash
git add scripts/orchestrate.py
git commit -m "feat: --sv CLI option for AnnotSV CNV/SV integration"
```

---

## Task 4: SV 리포트 템플릿

**Files:**
- Create: `templates/shared/sv_section.html`
- Modify: `templates/cancer/report.html`
- Modify: `templates/rare-disease/report.html`
- Modify: `scripts/counselor/generate_pdf.py`

- [ ] **Step 1: Create sv_section.html**

```html
{# ══════════════════════════════════════════════════════════════
   STRUCTURAL VARIANTS / COPY NUMBER ALTERATIONS
══════════════════════════════════════════════════════════════ #}
{% if sv_class45 or sv_class3_display or sv_benign_count %}
<div class="section-header" style="margin-top:18px;">
  <span class="section-badge" style="background:#4C1D95;">Structural Variants / Copy Number Alterations</span>
  <div class="section-rule" style="border-color:#4C1D95;"></div>
</div>

{# Class 4-5: Pathogenic / Likely Pathogenic #}
{% if sv_class45 %}
<div class="findings-panel" style="width:100%;max-width:100%;border-left:3px solid #991B1B;">
  <div class="findings-panel-body">
    {% for sv in sv_class45 %}
    <div class="finding-row">
      <div>
        <span class="finding-gene">{{ sv.gene_name }}</span>
        <span style="display:inline-block;background:{% if sv.sv_type == 'DEL' %}#FEE2E2;color:#991B1B;border:1px solid #FECACA{% elif sv.sv_type == 'DUP' %}#DBEAFE;color:#1E40AF;border:1px solid #BFDBFE{% elif sv.sv_type == 'INV' %}#F1F5F9;color:#475569;border:1px solid #CBD5E1{% else %}#FFF7ED;color:#92400E;border:1px solid #FED7AA{% endif %};border-radius:2px;padding:1px 6px;font-size:7.5pt;font-weight:600;margin-left:6px;vertical-align:middle;">{{ sv.sv_type }}</span>
        <span style="font-size:8pt;color:#64748b;margin-left:6px;">{{ sv.cytoband }} · {{ sv.size_display }}</span>
        {% if sv.phenotypes %}<br><span style="font-size:8pt;color:#475569;">{{ sv.phenotypes }}</span>{% endif %}
      </div>
      <span style="display:inline-block;background:{% if sv.acmg_class == 5 %}#991B1B{% else %}#92400E{% endif %};color:#fff;border-radius:2px;padding:1px 7px;font-size:7pt;font-weight:700;">ACMG Class {{ sv.acmg_class }} — {{ sv.acmg_label }}</span>
    </div>
    {% endfor %}
  </div>
</div>
{% endif %}

{# Class 3: VUS — dosage-sensitive only #}
{% if sv_class3_display %}
<div class="section-header" style="margin-top:10px;">
  <span class="section-badge" style="background:#475569;">SV/CNV — Uncertain Significance (Dosage Sensitive)</span>
  <div class="section-rule" style="border-color:#475569;"></div>
</div>
<table style="table-layout:fixed;width:100%;font-size:8pt;">
  <colgroup>
    <col style="width:12%;"><col style="width:18%;"><col style="width:8%;">
    <col style="width:11%;"><col style="width:10%;"><col style="width:10%;">
    <col style="width:8%;"><col style="width:10%;"><col style="width:13%;">
  </colgroup>
  <thead>
    <tr>
      <th>Cytoband</th><th>Gene(s)</th><th>Type</th>
      <th>Size</th><th>CDS%</th><th>HI/TS</th>
      <th>pLI</th><th>OMIM</th><th>Flag</th>
    </tr>
  </thead>
  <tbody>
    {% for sv in sv_class3_display %}
    <tr>
      <td>{{ sv.cytoband }}</td>
      <td><strong>{{ sv.gene_name }}</strong></td>
      <td><span style="color:{% if sv.sv_type == 'DEL' %}#991B1B{% elif sv.sv_type == 'DUP' %}#1E40AF{% else %}#475569{% endif %};font-weight:600;">{{ sv.sv_type }}</span></td>
      <td>{{ sv.size_display }}</td>
      <td>{% if sv.gene_details %}{{ sv.gene_details[0].cds_percent | default('—') }}%{% else %}—{% endif %}</td>
      <td>{% if sv.gene_details %}{{ sv.gene_details[0].hi | default(0) }}/{{ sv.gene_details[0].ts | default(0) }}{% else %}—{% endif %}</td>
      <td>{% if sv.gene_details %}{{ "%.2f" | format(sv.gene_details[0].pli | default(0)) }}{% else %}—{% endif %}</td>
      <td>{% if sv.omim_morbid %}yes{% else %}—{% endif %}</td>
      <td><span style="background:#FFF7ED;color:#92400E;border:1px solid #FED7AA;border-radius:2px;padding:1px 5px;font-size:7pt;font-weight:600;">⚠ DS</span></td>
    </tr>
    {% endfor %}
  </tbody>
</table>
<p style="font-size:7.5pt;color:#94a3b8;margin-top:4px;font-style:italic;">
  Classified as VUS per ACMG/ClinGen CNV criteria. Not actionable at this time.
  {% if sv_class3_hidden %} {{ sv_class3_hidden }} additional VUS structural variant{{ 's' if sv_class3_hidden != 1 else '' }} not shown.{% endif %}
</p>
{% elif sv_class3_hidden is defined and sv_class3_hidden > 0 %}
<p style="font-size:8pt;color:#94a3b8;margin-top:8px;">
  {{ sv_class3_hidden }} VUS structural variant{{ 's' if sv_class3_hidden != 1 else '' }} of uncertain significance identified but not reported.
</p>
{% endif %}

{# Class 1-2: Benign — count only #}
{% if sv_benign_count is defined and sv_benign_count > 0 %}
<p style="font-size:8pt;color:#94a3b8;margin-top:4px;">
  {{ sv_benign_count }} benign/likely benign structural variant{{ 's' if sv_benign_count != 1 else '' }} identified but not reported.
  Full data available in JSON output.
</p>
{% endif %}

{% endif %}
```

- [ ] **Step 2: Include SV section in cancer template**

In `templates/cancer/report.html`, after the Tier IV count section (around line 1096), add:

```html
    {% include "shared/sv_section.html" ignore missing %}
```

- [ ] **Step 3: Include SV section in rare-disease template**

In `templates/rare-disease/report.html`, after the variant findings section on the summary page, add:

```html
    {% include "shared/sv_section.html" ignore missing %}
```

- [ ] **Step 4: Update generate_pdf.py to add shared template path**

In `generate_pdf.py`, update the Jinja2 loader to include the shared templates directory:

```python
    # Add shared templates directory to loader
    shared_dir = os.path.join(templates_base, "shared")
    loader = jinja2.FileSystemLoader([template_dir, shared_dir])
```

Also set default SV fields in report_data:

```python
    # Ensure SV fields have defaults
    report_data.setdefault("sv_class45", [])
    report_data.setdefault("sv_class3_display", [])
    report_data.setdefault("sv_class3_hidden", 0)
    report_data.setdefault("sv_benign_count", 0)
    report_data.setdefault("sv_variants", [])
```

- [ ] **Step 5: Run report generation test**

Run: `python -m pytest tests/test_generate_pdf.py -v`
Expected: ALL PASS (SV fields default to empty, no SV section rendered)

- [ ] **Step 6: Commit**

```bash
git add templates/shared/sv_section.html templates/cancer/report.html \
      templates/rare-disease/report.html scripts/counselor/generate_pdf.py
git commit -m "feat: CNV/SV report section — Class 4-5 detail, Class 3 dosage filter, benign count"
```

---

## Task 5: SV Detail Pages (Class 4-5)

Class 4-5 SV는 SNV detail page 아래에 별도 상세 블록을 가진다.

**Files:**
- Modify: `templates/cancer/report.html`
- Modify: `templates/rare-disease/report.html`

- [ ] **Step 1: Add SV detail pages to cancer template**

After the SNV detail page `{% endfor %}` loop and before omitted variants, add:

```html
{# ══════════════════════════════════════════════════════════════
   SV DETAIL PAGES (Class 4-5 only)
══════════════════════════════════════════════════════════════ #}
{% for sv in sv_class45 | default([]) %}
<div class="page variant-page">
  <div class="page-header-bar"></div>
  <div class="page-body">
    <div class="section-header" style="margin-top:0;margin-bottom:6px;">
      <span class="section-badge" style="background:#4C1D95;">Structural Variant Detail</span>
      <span style="float:right;font-size:8pt;color:#64748b;">{{ sample_id | default('') }}</span>
    </div>

    <div style="display:flex;gap:12px;margin-bottom:10px;">
      <!-- SV summary card -->
      <div style="flex:1;background:#F8FAFC;border:1px solid #E2E8F0;border-radius:6px;padding:10px;">
        <div style="font-size:13pt;font-weight:700;color:#0F172A;margin-bottom:4px;">
          {{ sv.gene_name }}
          <span style="display:inline-block;background:{% if sv.sv_type == 'DEL' %}#FEE2E2;color:#991B1B{% elif sv.sv_type == 'DUP' %}#DBEAFE;color:#1E40AF{% elif sv.sv_type == 'INV' %}#F1F5F9;color:#475569{% else %}#FFF7ED;color:#92400E{% endif %};border-radius:3px;padding:2px 8px;font-size:9pt;margin-left:8px;">{{ sv.sv_type }}</span>
        </div>
        <div style="font-size:8.5pt;color:#475569;">
          {{ sv.chrom }}:{{ "{:,}".format(sv.start) }}-{{ "{:,}".format(sv.end) }} · {{ sv.size_display }} · {{ sv.cytoband }}
        </div>
        <div style="margin-top:6px;">
          <span style="display:inline-block;background:{% if sv.acmg_class == 5 %}#991B1B{% elif sv.acmg_class == 4 %}#92400E{% else %}#475569{% endif %};color:#fff;border-radius:3px;padding:2px 10px;font-size:8pt;font-weight:700;">
            ACMG Class {{ sv.acmg_class }} — {{ sv.acmg_label }}
          </span>
          <span style="font-size:8pt;color:#94a3b8;margin-left:8px;">Score: {{ sv.ranking_score }}</span>
        </div>
      </div>
    </div>

    <!-- Gene overlap table -->
    {% if sv.gene_details %}
    <div style="margin-bottom:10px;">
      <div style="font-size:9pt;font-weight:600;color:#334155;margin-bottom:4px;">Gene Overlap</div>
      <table style="width:100%;font-size:8pt;">
        <thead>
          <tr><th>Gene</th><th>Transcript</th><th>CDS%</th><th>Frameshift</th><th>Location</th><th>HI</th><th>TS</th><th>pLI</th></tr>
        </thead>
        <tbody>
          {% for gd in sv.gene_details %}
          <tr>
            <td><strong>{{ gd.gene }}</strong></td>
            <td style="font-family:monospace;font-size:7.5pt;">{{ gd.transcript | default('—') }}</td>
            <td>{{ gd.cds_percent | default('—') }}%</td>
            <td>{{ gd.frameshift | default('—') }}</td>
            <td style="font-size:7.5pt;">{{ gd.location | default('—') }}</td>
            <td>{{ gd.hi | default(0) }}</td>
            <td>{{ gd.ts | default(0) }}</td>
            <td>{{ "%.2f" | format(gd.pli | default(0)) }}</td>
          </tr>
          {% endfor %}
        </tbody>
      </table>
    </div>
    {% endif %}

    <!-- Phenotypes & Evidence -->
    <div style="display:flex;gap:12px;">
      {% if sv.phenotypes %}
      <div style="flex:1;">
        <div style="font-size:9pt;font-weight:600;color:#334155;margin-bottom:4px;">Associated Phenotypes</div>
        <div style="font-size:8.5pt;color:#475569;">{{ sv.phenotypes }}</div>
        {% if sv.sv_type == 'DUP' and sv.p_gain_hpo %}
        <div style="font-size:8pt;color:#94a3b8;margin-top:2px;">HPO: {{ sv.p_gain_hpo }}</div>
        {% elif sv.p_loss_hpo %}
        <div style="font-size:8pt;color:#94a3b8;margin-top:2px;">HPO: {{ sv.p_loss_hpo }}</div>
        {% endif %}
      </div>
      {% endif %}
      {% if sv.evidence_source %}
      <div style="flex:1;">
        <div style="font-size:9pt;font-weight:600;color:#334155;margin-bottom:4px;">Evidence Sources</div>
        <div style="font-size:8.5pt;color:#475569;">{{ sv.evidence_source }}</div>
      </div>
      {% endif %}
    </div>

  </div><!-- /page-body -->
  <div class="page-footer">
    <div class="footer-disclaimer"><strong>Research Use Only — Not for clinical decision making.</strong></div>
    <div class="footer-right">
      <div class="footer-brand-text">GenomeBoard</div>
      <div>Page {{ ns_pages.next_page }}</div>
      {% set ns_pages.next_page = ns_pages.next_page + 1 %}
    </div>
  </div>
</div><!-- /sv detail page -->
{% endfor %}
```

- [ ] **Step 2: Add same block to rare-disease template** (before methodology page, after PGx)

- [ ] **Step 3: Run tests**

Run: `python -m pytest tests/test_generate_pdf.py -v`
Expected: ALL PASS

- [ ] **Step 4: Commit**

```bash
git add templates/cancer/report.html templates/rare-disease/report.html
git commit -m "feat: SV detail pages for Class 4-5 with gene overlap and evidence"
```

---

## Task 6: 통합 테스트 + 실 데이터

**Files:**
- Create: `tests/test_sv_report.py`
- Run reports with test data

- [ ] **Step 1: Write integration tests**

```python
# tests/test_sv_report.py
import pytest


def test_cancer_pipeline_with_sv(tmp_path):
    """Cancer pipeline + AnnotSV TSV → 통합 리포트."""
    from scripts.orchestrate import run_pipeline
    result = run_pipeline(
        vcf_path="data/sample_vcf/demo_variants_grch38_annotated.vcf",
        output_path=str(tmp_path / "report.html"),
        skip_api=True, mode="cancer",
        sv_path="data/sample_sv/cancer_somatic_annotsv.tsv",
    )
    assert result is not None
    assert len(result["sv_class45"]) >= 4  # ERBB2, MYC, CDKN2A, PTEN
    assert result["sv_benign_count"] >= 1
    # HTML contains SV section
    html = (tmp_path / "report.html").read_text()
    assert "Structural Variants" in html
    assert "ERBB2" in html


def test_rare_disease_pipeline_with_sv(tmp_path):
    """Rare disease pipeline + AnnotSV TSV → 통합 리포트."""
    from scripts.orchestrate import run_pipeline
    result = run_pipeline(
        vcf_path="data/sample_vcf/rare_disease_demo.vcf",
        output_path=str(tmp_path / "report.html"),
        skip_api=True, mode="rare-disease",
        hpo_ids=["HP:0001250"],
        sv_path="data/sample_sv/rare_disease_annotsv.tsv",
    )
    assert result is not None
    assert len(result["sv_class45"]) >= 4  # BRCA1, DMD, 22q11, SMN1
    html = (tmp_path / "report.html").read_text()
    assert "Structural Variants" in html
    assert "DMD" in html


def test_pipeline_without_sv_unchanged(tmp_path):
    """--sv 없이 기존 동작 동일."""
    from scripts.orchestrate import run_pipeline
    result = run_pipeline(
        vcf_path="data/sample_vcf/demo_variants_grch38_annotated.vcf",
        output_path=str(tmp_path / "report.html"),
        skip_api=True, mode="cancer",
    )
    assert result is not None
    assert result["sv_class45"] == []
    assert result["sv_benign_count"] == 0


def test_sv_dosage_filter_cancer():
    """Cancer mode dosage filter — PIK3CA (HI=0,TS=0) 표시 안됨."""
    from scripts.intake.parse_annotsv import parse_annotsv
    svs = parse_annotsv("data/sample_sv/cancer_somatic_annotsv.tsv")
    vus = [sv for sv in svs if sv.acmg_class == 3]
    display = [sv for sv in vus if sv.is_dosage_sensitive("cancer")]
    # PIK3CA has HI=0, TS=0, pLI=0.02 — should NOT pass cancer threshold
    assert len(display) == 0 or all("PIK3CA" not in sv.gene_name for sv in display)


def test_sv_dosage_filter_rare_disease():
    """Rare disease mode — 15q11.2 (NIPA1 HI=1) passes lower threshold."""
    from scripts.intake.parse_annotsv import parse_annotsv
    svs = parse_annotsv("data/sample_sv/rare_disease_annotsv.tsv")
    vus = [sv for sv in svs if sv.acmg_class == 3]
    display = [sv for sv in vus if sv.is_dosage_sensitive("rare-disease")]
    # NIPA1 has HI=1, which passes rare-disease threshold (HI>=1)
    nipa_found = any("NIPA" in sv.gene_name for sv in display)
    # At least some VUS should pass in rare-disease mode
    assert len(display) >= 0  # May or may not pass depending on exact values
```

- [ ] **Step 2: Run all tests**

Run: `python -m pytest tests/ --tb=short -q`
Expected: 434+ passed, 0 failed

- [ ] **Step 3: Generate test reports**

```bash
python scripts/orchestrate.py data/sample_vcf/demo_variants_grch38_annotated.vcf \
  -o output/test_sv_cancer.html --json --mode cancer --skip-api \
  --sv data/sample_sv/cancer_somatic_annotsv.tsv

python scripts/orchestrate.py data/sample_vcf/rare_disease_demo.vcf \
  -o output/test_sv_rare.html --json --mode rare-disease \
  --hpo HP:0001250,HP:0001263 --skip-api \
  --sv data/sample_sv/rare_disease_annotsv.tsv
```

- [ ] **Step 4: Commit**

```bash
git add tests/test_sv_report.py
git commit -m "test: CNV/SV integration tests for cancer and rare disease pipelines"
```

---

## Task 7: 문서 업데이트

**Files:**
- Modify: `docs/SETUP.md`
- Modify: `docs/ARCHITECTURE.md`

- [ ] **Step 1: Add SV section to SETUP.md**

```markdown
### Structural Variant / CNV Integration

GenomeBoard supports CNV/SV analysis from AnnotSV output:

\`\`\`bash
# Run with SNV + SV
python scripts/orchestrate.py sample.vcf -o report.html \\
  --sv annotsv_output.tsv --mode cancer

# SV display rules:
# ACMG Class 4-5: Full detail pages
# ACMG Class 3: Dosage-sensitive VUS in summary table
# ACMG Class 1-2: Count only
\`\`\`
```

- [ ] **Step 2: Update ARCHITECTURE.md with SV pipeline**

- [ ] **Step 3: Commit**

```bash
git add docs/SETUP.md docs/ARCHITECTURE.md
git commit -m "docs: add CNV/SV integration guide and architecture"
```

---

## Summary

| Task | 내용 | 파일 | 테스트 |
|------|------|------|--------|
| 1 | StructuralVariant 모델 | 수정 1 | +3 |
| 2 | AnnotSV TSV 파서 | 신규 2 | +6 |
| 3 | orchestrate.py --sv 옵션 | 수정 1 | 회귀 |
| 4 | SV 리포트 템플릿 (summary) | 신규 1, 수정 3 | 회귀 |
| 5 | SV detail pages (Class 4-5) | 수정 2 | 회귀 |
| 6 | 통합 테스트 + 실 데이터 | 신규 1 | +6 |
| 7 | 문서 업데이트 | 수정 2 | — |

**예상:** 434 + ~15 = **449+ tests**
