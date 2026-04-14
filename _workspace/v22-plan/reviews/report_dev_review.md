# Report-dev Review — v2.2 Plan

**Reviewer:** report-dev
**Date:** 2026-04-14
**Plan:** `docs/superpowers/plans/2026-04-14-ai-board-v2.2.md`
**Focus:** render/HTML/UX, Jinja templates, A4 print, backward compat with `rerender_report.py`.

---

## Verdict

**Approve with required changes (2) + strong recommendations (5).** The plan is architecturally correct from the render layer's perspective, but two render-touching tasks (B3 references, B4 patient header) under-specify the data contract between LLM / orchestrator / renderer and will break `rerender_report.py` if implemented literally. Both are recoverable with small spec additions.

---

## Section-by-section critique

### A1 — Curated Treatments Module (render implications)

Not a render task on its face, but it creates a new `report_data["_curated_treatments"]` shape that the renderer must consume. The plan says "render the curated treatment row set" under A2 and B3, but **never specifies who hydrates `curated_id` → display row**. Two options:

- **(a) `runner.py` hydrates** — before calling render, resolve each `treatment_options[i].curated_id` into a full row and attach `treatment_options[i]._hydrated = CuratedTreatment`. Render just reads the attached row.
- **(b) `render.py` hydrates** — render receives both `opinion.treatment_options` and `report_data["_curated_treatments"]`, and joins on `curated_id` at render time.

**Recommendation:** option (a). `render.py` currently takes only `opinion` as input (see `render_board_opinion_html(opinion, language)` — it does **not** have access to `report_data`). Plumbing `report_data` into render.py is a non-trivial signature change that also breaks `rerender_report.py`'s contract (`rerender` already reconstructs `opinion` from a dict and calls `render_board_opinion_html(rebuilt, language)`). Hydrating in `runner.py` keeps render.py's signature stable and keeps the rerender path working for both old and new JSON dumps.

### A2 — Narrate-only prompts (render implications)

- The acceptance criterion "No curated treatment evidence for this variant set — see domain specialist notes" needs a **KO equivalent** for the rare-disease / bilingual path. Add: `"이 변이군에 대해 큐레이션된 치료 근거가 없습니다 — 전문의 소견 참조"`. Plan does not mention this.
- The `treatment_options[].curated_id` field must be retained in the serialised `CancerBoardOpinion` JSON dump (orchestrate `clinical_board` key) so that `rerender_report.py` can re-hydrate without re-running the curator. Plan does not say this explicitly.

### A3 / A4 — ACMG engine + ClinVar override

Rendering the `clinvar_override_reason` string is the only render touch. Plan says "suggest 'under the variant card with the override applied'". **Problem:** `render.py` (clinical_board module) does not render variant cards; variant cards live in `templates/cancer/report.html` (lines ~355–370, the `.variant-block` + `.variant-detail-header` structure). So the override reason must be piped into `report_data["variants"][i].clinvar_override_reason` and rendered from the Jinja template, not from `clinical_board/render.py`. The plan file-map does list `scripts/clinical_board/render.py` as touched for this — that line should be corrected to **`templates/cancer/report.html`** (the variant detail page), not `clinical_board/render.py`.

Recommended render location: inside the `.variant-detail-header` block on the variant-detail page, as an amber `.override-notice` strip with the full reason text and a monospace footnote linking to PMID 30224644. Needs `page-break-inside: avoid` because the reason text can exceed 1 line.

Add KO copy: `"ClinVar 오버라이드: ..."` prefix.

### B1 / B2 — Selector tightening + MMR carve-out

B1 says "Modify `scripts/clinical_board/render.py` (Tier III VUS table renderer must also apply the consequence gate if it builds its list directly from `report_data[\"variants\"]`)". I audited `clinical_board/render.py`: **there is no Tier III VUS table rendered from `render.py`.** The Tier III VUS table is emitted from `templates/cancer/report.html` Jinja (search pattern for Tier III VUS in the variant detail pages). This is a file-map inaccuracy — the consequence gate must be applied in **variant_selector.py upstream**, and no render changes are needed. Please correct the file map.

B2 is pure selector work — no render impact. ✓

### B3 — PMID references + superscript rendering (CRITICAL)

This is the weakest-spec'd task in the plan and the one most likely to regress the showcase if implemented literally. Problems:

**Problem 1 — LLM marker protocol is unspecified.** The plan says "render.py emits PMIDs as `<sup>[pmid]</sup>` inline at the point of the claim". But the claim text is free-form LLM prose. For the renderer to know *where* to place the superscript, the LLM must embed a marker in the narrative (e.g., `[ref:30224644]` or a curated_id like `[ref:a1b2c3d4]`). The plan does not define the marker syntax, nor say whether the marker is embedded in every free-text field (`therapeutic_implications`, `therapeutic_evidence`, `treatment_options[].resistance_notes`, `actionable_findings[]`, `clinical_actions[]`, `monitoring_plan[]`, `dissenting_opinions[]`, agent `findings[].finding`, agent `findings[].evidence`). Without this spec the renderer cannot function.

**Recommendation:** define a single canonical marker syntax — `[pmid:30224644]` — and instruct every free-text-emitting agent prompt (Chair + all domain specialists) to use that exact token. Render pass does a single regex substitution:

```python
_PMID_MARKER = re.compile(r"\[pmid:(\d+)\]")
def _superscriptify(text, ref_index):
    def _sub(m):
        pmid = m.group(1)
        n = ref_index.setdefault(pmid, len(ref_index) + 1)
        return f'<sup><a href="#ref-{pmid}">[{n}]</a></sup>'
    return _PMID_MARKER.sub(_sub, text)
```

An ordered `ref_index: dict[pmid → n]` accumulator threads through the whole opinion render so the numbered list at the bottom is globally consistent.

**Problem 2 — static HTML constraint.** Plan says "hover tooltips (no JS in current reports)" — correct, we must not add JS. But `<sup><a href="#ref-30224644">[1]</a></sup>` works fine in WeasyPrint PDF mode **only if** the target `<div id="ref-30224644">...</div>` is inside the same rendered document. It is not interactive in print (clicks don't matter on paper), so **keep it server-side-anchored HTML only; no `title=` tooltips either** because those only appear on screen, not in print.

**Problem 3 — `AgentOpinion.references` shape change breaks `rerender_report.py`.** The current `AgentOpinion.references: List[str]` is serialised in cached orchestrate JSON dumps as a list of strings. The plan changes it to `List[Reference]` where `Reference` is a dataclass. `rerender_report.py::_filter_dataclass_kwargs` will happily pass a `list[str]` through to the `AgentOpinion(references=...)` constructor — and then `_render_agent_opinions_section` does `", ".join(agent_op.references[:3])` expecting strings. After the type change, that line is `", ".join(Reference_instance, Reference_instance, ...)` which TypeErrors. **This will break re-rendering of every existing v2 showcase JSON**, including the three committed under `docs/showcase/`.

**Required fix:** add a coercion step in `rerender_report.py::_reconstruct_board_opinion` that upgrades `list[str]` → `list[Reference]` with `pmid=string, quote="", context="(legacy, no context)", source="legacy"`. And `_render_agent_opinions_section` must be rewritten to iterate `ref.pmid` instead of joining raw strings. Both changes must land in the same commit as the model change, not later. Plan does not flag this.

**Problem 4 — ordered reference numbering across cards.** Currently `_render_cancer_opinion` renders sections in sequence (Therapeutic Implications → Treatment Options → Actionable/Clinical → Monitoring/Dissent → Agent opinions). If `[pmid:X]` appears in multiple sections, the numbered list needs a single shared counter. The accumulator pattern above handles this, but the plan should explicitly require it as an acceptance criterion ("reference numbers are globally unique and appear in first-citation order").

**Problem 5 — references list placement.** Plan says "a numbered references list at the end of the AI Board section". The current `_render_cancer_opinion` return is a string inserted into the Jinja template at a specific page. Appending a reference list there is fine, but it must carry `page-break-inside: avoid` because a split reference list across an A4 boundary looks terrible, and it must be able to overflow naturally onto a new page if >12 entries. Please add this as a CSS acceptance criterion.

**Revised effort estimate:** plan says M (1.5–2 days). I estimate **M/L (2.5–3 days)** once marker-protocol spec + rerender_report coercion + KO/EN pairs + A4 overflow handling are included. Still feasible, just larger than stated.

### B4 — Patient demographics header (CRITICAL)

Clinical-advisor's open question #1: *"Is the existing `templates/cancer/report.html` layout flexible enough for a running patient header on every page, or does adding a header require a Jinja layout overhaul?"*

**Answer: partially flexible, but 'running header on every page' is not trivial.** Here's the ground truth:

1. The current `@page { size: A4; margin: 0 }` rule (line 850–853) gives no WeasyPrint margin boxes — there is no CSS `@top-center` / running-element setup. Every "header" on every page is a hand-authored `<div class="page-header-bar"></div>` at the top of each `<div class="page">` block.
2. Pages that currently exist and would each need a patient strip:
   - page 1: `#page-summary` (masthead line ~908)
   - Classification pages (Pathogenic/LP, VUS, Drug response, Risk, Benign) each wrap in `<div class="page">`
   - Variant detail pages (line 1680 area, `variant-page` class, already have a right-float `Sample ID` stub at line 1687)
   - SV detail pages (same pattern)
   - PGx pages
   - AI Clinical Board page (wraps `{{ clinical_board_html }}`)
   - Methodology page
3. That's ~7–10 distinct page divs. Manually inserting a patient strip into each is tedious but mechanical. **It is not a 1-day task.**

**Preferred approach (two alternatives):**

- **(a) Jinja macro + page header_bar upgrade.** Define `{% macro patient_strip(patient, sample_id, date, lang) %}` in `templates/cancer/report.html`, and replace every `<div class="page-header-bar"></div>` with `{{ patient_strip(patient, sample_id, date, language) }}`. Centralises markup, bounds the edit to one macro definition + N call sites. This is the cleanest and safest approach.
- **(b) WeasyPrint running element.** Use `@page { @top-center { content: element(patient-running); } }` + `.patient-strip { position: running(patient-running); }`. Prettier in principle but requires moving away from `@page { margin: 0 }` which will alter the left/right padding of *every* existing page — a high-risk CSS refactor. **Do not pick this for v2.2.**

Plan effort estimate S (1 day) is optimistic. Realistic: **S/M (1.5 days)** for approach (a), including a test that every `<div class="page">` in both cancer and rare-disease templates contains the macro call.

**Mandatory fields vs optional (per task description):**
- **Mandatory for display when any patient data is passed:** sample_id, report_date (already exist today). These already render with `{{ sample_id }}`.
- **Recommended mandatory when `patient` dict is present:** name, DOB, sex. If any of the three is missing, show `—` rather than omitting.
- **Optional (render only when present):** mrn, ordering_physician, facility, specimen_site, specimen_id, collected_date.
- **When patient dict is None / absent:** the entire patient block collapses to just `Sample ID: {{ sample_id }}` — the current behaviour. Plan acknowledges this ✓.

**PHI safety — cross-check with rerender_report.py.** If patient is passed via CLI and persisted into `report_data["patient"]`, then `orchestrate.py` dumps report_data to JSON, that JSON contains PHI. The plan says `reporting.persist_patient_metadata: false` by default — good. But it must **also** gate the JSON-dump step, not just the SQLite persistence step, because `docs/showcase/*.json` is committed to git. Add: "When `persist_patient_metadata=false`, the `patient` key is stripped from any `_workspace/`-destined JSON and from any JSON dump file located under `docs/showcase/`. It is only preserved in ephemeral report_data objects in-memory for the template render pass."

**rerender_report.py backward compat.** Plan's claim "When absent, the current Sample-ID-only header renders unchanged" needs to be enforced by `{% if patient %} ... {% else %} ...legacy masthead... {% endif %}` in the template, not by runtime Python. Old JSON dumps in `docs/showcase/` have no `patient` key, so the masthead path must still render via `sample_id`. The plan's B4 section says this but I want it in the acceptance criteria explicitly: *"The three existing showcase JSONs in `docs/showcase/` must re-render identically to their current committed HTML (modulo any other v2.2 changes) when passed through `rerender_report.py`."* Add this as a CI regression.

### C1 / C2 — Documentation / showcase regen

No render concerns. ✓

---

## Cross-cutting concerns

### A4 print fidelity

Current print CSS (lines 824–854):
- `@page { size: A4; margin: 0 }`
- `.variant-block { page-break-inside: avoid }`
- `.detail-grid { page-break-inside: avoid }`
- `table { page-break-inside: avoid }` — **problem for B3**: a long numbered references list will be table-like and `page-break-inside: avoid` will force it onto one page. Either render refs as `<ol>` (not `<table>`) and remove the avoid constraint, or use a dedicated `.references-list { page-break-inside: auto }` override.

**New elements v2.2 introduces and their print impact:**

| Element | Where | Risk | Mitigation |
|---|---|---|---|
| Patient strip | every `<div class="page">` | Adds ~14–18px of vertical space per page. Existing pages are sized for exact body height, so this **will** push the last line of dense pages (e.g., the variant-detail 2-per-page layout) below the A4 break. | Shrink `.page-body` padding top by the same amount; verify with WeasyPrint print preview; add snapshot test against the three showcase JSONs. |
| Curated treatments table | Cancer Board page | Adds a `source` column (OncoKB/CIViC/both) and potentially `pmids` — table could overflow horizontally. | Keep font 10px, table width 100%, `word-break: break-word` on the drug name, drop `resistance_notes` to its own row below the drug (two-row layout per entry) if needed. |
| References numbered list | end of Cancer Board section | Could push the Cancer Board onto page 2 of that section unexpectedly. | Give it `page-break-before: auto` and allow it to overflow naturally. Do **not** wrap in a `<table>`. |
| Override notes (A4) | variant-detail card | Adds variable-height text to already-compact variant blocks. Could push the second-variant-per-page below the fold. | Cap `clinvar_override_reason` display to ~200 chars with a `... (full reason in methodology section)` tail, and allow `.variant-block` height to be dynamic in this case (remove page-break-inside: avoid when override is present, or split to 1-variant-per-page). |

**Required verification before Phase B merges:** run the three existing showcase JSONs through `rerender_report.py` after each of B3 and B4 and diff the rendered PDF print-preview (not just HTML) against the current committed PDFs. Add this to the qa-engineer task list.

### Bilingual (KO/EN)

The plan does not enumerate every new UI string that needs a KO pair. Here is the missing set (all new in v2.2):

| EN | KO | Where |
|---|---|---|
| Patient | 환자 | B4 header |
| DOB | 생년월일 | B4 header |
| Sex | 성별 | B4 header |
| MRN | 의무기록번호 | B4 header (optional) |
| Ordering Physician | 의뢰 의사 | B4 header |
| Facility | 의뢰 기관 | B4 header |
| Specimen | 검체 | B4 header |
| Collected | 채취일 | B4 header |
| Curated Treatments | 큐레이션된 치료 근거 | A1/A2 table title |
| Source | 출처 | A1/A2 table column |
| References | 참고문헌 | B3 numbered list title |
| Override Notes | 오버라이드 근거 | A4 footnote block |
| Unverified source | 미확인 출처 | B3 when PMID has no curated backing |
| No curated treatment evidence... | 이 변이군에 대해 큐레이션된 치료 근거가 없습니다 — 전문의 소견 참조 | A2 empty state |

**Where translations live.** `render.py` already uses inline dicts keyed by language code (see `_NO_FINDINGS_MESSAGES`). Continue that pattern — do **not** introduce a gettext/po file infrastructure for 15 strings. For Jinja templates, use `{% if language == "ko" %}...{% else %}...{% endif %}` inline or add a `t(key, lang)` macro at the top of the template. Either is acceptable; pick one and apply consistently.

### Backward compatibility with `rerender_report.py`

Summary of breaking changes v2.2 introduces to the `report_data` / `opinion` shape and the required mitigations (all missing from the plan):

| Change | Breaking? | Mitigation required |
|---|---|---|
| `AgentOpinion.references: List[str]` → `List[Reference]` | **YES** — TypeError in render | Coerce in `_reconstruct_board_opinion`; update `_render_agent_opinions_section` iteration. |
| `CancerBoardOpinion.references` new field | No (default=[]) | `_filter_dataclass_kwargs` drops unknown keys; absent field defaults. ✓ |
| `CancerBoardOpinion.treatment_options[].curated_id` new | No if renderer tolerates missing key | `opt.get("curated_id", None)` — render path must not KeyError. |
| `ClassificationResult.clinvar_override_reason` new field | No (default="") | Template uses `{{ v.clinvar_override_reason \| default('') }}`. ✓ |
| `report_data["patient"]` new key | No | Template gates via `{% if patient %}`. ✓ |
| `report_data["_curated_treatments"]` new key | No | Renderer tolerates absent. ✓ |

**Required:** add a new task **B3.5** or fold into B3: "Update `rerender_report.py` coercion and add a regression test that re-renders all three `docs/showcase/*.json` files through `rerender_report.py` and diffs the output HTML against the committed `.html` file — must match byte-for-byte (modulo a date stamp)."

---

## Open questions answered (for clinical-advisor)

Re-stated from draft cover note §"For report-dev":

**Q1. Is the existing `templates/cancer/report.html` layout flexible enough for a running patient header on every page, or does adding a header require a Jinja layout overhaul?**
**A:** Partial. A Jinja macro + mechanical replacement of every `<div class="page-header-bar"></div>` is ~1.5 days of work, not 1 day. Avoid WeasyPrint running-elements approach — it requires changing `@page { margin: 0 }` which cascades to every page's layout.

**Q2. PMID superscript rendering — is there an existing reference-numbering utility in the codebase, or do we need a new one?**
**A:** None exists. Add one in `render.py` (`_superscriptify(text, ref_index)` pattern above). Keep it internal to clinical_board render; do not share with other templates — the AI Board is the only section that needs inline citations in v2.2.

**Q3. Should the curated Treatment Options table replace the current "Treatment Options" section entirely, or live alongside it during a transition period?**
**A:** Replace entirely. Side-by-side A/B would duplicate code paths, confuse clinicians, and double the A4 vertical budget. The plan's "narrate-only" architecture is the whole point; A/B defeats it.

**Q4. The `clinvar_override_reason` footnote (A4) — where should it render? Under the variant card, or as a separate "Override Notes" block at the bottom of the Classifications page?**
**A:** Under the variant detail card (inside `.variant-detail-header`). A separate block at the bottom loses the reason-variant coupling. BUT: the render-target is the Jinja template `templates/cancer/report.html`, not `scripts/clinical_board/render.py`. Please correct the file map in A4.

---

## Punch list (top → bottom priority)

1. **[CRITICAL / B3]** Define the LLM PMID marker protocol explicitly (`[pmid:NNNNNNNN]` recommended) in all cancer agent prompts. Without this, `render.py` has no way to place inline superscripts. Add to acceptance criteria.

2. **[CRITICAL / B3+rerender]** Add a forward-compat coercion in `rerender_report.py::_reconstruct_board_opinion` that upgrades `list[str]` → `list[Reference]`. Update `_render_agent_opinions_section` to iterate `ref.pmid` instead of joining raw strings. Bundle in the same commit as the model change. Add a regression test that re-renders all three existing `docs/showcase/*.json` files without error.

3. **[HIGH / B4]** Clarify render location: add a Jinja macro `patient_strip(...)` in `templates/cancer/report.html` (and `templates/rare-disease/report.html`), replace every `<div class="page-header-bar"></div>` with the macro call. Revise effort S→S/M (1.5 days). Add snapshot test that every `<div class="page">` in both templates contains the macro call.

4. **[HIGH / A4]** Correct the file map: `clinvar_override_reason` rendering lives in `templates/cancer/report.html` (variant detail page), **not** in `scripts/clinical_board/render.py`. Add `.override-notice` styling with `page-break-inside: avoid` and a length cap.

5. **[HIGH / B1]** Correct the file map: Tier III consequence gate lives in `variant_selector.py` — there is no Tier III VUS table in `clinical_board/render.py`. Remove `clinical_board/render.py` from B1's file list.

6. **[MEDIUM / A1+A2]** Specify that curated_id hydration happens in `runner.py` (not `render.py`), because `render_board_opinion_html(opinion, language)` does not receive `report_data`. Plumbing `report_data` into render would break `rerender_report.py`'s contract.

7. **[MEDIUM / PHI]** Extend `reporting.persist_patient_metadata=false` to also gate JSON dumps under `docs/showcase/` and `_workspace/`, not just SQLite persistence. Otherwise PHI leaks into git.

8. **[MEDIUM / bilingual]** Add the 15-string KO/EN table (above) to the plan's B3/B4 acceptance criteria. Pick the in-template `{% if language == "ko" %}` pattern (matches existing `_disclaimer_text`); do not introduce gettext.

9. **[MEDIUM / print CSS]** Verify that adding the patient strip does not push the last line of dense pages below the A4 break. Require `rerender_report.py` + print-preview diff against the three committed showcases before B4 merges.

10. **[LOW / B3 effort]** Revise B3 effort from M (1.5–2d) to M/L (2.5–3d) to cover marker protocol, coercion, KO/EN pairs, and A4 overflow.

11. **[LOW / table width]** Curated treatments table: if adding a `source` column, keep total columns ≤4 (Drug, Level, Source, Resistance Notes) at 10px font. Drop the curated_id from user-visible rendering.

12. **[LOW / CSS]** Do not wrap the references numbered list in a `<table>` — use `<ol>` so it can page-break naturally. The global `table { page-break-inside: avoid }` rule would otherwise force an oversized block onto a single page.

---

## Summary for team-lead

**Verdict:** approve with required changes. The two critical issues (B3 marker protocol + `rerender_report.py` coercion) are blocking — without them the showcase regeneration will crash and the feature will not be reproducible. Both are small (~half day total of spec work) but must land in the same commit as the model change.

**Top 3 punch list:**
1. **B3:** Define `[pmid:NNNNNNNN]` marker protocol in cancer agent prompts + add `rerender_report.py` coercion from `list[str]` to `list[Reference]`. Without both, existing showcases break.
2. **B4:** Render patient header via a Jinja macro replacing every `page-header-bar`. Effort S → S/M (1.5d). Do **not** use WeasyPrint running-element `@top-center` — it requires scrapping `@page { margin: 0 }`.
3. **File map corrections:** A4's override reason rendering is in `templates/cancer/report.html` (variant-detail card), not `clinical_board/render.py`. B1's Tier III consequence gate has no render-layer touch at all.

— report-dev
