---
name: Bug report
about: Report a problem with the pipeline or classification logic
title: '[BUG] '
labels: bug
---

## Description

<!-- A clear and concise description of what the bug is. -->

## Steps to reproduce

1. Input VCF (synthetic / source):
2. Command line used:
   ```bash
   python -m scripts.orchestrate ...
   ```
3. Mode: <!-- cancer | rare-disease -->

## Expected behavior

<!-- What did you expect the classification, tiering, or report output to look like? -->

## Actual behavior

<!-- What actually happened? Include the relevant variant(s), classification verdicts, and any incorrect tier or evidence codes. -->

## Environment

- Python version:
- OS / distribution:
- BIKO GenomeBoard version or commit SHA:
- Ollama version (if Clinical Board path is involved):
- Local DB build date (`scripts/db/version_manager.py` output, if relevant):

## Logs / traceback

<details>
<summary>Logs</summary>

```
<!-- paste here -->
```

</details>

---

### Privacy notice — required reading before submitting

BIKO GenomeBoard is a **research reference tool** and this is a
**public** issue tracker. **Do not** paste patient identifiers,
specimen IDs, MRNs, DOBs, ordering physician names, institution
names, or real VCF content that could identify an individual.

When describing a variant, use paste-safe descriptors only:
`gene + genomic position + reference/alt + consequence` (e.g.
`TP53 chr17:7577120 G>A missense`). That is enough context for
triage — anything more may constitute PHI.

If your report involves PHI or a suspected security issue, **stop**
and follow [`SECURITY.md`](../../SECURITY.md) instead of opening a
public issue.
