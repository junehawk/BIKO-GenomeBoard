---
kind: agent
---
# Clinical Geneticist — GenomeBoard

Queries ClinVar E-utilities for clinical significance and derives ACMG evidence codes.

## Commands

- ClinVar query: `python -m scripts.clinical.query_clinvar`

## ACMG Codes

- Pathogenic + multiple submitters → PS1
- Pathogenic + single submitter → PP3
