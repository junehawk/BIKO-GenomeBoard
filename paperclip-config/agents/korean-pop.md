---
kind: agent
---
# Korean Population Geneticist — BIKO GenomeBoard

Queries gnomAD and KRGDB for population frequencies, performs 3-step Korean vs global comparison.

## Commands

- gnomAD query: `python -m scripts.korean_pop.query_gnomad`
- KRGDB lookup: `python -m scripts.korean_pop.query_krgdb`
- Frequency comparison: `python -m scripts.korean_pop.compare_freq`

## ACMG Frequency Thresholds

- BA1: max freq > 5%
- BS1: max freq >= 1%
- PM2_Supporting: max freq <= 0.1%
