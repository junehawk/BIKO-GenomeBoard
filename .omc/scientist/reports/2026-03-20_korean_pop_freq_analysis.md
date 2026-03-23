# Korean Population Frequency Analysis

**Date:** 2026-03-20  
**Stage:** RESEARCH_STAGE:2  
**Dataset:** GenomeBoard demo VCF — 10 variants

## [OBJECTIVE]

Analyze Korean population allele frequencies for 10 GenomeBoard demo variants. Apply compare_freq.py thresholds (BA1>0.05, BS1>=0.01, PM2_Supporting<=0.001) to assign ACMG frequency codes. Flag Korean-specific enrichment (Korean/global >=5x) or depletion (Korean/global <=0.2x). Special focus on CYP2C19 (rs4244285) and HLA-B (rs2395029).

**Sources:** KRGDB, gnomAD v4 EAS/global, korean_pgx_table.json v2026-03-20 (CPIC 2022, PharmGKB, Yang et al. 2014)

## [DATA]

- 10 variants, 10 genes
- BA1 (common >5%): 3 variants (CYP2C19, HLA-B, APOE)
- BS1 (strong benign 1-5%): 1 variant (CFTR globally; Korean-context override to PM2_Supporting)
- PM2_Supporting (rare <=0.1%): 6 variants (TP53, BRCA2, ATM, MUTYH, PALB2, PTPN11)
- Korean-specific flags: 3 variants (V2 enrichment, V3 depletion critical, V4 enrichment critical)

