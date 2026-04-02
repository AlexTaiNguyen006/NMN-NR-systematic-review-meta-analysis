# NMN vs NR: Systematic Review and Network Meta-Analysis

Systematic review and network meta-analysis comparing nicotinamide mononucleotide (NMN) and nicotinamide riboside (NR) supplementation on cardiometabolic outcomes in human RCTs.

**PROSPERO Registration:** CRD420261330487

## Overview

This project includes 15 RCTs (5 NMN, 10 NR) in qualitative synthesis, of which 8 (4 NMN, 4 NR) contributed extractable data for quantitative network meta-analysis, yielding 73 effect sizes across 21 outcomes. Key analyses:

- Pairwise random-effects meta-analyses (NMN vs placebo, NR vs placebo)
- Indirect network meta-analysis comparisons (NMN vs NR via Bucher method)
- R/metafor REML validation
- Leave-one-out sensitivity analyses
- GRADE/CINeMA certainty of evidence ratings

## Repository Structure

```
data/
  extraction/      Curated data extraction files and RoB 2 assessments
  screening/       Search records from PubMed, Embase, Scopus, WoS, CENTRAL
src/               Analysis scripts (Python + R)
results/
  figures/         Forest plots, PRISMA flow, RoB figures, etc.
  tables/          Pairwise results, NMA results, supplementary tables
```

## Analysis Scripts

| Script | Purpose |
|--------|---------|
| `regenerate_figures.py` | Forest plots, network graph, RoB figures, GRADE heatmap |
| `sensitivity_analyses.py` | Leave-one-out and high-RoB exclusion analyses |
| `grade_cinema.py` | GRADE/CINeMA certainty of evidence assessment |
| `supplementary_tables.py` | Supplementary Tables S1–S9 |
| `prisma_flow.py` | PRISMA 2020 flow diagram |
| `metafor_analysis.R` | R/metafor REML validation of pairwise and indirect results |
| `data_verification.py` | Cross-checks extracted data against source papers |
| `fix_all_extractions.py` | Documents and applies data extraction corrections |
| `outcome_matrix.py` | Outcome reporting matrix heatmap |
| `publication_bias.py` | Funnel plots and Egger's regression test for publication bias |

## Requirements

- Python 3.11 with numpy, scipy, matplotlib
- R with metafor, readr, dplyr

See `environment.yml` for the conda environment specification.

## Reproducibility

All results (figures and tables) can be regenerated from the corrected source data in `data/extraction/nma_input_long.csv` using the scripts in `src/`.

## Citation

Citation metadata is provided in `CITATION.cff`.

For archival citation, use the Zenodo DOI associated with release tag `v1.0.0` once Zenodo records the release.
