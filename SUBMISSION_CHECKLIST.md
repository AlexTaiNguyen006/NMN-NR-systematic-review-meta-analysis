# Submission Readiness Checklist

**Paper:** Nicotinamide Mononucleotide vs Nicotinamide Riboside — Systematic Review & Network Meta-Analysis  
**Date verified:** 2026-04-01  

---

## 1. Manuscript (DOCX)

| Item | Status |
|------|--------|
| Title, Abstract, Intro, Methods, Results, Discussion, Conclusion | ✅ All present |
| Word count (~7,408 words) | ✅ Verified |
| References [1]–[33]: all cited and all present | ✅ Verified |
| Figure citations (Figures 1–6, S1–S6, 2A, 2B) | ✅ Verified |
| Table citations (Tables 1–2, S1–S9) | ✅ Verified |
| 10 embedded tables correctly structured | ✅ Verified |
| No broken cross-references or placeholders | ✅ Verified (except 4 author fields — see below) |

### Author fields still requiring completion

| Paragraph | Field | Current text |
|-----------|-------|-------------|
| Funding | Funding source | `[To be completed before submission]` |
| Conflicts of Interest | COI declaration | `[To be completed before submission]` |
| Data Availability | Repository URL | `[repository URL to be added before submission]` |
| Author Contributions | CRediT roles | `[To be completed before submission]` |

---

## 2. Data–Manuscript Number Consistency

| Check | Status |
|-------|--------|
| 28 pairwise MD/CI/p/I²/k values vs `pairwise_meta_analysis.csv` | ✅ All match |
| 14 indirect comparison MD/CI/p values vs `nma_results.csv` | ✅ All match |
| PRISMA flow counts (1687→562 dup→1125→874 excl→251 FT→15 incl→8 quant, 73 data points) | ✅ All match |
| LOO sensitivity: 56 iterations, 1 significance change (TG/Conze_2019) | ✅ Matches manuscript |
| GRADE/CINeMA: 14 outcomes, all Very Low certainty | ✅ Matches manuscript |
| Brakedal 2022 (dose=1000, duration=4.3 wk, Parkinson's) in DOCX + S2 CSV | ✅ Correct |
| Orr 2024 (dose=1000, duration=10 wk, MCI) in DOCX + S2 CSV | ✅ Correct |

---

## 3. Generated Outputs — All Regenerated from Corrected Source Data

### Figures (`results/figures/`)

| Figure | File(s) | Status |
|--------|---------|--------|
| 19 forest plots (pairwise) | `forest_*.png` | ✅ Regenerated |
| NMA forest (indirect comparisons) | `nma_forest_NMN_vs_NR_all.png/.pdf` | ✅ Regenerated |
| Network graph | `network_graph.png` | ✅ Regenerated |
| RoB 2 summary + traffic light | `rob2_summary.png/.pdf`, `rob2_traffic_light.png/.pdf` | ✅ Regenerated |
| PRISMA 2020 flow | `prisma_2020_flow.png/.pdf` | ✅ Regenerated |
| GRADE/CINeMA heatmap | `grade_cinema_heatmap.png/.pdf` | ✅ Regenerated |
| LOO indirect summary | `sensitivity_loo_indirect_summary.png/.pdf` | ✅ Regenerated |
| LOO pairwise forests (3) | `sensitivity_loo_forest_*.png/.pdf` | ✅ Regenerated |
| Outcome reporting matrix | `outcome_reporting_matrix.png` | ✅ Regenerated |

### Tables (`results/tables/`)

| Table | File | Status |
|-------|------|--------|
| S1: Search strategy | `supp_table_S1_search_strategy.csv` | ✅ Regenerated |
| S2: Study characteristics | `supp_table_S2_study_characteristics.csv` | ✅ Regenerated |
| S3: Full extraction | `supp_table_S3_full_extraction.csv` | ✅ Regenerated |
| S4: RoB 2 detailed | `supp_table_S4_rob2_detailed.csv` | ✅ Regenerated |
| S5: Pairwise all | `supp_table_S5_pairwise_all.csv` | ✅ Regenerated |
| S6: NMA all | `supp_table_S6_nma_all.csv` | ✅ Regenerated |
| S7: LOO summary | `supp_table_S7_loo_summary.csv` | ✅ Regenerated |
| S8: GRADE summary | `supp_table_S8_grade_summary.csv` | ✅ Regenerated |
| S9: Inter-rater agreement | `supp_table_S9_inter_rater_agreement.csv` | ✅ Regenerated |
| Core result tables | `pairwise_meta_analysis.csv`, `nma_results.csv`, `league_table.csv`, `grade_cinema_assessment.csv`, sensitivity CSVs | ✅ Regenerated |

---

## 4. Source Scripts

| Script | Purpose | Status |
|--------|---------|--------|
| `regenerate_figures.py` | All forest plots, network graph, RoB, GRADE heatmap | ✅ Runs clean |
| `supplementary_tables.py` | S1–S9 supplementary tables | ✅ Runs clean |
| `sensitivity_analyses.py` | LOO pairwise, LOO indirect, high-RoB exclusion | ✅ Runs clean |
| `grade_cinema.py` | GRADE/CINeMA assessment | ✅ Runs clean |
| `prisma_flow.py` | PRISMA 2020 flow diagram | ✅ Runs clean |
| `metafor_analysis.R` | R/metafor REML validation | Present (not re-run this session) |
| `data_verification.py` | Data integrity checks | Present |

---

## 5. Registration & Reporting

| Item | Status |
|------|--------|
| PROSPERO registration (CRD42025636869) | ✅ Noted as retrospective in manuscript |
| PRISMA 2020 adherence | ✅ Stated in Methods §2.1 |
| GRADE/CINeMA certainty ratings | ✅ All 14 outcomes rated Very Low |

---

## 6. Pre-Submission Actions Required

- [ ] Complete **Funding** declaration
- [ ] Complete **Conflicts of Interest** declaration
- [ ] Complete **Author Contributions** (CRediT taxonomy)
- [ ] Add **Data Availability** repository URL
- [ ] Confirm target journal formatting requirements
- [ ] Export final figures at journal-required resolution/format
- [ ] Prepare cover letter
- [ ] Confirm all co-author approvals

---

## 7. Workspace State

- `manuscript/` — Clean: 1 DOCX + 1 transcript
- `src/` — Clean: 7 production scripts + 1 R script + 1 utility
- `results/figures/` — 39 generated output files (PNG + PDF)
- `results/tables/` — 23 generated CSV files
- `archive/` — 3 DOCX backups + prior cleanup artifacts
- No temporary files remaining
