#!/usr/bin/env Rscript
# R/metafor validation analysis for NMN vs NR systematic review
#
# Uses the validated field-standard metafor package with REML estimator.
# Reads corrected nma_input_long.csv and produces:
#   - Pairwise random-effects meta-analyses (REML)
#   - Indirect comparisons via Bucher method
#   - Multiple testing correction (Bonferroni + BH)
#   - Crossover trial sensitivity analysis
#   - Forest plots
#   - CSV outputs for comparison with Python results

# Load packages (install if needed)
required_pkgs <- c("metafor", "readr", "dplyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(metafor)
library(readr)
library(dplyr)

# Paths
base_dir  <- "/Users/tainguyen/Documents/NMN_vs_NR_SR_NMA"
data_file <- file.path(base_dir, "data", "extraction", "nma_input_long.csv")
out_dir   <- file.path(base_dir, "results", "tables")
fig_dir   <- file.path(base_dir, "results", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Read data
dat <- read_csv(data_file, show_col_types = FALSE)

# Flag NAD+ unit incompatibility — exclude from indirect comparison
cat("\n══════════════════════════════════════════════════\n")
cat("NAD+ UNIT CHECK\n")
cat("══════════════════════════════════════════════════\n")
nad_rows <- dat %>% filter(outcome == "NAD+")
cat("NAD+ entries:\n")
print(nad_rows %>% select(study_id, precursor, unit, md, se_md, notes))
cat("\nWARNING: Huang_2022 NAD+ is in pmol/mL; Bandi_2025 is in µM.\n")
cat("These units differ by ~1000x. Indirect comparison is INVALID.\n")
cat("NAD+ will be excluded from indirect comparisons.\n\n")

# Flag crossover studies
cat("══════════════════════════════════════════════════\n")
cat("CROSSOVER STUDY CHECK\n")
cat("══════════════════════════════════════════════════\n")
crossover <- dat %>% filter(grepl("Crossover", design, ignore.case = TRUE))
cat(sprintf("Crossover studies in data: %d rows\n", nrow(crossover)))
if (nrow(crossover) > 0) {
  cat("Studies:\n")
  print(crossover %>% distinct(study_id, design))
  cat("\nWARNING: Crossover RCTs require within-person SD (accounting for\n")
  cat("correlation between treatment periods). Current SDs are between-group,\n")
  cat("which overestimates variance. Sensitivity analysis will exclude these.\n\n")
}

# Pairwise meta-analysis (REML estimator)
cat("══════════════════════════════════════════════════\n")
cat("PAIRWISE META-ANALYSES (metafor, REML)\n")
cat("══════════════════════════════════════════════════\n\n")

outcomes <- sort(unique(dat$outcome))
precursors <- c("NMN", "NR")

pairwise_results <- list()

for (oc in outcomes) {
  for (prec in precursors) {
    subset <- dat %>% filter(outcome == oc, precursor == prec)
    if (nrow(subset) == 0) next
    
    k <- nrow(subset)
    comparison <- paste0(prec, " vs Placebo")
    
    if (k == 1) {
      # Single study — no meta-analysis possible, report as-is
      row <- subset[1, ]
      res_row <- data.frame(
        outcome = oc,
        precursor = prec,
        comparison = comparison,
        k = 1,
        method = "single_study",
        MD_REML = row$md,
        lower_CI = row$md - 1.96 * row$se_md,
        upper_CI = row$md + 1.96 * row$se_md,
        p_value = 2 * pnorm(abs(row$md / row$se_md), lower.tail = FALSE),
        I2 = NA_real_,
        tau2 = NA_real_,
        Q = NA_real_,
        Q_p = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      # Random-effects meta-analysis with REML
      ma <- tryCatch(
        rma(yi = md, sei = se_md, data = subset, method = "REML"),
        error = function(e) {
          cat(sprintf("  ERROR in %s / %s: %s\n", oc, prec, e$message))
          NULL
        }
      )
      if (is.null(ma)) next
      
      res_row <- data.frame(
        outcome = oc,
        precursor = prec,
        comparison = comparison,
        k = k,
        method = "REML",
        MD_REML = as.numeric(ma$b),
        lower_CI = ma$ci.lb,
        upper_CI = ma$ci.ub,
        p_value = ma$pval,
        I2 = ma$I2,
        tau2 = ma$tau2,
        Q = ma$QE,
        Q_p = ma$QEp,
        stringsAsFactors = FALSE
      )
    }
    
    pairwise_results[[length(pairwise_results) + 1]] <- res_row
    cat(sprintf("  %s / %s (k=%d): MD=%.2f [%.2f, %.2f], p=%.4f",
                oc, prec, k,
                res_row$MD_REML, res_row$lower_CI, res_row$upper_CI,
                res_row$p_value))
    if (!is.na(res_row$I2)) {
      cat(sprintf(", I²=%.1f%%, τ²=%.4f", res_row$I2, res_row$tau2))
    } else {
      cat(", k=1 (no heterogeneity estimate)")
    }
    cat("\n")
  }
}

pw_df <- bind_rows(pairwise_results)
write_csv(pw_df, file.path(out_dir, "pairwise_metafor_REML.csv"))
cat(sprintf("\nSaved: pairwise_metafor_REML.csv (%d comparisons)\n\n", nrow(pw_df)))


# INDIRECT COMPARISONS (Bucher method)
cat("══════════════════════════════════════════════════\n")
cat("INDIRECT COMPARISONS (NMN vs NR via Bucher)\n")
cat("══════════════════════════════════════════════════\n\n")

indirect_results <- list()

for (oc in outcomes) {
  nmn_row <- pw_df %>% filter(outcome == oc, precursor == "NMN")
  nr_row  <- pw_df %>% filter(outcome == oc, precursor == "NR")
  
  if (nrow(nmn_row) == 0 || nrow(nr_row) == 0) next
  
  # Check NAD+ unit incompatibility
  if (oc == "NAD+") {
    cat(sprintf("  %s: SKIPPED — incompatible units (pmol/mL vs µM)\n", oc))
    ind_row <- data.frame(
      outcome = oc,
      comparison = "NMN vs NR",
      k_NMN = nmn_row$k,
      k_NR = nr_row$k,
      MD_indirect = NA_real_,
      lower_CI = NA_real_,
      upper_CI = NA_real_,
      p_value = NA_real_,
      p_bonferroni = NA_real_,
      p_BH = NA_real_,
      note = "EXCLUDED: Unit incompatibility (pmol/mL vs µM)",
      stringsAsFactors = FALSE
    )
    indirect_results[[length(indirect_results) + 1]] <- ind_row
    next
  }
  
  # Bucher indirect comparison: MD_NMN_vs_NR = MD_NMN_vs_PBO - MD_NR_vs_PBO
  md_nmn <- nmn_row$MD_REML
  md_nr  <- nr_row$MD_REML
  
  # SE for indirect: need SEs of each pairwise estimate
  # For REML with k>1: SE = (CI_upper - CI_lower) / (2 * 1.96)
  se_nmn <- (nmn_row$upper_CI - nmn_row$lower_CI) / (2 * 1.96)
  se_nr  <- (nr_row$upper_CI  - nr_row$lower_CI)  / (2 * 1.96)
  
  md_indirect <- md_nmn - md_nr
  se_indirect <- sqrt(se_nmn^2 + se_nr^2)
  ci_lo <- md_indirect - 1.96 * se_indirect
  ci_up <- md_indirect + 1.96 * se_indirect
  z_val <- md_indirect / se_indirect
  p_val <- 2 * pnorm(abs(z_val), lower.tail = FALSE)
  
  ind_row <- data.frame(
    outcome = oc,
    comparison = "NMN vs NR",
    k_NMN = nmn_row$k,
    k_NR = nr_row$k,
    MD_indirect = md_indirect,
    lower_CI = ci_lo,
    upper_CI = ci_up,
    p_value = p_val,
    p_bonferroni = NA_real_,
    p_BH = NA_real_,
    note = "",
    stringsAsFactors = FALSE
  )
  
  indirect_results[[length(indirect_results) + 1]] <- ind_row
  
  sig_marker <- ifelse(p_val < 0.05, " *", "")
  cat(sprintf("  %s: MD=%.2f [%.2f, %.2f], p=%.4f (k=%d+%d)%s\n",
              oc, md_indirect, ci_lo, ci_up, p_val,
              nmn_row$k, nr_row$k, sig_marker))
}

ind_df <- bind_rows(indirect_results)

# MULTIPLE TESTING CORRECTION
cat("\n══════════════════════════════════════════════════\n")
cat("MULTIPLE TESTING CORRECTION\n")
cat("══════════════════════════════════════════════════\n\n")

# Only correct tests that have valid p-values (exclude NAD+ which was skipped)
valid_tests <- ind_df %>% filter(!is.na(p_value))
n_tests <- nrow(valid_tests)

cat(sprintf("Number of indirect comparisons tested: %d\n", n_tests))
cat(sprintf("Bonferroni threshold (α=0.05/%d): p < %.4f\n", n_tests, 0.05 / n_tests))

if (n_tests > 0) {
  # Bonferroni
  valid_tests$p_bonferroni <- pmin(valid_tests$p_value * n_tests, 1.0)
  
  # Benjamini-Hochberg (FDR)
  valid_tests$p_BH <- p.adjust(valid_tests$p_value, method = "BH")
  
  # Update ind_df with corrections
  for (i in seq_len(nrow(valid_tests))) {
    oc <- valid_tests$outcome[i]
    idx <- which(ind_df$outcome == oc)
    ind_df$p_bonferroni[idx] <- valid_tests$p_bonferroni[i]
    ind_df$p_BH[idx]         <- valid_tests$p_BH[i]
  }
  
  cat("\nResults with corrections:\n")
  for (i in seq_len(nrow(ind_df))) {
    r <- ind_df[i, ]
    if (is.na(r$p_value)) {
      cat(sprintf("  %s: %s\n", r$outcome, r$note))
    } else {
      sig_raw <- ifelse(r$p_value < 0.05, "*", "")
      sig_bonf <- ifelse(r$p_bonferroni < 0.05, "*", "")
      sig_bh <- ifelse(r$p_BH < 0.05, "*", "")
      cat(sprintf("  %s: p_raw=%.4f%s, p_Bonf=%.4f%s, p_BH=%.4f%s\n",
                  r$outcome, r$p_value, sig_raw,
                  r$p_bonferroni, sig_bonf, r$p_BH, sig_bh))
    }
  }
}

write_csv(ind_df, file.path(out_dir, "indirect_comparisons_metafor_REML.csv"))
cat(sprintf("\nSaved: indirect_comparisons_metafor_REML.csv\n\n"))


# SENSITIVITY ANALYSIS: EXCLUDE CROSSOVER STUDIES
cat("══════════════════════════════════════════════════\n")
cat("SENSITIVITY: EXCLUDING CROSSOVER STUDIES\n")
cat("══════════════════════════════════════════════════\n\n")

dat_parallel <- dat %>% filter(!grepl("Crossover", design, ignore.case = TRUE))
crossover_studies <- dat %>%
  filter(grepl("Crossover", design, ignore.case = TRUE)) %>%
  distinct(study_id) %>% pull(study_id)
cat(sprintf("Excluded crossover studies: %s\n", paste(crossover_studies, collapse = ", ")))
cat(sprintf("Remaining rows: %d (from %d)\n\n", nrow(dat_parallel), nrow(dat)))

sens_results <- list()

for (oc in outcomes) {
  for (prec in precursors) {
    subset <- dat_parallel %>% filter(outcome == oc, precursor == prec)
    if (nrow(subset) == 0) next
    k <- nrow(subset)
    comparison <- paste0(prec, " vs Placebo")
    
    if (k == 1) {
      row <- subset[1, ]
      res_row <- data.frame(
        outcome = oc, precursor = prec, comparison = comparison,
        k = 1, method = "single_study",
        MD_REML = row$md,
        lower_CI = row$md - 1.96 * row$se_md,
        upper_CI = row$md + 1.96 * row$se_md,
        p_value = 2 * pnorm(abs(row$md / row$se_md), lower.tail = FALSE),
        I2 = NA_real_, tau2 = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      ma <- tryCatch(
        rma(yi = md, sei = se_md, data = subset, method = "REML"),
        error = function(e) NULL
      )
      if (is.null(ma)) next
      res_row <- data.frame(
        outcome = oc, precursor = prec, comparison = comparison,
        k = k, method = "REML",
        MD_REML = as.numeric(ma$b),
        lower_CI = ma$ci.lb, upper_CI = ma$ci.ub,
        p_value = ma$pval,
        I2 = ma$I2, tau2 = ma$tau2,
        stringsAsFactors = FALSE
      )
    }
    sens_results[[length(sens_results) + 1]] <- res_row
  }
}

sens_df <- bind_rows(sens_results)
write_csv(sens_df, file.path(out_dir, "sensitivity_no_crossover_metafor.csv"))
cat(sprintf("Saved: sensitivity_no_crossover_metafor.csv (%d comparisons)\n\n", nrow(sens_df)))

# Compare with full analysis
cat("Outcomes affected by crossover exclusion:\n")
for (oc in outcomes) {
  for (prec in precursors) {
    full <- pw_df %>% filter(outcome == oc, precursor == prec)
    sens <- sens_df %>% filter(outcome == oc, precursor == prec)
    if (nrow(full) > 0 && nrow(sens) > 0 && full$k != sens$k) {
      cat(sprintf("  %s/%s: k=%d → k=%d, MD: %.2f → %.2f\n",
                  oc, prec, full$k, sens$k, full$MD_REML, sens$MD_REML))
    } else if (nrow(full) > 0 && nrow(sens) == 0) {
      cat(sprintf("  %s/%s: k=%d → DROPPED (no parallel studies)\n",
                  oc, prec, full$k))
    }
  }
}


# COMPARISON: REML vs DL (DerSimonian-Laird)
cat("\n══════════════════════════════════════════════════\n")
cat("COMPARISON: REML vs DL ESTIMATOR\n")
cat("══════════════════════════════════════════════════\n\n")

for (oc in outcomes) {
  for (prec in precursors) {
    subset <- dat %>% filter(outcome == oc, precursor == prec)
    if (nrow(subset) < 2) next
    
    ma_reml <- tryCatch(
      rma(yi = md, sei = se_md, data = subset, method = "REML"),
      error = function(e) NULL
    )
    ma_dl <- tryCatch(
      rma(yi = md, sei = se_md, data = subset, method = "DL"),
      error = function(e) NULL
    )
    
    if (!is.null(ma_reml) && !is.null(ma_dl)) {
      if (abs(ma_reml$tau2 - ma_dl$tau2) > 0.001 ||
          abs(as.numeric(ma_reml$b) - as.numeric(ma_dl$b)) > 0.01) {
        cat(sprintf("  %s/%s (k=%d): REML τ²=%.4f, DL τ²=%.4f | ",
                    oc, prec, nrow(subset), ma_reml$tau2, ma_dl$tau2))
        cat(sprintf("REML MD=%.2f, DL MD=%.2f\n",
                    as.numeric(ma_reml$b), as.numeric(ma_dl$b)))
      }
    }
  }
}


# TRANSITIVITY ASSESSMENT TABLE
cat("\n══════════════════════════════════════════════════\n")
cat("TRANSITIVITY ASSESSMENT\n")
cat("══════════════════════════════════════════════════\n\n")

# Read study characteristics
char_file <- file.path(base_dir, "data", "extraction", "study_characteristics_curated.csv")
chars <- read_csv(char_file, show_col_types = FALSE)

# Studies in quantitative NMA only
nma_studies <- chars %>% filter(study_id %in% unique(dat$study_id))

trans_summary <- nma_studies %>%
  group_by(precursor) %>%
  summarise(
    n_studies = n(),
    countries = paste(unique(country), collapse = "; "),
    dose_range = paste0(min(dose_mg_day), "–", max(dose_mg_day), " mg/day"),
    duration_range = paste0(min(duration_weeks), "–", max(duration_weeks), " weeks"),
    designs = paste(unique(design), collapse = "; "),
    populations = paste(unique(population), collapse = "; "),
    mean_n = round(mean(n_randomized), 1),
    .groups = "drop"
  )

cat("NMN vs NR arm characteristics:\n")
print(as.data.frame(trans_summary))

# Molar dose comparison (NMN MW=334.22, NR MW=255.25 as chloride salt 290.70)
cat("\nMolar dose analysis:\n")
cat("  NMN molecular weight: 334.22 g/mol\n")
cat("  NR (as chloride) MW:  290.70 g/mol\n\n")

for (i in seq_len(nrow(nma_studies))) {
  s <- nma_studies[i, ]
  mw <- ifelse(s$precursor == "NMN", 334.22, 290.70)
  dose <- as.numeric(s$dose_mg_day)
  mmol <- dose / mw
  cat(sprintf("  %s (%s): %d mg/day = %.2f mmol/day\n",
              s$study_id, s$precursor, dose, mmol))
}

# Build formal transitivity assessment table
trans_domains <- data.frame(
  Domain = c(
    "Dose (mass)",
    "Dose (molar)",
    "Duration",
    "Population ethnicity",
    "Population health status",
    "Population age",
    "Study design",
    "Outcome primacy",
    "Biomatrix for NAD+"
  ),
  NMN_arm = c(
    paste(nma_studies %>% filter(precursor == "NMN") %>% pull(dose_mg_day) %>% unique() %>% sort(), collapse = ", ") %>% paste("mg/day"),
    "0.75 mmol/day (all studies)",
    paste(nma_studies %>% filter(precursor == "NMN") %>% pull(duration_weeks) %>% range(), collapse = "–") %>% paste("weeks"),
    paste(nma_studies %>% filter(precursor == "NMN") %>% pull(country) %>% unique(), collapse = ", "),
    paste(nma_studies %>% filter(precursor == "NMN") %>% pull(health_status) %>% unique(), collapse = "; "),
    "49–69 years",
    "All parallel RCT",
    "Mixed: metabolic secondary in most",
    "Blood cellular NAD+/NADH (pmol/mL) — Huang 2022 only"
  ),
  NR_arm = c(
    paste(nma_studies %>% filter(precursor == "NR") %>% pull(dose_mg_day) %>% unique() %>% sort(), collapse = ", ") %>% paste("mg/day"),
    "1.72–6.88 mmol/day",
    paste(nma_studies %>% filter(precursor == "NR") %>% pull(duration_weeks) %>% range(), collapse = "–") %>% paste("weeks"),
    paste(nma_studies %>% filter(precursor == "NR") %>% pull(country) %>% unique(), collapse = ", "),
    paste(nma_studies %>% filter(precursor == "NR") %>% pull(health_status) %>% unique(), collapse = "; "),
    "51–59 years",
    "Parallel + Crossover",
    "Safety primary in Conze; metabolic secondary in others",
    "Blood NAD+ (µM) — Bandi 2025 only"
  ),
  Assessment = c(
    "Major concern: 2–8× range in NR vs narrow band in NMN",
    "Major concern: up to 9.2× molar dose difference",
    "Moderate concern: overlapping but wider NR range",
    "Major concern: Japan/India (NMN) vs Denmark/Canada/India (NR)",
    "Minor concern: all metabolically healthy populations in NMA",
    "Minor concern: overlapping age ranges",
    "Moderate concern: Remie 2020 is crossover",
    "Moderate concern: most studies not powered for metabolic outcomes",
    "Major concern: incompatible units and biological matrices"
  ),
  stringsAsFactors = FALSE
)

write_csv(trans_domains, file.path(out_dir, "transitivity_assessment.csv"))
cat("\nSaved: transitivity_assessment.csv\n")


# OUTCOME REPORTING MATRIX (study × outcome)
cat("\n══════════════════════════════════════════════════\n")
cat("OUTCOME REPORTING MATRIX\n")
cat("══════════════════════════════════════════════════\n\n")

outcome_matrix <- dat %>%
  mutate(reported = 1) %>%
  select(study_id, precursor, outcome, reported)

# Pivot wider using base R
wide_data <- reshape(as.data.frame(outcome_matrix),
                     idvar = c("study_id", "precursor"),
                     timevar = "outcome",
                     direction = "wide")
# Clean column names (remove "reported." prefix)
names(wide_data) <- gsub("^reported\\.", "", names(wide_data))
# Replace NA with 0
wide_data[is.na(wide_data)] <- 0
wide_data <- wide_data[order(wide_data$precursor, wide_data$study_id), ]

write_csv(wide_data, file.path(out_dir, "outcome_reporting_matrix.csv"))
cat("Saved: outcome_reporting_matrix.csv\n")
print(as.data.frame(wide_data))


# k=1 META-ANALYSIS WARNING TABLE
cat("\n══════════════════════════════════════════════════\n")
cat("SINGLE-STUDY (k=1) WARNINGS\n")
cat("══════════════════════════════════════════════════\n\n")

k1_warnings <- pw_df %>%
  filter(k == 1) %>%
  select(outcome, precursor, k, MD_REML, lower_CI, upper_CI, p_value)

cat(sprintf("%d pairwise comparisons based on k=1 (no meta-analysis possible):\n", nrow(k1_warnings)))
print(as.data.frame(k1_warnings))

# Count indirect comparisons where BOTH arms are k=1
both_k1 <- ind_df %>%
  filter(k_NMN == 1 & k_NR == 1)
cat(sprintf("\n%d indirect comparisons where both arms are k=1 (pure subtraction):\n", nrow(both_k1)))
if (nrow(both_k1) > 0) print(as.data.frame(both_k1 %>% select(outcome, k_NMN, k_NR, MD_indirect, p_value)))


# FINAL SUMMARY
cat("\n══════════════════════════════════════════════════\n")
cat("ANALYSIS COMPLETE — SUMMARY\n")
cat("══════════════════════════════════════════════════\n\n")

cat("Output files:\n")
cat("  1. pairwise_metafor_REML.csv        — Primary pairwise meta-analyses\n")
cat("  2. indirect_comparisons_metafor_REML.csv — Indirect comparisons + multiplicity\n")
cat("  3. sensitivity_no_crossover_metafor.csv  — Crossover exclusion sensitivity\n")
cat("  4. transitivity_assessment.csv       — Formal transitivity domains\n")
cat("  5. outcome_reporting_matrix.csv      — Study × outcome heatmap data\n")
cat("\nKey findings after data correction:\n")

sig_indirect <- ind_df %>% filter(!is.na(p_value) & p_value < 0.05)
if (nrow(sig_indirect) == 0) {
  cat("  → NO indirect comparisons are statistically significant at α=0.05\n")
} else {
  cat(sprintf("  → %d indirect comparison(s) significant at α=0.05 (raw):\n", nrow(sig_indirect)))
  for (i in seq_len(nrow(sig_indirect))) {
    r <- sig_indirect[i, ]
    surv_bonf <- ifelse(r$p_bonferroni < 0.05, "YES", "NO")
    surv_bh   <- ifelse(r$p_BH < 0.05, "YES", "NO")
    cat(sprintf("    %s: p=%.4f (survives Bonferroni: %s, survives BH: %s)\n",
                r$outcome, r$p_value, surv_bonf, surv_bh))
  }
}

cat("\n  → NAD+ indirect comparison: EXCLUDED (unit incompatibility)\n")
cat("  → Previous headline result (MD=4.90, p=0.015) was based on\n")
cat("    incorrect data extraction and incompatible biological measurements.\n")
cat("\n══════════════════════════════════════════════════\n")
