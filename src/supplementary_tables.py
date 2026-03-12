#!/usr/bin/env python3
"""
Generate supplementary tables for manuscript.

Outputs:
  results/tables/supp_table_S1_search_strategy.csv
  results/tables/supp_table_S2_study_characteristics.csv
  results/tables/supp_table_S3_full_extraction.csv
  results/tables/supp_table_S4_rob2_detailed.csv
  results/tables/supp_table_S5_pairwise_all.csv
  results/tables/supp_table_S6_nma_all.csv
  results/tables/supp_table_S7_loo_summary.csv
  results/tables/supp_table_S8_grade_summary.csv
  results/tables/supp_table_S9_inter_rater_agreement.csv
"""

import csv, os
from collections import defaultdict

BASE = "/Users/tainguyen/Documents/NMN_vs_NR_SR_NMA"
DATA_DIR = os.path.join(BASE, "data", "extraction")
RES = os.path.join(BASE, "results", "tables")
os.makedirs(RES, exist_ok=True)


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved: {os.path.basename(path)}")


# -----------------------------------------------------------------------------
# S1: Search strategy summary
# -----------------------------------------------------------------------------

def make_S1():
    rows = [
        {"database": "PubMed", "records": 286,
         "search_query": '("nicotinamide mononucleotide" OR "NMN" OR "nicotinamide riboside" OR "NR"[tiab] OR "NIAGEN" OR "NAD+ precursor" OR "NAD precursor") AND ("randomized controlled trial"[pt] OR "randomized"[tiab] OR "placebo"[tiab] OR "RCT"[tiab]) AND ("insulin"[tiab] OR "glucose"[tiab] OR "HOMA"[tiab] OR "HbA1c"[tiab] OR "glycated hemoglobin"[tiab] OR "lipid"[tiab] OR "cholesterol"[tiab] OR "triglyceride"[tiab] OR "body weight"[tiab] OR "body composition"[tiab] OR "BMI"[tiab] OR "blood pressure"[tiab] OR "energy expenditure"[tiab] OR "NAD+"[tiab] OR "metabolic"[tiab] OR "fasting"[tiab] OR "C-reactive protein"[tiab] OR "interleukin"[tiab] OR "liver enzyme"[tiab] OR "ALT"[tiab] OR "AST"[tiab]); Filter: Humans, Clinical Trial or RCT',
         "filters": "Humans; Clinical Trial or RCT publication type",
         "date_searched": "2025-03-03"},
        {"database": "Embase (Ovid)", "records": 375,
         "search_query": '1: nicotinamide mononucleotide.mp. OR NMN.ti,ab. OR beta-nicotinamide mononucleotide.mp. OR MIB-626.mp. 2: nicotinamide riboside.mp. OR NIAGEN.mp. OR NR.ti,ab. 3: NAD+ precursor*.mp. OR NAD precursor*.mp. 4: 1 OR 2 OR 3; 5: insulin.ti,ab. OR glucose.ti,ab. OR HOMA*.ti,ab. OR HbA1c.ti,ab. OR glycated hemoglobin.ti,ab. OR lipid*.ti,ab. OR cholesterol.ti,ab. OR triglyceride*.ti,ab. OR body weight.ti,ab. OR body mass index.ti,ab. OR BMI.ti,ab. OR body composition.ti,ab. OR body fat.ti,ab. OR lean mass.ti,ab. OR energy expenditure.ti,ab. OR blood pressure.ti,ab. OR NAD+.ti,ab. OR C-reactive protein.ti,ab. OR CRP.ti,ab. OR interleukin*.ti,ab. OR ALT.ti,ab. OR AST.ti,ab. OR metabolic.ti,ab. 6: randomized controlled trial/ OR randomized.ti,ab. OR placebo.ti,ab. OR double-blind*.ti,ab. OR clinical trial/ OR crossover.ti,ab. 7: 4 AND 5 AND 6; 8: limit 7 to human; 9: remove duplicates from 8',
         "filters": "Humans; remove duplicates",
         "date_searched": "2025-03-03"},
        {"database": "Scopus", "records": 398,
         "search_query": 'TITLE-ABS-KEY(("nicotinamide mononucleotide" OR "NMN" OR "beta-nicotinamide mononucleotide" OR "MIB-626" OR "nicotinamide riboside" OR "NIAGEN" OR "NAD+ precursor" OR "NAD precursor") AND ("randomized" OR "randomised" OR "placebo" OR "RCT" OR "clinical trial") AND ("insulin" OR "glucose" OR "HOMA" OR "HbA1c" OR "lipid" OR "cholesterol" OR "triglyceride" OR "body weight" OR "body composition" OR "BMI" OR "blood pressure" OR "energy expenditure" OR "NAD+" OR "C-reactive protein" OR "metabolic" OR "ALT" OR "AST")); Limit: Document type = Article OR Review; Source type = Journal',
         "filters": "Article or Review; Journal source",
         "date_searched": "2025-03-03"},
        {"database": "Web of Science", "records": 319,
         "search_query": 'TS=(("nicotinamide mononucleotide" OR "NMN" OR "beta-nicotinamide mononucleotide" OR "MIB-626" OR "nicotinamide riboside" OR "NIAGEN" OR "NAD+ precursor" OR "NAD precursor") AND ("randomized" OR "randomised" OR "placebo" OR "RCT" OR "clinical trial") AND ("insulin" OR "glucose" OR "HOMA" OR "HbA1c" OR "lipid" OR "cholesterol" OR "triglyceride" OR "body weight" OR "body composition" OR "BMI" OR "blood pressure" OR "energy expenditure" OR "NAD+" OR "C-reactive protein" OR "metabolic" OR "ALT" OR "AST")); Refined by: Document Types = Article',
         "filters": "Article document type",
         "date_searched": "2025-03-03"},
        {"database": "Cochrane CENTRAL", "records": 309,
         "search_query": '#1: "nicotinamide mononucleotide" OR "NMN" OR "beta-nicotinamide mononucleotide" OR "MIB-626"; #2: "nicotinamide riboside" OR "NR" OR "NIAGEN"; #3: "NAD+ precursor" OR "NAD precursor"; #4: #1 OR #2 OR #3; #5: "insulin" OR "glucose" OR "HOMA-IR" OR "HbA1c" OR "glycated hemoglobin" OR "lipid" OR "cholesterol" OR "triglyceride" OR "body weight" OR "body mass index" OR "BMI" OR "body composition" OR "body fat" OR "lean mass" OR "energy expenditure" OR "blood pressure" OR "NAD+" OR "C-reactive protein" OR "CRP" OR "interleukin" OR "ALT" OR "AST" OR "metabolic"; #6: #4 AND #5; #7: Limit #6 to Trials',
         "filters": "Trials",
         "date_searched": "2025-03-03"},
        {"database": "TOTAL", "records": 1687,
         "search_query": "",
         "filters": "",
         "date_searched": ""},
        {"database": "After deduplication", "records": 1125,
         "search_query": "Duplicates removed: 562 (DOI, PMID, title fuzzy matching)",
         "filters": "",
         "date_searched": ""},
    ]
    write_csv(
        os.path.join(RES, "supp_table_S1_search_strategy.csv"),
        rows,
        ["database", "records", "search_query", "filters", "date_searched"],
    )


# -----------------------------------------------------------------------------
# S2: Study characteristics (from curated CSV)
# -----------------------------------------------------------------------------

def make_S2():
    src = os.path.join(DATA_DIR, "study_characteristics_curated.csv")
    with open(src) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # Select key columns for the supplementary table
    out_rows = []
    for r in rows:
        out_rows.append({
            "Study": r["study_id"].replace("_", " "),
            "Precursor": r["precursor"],
            "Country": r["country"],
            "Design": r["design"],
            "Population": r["population"],
            "N_randomized": r["n_randomized"],
            "N_intervention": r["n_analyzed_intervention"],
            "N_control": r["n_analyzed_control"],
            "Dose_mg_day": r["dose_mg_day"],
            "Duration_weeks": r["duration_weeks"],
            "Age_intervention": r["age_mean_intervention"],
            "Age_control": r["age_mean_control"],
            "Pct_female": r["pct_female"],
            "BMI_mean": r["bmi_mean"],
            "Health_status": r["health_status"],
            "Registration": r["registration"],
            "Funding": r["funding"],
            "RoB_concerns": r["risk_of_bias_concerns"],
        })

    write_csv(
        os.path.join(RES, "supp_table_S2_study_characteristics.csv"),
        out_rows,
        list(out_rows[0].keys()),
    )


# -----------------------------------------------------------------------------
# S3: Full data extraction (from nma_input_long.csv)
# -----------------------------------------------------------------------------

def make_S3():
    src = os.path.join(DATA_DIR, "nma_input_long.csv")
    with open(src) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    out_rows = []
    for r in rows:
        out_rows.append({
            "Study": r["study_id"].replace("_", " "),
            "Precursor": r["precursor"],
            "Outcome": r["outcome"],
            "N_treatment": r["n_treatment"],
            "N_control": r["n_control"],
            "Treatment_mean": r["treat_mean"],
            "Treatment_SD": r["treat_sd"],
            "Control_mean": r["ctrl_mean"],
            "Control_SD": r["ctrl_sd"],
            "MD": r["md"],
            "SE_MD": r["se_md"],
            "Unit": r["unit"],
            "Design": r["design"],
            "Original_dispersion": r["original_dispersion"],
            "Original_unit": r["original_unit"],
            "Notes": r.get("notes", ""),
        })

    write_csv(
        os.path.join(RES, "supp_table_S3_full_extraction.csv"),
        out_rows,
        list(out_rows[0].keys()),
    )


# -----------------------------------------------------------------------------
# S4: Detailed RoB 2 assessment
# -----------------------------------------------------------------------------

def make_S4():
    src = os.path.join(DATA_DIR, "rob2_assessment.csv")
    with open(src) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    out_rows = []
    for r in rows:
        out_rows.append({
            "Study": r["study_id"].replace("_", " "),
            "D1_Randomization": r["D1_randomization"],
            "D1_Justification": r["D1_justification"],
            "D2_Deviations": r["D2_deviations"],
            "D2_Justification": r["D2_justification"],
            "D3_Missing_data": r["D3_missing_data"],
            "D3_Justification": r["D3_justification"],
            "D4_Measurement": r["D4_measurement"],
            "D4_Justification": r["D4_justification"],
            "D5_Reporting": r["D5_reporting"],
            "D5_Justification": r["D5_justification"],
            "Overall": r["Overall"],
        })

    write_csv(
        os.path.join(RES, "supp_table_S4_rob2_detailed.csv"),
        out_rows,
        list(out_rows[0].keys()),
    )


# -----------------------------------------------------------------------------
# S5: All pairwise meta-analyses
# -----------------------------------------------------------------------------

def make_S5():
    src = os.path.join(RES, "pairwise_meta_analysis.csv")
    if not os.path.exists(src):
        print("  Skipped S5: pairwise_meta_analysis.csv not found")
        return
    with open(src) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    # Just copy with cleaner naming
    write_csv(
        os.path.join(RES, "supp_table_S5_pairwise_all.csv"),
        rows,
        list(rows[0].keys()),
    )


# -----------------------------------------------------------------------------
# S6: All NMA results
# -----------------------------------------------------------------------------

def make_S6():
    src = os.path.join(RES, "nma_results.csv")
    if not os.path.exists(src):
        print("  Skipped S6: nma_results.csv not found")
        return
    with open(src) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    write_csv(
        os.path.join(RES, "supp_table_S6_nma_all.csv"),
        rows,
        list(rows[0].keys()),
    )


# -----------------------------------------------------------------------------
# S7: LOO sensitivity summary
# -----------------------------------------------------------------------------

def make_S7():
    pw_src = os.path.join(RES, "sensitivity_loo_pairwise.csv")
    ind_src = os.path.join(RES, "sensitivity_loo_indirect.csv")

    summary_rows = []

    # Summarize pairwise LOO
    if os.path.exists(pw_src):
        with open(pw_src) as f:
            pw_rows = list(csv.DictReader(f))

        # Group by outcome + comparison
        groups = defaultdict(list)
        for r in pw_rows:
            key = (r["outcome"], r["comparison"])
            groups[key].append(r)

        for (oc, comp), rows in sorted(groups.items()):
            n_dir = sum(1 for r in rows if r["direction_change"] == "Yes")
            n_sig = sum(1 for r in rows if r["significance_change"] == "Yes")
            loo_mds = [float(r["loo_MD"]) for r in rows]
            summary_rows.append({
                "type": "Pairwise",
                "outcome": oc,
                "comparison": comp,
                "n_loo": len(rows),
                "full_MD": rows[0]["full_MD"],
                "full_p": rows[0]["full_p"],
                "loo_MD_range": f"{min(loo_mds):.4f} to {max(loo_mds):.4f}",
                "direction_changes": n_dir,
                "significance_changes": n_sig,
                "robust": "Yes" if n_dir == 0 and n_sig == 0 else "No",
            })

    # Summarize indirect LOO
    # Also load nma_results.csv to find all indirect outcomes (including
    # those with k=1 each arm where LOO is not possible)
    nma_src = os.path.join(RES, "nma_results.csv")
    all_indirect_outcomes = {}
    if os.path.exists(nma_src):
        with open(nma_src) as f:
            for r in csv.DictReader(f):
                if r["type"] == "indirect":
                    all_indirect_outcomes[r["outcome"]] = r

    if os.path.exists(ind_src):
        with open(ind_src) as f:
            ind_rows = list(csv.DictReader(f))

        groups = defaultdict(list)
        for r in ind_rows:
            groups[r["outcome"]].append(r)

        # Include all indirect outcomes, even those with no LOO entries
        for oc in sorted(all_indirect_outcomes.keys()):
            rows = groups.get(oc, [])
            nma_row = all_indirect_outcomes[oc]
            if rows:
                n_sig = sum(1 for r in rows if r["significance_change"] == "Yes")
                loo_mds = [float(r["loo_MD"]) for r in rows]
                note = ""
                # Flag narrow ranges where only one arm has k>=2
                k_str = nma_row.get("k", "")
                if "+" in k_str:
                    parts = k_str.split("+")
                    k_a, k_b = int(parts[0]), int(parts[1])
                    if k_a == 1 or k_b == 1:
                        note = f"Only {'NMN' if k_b == 1 else 'NR'} arm (k={'%d' % max(k_a,k_b)}) varied; other arm k=1"
                summary_rows.append({
                    "type": "Indirect (NMN vs NR)",
                    "outcome": oc,
                    "comparison": "NMN vs NR",
                    "n_loo": len(rows),
                    "full_MD": rows[0]["full_MD"],
                    "full_p": rows[0]["full_p"],
                    "loo_MD_range": f"{min(loo_mds):.4f} to {max(loo_mds):.4f}",
                    "direction_changes": "",
                    "significance_changes": n_sig,
                    "robust": "Yes" if n_sig == 0 else "No",
                    "note": note,
                })
            else:
                # Outcome has indirect comparison but no LOO possible
                summary_rows.append({
                    "type": "Indirect (NMN vs NR)",
                    "outcome": oc,
                    "comparison": "NMN vs NR",
                    "n_loo": 0,
                    "full_MD": f"{float(nma_row['MD']):.4f}",
                    "full_p": f"{float(nma_row['p_value']):.4f}",
                    "loo_MD_range": "—",
                    "direction_changes": "",
                    "significance_changes": "",
                    "robust": "—",
                    "note": "k=1 each arm; LOO not possible",
                })

    # Ensure note column exists in pairwise rows too
    for r in summary_rows:
        if "note" not in r:
            r["note"] = ""

    if summary_rows:
        fieldnames = list(summary_rows[0].keys())
        if "note" not in fieldnames:
            fieldnames.append("note")
        write_csv(
            os.path.join(RES, "supp_table_S7_loo_summary.csv"),
            summary_rows,
            fieldnames,
        )


# -----------------------------------------------------------------------------
# S8: GRADE/CINeMA summary
# -----------------------------------------------------------------------------

def make_S8():
    src = os.path.join(RES, "grade_cinema_assessment.csv")
    if not os.path.exists(src):
        print("  Skipped S8: grade_cinema_assessment.csv not found")
        return

    with open(src) as f:
        rows = list(csv.DictReader(f))

    out_rows = []
    for r in rows:
        out_rows.append({
            "Outcome": r["outcome_label"],
            "k_NMN": r["k_NMN"],
            "k_NR": r["k_NR"],
            "MD_95CI": f"{r['MD']} [{r['CI_lower']}, {r['CI_upper']}]",
            "p_value": r["p_value"],
            "Within_study_bias": r["Within_study_bias_rating"],
            "Reporting_bias": r["Reporting_bias_rating"],
            "Indirectness": r["Indirectness_rating"],
            "Imprecision": r["Imprecision_rating"],
            "Heterogeneity": r["Heterogeneity_rating"],
            "Incoherence": r["Incoherence_rating"],
            "Overall_certainty": r["overall_certainty"],
            "Downgrades": r["n_downgrades"],
        })

    write_csv(
        os.path.join(RES, "supp_table_S8_grade_summary.csv"),
        out_rows,
        list(out_rows[0].keys()),
    )


# -----------------------------------------------------------------------------
# S9: Inter-rater agreement
# -----------------------------------------------------------------------------

def make_S9():
    src = os.path.join(DATA_DIR, "inter_rater_agreement.csv")
    if not os.path.exists(src):
        print("  Skipped S9: inter_rater_agreement.csv not found")
        return

    with open(src) as f:
        rows = list(csv.DictReader(f))

    out_rows = []
    for r in rows:
        out_rows.append({
            "Phase": r["phase"],
            "Items_assessed": r["n_items"],
            "Agreements": r["agreements"],
            "Percent_agreement": r["pct_agreement"],
            "Cohen_kappa": r["kappa"],
            "Disagreement_resolution": r["resolution"],
        })

    write_csv(
        os.path.join(RES, "supp_table_S9_inter_rater_agreement.csv"),
        out_rows,
        list(out_rows[0].keys()),
    )


# -----------------------------------------------------------------------------
# Main execution
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("  GENERATING SUPPLEMENTARY TABLES")
    print("=" * 60)

    make_S1(); print()
    make_S2(); print()
    make_S3(); print()
    make_S4(); print()
    make_S5(); print()
    make_S6(); print()
    make_S7(); print()
    make_S8(); print()
    make_S9(); print()

    print("=" * 60)
    print("  All supplementary tables generated.")
    print("=" * 60)
