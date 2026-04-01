#!/usr/bin/env python3
"""
Data Verification Report: Cross-check extracted data against source papers.

Compares nma_input_long.csv values to verified values from published fulltexts.
Outputs a verification report documenting all discrepancies.
"""

import csv, math, os

BASE = "/Users/tainguyen/Documents/NMN_vs_NR_SR_NMA"
DATA = os.path.join(BASE, "data", "extraction", "nma_input_long.csv")
RES = os.path.join(BASE, "results", "tables")
os.makedirs(RES, exist_ok=True)


# Verified values from published papers (manually checked against full texts)
# Format: (study_id, outcome) -> {
#   "paper_treat_mean", "paper_treat_sd_or_sem", "paper_ctrl_mean",
#   "paper_ctrl_sd_or_sem", "paper_unit", "paper_dispersion_type",
#   "paper_table", "notes"
# }

VERIFIED = {
    # Huang 2022 (Frontiers in Aging, Table 2–8, 10–13) ───────────
    # Paper reports Day 60 end-of-study values; dispersion is SD unless noted
    ("Huang_2022", "NAD+"): {
        "paper_treat_mean": 9.07, "paper_treat_disp": 5.65,
        "paper_ctrl_mean": 8.14, "paper_ctrl_disp": 4.86,
        "paper_unit": "pmol/mL", "disp_type": "SD",
        "paper_table": "Table 2",
        "paper_p_between": 0.40,  # change-from-baseline between-group p
        "notes": "Blood cellular NAD+/NADH. NOT ng/mL. Paper p=0.40 (NS). "
                 "Change from baseline: Uthever 2.50±8.21 (p=0.10), "
                 "Placebo 1.01±5.35 (p=0.30). Between-group change p=0.40.",
        "severity": "CRITICAL",
    },
    ("Huang_2022", "HOMA-IR"): {
        "paper_treat_mean": 1.79, "paper_treat_disp": 0.94,
        "paper_ctrl_mean": 2.39, "paper_ctrl_disp": 2.08,
        "paper_unit": "", "disp_type": "SD",
        "paper_table": "Table 6",
        "paper_p_between": 0.18,
        "notes": "Extraction has 3.6 and 4.6; paper shows 1.79 and 2.39. "
                 "Between-group change p=0.18 (NS).",
        "severity": "CRITICAL",
    },
    ("Huang_2022", "FBG"): {
        "paper_treat_mean": 91.53, "paper_treat_disp": 21.41,
        "paper_ctrl_mean": 108.48, "paper_ctrl_disp": 47.12,
        "paper_unit": "mg/dL", "disp_type": "SD",
        "paper_table": "Table 7",
        "paper_p_between": 0.21,
        "notes": "Paper shows 91.53±21.41 vs 108.48±47.12; extraction has "
                 "96.5±10.2 vs 105.4±15.5. Neither values nor SDs match.",
        "severity": "CRITICAL",
    },
    ("Huang_2022", "fasting_insulin"): {
        "paper_treat_mean": 14.06, "paper_treat_disp": 7.52,
        "paper_ctrl_mean": 19.08, "paper_ctrl_disp": 18.09,
        "paper_unit": "µU/mL", "disp_type": "SD",
        "paper_table": "Table 8",
        "paper_p_between": 0.22,
        "notes": "Extraction SDs wildly wrong: 3.3 vs 7.52 (treat), "
                 "4.1 vs 18.09 (ctrl).",
        "severity": "CRITICAL",
    },
    ("Huang_2022", "SBP"): {
        "paper_treat_mean": 124.19, "paper_treat_disp": 12.03,
        "paper_ctrl_mean": 126.26, "paper_ctrl_disp": 13.51,
        "paper_unit": "mmHg", "disp_type": "SD",
        "paper_table": "Table 5",
        "paper_p_between": 0.42,
        "notes": "Extraction: 128.3 vs 133.2; paper: 124.19 vs 126.26.",
        "severity": "CRITICAL",
    },
    ("Huang_2022", "DBP"): {
        "paper_treat_mean": 76.81, "paper_treat_disp": 8.30,
        "paper_ctrl_mean": 77.84, "paper_ctrl_disp": 9.01,
        "paper_unit": "mmHg", "disp_type": "SD",
        "paper_table": "Table 5",
        "paper_p_between": 0.79,
        "notes": "Extraction: 78.1 vs 82.5; paper: 76.81 vs 77.84.",
        "severity": "CRITICAL",
    },
    ("Huang_2022", "TC"): {
        "paper_treat_mean": 172.86, "paper_treat_disp": 31.89,
        "paper_ctrl_mean": 179.79, "paper_ctrl_disp": 39.40,
        "paper_unit": "mg/dL", "disp_type": "SD",
        "paper_table": "Table 13",
        "paper_p_between": 0.90,
        "notes": "Extraction: 204.1 vs 214.1; paper: 172.86 vs 179.79.",
        "severity": "CRITICAL",
    },
    ("Huang_2022", "TG"): {
        "paper_treat_mean": 132.48, "paper_treat_disp": 60.41,
        "paper_ctrl_mean": 167.91, "paper_ctrl_disp": 198.89,
        "paper_unit": "mg/dL", "disp_type": "SD",
        "paper_table": "Table 13",
        "paper_p_between": 0.75,
        "notes": "Extraction: 165.4±42.2 vs 183.1±46.1; paper: "
                 "132.48±60.41 vs 167.91±198.89. Ctrl SD off by 4x.",
        "severity": "CRITICAL",
    },
    ("Huang_2022", "LDL"): {
        "paper_treat_mean": 114.45, "paper_treat_disp": 30.75,
        "paper_ctrl_mean": 119.42, "paper_ctrl_disp": 37.76,
        "paper_unit": "mg/dL", "disp_type": "SD",
        "paper_table": "Table 13",
        "paper_p_between": 0.39,
        "notes": "Extraction: 100.9 vs 110.2; paper: 114.45 vs 119.42.",
        "severity": "CRITICAL",
    },
    ("Huang_2022", "ALT"): {
        "paper_treat_mean": 17.20, "paper_treat_disp": 6.58,
        "paper_ctrl_mean": 19.45, "paper_ctrl_disp": 8.77,
        "paper_unit": "U/L", "disp_type": "SD",
        "paper_table": "Table 11",
        "paper_p_between": 0.06,
        "notes": "Extraction: 19.2±3.3 vs 22.3±6.1; paper: 17.20±6.58 "
                 "vs 19.45±8.77.",
        "severity": "MAJOR",
    },
    ("Huang_2022", "AST"): {
        "paper_treat_mean": 18.31, "paper_treat_disp": 4.09,
        "paper_ctrl_mean": 18.97, "paper_ctrl_disp": 4.11,
        "paper_unit": "U/L", "disp_type": "SD",
        "paper_table": "Table 11",
        "paper_p_between": 0.53,
        "notes": "Extraction: 22.6±4.2 vs 24.8±4.8; paper: 18.31±4.09 "
                 "vs 18.97±4.11.",
        "severity": "MAJOR",
    },

    # Bandi 2025 (Table 2) ────────────────────────────────────────
    ("Bandi_2025", "NAD+"): {
        "paper_treat_mean": 31.9, "paper_treat_disp": 6.5,
        "paper_ctrl_mean": 25.8, "paper_ctrl_disp": 7.1,
        "paper_unit": "µM", "disp_type": "SD",
        "paper_table": "Table 2",
        "notes": "NR-500 Day 60 values. Unit is µM (micromolar). "
                 "Extraction lists unit as 'uM' which matches. "
                 "Values match paper. CORRECT extraction.",
        "severity": "OK_BUT_UNIT_INCOMPATIBLE",
    },

    # Yoshino 2021 (Science, Table 1) ─────────────────────────────
    # Paper reports values as SEM; glucose in mmol/L, TG/HDL in mmol/L
    ("Yoshino_2021", "FBG"): {
        "paper_treat_mean_mmol": 5.6, "paper_treat_sem": 0.2,
        "paper_ctrl_mean_mmol": 5.6, "paper_ctrl_sem": 0.2,
        "paper_unit": "mmol/L (convert ×18.018 to mg/dL)",
        "disp_type": "SEM",
        "paper_table": "Table 1",
        "notes": "After values: NMN 5.6±0.2, Placebo 5.6±0.2 mmol/L (SEM). "
                 "Both groups = 100.9 mg/dL. MD should be ~0, not 2.0. "
                 "Extraction has 101.0 vs 99.0 — may be from supplementary "
                 "data with more decimal precision, but p-values show NS.",
        "severity": "NEEDS_VERIFICATION",
    },
}


def read_extraction():
    """Read the current extraction data."""
    rows = []
    with open(DATA) as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def verify_all():
    """Cross-check extracted values against verified paper values."""
    extraction = read_extraction()
    
    report_rows = []
    
    # Build lookup from extraction
    ext_lookup = {}
    for r in extraction:
        key = (r["study_id"], r["outcome"])
        ext_lookup[key] = r
    
    print("=" * 80)
    print("DATA VERIFICATION REPORT")
    print("=" * 80)
    print()
    
    critical_count = 0
    major_count = 0
    ok_count = 0
    
    for (study, outcome), verified in sorted(VERIFIED.items()):
        ext = ext_lookup.get((study, outcome))
        if ext is None:
            print(f"  WARNING: {study} / {outcome} not found in extraction")
            continue
        
        sev = verified.get("severity", "UNKNOWN")
        
        row = {
            "study_id": study,
            "outcome": outcome,
            "severity": sev,
            "paper_table": verified.get("paper_table", ""),
            "paper_unit": verified.get("paper_unit", ""),
            "extraction_unit": ext.get("unit", ""),
        }
        
        # Compare means if available
        if "paper_treat_mean" in verified:
            ext_t = float(ext["treat_mean"])
            ext_c = float(ext["ctrl_mean"])
            paper_t = verified["paper_treat_mean"]
            paper_c = verified["paper_ctrl_mean"]
            
            row["ext_treat_mean"] = ext_t
            row["paper_treat_mean"] = paper_t
            row["treat_mean_match"] = abs(ext_t - paper_t) < 0.5
            
            row["ext_ctrl_mean"] = ext_c
            row["paper_ctrl_mean"] = paper_c
            row["ctrl_mean_match"] = abs(ext_c - paper_c) < 0.5
            
            # Compare SDs
            ext_t_sd = float(ext["treat_sd"])
            ext_c_sd = float(ext["ctrl_sd"])
            paper_t_disp = verified["paper_treat_disp"]
            paper_c_disp = verified["paper_ctrl_disp"]
            
            row["ext_treat_sd"] = ext_t_sd
            row["paper_treat_disp"] = paper_t_disp
            row["ext_ctrl_sd"] = ext_c_sd
            row["paper_ctrl_disp"] = paper_c_disp
            
        row["notes"] = verified.get("notes", "")
        
        # Print summary
        if sev == "CRITICAL":
            critical_count += 1
            marker = "❌ CRITICAL"
        elif sev == "MAJOR":
            major_count += 1
            marker = "⚠️  MAJOR"
        elif "OK" in sev:
            ok_count += 1
            marker = "✅ OK"
        else:
            marker = "❓ " + sev
        
        print(f"  {marker}: {study} / {outcome}")
        if "paper_treat_mean" in verified:
            ext_t = float(ext["treat_mean"])
            ext_c = float(ext["ctrl_mean"])
            paper_t = verified["paper_treat_mean"]
            paper_c = verified["paper_ctrl_mean"]
            print(f"    Extracted:  treat={ext_t}, ctrl={ext_c}, "
                  f"unit={ext.get('unit', 'N/A')}")
            print(f"    Paper:      treat={paper_t}, ctrl={paper_c}, "
                  f"unit={verified['paper_unit']}")
            if not (abs(ext_t - paper_t) < 0.5 and abs(ext_c - paper_c) < 0.5):
                print(f"    *** VALUES DO NOT MATCH ***")
        
        if verified.get("notes"):
            print(f"    Notes: {verified['notes']}")
        print()
        
        report_rows.append(row)
    
    # Summary
    print("=" * 80)
    print(f"SUMMARY: {critical_count} CRITICAL, {major_count} MAJOR, "
          f"{ok_count} OK/compatible")
    print("=" * 80)
    print()
    
    # NAD+ unit incompatibility analysis
    print("=" * 80)
    print("NAD+ UNIT INCOMPATIBILITY ANALYSIS")
    print("=" * 80)
    print()
    print("  Huang 2022 (NMN arm): NAD+/NADH in pmol/mL")
    print("    Paper values: Uthever Day 60 = 9.07 ± 5.65 pmol/mL")
    print("    Placebo Day 60 = 8.14 ± 4.86 pmol/mL")
    print("    Between-group change from baseline: p = 0.40 (NOT SIGNIFICANT)")
    print()
    print("  Bandi 2025 (NR arm): NAD+ in µM (micromolar)")
    print("    Paper values: NR-500 Day 60 = 31.9 ± 6.5 µM")
    print("    Placebo Day 60 = 25.8 ± 7.1 µM")
    print()
    print("  CONVERSION:")
    print("    1 µM = 1,000 pmol/mL")
    print("    Huang MD in µM: (9.07 - 8.14) / 1000 = 0.00093 µM")
    print("    Bandi MD in µM: 31.9 - 25.8 = 6.1 µM")
    print()
    print("  The current indirect comparison computes 11.0 - 6.1 = 4.9")
    print("  using incompatible units (the extraction says ng/mL for Huang,")
    print("  but the paper says pmol/mL, and neither is µM).")
    print()
    print("  Even if converted to common units, the measurements are from")
    print("  different biological matrices (blood cellular vs blood/serum)")
    print("  using different assays — making comparison scientifically invalid.")
    print()
    print("  CONCLUSION: The NAD+ indirect comparison CANNOT be performed.")
    print("  The only statistically significant result in the paper does not exist.")
    print()
    
    # Huang 2022 comprehensive error summary
    print("=" * 80)
    print("HUANG 2022 EXTRACTION ERROR SUMMARY")
    print("=" * 80)
    print()
    print("  ALL 11 outcomes extracted from Huang 2022 have values that DO NOT")
    print("  match the published paper tables. The extracted values (treat_mean,")
    print("  ctrl_mean, treat_sd, ctrl_sd) cannot be found in any table of the")
    print("  published Frontiers in Aging article.")
    print()
    print("  This affects the following pairwise meta-analyses:")
    print("    - FBG (NMN arm, 1 of 4 studies)")
    print("    - HOMA-IR (NMN arm, sole study → entire indirect comparison)")
    print("    - Fasting insulin (NMN arm, 1 of 2 studies)")  
    print("    - TC (NMN arm, 1 of 3 studies)")
    print("    - TG (NMN arm, 1 of 3 studies)")
    print("    - LDL (NMN arm, 1 of 3 studies)")
    print("    - SBP (NMN arm, 1 of 4 studies)")
    print("    - DBP (NMN arm, 1 of 4 studies)")
    print("    - ALT (NMN arm, 1 of 3 studies)")
    print("    - AST (NMN arm, 1 of 2 studies)")
    print("    - NAD+ (NMN arm, sole study → entire indirect comparison)")
    print()
    print("  REQUIRED ACTION: Complete re-extraction from Huang 2022 original")
    print("  paper (Frontiers in Aging 3:851698, doi:10.3389/fragi.2022.851698)")
    print("  using the verified values documented in this script.")
    
    # Write CSV report
    fieldnames = [
        "study_id", "outcome", "severity", "paper_table", "paper_unit",
        "extraction_unit", "ext_treat_mean", "paper_treat_mean",
        "treat_mean_match", "ext_ctrl_mean", "paper_ctrl_mean",
        "ctrl_mean_match", "ext_treat_sd", "paper_treat_disp",
        "ext_ctrl_sd", "paper_ctrl_disp", "notes",
    ]
    out_path = os.path.join(RES, "data_verification_report.csv")
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(report_rows)
    print(f"\n  Report saved: {out_path}")


if __name__ == "__main__":
    verify_all()
