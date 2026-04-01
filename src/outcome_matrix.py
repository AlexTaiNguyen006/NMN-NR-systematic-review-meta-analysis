#!/usr/bin/env python3
"""
Outcome reporting matrix heatmap: which studies report which outcomes.
Produces a study × outcome grid showing data availability and sparsity.
"""

import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

BASE = "/Users/tainguyen/Documents/NMN_vs_NR_SR_NMA"
DATA = os.path.join(BASE, "data", "extraction", "nma_input_long.csv")
FIGS = os.path.join(BASE, "results", "figures")
os.makedirs(FIGS, exist_ok=True)

# --- Style ---
FONT_FAMILY = "Arial"
DPI = 300
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": [FONT_FAMILY, "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "figure.dpi": DPI,
    "savefig.dpi": DPI,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.2,
})

COL_NMN = "#2196F3"
COL_NR  = "#4CAF50"

OUTCOME_ORDER = [
    "FBG", "HbA1c", "HOMA-IR", "fasting_insulin",
    "TC", "LDL", "HDL", "TG",
    "body_weight", "BMI", "fat_mass", "body_fat_pct", "skeletal_muscle",
    "SBP", "DBP",
    "ALT", "AST",
    "NAD+", "IL-6",
    "FFM_pct", "sleeping_metabolic_rate",
]

OUTCOME_LABELS = {
    "FBG": "FBG", "HbA1c": "HbA1c", "HOMA-IR": "HOMA-IR",
    "fasting_insulin": "Fasting Insulin",
    "TC": "Total Cholesterol", "LDL": "LDL-C", "HDL": "HDL-C",
    "TG": "Triglycerides",
    "body_weight": "Body Weight", "BMI": "BMI",
    "fat_mass": "Fat Mass", "body_fat_pct": "Body Fat %",
    "skeletal_muscle": "Skeletal Muscle",
    "SBP": "Systolic BP", "DBP": "Diastolic BP",
    "ALT": "ALT", "AST": "AST",
    "NAD+": "Blood NAD+", "IL-6": "IL-6",
    "FFM_pct": "Fat-Free Mass %",
    "sleeping_metabolic_rate": "Sleeping MR",
}


def main():
    # Read data
    rows = []
    with open(DATA) as f:
        for r in csv.DictReader(f):
            rows.append(r)
    
    # Get unique studies and outcomes present in data
    studies = []
    study_precursors = {}
    for r in rows:
        sid = r["study_id"]
        if sid not in study_precursors:
            studies.append(sid)
            study_precursors[sid] = r["precursor"]
    
    # Sort: NMN first, then NR, alphabetically within each
    nmn_studies = sorted([s for s in studies if study_precursors[s] == "NMN"])
    nr_studies = sorted([s for s in studies if study_precursors[s] == "NR"])
    studies_sorted = nmn_studies + nr_studies
    
    # Get outcomes present in data, ordered
    outcomes_in_data = set(r["outcome"] for r in rows)
    outcomes_sorted = [o for o in OUTCOME_ORDER if o in outcomes_in_data]
    # Add any not in our order list
    for o in sorted(outcomes_in_data):
        if o not in outcomes_sorted:
            outcomes_sorted.append(o)
    
    # Build matrix
    n_studies = len(studies_sorted)
    n_outcomes = len(outcomes_sorted)
    matrix = np.zeros((n_studies, n_outcomes))
    
    # Build lookup
    reported = set()
    for r in rows:
        reported.add((r["study_id"], r["outcome"]))
    
    for i, s in enumerate(studies_sorted):
        for j, o in enumerate(outcomes_sorted):
            if (s, o) in reported:
                matrix[i, j] = 1 if study_precursors[s] == "NMN" else 2
    
    # --- Plot ---
    fig_h = max(4, 0.4 * n_studies + 1.5)
    fig_w = max(8, 0.45 * n_outcomes + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    
    # Custom colormap: 0=white, 1=blue(NMN), 2=green(NR)
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#FFFFFF", COL_NMN, COL_NR])
    
    ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=0, vmax=2,
              interpolation="nearest")
    
    # Grid lines
    for i in range(n_studies + 1):
        ax.axhline(i - 0.5, color="#cccccc", linewidth=0.5)
    for j in range(n_outcomes + 1):
        ax.axvline(j - 0.5, color="#cccccc", linewidth=0.5)
    
    # Separator between NMN and NR groups
    sep_y = len(nmn_studies) - 0.5
    ax.axhline(sep_y, color="#333333", linewidth=2)
    
    # Labels
    study_labels = [s.replace("_", " ") for s in studies_sorted]
    outcome_labels = [OUTCOME_LABELS.get(o, o) for o in outcomes_sorted]
    
    ax.set_yticks(range(n_studies))
    ax.set_yticklabels(study_labels, fontsize=8)
    ax.set_xticks(range(n_outcomes))
    ax.set_xticklabels(outcome_labels, rotation=45, ha="right", fontsize=7.5)
    
    # Color y-axis labels by precursor
    for i, s in enumerate(studies_sorted):
        color = COL_NMN if study_precursors[s] == "NMN" else COL_NR
        ax.get_yticklabels()[i].set_color(color)
        ax.get_yticklabels()[i].set_fontweight("bold")
    
    ax.set_title("Outcome Reporting Matrix: Study × Outcome Availability",
                 fontsize=11, fontweight="bold", pad=12)
    
    # Annotate column counts at bottom
    for j in range(n_outcomes):
        col_sum = int(np.sum(matrix[:, j] > 0))
        ax.text(j, n_studies - 0.1, f"k={col_sum}", ha="center", va="top",
                fontsize=6.5, color="#555555")
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COL_NMN, edgecolor="#333", label="NMN study"),
        mpatches.Patch(facecolor=COL_NR, edgecolor="#333", label="NR study"),
        mpatches.Patch(facecolor="white", edgecolor="#ccc", label="Not reported"),
    ]
    ax.legend(handles=legend_elements, loc="upper right",
              bbox_to_anchor=(1.0, -0.12), ncol=3, fontsize=8,
              frameon=True, edgecolor="#cccccc")
    
    # Category annotations along top
    categories = [
        (0, 3, "Glycemic"),
        (4, 7, "Lipids"),
        (8, 12, "Body Composition"),
        (13, 14, "BP"),
        (15, 16, "Liver"),
    ]
    for start, end, label in categories:
        if end < n_outcomes:
            mid = (start + end) / 2
            ax.text(mid, -1.2, label, ha="center", va="bottom",
                    fontsize=7, fontweight="bold", color="#555555")
    
    plt.tight_layout()
    out_path = os.path.join(FIGS, "outcome_reporting_matrix.png")
    plt.savefig(out_path, dpi=DPI)
    plt.close()
    print(f"Saved: {out_path}")
    
    # Print sparsity summary
    total_cells = n_studies * n_outcomes
    filled = int(np.sum(matrix > 0))
    print(f"\nMatrix dimensions: {n_studies} studies × {n_outcomes} outcomes")
    print(f"Cells filled: {filled}/{total_cells} ({100*filled/total_cells:.1f}%)")
    print(f"Sparsity: {100*(1-filled/total_cells):.1f}%")
    
    # Count outcomes with indirect comparison possible (both NMN and NR)
    indirect_possible = 0
    for j, o in enumerate(outcomes_sorted):
        has_nmn = any(matrix[i, j] == 1 for i in range(n_studies))
        has_nr = any(matrix[i, j] == 2 for i in range(n_studies))
        if has_nmn and has_nr:
            indirect_possible += 1
    print(f"Outcomes with NMN + NR data (indirect possible): {indirect_possible}/{n_outcomes}")


if __name__ == "__main__":
    main()
