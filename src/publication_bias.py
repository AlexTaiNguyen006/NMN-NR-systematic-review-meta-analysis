#!/usr/bin/env python3
"""
Publication bias assessment for NMN vs NR NMA.

Generates funnel plots for all pairwise comparisons with k >= 3 and
runs Egger's regression test where applicable.

NOTE: No pairwise comparison exceeds k = 4. The Cochrane Handbook
recommends k >= 10 for reliable funnel plot interpretation and Egger's
test. These outputs are provided for transparency and completeness;
they should not be over-interpreted.

Outputs:
  results/figures/funnel_<outcome>_<precursor>_vs_PBO.png
  results/tables/publication_bias_assessment.csv
"""

import csv, math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, linregress

# Paths
BASE = "/Users/tainguyen/Documents/NMN_vs_NR_SR_NMA"
DATA = os.path.join(BASE, "data", "extraction", "nma_input_long.csv")
FIGS = os.path.join(BASE, "results", "figures")
RES  = os.path.join(BASE, "results", "tables")
os.makedirs(FIGS, exist_ok=True)
os.makedirs(RES, exist_ok=True)

# Style (matches regenerate_figures.py)
FONT_FAMILY     = "Arial"
TITLE_SIZE      = 12
AXIS_LABEL_SIZE = 10
TICK_LABEL_SIZE = 9
ANNOT_SIZE      = 8.5
DPI             = 300

COL_NMN  = "#2196F3"
COL_NR   = "#4CAF50"
COL_FILL = "#E3F2FD"

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   [FONT_FAMILY, "Helvetica", "DejaVu Sans"],
    "font.size":         TICK_LABEL_SIZE,
    "axes.titlesize":    TITLE_SIZE,
    "axes.titleweight":  "bold",
    "axes.labelsize":    AXIS_LABEL_SIZE,
    "xtick.labelsize":   TICK_LABEL_SIZE,
    "ytick.labelsize":   TICK_LABEL_SIZE,
    "figure.dpi":        DPI,
    "savefig.dpi":       DPI,
    "savefig.bbox":      "tight",
    "savefig.pad_inches": 0.15,
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

OUTCOME_LABELS = {
    "FBG": "Fasting Blood Glucose", "HbA1c": "HbA1c", "HOMA-IR": "HOMA-IR",
    "fasting_insulin": "Fasting Insulin", "TC": "Total Cholesterol",
    "LDL": "LDL Cholesterol", "HDL": "HDL Cholesterol", "TG": "Triglycerides",
    "body_weight": "Body Weight", "BMI": "BMI", "SBP": "Systolic BP",
    "DBP": "Diastolic BP", "ALT": "ALT", "AST": "AST", "IL-6": "IL-6",
    "NAD+": "Blood NAD+",
}

OUTCOME_UNITS = {
    "FBG": "mg/dL", "HbA1c": "%", "HOMA-IR": "", "fasting_insulin": "µU/mL",
    "TC": "mg/dL", "LDL": "mg/dL", "HDL": "mg/dL", "TG": "mg/dL",
    "body_weight": "kg", "BMI": "kg/m²", "SBP": "mmHg", "DBP": "mmHg",
    "ALT": "U/L", "AST": "U/L", "IL-6": "pg/mL",
}

MIN_K = 3  # minimum studies for funnel plot


def read_data(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            r["md"]    = float(r["md"])
            r["se_md"] = float(r["se_md"])
            rows.append(r)
    return rows


def pairwise_ma(studies):
    """DerSimonian-Laird random effects, returns summary estimate."""
    k = len(studies)
    tes = [s["md"] for s in studies]
    ses = [s["se_md"] for s in studies]
    ws  = [1.0 / (se**2) for se in ses]
    W   = sum(ws)
    te_f = sum(w * te for w, te in zip(ws, tes)) / W
    Q    = sum(w * (te - te_f)**2 for w, te in zip(ws, tes))
    df   = k - 1
    C    = W - sum(w**2 for w in ws) / W
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    ws_r = [1.0 / (se**2 + tau2) for se in ses]
    W_r  = sum(ws_r)
    te_r = sum(w * te for w, te in zip(ws_r, tes)) / W_r
    se_r = 1.0 / math.sqrt(W_r)
    return te_r, se_r, tau2


def egger_test(tes, ses):
    """
    Egger's regression test for funnel plot asymmetry.
    Regresses standardized effect (te/se) on precision (1/se).
    Returns intercept, SE of intercept, t-statistic, p-value.
    Returns None if k < 3.
    """
    k = len(tes)
    if k < 3:
        return None
    precision = [1.0 / se for se in ses]
    std_effect = [te / se for te, se in zip(tes, ses)]
    slope, intercept, r_value, p_value, std_err = linregress(precision, std_effect)
    # For Egger's, the intercept test is what matters
    # The linregress p_value is for the slope; we need the intercept test
    # Recompute using the intercept's SE
    n = len(precision)
    x = np.array(precision)
    y = np.array(std_effect)
    y_pred = slope * x + intercept
    residuals = y - y_pred
    mse = np.sum(residuals**2) / (n - 2) if n > 2 else np.nan
    x_mean = np.mean(x)
    ss_x = np.sum((x - x_mean)**2)
    se_intercept = math.sqrt(mse * (1.0/n + x_mean**2 / ss_x)) if ss_x > 0 and n > 2 else np.nan
    if np.isnan(se_intercept) or se_intercept == 0:
        return None
    from scipy.stats import t as t_dist
    t_stat = intercept / se_intercept
    df = n - 2
    p_val = 2 * t_dist.sf(abs(t_stat), df)
    return {
        "intercept": intercept,
        "se_intercept": se_intercept,
        "t_stat": t_stat,
        "df": df,
        "p_value": p_val,
    }


def draw_funnel(tes, ses, summary_te, outcome, precursor, k, egger_result, filename):
    """Standard funnel plot: x = effect size (MD), y = SE (inverted)."""
    fig, ax = plt.subplots(figsize=(6, 5))

    color = COL_NMN if precursor == "NMN" else COL_NR

    # Pseudo-CI region (triangular funnel)
    se_max = max(ses) * 1.3
    se_range = np.linspace(0, se_max, 200)
    ax.fill_betweenx(se_range,
                     summary_te - 1.96 * se_range,
                     summary_te + 1.96 * se_range,
                     alpha=0.12, color=color, label="Pseudo 95% CI")

    # Individual studies
    ax.scatter(tes, ses, s=60, c=color, edgecolors="white", linewidths=0.8,
              zorder=5)

    # Summary line
    ax.axvline(summary_te, color=color, linestyle="--", linewidth=1, alpha=0.7,
              label=f"Summary MD = {summary_te:.2f}")

    # Zero line
    ax.axvline(0, color="#888888", linestyle=":", linewidth=0.8, alpha=0.5)

    # Invert y-axis (larger SE at bottom)
    ax.invert_yaxis()
    ax.set_ylim(se_max, 0)

    unit = OUTCOME_UNITS.get(outcome, "")
    unit_str = f" ({unit})" if unit else ""
    pretty = OUTCOME_LABELS.get(outcome, outcome)
    ax.set_xlabel(f"Mean Difference{unit_str}", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel("Standard Error", fontsize=AXIS_LABEL_SIZE)
    ax.set_title(f"Funnel Plot: {pretty}\n{precursor} vs Placebo (k={k})",
                fontsize=TITLE_SIZE, fontweight="bold")

    # Caveat annotation
    caveat = f"Caution: k = {k} (recommended ≥ 10 for interpretation)"
    if egger_result:
        caveat += f"\nEgger's test: intercept = {egger_result['intercept']:.2f}, p = {egger_result['p_value']:.3f}"
    ax.annotate(caveat, xy=(0.02, 0.98), xycoords="axes fraction",
               fontsize=7.5, va="top", ha="left", style="italic",
               color="#666666",
               bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow",
                        edgecolor="#cccccc", alpha=0.9))

    ax.legend(fontsize=8, loc="lower right")

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGS, f"{filename}.{ext}"),
                   dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def main():
    data = read_data(DATA)

    # Group by outcome + precursor
    from collections import defaultdict
    groups = defaultdict(list)
    for r in data:
        groups[(r["outcome"], r["precursor"])].append(r)

    results = []

    for (outcome, precursor), studies in sorted(groups.items()):
        k = len(studies)
        if k < MIN_K:
            results.append({
                "outcome": outcome,
                "precursor": precursor,
                "k": k,
                "funnel_plot": "Not generated (k < 3)",
                "egger_intercept": "",
                "egger_se": "",
                "egger_t": "",
                "egger_df": "",
                "egger_p": "",
                "interpretation": f"k = {k}: insufficient for publication bias assessment",
            })
            continue

        tes = [s["md"] for s in studies]
        ses = [s["se_md"] for s in studies]
        summary_te, summary_se, tau2 = pairwise_ma(studies)

        # Egger's test
        egger = egger_test(tes, ses)

        # Funnel plot
        fname = f"funnel_{outcome}_{precursor}_vs_PBO"
        draw_funnel(tes, ses, summary_te, outcome, precursor, k, egger, fname)

        # Interpretation
        if egger:
            if egger["p_value"] < 0.10:
                interp = f"Egger's test significant (p = {egger['p_value']:.3f}), but k = {k} is far below recommended minimum (k >= 10); result unreliable"
            else:
                interp = f"Egger's test non-significant (p = {egger['p_value']:.3f}); k = {k} provides insufficient power to detect asymmetry"
        else:
            interp = f"Egger's test not computable (k = {k})"

        results.append({
            "outcome": outcome,
            "precursor": precursor,
            "k": k,
            "funnel_plot": f"{fname}.png",
            "egger_intercept": f"{egger['intercept']:.4f}" if egger else "N/A",
            "egger_se": f"{egger['se_intercept']:.4f}" if egger else "N/A",
            "egger_t": f"{egger['t_stat']:.4f}" if egger else "N/A",
            "egger_df": str(egger["df"]) if egger else "N/A",
            "egger_p": f"{egger['p_value']:.4f}" if egger else "N/A",
            "interpretation": interp,
        })
        print(f"  {outcome:20s} {precursor} vs PBO  k={k}  "
              f"Egger p={egger['p_value']:.3f}" if egger else
              f"  {outcome:20s} {precursor} vs PBO  k={k}  Egger: N/A")

    # Write CSV
    outpath = os.path.join(RES, "publication_bias_assessment.csv")
    fields = ["outcome", "precursor", "k", "funnel_plot",
              "egger_intercept", "egger_se", "egger_t", "egger_df",
              "egger_p", "interpretation"]
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    # Also write as Supplementary Table S10
    s10_path = os.path.join(RES, "supp_table_S10_publication_bias.csv")
    import shutil
    shutil.copy2(outpath, s10_path)

    n_plots = sum(1 for r in results if r["funnel_plot"].endswith(".png"))
    print(f"\nGenerated {n_plots} funnel plots")
    print(f"Saved: {outpath}")
    print(f"Saved: {s10_path}")
    print(f"\nNOTE: All comparisons have k <= 4. Funnel plots and Egger's")
    print(f"test are provided for transparency but are not reliably")
    print(f"interpretable below k = 10 (Cochrane Handbook, §13.3.5.4).")


if __name__ == "__main__":
    main()
