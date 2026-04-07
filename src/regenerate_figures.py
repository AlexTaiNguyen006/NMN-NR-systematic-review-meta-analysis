#!/usr/bin/env python3
# regenerate_figures.py -- generates all forest plots, network graph,
# RoB summary/traffic-light, and GRADE heatmap from nma_input_long.csv.

import csv, math, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import numpy as np
from collections import defaultdict
from scipy.stats import norm

# Global style configuration

# Paths
BASE = "/Users/tainguyen/Documents/NMN_vs_NR_SR_NMA"
DATA = os.path.join(BASE, "data", "extraction", "nma_input_long.csv")
ROB  = os.path.join(BASE, "data", "extraction", "rob2_assessment.csv")
FIGS = os.path.join(BASE, "results", "figures"); os.makedirs(FIGS, exist_ok=True)
RES  = os.path.join(BASE, "results", "tables"); os.makedirs(RES, exist_ok=True)

# Typography
FONT_FAMILY     = "Arial"
TITLE_SIZE      = 12
AXIS_LABEL_SIZE = 10
TICK_LABEL_SIZE = 9
ANNOT_SIZE      = 8.5   # annotations inside plots (CI text, weights, %)
LEGEND_SIZE     = 9
SMALL_SIZE      = 8     # minor annotations

# Colours
COL_NMN        = "#2196F3"   # blue
COL_NR         = "#4CAF50"   # green
COL_PBO        = "#9E9E9E"   # grey
COL_DIAMOND    = "#D32F2F"   # red (summary)
COL_SIG        = "#D32F2F"   # significant markers
COL_NONSIG     = "#1565C0"   # non-sig markers
COL_WHISKER    = "#333333"
COL_ZERO_LINE  = "#888888"

# RoB 2 Cochrane palette
ROB_LOW        = "#4CAF50"
ROB_SOME       = "#FFC107"
ROB_HIGH       = "#F44336"
ROB_COLORS     = {"Low": ROB_LOW, "Some concerns": ROB_SOME, "High": ROB_HIGH}
ROB_SYMBOLS    = {"Low": "+", "Some concerns": "\u2212", "High": "\u00D7"}

# Figure widths
FIG_W_FOREST   = 11     # forest plot total width (inches)
FIG_W_NMA      = 11     # NMA summary forest
FIG_W_NETWORK  = 6.5
FIG_W_ROB_TL   = 10     # RoB traffic-light
FIG_W_ROB_SUM  = 8.5    # RoB summary bar
DPI            = 300

# Apply global matplotlib rcParams
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   [FONT_FAMILY, "Helvetica", "DejaVu Sans"],
    "font.size":         TICK_LABEL_SIZE,
    "axes.titlesize":    TITLE_SIZE,
    "axes.titleweight":  "bold",
    "axes.labelsize":    AXIS_LABEL_SIZE,
    "xtick.labelsize":   TICK_LABEL_SIZE,
    "ytick.labelsize":   TICK_LABEL_SIZE,
    "legend.fontsize":   LEGEND_SIZE,
    "figure.dpi":        DPI,
    "savefig.dpi":       DPI,
    "savefig.bbox":      "tight",
    "savefig.pad_inches": 0.15,
    "axes.spines.top":   False,
    "axes.spines.right": False,
})


# Data loading and meta-analysis functions (unchanged logic)

def read_data(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            r["md"]    = float(r["md"])
            r["se_md"] = float(r["se_md"])
            r["n_treatment"] = int(r["n_treatment"])
            r["n_control"]   = int(r["n_control"])
            r["treat_mean"]  = float(r["treat_mean"])
            r["treat_sd"]    = float(r["treat_sd"])
            r["ctrl_mean"]   = float(r["ctrl_mean"])
            r["ctrl_sd"]     = float(r["ctrl_sd"])
            rows.append(r)
    return rows


def _norm_unit(u):
    if u is None:
        return ""
    x = u.strip().lower().replace("μ", "u")
    x = x.replace("µ", "u")
    x = x.replace("umol/l", "um")
    x = x.replace("u mol/l", "um")
    x = x.replace("pmol/ml", "pmol/ml")
    return x


def indirect_comparable(outcome_rows):
    """
    Return (is_comparable, reason) for NMN-vs-NR indirect comparison.
    Excludes outcomes explicitly flagged as incompatible in extraction notes,
    or with mismatched units across NMN and NR arms.
    """
    nmn_rows = [r for r in outcome_rows if r["precursor"] == "NMN"]
    nr_rows = [r for r in outcome_rows if r["precursor"] == "NR"]
    if not nmn_rows or not nr_rows:
        return False, "missing one precursor arm"

    for r in outcome_rows:
        notes = (r.get("notes") or "").upper()
        if "INCOMPATIBLE UNIT" in notes:
            return False, "flagged as incompatible unit in extraction notes"

    nmn_units = sorted({_norm_unit(r.get("unit", "")) for r in nmn_rows})
    nr_units = sorted({_norm_unit(r.get("unit", "")) for r in nr_rows})
    if nmn_units != nr_units:
        return False, f"unit mismatch across arms: NMN={nmn_units}, NR={nr_units}"

    return True, ""


def pairwise_ma(studies):
    k = len(studies)
    if k == 0:
        return None
    tes = [s["md"] for s in studies]
    ses = [s["se_md"] for s in studies]
    ws  = [1.0 / (se**2) for se in ses]
    W   = sum(ws)

    te_f = sum(w * te for w, te in zip(ws, tes)) / W
    se_f = 1.0 / math.sqrt(W)
    Q    = sum(w * (te - te_f)**2 for w, te in zip(ws, tes))
    df   = k - 1
    I2   = max(0.0, (Q - df) / Q) if Q > 0 and df > 0 else 0.0
    C    = W - sum(w**2 for w in ws) / W
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0

    ws_r = [1.0 / (se**2 + tau2) for se in ses]
    W_r  = sum(ws_r)
    te_r = sum(w * te for w, te in zip(ws_r, tes)) / W_r
    se_r = 1.0 / math.sqrt(W_r)

    ci = lambda t, s: (t - 1.96 * s, t + 1.96 * s)
    lo_f, up_f = ci(te_f, se_f)
    lo_r, up_r = ci(te_r, se_r)
    p_f = 2 * (1 - norm.cdf(abs(te_f / se_f))) if se_f > 0 else 1.0
    p_r = 2 * (1 - norm.cdf(abs(te_r / se_r))) if se_r > 0 else 1.0

    return {
        "k": k, "Q": Q, "df": df, "I2": I2, "tau2": tau2,
        "te_f": te_f, "se_f": se_f, "lo_f": lo_f, "up_f": up_f, "p_f": p_f,
        "te_r": te_r, "se_r": se_r, "lo_r": lo_r, "up_r": up_r, "p_r": p_r,
        "studies": studies, "tes": tes, "ses": ses, "ws_r": ws_r, "ws": ws,
    }


def indirect_comparison(ma_a, ma_b):
    te = ma_a["te_r"] - ma_b["te_r"]
    se = math.sqrt(ma_a["se_r"]**2 + ma_b["se_r"]**2)
    lo = te - 1.96 * se
    up = te + 1.96 * se
    p  = 2 * (1 - norm.cdf(abs(te / se))) if se > 0 else 1.0
    return {"te": te, "se": se, "lo": lo, "up": up, "p": p}


# PRETTY OUTCOME NAMES

OUTCOME_LABELS = {
    "FBG": "Fasting Blood Glucose",
    "HbA1c": "HbA1c",
    "HOMA-IR": "HOMA-IR",
    "fasting_insulin": "Fasting Insulin",
    "TC": "Total Cholesterol",
    "LDL": "LDL Cholesterol",
    "HDL": "HDL Cholesterol",
    "TG": "Triglycerides",
    "body_weight": "Body Weight",
    "BMI": "BMI",
    "body_fat_pct": "Body Fat %",
    "lean_mass": "Lean Mass",
    "REE": "Resting Energy Expenditure",
    "NAD+": "Blood NAD+",
    "SBP": "Systolic BP",
    "DBP": "Diastolic BP",
    "ALT": "ALT",
    "AST": "AST",
    "CRP": "CRP / hs-CRP",
    "IL-6": "IL-6",
}

OUTCOME_UNITS = {
    "FBG": "mg/dL", "HbA1c": "%", "HOMA-IR": "",
    "fasting_insulin": "\u00B5U/mL", "TC": "mg/dL", "LDL": "mg/dL",
    "HDL": "mg/dL", "TG": "mg/dL", "body_weight": "kg",
    "BMI": "kg/m\u00B2", "body_fat_pct": "%", "lean_mass": "kg",
    "REE": "kcal/d", "NAD+": "ng/mL", "SBP": "mmHg", "DBP": "mmHg",
    "ALT": "U/L", "AST": "U/L", "CRP": "mg/L", "IL-6": "pg/mL",
}


def pretty_outcome(oc):
    return OUTCOME_LABELS.get(oc, oc)


def outcome_unit(oc):
    return OUTCOME_UNITS.get(oc, "")


# FOREST PLOT -- 3-column layout: labels | plot | stats

def _draw_hline(axes_list, y, color="#cccccc", lw=0.5):
    """Draw a horizontal separator across all panels at the same y."""
    for ax in axes_list:
        ax.axhline(y, color=color, linewidth=lw, zorder=0)


def forest_plot(ma, outcome, comparison_label, filename):
    """
    Publication-quality forest plot with separate columns for:
      Left   - study label (single line with n)
      Centre - CI whiskers + square
      Right  - MD [95% CI]  Weight%
    """
    if ma is None or ma["k"] == 0:
        return

    k       = ma["k"]
    studies = ma["studies"]
    tes     = ma["tes"]
    ses     = ma["ses"]
    ws_r    = ma["ws_r"]
    max_w   = max(ws_r)
    total_w = sum(ws_r)

    ROW_H   = 0.50     # height per study row (inches)
    HEADER  = 1.0      # extra space at top for title + header row
    FOOTER  = 1.1      # space at bottom for summary diamond + x-label
    fig_h   = HEADER + k * ROW_H + FOOTER

    fig = plt.figure(figsize=(FIG_W_FOREST, max(3.5, fig_h)))
    gs = gridspec.GridSpec(1, 3, width_ratios=[0.28, 0.40, 0.32],
                           figure=fig, wspace=0.04)

    ax_lab  = fig.add_subplot(gs[0])
    ax_for  = fig.add_subplot(gs[1])
    ax_stat = fig.add_subplot(gs[2])

    y_positions = list(range(k, 0, -1))   # top study = k, bottom = 1
    y_summary   = 0

    YLIM = (-0.8, k + 1.0)
    for ax in (ax_lab, ax_for, ax_stat):
        ax.set_ylim(*YLIM)
    HEADER_Y = k + 0.60   # y-coord for column headers
    HLINE_Y  = k + 0.35   # separator under header
    SUMLINE_Y = 0.55      # separator above summary

    # CI range for x-axis scaling
    all_lo = [te - 1.96 * se for te, se in zip(tes, ses)]
    all_up = [te + 1.96 * se for te, se in zip(tes, ses)]
    all_lo.append(ma["lo_r"]); all_up.append(ma["up_r"])
    x_min = min(all_lo); x_max = max(all_up)
    x_pad = (x_max - x_min) * 0.18
    x_left = x_min - x_pad; x_right = x_max + x_pad

    # Shared decorations
    all_axes = [ax_lab, ax_for, ax_stat]
    _draw_hline(all_axes, HLINE_Y)
    _draw_hline(all_axes, SUMLINE_Y)

    for ax in all_axes:
        ax.set_yticks([])
        ax.tick_params(left=False, bottom=False)
        for spine in ax.spines.values():
            spine.set_visible(False)

    # LEFT: study labels (single line)
    ax_lab.set_xlim(0, 1)
    ax_lab.text(0.02, HEADER_Y, "Study", fontsize=ANNOT_SIZE, fontweight="bold",
                va="center", ha="left")

    for i, (s, y) in enumerate(zip(studies, y_positions)):
        sid = s["study_id"].replace("_", " ")
        n_total = s['n_treatment'] + s['n_control']
        label = f"{sid}  (n={n_total})"
        ax_lab.text(0.02, y, label, fontsize=ANNOT_SIZE, va="center", ha="left")

    ax_lab.text(0.02, y_summary,
                f"Overall  (k={k}, I\u00B2={ma['I2']*100:.0f}%)",
                fontsize=ANNOT_SIZE, va="center", ha="left", fontweight="bold",
                color=COL_DIAMOND)

    # CENTRE: forest plot
    ax_for.set_xlim(x_left, x_right)
    ax_for.spines["bottom"].set_visible(True)
    ax_for.tick_params(bottom=True)
    ax_for.axvline(0, color=COL_ZERO_LINE, linestyle="--", linewidth=0.7, zorder=1)

    unit = outcome_unit(outcome)
    xlabel = f"Mean Difference ({unit})" if unit else "Mean Difference"
    ax_for.set_xlabel(xlabel, fontsize=AXIS_LABEL_SIZE, labelpad=8)

    for i, (s, te, se, y) in enumerate(zip(studies, tes, ses, y_positions)):
        lo = te - 1.96 * se; up = te + 1.96 * se
        sz = max(4, 18 * (ws_r[i] / max_w))
        ax_for.plot([lo, up], [y, y], color=COL_WHISKER, linewidth=1.0, zorder=2)
        ax_for.plot(te, y, "s",
                    color=COL_NMN if "NMN" in comparison_label else COL_NR,
                    markersize=sz, zorder=4,
                    markeredgecolor="white", markeredgewidth=0.4)

    # Summary diamond
    te_r = ma["te_r"]; lo_r = ma["lo_r"]; up_r = ma["up_r"]
    dh = 0.25
    diamond_x = [lo_r, te_r, up_r, te_r, lo_r]
    diamond_y = [y_summary, y_summary + dh, y_summary, y_summary - dh, y_summary]
    ax_for.fill(diamond_x, diamond_y, color=COL_DIAMOND, alpha=0.70, zorder=5)
    ax_for.plot(diamond_x, diamond_y, color=COL_DIAMOND, linewidth=0.8, zorder=5)

    # RIGHT: statistics column
    ax_stat.set_xlim(0, 1)

    ax_stat.text(0.02, HEADER_Y, "MD [95% CI]", fontsize=ANNOT_SIZE,
                 fontweight="bold", va="center", ha="left")
    ax_stat.text(0.62, HEADER_Y, "Weight", fontsize=ANNOT_SIZE,
                 fontweight="bold", va="center", ha="center")
    ax_stat.text(0.84, HEADER_Y, "p", fontsize=ANNOT_SIZE,
                 fontweight="bold", va="center", ha="center")

    for i, (te, se, y) in enumerate(zip(tes, ses, y_positions)):
        lo = te - 1.96 * se; up = te + 1.96 * se
        w_pct = ws_r[i] / total_w * 100
        ci_txt = f"{te:.2f} [{lo:.2f}, {up:.2f}]"
        ax_stat.text(0.02, y, ci_txt,
                     fontsize=ANNOT_SIZE, va="center", ha="left")
        ax_stat.text(0.62, y, f"{w_pct:.1f}%",
                     fontsize=ANNOT_SIZE, va="center", ha="center")

    p_str = f"p = {ma['p_r']:.3f}" if ma["p_r"] >= 0.001 else "p < 0.001"
    ax_stat.text(0.02, y_summary,
                 f"{te_r:.2f} [{lo_r:.2f}, {up_r:.2f}]",
                 fontsize=ANNOT_SIZE, va="center", ha="left",
                 fontweight="bold", color=COL_DIAMOND)
    ax_stat.text(0.84, y_summary, p_str,
                 fontsize=ANNOT_SIZE, va="center", ha="center",
                 fontweight="bold", color=COL_DIAMOND)

    # Title
    pretty = pretty_outcome(outcome)
    fig.suptitle(f"{pretty}: {comparison_label}",
                 fontsize=TITLE_SIZE, fontweight="bold", y=0.98)

    fig.savefig(os.path.join(FIGS, filename), dpi=DPI, bbox_inches="tight")
    fig.savefig(os.path.join(FIGS, filename.replace(".png", ".pdf")), bbox_inches="tight")
    plt.close(fig)


# NMA SUMMARY FOREST  (NMN vs NR indirect, all outcomes)

def nma_summary_forest(nma_results_dict, filename_base):
    ocs = sorted(nma_results_dict.keys())
    n = len(ocs)
    if n == 0:
        return

    ROW_H  = 0.48
    HEADER = 1.0
    FOOTER = 1.3
    fig_h  = HEADER + n * ROW_H + FOOTER

    fig = plt.figure(figsize=(FIG_W_NMA, max(4.5, fig_h)))
    gs = gridspec.GridSpec(1, 3, width_ratios=[0.28, 0.40, 0.32],
                           figure=fig, wspace=0.04)

    ax_lab  = fig.add_subplot(gs[0])
    ax_for  = fig.add_subplot(gs[1])
    ax_stat = fig.add_subplot(gs[2])

    y_positions = list(range(n - 1, -1, -1))   # top = n-1, bottom = 0

    YLIM = (-1.0, n + 0.2)
    HEADER_Y = n - 0.15
    HLINE_Y  = n - 0.40

    # Compute x-limits
    all_lo, all_up = [], []
    for oc in ocs:
        v = nma_results_dict[oc]["NMN vs NR"]
        all_lo.append(v["lo"]); all_up.append(v["up"])
    x_min = min(all_lo); x_max = max(all_up)
    x_pad = (x_max - x_min) * 0.18
    x_left = x_min - x_pad; x_right = x_max + x_pad

    all_axes = [ax_lab, ax_for, ax_stat]
    for ax in all_axes:
        ax.set_ylim(*YLIM)
        ax.set_yticks([])
        ax.tick_params(left=False, bottom=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
    _draw_hline(all_axes, HLINE_Y)

    # LEFT: outcome labels
    ax_lab.set_xlim(0, 1)
    ax_lab.text(0.02, HEADER_Y, "Outcome", fontsize=ANNOT_SIZE, fontweight="bold",
                va="center", ha="left")
    for oc, y in zip(ocs, y_positions):
        ax_lab.text(0.02, y, pretty_outcome(oc), fontsize=ANNOT_SIZE,
                    va="center", ha="left")

    # CENTRE: forest
    ax_for.set_xlim(x_left, x_right)
    ax_for.spines["bottom"].set_visible(True)
    ax_for.tick_params(bottom=True)
    ax_for.axvline(0, color=COL_ZERO_LINE, linestyle="--", linewidth=0.7, zorder=1)
    ax_for.set_xlabel("Mean Difference (MD) - NMN vs NR (indirect)",
                       fontsize=AXIS_LABEL_SIZE, labelpad=8)

    for oc, y in zip(ocs, y_positions):
        v = nma_results_dict[oc]["NMN vs NR"]
        color = COL_SIG if v["p"] < 0.05 else COL_NONSIG
        ax_for.plot([v["lo"], v["up"]], [y, y], color=COL_WHISKER, linewidth=1.2, zorder=2)
        ax_for.plot(v["te"], y, "D", color=color, markersize=8, zorder=4,
                    markeredgecolor="white", markeredgewidth=0.4)

    # Favours label centred below x-axis
    ax_for.text(0.5, -0.14, "\u2190 Favours NMN          Favours NR \u2192",
                transform=ax_for.transAxes, fontsize=SMALL_SIZE, ha="center",
                style="italic", color="#888888")

    # RIGHT: stats
    ax_stat.set_xlim(0, 1)
    ax_stat.text(0.04, HEADER_Y, "MD [95% CI]", fontsize=ANNOT_SIZE,
                 fontweight="bold", va="center", ha="left")
    ax_stat.text(0.92, HEADER_Y, "p-value", fontsize=ANNOT_SIZE,
                 fontweight="bold", va="center", ha="right")

    for oc, y in zip(ocs, y_positions):
        v = nma_results_dict[oc]["NMN vs NR"]
        sig = v["p"] < 0.05
        col = COL_SIG if sig else "#333333"
        ci_txt = f"{v['te']:+7.2f}  [{v['lo']:+7.2f}, {v['up']:+7.2f}]"
        p_txt  = f"{v['p']:.3f}" if v["p"] >= 0.001 else "<0.001"
        ax_stat.text(0.04, y, ci_txt, fontsize=ANNOT_SIZE, va="center",
                     ha="left", color=col)
        ax_stat.text(0.82, y, p_txt, fontsize=ANNOT_SIZE, va="center",
                     ha="right", color=col)
        if sig:
            ax_stat.text(0.84, y, "*", fontsize=ANNOT_SIZE, va="center",
                         ha="left", fontweight="bold", color=col)

    fig.suptitle("NMN vs NR: Indirect Comparisons via Network Meta-Analysis",
                 fontsize=TITLE_SIZE, fontweight="bold", y=0.98)

    fig.savefig(os.path.join(FIGS, filename_base + ".png"), dpi=DPI, bbox_inches="tight")
    fig.savefig(os.path.join(FIGS, filename_base + ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print("  NMA summary forest saved.")


# NETWORK GRAPH

def network_graph(n_nmn, n_nr, filename):
    fig, ax = plt.subplots(figsize=(FIG_W_NETWORK, 5.0))

    # Equilateral-ish triangle: PBO at bottom-centre, NMN top-left, NR top-right
    pos = {"Placebo": (0.0, 0.0), "NMN": (-1.5, 2.0), "NR": (1.5, 2.0)}
    colors = {"Placebo": COL_PBO, "NMN": COL_NMN, "NR": COL_NR}
    node_k = {"Placebo": n_nmn + n_nr, "NMN": n_nmn, "NR": n_nr}

    # Edges
    lw_nmn = max(2.0, n_nmn * 0.8)
    lw_nr  = max(2.0, n_nr * 0.8)

    # NMN ↔ Placebo (solid, direct)
    ax.plot([pos["NMN"][0], pos["Placebo"][0]],
            [pos["NMN"][1], pos["Placebo"][1]],
            "-", color="#555555", linewidth=lw_nmn, alpha=0.60, zorder=1)
    # NR ↔ Placebo (solid, direct)
    ax.plot([pos["NR"][0], pos["Placebo"][0]],
            [pos["NR"][1], pos["Placebo"][1]],
            "-", color="#555555", linewidth=lw_nr, alpha=0.60, zorder=1)

    # Edge labels -- compute real angle from coordinates
    import math as _m
    def _edge_angle(x1, y1, x2, y2):
        return _m.degrees(_m.atan2(y2 - y1, x2 - x1))

    ang_nmn = _edge_angle(pos["Placebo"][0], pos["Placebo"][1],
                          pos["NMN"][0], pos["NMN"][1])
    ang_nr  = _edge_angle(pos["Placebo"][0], pos["Placebo"][1],
                          pos["NR"][0], pos["NR"][1])

    mid_nmn = (0.5*(pos["NMN"][0]+pos["Placebo"][0]) - 0.18,
               0.5*(pos["NMN"][1]+pos["Placebo"][1]))
    mid_nr  = (0.5*(pos["NR"][0]+pos["Placebo"][0]) + 0.18,
               0.5*(pos["NR"][1]+pos["Placebo"][1]))

    ax.text(*mid_nmn, f"k = {n_nmn}", fontsize=ANNOT_SIZE, ha="center",
            va="center", color="#555555", rotation=ang_nmn,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))
    ax.text(*mid_nr, f"k = {n_nr}", fontsize=ANNOT_SIZE, ha="center",
            va="center", color="#555555", rotation=ang_nr,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))

    # NMN ↔ NR (dashed, indirect)
    ax.plot([pos["NMN"][0], pos["NR"][0]],
            [pos["NMN"][1], pos["NR"][1]],
            "--", color="#bbbbbb", linewidth=1.0, alpha=0.6, zorder=1)
    ax.text(0.0, 2.22, "indirect", fontsize=SMALL_SIZE, ha="center",
            color="#999999", style="italic")

    # Nodes
    for node, (x, y) in pos.items():
        sz = 700 + node_k[node] * 80
        ax.scatter(x, y, s=sz, c=colors[node], zorder=5,
                   edgecolors="white", linewidths=2.5)
        # Label BELOW node on a single line
        ax.text(x, y - 0.50, f"{node}  (k = {node_k[node]})",
                ha="center", va="top", fontsize=TICK_LABEL_SIZE,
                fontweight="bold")

    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-0.9, 2.7)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Network Geometry", fontsize=TITLE_SIZE, fontweight="bold", pad=14)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, filename), dpi=DPI, bbox_inches="tight")
    fig.savefig(os.path.join(FIGS, filename.replace(".png", ".pdf")), bbox_inches="tight")
    plt.close(fig)
    print("  Network graph saved.")


# ROB 2 FIGURES

ROB_DOMAIN_KEYS = [
    "D1_randomization", "D2_deviations", "D3_missing_data",
    "D4_measurement", "D5_reporting", "Overall",
]
ROB_DOMAIN_LABELS = [
    "D1: Randomization\nprocess",
    "D2: Deviations from\nintended interventions",
    "D3: Missing\noutcome data",
    "D4: Measurement\nof the outcome",
    "D5: Selection of\nthe reported result",
    "Overall",
]

NMN_STUDY_IDS = {
    "Yoshino_2021", "Huang_2022", "Igarashi_2022",
    "Katayoshi_2023", "Morifuji_2024",
}


def _load_rob():
    rows = []
    with open(ROB) as f:
        for r in csv.DictReader(f):
            rows.append(r)
    nmn = sorted([s for s in rows if s["study_id"] in NMN_STUDY_IDS],
                 key=lambda x: x["study_id"])
    nr  = sorted([s for s in rows if s["study_id"] not in NMN_STUDY_IDS],
                 key=lambda x: x["study_id"])
    return nmn + nr


def rob2_traffic_light(filename_base):
    studies = _load_rob()
    n = len(studies)
    nd = len(ROB_DOMAIN_KEYS)

    # Figure sizing -- no aspect="equal"; use Ellipse to keep circles round
    col_w_in = 1.45   # inches per domain column (wide enough for labels)
    row_h_in = 0.52   # inches per study row
    margin_left   = 2.6
    margin_right  = 0.3
    margin_top    = 1.8   # suptitle + multiline domain labels
    margin_bottom = 0.8   # legend
    plot_w = nd * col_w_in
    plot_h = n  * row_h_in
    fig_w = margin_left + plot_w + margin_right
    fig_h = margin_top  + plot_h + margin_bottom

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    # Position axes explicitly so we control margins
    ax_left   = margin_left / fig_w
    ax_bottom = margin_bottom / fig_h
    ax_w      = plot_w / fig_w
    ax_h      = plot_h / fig_h
    ax.set_position([ax_left, ax_bottom, ax_w, ax_h])

    # Data coordinates: each column spans [0,1) in x per domain, each row [0,1) in y
    ax.set_xlim(0, nd)
    ax.set_ylim(0, n)

    # Study labels
    study_labels = []
    for s in studies:
        sid = s["study_id"].replace("_", " ")
        tag = "NMN" if s["study_id"] in NMN_STUDY_IDS else "NR"
        study_labels.append(f"{sid}  [{tag}]")

    # Compute Ellipse radii so circles LOOK round on screen
    # Data-unit per inch in each direction:
    du_per_in_x = nd / plot_w      # data units per inch, x
    du_per_in_y = n  / plot_h      # data units per inch, y
    visual_r_in = 0.19             # desired visual radius in inches
    rx = visual_r_in * du_per_in_x
    ry = visual_r_in * du_per_in_y

    from matplotlib.patches import Ellipse
    for i, study in enumerate(studies):
        y = n - 1 - i
        for j, dk in enumerate(ROB_DOMAIN_KEYS):
            val = study[dk].strip()
            color = ROB_COLORS.get(val, "#cccccc")
            sym   = ROB_SYMBOLS.get(val, "?")
            ell = Ellipse((j + 0.5, y + 0.5), width=2*rx, height=2*ry,
                          facecolor=color, edgecolor="white", lw=1.2, zorder=2)
            ax.add_patch(ell)
            ax.text(j + 0.5, y + 0.5, sym,
                    ha="center", va="center",
                    fontsize=10, fontweight="bold", color="white", zorder=3)

    # Axes ticks & labels
    ax.set_xticks([j + 0.5 for j in range(nd)])
    ax.set_xticklabels(ROB_DOMAIN_LABELS, fontsize=ANNOT_SIZE, ha="center",
                        linespacing=1.15)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_yticks([i + 0.5 for i in range(n)])
    ax.set_yticklabels([study_labels[n - 1 - i] for i in range(n)],
                        fontsize=ANNOT_SIZE)
    ax.tick_params(left=False, bottom=False, top=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Legend
    legend_patches = [
        mpatches.Patch(color=ROB_LOW, label="Low risk"),
        mpatches.Patch(color=ROB_SOME, label="Some concerns"),
        mpatches.Patch(color=ROB_HIGH, label="High risk"),
    ]
    ax.legend(handles=legend_patches, loc="lower center",
              bbox_to_anchor=(0.5, -0.07), ncol=3, fontsize=LEGEND_SIZE,
              frameon=False, handletextpad=0.4, columnspacing=1.5)

    # Title via suptitle so it sits above domain labels
    fig.suptitle("Risk of Bias 2 - Traffic-Light Plot",
                 fontsize=TITLE_SIZE, fontweight="bold",
                 x=ax_left + ax_w/2, y=1 - 0.02)

    fig.savefig(os.path.join(FIGS, filename_base + ".png"), dpi=DPI, bbox_inches="tight")
    fig.savefig(os.path.join(FIGS, filename_base + ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print("  RoB 2 traffic-light saved.")


def rob2_summary_bar(filename_base):
    studies = _load_rob()
    n = len(studies)
    nd = len(ROB_DOMAIN_KEYS)
    cats   = ["Low", "Some concerns", "High"]
    colors = [ROB_LOW, ROB_SOME, ROB_HIGH]

    fig, ax = plt.subplots(figsize=(FIG_W_ROB_SUM, max(3.5, nd * 0.55 + 1.8)))

    y_pos = np.arange(nd)
    bar_h = 0.55

    left = np.zeros(nd)
    for cat, col in zip(cats, colors):
        pcts = []
        for dk in ROB_DOMAIN_KEYS:
            cnt = sum(1 for s in studies if s[dk].strip() == cat)
            pcts.append(cnt / n * 100)
        pcts = np.array(pcts)
        ax.barh(y_pos, pcts, bar_h, left=left, color=col,
                edgecolor="white", linewidth=0.6)
        for idx, (val, l) in enumerate(zip(pcts, left)):
            if val >= 8:
                ax.text(l + val / 2, idx, f"{val:.0f}%", ha="center", va="center",
                        fontsize=SMALL_SIZE, color="white", fontweight="bold")
        left += pcts

    ax.set_yticks(y_pos)
    ax.set_yticklabels(ROB_DOMAIN_LABELS, fontsize=ANNOT_SIZE)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Proportion of studies (%)", fontsize=AXIS_LABEL_SIZE)
    ax.invert_yaxis()
    ax.tick_params(left=False)

    legend_patches = [mpatches.Patch(facecolor=c, label=l)
                      for c, l in zip(colors, cats)]
    ax.legend(handles=legend_patches, loc="lower center",
              bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=LEGEND_SIZE,
              frameon=False, handletextpad=0.4, columnspacing=1.5)

    ax.set_title("Risk of Bias 2 - Summary",
                 fontsize=TITLE_SIZE, fontweight="bold", pad=10)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, filename_base + ".png"), dpi=DPI, bbox_inches="tight")
    fig.savefig(os.path.join(FIGS, filename_base + ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print("  RoB 2 summary bar saved.")


# Main execution

if __name__ == "__main__":
    print("  Regenerating all figures with unified style")

    data = read_data(DATA)
    by_outcome = defaultdict(list)
    for r in data:
        by_outcome[r["outcome"]].append(r)

    # Pairwise meta-analyses + forest plots
    print("\n[1/5] Pairwise forest plots ...")
    all_pairwise = {}
    count = 0
    for oc in sorted(by_outcome.keys()):
        all_pairwise[oc] = {}
        for prec in ["NMN", "NR"]:
            arm = [s for s in by_outcome[oc] if s["precursor"] == prec]
            if len(arm) < 1:
                continue
            ma = pairwise_ma(arm)
            if ma is None:
                continue
            all_pairwise[oc][prec] = ma
            if ma["k"] >= 2:
                comp_label = f"{prec} vs Placebo"
                fname = f"forest_{oc}_{prec}_vs_PBO.png"
                forest_plot(ma, oc, comp_label, fname)
                count += 1
                print(f"       {fname}")
    print(f"       -> {count} forest plots")

    # NMA indirect comparisons
    print("\n[2/5] NMA indirect comparisons ...")
    nma_results_dict = {}
    excluded_indirect = []
    for oc in sorted(all_pairwise.keys()):
        pw = all_pairwise[oc]
        if "NMN" not in pw or "NR" not in pw:
            continue
        ok, reason = indirect_comparable(by_outcome[oc])
        if not ok:
            excluded_indirect.append({"outcome": oc, "reason": reason})
            continue
        ma_nmn = pw["NMN"]; ma_nr = pw["NR"]
        indirect = indirect_comparison(ma_nmn, ma_nr)
        nma_results_dict[oc] = {
            "NMN vs Placebo": {"te": ma_nmn["te_r"], "lo": ma_nmn["lo_r"],
                               "up": ma_nmn["up_r"], "p": ma_nmn["p_r"],
                               "k": ma_nmn["k"], "I2": ma_nmn["I2"],
                               "tau2": ma_nmn["tau2"]},
            "NR vs Placebo":  {"te": ma_nr["te_r"], "lo": ma_nr["lo_r"],
                               "up": ma_nr["up_r"], "p": ma_nr["p_r"],
                               "k": ma_nr["k"], "I2": ma_nr["I2"],
                               "tau2": ma_nr["tau2"]},
            "NMN vs NR":      {"te": indirect["te"], "lo": indirect["lo"],
                               "up": indirect["up"], "p": indirect["p"]},
        }
    print(f"       -> {len(nma_results_dict)} outcomes with indirect comparison")
    if excluded_indirect:
        ex_path = os.path.join(RES, "indirect_exclusion_reasons.csv")
        with open(ex_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["outcome", "reason"])
            w.writeheader(); w.writerows(excluded_indirect)
        print(f"       -> {len(excluded_indirect)} outcomes excluded from indirect comparison")

    # Write result CSVs
    pw_rows = []
    for oc in sorted(all_pairwise.keys()):
        for prec in ["NMN", "NR"]:
            if prec not in all_pairwise[oc]:
                continue
            ma = all_pairwise[oc][prec]
            k1_flag = " [k=1: single study, no meta-analysis]" if ma["k"] == 1 else ""
            pw_rows.append({
                "outcome": oc, "precursor": prec,
                "comparison": f"{prec} vs Placebo",
                "k": ma["k"],
                "MD_random": round(ma["te_r"], 4),
                "lower_CI": round(ma["lo_r"], 4),
                "upper_CI": round(ma["up_r"], 4),
                "p_random": round(ma["p_r"], 4),
                "I2_pct": round(ma["I2"] * 100, 1),
                "tau2": round(ma["tau2"], 4),
                "MD_fixed": round(ma["te_f"], 4),
                "lower_fixed": round(ma["lo_f"], 4),
                "upper_fixed": round(ma["up_f"], 4),
                "p_fixed": round(ma["p_f"], 4),
                "note": k1_flag.strip(),
            })
    pw_path = os.path.join(RES, "pairwise_meta_analysis.csv")
    with open(pw_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pw_rows[0].keys()))
        w.writeheader(); w.writerows(pw_rows)

    nma_rows = []
    for oc in sorted(nma_results_dict.keys()):
        d = nma_results_dict[oc]
        for comp_key, comp_label in [("NMN vs Placebo", "NMN vs Placebo"),
                                      ("NR vs Placebo", "NR vs Placebo"),
                                      ("NMN vs NR", "NMN vs NR")]:
            v = d[comp_key]
            is_indirect = comp_key == "NMN vs NR"
            k_str = f"{d['NMN vs Placebo']['k']}+{d['NR vs Placebo']['k']}" if is_indirect else str(v.get("k", ""))
            k_nmn = d["NMN vs Placebo"]["k"]
            k_nr  = d["NR vs Placebo"]["k"]
            note = ""
            if is_indirect and (k_nmn == 1 or k_nr == 1):
                parts = []
                if k_nmn == 1: parts.append("NMN arm k=1")
                if k_nr == 1: parts.append("NR arm k=1")
                note = f"Caution: {', '.join(parts)}"
            elif not is_indirect and v.get("k") == 1:
                note = "Single study - no meta-analysis"
            nma_rows.append({
                "outcome": oc, "comparison": comp_label,
                "MD": round(v["te"], 4), "lower_CI": round(v["lo"], 4),
                "upper_CI": round(v["up"], 4), "p_value": round(v["p"], 4),
                "k": k_str, "type": "indirect" if is_indirect else "direct",
                "note": note,
            })
    nma_path = os.path.join(RES, "nma_results.csv")
    with open(nma_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(nma_rows[0].keys()))
        w.writeheader(); w.writerows(nma_rows)

    # League table
    lt_rows = []
    for oc in sorted(nma_results_dict.keys()):
        d = nma_results_dict[oc]
        for comp_key in ["NMN vs Placebo", "NR vs Placebo", "NMN vs NR"]:
            v = d[comp_key]
            is_indirect = comp_key == "NMN vs NR"
            lt_rows.append({
                "outcome": oc, "comparison": comp_key,
                "MD_random": round(v["te"], 4),
                "lower_95CI": round(v["lo"], 4),
                "upper_95CI": round(v["up"], 4),
                "p_value": round(v["p"], 4),
                "n_studies": f"{d['NMN vs Placebo']['k']}+{d['NR vs Placebo']['k']}" if is_indirect else v.get("k", ""),
                "I2_pct": "" if is_indirect else round(v.get("I2", 0) * 100, 1),
                "tau2": "" if is_indirect else round(v.get("tau2", 0), 4),
            })
    lt_path = os.path.join(RES, "league_table.csv")
    with open(lt_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(lt_rows[0].keys()))
        w.writeheader(); w.writerows(lt_rows)

    n_k1 = sum(1 for r in pw_rows if r["k"] == 1)
    print(f"\n  Result tables saved ({len(pw_rows)} pairwise, {len(nma_rows)} NMA, {n_k1} single-study warnings)")

    # NMA summary forest
    print("\n[3/5] NMA summary forest plot ...")
    nma_summary_forest(nma_results_dict, "nma_forest_NMN_vs_NR_all")

    # Network graph
    print("\n[4/5] Network graph ...")
    study_prec = set()
    for r in data:
        study_prec.add((r["study_id"], r["precursor"]))
    n_nmn = len(set(s for s, p in study_prec if p == "NMN"))
    n_nr  = len(set(s for s, p in study_prec if p == "NR"))
    network_graph(n_nmn, n_nr, "network_graph.png")

    # RoB 2 figures
    print("\n[5/5] RoB 2 figures ...")
    rob2_traffic_light("rob2_traffic_light")
    rob2_summary_bar("rob2_summary")

    # Done
    print("\nAll figures regenerated.")
    print(f"  Output: {FIGS}/")
