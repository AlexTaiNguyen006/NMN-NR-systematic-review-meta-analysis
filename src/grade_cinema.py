#!/usr/bin/env python3
"""
GRADE / CINeMA certainty-of-evidence assessment for NMN vs NR NMA.

Evaluates each indirect NMN vs NR comparison across the CINeMA framework:
  1. Within-study bias (from RoB 2)
  2. Reporting bias
  3. Indirectness
  4. Imprecision
  5. Heterogeneity
  6. Incoherence (not applicable — no direct NMN vs NR evidence)

Outputs:
  results/tables/grade_cinema_assessment.csv
  results/figures/grade_cinema_heatmap.png
"""

import csv, math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from collections import defaultdict
from scipy.stats import norm

BASE = "/Users/tainguyen/Documents/NMN_vs_NR_SR_NMA"
DATA = os.path.join(BASE, "data", "extraction", "nma_input_long.csv")
ROB  = os.path.join(BASE, "data", "extraction", "rob2_assessment.csv")
FIGS = os.path.join(BASE, "results", "figures")
RES  = os.path.join(BASE, "results", "tables")

# Style
FONT_FAMILY = "Arial"
DPI = 300
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": [FONT_FAMILY, "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "figure.dpi": DPI,
    "savefig.dpi": DPI,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.15,
})

OUTCOME_LABELS = {
    "FBG": "Fasting Blood Glucose", "HbA1c": "HbA1c", "HOMA-IR": "HOMA-IR",
    "fasting_insulin": "Fasting Insulin", "TC": "Total Cholesterol",
    "LDL": "LDL Cholesterol", "HDL": "HDL Cholesterol", "TG": "Triglycerides",
    "body_weight": "Body Weight", "BMI": "BMI", "NAD+": "Blood NAD+",
    "SBP": "Systolic BP", "DBP": "Diastolic BP", "ALT": "ALT", "AST": "AST",
}


def pretty(oc):
    return OUTCOME_LABELS.get(oc, oc)


# Data loading and meta-analysis (same as main script)

def read_data(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            r["md"] = float(r["md"])
            r["se_md"] = float(r["se_md"])
            r["n_treatment"] = int(r["n_treatment"])
            r["n_control"] = int(r["n_control"])
            rows.append(r)
    return rows


def _norm_unit(u):
    if u is None:
        return ""
    x = u.strip().lower().replace("μ", "u")
    x = x.replace("µ", "u")
    x = x.replace("umol/l", "um")
    x = x.replace("u mol/l", "um")
    return x


def indirect_comparable(outcome_rows):
    nmn_rows = [r for r in outcome_rows if r["precursor"] == "NMN"]
    nr_rows = [r for r in outcome_rows if r["precursor"] == "NR"]
    if not nmn_rows or not nr_rows:
        return False

    for r in outcome_rows:
        notes = (r.get("notes") or "").upper()
        if "INCOMPATIBLE UNIT" in notes:
            return False

    nmn_units = sorted({_norm_unit(r.get("unit", "")) for r in nmn_rows})
    nr_units = sorted({_norm_unit(r.get("unit", "")) for r in nr_rows})
    return nmn_units == nr_units


def pairwise_ma(studies):
    k = len(studies)
    if k == 0:
        return None
    tes = [s["md"] for s in studies]
    ses = [s["se_md"] for s in studies]
    ws = [1.0 / (se**2) for se in ses]
    W = sum(ws)
    te_f = sum(w * te for w, te in zip(ws, tes)) / W
    se_f = 1.0 / math.sqrt(W)
    Q = sum(w * (te - te_f)**2 for w, te in zip(ws, tes))
    df = k - 1
    I2 = max(0.0, (Q - df) / Q) if Q > 0 and df > 0 else 0.0
    C = W - sum(w**2 for w in ws) / W
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    ws_r = [1.0 / (se**2 + tau2) for se in ses]
    W_r = sum(ws_r)
    te_r = sum(w * te for w, te in zip(ws_r, tes)) / W_r
    se_r = 1.0 / math.sqrt(W_r)
    lo_r = te_r - 1.96 * se_r
    up_r = te_r + 1.96 * se_r
    p_r = 2 * (1 - norm.cdf(abs(te_r / se_r))) if se_r > 0 else 1.0
    return {"k": k, "Q": Q, "df": df, "I2": I2, "tau2": tau2,
            "te_r": te_r, "se_r": se_r, "lo_r": lo_r, "up_r": up_r, "p_r": p_r,
            "studies": studies}


def indirect_comparison(ma_a, ma_b):
    te = ma_a["te_r"] - ma_b["te_r"]
    se = math.sqrt(ma_a["se_r"]**2 + ma_b["se_r"]**2)
    lo = te - 1.96 * se
    up = te + 1.96 * se
    p = 2 * (1 - norm.cdf(abs(te / se))) if se > 0 else 1.0
    return {"te": te, "se": se, "lo": lo, "up": up, "p": p}


# RoB data

def load_rob():
    rob = {}
    with open(ROB) as f:
        for r in csv.DictReader(f):
            rob[r["study_id"]] = r["Overall"].strip()
    return rob


# CINeMA domain assessments

def assess_within_study_bias(studies_nmn, studies_nr, rob_dict):
    """
    Domain 1: Within-study bias (from RoB 2 overall judgments).
    All studies in our NMA are 'Some concerns' (the 2 High-RoB studies
    were excluded a priori). Rating:
      - All Low → Low concern
      - Mostly Low, some Some concerns → Some concerns
      - Any High → Major concerns
    """
    judgments = []
    for arm in [studies_nmn, studies_nr]:
        for s in arm:
            sid = s["study_id"]
            j = rob_dict.get(sid, "Some concerns")
            judgments.append(j)

    n_high = sum(1 for j in judgments if j == "High")
    n_some = sum(1 for j in judgments if j == "Some concerns")
    n_low  = sum(1 for j in judgments if j == "Low")

    if n_high > 0:
        return "Major concerns", f"{n_high} High, {n_some} Some concerns, {n_low} Low"
    elif n_some > n_low:
        return "Some concerns", f"{n_some} Some concerns, {n_low} Low"
    elif n_some > 0:
        return "Low concern", f"{n_low} Low, {n_some} Some concerns"
    else:
        return "No concerns", f"All {n_low} Low"


def assess_reporting_bias(studies_nmn, studies_nr):
    """
    Domain 2: Reporting bias.
    Concerns:
      - Small number of studies per arm (potential for small-study effects)
      - Most studies industry-funded or authored
      - Several outcomes are secondary/exploratory
    → Some concerns for all outcomes (systematic issue)
    """
    k_total = len(studies_nmn) + len(studies_nr)
    if k_total <= 2:
        return "Some concerns", f"k={k_total}; too few studies to assess funnel asymmetry; most industry-funded"
    elif k_total <= 4:
        return "Some concerns", f"k={k_total}; insufficient studies for funnel plot; most industry-funded"
    else:
        return "Some concerns", f"k={k_total}; industry funding prevalent; metabolic outcomes often secondary"


def assess_indirectness(outcome, studies_nmn, studies_nr):
    """
    Domain 3: Indirectness (transitivity assumption).
    Key concern: NMN and NR study populations may differ systematically.
    NMN studies: mostly healthy Japanese adults, younger, lower BMI
    NR studies: mostly Western, overweight/obese, some disease-specific
    → Some to Major concerns depending on population heterogeneity
    """
    # Check population mix
    nmn_ids = {s["study_id"] for s in studies_nmn}
    nr_ids  = {s["study_id"] for s in studies_nr}

    # NMN: Yoshino (US, obese women), Huang (India, healthy overweight),
    #       Katayoshi (Japan, healthy), Morifuji (Japan, healthy elderly)
    # NR:  Dollerup (Denmark, obese men), Conze (Canada, healthy overweight),
    #       Remie (Netherlands, obese), Bandi (India, healthy aging)
    # → Different populations, doses, durations — systematic concern for transitivity

    note = ("Populations differ: NMN trials mostly Asian (3/4), "
            "NR trials mostly Western (3/4); BMI ranges differ; "
            "doses 250mg NMN vs 500-2000mg NR; "
            "indirect comparison assumes transitivity")

    return "Some concerns", note


def assess_imprecision(indirect_result, outcome):
    """
    Domain 4: Imprecision.
    Based on CI width and whether CI includes clinically important effects in both directions.
    """
    te = indirect_result["te"]
    lo = indirect_result["lo"]
    up = indirect_result["up"]
    ci_width = up - lo

    # Check if CI crosses zero
    crosses_zero = (lo < 0 < up)

    # Assess based on CI width relative to the effect scale
    # and whether p < 0.05
    sig = indirect_result["p"] < 0.05

    if sig:
        return "Low concern", f"p={indirect_result['p']:.4f}; CI [{lo:.2f}, {up:.2f}] excludes null"
    elif crosses_zero and ci_width > 20:
        return "Major concerns", f"Very wide CI [{lo:.2f}, {up:.2f}]; width={ci_width:.1f}; crosses null"
    elif crosses_zero:
        return "Some concerns", f"CI [{lo:.2f}, {up:.2f}] crosses null; width={ci_width:.1f}"
    else:
        return "Low concern", f"CI [{lo:.2f}, {up:.2f}] narrow and excludes null"


def assess_heterogeneity(ma_nmn, ma_nr):
    """
    Domain 5: Heterogeneity.
    Based on I² and tau² from each pairwise arm.
    """
    i2_max = max(ma_nmn["I2"], ma_nr["I2"])
    k_nmn = ma_nmn["k"]
    k_nr  = ma_nr["k"]

    details = (f"NMN arm: I²={ma_nmn['I2']*100:.0f}% (k={k_nmn}); "
               f"NR arm: I²={ma_nr['I2']*100:.0f}% (k={k_nr})")

    if i2_max == 0 and k_nmn == 1 and k_nr == 1:
        return "Some concerns", details + "; single study per arm — heterogeneity not estimable"
    elif i2_max < 0.25:
        return "No concerns", details
    elif i2_max < 0.50:
        return "Low concern", details
    elif i2_max < 0.75:
        return "Some concerns", details
    else:
        return "Major concerns", details + "; substantial heterogeneity"


def assess_incoherence():
    """
    Domain 6: Incoherence.
    Not applicable — no direct NMN vs NR evidence exists.
    Cannot assess consistency between direct and indirect evidence.
    """
    return "Not applicable", "No direct NMN vs NR comparisons; incoherence cannot be evaluated"


def overall_certainty(ratings):
    """
    Derive overall GRADE certainty from individual domain ratings.
    Start at HIGH (for RCT evidence), downgrade for each concern.
    """
    level = 4  # HIGH
    downgrades = 0
    for domain, (rating, _) in ratings.items():
        if domain == "Incoherence":
            continue  # N/A
        if rating == "Major concerns":
            downgrades += 2
        elif rating == "Some concerns":
            downgrades += 1
        elif rating == "Low concern":
            downgrades += 0  # borderline, no downgrade

    level = max(1, level - downgrades)
    labels = {4: "High", 3: "Moderate", 2: "Low", 1: "Very low"}
    return labels[level], downgrades


# Main assessment

def run_assessment():
    data = read_data(DATA)
    rob_dict = load_rob()

    by_outcome = defaultdict(list)
    for r in data:
        by_outcome[r["outcome"]].append(r)

    results = []
    for oc in sorted(by_outcome.keys()):
        nmn_arm = [s for s in by_outcome[oc] if s["precursor"] == "NMN"]
        nr_arm  = [s for s in by_outcome[oc] if s["precursor"] == "NR"]
        if not nmn_arm or not nr_arm:
            continue
        if not indirect_comparable(by_outcome[oc]):
            continue

        ma_nmn = pairwise_ma(nmn_arm)
        ma_nr  = pairwise_ma(nr_arm)
        ind    = indirect_comparison(ma_nmn, ma_nr)

        # Assess each domain
        ratings = {}
        ratings["Within-study bias"] = assess_within_study_bias(nmn_arm, nr_arm, rob_dict)
        ratings["Reporting bias"]    = assess_reporting_bias(nmn_arm, nr_arm)
        ratings["Indirectness"]      = assess_indirectness(oc, nmn_arm, nr_arm)
        ratings["Imprecision"]       = assess_imprecision(ind, oc)
        ratings["Heterogeneity"]     = assess_heterogeneity(ma_nmn, ma_nr)
        ratings["Incoherence"]       = assess_incoherence()

        overall, n_downgrades = overall_certainty(ratings)

        row = {
            "outcome": oc,
            "outcome_label": pretty(oc),
            "k_NMN": ma_nmn["k"],
            "k_NR": ma_nr["k"],
            "MD": f"{ind['te']:.2f}",
            "CI_lower": f"{ind['lo']:.2f}",
            "CI_upper": f"{ind['up']:.2f}",
            "p_value": f"{ind['p']:.4f}",
        }
        for domain, (rating, justification) in ratings.items():
            domain_key = domain.replace(" ", "_").replace("-", "_")
            row[f"{domain_key}_rating"] = rating
            row[f"{domain_key}_justification"] = justification

        row["overall_certainty"] = overall
        row["n_downgrades"] = n_downgrades
        results.append(row)

    return results


def write_csv(results):
    if not results:
        return
    path = os.path.join(RES, "grade_cinema_assessment.csv")
    fieldnames = list(results[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(results)
    print(f"  Saved: {path}")


# Heatmap figure

RATING_COLORS = {
    "No concerns":    "#4CAF50",   # green
    "Low concern":    "#8BC34A",   # light green
    "Some concerns":  "#FFC107",   # amber
    "Major concerns": "#F44336",   # red
    "Not applicable": "#E0E0E0",   # grey
}

CERTAINTY_COLORS = {
    "High":     "#4CAF50",
    "Moderate": "#8BC34A",
    "Low":      "#FFC107",
    "Very low": "#F44336",
}

DOMAINS = [
    "Within-study bias", "Reporting bias", "Indirectness",
    "Imprecision", "Heterogeneity", "Incoherence",
]


def make_heatmap(results):
    n = len(results)
    nd = len(DOMAINS) + 1  # +1 for overall certainty

    fig_h = max(5.5, 1.8 + n * 0.46)
    fig_w = max(9.5, 2.5 + nd * 1.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Layout parameters
    col_w = 1.35
    row_h = 0.44
    margin_left = 2.8
    margin_top  = 2.0

    ax.set_xlim(0, margin_left + nd * col_w + 0.5)
    ax.set_ylim(0, margin_top + n * row_h + 0.5)
    ax.axis("off")

    # Column headers
    headers = DOMAINS + ["Overall\nCertainty"]
    for j, header in enumerate(headers):
        x = margin_left + j * col_w + col_w / 2
        y = margin_top + n * row_h + 0.15
        fw = "bold" if j == len(DOMAINS) else "normal"
        ax.text(x, y, header, ha="center", va="bottom", fontsize=7.5,
                fontweight=fw, linespacing=1.15)

    # Rows
    from matplotlib.patches import FancyBboxPatch
    for i, row in enumerate(results):
        y = margin_top + (n - 1 - i) * row_h + row_h / 2

        # Outcome label
        ax.text(margin_left - 0.15, y, pretty(row["outcome"]),
                ha="right", va="center", fontsize=8)

        # Domain cells
        for j, domain in enumerate(DOMAINS):
            domain_key = domain.replace(" ", "_").replace("-", "_")
            rating = row[f"{domain_key}_rating"]
            color = RATING_COLORS.get(rating, "#E0E0E0")

            x = margin_left + j * col_w
            rect = FancyBboxPatch(
                (x + 0.05, y - row_h/2 + 0.03), col_w - 0.1, row_h - 0.06,
                boxstyle="round,pad=0.02",
                facecolor=color, edgecolor="white", linewidth=0.8, zorder=2
            )
            ax.add_patch(rect)

            # Short label
            short = {
                "No concerns": "Low", "Low concern": "Low",
                "Some concerns": "Some", "Major concerns": "Major",
                "Not applicable": "N/A"
            }.get(rating, rating)
            text_color = "white" if rating in ("Major concerns",) else "#333333"
            ax.text(x + col_w/2, y, short, ha="center", va="center",
                    fontsize=7.5, fontweight="bold", color=text_color, zorder=3)

        # Overall certainty cell
        j = len(DOMAINS)
        x = margin_left + j * col_w
        overall = row["overall_certainty"]
        color = CERTAINTY_COLORS.get(overall, "#E0E0E0")
        rect = FancyBboxPatch(
            (x + 0.05, y - row_h/2 + 0.03), col_w - 0.1, row_h - 0.06,
            boxstyle="round,pad=0.02",
            facecolor=color, edgecolor="#333333", linewidth=1.0, zorder=2
        )
        ax.add_patch(rect)
        text_color = "white" if overall in ("Very low",) else "#333333"
        # Combined label with GRADE symbols (ASCII-safe)
        symbols = {"High": "++++", "Moderate": "+++-",
                   "Low": "++--", "Very low": "+---"}
        sym = symbols.get(overall, "")
        combined = f"{overall}\n({sym})" if sym else overall
        ax.text(x + col_w/2, y, combined, ha="center", va="center",
                fontsize=6.5, fontweight="bold", color=text_color, zorder=3,
                linespacing=1.3)

    # Legend
    legend_items = [
        ("No concerns / Low concern", "#4CAF50"),
        ("Some concerns", "#FFC107"),
        ("Major concerns", "#F44336"),
        ("Not applicable", "#E0E0E0"),
    ]
    legend_y = 0.6
    legend_x = margin_left
    for label, color in legend_items:
        rect = FancyBboxPatch(
            (legend_x, legend_y - 0.12), 0.3, 0.24,
            boxstyle="round,pad=0.02",
            facecolor=color, edgecolor="white", linewidth=0.5
        )
        ax.add_patch(rect)
        ax.text(legend_x + 0.45, legend_y, label,
                ha="left", va="center", fontsize=7)
        legend_x += 2.3

    fig.suptitle("GRADE/CINeMA Certainty of Evidence: NMN vs NR (Indirect Comparisons)",
                 fontsize=11, fontweight="bold", y=0.98)

    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(FIGS, f"grade_cinema_heatmap.{ext}"),
                    dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("  GRADE/CINeMA heatmap saved.")


# Summary table for manuscript

def print_summary(results):
    print("\n" + "=" * 80)
    print("GRADE/CINeMA SUMMARY — NMN vs NR Indirect Comparisons")
    print("=" * 80)
    print(f"{'Outcome':<22} {'MD [95% CI]':<32} {'Certainty':<12} {'Downgrades'}")
    print("-" * 80)
    for r in results:
        ci = f"{r['MD']} [{r['CI_lower']}, {r['CI_upper']}]"
        print(f"{pretty(r['outcome']):<22} {ci:<32} {r['overall_certainty']:<12} "
              f"{r['n_downgrades']}")
    print("-" * 80)

    # Count by certainty
    counts = defaultdict(int)
    for r in results:
        counts[r["overall_certainty"]] += 1
    print("\nCertainty distribution:")
    for level in ["High", "Moderate", "Low", "Very low"]:
        if counts[level] > 0:
            print(f"  {level}: {counts[level]}/{len(results)} outcomes")


if __name__ == "__main__":
    print("=" * 60)
    print("  GRADE / CINeMA CERTAINTY ASSESSMENT")
    print("=" * 60)

    results = run_assessment()
    write_csv(results)
    make_heatmap(results)
    print_summary(results)

    print("\n" + "=" * 60)
    print("  Assessment complete.")
    print("=" * 60)
