#!/usr/bin/env python3
# prisma_flow.py
# Draws the PRISMA 2020 flow diagram. Exclusion boxes are positioned
# dynamically to avoid overlap. Outputs png + pdf.

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import os

BASE = "/Users/tainguyen/Documents/NMN_vs_NR_SR_NMA"
FIGS = os.path.join(BASE, "results", "figures")
os.makedirs(FIGS, exist_ok=True)

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

# Colours
C_ID     = "#E3F2FD"   # identification
C_SCR    = "#E8F5E9"
C_INCL   = "#FFF3E0"
C_EXCL   = "#FFEBEE"
C_BORDER = "#444444"
C_ARROW  = "#555555"

# layout constants
FIG_W    = 13.0
LABEL_X  = 0.55
CX_MAIN  = 5.2                # main column
CX_EXCL  = 10.5               # exclusion column
WM       = 4.4
WE       = 3.6
MIN_GAP  = 0.35
DEFAULT_GAP = 0.50


class Box:
    """Rounded box w/ stacked text lines."""
    def __init__(self, ax, cx, cy, w, h, facecolor, texts, fontsize=8.5,
                 linewidth=0.8):
        self.cx, self.cy = cx, cy
        self.w, self.h = w, h
        self.top    = cy + h / 2
        self.bottom = cy - h / 2
        self.left   = cx - w / 2
        self.right  = cx + w / 2

        patch = FancyBboxPatch(
            (self.left, self.bottom), w, h,
            boxstyle="round,pad=0.06",
            facecolor=facecolor, edgecolor=C_BORDER,
            linewidth=linewidth, zorder=2,
        )
        ax.add_patch(patch)

        n_real = sum(1 for t, _ in texts if t)
        n_spacer = sum(1 for t, _ in texts if not t)
        if n_real == 0:
            return
        total_wt = n_real + n_spacer * 0.4
        usable = h * 0.82
        step = usable / max(total_wt - 1, 1) if total_wt > 1 else 0
        cur_y = cy + usable / 2

        for txt, kwargs in texts:
            if not txt:
                cur_y -= step * 0.4
                continue
            fs = kwargs.pop("fontsize", fontsize)
            ax.text(cx, cur_y, txt, ha="center", va="center",
                    fontsize=fs, zorder=3, **kwargs)
            cur_y -= step


def _arrow(ax, x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=C_ARROW,
                                lw=1.2, shrinkA=0, shrinkB=0), zorder=1)


def _line(ax, x1, y1, x2, y2):
    ax.plot([x1, x2], [y1, y2], color=C_ARROW, lw=1.2, zorder=1)


def _section_label(ax, y, text):
    ax.text(LABEL_X, y, text, ha="center", va="center", fontsize=10.5,
            fontweight="bold", rotation=90, color="#555555", zorder=3)


def _section_line(ax, y):
    ax.plot([1.2, FIG_W - 0.3], [y, y], color="#cccccc", lw=0.5,
            ls="--", zorder=0)


# compute positions top-down, then draw everything

def main():
    # row specs: (id, main_h, excl_h, excl_cy_offset_below_main)
    # None = no exclusion box for that row
    rows = [
        ("db",    1.20, None, None),
        ("dedup", 0.55, 0.70, 0.00),   # L-branch to duplicates
        ("scr",   0.55, 1.50, 0.00),
        ("ret",   0.55, 0.50, 0.00),
        ("ft",    0.70, 2.80, 0.35),
        ("met",   0.55, 1.55, 0.35),
        ("qual",  0.70, 1.20, 0.10),
        ("nma",   1.00, None, None),
    ]

    # Pass 1: compute cy values
    positions = {}        # id → (cy_main, cy_excl | None)
    cursor = 15.0         # y of top edge of first box
    prev_excl_bottom = None

    for i, (rid, mh, eh, eoff) in enumerate(rows):
        if i == 0:
            cy = cursor - mh / 2
            positions[rid] = (cy, None)
            cursor = cy - mh / 2
            continue

        # dedup has special L-branch logic
        if rid == "dedup":
            gap = DEFAULT_GAP
            mid_y = cursor - gap / 2
            cy = mid_y - gap / 2 - mh / 2
            excl_cy = mid_y
            positions[rid] = (cy, excl_cy)
            prev_excl_bottom = excl_cy - eh / 2
            cursor = cy - mh / 2
            continue

        gap = DEFAULT_GAP

        # Enlarge gap if the right-column box would overlap the previous one
        if eh is not None and prev_excl_bottom is not None:
            proposed_cy = cursor - gap - mh / 2
            proposed_excl_top = proposed_cy - eoff + eh / 2
            need_below = prev_excl_bottom - MIN_GAP
            if proposed_excl_top > need_below:
                gap += (proposed_excl_top - need_below)

        cy = cursor - gap - mh / 2
        excl_cy = (cy - eoff) if eh is not None else None
        positions[rid] = (cy, excl_cy)

        if eh is not None:
            prev_excl_bottom = excl_cy - eh / 2

        cursor = cy - mh / 2

    # Figure height from extent of all boxes
    all_tops = []
    all_bottoms = []
    for rid, (cy, ecy) in positions.items():
        mh = next(r[1] for r in rows if r[0] == rid)
        all_tops.append(cy + mh / 2)
        all_bottoms.append(cy - mh / 2)
        if ecy is not None:
            eh = next(r[2] for r in rows if r[0] == rid)
            all_tops.append(ecy + eh / 2)
            all_bottoms.append(ecy - eh / 2)

    highest = max(all_tops)
    lowest  = min(all_bottoms)
    pad_top = 1.6   # room for title
    pad_bot = 1.6   # room for notes
    FIG_H = (highest + pad_top) - (lowest - pad_bot)
    y_shift = -(lowest - pad_bot)

    # Shift everything so lowest point sits at pad_bot
    for rid in positions:
        cy, ecy = positions[rid]
        positions[rid] = (cy + y_shift,
                          (ecy + y_shift) if ecy is not None else None)

    # Pass 2: draw
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(0, FIG_H)
    ax.axis("off")

    # IDENTIFICATION -
    cy_db, _ = positions["db"]
    b_db = Box(ax, CX_MAIN, cy_db, WM, 1.20, C_ID, [
        ("Records identified from databases", {"fontweight": "bold", "fontsize": 9}),
        ("(n = 1,687)", {"fontsize": 8.5}),
        ("", {}),
        ("PubMed 286  |  Embase 375  |  Scopus 398", {"fontsize": 7.5, "color": "#555"}),
        ("Web of Science 319  |  CENTRAL 309", {"fontsize": 7.5, "color": "#555"}),
    ])

    cy_dedup, cy_dup = positions["dedup"]
    b_dedup = Box(ax, CX_MAIN, cy_dedup, WM, 0.55, C_ID, [
        ("Records after deduplication", {"fontsize": 8.5}),
        ("(n = 1,125)", {"fontsize": 8.5}),
    ])
    b_dup = Box(ax, CX_EXCL, cy_dup, WE, 0.70, C_EXCL, [
        ("Duplicate records removed", {"fontweight": "bold", "fontsize": 8.5}),
        ("(n = 562)", {"fontsize": 8.5}),
    ])
    _line(ax, CX_MAIN, b_db.bottom, CX_MAIN, cy_dup)
    _arrow(ax, CX_MAIN, cy_dup, CX_MAIN, b_dedup.top)
    _arrow(ax, CX_MAIN, cy_dup, b_dup.left, cy_dup)

    sec1 = b_dedup.bottom - DEFAULT_GAP * 0.4
    _section_line(ax, sec1)

    # SCREENING -
    cy_scr, cy_excl_ta = positions["scr"]
    b_scr = Box(ax, CX_MAIN, cy_scr, WM, 0.55, C_SCR, [
        ("Records screened (title / abstract)", {"fontweight": "bold", "fontsize": 8.5}),
        ("(n = 1,125)", {"fontsize": 8.5}),
    ])
    b_excl_ta = Box(ax, CX_EXCL, cy_excl_ta, WE, 1.50, C_EXCL, [
        ("Records excluded (n = 874)", {"fontweight": "bold", "fontsize": 8.5}),
        ("", {}),
        ("No NMN/NR intervention (n = 770)", {"fontsize": 7.5}),
        ("Not an RCT (n = 39)", {"fontsize": 7.5}),
        ("Review / protocol (n = 32)", {"fontsize": 7.5}),
        ("Animal / in vitro (n = 27)", {"fontsize": 7.5}),
        ("Pediatric / T1DM (n = 6)", {"fontsize": 7.5}),
    ])
    _arrow(ax, CX_MAIN, b_dedup.bottom, CX_MAIN, b_scr.top)
    _arrow(ax, b_scr.right, b_scr.cy, b_excl_ta.left, b_scr.cy)

    cy_ret, cy_notret = positions["ret"]
    b_ret = Box(ax, CX_MAIN, cy_ret, WM, 0.55, C_SCR, [
        ("Reports sought for retrieval", {"fontweight": "bold", "fontsize": 8.5}),
        ("(n = 251)", {"fontsize": 8.5}),
    ])
    b_notret = Box(ax, CX_EXCL, cy_notret, WE, 0.50, C_EXCL, [
        ("Reports not retrieved (n = 0)", {"fontsize": 8.5}),
    ])
    _arrow(ax, CX_MAIN, b_scr.bottom, CX_MAIN, b_ret.top)
    _arrow(ax, b_ret.right, b_ret.cy, b_notret.left, b_ret.cy)

    sec2 = b_ret.bottom - DEFAULT_GAP * 0.4
    _section_line(ax, sec2)

    # ELIGIBILITY -
    cy_ft, cy_excl_ft = positions["ft"]
    b_ft = Box(ax, CX_MAIN, cy_ft, WM, 0.70, C_SCR, [
        ("Reports assessed for eligibility", {"fontweight": "bold", "fontsize": 8.5}),
        ("(full-text review)", {"fontsize": 8}),
        ("(n = 251)", {"fontsize": 8.5}),
    ])
    b_excl_ft = Box(ax, CX_EXCL, cy_excl_ft, WE, 2.80, C_EXCL, [
        ("Reports excluded (n = 217)", {"fontweight": "bold", "fontsize": 8.5}),
        ("", {}),
        ("Title/abstract screening (n = 133)", {"fontweight": "bold", "fontsize": 7.5}),
        ("  Not RCT / no placebo arm", {"fontsize": 7, "color": "#555"}),
        ("  Non-metabolic primary outcome", {"fontsize": 7, "color": "#555"}),
        ("  Non-human / in vitro", {"fontsize": 7, "color": "#555"}),
        ("", {}),
        ("Re-screening of uncertain (n = 65)", {"fontweight": "bold", "fontsize": 7.5}),
        ("  Protocol / commentary only", {"fontsize": 7, "color": "#555"}),
        ("  Duplicate cohort data", {"fontsize": 7, "color": "#555"}),
        ("", {}),
        ("Manual full-text review (n = 19)", {"fontweight": "bold", "fontsize": 7.5}),
        ("  No numeric endpoint data", {"fontsize": 7, "color": "#555"}),
        ("  Combination supplement", {"fontsize": 7, "color": "#555"}),
        ("  Population outside PICO", {"fontsize": 7, "color": "#555"}),
    ])
    _arrow(ax, CX_MAIN, b_ret.bottom, CX_MAIN, b_ft.top)
    _arrow(ax, b_ft.right, b_ft.cy, b_excl_ft.left, b_ft.cy)

    cy_met, cy_excl_qual = positions["met"]
    b_met = Box(ax, CX_MAIN, cy_met, WM, 0.55, C_SCR, [
        ("Reports meeting inclusion criteria", {"fontweight": "bold", "fontsize": 8.5}),
        ("(n = 34; NMN = 14, NR = 20)", {"fontsize": 8.5}),
    ])
    b_excl_qual = Box(ax, CX_EXCL, cy_excl_qual, WE, 1.55, C_EXCL, [
        ("Reports excluded from", {"fontweight": "bold", "fontsize": 8.5}),
        ("qualitative synthesis (n = 19)", {"fontweight": "bold", "fontsize": 8.5}),
        ("", {}),
        ("Trial registration only (n = 9)", {"fontsize": 7.5}),
        ("No extractable data / image PDF (n = 3)", {"fontsize": 7.5}),
        ("Combination supplement (n = 2)", {"fontsize": 7.5}),
        ("Secondary/duplicate publication (n = 5)", {"fontsize": 7.5}),
    ])
    _arrow(ax, CX_MAIN, b_ft.bottom, CX_MAIN, b_met.top)
    _arrow(ax, b_met.right, b_met.cy, b_excl_qual.left, b_met.cy)

    sec3 = b_met.bottom - DEFAULT_GAP * 0.4
    _section_line(ax, sec3)

    # INCLUDED -
    cy_qual, cy_excl_nma = positions["qual"]
    b_qual = Box(ax, CX_MAIN, cy_qual, WM, 0.70, C_INCL, [
        ("Studies included in", {"fontweight": "bold", "fontsize": 9}),
        ("qualitative synthesis (n = 15)", {"fontweight": "bold", "fontsize": 9}),
        ("NMN = 5  |  NR = 10", {"fontsize": 8.5, "color": "#555"}),
    ])
    b_excl_nma = Box(ax, CX_EXCL, cy_excl_nma, WE, 1.20, C_EXCL, [
        ("Studies excluded from", {"fontweight": "bold", "fontsize": 8.5}),
        ("quantitative NMA (n = 7)", {"fontweight": "bold", "fontsize": 8.5}),
        ("", {}),
        ("No metabolic endpoint data (n = 5)", {"fontsize": 7.5}),
        ("High RoB + insufficient data (n = 2)", {"fontsize": 7.5}),
    ])
    _arrow(ax, CX_MAIN, b_met.bottom, CX_MAIN, b_qual.top)
    _arrow(ax, b_qual.right, b_qual.cy, b_excl_nma.left, b_qual.cy)

    cy_nma, _ = positions["nma"]
    b_nma = Box(ax, CX_MAIN, cy_nma, WM, 1.0, C_INCL, [
        ("Studies included in quantitative", {"fontweight": "bold", "fontsize": 9.5}),
        ("synthesis - NMA (n = 8)", {"fontweight": "bold", "fontsize": 9.5}),
        ("NMN = 4  |  NR = 4", {"fontsize": 9, "color": "#555"}),
        ("14 outcomes  |  73 data points", {"fontsize": 8, "color": "#777"}),
    ])
    _arrow(ax, CX_MAIN, b_qual.bottom, CX_MAIN, b_nma.top)

    # Section labels (centred vertically in each section) -
    _section_label(ax, (b_db.top + sec1) / 2, "Identification")
    _section_label(ax, (sec1 + sec2) / 2, "Screening")
    _section_label(ax, (sec2 + sec3) / 2, "Eligibility")
    _section_label(ax, (sec3 + b_nma.bottom) / 2, "Included")

    # Bottom notes -
    note_y = b_nma.bottom - 0.70
    ax.text(FIG_W / 2, note_y,
            "Search: January 2018 - May 2025  |  "
            "Databases: PubMed, Embase, Scopus, Web of Science, Cochrane CENTRAL",
            ha="center", va="center", fontsize=7.5, color="#888888",
            style="italic")
    ax.text(FIG_W / 2, note_y - 0.40,
            "PRISMA 2020 flow diagram (Page et al., BMJ 2021;372:n71)",
            ha="center", va="center", fontsize=7.5, color="#888888",
            style="italic")

    fig.suptitle("PRISMA 2020 Flow Diagram",
                 fontsize=14, fontweight="bold", y=0.985)

    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(FIGS, f"prisma_2020_flow.{ext}"),
                    dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("PRISMA flow diagram saved.")


if __name__ == "__main__":
    main()
