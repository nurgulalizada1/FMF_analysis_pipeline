#!/usr/bin/env python3
"""
07_external_figures.py -- Figure 6: cohort allele frequencies against the
population reference (stage 06 output).

Panel A  enrichment of each variant relative to the gnomAD ancestry group with
         the highest frequency, log scale, ordered by enrichment
Panel B  heatmap of allele frequency across all ten gnomAD ancestry groups
         plus the cohort, log scale
"""

# Make the package importable when the script is run without installing it.
import pathlib as _pathlib
import sys as _sys
_sys.path.insert(0, str(_pathlib.Path(__file__).resolve().parents[1] / "src"))

from mefv import config

config.pin_build_timestamp()

import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

FIG = str(config.FIGURES)
TAB = str(config.TABLES)
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "figure.dpi": 150,
    "savefig.bbox": "tight", "pdf.fonttype": 42,
})
CM = 1 / 2.54
NAVY, TEAL, AMBER, GREY = "#1F3864", "#1D9E75", "#BA7517", "#888780"

t = pd.read_csv(f"{TAB}/table6_gnomad_comparison.csv", sep=";", encoding="utf-8-sig")
t = t[t.cohort_af > 0].copy()

# Classification for colour: the pathogenic and likely pathogenic tier of the
# variant dictionary, the same classification every other stage uses. (Earlier
# versions read a separate list from the study parameters, which disagreed
# with the dictionary; see docs/development_status.md.)
_dic = pd.read_csv(config.DICTIONARY, sep=";", dtype=str, keep_default_na=False,
                   encoding="utf-8-sig")
PLP = set(_dic.loc[_dic.tier.isin(["P", "LP"]), "variant"])
BENIGN_POP = set(t.loc[t.BA1 == "yes", "variant"])

def tier(v):
    if v in PLP:
        return "Pathogenic / Likely pathogenic"
    if v in BENIGN_POP:
        return "BA1 met (population frequency)"
    return "Uncertain significance"

COLOUR = {"Pathogenic / Likely pathogenic": NAVY,
          "BA1 met (population frequency)": AMBER,
          "Uncertain significance": GREY}
t["tier"] = t.variant.map(tier)

# ---------------------------------------------------------------- Panel A ---
a = t[np.isfinite(t.x_vs_highest)].sort_values("x_vs_highest")
fig, ax = plt.subplots(figsize=(9 * CM, 10 * CM))
y = np.arange(len(a))
ax.barh(y, a.x_vs_highest, height=0.66, zorder=2,
        color=[COLOUR[x] for x in a.tier], alpha=0.9)
ax.axvline(1, color="0.35", lw=0.8, ls="--", zorder=3)
ax.set_xscale("log")
ax.set_yticks(y)
ax.set_yticklabels(a.variant)
ax.set_xlabel("Cohort allele frequency ÷ highest gnomAD ancestry group")
ax.grid(axis="x", lw=0.4, color="0.9", zorder=0)
ax.set_axisbelow(True)
ax.set_xlim(0.02, 200)
for yi, (v, g) in enumerate(zip(a.x_vs_highest, a.highest_group)):
    ax.text(v * 1.15, yi, g, va="center", fontsize=5.5, color="0.35")
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in COLOUR.values()]
ax.legend(handles, COLOUR.keys(), frameon=False, loc="lower right", fontsize=6)
for ext in ("pdf", "png"):
    fig.savefig(f"{FIG}/figure6a_enrichment.{ext}", dpi=300)
plt.close(fig)

# ---------------------------------------------------------------- Panel B ---
GRP = ["gnomad_mid", "gnomad_nfe", "gnomad_fin", "gnomad_asj", "gnomad_sas",
       "gnomad_eas", "gnomad_afr", "gnomad_amr", "gnomad_rmi", "gnomad_ami"]
NAMES = ["Middle Eastern", "European (non-Finnish)", "European (Finnish)",
         "Ashkenazi Jewish", "South Asian", "East Asian", "African",
         "Admixed American", "Remaining", "Amish"]
b = t.sort_values("cohort_af", ascending=False)
M = np.vstack([b.cohort_af.values, *[b[g].values for g in GRP]]).T
M = np.where(M <= 0, np.nan, M)

fig, ax = plt.subplots(figsize=(13.5 * CM, 11.5 * CM))
im = ax.imshow(M, cmap="YlOrBr", norm=LogNorm(vmin=1e-5, vmax=0.4), aspect="auto")
ax.set_xticks(range(len(GRP) + 1))
_label = config.load_study_parameters().get("cohort_label", "COHORT")
ax.set_xticklabels([_label] + [n.replace("\n", " ") for n in NAMES],
                   rotation=55, ha="right", fontsize=6.2)
ax.set_yticks(range(len(b)))
ax.set_yticklabels(b.variant, fontsize=6.5)
ax.axvline(0.5, color="0.15", lw=1.4)
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(length=0)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if np.isnan(M[i, j]):
            ax.text(j, i, "·", ha="center", va="center", fontsize=6, color="0.6")
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
cb.set_label("Allele frequency", fontsize=7)
cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=6, length=2)
ax.set_title("Allele frequency: cohort vs gnomAD v4.1.1 ancestry groups",
             fontsize=8, pad=8)
for ext in ("pdf", "png"):
    fig.savefig(f"{FIG}/figure6b_ancestry_heatmap.{ext}", dpi=300)
plt.close(fig)
print("  figure6a_enrichment\n  figure6b_ancestry_heatmap")
