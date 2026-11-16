#!/usr/bin/env python3
"""
08_functional_and_panel.py

A. Functional concordance
   Do variants classified as hypermorphic in the SpeckSeq assay (Bronnec et al.,
   as recorded in Infevers) show greater enrichment in this clinically
   ascertained cohort than variants classified as not hypermorphic?
   These are two independent lines of evidence: a laboratory functional assay
   and allele frequency in referred patients.

B. Diagnostic panel performance
   If a targeted panel tested only the n most frequent variants, what fraction
   of patients carrying at least one variant would it detect? Computed under
   two definitions of a reportable finding.

Outputs -> tables/ and figures/
"""

# Make the package importable when the script is run without installing it.
import pathlib as _pathlib
import sys as _sys
_sys.path.insert(0, str(_pathlib.Path(__file__).resolve().parents[1] / "src"))

from mefv import config
from mefv.utils import _leading_comment_lines as _banner

config.pin_build_timestamp()

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mefv import DOSE
from mefv.config import AMBER, GREY, NAVY, TEAL

SYNONYMOUS = config.variant_groups()["synonymous"]

VAR = str(config.VARIANTS)
PAT = str(config.PATIENTS)
DICT = str(config.DICTIONARY)
TAB = str(config.TABLES)
FIG = str(config.FIGURES)
for p in (TAB, FIG):
    os.makedirs(p, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "figure.dpi": 150,
    "savefig.bbox": "tight", "pdf.fonttype": 42,
})
CM = 1 / 2.54


def read(p):
    d = pd.read_csv(p, sep=";", encoding="utf-8-sig", dtype=str, keep_default_na=False, skiprows=_banner(p))
    return d.loc[:, [c for c in d.columns if not c.startswith("Unnamed")]]


PANEL_EXONS = set(config.study_design().panel_regions)   # tested in every individual
var, pat, dic = read(VAR), read(PAT), read(DICT)
var = var[var.exon.isin(PANEL_EXONS)].copy()
N = len(pat)
gn = pd.read_csv(f"{TAB}/table6_gnomad_comparison.csv", sep=";", encoding="utf-8-sig")

# ------------------------------------------------- A. functional concordance --
SIMPLE = {"Strongly hypermorphic (FMF-like)": "hypermorphic",
          "Strongly hypermorphic (PAAND-like)": "hypermorphic",
          "Moderately hypermorphic (FMF-like)": "hypermorphic",
          "Moderately hypermorphic (PAAND-like)": "hypermorphic",
          "Not hypermorphic": "not hypermorphic",
          "Uncertain": "functionally uncertain"}
speck = dict(zip(dic.variant, dic.speckseq))
a = gn.copy()
a["speckseq"] = a.variant.map(speck).replace("", np.nan)
a["functional_group"] = a.speckseq.map(SIMPLE)
a = a[a.functional_group.notna() & np.isfinite(a.x_vs_highest)]

summary = (a.groupby("functional_group").x_vs_highest
           .agg(n="count", median="median", q1=lambda s: s.quantile(.25),
                q3=lambda s: s.quantile(.75), min="min", max="max").round(3))

h = a.loc[a.functional_group == "hypermorphic", "x_vs_highest"]
nh = a.loc[a.functional_group == "not hypermorphic", "x_vs_highest"]
u, p_mw = stats.mannwhitneyu(h, nh, alternative="greater")

a[["variant", "speckseq", "functional_group", "cohort_af", "x_vs_highest",
   "x_vs_mid", "grpmax_filtering_af", "clinvar", "clinvar_stars", "revel"]].to_csv(
    f"{TAB}/table8_functional_concordance.csv", sep=";", index=False, encoding="utf-8-sig")

# figure
order = ["hypermorphic", "functionally uncertain", "not hypermorphic"]
cols = {"hypermorphic": NAVY, "functionally uncertain": GREY, "not hypermorphic": AMBER}
fig, ax = plt.subplots(figsize=(9 * CM, 7 * CM))
rng = np.random.default_rng(config.SEED_FIGURES)
for i, g in enumerate(order):
    v = a.loc[a.functional_group == g, "x_vs_highest"].values
    ax.scatter(rng.normal(i, 0.07, len(v)), v, s=26, color=cols[g],
               edgecolor="white", lw=0.5, zorder=3)
    ax.plot([i - 0.26, i + 0.26], [np.median(v)] * 2, color="0.2", lw=1.4, zorder=4)
ax.axhline(1, color="0.45", lw=0.8, ls="--", zorder=1)
ax.set_yscale("log")
ax.set_xticks(range(3))
ax.set_xticklabels([g.replace(" ", "\n") for g in order])
ax.set_ylabel("Cohort AF ÷ highest gnomAD group")
ax.grid(axis="y", lw=0.4, color="0.9", zorder=0)
ax.set_axisbelow(True)
ax.set_title(f"SpeckSeq functional class vs cohort enrichment\n"
             f"hypermorphic vs not hypermorphic: Mann–Whitney p = {p_mw:.3f}",
             fontsize=7.5)
for ext in ("pdf", "png"):
    fig.savefig(f"{FIG}/figure7_functional_concordance.{ext}", dpi=300)
plt.close(fig)

# ------------------------------------------------- B. panel performance -------
TIER = dict(zip(dic.variant, dic.tier))
PLP = {v for v, t in TIER.items() if t in ("P", "LP")}
NONSYN = set(var.variant) - SYNONYMOUS

def cumulative(target_variants, label):
    sub = var[var.variant.isin(target_variants)]
    order_ = (sub.groupby("variant").patient_id.nunique()
              .sort_values(ascending=False).index.tolist())
    covered, rows = set(), []
    denom = sub.patient_id.nunique()
    for i, v in enumerate(order_, 1):
        covered |= set(sub.loc[sub.variant == v, "patient_id"])
        rows.append(dict(definition=label, n_variants_in_panel=i, variant_added=v,
                         patients_detected=len(covered),
                         pct_of_positive_patients=round(len(covered) / denom * 100, 1),
                         pct_of_whole_cohort=round(len(covered) / N * 100, 1)))
    return pd.DataFrame(rows)


p1 = cumulative(PLP, "Pathogenic / Likely pathogenic only")
p2 = cumulative(NONSYN, "All non-synonymous variants")
panel = pd.concat([p1, p2], ignore_index=True)
panel.to_csv(f"{TAB}/table9_panel_performance.csv", sep=";", index=False, encoding="utf-8-sig")

fig, ax = plt.subplots(figsize=(9.5 * CM, 7 * CM))
for df, col, lab in [(p1, NAVY, "Pathogenic / Likely pathogenic"),
                     (p2, TEAL, "All non-synonymous")]:
    ax.step(df.n_variants_in_panel, df.pct_of_positive_patients, where="post",
            color=col, lw=1.6, label=lab)
    ax.scatter(df.n_variants_in_panel, df.pct_of_positive_patients, s=14,
               color=col, zorder=3)
for thr in (90, 95, 99):
    ax.axhline(thr, color="0.8", lw=0.5, ls=":")
    ax.text(0.4, thr + 0.6, f"{thr}%", fontsize=5.5, color="0.5")
ax.set_xlabel("Number of variants included in the panel")
ax.set_ylabel("Variant-positive patients detected (%)")
ax.set_ylim(0, 103)
ax.set_xlim(0.5, max(len(p1), len(p2)) + 0.5)
ax.legend(frameon=False, loc="lower right")
ax.grid(axis="y", lw=0.4, color="0.92")
ax.set_axisbelow(True)
for ext in ("pdf", "png"):
    fig.savefig(f"{FIG}/figure8_panel_performance.{ext}", dpi=300)
plt.close(fig)

# ------------------------------------------------------------------ print ----
pd.set_option("display.width", 220)
print("=== A. SpeckSeq functional class vs cohort enrichment ===")
print(summary.to_string())
print(f"\n  Mann–Whitney U (hypermorphic > not hypermorphic): U = {u:.1f}, p = {p_mw:.4f}")
print("\n  notable discordances:")
disc = a[((a.functional_group == "not hypermorphic") & (a.x_vs_highest > 5)) |
         ((a.functional_group == "hypermorphic") & (a.x_vs_highest < 1)) |
         ((a.functional_group == "functionally uncertain") & (a.x_vs_highest > 20))]
print(disc[["variant", "speckseq", "x_vs_highest", "clinvar"]].to_string(index=False))

print("\n=== B. panel performance ===")
for lab, df in [("P/LP only", p1), ("all non-synonymous", p2)]:
    tot = df.patients_detected.max()
    print(f"\n  {lab}: {len(df)} variants, {tot} positive patients "
          f"({tot / N * 100:.1f}% of cohort)")
    for k in (1, 2, 3, 5, 7, 10):
        if k <= len(df):
            r = df.iloc[k - 1]
            print(f"    top {k:2d} variants -> {r.patients_detected:3d} patients "
                  f"({r.pct_of_positive_patients:5.1f}% of positives, "
                  f"{r.pct_of_whole_cohort:4.1f}% of cohort)")

print("\nWritten: table8_functional_concordance.csv, table9_panel_performance.csv,"
      " figure7_functional_concordance, figure8_panel_performance")
