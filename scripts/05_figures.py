#!/usr/bin/env python3
"""
05_figures.py -- publication figures.

Figure 1  variant frequencies with 95% Wilson CIs
Figure 2  genotype matrix of the primary variant pair
Figure 3  linkage disequilibrium heatmap (|D'| upper, r2 lower)
Figure 4  variant positions along pyrin (lollipop)
Figure 5  allele burden by era and assay

Domain boundaries in Figure 4 are approximate and must be verified against
UniProt O15553 before submission.
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
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mefv.eligibility import (
    STATUS_EVALUABLE, cohort_frequency, known_regions, panel_eligible, resolve_region,
    variant_regions,
)
from mefv import allele_dose, genotype_vector, ld_stats, load, wilson_ci
from mefv.config import GREY, NAVY

_GROUPS = config.variant_groups()
SYNONYMOUS = _GROUPS["synonymous"]
BENIGN_COMMON = _GROUPS["benign_common"]
_A, _B = config.linkage_plan()["primary_pair"]

PAT = str(config.PATIENTS)
VAR = str(config.VARIANTS)
FIG = str(config.FIGURES)
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "figure.dpi": 150, "savefig.bbox": "tight", "pdf.fonttype": 42,
})
CM = 1 / 2.54


def no_data(ax, message):
    """A figure with nothing to draw says so, instead of failing or drawing an empty frame."""
    ax.text(0.5, 0.5, message, transform=ax.transAxes, ha="center", va="center",
            fontsize=7, color="0.3")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIG}/{name}.{ext}", dpi=300)
    plt.close(fig)
    print(f"  {name}")


_DESIGN = config.study_design()
ERA_EARLY, ERA_LATE = _DESIGN.primary.label, _DESIGN.sensitivity.label
PANEL_EXONS = set(_DESIGN.panel_regions)   # regions tested in every individual
ASSAY_PANEL, ASSAY_EXT = _DESIGN.panel_assay, _DESIGN.extended_assay

pat, var = load(PAT, VAR, _DESIGN)
# Figure 1 uses per-variant eligible denominators (src/mefv/eligibility.py). The other
# figures compare variants across the panel regions, so they use only individuals whose
# testing record lists every panel region.
pat_all, var_all = pat, var
pat = pat[panel_eligible(pat, _DESIGN.panel_regions)].copy()
var = var[var.patient_id.isin(set(pat.patient_id))].copy()
var_panel = var[var.exon.isin(PANEL_EXONS)]
N = len(pat)
ids_early = list(pat[pat.era == ERA_EARLY].patient_id)
N_early = len(ids_early)
gv = lambda v, ids: genotype_vector(var, v, ids)

print("Writing figures to", FIG)

# ---------------------------------------------------------------- Fig 1 ----
rows = []
_obs = variant_regions(var_all)
_vocab = known_regions(pat_all) | PANEL_EXONS
_panel_ns = set(var_all[var_all.exon.isin(PANEL_EXONS)
                        & ~var_all.variant.isin(SYNONYMOUS)].variant)
for v in sorted(_panel_ns):
    region, src = resolve_region(v, _obs, None, _vocab)
    c = cohort_frequency(pat_all, var_all, v, region, src)
    if c["eligibility_status"] != STATUS_EVALUABLE:
        continue                      # undefined frequency: not plotted, never plotted as 0
    lo, hi = wilson_ci(c["alleles"], c["allele_denominator"])
    rows.append((v, c["allele_freq"], lo, hi, c["n_eligible"]))
f1 = pd.DataFrame(rows, columns=["v", "af", "lo", "hi", "n_eligible"]).sort_values("af")
_n_lo, _n_hi = (int(f1.n_eligible.min()), int(f1.n_eligible.max())) if len(f1) else (0, 0)
f1 = f1[f1.af > 0.002]

fig, ax = plt.subplots(figsize=(9 * CM, 7 * CM))
if f1.empty:
    no_data(ax, "No variant has an evaluable allele frequency\nabove the plotting threshold")
else:
    y = np.arange(len(f1))
    ax.barh(y, f1.af, height=0.62, color=NAVY, alpha=0.85, zorder=2)
    ax.errorbar(f1.af, y, xerr=[f1.af - f1.lo, f1.hi - f1.af], fmt="none",
                ecolor="0.25", elinewidth=0.7, capsize=1.8, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(f1.v)
    ax.set_xlabel("Allele frequency (eligible n = "
                  + (str(_n_lo) if _n_lo == _n_hi else f"{_n_lo}–{_n_hi}") + ")")
    ax.set_xlim(0, max(f1.hi) * 1.14)
    ax.grid(axis="x", lw=0.4, color="0.88", zorder=0); ax.set_axisbelow(True)
    for yi, (a, h) in enumerate(zip(f1.af, f1.hi)):
        ax.text(h + max(f1.hi) * 0.013, yi, f"{a:.3f}", va="center", fontsize=6, color="0.3")
save(fig, "figure1_allele_frequencies")

# ---------------------------------------------------------------- Fig 2 ----
gr, gm = gv(_A, ids_early), gv(_B, ids_early)
M = np.zeros((3, 3), int)
for a, b in zip(gr, gm):
    M[2 - a, b] += 1
fig, ax = plt.subplots(figsize=(8 * CM, 7 * CM))
ax.imshow(M, cmap="Blues", vmin=0, vmax=M.max())
for i in range(3):
    for j in range(3):
        ax.text(j, i, str(M[i, j]), ha="center", va="center", fontsize=10,
                color="white" if M[i, j] > M.max() * 0.55 else "0.15",
                fontweight="bold" if M[i, j] == 0 else "normal")
ax.set_xticks(range(3)); ax.set_xticklabels(["Wild type", "Heterozygous", "Homozygous"])
ax.set_yticks(range(3)); ax.set_yticklabels(["Homozygous", "Heterozygous", "Wild type"])
ax.set_xlabel(_B); ax.set_ylabel(_A)
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(length=0)
r = ld_stats(gr, gm)
ax.set_title(f"n = {N_early}    D' = {r['Dprime']:.2f}    r² = {r['r2']:.2f}", pad=8, fontsize=8)
save(fig, "figure2_primary_pair_matrix")

# ---------------------------------------------------------------- Fig 3 ----
cnt = var_panel[var_panel.patient_id.isin(ids_early)].groupby("variant").patient_id.nunique()
common = sorted(cnt[cnt >= 10].index,
                key=lambda v: _DESIGN.region_rank(var[var.variant == v].exon.iloc[0]) * 1000 +
                int("".join(c for c in v if c.isdigit()) or 0))
G = {v: gv(v, ids_early) for v in common}
n = len(common)
Dm, Rm = np.full((n, n), np.nan), np.full((n, n), np.nan)
for i in range(n):
    for j in range(i + 1, n):
        s = ld_stats(G[common[i]], G[common[j]])
        Dm[i, j], Rm[j, i] = abs(s["Dprime"]), s["r2"]
comb = np.where(np.isnan(Dm), Rm, Dm)

fig, ax = plt.subplots(figsize=(11 * CM, 10 * CM))
if n < 2:
    no_data(ax, "Fewer than two variants reach the carrier threshold:\nno pairwise LD to show")
else:
    im = ax.imshow(comb, cmap="YlOrBr", vmin=0, vmax=1)
    for i in range(n):
        for j in range(n):
            if not np.isnan(comb[i, j]):
                ax.text(j, i, f"{comb[i,j]:.2f}", ha="center", va="center", fontsize=5.5,
                        color="white" if comb[i, j] > 0.65 else "0.2")
    ax.set_xticks(range(n)); ax.set_xticklabels(common, rotation=90)
    ax.set_yticks(range(n)); ax.set_yticklabels(common)
    ax.plot([-0.5, n - 0.5], [-0.5, n - 0.5], color="0.35", lw=0.8)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("Upper triangle |D'|   ·   lower triangle r²", fontsize=8, pad=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=6, length=2)
save(fig, "figure3_ld_heatmap")

# ---------------------------------------------------------------- Fig 4 ----
DOMAINS = [("PYD", 1, 92, "#CECBF6"), ("bZIP", 266, 280, "#9FE1CB"),
           ("B-box", 370, 412, "#FAC775"), ("Coiled-coil", 413, 585, "#D3D1C7"),
           ("B30.2 / SPRY", 586, 776, "#F5C4B3")]
PROT = 781
pc = []
for v, g in var_panel.groupby("variant"):
    d = "".join(c for c in v if c.isdigit())
    if d and int(d) <= PROT:
        pc.append((int(d), len(g), v, v in SYNONYMOUS))
pc.sort(key=lambda x: x[0])
mx = max((c for _, c, _, _ in pc), default=1)

fig, ax = plt.subplots(figsize=(17 * CM, 7.5 * CM))
for name, s, e, col in DOMAINS:
    ax.add_patch(Rectangle((s, -0.7), e - s, 1.4, facecolor=col, edgecolor="0.4", lw=0.4, zorder=1))
    ax.text((s + e) / 2, 0, name, ha="center", va="center", fontsize=6, color="0.15")
ax.plot([1, PROT], [0, 0], color="0.45", lw=0.7, zorder=0)
labelled = []
for codon, cnt_, name, syn in pc:
    h = np.log10(cnt_ + 1) * 3.0
    ax.plot([codon, codon], [0.7, 0.7 + h], color="0.55", lw=0.6, zorder=2)
    ax.scatter(codon, 0.7 + h, s=14 + 26 * np.sqrt(cnt_ / mx), zorder=3,
               facecolor=GREY if syn else NAVY, edgecolor="white", lw=0.5,
               marker="s" if syn else "o")
    if cnt_ >= 20 and not any(abs(codon - c0) < 10 for c0, _ in labelled):
        tier = 0
        while any(abs(codon - c0) < 26 and t0 == tier for c0, t0 in labelled):
            tier += 1
        labelled.append((codon, tier))
        ax.annotate(name, xy=(codon, 0.7 + h), xytext=(codon, 10.6 + tier * 2.6),
                    ha="center", va="bottom", fontsize=6, rotation=90, color="0.2",
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="0.7", shrinkA=0, shrinkB=2))
ax.set_xlim(-25, PROT + 25); ax.set_ylim(-1.4, 20); ax.set_yticks([])
ax.set_xlabel("Amino acid position (pyrin, 781 aa)")
ax.spines["left"].set_visible(False)
h1 = ax.scatter([], [], s=26, facecolor=NAVY, edgecolor="white", lw=0.5, label="Non-synonymous")
h2 = ax.scatter([], [], s=26, facecolor=GREY, edgecolor="white", lw=0.5, marker="s",
                label="Benign polymorphism")
ax.legend(handles=[h1, h2], loc="lower center", bbox_to_anchor=(0.5, -0.36), frameon=False, ncol=2)
save(fig, "figure4_protein_lollipop")

# ---------------------------------------------------------------- Fig 5 ----
b = var_panel[(~var_panel.variant.isin(SYNONYMOUS)) & (~var_panel.variant.isin(BENIGN_COMMON))].copy()
b["al"] = allele_dose(b.zygosity)
pat["burden"] = pat.patient_id.map(b.groupby("patient_id").al.sum()).fillna(0).astype(int)
cats, cols = ["0", "1", "2", "≥3"], ["#D3D1C7", "#85B7EB", "#378ADD", "#1F3864"]
groups = [(ERA_EARLY, pat[pat.era == ERA_EARLY]), (ERA_LATE, pat[pat.era == ERA_LATE]),
          (ASSAY_PANEL, pat[pat.assay == ASSAY_PANEL]), (ASSAY_EXT, pat[pat.assay == ASSAY_EXT]),
          ("All", pat)]
fig, ax = plt.subplots(figsize=(9.5 * CM, 6.5 * CM))
ypos = np.arange(len(groups)); left = np.zeros(len(groups))
for k, (cat, col) in enumerate(zip(cats, cols)):
    vals = np.array([(((s.burden >= 3) if cat == "≥3" else (s.burden == int(cat))).sum()
                      / len(s) * 100) if len(s) else np.nan   # an empty group has no percentage
                     for _, s in groups])
    ax.barh(ypos, vals, left=left, height=0.55, color=col, label=f"{cat} alleles",
            edgecolor="white", lw=0.6)
    for yi, (v, l) in enumerate(zip(vals, left)):
        if v > 6:
            ax.text(l + v / 2, yi, f"{v:.0f}%", ha="center", va="center",
                    fontsize=6.5, color="white" if k >= 2 else "0.15")
    left += vals
ax.set_yticks(ypos)
ax.set_yticklabels([f"{g}\n(n={len(s)})" for g, s in groups])
ax.set_xlabel("Patients tested for every panel region (%)"); ax.set_xlim(0, 100)
ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.18))
save(fig, "figure5_allele_burden")
print("done")
