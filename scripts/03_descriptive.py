#!/usr/bin/env python3
"""
03_descriptive.py -- cohort description, variant frequencies, HWE, allele burden.

Tables written to tables/:
    table1_cohort.csv
    table2_variant_frequencies.csv
    table3_allele_burden.csv
    supp_hwe.csv
"""

# Make the package importable when the script is run without installing it.
import pathlib as _pathlib
import sys as _sys
_sys.path.insert(0, str(_pathlib.Path(__file__).resolve().parents[1] / "src"))

from mefv import config

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mefv import DOSE, fmt_ci, fmt_p, hwe_exact, load, wilson_ci
from mefv.utils import _leading_comment_lines as _banner

_DESIGN = config.study_design()
ERA_EARLY, ERA_LATE = _DESIGN.primary.label, _DESIGN.sensitivity.label
PANEL_EXONS = set(_DESIGN.panel_regions)   # regions tested in every individual
ASSAY_PANEL, ASSAY_EXT = _DESIGN.panel_assay, _DESIGN.extended_assay

_GROUPS = config.variant_groups()
SYNONYMOUS = _GROUPS["synonymous"]
BLOCK_SYN = _GROUPS["block_synonymous"]
BENIGN_COMMON = _GROUPS["benign_common"]
UNVERIFIED = _GROUPS["unverified"]
_LINKED = config.linkage_plan()["primary_pair"][0]

PAT = str(config.PATIENTS)
VAR = str(config.VARIANTS)
OUT = str(config.TABLES)
os.makedirs(OUT, exist_ok=True)

pat, var = load(PAT, VAR, _DESIGN)
N = len(pat)
pat_early = pat[pat.era == ERA_EARLY]
N_early = len(pat_early)
var_panel = var[var.exon.isin(PANEL_EXONS)].copy()
outside = var[~var.exon.isin(PANEL_EXONS)]

# --------------------------------------------------------------- TABLE 1 ----
rows = []


def add(section, item, a, b, t):
    rows.append(dict(section=section, item=item, era_early=a, era_late=b, total=t))


nA, nB = N_early, N - N_early


def pct(k, n):
    return f"{k} ({k / n * 100:.1f})"


add("Cohort", "Patients, n", nA, nB, N)
add("Cohort", "Period",
    f"{pat_early.year.min()}–{pat_early.year.max()}",
    f"{pat[pat.era == ERA_LATE].year.min()}–{pat[pat.era == ERA_LATE].year.max()}",
    f"{pat.year.min()}–{pat.year.max()}")

for sx in ["male", "female"]:
    a = int(((pat.era == ERA_EARLY) & (pat.sex == sx)).sum())
    b = int(((pat.era == ERA_LATE) & (pat.sex == sx)).sum())
    add("Sex", f"{sx}, n (%)", pct(a, nA), pct(b, nB), pct(a + b, N))

for asy in [ASSAY_PANEL, ASSAY_EXT]:
    a = int(((pat.era == ERA_EARLY) & (pat.assay == asy)).sum())
    b = int(((pat.era == ERA_LATE) & (pat.assay == asy)).sum())
    add("Assay", f"{asy}, n (%)", pct(a, nA), pct(b, nB), pct(a + b, N))

_inferred = pat.assay_source == _DESIGN.inferred_source_label
a = int((_inferred & (pat.era == ERA_EARLY)).sum())
b = int((_inferred & (pat.era == ERA_LATE)).sum())
add("Assay", "assay inferred from laboratory statement, n", a, b, a + b)

for lab, m in [("Variant detected", pat.variant_detected == "yes"),
               ("No variant detected", pat.variant_detected == "no")]:
    a = int((m & (pat.era == ERA_EARLY)).sum())
    b = int((m & (pat.era == ERA_LATE)).sum())
    add("Detection", f"{lab}, n (%)", pct(a, nA), pct(b, nB), pct(a + b, N))

syn_pat = set(var[var.variant.isin(SYNONYMOUS)].patient_id)
a = int(pat[(pat.era == ERA_EARLY) & pat.patient_id.isin(syn_pat)].shape[0])
b = int(pat[(pat.era == ERA_LATE) & pat.patient_id.isin(syn_pat)].shape[0])
add("Record completeness", "≥1 benign polymorphism recorded, n (%)",
    pct(a, nA), pct(b, nB), pct(a + b, N))

# benign recording by assay within the late era -- shows it is not assay-driven
late = pat[pat.era == ERA_LATE]
for asy in [ASSAY_PANEL, ASSAY_EXT]:
    sub = late[late.assay == asy]
    if len(sub):
        k = int(sub.patient_id.isin(syn_pat).sum())
        add("Record completeness", f"  of which {asy}, n (%)", "—", pct(k, len(sub)), "")

vpp = pat.groupby("era").n_variants.agg(["median", "min", "max"])
add("Variants recorded", "Per patient, median (range)",
    f"{vpp.loc[ERA_EARLY,'median']:.0f} ({vpp.loc[ERA_EARLY,'min']}–{vpp.loc[ERA_EARLY,'max']})",
    f"{vpp.loc[ERA_LATE,'median']:.0f} ({vpp.loc[ERA_LATE,'min']}–{vpp.loc[ERA_LATE,'max']})",
    f"{pat.n_variants.median():.0f} ({pat.n_variants.min()}–{pat.n_variants.max()})")

t1 = pd.DataFrame(rows)
t1.to_csv(f"{OUT}/table1_cohort.csv", sep=";", index=False, encoding="utf-8-sig")


# --------------------------------------------------------------- TABLE 2 ----
def freq_table(v_sub, denom, label):
    out = []
    for name, g in v_sub.groupby("variant"):
        het = int((g.zygosity == "Het").sum())
        hom = int((g.zygosity == "Hom").sum())
        wt = denom - het - hom
        if wt < 0:
            continue
        alleles = 2 * hom + het
        lo, hi = wilson_ci(alleles, 2 * denom)
        clo, chi = wilson_ci(het + hom, denom)
        out.append(dict(
            cohort=label, variant=name, exon=g.exon.iloc[0],
            n_tested=denom, wildtype=wt, heterozygous=het, homozygous=hom,
            carriers=het + hom, carrier_freq=round((het + hom) / denom, 4),
            carrier_95CI=fmt_ci(clo, chi), alleles=alleles,
            allele_freq=round(alleles / (2 * denom), 4), allele_95CI=fmt_ci(lo, hi),
            hwe_p=fmt_p(hwe_exact(hom, het, wt)),
            variant_type="synonymous" if name in SYNONYMOUS else "non-synonymous",
            verification="unverified" if name in UNVERIFIED else "verified"))
    df = pd.DataFrame(out)
    return df.sort_values("alleles", ascending=False) if len(df) else df


ns = var_panel[~var_panel.variant.isin(SYNONYMOUS)]
t2_ns = freq_table(ns, N, f"non-synonymous, full cohort (n={N})")

syn_early = var_panel[(var_panel.variant.isin(SYNONYMOUS)) & (var_panel.era == ERA_EARLY)]
t2_syn = freq_table(syn_early, N_early, f"benign polymorphisms, {ERA_EARLY} (n={N_early})")

t2 = pd.concat([t2_ns, t2_syn], ignore_index=True)
t2.to_csv(f"{OUT}/table2_variant_frequencies.csv", sep=";", index=False, encoding="utf-8-sig")
t2[["cohort", "variant", "exon", "variant_type", "wildtype", "heterozygous",
    "homozygous", "allele_freq", "hwe_p"]].to_csv(
    f"{OUT}/supp_hwe.csv", sep=";", index=False, encoding="utf-8-sig")


# --------------------------------------------------------------- TABLE 3 ----
# Two burden definitions.
#   burden_broad : every non-synonymous variant except the common benign one -- the count a
#                  laboratory would produce if it reported all coding changes.
#   burden_plp   : only variants currently classified Pathogenic or Likely
#                  pathogenic in Infevers. This is the clinically meaningful
#                  count and is the primary definition.
DICT = pd.read_csv(str(config.DICTIONARY), sep=";", encoding="utf-8-sig",
                   dtype=str, keep_default_na=False,
                   skiprows=_banner(config.DICTIONARY))
TIER = dict(zip(DICT.variant, DICT.tier))
PLP = {v for v, t in TIER.items() if t in ("P", "LP")}

for name, keep in [("burden_broad", lambda v: v not in SYNONYMOUS and v not in BENIGN_COMMON),
                   ("burden_plp", lambda v: v in PLP)]:
    src = var_panel[var_panel.variant.map(keep)].copy()
    src["alleles"] = src.zygosity.map(DOSE).fillna(1)
    pat[name] = pat.patient_id.map(src.groupby("patient_id").alleles.sum()).fillna(0).astype(int)
pat["burden"] = pat.burden_plp

rows3 = []
groups = [(ERA_EARLY, pat[pat.era == ERA_EARLY]), (ERA_LATE, pat[pat.era == ERA_LATE]),
          (ASSAY_PANEL, pat[pat.assay == ASSAY_PANEL]), (ASSAY_EXT, pat[pat.assay == ASSAY_EXT]),
          ("total", pat)]
for definition in ["burden_plp", "burden_broad"]:
    for label, sub in groups:
        n = len(sub)
        b = sub[definition]
        for cat, m in [("0 alleles", b == 0), ("1 allele", b == 1),
                       ("2 alleles", b == 2), ("≥3 alleles", b >= 3),
                       ("Detection rate (≥1 allele)", b >= 1)]:
            k = int(m.sum())
            lo, hi = wilson_ci(k, n)
            rows3.append(dict(definition=definition, group=label, category=cat,
                              n=k, total=n, percent=round(k / n * 100, 1),
                              ci_95=fmt_ci(lo, hi, 3)))
t3 = pd.DataFrame(rows3)
t3.to_csv(f"{OUT}/table3_allele_burden.csv", sep=";", index=False, encoding="utf-8-sig")

# ------------------------------------------------------------- statistics ---
chi2_sex, p_sex, dof_sex, _ = stats.chi2_contingency(
    pd.crosstab(pat.sex, pat.burden.clip(upper=3)))
_, p_era, _, _ = stats.chi2_contingency(
    pd.crosstab(pat.era, pat.burden >= 1))
_, p_assay, _, _ = stats.chi2_contingency(
    pd.crosstab(pat.assay, pat.burden >= 1))

# Added value of the extended assay: variants found outside the panel regions
ext_ids = set(pat[pat.assay == ASSAY_EXT].patient_id)
extra = outside[outside.patient_id.isin(ext_ids)]

# ------------------------------------------------------------------ print ---
print(f"Patients {N}  |  variant rows {len(var)}  |  within panel exons {len(var_panel)}")
print(f"Assay: {ASSAY_PANEL} {int((pat.assay == ASSAY_PANEL).sum())}, "
      f"{ASSAY_EXT} {int((pat.assay == ASSAY_EXT).sum())}")

print("\n=== TABLE 1 — cohort ===")
print(t1.to_string(index=False))

cols = ["variant", "exon", "n_tested", "heterozygous", "homozygous",
        "allele_freq", "allele_95CI", "hwe_p", "verification"]
print(f"\n=== TABLE 2A — non-synonymous, N={N} ===")
print(t2_ns[cols].to_string(index=False))
print(f"\n=== TABLE 2B — benign polymorphisms, {ERA_EARLY}, N={N_early} ===")
print(t2_syn[cols].to_string(index=False))

print("\n=== HWE by stratum ===")
for lab, s in [("benign, outside the block (control)", t2_syn[~t2_syn.variant.isin(BLOCK_SYN)]),
               (f"benign, block (linked to {_LINKED})", t2_syn[t2_syn.variant.isin(BLOCK_SYN)]),
               ("non-synonymous", t2_ns)]:
    ps = pd.to_numeric(s.hwe_p, errors="coerce").dropna()
    if len(ps):
        print(f"  {lab:38s} n={len(ps):2d}  p>0.05 in {int((ps>0.05).sum())}/{len(ps)}"
              f"  median p={ps.median():.3g}")

print("\n=== TABLE 3 — allele burden ===")
print("  burden_plp   : Pathogenic / Likely pathogenic only (Infevers, primary)")
print("  burden_broad : all non-synonymous except the common benign variant (sensitivity)")
print(t3[t3.definition == "burden_plp"].drop(columns="definition").to_string(index=False))
print()
print(t3[(t3.definition == "burden_broad") & (t3.group == "total")]
      .drop(columns="definition").to_string(index=False))

print("\n=== Variant classification tiers (source annotations) ===")
tier_tab = (var_panel.assign(tier=var_panel.variant.map(TIER))
            .groupby("tier").agg(variants=("variant", "nunique"),
                                 observations=("variant", "size")))
print(tier_tab.to_string())

print(f"\nSex × allele burden    : chi2={chi2_sex:.3f}, df={dof_sex}, p={p_sex:.3f}")
print(f"Era × detection rate   : p={p_era:.3f}")
print(f"Assay × detection rate : p={p_assay:.3f}")

# Counts only: the console log is not a place for record-level rows.
print(f"\n=== Variants outside the panel regions ({ASSAY_EXT} only, "
      f"n={len(ext_ids)} individuals) ===")
print(extra.groupby(["variant", "exon"]).size().rename("observations").to_string()
      if len(extra) else "  none")

print(f"\nWritten to {OUT}/")
