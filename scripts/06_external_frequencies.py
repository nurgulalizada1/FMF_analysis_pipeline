#!/usr/bin/env python3
"""
06_external_frequencies.py -- cohort allele frequencies against a population
reference, across all genetic ancestry groups rather than a single comparator.

The reference table, its database, release and retrieval date are read from
the file named by MEFV_POPULATION_FREQUENCIES and validated by
mefv.reference; its provenance is printed with the results.

  A. enrichment of each variant relative to every ancestry group
  B. population affinity: how closely the cohort allele-frequency profile
     matches each ancestry group, computed on non-enriched variants only.
     Variants enriched by clinical ascertainment are excluded, otherwise
     they would swamp any genuine ancestry signal.
  C. sex-specific frequencies as a negative control (MEFV is autosomal)
"""

# Make the package importable when the script is run without installing it.
import pathlib as _pathlib
import sys as _sys
_sys.path.insert(0, str(_pathlib.Path(__file__).resolve().parents[1] / "src"))

from mefv import config
from mefv.reference import ReferenceDataError, load_population_frequencies
from mefv.utils import _leading_comment_lines as _banner

import os
import numpy as np
import pandas as pd

VAR = str(config.VARIANTS)
PAT = str(config.PATIENTS)
OUT = str(config.TABLES)
os.makedirs(OUT, exist_ok=True)

GROUPS = ["mid", "nfe", "fin", "asj", "sas", "eas", "afr", "amr", "rmi", "ami"]
LABEL = {"mid": "Middle Eastern", "nfe": "European (non-Finnish)",
         "fin": "European (Finnish)", "asj": "Ashkenazi Jewish",
         "sas": "South Asian", "eas": "East Asian",
         "afr": "African/African American", "amr": "Admixed American",
         "rmi": "Remaining", "ami": "Amish"}

# Population reference frequencies, transcribed from the external database.
# The table lives outside the repository because the set of variants it lists
# is the study's observed spectrum. mefv.reference validates its columns,
# values and provenance header; any problem stops the stage.
_PF = config.POPULATION_FREQUENCIES
_PARAMS = config.load_study_parameters()
_EXPECTED = (_PARAMS.get("reference_data") or {}).get("population_frequencies_version")
try:
    _pf, _prov = load_population_frequencies(_PF, expected_version=_EXPECTED)
except ReferenceDataError as exc:
    raise SystemExit(f"stage 06: {exc}") from None
G = {r.variant: {k: v for k, v in r._asdict().items() if k not in ("Index", "variant")}
     for r in _pf.itertuples()}
print(f"population reference  <- {_PF} ({len(G)} variants)")
for _k in ("source", "version", "retrieved"):
    print(f"reference {_k:9s}: {_prov[_k]}")

# Variants used for the population-affinity comparison. Cohort-specific, so the
# list is a study parameter rather than a literal here.
NON_ENRICHED = list(_PARAMS["affinity_comparison_variants"])


def read(p):
    d = pd.read_csv(p, sep=";", encoding="utf-8-sig", dtype=str, keep_default_na=False, skiprows=_banner(p))
    return d.loc[:, [c for c in d.columns if not c.startswith("Unnamed")]]


var, pat = read(VAR), read(PAT)
N = len(pat)


def cohort_af(v):
    s = var[var.variant == v]
    het = int((s.zygosity == "Het").sum())
    hom = int((s.zygosity == "Hom").sum())
    return (2 * hom + het) / (2 * N), het + hom


rows = []
for v, g in G.items():
    af, carriers = cohort_af(v)
    r = dict(variant=v, carriers=carriers, cohort_af=round(af, 5))
    for k in GROUPS:
        r[f"gnomad_{k}"] = g[k]
    r["gnomad_global"] = g["glob"]
    best = max(GROUPS, key=lambda k: g[k])
    r["highest_group"] = LABEL[best]
    r["highest_af"] = g[best]
    r["x_vs_highest"] = round(af / g[best], 2) if g[best] > 0 else np.inf
    r["x_vs_mid"] = round(af / g["mid"], 2) if g["mid"] > 0 else np.inf
    r["x_vs_global"] = round(af / g["glob"], 1) if g["glob"] > 0 else np.inf
    r.update(grpmax_filtering_af=g["grpmax"], gnomad_n_hom=g["hom"],
             gnomad_af_XX=g["xx"], gnomad_af_XY=g["xy"],
             BA1="yes" if g["grpmax"] > 0.05 else "no",
             BS1="yes" if g["grpmax"] > 0.01 else "no",
             revel=g["revel"], cadd=g["cadd"],
             clinvar=g["cv"], clinvar_stars=g["stars"])
    rows.append(r)

t = pd.DataFrame(rows).sort_values("x_vs_highest", ascending=False)
t.to_csv(f"{OUT}/table6_gnomad_comparison.csv", sep=";", index=False, encoding="utf-8-sig")

sub = [v for v in NON_ENRICHED if v in G]
coh = np.array([cohort_af(v)[0] for v in sub])
aff = []
for k in GROUPS:
    grp = np.array([G[v][k] for v in sub])
    rs = pd.Series(coh).corr(pd.Series(grp), method="spearman")
    rp = pd.Series(np.log10(coh + 1e-6)).corr(pd.Series(np.log10(grp + 1e-6)))
    ratio = np.abs(np.log10((coh + 1e-6) / (grp + 1e-6))).mean()
    aff.append(dict(group=LABEL[k], spearman_rho=round(rs, 3),
                    pearson_log10=round(rp, 3), mean_abs_log_ratio=round(ratio, 3)))
a = pd.DataFrame(aff).sort_values("mean_abs_log_ratio")
a.to_csv(f"{OUT}/table7_population_affinity.csv", sep=";", index=False, encoding="utf-8-sig")

pd.set_option("display.width", 260)
print(f"Cohort n = {N}\n")
print("=== A. enrichment relative to the HIGHEST gnomAD ancestry group ===")
print(t[["variant", "carriers", "cohort_af", "highest_group", "highest_af",
         "x_vs_highest", "x_vs_mid", "x_vs_global"]].to_string(index=False))

print("\n=== ACMG population criteria ===")
print(t[["variant", "grpmax_filtering_af", "gnomad_n_hom", "BA1", "BS1",
         "clinvar", "clinvar_stars", "revel"]].to_string(index=False))

print("\n=== C. sex-specific gnomAD frequencies (autosomal — expect equality) ===")
s = t[["variant", "gnomad_af_XX", "gnomad_af_XY"]].copy()
s["ratio_XX_XY"] = (s.gnomad_af_XX / s.gnomad_af_XY.replace(0, np.nan)).round(3)
print(f"  median XX/XY ratio: {s.ratio_XX_XY.median():.3f}"
      f"   range {s.ratio_XX_XY.min():.2f}–{s.ratio_XX_XY.max():.2f}")

print(f"\n=== B. population affinity, non-enriched variants only (n={len(sub)}) ===")
print(f"  variants: {', '.join(sub)}")
print(a.to_string(index=False))
print(f"\nWritten to {OUT}/")
