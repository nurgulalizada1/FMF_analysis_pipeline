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
from mefv.comparison import (
    BA1_THRESHOLD, BS1_THRESHOLD, enrichment, highest_group, population_flag,
)
from mefv.eligibility import (
    cohort_frequency, known_regions, primary_period_scope, resolve_region, variant_regions,
)
from mefv.reference import ReferenceDataError, load_population_frequencies
from mefv.utils import _leading_comment_lines as _banner
from mefv.utils import require_released

import os
import numpy as np
import pandas as pd

VAR = str(config.VARIANTS)
PAT = str(config.PATIENTS)
require_released(PAT)   # refuse tables that a blocked or diagnostic-only stage 01 did not release
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

# Cohort frequency of each reference variant over its ELIGIBLE population: the
# individuals whose testing record lists the variant's region
# (src/mefv/eligibility.py). The earlier version divided by the whole cohort, which
# understates the frequency of any variant in a region that only some individuals
# were tested for, and reported 0 for a variant whose region nobody was recorded as
# tested for. A frequency with no eligible individual is NA, never 0.
OBS_REGIONS = variant_regions(var)
VOCAB = known_regions(pat)
_dict_path = config.DICTIONARY
DICT_REGION = {}
if _dict_path.exists():
    _d = read(_dict_path)
    DICT_REGION = dict(zip(_d.variant, _d.exon_infevers))


_DESIGN = config.study_design()
_PERIOD_RESTRICTED = config.variant_groups()["synonymous"]
_PRIMARY = primary_period_scope(pat, _DESIGN)


def cohort_row(v):
    region, source = resolve_region(v, OBS_REGIONS, DICT_REGION, VOCAB)
    # Variants recorded only in the primary period are evaluated within it.
    if v in _PERIOD_RESTRICTED:
        return cohort_frequency(pat, var, v, region, source, _PRIMARY, _DESIGN.primary.label)
    return cohort_frequency(pat, var, v, region, source)


COHORT = {v: cohort_row(v) for v in G}

rows = []
for v, g in G.items():
    c = COHORT[v]
    af = c["allele_freq"]                       # NaN when not evaluable
    r = dict(variant=v, region=c["region"], region_source=c["region_source"],
             scope=c["scope"], n_cohort=c["n_cohort"], n_in_scope=c["n_in_scope"],
             n_tested=c["n_tested"], n_not_tested=c["n_not_tested"],
             n_undetermined=c["n_undetermined"], n_eligible=c["n_eligible"],
             allele_denominator=c["allele_denominator"],
             n_detected_ineligible=c["n_detected_ineligible"],
             eligibility_status=c["eligibility_status"],
             carriers=c["carriers"], cohort_af=round(af, 5))
    for k in GROUPS:
        r[f"gnomad_{k}"] = g[k]
    r["gnomad_global"] = g["glob"]
    best, best_af, complete = highest_group({k: g[k] for k in GROUPS}, GROUPS)
    r["highest_group"] = LABEL[best] if best else ""
    r["highest_af"] = best_af
    r["highest_group_complete"] = complete
    for col, ref, dp in (("x_vs_highest", best_af, 2), ("x_vs_mid", g["mid"], 2),
                         ("x_vs_global", g["glob"], 1)):
        val, status = enrichment(af, ref)
        r[col] = round(val, dp) if not np.isnan(val) else np.nan
        r[f"{col}_status"] = status
    r.update(grpmax_filtering_af=g["grpmax"], gnomad_n_hom=g["hom"],
             gnomad_af_XX=g["xx"], gnomad_af_XY=g["xy"],
             BA1=population_flag(g["grpmax"], BA1_THRESHOLD),
             BS1=population_flag(g["grpmax"], BS1_THRESHOLD),
             revel=g["revel"], cadd=g["cadd"],
             clinvar=g["cv"], clinvar_stars=g["stars"],
             reference_source=_prov["source"], reference_version=_prov["version"],
             reference_retrieved=_prov["retrieved"])
    rows.append(r)

t = pd.DataFrame(rows).sort_values("x_vs_highest", ascending=False, na_position="last")
t.to_csv(f"{OUT}/table6_gnomad_comparison.csv", sep=";", index=False, encoding="utf-8-sig")

# B. population affinity: only variants with a defined cohort frequency and a
# reference value in every group enter, and the number used is reported. The
# 1e-6 offset in the log terms is an existing choice that keeps log10 defined at
# a cohort frequency of exactly 0; it is not validated (docs/statistics.md).
sub = [v for v in NON_ENRICHED if v in G and not np.isnan(COHORT[v]["allele_freq"])]
aff = []
for k in GROUPS:
    use = [v for v in sub if not np.isnan(G[v][k])]
    coh = np.array([COHORT[v]["allele_freq"] for v in use])
    grp = np.array([G[v][k] for v in use])
    if len(use) >= 2:
        rs = pd.Series(coh).corr(pd.Series(grp), method="spearman")
        rp = pd.Series(np.log10(coh + 1e-6)).corr(pd.Series(np.log10(grp + 1e-6)))
        ratio = np.abs(np.log10((coh + 1e-6) / (grp + 1e-6))).mean()
    else:
        rs = rp = ratio = np.nan                  # too few variants to compare: undefined
    aff.append(dict(group=LABEL[k], n_variants_compared=len(use),
                    n_variants_requested=len([v for v in NON_ENRICHED if v in G]),
                    spearman_rho=round(rs, 3), pearson_log10=round(rp, 3),
                    mean_abs_log_ratio=round(ratio, 3)))
a = pd.DataFrame(aff).sort_values("mean_abs_log_ratio", na_position="last")
a.to_csv(f"{OUT}/table7_population_affinity.csv", sep=";", index=False, encoding="utf-8-sig")

pd.set_option("display.width", 260)
print(f"Cohort n = {N} (total); each frequency uses its own eligible population "
      f"(columns n_eligible, eligibility_status)\n")
print("=== A. enrichment relative to the HIGHEST gnomAD ancestry group ===")
print(t[["variant", "n_eligible", "carriers", "cohort_af", "highest_group", "highest_af",
         "x_vs_highest", "x_vs_mid", "x_vs_global"]].to_string(index=False))
_ne = t[t.x_vs_highest_status != "computed"]
print(f"\n  enrichment not computed for {len(_ne)} variant(s): "
      + ", ".join(f"{k}={v}" for k, v in _ne.x_vs_highest_status.value_counts().items()))

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
