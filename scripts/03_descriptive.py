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
from mefv import DataIntegrityError, allele_dose, fmt_ci, fmt_p, hwe_exact, load, wilson_ci
from mefv.eligibility import (
    STATUS_EVALUABLE, cohort_frequency, known_regions, panel_eligible, region_sets,
    resolve_region, variant_regions,
)
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

# Eligibility (src/mefv/eligibility.py, docs/eligibility.md): a frequency is
# computed over the individuals whose testing record lists the variant's region.
# The whole-cohort count N is used only to describe the cohort.
OBS_REGIONS = variant_regions(var)
VOCAB = known_regions(pat) | PANEL_EXONS
PANEL_OK = panel_eligible(pat, _DESIGN.panel_regions)   # tested for every panel region
UNUSABLE = region_sets(pat).isna()                      # testing record missing/malformed

# --------------------------------------------------------------- TABLE 1 ----
rows = []


def add(section, item, a, b, t):
    rows.append(dict(section=section, item=item, era_early=a, era_late=b, total=t))


nA, nB = N_early, N - N_early


def pct(k, n):
    return f"{k} ({k / n * 100:.1f})" if n else f"{k} (n/a)"      # no denominator: no percentage


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
               ("No variant detected", pat.variant_detected == "no"),
               ("Detection status undetermined", pat.variant_detected == "undetermined")]:
    a = int((m & (pat.era == ERA_EARLY)).sum())
    b = int((m & (pat.era == ERA_LATE)).sum())
    add("Detection", f"{lab}, n (%)", pct(a, nA), pct(b, nB), pct(a + b, N))

for lab, m in [("Tested for every panel region (eligible for panel analyses)", PANEL_OK),
               ("Usable record, not every panel region listed", ~PANEL_OK & ~UNUSABLE),
               ("Testing record missing or unusable (undetermined)", UNUSABLE)]:
    a = int((m & (pat.era == ERA_EARLY)).sum())
    b = int((m & (pat.era == ERA_LATE)).sum())
    add("Eligibility", f"{lab}, n (%)", pct(a, nA), pct(b, nB), pct(a + b, N))

_vr = config.VALIDATION_REPORT
if _vr.exists():
    import json as _json
    _v = _json.loads(_vr.read_text())
    add("Data integrity", "Stage 01 validation status", "", "", _v["status"])
    for _d in _v.get("owner_dispositions", []):
        add("Data integrity",
            f"Owner disposition {_d['rule']}: {_d['decision']} ({_d['records']} records, "
            "unresolved in the data; rows counted exactly as recorded)", "", "", "recorded")

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

def _median_range(x):
    """median (min–max) of the declared n_variants; blank when there is no declared value."""
    x = pd.to_numeric(x, errors="coerce").dropna()
    return f"{x.median():.0f} ({int(x.min())}–{int(x.max())})" if len(x) else ""


add("Variants recorded", "Per patient, median (range)",
    _median_range(pat.loc[pat.era == ERA_EARLY, "n_variants"]),
    _median_range(pat.loc[pat.era == ERA_LATE, "n_variants"]),
    _median_range(pat.n_variants))

t1 = pd.DataFrame(rows)
t1.to_csv(f"{OUT}/table1_cohort.csv", sep=";", index=False, encoding="utf-8-sig")


# --------------------------------------------------------------- TABLE 2 ----
FREQ_COLUMNS = [
    "cohort", "variant", "exon", "n_cohort", "n_in_scope", "n_tested", "n_not_tested",
    "n_undetermined", "n_eligible", "allele_denominator", "wildtype", "heterozygous",
    "homozygous", "carriers", "carrier_freq", "carrier_95CI", "alleles", "allele_freq",
    "allele_95CI", "hwe_p", "n_detected_ineligible", "eligibility_status", "scope",
    "region_source", "variant_type", "verification",
]


def freq_table(variants, label, scope=None, scope_label="whole cohort"):
    """One row per variant, over its eligible population (see eligibility.py)."""
    out = []
    for name in sorted(variants):
        region, source = resolve_region(name, OBS_REGIONS, None, VOCAB)
        c = cohort_frequency(pat, var, name, region, source, scope, scope_label)
        ok = c["eligibility_status"] == STATUS_EVALUABLE
        lo, hi = wilson_ci(c["alleles"], c["allele_denominator"]) if ok else (np.nan, np.nan)
        clo, chi = wilson_ci(c["carriers"], c["n_eligible"]) if ok else (np.nan, np.nan)
        out.append(dict(
            cohort=label, variant=name, exon=c["region"],
            n_cohort=c["n_cohort"], n_in_scope=c["n_in_scope"], n_tested=c["n_tested"],
            n_not_tested=c["n_not_tested"], n_undetermined=c["n_undetermined"],
            n_eligible=c["n_eligible"], allele_denominator=c["allele_denominator"],
            wildtype=c["wildtype"], heterozygous=c["heterozygous"], homozygous=c["homozygous"],
            carriers=c["carriers"],
            carrier_freq=round(c["carrier_freq"], 4), carrier_95CI=fmt_ci(clo, chi),
            alleles=c["alleles"],
            allele_freq=round(c["allele_freq"], 4), allele_95CI=fmt_ci(lo, hi),
            hwe_p=fmt_p(hwe_exact(c["homozygous"], c["heterozygous"], int(c["wildtype"])))
                  if ok else "",
            n_detected_ineligible=c["n_detected_ineligible"],
            eligibility_status=c["eligibility_status"], scope=c["scope"],
            region_source=c["region_source"],
            variant_type="synonymous" if name in SYNONYMOUS else "non-synonymous",
            verification="unverified" if name in UNVERIFIED else "verified"))
    df = pd.DataFrame(out, columns=FREQ_COLUMNS)    # fixed columns: also valid with no variant
    return df.sort_values("alleles", ascending=False) if len(df) else df


ns_variants = set(var_panel[~var_panel.variant.isin(SYNONYMOUS)].variant)
t2_ns = freq_table(ns_variants, f"non-synonymous, full cohort (n={N})")

# Benign polymorphisms were recorded completely only in the primary period (the
# configured period restriction, see docs/limitations.md): their denominator is
# the eligible individuals of that period.
syn_variants = set(var_panel[(var_panel.variant.isin(SYNONYMOUS))
                             & (var_panel.era == ERA_EARLY)].variant)
t2_syn = freq_table(syn_variants, f"benign polymorphisms, {ERA_EARLY} (n={N_early})",
                    scope=(pat.era == ERA_EARLY), scope_label=ERA_EARLY)

# Period-restricted variants (the configured metadata is variant_groups.synonymous:
# recorded completely only in the primary period). Report, for each, the metadata field
# that restricts it, the window used, and whether the records respect that window.
# The metadata holds no other per-variant recording dates, so nothing else is restricted.
_scope_rows = []
for _v in sorted(SYNONYMOUS):
    _r = var_panel[var_panel.variant == _v]
    _n_in, _n_out = int((_r.era == ERA_EARLY).sum()), int((_r.era != ERA_EARLY).sum())
    _scope_rows.append(dict(
        variant=_v, restriction_metadata="variant_groups.synonymous", scope=ERA_EARLY,
        rows_in_scope=_n_in, rows_outside_scope=_n_out,
        status=("not_observed_in_panel_regions" if not len(_r) else
                "records_outside_scope" if _n_out else "consistent")))
pd.DataFrame(_scope_rows).to_csv(f"{OUT}/supp_period_scope.csv", sep=";", index=False,
                                 encoding="utf-8-sig")

t2 = pd.concat([t2_ns, t2_syn], ignore_index=True)
t2.to_csv(f"{OUT}/table2_variant_frequencies.csv", sep=";", index=False, encoding="utf-8-sig")
t2[["cohort", "variant", "exon", "variant_type", "n_eligible", "wildtype", "heterozygous",
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
    src = var_panel[var_panel.variant.map(keep).astype(bool)].copy()   # astype: empty map is object
    src["alleles"] = allele_dose(src.zygosity)
    pat[name] = pat.patient_id.map(src.groupby("patient_id").alleles.sum()).fillna(0).astype(int)
pat["burden"] = pat.burden_plp

# Burden counts variants across all panel regions, so an individual enters only if
# the record lists every panel region. Others are excluded and counted.
pat_all = pat
pat = pat[PANEL_OK].copy()
rows3 = []
groups = [(ERA_EARLY, pat[pat.era == ERA_EARLY]), (ERA_LATE, pat[pat.era == ERA_LATE]),
          (ASSAY_PANEL, pat[pat.assay == ASSAY_PANEL]), (ASSAY_EXT, pat[pat.assay == ASSAY_EXT]),
          ("total", pat)]
_in_group = {ERA_EARLY: pat_all.era == ERA_EARLY, ERA_LATE: pat_all.era == ERA_LATE,
             ASSAY_PANEL: pat_all.assay == ASSAY_PANEL, ASSAY_EXT: pat_all.assay == ASSAY_EXT,
             "total": pd.Series(True, index=pat_all.index)}
for definition in ["burden_plp", "burden_broad"]:
    for label, sub in groups:
        n = len(sub)
        n_excluded = int(_in_group[label].sum()) - n
        b = sub[definition]
        for cat, m in [("0 alleles", b == 0), ("1 allele", b == 1),
                       ("2 alleles", b == 2), ("≥3 alleles", b >= 3),
                       ("Detection rate (≥1 allele)", b >= 1)]:
            k = int(m.sum())
            lo, hi = wilson_ci(k, n)
            rows3.append(dict(definition=definition, group=label, category=cat,
                              n=k, total=n, percent=round(k / n * 100, 1) if n else np.nan,
                              ci_95=fmt_ci(lo, hi, 3),
                              n_excluded_not_eligible=n_excluded))
t3 = pd.DataFrame(rows3)
t3.to_csv(f"{OUT}/table3_allele_burden.csv", sep=";", index=False, encoding="utf-8-sig")

# ------------------------------------------------------------- statistics ---
def chi2_or_nan(table):
    """Chi-square test of independence, or (nan, nan, nan) where it is not defined: a table
    with fewer than two rows or columns, or with an empty row or column (no individuals).
    Never a p value of 1 for a table that cannot be tested."""
    t = np.asarray(table)
    if t.ndim != 2 or min(t.shape) < 2 or (t.sum(axis=0) == 0).any() or (t.sum(axis=1) == 0).any():
        return np.nan, np.nan, np.nan
    chi2, p, dof, _ = stats.chi2_contingency(table)
    return chi2, p, dof


chi2_sex, p_sex, dof_sex = chi2_or_nan(pd.crosstab(pat.sex, pat.burden.clip(upper=3)))
_, p_era, _ = chi2_or_nan(pd.crosstab(pat.era, pat.burden >= 1))
_, p_assay, _ = chi2_or_nan(pd.crosstab(pat.assay, pat.burden >= 1))

# Added value of the extended assay: variants found outside the panel regions
ext_ids = set(pat[pat.assay == ASSAY_EXT].patient_id)
extra = outside[outside.patient_id.isin(ext_ids)]

# ------------------------------------------------------------------ print ---
print(f"Patients {N}  |  variant rows {len(var)}  |  within panel exons {len(var_panel)}")
print(f"Assay: {ASSAY_PANEL} {int((pat.assay == ASSAY_PANEL).sum())}, "
      f"{ASSAY_EXT} {int((pat.assay == ASSAY_EXT).sum())}")

print("\n=== TABLE 1 — cohort ===")
print(t1.to_string(index=False))

cols = ["variant", "exon", "n_eligible", "heterozygous", "homozygous",
        "allele_freq", "allele_95CI", "hwe_p", "verification"]
print(f"\n=== TABLE 2A — non-synonymous, N={N} ===")
print(t2_ns[cols].to_string(index=False))
print(f"\n=== TABLE 2B — benign polymorphisms, {ERA_EARLY}, N={N_early} (cohort in period) ===")
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

print(f"\nSex × allele burden    : chi2={chi2_sex:.3f}, df={dof_sex}, p={p_sex:.3f}"
      + ("   (not defined: too little data)" if np.isnan(p_sex) else ""))
print(f"Era × detection rate   : p={p_era:.3f}")
print(f"Assay × detection rate : p={p_assay:.3f}")

# Counts only: the console log is not a place for record-level rows.
print(f"\n=== Variants outside the panel regions ({ASSAY_EXT} only, "
      f"n={len(ext_ids)} individuals) ===")
print(extra.groupby(["variant", "exon"]).size().rename("observations").to_string()
      if len(extra) else "  none")

print(f"\nWritten to {OUT}/")
