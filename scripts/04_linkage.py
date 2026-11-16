#!/usr/bin/env python3
"""
04_linkage.py -- allelic architecture of the MEFV locus.

  A. the primary variant pair, stratified by recording era and by assay
  B. all pairwise associations with Benjamini-Hochberg FDR
  C. candidate complex alleles (perfect co-occurrence)
  D. benign-polymorphism haplotype blocks and their link to the ascertained locus
  E. patients carrying >=3 FMF-associated alleles

Phase was not determined experimentally. D' and haplotype frequencies are
population-level estimates from unphased genotypes: they are compatible with,
but do not prove, a cis configuration in any individual patient.
"""

# Make the package importable when the script is run without installing it.
import pathlib as _pathlib
import sys as _sys
_sys.path.insert(0, str(_pathlib.Path(__file__).resolve().parents[1] / "src"))

from mefv import config

import os
import sys
import itertools
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mefv import DOSE, bh_fdr, genotype_vector, ld_stats, load

_GROUPS = config.variant_groups()
SYNONYMOUS = _GROUPS["synonymous"]
BENIGN_COMMON = _GROUPS["benign_common"]
_PLAN = config.linkage_plan()

PAT = str(config.PATIENTS)
VAR = str(config.VARIANTS)
OUT = str(config.TABLES)
os.makedirs(OUT, exist_ok=True)
_DESIGN = config.study_design()
ERA_EARLY, ERA_LATE = _DESIGN.primary.label, _DESIGN.sensitivity.label
PANEL_EXONS = set(_DESIGN.panel_regions)   # regions tested in every individual
ASSAY_PANEL, ASSAY_EXT = _DESIGN.panel_assay, _DESIGN.extended_assay
RNG = np.random.default_rng(config.SEED_BOOTSTRAP)

pat, var = load(PAT, VAR, _DESIGN)
var = var[var.exon.isin(PANEL_EXONS)].copy()

ids_all = list(pat.patient_id)
ids_early = list(pat[pat.era == ERA_EARLY].patient_id)
ids_late = list(pat[pat.era == ERA_LATE].patient_id)
ids_panel_assay = list(pat[pat.assay == ASSAY_PANEL].patient_id)
ids_ext_assay = list(pat[pat.assay == ASSAY_EXT].patient_id)


def gv(v, ids):
    return genotype_vector(var, v, ids)


def bootstrap_dprime(gA, gB, n_boot=config.BOOTSTRAP_REPLICATES):
    g = np.array(list(zip(gA, gB)))
    n = len(g)
    vals = []
    for _ in range(n_boot):
        s = g[RNG.integers(0, n, n)]
        r = ld_stats(s[:, 0].tolist(), s[:, 1].tolist())
        if not np.isnan(r["Dprime"]):
            vals.append(r["Dprime"])
    return tuple(np.percentile(vals, [2.5, 97.5])) if vals else (np.nan, np.nan)


# ------------------------------------------- A. the primary variant pair ----
def pair_report(vA, vB, ids, label, boot=True):
    gA, gB = gv(vA, ids), gv(vB, ids)
    r = ld_stats(gA, gB)
    a = sum(1 for x, y in zip(gA, gB) if x > 0 and y > 0)
    b = sum(1 for x, y in zip(gA, gB) if x > 0 and y == 0)
    c = sum(1 for x, y in zip(gA, gB) if x == 0 and y > 0)
    d = sum(1 for x, y in zip(gA, gB) if x == 0 and y == 0)
    _, p = stats.fisher_exact([[a, b], [c, d]])
    lo, hi = bootstrap_dprime(gA, gB) if boot else (np.nan, np.nan)
    frac = r["h"]["11"] / r["pB"] if r["pB"] > 0 else np.nan
    return dict(stratum=label, n=len(ids), freq_A=round(r["pA"], 4),
                freq_B=round(r["pB"], 4),
                hap_AB=round(r["h"]["11"], 4), hap_Ab=round(r["h"]["10"], 4),
                hap_aB=round(r["h"]["01"], 4), hap_ab=round(r["h"]["00"], 4),
                Dprime=round(r["Dprime"], 4),
                Dprime_95CI=f"{lo:.3f}–{hi:.3f}" if not np.isnan(lo) else "",
                r2=round(r["r2"], 4), carriers_both=a, A_only=b, B_only=c,
                neither=d, fisher_p=f"{p:.2e}",
                pct_B_alleles_on_A_background=round(frac * 100, 1))


strata = [(f"{ERA_EARLY} (primary)", ids_early),
          (f"{ERA_LATE} (sensitivity)", ids_late),
          (ASSAY_PANEL, ids_panel_assay), (ASSAY_EXT, ids_ext_assay),
          ("full cohort", ids_all)]
_A, _B = _PLAN["primary_pair"]
t4 = pd.DataFrame([pair_report(_A, _B, ids, lab) for lab, ids in strata])
t4.insert(0, "variant_pair", f"{_A} / {_B}")
t4.to_csv(f"{OUT}/table4_primary_pair.csv", sep=";", index=False, encoding="utf-8-sig")


def crosstab(vA, vB, ids):
    lbl = {0: "wild type", 1: "heterozygous", 2: "homozygous"}
    ct = pd.crosstab(pd.Series([lbl[x] for x in gv(vA, ids)], name=vA),
                     pd.Series([lbl[x] for x in gv(vB, ids)], name=vB))
    order = ["homozygous", "heterozygous", "wild type"]
    return ct.reindex(index=[o for o in order if o in ct.index],
                      columns=[o for o in reversed(order) if o in ct.columns]).fillna(0).astype(int)


ct_early = crosstab(_A, _B, ids_early)
ct_early.to_csv(f"{OUT}/table4b_crosstab_primary.csv", sep=";", encoding="utf-8-sig")


# ------------------------------------------- B. pairwise associations -------
def pairwise(ids, variants, label, min_carriers=2):
    rows = []
    for vA, vB in itertools.combinations(sorted(variants), 2):
        gA, gB = gv(vA, ids), gv(vB, ids)
        a = sum(1 for x, y in zip(gA, gB) if x > 0 and y > 0)
        b = sum(1 for x, y in zip(gA, gB) if x > 0 and y == 0)
        c = sum(1 for x, y in zip(gA, gB) if x == 0 and y > 0)
        d = sum(1 for x, y in zip(gA, gB) if x == 0 and y == 0)
        if (a + b) < min_carriers or (a + c) < min_carriers:
            continue
        _, p = stats.fisher_exact([[a, b], [c, d]])
        r = ld_stats(gA, gB)
        rows.append(dict(cohort=label, variant_A=vA, variant_B=vB,
                         carriers_A=a + b, carriers_B=a + c, both=a,
                         Dprime=round(r["Dprime"], 4) if not np.isnan(r["Dprime"]) else "",
                         r2=round(r["r2"], 4) if not np.isnan(r["r2"]) else "",
                         fisher_p_raw=p,
                         perfect_cooccurrence=(a == a + b == a + c)))
    df = pd.DataFrame(rows)
    if len(df):
        df["fdr_q"] = bh_fdr(df.fisher_p_raw.values)
        df["fisher_p"] = df.fisher_p_raw.map(lambda x: f"{x:.2e}")
        df["fdr_q"] = df.fdr_q.map(lambda x: f"{x:.2e}")
        df = df.drop(columns="fisher_p_raw").sort_values("both", ascending=False)
    return df


nonsyn = sorted(set(var[~var.variant.isin(SYNONYMOUS)].variant))
pw_ns = pairwise(ids_all, nonsyn, f"non-synonymous, full cohort (n={len(ids_all)})")
syn = sorted(set(var[var.variant.isin(SYNONYMOUS)].variant))
pw_syn = pairwise(ids_early, syn, f"benign polymorphisms, {ERA_EARLY} (n={len(ids_early)})")
pd.concat([pw_ns, pw_syn], ignore_index=True).to_csv(
    f"{OUT}/supp_pairwise_associations.csv", sep=";", index=False, encoding="utf-8-sig")

cands = pw_ns[pw_ns.perfect_cooccurrence] if len(pw_ns) else pd.DataFrame()
if len(cands):
    cands.to_csv(f"{OUT}/table5_complex_alleles.csv", sep=";", index=False, encoding="utf-8-sig")


# ------------------------------------------- D. haplotype blocks ------------
blocks = _PLAN["haplotype_blocks"]
rows = []
for name, members in blocks.items():
    for vA, vB in itertools.combinations(members, 2):
        gA, gB = gv(vA, ids_early), gv(vB, ids_early)
        r = ld_stats(gA, gB)
        rows.append(dict(comparison=f"within {name}", variant_A=vA, variant_B=vB,
                         discordant_genotypes=sum(1 for x, y in zip(gA, gB) if x != y),
                         n=len(ids_early), Dprime=round(r["Dprime"], 4),
                         r2=round(r["r2"], 4)))
for vA, _partners in _PLAN["cross_locus_pairs"].items():
    for vB in _partners:
        r = ld_stats(gv(vA, ids_early), gv(vB, ids_early))
        rows.append(dict(comparison="benign block vs ascertained locus",
                         variant_A=vA, variant_B=vB, discordant_genotypes="",
                         n=len(ids_early), Dprime=round(r["Dprime"], 4),
                         r2=round(r["r2"], 4)))
t_blocks = pd.DataFrame(rows)
t_blocks.to_csv(f"{OUT}/supp_haplotype_blocks.csv", sep=";", index=False, encoding="utf-8-sig")


# ------------------------------------------- E. high allele burden ----------
b = var[(~var.variant.isin(SYNONYMOUS)) & (~var.variant.isin(BENIGN_COMMON))].copy()
b["alleles"] = b.zygosity.map(DOSE).fillna(1)
burden = b.groupby("patient_id").alleles.sum()
high = burden[burden >= 3].sort_values(ascending=False)
t_high = pd.DataFrame([
    dict(patient_id=pid, year=int(b[b.patient_id == pid].year.iloc[0]),
         assay=b[b.patient_id == pid].assay.iloc[0], total_alleles=int(high[pid]),
         genotype="; ".join(f"{r.variant}({r.zygosity})"
                            for _, r in b[b.patient_id == pid].iterrows()))
    for pid in high.index])
t_high.to_csv(f"{OUT}/supp_high_allele_burden.csv", sep=";", index=False, encoding="utf-8-sig")


# ------------------------------------------------------------------ print ---
pd.set_option("display.width", 220)
print(f"=== A. {_A} × {_B} ===")
print(t4[["stratum", "n", "freq_A", "freq_B", "hap_AB", "hap_aB", "Dprime",
          "Dprime_95CI", "r2", "fisher_p", "pct_B_alleles_on_A_background"]].to_string(index=False))

print(f"\n--- genotype cross-tabulation, {ERA_EARLY} (n={len(ids_early)}) ---")
print(ct_early.to_string())

cols = ["variant_A", "variant_B", "carriers_A", "carriers_B", "both",
        "Dprime", "r2", "fisher_p", "fdr_q", "perfect_cooccurrence"]
print("\n=== B. strongest pairwise associations (non-synonymous) ===")
print(pw_ns.head(10)[cols].to_string(index=False))

print("\n=== C. candidate complex alleles (perfect co-occurrence) ===")
print(cands[cols].to_string(index=False) if len(cands) else "  none")

print(f"\n=== D. haplotype blocks, {ERA_EARLY} ===")
print(t_blocks.to_string(index=False))

print(f"\n=== E. patients with ≥3 FMF-associated alleles (n={len(t_high)}) ===")
print(t_high.to_string(index=False))
print(f"\nWritten to {OUT}/")
