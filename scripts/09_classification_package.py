#!/usr/bin/env python3
"""
09_acmg_package.py -- assemble the evidence package for ACMG classification
and write two blank assessment forms for independent assessors.

Evidence assembled per variant: HGVS, exon, cohort counts and frequency,
gnomAD frequencies (global, grpmax filtering AF, Middle Eastern, highest group),
homozygote count, enrichment ratio, ClinVar classification and review status,
Infevers classification, SpeckSeq functional class, REVEL, CADD.

Assessors fill the ACMG code columns independently and blinded to each other.
Agreement is quantified with Cohen's kappa in 10_acmg_concordance.py.
"""

# Make the package importable when the script is run without installing it.
import pathlib as _pathlib
import sys as _sys
_sys.path.insert(0, str(_pathlib.Path(__file__).resolve().parents[1] / "src"))

from mefv import config
from mefv.utils import _leading_comment_lines as _banner

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

DICT = str(config.DICTIONARY)
VAR = str(config.VARIANTS)
PAT = str(config.PATIENTS)
TAB = str(config.TABLES)
OUT = str(config.CLASSIFICATION)
os.makedirs(OUT, exist_ok=True)


def read(p, sep=";"):
    d = pd.read_csv(p, sep=sep, encoding="utf-8-sig", dtype=str, keep_default_na=False, skiprows=_banner(p))
    return d.loc[:, [c for c in d.columns if not c.startswith("Unnamed")]]


dic, var, pat = read(DICT), read(VAR), read(PAT)
gn = pd.read_csv(f"{TAB}/table6_gnomad_comparison.csv", sep=";", encoding="utf-8-sig")
N = len(pat)

G = gn.set_index("variant")
rows = []
for _, d in dic.iterrows():
    v = d.variant
    s = var[var.variant == v]
    het = int((s.zygosity == "Het").sum())
    hom = int((s.zygosity == "Hom").sum())
    g = G.loc[v] if v in G.index else None

    def gv(col, dp=None):
        if g is None or pd.isna(g.get(col)):
            return ""
        x = g.get(col)
        return round(float(x), dp) if dp is not None and x != "" else x

    rows.append(dict(
        variant=v,
        hgvs_c=f"NM_000243.3:{d.hgvs_c}" if d.hgvs_c else "",
        hgvs_p=d.hgvs_p, exon=d.exon_infevers, rsid=d.rsid,
        cohort_het=het, cohort_hom=hom, cohort_carriers=het + hom,
        cohort_af=round((2 * hom + het) / (2 * N), 5),
        gnomad_global_af=gv("gnomad_global"),
        gnomad_grpmax_filtering_af=gv("grpmax_filtering_af"),
        gnomad_middle_eastern_af=gv("gnomad_mid"),
        gnomad_highest_group=gv("highest_group"),
        gnomad_highest_af=gv("highest_af"),
        gnomad_n_homozygotes=gv("gnomad_n_hom"),
        enrichment_vs_highest=gv("x_vs_highest"),
        BA1_flag=gv("BA1"), BS1_flag=gv("BS1"),
        clinvar_class=gv("clinvar"), clinvar_stars=gv("clinvar_stars"),
        infevers_class=d.infevers_class,
        infevers_previous=d.previous_class,
        speckseq_functional=d.speckseq,
        other_functional=d.other_functional,
        revel=gv("revel"), cadd=gv("cadd"),
        infevers_comment=d.comment,
    ))

pkg = pd.DataFrame(rows)
pkg["_sort"] = pd.to_numeric(pkg.cohort_carriers, errors="coerce")
pkg = pkg.sort_values("_sort", ascending=False).drop(columns="_sort")
pkg.to_csv(f"{OUT}/acmg_evidence_package.csv", sep=";", index=False, encoding="utf-8-sig")
pkg.to_csv(f"{OUT}/acmg_evidence_package.tsv", sep="\t", index=False,
           encoding="utf-8", lineterminator="\n")

# ---------------------------------------------------------- blank forms -----
CODES = ["PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5",
         "PM6", "PP1", "PP2", "PP3", "PP4", "PP5",
         "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5",
         "BP6", "BP7"]
NOT_APPLICABLE = {"PS2", "PM3", "PM6", "PP1", "PP4", "PP5", "BS4", "BP2", "BP5", "BP6"}

for assessor in ("assessor_1", "assessor_2"):
    f = pkg[["variant", "hgvs_c", "hgvs_p", "exon", "cohort_carriers", "cohort_af",
             "gnomad_grpmax_filtering_af", "gnomad_n_homozygotes",
             "enrichment_vs_highest", "clinvar_class", "clinvar_stars",
             "infevers_class", "speckseq_functional", "revel", "cadd"]].copy()
    for c in CODES:
        f[c] = "n/a" if c in NOT_APPLICABLE else ""
    f["strength_modifiers"] = ""
    f["final_classification"] = ""
    f["rationale"] = ""
    f["date"] = ""
    f.to_csv(f"{OUT}/acmg_form_{assessor}.tsv", sep="\t", index=False,
             encoding="utf-8", lineterminator="\n")

print(f"variants in package : {len(pkg)}")
print(f"with gnomAD data    : {(pkg.gnomad_global_af != '').sum()}")
print(f"with SpeckSeq data  : {(pkg.speckseq_functional != '').sum()}")
print(f"BA1 flagged         : {(pkg.BA1_flag == 'yes').sum()}")
print(f"BS1 flagged         : {(pkg.BS1_flag == 'yes').sum()}")
print(f"\nwritten to {OUT}/")
