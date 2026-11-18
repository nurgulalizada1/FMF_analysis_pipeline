#!/usr/bin/env python3
"""
09_acmg_package.py -- assemble the evidence package for ACMG classification
and write two blank assessment forms for independent assessors.

Evidence assembled per variant: HGVS, exon, cohort counts and frequency,
gnomAD frequencies (global, grpmax filtering AF, Middle Eastern, highest group),
homozygote count, enrichment ratio, ClinVar classification and review status,
Infevers classification, SpeckSeq functional class, REVEL, CADD.

Assessors fill the ACMG code columns independently and blinded to each other.
Agreement is quantified with Cohen's kappa; that script is not part of this
repository.

The package assigns no ACMG code and no classification. Cohort frequencies use
each variant's eligible population (see src/mefv/eligibility.py); BA1/BS1 are
"not_evaluable" when the reference has no grpmax value. The thresholds and the
pre-filled "n/a" codes are provisional rules whose provenance is not documented;
they are written to acmg_rule_provenance.tsv beside the package.
"""

# Make the package importable when the script is run without installing it.
import pathlib as _pathlib
import sys as _sys
_sys.path.insert(0, str(_pathlib.Path(__file__).resolve().parents[1] / "src"))

from mefv import config
from mefv.acmg import CODES, NOT_APPLICABLE, provenance_table
from mefv.eligibility import (
    cohort_frequency, known_regions, primary_period_scope, resolve_region, variant_regions,
)
from mefv.utils import _leading_comment_lines as _banner
from mefv.utils import require_released

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

DICT = str(config.DICTIONARY)
VAR = str(config.VARIANTS)
PAT = str(config.PATIENTS)
require_released(PAT)   # refuse tables that a blocked or diagnostic-only stage 01 did not release
TAB = str(config.TABLES)
OUT = str(config.CLASSIFICATION)
os.makedirs(OUT, exist_ok=True)


def read(p, sep=";"):
    d = pd.read_csv(p, sep=sep, encoding="utf-8-sig", dtype=str, keep_default_na=False, skiprows=_banner(p))
    return d.loc[:, [c for c in d.columns if not c.startswith("Unnamed")]]


dic, var, pat = read(DICT), read(VAR), read(PAT)
gn = pd.read_csv(f"{TAB}/table6_gnomad_comparison.csv", sep=";", encoding="utf-8-sig")
N = len(pat)
OBS_REGIONS = variant_regions(var)
VOCAB = known_regions(pat)
DICT_REGION = dict(zip(dic.variant, dic.exon_infevers))
_DESIGN = config.study_design()
_PERIOD_RESTRICTED = config.variant_groups()["synonymous"]
_PRIMARY = primary_period_scope(pat, _DESIGN)

G = gn.set_index("variant")
rows = []
for _, d in dic.iterrows():
    v = d.variant
    region, source = resolve_region(v, OBS_REGIONS, DICT_REGION, VOCAB)
    if v in _PERIOD_RESTRICTED:     # recorded only in the primary period: same time window
        c = cohort_frequency(pat, var, v, region, source, _PRIMARY, _DESIGN.primary.label)
    else:
        c = cohort_frequency(pat, var, v, region, source)
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
        cohort_het=c["heterozygous"], cohort_hom=c["homozygous"], cohort_carriers=c["carriers"],
        cohort_scope=c["scope"], cohort_n_total=c["n_cohort"], cohort_n_in_scope=c["n_in_scope"],
        cohort_n_tested=c["n_tested"],
        cohort_n_not_tested=c["n_not_tested"], cohort_n_undetermined=c["n_undetermined"],
        cohort_n_eligible=c["n_eligible"], cohort_allele_denominator=c["allele_denominator"],
        cohort_n_detected_ineligible=c["n_detected_ineligible"],
        cohort_eligibility_status=c["eligibility_status"],
        cohort_af=round(c["allele_freq"], 5) if c["allele_freq"] == c["allele_freq"] else "",
        gnomad_global_af=gv("gnomad_global"),
        gnomad_grpmax_filtering_af=gv("grpmax_filtering_af"),
        gnomad_middle_eastern_af=gv("gnomad_mid"),
        gnomad_highest_group=gv("highest_group"),
        gnomad_highest_af=gv("highest_af"),
        gnomad_n_homozygotes=gv("gnomad_n_hom"),
        enrichment_vs_highest=gv("x_vs_highest"),
        enrichment_status=gv("x_vs_highest_status"),
        reference_source=gv("reference_source"), reference_version=gv("reference_version"),
        reference_retrieved=gv("reference_retrieved"),
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
# The thresholds, the code list and the pre-filled "n/a" set live in mefv.acmg with
# their provenance status; they are provisional, not validated, rules.
provenance_table().to_csv(f"{OUT}/acmg_rule_provenance.tsv", sep="\t", index=False,
                          encoding="utf-8", lineterminator="\n")

for assessor in ("assessor_1", "assessor_2"):
    f = pkg[["variant", "hgvs_c", "hgvs_p", "exon", "cohort_carriers", "cohort_n_eligible",
             "cohort_eligibility_status", "cohort_af",
             "gnomad_grpmax_filtering_af", "gnomad_n_homozygotes",
             "enrichment_vs_highest", "enrichment_status", "BA1_flag", "BS1_flag",
             "reference_source", "reference_version", "clinvar_class", "clinvar_stars",
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
print(f"BA1 flagged         : {(pkg.BA1_flag == 'yes').sum()}"
      f"  (not evaluable: {(pkg.BA1_flag == 'not_evaluable').sum()})")
print(f"BS1 flagged         : {(pkg.BS1_flag == 'yes').sum()}"
      f"  (not evaluable: {(pkg.BS1_flag == 'not_evaluable').sum()})")
print(f"cohort frequency not evaluable: "
      f"{(pkg.cohort_eligibility_status != 'evaluable').sum()}")
print(f"\nwritten to {OUT}/")
