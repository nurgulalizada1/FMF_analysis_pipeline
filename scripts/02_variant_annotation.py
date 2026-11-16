#!/usr/bin/env python3
"""
02_variant_annotation.py -- build the variant annotation dictionary.

The dictionary has two parts. This stage generates the twelve columns that are
transcribed from the external variant database, then joins the reviewed
overlay that carries functional-evidence fields and manual scientific
corrections. The overlay is read, never written: a routine re-run of this
stage cannot overwrite a reviewed value. A reviewed value that deviates from
its external source must be recorded in the annotation decision log, or the
merge fails. See docs/annotation_schema.md.

The source annotations (MEFV_SOURCE_ANNOTATIONS) follow the record structure of
the Infevers database: name, HGVS descriptions on NM_000243.3, exon, rsID,
current classification and, where it changed at the database's 2025
re-evaluation, the previous one. When and from which release they were
transcribed is recorded with the private source file, not here.
"""

# Make the package importable when the script is run without installing it.
import pathlib as _pathlib
import sys as _sys
_sys.path.insert(0, str(_pathlib.Path(__file__).resolve().parents[1] / "src"))

from mefv import annotation, config
from mefv.utils import _leading_comment_lines as _banner

import os
import numpy as np
import pandas as pd

# The merged dictionary is written where every later stage reads it from.
DICTIONARY = config.DICTIONARY
DICTIONARY.parent.mkdir(parents=True, exist_ok=True)

SOURCE = config.SOURCE_ANNOTATIONS
if not SOURCE.exists():
    raise SystemExit(
        f"source annotation table not found at {SOURCE}.\n"
        "This table is the transcription of the external variant database for the "
        "variants observed in the study. It is not kept in the repository, because "
        "the set of variants it lists is itself a study result. Copy "
        "data/external/source_annotations.example.csv and fill it in, or set "
        "MEFV_SOURCE_ANNOTATIONS."
    )
d = pd.read_csv(SOURCE, sep=";", dtype=str, keep_default_na=False, encoding="utf-8-sig", skiprows=_banner(SOURCE))
print(f"source annotations   <- {SOURCE} ({len(d)} variants)")


d.insert(2, "hgvs_full", "NM_000243.3:" + d.hgvs_c)

TIER = {"Pathogenic": "P", "Likely pathogenic": "LP",
        "Uncertain significance": "VUS",
        "Uncertain significance (provisional)": "VUS",
        "Likely benign": "LB", "Benign": "B", "not in Infevers": "unclassified"}
d["tier"] = d.infevers_class.map(TIER)
d["reclassified_2025"] = d.previous_class.str.len() > 0

# ---- generated annotations -------------------------------------------------
generated_path = config.ANNOTATIONS_GENERATED
generated_path.parent.mkdir(parents=True, exist_ok=True)
d.to_csv(generated_path, sep=";", index=False, encoding="utf-8-sig")
print(f"generated annotations -> {generated_path} ({len(d)} variants, {len(d.columns)} columns)")

# ---- join the reviewed overlay ---------------------------------------------
overlay = annotation.load_overlay(config.ANNOTATION_OVERLAY)
log = annotation.load_decision_log(config.ANNOTATION_DECISION_LOG)
merged = annotation.merge(d, overlay, log)
annotation.validate(merged)
print(f"reviewed overlay     <- {config.ANNOTATION_OVERLAY} ({len(overlay)} variants)")
print(f"decision log         <- {config.ANNOTATION_DECISION_LOG} ({len(log)} recorded decision(s))")
for r in log.itertuples():
    print(f"    {r.variant} / {r.column}: {r.previous_value!r} -> {r.accepted_value!r} "
          f"[{r.review_status}]")
_open = [r for r in log.itertuples() if r.review_status != "accepted"]
if _open:
    print(f"NOTE: {len(_open)} recorded decision(s) are provisional or awaiting a second "
          "review. Results that depend on them are provisional too.")

merged.to_csv(DICTIONARY, sep=";", index=False, encoding="utf-8-sig")
merged.to_csv(DICTIONARY.with_suffix(".tsv"), sep="\t", index=False, encoding="utf-8")
print(f"merged dictionary    -> {DICTIONARY} "
      f"({len(merged)} variants, {len(merged.columns)} columns)")

d = merged
print(f"\n{len(d)} variants")
print(d.tier.value_counts().to_string())
print(f"\nreclassified since 2023: {int(d.reclassified_2025.sum())}")
print("\n=== reclassified ===")
print(d[d.reclassified_2025][["variant", "previous_class", "infevers_class"]].to_string(index=False))

# codon vs c. position sanity check
import re
print("\n=== codon consistency check ===")
bad = []
for _, r in d.iterrows():
    m = re.search(r"c\.(\d+)", r.hgvs_c)
    p = re.search(r"\((\D+?)(\d+)", r.hgvs_p)
    if m and p:
        codon_from_c = (int(m.group(1)) + 2) // 3
        if codon_from_c != int(p.group(2)):
            bad.append((r.variant, r.hgvs_c, r.hgvs_p, codon_from_c))
print("  mismatches:", bad if bad else "none — all c. positions match their codon")
