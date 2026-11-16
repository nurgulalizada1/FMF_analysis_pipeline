#!/usr/bin/env python3
"""
01_build_dataset.py -- assemble the two analysis tables from the laboratory export.

Joins the assay (Method) column from the method table onto the patient and
variant tables, restores any record present in the method table but missing
from the patient table, and applies confirmed nomenclature corrections.

Assay assignment, in order of precedence
    1. a record-level override confirmed by the laboratory
       (``assay_overrides`` in the study parameters);
    2. the value recorded in the export;
    3. if ``study_design.assay_inference`` is set and the record's year is no
       later than its ``until_year``: the assay named there. This encodes a
       documented laboratory statement, not a data-derived fact, and is
       labelled with ``source_label`` in the ``assay_source`` column so that it
       can be reported as such.
    Anything else stops the stage.

Patients in whom no variant was detected are retained in patients.csv with
n_variants = 0. They were tested and are part of the denominator.

Outputs -> data/processed/{patients.csv, variants_long.csv}
"""

# Make the package importable when the script is run without installing it.
import pathlib as _pathlib
import sys as _sys
_sys.path.insert(0, str(_pathlib.Path(__file__).resolve().parents[1] / "src"))

from mefv import config
from mefv.utils import _leading_comment_lines as _banner

import os
import pandas as pd

# Inputs are the raw laboratory export; outputs are the processed tables that
# every later stage consumes. The original script reused MEFV_PATIENTS and
# MEFV_VARIANTS for both roles, so setting either variable globally silently
# pointed this stage at its own output. The two roles now have distinct names.
AUG = str(config.RAW_AUG)
VAR = str(config.RAW_VARIANTS)
PAT = str(config.RAW_PATIENTS)
OUT = str(config.PROCESSED_DIR)

_PARAMS = config.load_study_parameters()
DESIGN = config.parse_study_design(_PARAMS)

# Records whose assay the laboratory corrected or confirmed individually. The
# identifiers are pseudonymous but the laboratory holds the key that links
# them to individuals, so they are not published here. They are read from the
# local study-parameters file; see study_parameters.example.json.
MANUAL = dict(_PARAMS["assay_overrides"])

# Nomenclature corrections confirmed against the external variant database and
# with the laboratory. Each maps a transcription error in the records to the
# correct variant label. The mapping names study variants, so it is read from
# the local study-parameters file together with its justification.
VARIANT_FIX = dict(_PARAMS["nomenclature_fixes"])
os.makedirs(OUT, exist_ok=True)


def read(path):
    d = pd.read_csv(path, sep=";", encoding="utf-8-sig", dtype=str,
                    keep_default_na=False, engine="python", skiprows=_banner(path))
    return d.loc[:, [c for c in d.columns if not c.startswith("Unnamed")]]


# The identifier scheme is a study parameter: the prefix and zero-padding of the
# record identifiers. Identifiers in the export are inconsistent in case and
# padding, so each is normalised to PREFIX + zero-padded number.
ID_PREFIX = str(_PARAMS["id_prefix"])
ID_WIDTH = int(_PARAMS["id_width"])


def canonical_id(x):
    s = str(x).strip().lower().replace(ID_PREFIX.lower(), "")
    return f"{ID_PREFIX}{int(s):0{ID_WIDTH}d}" if s.isdigit() else None


aug = read(AUG)
var = read(VAR)
pat = read(PAT)
pat = pat[pat.patient_id.str.strip() != ""].copy()

aug["patient_id"] = aug.Patient_ID.map(canonical_id)
aug["year"] = aug.YEAR.astype(int)

# ---- resolve Method -------------------------------------------------------
def resolve(row):
    if row.patient_id in MANUAL:
        return MANUAL[row.patient_id]
    m = row.Method.strip()
    if m:
        return m
    if _inferable(row.year):
        return DESIGN.inferred_assay
    return "not_recorded"


def _inferable(year):
    return DESIGN.inference_until_year is not None and year <= DESIGN.inference_until_year


aug["assay"] = aug.apply(resolve, axis=1)
# The provenance label follows the same precedence as resolve() above, so it
# always says where the assay value came from. (The original script checked
# the export field first, which labelled an overridden record "recorded".)
def provenance(pid, method, year):
    if pid in MANUAL:
        return "laboratory_confirmed"
    if method.strip():
        return "recorded"
    if _inferable(year):
        return DESIGN.inferred_source_label
    return "missing"


aug["assay_source"] = [provenance(pid, m, y)
                       for pid, m, y in zip(aug.patient_id, aug.Method, aug.year)]

meth = aug.set_index("patient_id")[["assay", "assay_source"]]

# ---- restore records present in the augmented export only ----------------
missing = sorted(set(aug.patient_id) - set(pat.patient_id))
for pid in missing:
    r = aug[aug.patient_id == pid].iloc[0]
    pat = pd.concat([pat, pd.DataFrame([dict(
        patient_id=pid, year_requested=r.YEAR, sex=r.sex, n_variants="0",
        variant_detected="no", exons_tested=",".join(DESIGN.panel_regions),
        coverage_source="lab_confirmed_baseline_panel")])], ignore_index=True)

# ---- confirmed nomenclature corrections -----------------------------------
RENAME = VARIANT_FIX
# Count before renaming. The original script counted afterwards, when no
# mistyped label remained, so it always reported zero corrections.
n_fixed = int(var.variant.isin(RENAME).sum())
for old, new in RENAME.items():
    n = int((var.variant == old).sum())
    if n:
        var.loc[var.variant == old, "variant"] = new
        print(f"renamed {old} -> {new} ({n} row(s))")

# ---- attach ---------------------------------------------------------------
var["variant"] = var.variant.replace(VARIANT_FIX)

for df in (pat, var):
    df["assay"] = df.patient_id.map(meth.assay)
    df["assay_source"] = df.patient_id.map(meth.assay_source)

pat = pat.sort_values("patient_id").reset_index(drop=True)
pat.to_csv(f"{OUT}/patients.csv", sep=";", index=False, encoding="utf-8-sig")
var.to_csv(f"{OUT}/variants_long.csv", sep=";", index=False, encoding="utf-8-sig")

# ---- report ---------------------------------------------------------------
print(f"variant labels fixed: {n_fixed} ({len(RENAME)} documented correction rule(s))")
print(f"patients          : {len(pat)}")
print(f"variant rows      : {len(var)}")
print(f"restored patients : {missing if missing else 'none'}")
print(f"\nassay:\n{pat.assay.value_counts().to_string()}")
print(f"\nassay_source:\n{pat.assay_source.value_counts().to_string()}")
pat["year"] = pat.year_requested.astype(int)
print(f"\nassay x year:\n{pd.crosstab(pat.year, pat.assay).to_string()}")
print(f"\nvariant detected:\n{pat.variant_detected.value_counts().to_string()}")

_expected = _PARAMS["expected_cohort_size"]
assert len(pat) == _expected, f"patient count should be {_expected}, got {len(pat)}"
assert pat.assay.notna().all(), "unassigned assay"
assert (pat.assay != "not_recorded").all(), "assay still missing"
assert set(pat.assay) <= set(DESIGN.assays), (
    f"assay labels {sorted(set(pat.assay) - set(DESIGN.assays))} are not among the "
    f"configured study_design.assays {list(DESIGN.assays)}")
assert set(var.patient_id) <= set(pat.patient_id), "orphan variant rows"
print("\nchecks passed")
