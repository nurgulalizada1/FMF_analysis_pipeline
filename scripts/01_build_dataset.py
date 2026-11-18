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

Validation (src/mefv/validation.py, rules in docs/validation.md)
    The three raw tables are checked before anything is built. Every finding
    has a stable rule identifier and is reported as counts only; no identifier,
    label or cell value is printed or stored. Nothing is repaired, normalised,
    merged or dropped (apart from rows that are entirely empty, which are
    counted).
      * fatal findings stop the stage;
      * exact duplicates, conflicting genotypes and labels that differ by case
        or whitespace (V5, V6a, V6b, V7) withhold the outputs until the owner
        has recorded an explicit disposition (``--owner-disposition FILE``,
        stored outside the repository). Naming a rule with
        ``--accept-reviewed-findings`` only marks it as seen: the status is
        ``acknowledged_unresolved`` and no final output is written;
      * ``--diagnostic-only`` runs the checks, writes only the validation
        report and exits 0 (no blocking finding) or 2.

Outputs -> data/processed/{patients.csv, variants_long.csv, validation_report.json}
"""

# Make the package importable when the script is run without installing it.
import pathlib as _pathlib
import sys as _sys
_sys.path.insert(0, str(_pathlib.Path(__file__).resolve().parents[1] / "src"))

from mefv import config
from mefv.utils import _leading_comment_lines as _banner

import argparse
import json
import os
import sys

import pandas as pd

from mefv.disclosure import suppress_frame, suppressed_counts
from mefv.disclosure import text as _n
from mefv.eligibility import RESTORED_CONFIRMED, RESTORED_UNCONFIRMED
from mefv.validation import (
    ACKNOWLEDGEABLE,
    DataIntegrityError,
    digest_inputs,
    digest_config,
    load_dispositions,
    ruleset_sha256,
    canonical_id as _canonical_id,
    validate_built_tables,
    validate_raw_inputs,
)

_ap = argparse.ArgumentParser(description="Build the analysis tables from the raw export.")
_ap.add_argument("--diagnostic-only", action="store_true",
                 help="validate the raw tables and write only the validation report")
_ap.add_argument("--accept-reviewed-findings", nargs="*", default=None, metavar="RULE",
                 help=f"mark blocking findings ({', '.join(ACKNOWLEDGEABLE)}) as seen; this does "
                      "NOT release the run (see --owner-disposition); also read from "
                      "MEFV_ACCEPT_REVIEWED_FINDINGS (comma-separated)")
_ap.add_argument("--owner-disposition", default=None, metavar="FILE",
                 help="private JSON file with the owner's explicit disposition of each blocking "
                      "finding (outside the repository); also MEFV_OWNER_DISPOSITION")
_args = _ap.parse_args()
_accepted = set(_args.accept_reviewed_findings or [])
_accepted |= {r.strip() for r in os.environ.get("MEFV_ACCEPT_REVIEWED_FINDINGS", "").split(",")
              if r.strip()}
_DISPOSITION_PATH = _args.owner_disposition or os.environ.get("MEFV_OWNER_DISPOSITION") or None

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
    return _canonical_id(x, ID_PREFIX, ID_WIDTH)


aug = read(AUG)
var = read(VAR)
pat = read(PAT)

# ---- validate the raw tables before building anything --------------------
REPORT_PATH = config.VALIDATION_REPORT


def _write_report(rep):
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(json.dumps(rep.to_dict(), indent=1, sort_keys=True) + "\n")


def _withdraw_earlier_outputs(remove_report=False):
    """A stopped run must not leave the previous run's tables looking current.

    patients.csv and variants_long.csv of an earlier run are stage-01 products
    that this run supersedes; they are removed so nothing downstream can read
    them. (Later stages also refuse tables whose report did not release them.)
    """
    gone = 0
    for f in (config.PATIENTS, config.VARIANTS) + ((REPORT_PATH,) if remove_report else ()):
        if f.exists():
            f.unlink()
            gone += 1
    if gone:
        print(f"removed {gone} file(s) from an earlier run: no output of that run is current",
              file=sys.stderr)


def _stop(code, message, remove_report=False):
    _withdraw_earlier_outputs(remove_report)
    print(message, file=sys.stderr)
    raise SystemExit(code)


def _load_owner_disposition():
    """Entries and digest of the owner's disposition file, or none. Fails closed, value-free."""
    if not _DISPOSITION_PATH:
        return None, "", ""
    path = _pathlib.Path(_DISPOSITION_PATH).expanduser().resolve()
    if path.is_relative_to(_pathlib.Path(__file__).resolve().parents[1]):
        _stop(2, "stage 01 stopped: the owner disposition file must be stored outside the "
                 "repository", True)
    if not path.is_file():
        _stop(2, "stage 01 stopped: the owner disposition file was not found", True)
    try:
        return load_dispositions(str(path))
    except DataIntegrityError as exc:
        _stop(2, f"stage 01 stopped: {exc}", True)


_dispositions, _disposition_digest, _disposition_binding = _load_owner_disposition()


_report = validate_raw_inputs(aug, var, pat, design=DESIGN, id_prefix=ID_PREFIX,
                              id_width=ID_WIDTH, nomenclature_fixes=VARIANT_FIX)
# The release of a blocked run is bound to these exact files (see docs/validation.md).
_report.input_digests, _report.input_sha256 = digest_inputs(
    {"method": AUG, "variants": VAR, "patients": PAT})
_report.config_sha256 = digest_config(_PARAMS)
print(_report.render())
print("digests an owner disposition must carry (copy them from validation_report.json):")
print(f"  input_sha256   {_report.input_sha256}")
print(f"  config_sha256  {_report.config_sha256}")
print(f"  ruleset_sha256 {ruleset_sha256()}")
if _args.diagnostic_only:
    # Diagnostics only: acknowledgements and dispositions do not apply and nothing else is
    # written. The report says "not released" and replaces any earlier one, so tables from an
    # earlier run cannot be analysed after it; re-run stage 01 to regenerate them.
    _report.diagnostic_only = True
    _write_report(_report)
    print(f"\ndiagnostic-only run; validation report -> {REPORT_PATH}")
    raise SystemExit(0 if _report.status in ("passed", "passed_with_findings") else 2)
# (status "owner_disposition_recorded" cannot occur here: dispositions are not applied)
try:
    _report.enforce(_accepted, _dispositions, _disposition_digest, _disposition_binding)
except DataIntegrityError as exc:
    _write_report(_report)
    _stop(2, f"stage 01 stopped: {exc}")

# Rows that are entirely empty (trailing separators) carry no record. They were
# dropped silently before; they are now counted in the report (P1a, V1a, M3).
def _drop_empty(df, id_col):
    empty = df.apply(lambda c: c.astype(str).str.strip() == "").all(axis=1)
    return df[~(empty & (df[id_col].astype(str).str.strip() == ""))].copy()


aug, var, pat = _drop_empty(aug, "Patient_ID"), _drop_empty(var, "patient_id"), _drop_empty(pat, "patient_id")

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
# An individual who is in the method table but not in the patients table has no
# testing record. Whether they were tested, and for which regions, is not
# something the pipeline can infer: until the data owner confirms it
# (study_parameters: restored_individuals.coverage_confirmed = true) they are
# restored with an empty ``exons_tested`` and an undetermined detection
# status, so that every eligibility-aware analysis leaves them out of its
# denominators instead of counting them as tested-and-negative.
_RESTORED = _PARAMS.get("restored_individuals") or {}
_RESTORED_CONFIRMED = _RESTORED.get("coverage_confirmed", False) is True
missing = sorted(set(aug.patient_id) - set(pat.patient_id))
for pid in missing:
    r = aug[aug.patient_id == pid].iloc[0]
    if _RESTORED_CONFIRMED:
        extra = dict(n_variants="0", variant_detected="no",
                     exons_tested=",".join(DESIGN.panel_regions),
                     coverage_source=RESTORED_CONFIRMED)
    else:
        extra = dict(n_variants="", variant_detected="undetermined", exons_tested="",
                     coverage_source=RESTORED_UNCONFIRMED)
    pat = pd.concat([pat, pd.DataFrame([dict(
        patient_id=pid, year_requested=r.YEAR, sex=r.sex, **extra)])], ignore_index=True)

# ---- confirmed nomenclature corrections -----------------------------------
RENAME = VARIANT_FIX
# Count before renaming. The original script counted afterwards, when no
# mistyped label remained, so it always reported zero corrections.
n_fixed = int(var.variant.isin(RENAME).sum())
_n_rules_applied = 0
for old, new in RENAME.items():
    n = int((var.variant == old).sum())
    if n:
        var.loc[var.variant == old, "variant"] = new
        _n_rules_applied += 1
print(f"nomenclature corrections applied: {_n_rules_applied} rule(s), {_n(n_fixed)} row(s)")

# ---- attach ---------------------------------------------------------------
var["variant"] = var.variant.replace(VARIANT_FIX)

for df in (pat, var):
    df["assay"] = df.patient_id.map(meth.assay)
    df["assay_source"] = df.patient_id.map(meth.assay_source)

pat = pat.sort_values("patient_id").reset_index(drop=True)

# ---- post-assembly checks (formerly assert statements) -------------------
validate_built_tables(pat, DESIGN.assays, _PARAMS["expected_cohort_size"], _report)
if not set(var.patient_id) <= set(pat.patient_id):
    _report.add("X1", records=int((~var.patient_id.isin(set(pat.patient_id))).sum()))
try:
    _report.enforce(_accepted, _dispositions, _disposition_digest, _disposition_binding)
except DataIntegrityError as exc:
    _write_report(_report)
    _stop(2, f"stage 01 stopped: {exc}")
_write_report(_report)

pat.to_csv(f"{OUT}/patients.csv", sep=";", index=False, encoding="utf-8-sig")
var.to_csv(f"{OUT}/variants_long.csv", sep=";", index=False, encoding="utf-8-sig")

# ---- report ---------------------------------------------------------------
print(f"variant labels fixed: {_n(n_fixed)} ({len(RENAME)} documented correction rule(s))")
print("(counts below 5 are shown as '<5')")
print(f"patients          : {_n(len(pat))}")
print(f"variant rows      : {_n(len(var))}")
print(f"restored patients : {_n(len(missing))} "
      f"({'coverage confirmed by owner' if _RESTORED_CONFIRMED else 'testing record unconfirmed -> undetermined'})")
print(f"\nassay:\n{suppressed_counts(pat.assay.value_counts()).to_string()}")
print(f"\nassay_source:\n{suppressed_counts(pat.assay_source.value_counts()).to_string()}")
pat["year"] = pat.year_requested.astype(int)
print(f"\nassay x year:\n{suppress_frame(pd.crosstab(pat.year, pat.assay)).to_string()}")
print(f"\nvariant detected:\n{suppressed_counts(pat.variant_detected.value_counts()).to_string()}")

print(f"\nvalidation: {_report.status} (report: {REPORT_PATH.name})")
