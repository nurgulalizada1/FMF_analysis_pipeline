"""Input/output helpers and genotype encoding.

The analysis uses a two-table design:

* ``patients``       one row per tested individual; this is the denominator
                     for every rate reported by the pipeline.
* ``variants_long``  one row per observed variant in one patient.

Individuals without an observed variant have no row in ``variants_long`` and
are therefore encoded as homozygous reference. This is only valid for regions
that were tested in every individual - the configured panel regions - which is
why analyses that rely on it restrict themselves to those regions; see
docs/limitations.md.
"""
from __future__ import annotations

import warnings
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from mefv.config import StudyDesign


def _leading_comment_lines(path: str | Path) -> int:
    """Count the ``#`` banner lines at the top of a file, if any.

    Synthetic data files carry a banner marking them as artificial. Only
    leading lines are considered, so a ``#`` inside the data is untouched.
    """
    n = 0
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            if line.lstrip().startswith("#"):
                n += 1
            else:
                break
    return n


def read_csv(path: str | Path) -> pd.DataFrame:
    d = pd.read_csv(path, sep=";", encoding="utf-8-sig", dtype=str,
                    keep_default_na=False, engine="python",
                    skiprows=_leading_comment_lines(path))
    return d.loc[:, [c for c in d.columns if not c.startswith("Unnamed")]]


def load(pat_path: str | Path, var_path: str | Path,
         design: StudyDesign | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the two analysis tables and label each row with its study period.

    ``design`` defaults to the one in the study-parameters file. A year outside
    both configured periods raises rather than being assigned silently.
    """
    if design is None:
        from mefv.config import study_design
        design = study_design()
    pat = read_csv(pat_path)
    var = read_csv(var_path)
    pat = pat[pat.patient_id.str.strip() != ""].copy()
    pat["year"] = pat.year_requested.astype(int)
    var["year"] = var.year_requested.astype(int)
    pat["n_variants"] = pat.n_variants.astype(int)
    for d in (pat, var):
        d["era"] = [design.period_of(y) for y in d.year]
    return pat, var


DOSE = {"Het": 1, "Hom": 2}


def genotype_vector(
    var: pd.DataFrame, variant: str, patient_ids: Sequence[str]
) -> list[int]:
    """Copies of the variant allele per patient, in the order of patient_ids.

    Returns 0, 1 or 2 for every identifier in ``patient_ids``. A patient with
    no row for this variant is encoded as 0.

    An unrecognised zygosity label is counted as one copy, which reproduces the
    original behaviour of the analysis, but a warning is now raised rather than
    the value being substituted silently.
    """
    sub = var[var.variant == variant]
    unknown = set(sub.zygosity) - set(DOSE)
    if unknown:
        warnings.warn(
            f"{variant}: unrecognised zygosity label(s) {sorted(unknown)} "
            "counted as one allele copy",
            stacklevel=2,
        )
    d = dict(zip(sub.patient_id, sub.zygosity.map(DOSE).fillna(1)))
    return [int(d.get(p, 0)) for p in patient_ids]


