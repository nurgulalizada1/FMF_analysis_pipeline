"""Population reference data: schema, provenance and validation.

Stage 06 compares cohort allele frequencies with a population reference table
transcribed from an external database. This module defines what that table
must look like and refuses one that does not, so that a malformed, incomplete
or undocumented reference file stops the pipeline instead of silently
changing a result.

The file is semicolon-separated and opens with ``#`` header lines. Three of
them are required, each as ``# key : value``:

* ``source``    - the database the values were transcribed from;
* ``version``   - the database release;
* ``retrieved`` - when the values were obtained.

Further ``# key : value`` lines (for example ``verified`` or ``open``) are kept
as provenance; a line beginning with ``#`` and white space continues the
previous value. If the study parameters name an expected version
(``reference_data.population_frequencies_version``), the header must state
exactly that version, so a different release cannot be substituted unnoticed.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from mefv.disclosure import text as _n

ANCESTRY_GROUPS = ("mid", "nfe", "fin", "asj", "sas", "eas", "afr", "amr", "rmi", "ami")

#: Required columns, in order.
COLUMNS = ("variant", *ANCESTRY_GROUPS, "glob", "grpmax", "hom", "xx", "xy",
           "revel", "cadd", "cv", "stars")

#: Columns holding allele frequencies, which must lie in [0, 1].
FREQUENCY_COLUMNS = (*ANCESTRY_GROUPS, "glob", "grpmax", "xx", "xy")

#: Columns that must be numeric when not empty.
NUMERIC_COLUMNS = (*FREQUENCY_COLUMNS, "hom", "revel", "cadd", "stars")

REQUIRED_PROVENANCE = ("source", "version", "retrieved")

_KEY = re.compile(r"^#\s([a-z_]+)\s*:\s*(.*)$")
_CONTINUATION = re.compile(r"^#\s{2,}(\S.*)$")


class ReferenceDataError(ValueError):
    """The reference table is missing, malformed or undocumented."""


def read_provenance(path: str | Path) -> dict[str, str]:
    """The ``# key : value`` provenance lines at the top of the file."""
    out: dict[str, str] = {}
    last = None
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.lstrip().startswith("#"):
                break
            m = _KEY.match(line)
            if m:
                last = m.group(1)
                out[last] = m.group(2).strip()
                continue
            c = _CONTINUATION.match(line)
            if c and last is not None:
                out[last] = f"{out[last]} {c.group(1).strip()}".strip()
            else:
                last = None
    return out


def load_population_frequencies(
    path: str | Path, expected_version: str | None = None
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Read and validate the population reference table.

    Returns the table, with numeric columns converted, and its provenance.
    Raises ``ReferenceDataError`` on any departure from the schema.
    """
    path = Path(path)
    if not path.exists():
        raise ReferenceDataError(
            f"population reference table not found at {path}. Copy "
            "data/external/population_frequencies.example.csv, fill it in and "
            "record its provenance, or set MEFV_POPULATION_FREQUENCIES."
        )
    provenance = read_provenance(path)
    missing = [k for k in REQUIRED_PROVENANCE if not provenance.get(k)]
    if missing:
        raise ReferenceDataError(
            f"{path}: provenance header lines missing or empty: {missing}. "
            "Each must be written as '# key : value' above the column names."
        )
    if expected_version is not None and provenance["version"] != expected_version:
        raise ReferenceDataError(
            f"{path}: header states version {provenance['version']!r}, but the study "
            f"parameters expect {expected_version!r}. A different reference release "
            "is not substituted silently."
        )

    n_header = 0
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            if line.lstrip().startswith("#"):
                n_header += 1
            else:
                break
    df = pd.read_csv(path, sep=";", encoding="utf-8-sig", dtype=str,
                     keep_default_na=False, skiprows=n_header)
    if tuple(df.columns) != COLUMNS:
        raise ReferenceDataError(
            f"{path}: columns must be {list(COLUMNS)}, found {list(df.columns)}")
    if (df.variant.str.strip() == "").any():
        raise ReferenceDataError(f"{path}: empty variant label")
    dup = sorted(df.variant[df.variant.duplicated()].unique())
    if dup:
        raise ReferenceDataError(f"{path}: {_n(len(dup))} duplicated variant label(s)")

    for c in NUMERIC_COLUMNS:
        raw = df[c].str.strip()
        num = pd.to_numeric(raw, errors="coerce")
        bad = raw.ne("") & num.isna()
        if bad.any():
            raise ReferenceDataError(
                f"{path}: column {c!r} has {_n(int(bad.sum()))} non-numeric value(s)")
        df[c] = num
    for c in FREQUENCY_COLUMNS:
        out_of_range = df[c].notna() & ~df[c].between(0, 1)
        if out_of_range.any():
            raise ReferenceDataError(
                f"{path}: column {c!r} has {_n(int(out_of_range.sum()))} frequency value(s) "
                "outside [0, 1]")
    return df, provenance
