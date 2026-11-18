"""Small-count suppression for diagnostic output.

A count of one to four individuals or records can single out a person, so
diagnostic output (validation reports, scanner reports, error messages and the
run log) never states one. A count below :data:`SMALL_COUNT_LIMIT` is shown as
``"<5"``. Zero is shown as zero: it identifies nobody.

Only what a reader *sees* is suppressed. The software itself keeps exact counts
to make its decisions (for example whether a stored owner disposition matches
the data), and the study's analytical result tables are not touched by this
module; whether they need cell suppression before sharing is a disclosure
decision for the data owner (docs/owner_decisions.md).
"""
from __future__ import annotations

import numbers

import pandas as pd

SMALL_COUNT_LIMIT = 5
SUPPRESSED = f"<{SMALL_COUNT_LIMIT}"


def is_small(n: object) -> bool:
    """True for a count of 1 to ``SMALL_COUNT_LIMIT - 1``."""
    try:
        return 0 < int(n) < SMALL_COUNT_LIMIT
    except (TypeError, ValueError):
        return False


def shown(n: object) -> int | str:
    """The count as it may be displayed: the integer, or ``"<5"``."""
    return SUPPRESSED if is_small(n) else int(n)


def text(n: object) -> str:
    """:func:`shown` as text, for messages."""
    return str(shown(n))


def _cell(x: object) -> object:
    if isinstance(x, bool) or not isinstance(x, numbers.Real) or pd.isna(x):
        return x
    return shown(x) if float(x) == int(x) else x


def suppress_frame(df: pd.DataFrame) -> pd.DataFrame:
    """A copy of a table of counts with every small whole number replaced by ``"<5"``."""
    return df.astype(object).map(_cell)


def suppressed_counts(series: pd.Series) -> pd.Series:
    """A ``value_counts()`` result with small counts replaced by ``"<5"``."""
    return series.astype(object).map(shown)
