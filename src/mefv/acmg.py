"""ACMG/AMP evidence-code configuration and its provenance.

The population-frequency thresholds and the list of codes treated as not
applicable were carried over from the original analysis. They are **provisional
expert-judgment rules**: this repository does not document who set them, from
which guideline version, or on what evidence. They are kept so that results do
not change silently, written beside every evidence package so that they cannot
be mistaken for validated rules, and listed as owner decisions. Nothing here
classifies a variant.
"""

from __future__ import annotations

import pandas as pd

from mefv.comparison import BA1_THRESHOLD, BS1_THRESHOLD

#: Codes of the ACMG/AMP evidence framework, in the order used on the forms.
CODES = (
    "PVS1",
    "PS1",
    "PS2",
    "PS3",
    "PS4",
    "PM1",
    "PM2",
    "PM3",
    "PM4",
    "PM5",
    "PM6",
    "PP1",
    "PP2",
    "PP3",
    "PP4",
    "PP5",
    "BA1",
    "BS1",
    "BS2",
    "BS3",
    "BS4",
    "BP1",
    "BP2",
    "BP3",
    "BP4",
    "BP5",
    "BP6",
    "BP7",
)

#: Codes pre-filled "n/a" on the assessor forms because the data this pipeline
#: holds cannot inform them (no segregation, de novo, phenotype or allelic data).
NOT_APPLICABLE = frozenset({"PS2", "PM3", "PM6", "PP1", "PP4", "PP5", "BS4", "BP2", "BP5", "BP6"})

STATUS = "provisional; carried over from the original analysis; provenance not documented"

RULE_PROVENANCE = [
    {
        "item": "BA1_threshold",
        "value": BA1_THRESHOLD,
        "rule": "BA1 flagged when reference grpmax filtering AF is strictly above this value",
        "evidence_status": STATUS,
        "owner_decision": "confirm or replace; state the guideline/VCEP source",
    },
    {
        "item": "BS1_threshold",
        "value": BS1_THRESHOLD,
        "rule": "BS1 flagged when reference grpmax filtering AF is strictly above this value",
        "evidence_status": STATUS,
        "owner_decision": "confirm or replace; BS1 thresholds are disease-specific",
    },
    {
        "item": "missing_grpmax",
        "value": "not_evaluable",
        "rule": "a variant without a reference grpmax is BA1/BS1 not_evaluable, never 'no'",
        "evidence_status": "conservative rule introduced during remediation",
        "owner_decision": "none required",
    },
    {
        "item": "not_applicable_codes",
        "value": " ".join(sorted(NOT_APPLICABLE)),
        "rule": "codes pre-filled 'n/a' on the assessor forms",
        "evidence_status": STATUS,
        "owner_decision": "confirm the list or leave all codes to the assessors",
    },
    {
        "item": "classification",
        "value": "none",
        "rule": "the pipeline assembles evidence and blank forms only; it assigns no ACMG code "
        "and no pathogenicity classification",
        "evidence_status": "by design",
        "owner_decision": "none required",
    },
]


def provenance_table() -> pd.DataFrame:
    return pd.DataFrame(RULE_PROVENANCE)
