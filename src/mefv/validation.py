"""Structured validation of the raw input tables.

Every check has a stable rule identifier (``P2``, ``V5``, ``X8`` ...), a
severity and a one-line statement of the rule. The checks never repair,
normalise, merge or drop anything: they count what they find and report it.

Severities
    fatal     structural defect; analysis cannot continue and the finding
              cannot be waived;
    blocking  ambiguous integrity problem (exact duplicates, conflicting
              genotypes, labels differing by case or whitespace); final
              outputs are withheld until the data owner has recorded an
              explicit disposition for each finding (see "Owner disposition").
              Acknowledging a rule by identifier does NOT release the run;
    warning   reported, does not stop the run;
    info      reported for completeness.

Owner disposition
    A blocking finding is a scientific conflict that the software cannot
    resolve. Naming a rule on the command line ("acknowledging" it) only marks
    it as seen: the status becomes ``acknowledged_unresolved`` and no final
    output is written. The run is released only by a disposition file, a
    private JSON record kept outside the repository, with one entry per
    blocking rule: the rule, the number of records reviewed (an integer of 5 or more equal to
    the current count, or "<5" for a finding with fewer than 5 records; a changed
    export invalidates it), the decision, who
    decided, when, and why. The only decision supported is
    ``use_as_recorded``: every row is counted exactly as it is in the export,
    nothing is removed, merged or selected. The report records the rule, the
    count, the decision and the SHA-256 of the file, never the names or the
    free text.

Privacy
    A report holds rule identifiers and counts only. It never contains an
    identifier, a variant label, a region name or any other cell value, and
    neither do the exception messages built from it. Invalid categorical values
    are described by their *shape* (blank, case/whitespace variant of a valid
    label, unrecognised), never echoed.

Scope
    This module checks structure and internal consistency. It cannot decide
    which of two conflicting records is correct, whether the declared
    ``n_variants`` or the observed rows are authoritative, or whether an
    individual who has no variant row was tested for it. Those are study
    decisions; see docs/validation.md.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import re
from collections import Counter
from collections.abc import Collection, Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import pandas as pd

from mefv.config import StudyDesign
from mefv.disclosure import SMALL_COUNT_LIMIT, SUPPRESSED, is_small, shown
from mefv.disclosure import text as _n

FATAL, BLOCKING, WARNING, INFO = "fatal", "blocking", "warning", "info"
SEVERITIES = (FATAL, BLOCKING, WARNING, INFO)

VALID_ZYGOSITY = ("Het", "Hom")
VALID_SEX = ("male", "female")
VALID_DETECTED = ("yes", "no")
_REGION_TOKEN = re.compile(r"^[A-Za-z0-9_.\-]+$")
_UINT = re.compile(r"^[0-9]+$")

REQUIRED_COLUMNS = {
    "method": ("Patient_ID", "YEAR", "sex", "Method"),
    "variants": ("patient_id", "variant", "exon", "zygosity", "year_requested"),
    "patients": (
        "patient_id",
        "year_requested",
        "sex",
        "n_variants",
        "variant_detected",
        "exons_tested",
        "coverage_source",
    ),
}


@dataclass(frozen=True)
class Rule:
    id: str
    severity: str
    table: str
    statement: str


def _catalogue() -> dict[str, Rule]:
    rows = [
        # ---- schema
        ("S1", FATAL, "any", "every required column is present"),
        # ---- patients table
        ("P1a", INFO, "patients", "rows that are entirely empty are removed and counted"),
        ("P1b", FATAL, "patients", "a row with content has a patient_id"),
        (
            "P2",
            FATAL,
            "patients",
            "patient_id is in canonical form (configured prefix + zero-padded number)",
        ),
        (
            "P3",
            FATAL,
            "patients",
            "distinct raw identifiers do not collapse to one canonical identifier",
        ),
        ("P4a", FATAL, "patients", "no individual has identical duplicate rows"),
        ("P4b", FATAL, "patients", "no individual has more than one distinct row"),
        ("P5a", FATAL, "patients", "year_requested is a non-negative integer"),
        ("P5b", FATAL, "patients", "year_requested lies in a configured study period"),
        ("P6", WARNING, "patients", "sex is one of the permitted labels"),
        ("P7", WARNING, "patients", "variant_detected is one of the permitted labels"),
        ("P8", FATAL, "patients", "n_variants is a non-negative integer"),
        # ---- exons_tested
        ("E1", WARNING, "patients", "exons_tested is not blank"),
        (
            "E2",
            WARNING,
            "patients",
            "exons_tested is comma-separated, with no whitespace or empty token",
        ),
        ("E3", WARNING, "patients", "each exons_tested token uses only the permitted characters"),
        ("E4a", WARNING, "patients", "each exons_tested token is a known region"),
        ("E5", WARNING, "patients", "exons_tested contains every configured panel region"),
        ("E6", INFO, "patients", "number of distinct exons_tested region sets"),
        # ---- variants table
        ("V1a", INFO, "variants", "rows that are entirely empty are removed and counted"),
        ("V1b", FATAL, "variants", "a row with content has a patient_id"),
        ("V10", FATAL, "variants", "patient_id is in canonical form"),
        ("V2", FATAL, "variants", "variant label is not blank"),
        ("V3", FATAL, "variants", "exon/region is not blank"),
        ("V4", FATAL, "variants", "zygosity is exactly 'Het' or 'Hom'"),
        ("V5", BLOCKING, "variants", "no row is an exact duplicate of another"),
        (
            "V6a",
            BLOCKING,
            "variants",
            "an individual has one row per variant (same genotype repeated)",
        ),
        (
            "V6b",
            BLOCKING,
            "variants",
            "an individual has one genotype per variant (conflicting genotypes)",
        ),
        (
            "V7",
            BLOCKING,
            "variants",
            "one spelling per variant label (case/whitespace); never merged automatically",
        ),
        ("V8", INFO, "variants", "rows touched by documented nomenclature corrections"),
        ("V9", WARNING, "variants", "a variant is reported in one region only"),
        ("V11", FATAL, "variants", "year_requested is a non-negative integer"),
        ("V12", FATAL, "variants", "year_requested lies in a configured study period"),
        # ---- method table
        ("M1", FATAL, "method", "Patient_ID can be parsed to the configured identifier scheme"),
        ("M2", FATAL, "method", "YEAR is a non-negative integer"),
        ("M3", INFO, "method", "rows that are entirely empty are removed and counted"),
        # ---- cross-table
        (
            "X1",
            FATAL,
            "cross",
            "every variant row belongs to an individual in the patients or method table",
        ),
        ("X2", INFO, "cross", "variant rows of individuals present only in the method table"),
        ("X2b", INFO, "cross", "individuals present only in the method table (restored by design)"),
        ("X3", FATAL, "cross", "every patient has a method record"),
        ("X4", FATAL, "cross", "one method row per individual"),
        ("X5a", WARNING, "cross", "declared n_variants does not exceed the observed variant rows"),
        ("X5b", WARNING, "cross", "declared n_variants is not below the observed variant rows"),
        (
            "X5c",
            INFO,
            "cross",
            "declared n_variants equals the distinct individual+variant pairs but not the rows",
        ),
        ("X6", WARNING, "cross", "variant_detected agrees with the declared n_variants"),
        ("X7", WARNING, "cross", "variant_detected agrees with the observed variant rows"),
        ("X8", WARNING, "cross", "a variant's region is listed in the individual's exons_tested"),
        (
            "X9",
            WARNING,
            "cross",
            "individuals with variant rows have an interpretable exons_tested",
        ),
        ("X10", WARNING, "cross", "year agrees between the patients and method tables"),
        ("X11", WARNING, "cross", "sex agrees between the patients and method tables"),
        # ---- built tables
        ("A1", FATAL, "built", "every individual has a resolved assay"),
        ("A2", FATAL, "built", "every assay label is one of the configured assays"),
        ("C1", FATAL, "built", "the cohort has the configured number of individuals"),
    ]
    return {r[0]: Rule(*r) for r in rows}


RULES: dict[str, Rule] = _catalogue()
#: Version of the report layout and of the rule semantics. Raise it whenever a rule is
#: changed in a way its statement does not show.
RULESET_VERSION = 2


def ruleset_sha256() -> str:
    """Digest of the rule catalogue (identifiers, severities, statements) and its version.

    An owner disposition is bound to it: if a rule or a severity changes, a disposition
    recorded under the old rules no longer applies.
    """
    h = hashlib.sha256(f"ruleset-version={RULESET_VERSION}\n".encode())
    for r in RULES.values():
        h.update(f"{r.id}\0{r.severity}\0{r.table}\0{r.statement}\n".encode())
    return h.hexdigest()


#: Rules whose findings stop the run until the owner records a disposition. They
#: may be named (acknowledged) on the command line; that alone releases nothing.
ACKNOWLEDGEABLE = tuple(r.id for r in RULES.values() if r.severity == BLOCKING)
#: Dispositions the software can honour. It never repairs, merges or deletes.
DISPOSITIONS = ("use_as_recorded",)
_DISPOSITION_FIELDS = ("rule", "records", "decision", "decided_by", "decided_on", "basis")
_ISO_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


class DispositionBinding(NamedTuple):
    """What an owner disposition file was made for. A change in any part voids it."""

    input_sha256: str
    ruleset_sha256: str
    config_sha256: str


@dataclass
class Finding:
    rule_id: str
    records: int = 0
    individuals: int = 0
    detail: dict[str, int] = field(default_factory=dict)

    @property
    def rule(self) -> Rule:
        return RULES[self.rule_id]

    def to_dict(self) -> dict:
        return {
            "rule": self.rule_id,
            "severity": self.rule.severity,
            "table": self.rule.table,
            "statement": self.rule.statement,
            "records": shown(self.records),
            "individuals": shown(self.individuals),
            "detail": {k: shown(v) for k, v in sorted(self.detail.items())},
        }


class DataIntegrityError(ValueError):
    """Input data failed validation. The message contains rule IDs and counts only."""

    def __init__(self, message: str, report: ValidationReport | None = None):
        super().__init__(message)
        self.report = report


@dataclass
class ValidationReport:
    tables: dict[str, int] = field(default_factory=dict)
    findings: list[Finding] = field(default_factory=list)
    acknowledged: tuple[str, ...] = ()
    #: (rule, records, decision) entries of a loaded owner disposition
    dispositions: tuple[tuple[str, int, str], ...] = ()
    disposition_sha256: str = ""
    #: SHA-256 of the input files this run read (set by the caller), per file and combined
    input_digests: dict[str, str] = field(default_factory=dict)
    input_sha256: str = ""
    #: the input digest the owner's disposition was made for
    disposition_input_sha256: str = ""
    #: the rule-catalogue digest the owner's disposition was made for
    disposition_ruleset_sha256: str = ""
    #: digest of the effective validation settings (set by the caller)
    config_sha256: str = ""
    disposition_config_sha256: str = ""
    #: a diagnostic-only run produces no analysis tables, so it never releases any
    diagnostic_only: bool = False

    def add(
        self,
        rule_id: str,
        records: int = 0,
        individuals: int = 0,
        detail: Mapping[str, int] | None = None,
        always: bool = False,
    ) -> None:
        """Record a finding. Rules with nothing to report are omitted unless ``always``."""
        detail = {k: int(v) for k, v in (detail or {}).items() if int(v)}
        if records or individuals or detail or always:
            self.findings.append(Finding(rule_id, int(records), int(individuals), detail))

    def by_severity(self, severity: str) -> list[Finding]:
        return [f for f in self.findings if f.rule.severity == severity]

    @property
    def fatal(self) -> list[Finding]:
        return [f for f in self.by_severity(FATAL) if f.records or f.individuals]

    def _covered(self, f: Finding) -> bool:
        """An owner disposition exists for this rule and this record count, and was made for
        exactly the input files, settings and rules in force now."""
        return self._bound() and any(
            r == f.rule_id and (n == f.records or (n == SUPPRESSED and is_small(f.records)))
            for r, n, _ in self.dispositions
        )

    def _input_matches(self) -> bool:
        return bool(self.input_sha256) and self.input_sha256 == self.disposition_input_sha256

    def _bound(self) -> bool:
        """The disposition was made for these input files, these settings and these rules."""
        return (
            self._input_matches()
            and bool(self.config_sha256)
            and self.config_sha256 == self.disposition_config_sha256
            and self.disposition_ruleset_sha256 == ruleset_sha256()
        )

    @property
    def blocking(self) -> list[Finding]:
        """Blocking findings with no matching owner disposition (acknowledgement is irrelevant)."""
        return [
            f
            for f in self.by_severity(BLOCKING)
            if (f.records or f.individuals) and not self._covered(f)
        ]

    @property
    def acknowledged_unresolved(self) -> list[str]:
        """Rules that were acknowledged but still have an undisposed blocking finding.

        Acknowledging a rule marks it as seen. It does not resolve the data
        problem and does not release the run.
        """
        return sorted({f.rule_id for f in self.blocking if f.rule_id in self.acknowledged})

    @property
    def dispositioned(self) -> list[Finding]:
        return [
            f
            for f in self.by_severity(BLOCKING)
            if (f.records or f.individuals) and self._covered(f)
        ]

    @property
    def releases_outputs(self) -> bool:
        return not self.fatal and not self.blocking

    @property
    def status(self) -> str:
        if self.fatal:
            return "failed"
        if self.blocking:
            only_ack = all(f.rule_id in self.acknowledged for f in self.blocking)
            return "acknowledged_unresolved" if only_ack else "blocked"
        if self.dispositioned:
            return "owner_disposition_recorded"
        return (
            "passed_with_findings"
            if any(
                f.rule.severity == WARNING and (f.records or f.individuals) for f in self.findings
            )
            else "passed"
        )

    def to_dict(self) -> dict:
        order = {rid: i for i, rid in enumerate(RULES)}
        return {
            "report_version": RULESET_VERSION,
            "ruleset_sha256": ruleset_sha256(),
            "small_counts": f"counts below 5 are shown as '{SUPPRESSED}'",
            "status": self.status,
            "tables": {k: shown(v) for k, v in sorted(self.tables.items())},
            "acknowledged_rules": sorted(self.acknowledged),
            "acknowledged_unresolved_rules": self.acknowledged_unresolved,
            "owner_dispositions": [
                {"rule": r, "records": shown(n) if isinstance(n, int) else n, "decision": d}
                for r, n, d in sorted(self.dispositions, key=lambda t: (t[0], str(t[1]), t[2]))
            ],
            "owner_disposition_file_sha256": self.disposition_sha256,
            "inputs_sha256": dict(sorted(self.input_digests.items())),
            "input_sha256": self.input_sha256,
            "config_sha256": self.config_sha256,
            "outputs_released": self.releases_outputs and not self.diagnostic_only,
            "diagnostic_only": self.diagnostic_only,
            "findings": [
                f.to_dict() for f in sorted(self.findings, key=lambda f: order[f.rule_id])
            ],
            "note": "Counts only, small counts suppressed. No identifier, label or cell value is "
            "recorded.",
        }

    def render(self, *, only_findings: bool = True) -> str:
        lines = [f"validation status: {self.status}"]
        order = {rid: i for i, rid in enumerate(RULES)}
        for f in sorted(self.findings, key=lambda f: order[f.rule_id]):
            if only_findings and not (f.records or f.individuals or f.detail):
                continue
            det = "; ".join(f"{k}={_n(v)}" for k, v in sorted(f.detail.items()))
            ack = (" [owner disposition recorded]" if self._covered(f)
                   else " [acknowledged, NOT resolved]" if f.rule_id in self.acknowledged
                   else "")
            lines.append(
                f"  {f.rule_id:4} {f.rule.severity:8} records={_n(f.records):<6} "
                f"individuals={_n(f.individuals):<5} {f.rule.statement}"
                + (f" ({det})" if det else "")
                + ack
            )
        return "\n".join(lines)

    def enforce(
        self,
        accepted: Collection[str] = (),
        dispositions: Collection[Mapping[str, object]] | None = None,
        disposition_sha256: str = "",
        binding: DispositionBinding | str = "",
    ) -> None:
        """Raise unless there is no fatal finding and every blocking finding has a disposition.

        ``accepted`` names blocking rules the owner has looked at. It never
        releases the run. ``dispositions`` are the validated entries of the
        owner's disposition file (see :func:`load_dispositions`). An entry covers a
        finding only if it names the same rule and the same record count AND the
        file was made for the same input files, the same effective settings and the
        same rule catalogue (``binding``; see :class:`DispositionBinding`): a changed
        export, configuration or rule set never inherits an old disposition, even
        when the counts happen to agree.
        Fatal rules can be neither acknowledged nor disposed of.
        """
        bad = sorted(set(accepted) - set(ACKNOWLEDGEABLE))
        if bad:
            raise DataIntegrityError(
                f"cannot acknowledge rule(s) {bad}: only {list(ACKNOWLEDGEABLE)} can be "
                "acknowledged; other rules are structural or informational.",
                self,
            )
        self.acknowledged = tuple(sorted(accepted))
        self.dispositions = tuple(
            (str(d["rule"]), d["records"], str(d["decision"])) for d in (dispositions or ())
        )
        self.disposition_sha256 = disposition_sha256
        if isinstance(binding, DispositionBinding):
            self.disposition_input_sha256 = binding.input_sha256
            self.disposition_ruleset_sha256 = binding.ruleset_sha256
            self.disposition_config_sha256 = binding.config_sha256
        else:  # a bare input digest binds nothing else, so it can never cover a finding
            self.disposition_input_sha256 = binding
            self.disposition_ruleset_sha256 = self.disposition_config_sha256 = ""
        if self.fatal:
            raise DataIntegrityError(
                "input validation failed ("
                + ", ".join(f"{f.rule_id}: {_n(f.records)} record(s)" for f in self.fatal)
                + "). See the validation report; no analytical output was produced.",
                self,
            )
        if self.blocking:
            stale = sorted(
                {r for r, _, _ in self.dispositions}
                & {f.rule_id for f in self.blocking}
            )
            other_input = bool(self.dispositions) and not self._bound()
            raise DataIntegrityError(
                "final outputs are withheld: unresolved blocking finding(s) "
                + ", ".join(f"{f.rule_id}: {_n(f.records)} record(s)" for f in self.blocking)
                + (
                    "; acknowledged but not resolved: "
                    + ", ".join(self.acknowledged_unresolved)
                    if self.acknowledged_unresolved
                    else ""
                )
                + (
                    "; the disposition was made for different input files, settings or validation "
                    "rules (the export, the study parameters or the rule catalogue changed since "
                    "it was recorded)"
                    if other_input
                    else ""
                )
                + (
                    "; a disposition exists for "
                    + ", ".join(stale)
                    + " but for a different record count (the data changed)"
                    if stale and not other_input
                    else ""
                )
                + ". Acknowledging does not release the run. Correct the export, or record an "
                "explicit owner disposition (--owner-disposition FILE, see docs/validation.md).",
                self,
            )


def digest_inputs(paths: Mapping[str, str | Path]) -> tuple[dict[str, str], str]:
    """SHA-256 of each input file, and one combined digest over names and digests."""
    each = {k: hashlib.sha256(Path(p).read_bytes()).hexdigest() for k, p in sorted(paths.items())}
    h = hashlib.sha256()
    for k, d in each.items():
        h.update(f"{k}\0{d}\n".encode())
    return each, h.hexdigest()


def digest_config(settings: Mapping[str, object]) -> str:
    """Digest of the effective study settings, independent of formatting and key order.

    Two parameter files with the same values give the same digest; any changed value
    (identifier scheme, periods, panel regions, nomenclature corrections, expected cohort
    size ...) gives a different one. An owner disposition is bound to it.
    """
    canonical = json.dumps(settings, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _no_duplicate_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    keys = [k for k, _ in pairs]
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate key")
    return dict(pairs)


def _reject_constant(_name: str) -> None:
    raise ValueError("non-finite number")


def _is_calendar_date(value: object) -> bool:
    if not (isinstance(value, str) and _ISO_DATE.match(value)):
        return False
    try:
        datetime.date.fromisoformat(value)
    except ValueError:
        return False
    return True


_TOP_LEVEL = ("input_sha256", "ruleset_sha256", "config_sha256", "dispositions")


def load_dispositions(path: str) -> tuple[list[dict[str, object]], str, DispositionBinding]:
    """Read and validate an owner disposition file.

    Returns the entries, the SHA-256 of the file and the :class:`DispositionBinding`
    it carries. The file is JSON::

        {"input_sha256": "...", "ruleset_sha256": "...", "config_sha256": "...",
         "dispositions": [{"rule", "records", "decision", "decided_by", "decided_on",
                           "basis"}, ...]}

    The three digests are copied from the validation report of the run being
    reviewed. A repeated key, an unknown key, a non-finite number, a date that is not
    a real calendar date, or a malformed digest is refused. ``records`` is the count
    shown in the report: a positive integer, or the string ``"<5"`` for a finding
    whose count is suppressed (the input digest still pins the exact data). Errors
    name the entry number and the field, never a value.
    """
    raw = Path(path).read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    try:
        doc = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=_no_duplicate_keys,
            parse_constant=_reject_constant,
        )
    except UnicodeDecodeError:
        raise DataIntegrityError("the owner disposition file is not valid UTF-8 JSON") from None
    except ValueError as exc:
        why = " (a key is repeated)" if "duplicate key" in str(exc) else ""
        raise DataIntegrityError(f"the owner disposition file is not valid JSON{why}") from None
    if not isinstance(doc, dict):
        raise DataIntegrityError("the owner disposition file must contain a JSON object")
    if set(doc) - set(_TOP_LEVEL):
        raise DataIntegrityError("the owner disposition file has an unknown top-level key")
    digests = []
    for k in _TOP_LEVEL[:3]:
        v = doc.get(k)
        if not (isinstance(v, str) and re.fullmatch(r"[0-9a-f]{64}", v)):
            raise DataIntegrityError(
                f"the owner disposition file needs '{k}': a SHA-256 digest copied from the "
                "validation report of the run that was reviewed"
            )
        digests.append(v)
    entries = doc.get("dispositions")
    if not isinstance(entries, list) or not entries:
        raise DataIntegrityError(
            "the owner disposition file needs a non-empty list under 'dispositions'"
        )
    out: list[dict[str, object]] = []
    seen: set[str] = set()
    for n, e in enumerate(entries, 1):
        if not isinstance(e, dict):
            raise DataIntegrityError(f"owner disposition entry {n}: not an object")
        missing = [k for k in _DISPOSITION_FIELDS if k not in e]
        if missing:
            raise DataIntegrityError(f"owner disposition entry {n}: missing field(s) {missing}")
        if set(e) - set(_DISPOSITION_FIELDS):
            raise DataIntegrityError(f"owner disposition entry {n}: unknown field")
        if e["rule"] not in ACKNOWLEDGEABLE:
            raise DataIntegrityError(
                f"owner disposition entry {n}: 'rule' must be one of {list(ACKNOWLEDGEABLE)}"
            )
        if e["rule"] in seen:
            raise DataIntegrityError(f"owner disposition entry {n}: rule listed more than once")
        seen.add(e["rule"])
        rec = e["records"]
        exact = isinstance(rec, int) and not isinstance(rec, bool) and rec >= SMALL_COUNT_LIMIT
        if not (exact or rec == SUPPRESSED):
            raise DataIntegrityError(
                f"owner disposition entry {n}: 'records' must be the count shown in the report: "
                f"an integer of {SMALL_COUNT_LIMIT} or more, or the text '{SUPPRESSED}'"
            )
        if e["decision"] not in DISPOSITIONS:
            raise DataIntegrityError(
                f"owner disposition entry {n}: 'decision' must be one of {list(DISPOSITIONS)}"
            )
        for k in ("decided_by", "basis"):
            if not isinstance(e[k], str) or not e[k].strip():
                raise DataIntegrityError(
                    f"owner disposition entry {n}: '{k}' must be non-empty text"
                )
        if not _is_calendar_date(e["decided_on"]):
            raise DataIntegrityError(
                f"owner disposition entry {n}: 'decided_on' must be a real date, YYYY-MM-DD"
            )
        out.append(e)
    return out, digest, DispositionBinding(*digests)


# ---------------------------------------------------------------- helpers ---
def canonical_id(value: object, prefix: str, width: int) -> str | None:
    """PREFIX + zero-padded number, or None when the value is not of that form.

    Accepted: optional surrounding whitespace, the prefix in any letter case (at most once,
    at the start), then ASCII digits only. Anything else is None, never repaired: an
    embedded or repeated prefix, a separator, a sign or a missing value (None, NaN, NA)
    does not become a number, so two different raw identifiers cannot be forced onto the
    same canonical one by stripping characters out of the middle. Numbers longer than
    ``width`` are kept whole, never truncated.
    """
    if _is_missing(value):
        return None
    s = str(value).strip().lower()
    if prefix and s.startswith(prefix.lower()):
        s = s[len(prefix):]
    return f"{prefix}{int(s):0{width}d}" if (s.isascii() and s.isdigit()) else None


def _is_missing(value: object) -> bool:
    """None, NaN or pd.NA (a scalar only; a container is never 'missing')."""
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def shape_of(value: str, valid: Iterable[str]) -> str:
    """Classify a categorical value without echoing it."""
    valid = tuple(valid)
    if value in valid:
        return "valid"
    if value.strip() == "":
        return "blank"
    if value.strip().casefold() in {v.casefold() for v in valid}:
        return "case_or_whitespace_variant_of_valid"
    return "unrecognised"


def is_uint(series: pd.Series) -> pd.Series:
    """True where the (stripped) text is a non-negative integer."""
    return series.astype(str).str.strip().str.fullmatch(_UINT)


def parse_uint_column(series: pd.Series, what: str, allow_blank: bool = False) -> pd.Series:
    """Integer column from text; raises a value-free error if any entry is not an integer.

    With ``allow_blank`` an empty entry is read as missing (nullable ``Int64``
    result). It is never read as zero.
    """
    text = series.astype(str).str.strip()
    ok = is_uint(series)
    blank = text == ""
    if allow_blank:
        ok = ok | blank
    if not ok.all():
        raise DataIntegrityError(
            f"{what}: {_n(int((~ok).sum()))} entr(ies) are not non-negative integers."
        )
    if allow_blank:
        return pd.to_numeric(text.where(~blank)).astype("Int64")
    return text.astype(int)


def parse_region_set(value: str) -> frozenset[str] | None:
    """The set of regions in an ``exons_tested`` value, or None if it is unusable.

    None is returned for a blank or malformed value (wrong separator,
    whitespace, empty token, unsupported characters). It is never read as "all
    regions tested": eligibility for such an individual is undetermined.
    """
    if _is_missing(value) or str(value).strip() == "":
        return None
    tokens = str(value).split(",")
    if any(t == "" or not _REGION_TOKEN.match(t) for t in tokens):
        return None
    return frozenset(tokens)


def _blank(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip() == ""


def _empty_rows(df: pd.DataFrame) -> pd.Series:
    return df.apply(_blank).all(axis=1) if len(df.columns) else pd.Series(False, index=df.index)


def missing_columns(
    method: pd.DataFrame, variants: pd.DataFrame, patients: pd.DataFrame
) -> dict[str, list[str]]:
    tables = {"method": method, "variants": variants, "patients": patients}
    return {
        t: [c for c in cols if c not in tables[t].columns]
        for t, cols in REQUIRED_COLUMNS.items()
        if any(c not in tables[t].columns for c in cols)
    }


# -------------------------------------------------------------- validation --
def validate_raw_inputs(
    method: pd.DataFrame,
    variants: pd.DataFrame,
    patients: pd.DataFrame,
    *,
    design: StudyDesign,
    id_prefix: str,
    id_width: int,
    nomenclature_fixes: Mapping[str, str] | None = None,
) -> ValidationReport:
    """Check the three raw tables. Nothing is changed; findings are counts only."""
    rep = ValidationReport(
        tables={
            "method_rows": len(method),
            "variants_rows": len(variants),
            "patients_rows": len(patients),
        }
    )
    gaps = missing_columns(method, variants, patients)
    if gaps:
        rep.add(
            "S1",
            records=sum(len(v) for v in gaps.values()),
            detail={f"{t}_missing_columns": len(v) for t, v in gaps.items()},
        )
        return rep

    fixes = dict(nomenclature_fixes or {})
    periods = (design.primary, design.sensitivity)
    canon = lambda x: canonical_id(x, id_prefix, id_width)  # noqa: E731

    def in_period(years: pd.Series) -> pd.Series:
        return years.map(lambda y: any(p.contains(int(y)) for p in periods)).astype(bool)

    # ----------------------------------------------------------- patients --
    p = patients.astype(str)
    empty = _empty_rows(p) & _blank(p.patient_id)
    rep.add("P1a", records=int(empty.sum()))
    rep.add("P1b", records=int((_blank(p.patient_id) & ~empty).sum()))
    p = p[~_blank(p.patient_id)]

    cid = p.patient_id.map(canon)
    rep.add("P2", records=int((cid != p.patient_id).sum()))
    ok = cid.notna()
    raw_per_canon = p[ok].groupby(cid[ok]).patient_id.nunique()
    clash = raw_per_canon.index[raw_per_canon > 1]
    rep.add("P3", records=int(cid.isin(clash).sum()), individuals=len(clash))

    size = p.groupby("patient_id").size()
    distinct = p.drop_duplicates().groupby("patient_id").size()
    exact_ids = size.index[size > distinct.reindex(size.index)]
    conflict_ids = distinct.index[distinct > 1]
    rep.add("P4a", records=int(p.patient_id.isin(exact_ids).sum()), individuals=len(exact_ids))
    rep.add(
        "P4b", records=int(p.patient_id.isin(conflict_ids).sum()), individuals=len(conflict_ids)
    )

    yr_ok = is_uint(p.year_requested)
    rep.add("P5a", records=int((~yr_ok).sum()))
    rep.add("P5b", records=int((~in_period(p.year_requested[yr_ok].str.strip().astype(int))).sum()))
    for rid, col, valid in (("P6", "sex", VALID_SEX), ("P7", "variant_detected", VALID_DETECTED)):
        c = Counter(shape_of(v, valid) for v in p[col])
        c.pop("valid", None)
        rep.add(rid, records=sum(c.values()), detail=c)
    nv_ok = is_uint(p.n_variants)
    rep.add("P8", records=int((~nv_ok).sum()))

    # ------------------------------------------------------- exons_tested --
    # Known regions: the configured panel plus every region observed in the
    # variant table. No closed vocabulary of extended-only regions is
    # configured, so "known" cannot be stricter than this.
    v_exon = variants.exon.astype(str)
    known = set(design.panel_regions) | set(v_exon[~_blank(v_exon)])
    sets: Counter = Counter()
    n_blank = n_sep = n_chr = n_unknown = n_nopanel = 0
    for s in p.exons_tested:
        if s.strip() == "":
            n_blank += 1
            continue
        toks = s.split(",")
        sep_bad = bool(re.search(r"[;|/\s]", s)) or any(t == "" for t in toks)
        chr_bad = (not sep_bad) and any(not _REGION_TOKEN.match(t) for t in toks)
        n_sep += sep_bad
        n_chr += chr_bad
        if not sep_bad and not chr_bad:
            n_unknown += any(t not in known for t in toks)
            n_nopanel += not set(design.panel_regions) <= set(toks)
            sets[frozenset(toks)] += 1
    rep.add("E1", records=n_blank)
    rep.add("E2", records=n_sep)
    rep.add("E3", records=n_chr)
    rep.add("E4a", records=n_unknown)
    rep.add("E5", records=n_nopanel)
    e6 = {"distinct_region_sets": len(sets)}
    for n_regions, n in Counter(len(k) for k in sets.elements()).items():
        e6[f"individuals_with_{n_regions}_regions"] = n
    rep.add("E6", detail=e6)

    # ----------------------------------------------------------- variants --
    v = variants.astype(str)
    v_empty = _empty_rows(v) & _blank(v.patient_id)
    rep.add("V1a", records=int(v_empty.sum()))
    rep.add("V1b", records=int((_blank(v.patient_id) & ~v_empty).sum()))
    v = v[~_blank(v.patient_id)].copy()
    rep.add("V10", records=int((v.patient_id.map(canon) != v.patient_id).sum()))
    rep.add("V2", records=int(_blank(v.variant).sum()))
    rep.add("V3", records=int(_blank(v.exon).sum()))
    z = Counter(shape_of(x, VALID_ZYGOSITY) for x in v.zygosity)
    z.pop("valid", None)
    rep.add("V4", records=sum(z.values()), detail=z)
    vy_ok = is_uint(v.year_requested)
    rep.add("V11", records=int((~vy_ok).sum()))
    rep.add("V12", records=int((~in_period(v.year_requested[vy_ok].str.strip().astype(int))).sum()))

    dup_rows = v.duplicated(keep=False)
    rep.add(
        "V5",
        records=int(dup_rows.sum()),
        individuals=int(v[dup_rows].patient_id.nunique()),
        detail={"surplus_rows": int(v.duplicated().sum())},
    )
    v["variant_fixed"] = v.variant.replace(fixes)
    rep.add("V8", records=int(v.variant.isin(fixes).sum()))
    # Genotype multiplicity is judged on the labels the pipeline will count,
    # i.e. after the documented nomenclature corrections.
    pair = v.groupby(["patient_id", "variant_fixed"]).agg(
        n=("zygosity", "size"), nz=("zygosity", "nunique")
    )
    multi = pair[pair.n > 1].reset_index()
    same, conf = multi[multi.nz == 1], multi[multi.nz > 1]
    rep.add(
        "V6a",
        records=int(same.n.sum()),
        individuals=same.patient_id.nunique(),
        detail={"individual_variant_pairs": len(same)},
    )
    rep.add(
        "V6b",
        records=int(conf.n.sum()),
        individuals=conf.patient_id.nunique(),
        detail={"individual_variant_pairs": len(conf)},
    )
    labels = v.variant[~_blank(v.variant)]
    folded = labels.str.strip().str.casefold()
    spell = labels.groupby(folded).nunique()
    split = spell.index[spell > 1]
    rep.add("V7", records=int(folded.isin(split).sum()), detail={"label_groups": len(split)})
    with_region = v[~_blank(v.exon)]
    regions_per_variant = with_region.groupby("variant_fixed").exon.nunique()
    multi_region = regions_per_variant.index[regions_per_variant > 1]
    rep.add(
        "V9",
        records=int(with_region.variant_fixed.isin(multi_region).sum()),
        detail={"variant_labels": len(multi_region)},
    )

    # -------------------------------------------------------------- method --
    m = method.astype(str)
    m_empty = _empty_rows(m) & _blank(m.Patient_ID)
    rep.add("M3", records=int(m_empty.sum()))
    m = m[~m_empty].copy()
    m["cid"] = m.Patient_ID.map(canon)
    rep.add("M1", records=int(m.cid.isna().sum()))
    rep.add("M2", records=int((~is_uint(m.YEAR)).sum()))
    mm = m[m.cid.notna()]
    dup_m = mm.cid.duplicated(keep=False)
    rep.add("X4", records=int(dup_m.sum()), individuals=int(mm[dup_m].cid.nunique()))

    # --------------------------------------------------------- cross-table --
    pat_ids, meth_ids = set(p.patient_id), set(mm.cid)
    only_meth = meth_ids - pat_ids
    vs = v[(~_blank(v.patient_id)) & (v.patient_id.map(canon) == v.patient_id)]
    orphan = vs[~vs.patient_id.isin(pat_ids | meth_ids)]
    rep.add("X1", records=len(orphan), individuals=orphan.patient_id.nunique())
    from_meth = vs[vs.patient_id.isin(only_meth)]
    rep.add("X2", records=len(from_meth), individuals=from_meth.patient_id.nunique())
    rep.add("X2b", records=len(only_meth), individuals=len(only_meth))
    no_method = pat_ids - meth_ids
    rep.add("X3", records=len(no_method), individuals=len(no_method))

    # Declared vs observed, for individuals with a single, parseable patient row.
    p1 = p[nv_ok & ~p.patient_id.isin(set(exact_ids) | set(conflict_ids))].drop_duplicates(
        "patient_id"
    )
    p1 = p1.set_index("patient_id")
    declared = p1.n_variants.str.strip().astype(int)
    rows_by = v.groupby("patient_id").size()
    pairs_by = v.drop_duplicates(["patient_id", "variant_fixed"]).groupby("patient_id").size()
    observed = declared.index.to_series().map(rows_by).fillna(0).astype(int)
    pairs = declared.index.to_series().map(pairs_by).fillna(0).astype(int)
    above, below = declared > observed, declared < observed
    rep.add(
        "X5a",
        records=int(above.sum()),
        individuals=int(above.sum()),
        detail={
            "declared_total": int(declared[above].sum()),
            "observed_rows_total": int(observed[above].sum()),
        },
    )
    rep.add(
        "X5b",
        records=int(below.sum()),
        individuals=int(below.sum()),
        detail={
            "declared_total": int(declared[below].sum()),
            "observed_rows_total": int(observed[below].sum()),
        },
    )
    rep.add("X5c", records=int(((declared == pairs) & (declared != observed)).sum()))
    flag = p1.variant_detected
    rep.add(
        "X6",
        records=int(
            (((flag == "yes") & (declared == 0)) | ((flag == "no") & (declared > 0))).sum()
        ),
    )
    rep.add(
        "X7",
        records=int(
            (((flag == "yes") & (observed == 0)) | ((flag == "no") & (observed > 0))).sum()
        ),
    )

    regions = p1.exons_tested.map(parse_region_set)
    vv = v[v.patient_id.isin(p1.index) & ~_blank(v.exon)]
    tested = vv.patient_id.map(regions)
    unusable = tested.isna()
    contra = ~unusable & pd.Series(
        [e not in t if t is not None else False for e, t in zip(vv.exon, tested)], index=vv.index
    )
    rep.add("X8", records=int(contra.sum()), individuals=vv[contra].patient_id.nunique())
    rep.add("X9", records=int(unusable.sum()), individuals=vv[unusable].patient_id.nunique())

    both = p1.index.intersection(mm.drop_duplicates("cid").set_index("cid").index)
    ms = mm.drop_duplicates("cid").set_index("cid")
    yd = p1.loc[both, "year_requested"].str.strip() != ms.loc[both, "YEAR"].str.strip()
    sd = p1.loc[both, "sex"].str.strip() != ms.loc[both, "sex"].str.strip()
    rep.add("X10", records=int(yd.sum()), individuals=int(yd.sum()))
    rep.add("X11", records=int(sd.sum()), individuals=int(sd.sum()))
    return rep


def validate_built_tables(
    pat: pd.DataFrame, assays: Collection[str], expected_cohort_size: int, report: ValidationReport
) -> ValidationReport:
    """Post-assembly checks that replace the former ``assert`` statements."""
    a = pat.assay.astype(object)
    unresolved = a.isna() | (a.astype(str).str.strip().isin(("", "not_recorded")))
    report.add("A1", records=int(unresolved.sum()))
    report.add("A2", records=int((~unresolved & ~a.isin(set(assays))).sum()))
    if len(pat) != expected_cohort_size:
        report.add("C1", records=1, detail={"expected": expected_cohort_size, "found": len(pat)})
    return report
