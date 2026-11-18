"""Repository privacy scanner: files that should not be published.

Two layers, both standard-library only:

1. **General rules**, built in: sensitive file locations and types, identifier-
   like and contact-like content, tables that look like patient or variant
   exports, and credentials.
2. **A private denylist**, kept *outside* the repository, for the strings that
   only the data owner knows are sensitive (the real identifier scheme, names,
   laboratory file names). Entries may be literals, regular expressions or
   SHA-256 digests of tokens, so the list itself need not hold the plain text.
   The scanner refuses a denylist that lives inside the repository, and works
   normally (reporting that it ran without one) when none is available.

What a report contains
    file path, rule identifier and a count, nothing else. A matched value, a
    line or a snippet is never printed, logged or written to a report, and the
    denylist entries are reported only by their position (``entry 3``).

What this does not do
    Pattern scanning finds what its patterns describe. It cannot recognise
    arbitrary personal data, it does not read binary files, it does not inspect
    Git history, and a clean result is **not** a guarantee that nothing
    sensitive is present. It is a safety net behind careful review, not a
    substitute for it. See docs/privacy_policy.md.
"""
from __future__ import annotations

import fnmatch
import hashlib
import os
import re
import subprocess
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from mefv.disclosure import shown
from mefv.disclosure import text as _n

DENYLIST_ENV = "MEFV_PRIVATE_DENYLIST"
MAX_BYTES = 20_000_000
#: Paths (fnmatch, relative to the repository root) of documented placeholder examples that
#: are allowed to look like tables. Nothing else is exempt, whatever a file says about itself.
SYNTHETIC_EXAMPLES = ("data/external/*.example.csv", "data/reviewed/*.example.csv")

# ---------------------------------------------------------------- rules ------
#: rule id -> description
RULES: dict[str, str] = {
    "PATH-DATA": "file under data/ or results/ other than the documented placeholders",
    "PATH-TYPE": "file type that normally carries laboratory or patient data "
                 "(spreadsheet, archive, bundle, patch, database, image)",
    "PATH-NAME": "file name suggesting an export, patient file, credential or private list",
    "PATH-SYMLINK": "symbolic link (never followed; its target is not scanned, so it is reported)",
    "CONTENT-EMAIL": "e-mail address",
    "CONTENT-PHONE": "international telephone number",
    "CONTENT-LONGNUM": "long digit sequence that could be a national or medical record number",
    "CONTENT-DOB": "date-of-birth field with a value",
    "CONTENT-TABLE": "delimited table with a patient-identifier header and data rows",
    "CONTENT-GENOTYPE": "several rows with a Het/Hom genotype column (variant export)",
    "CONTENT-CREDENTIAL": "private key or access-token pattern",
    "DENYLIST": "match against the owner's private denylist",
}

_ALLOWED_DATA = ("data/README.md",)
_DATA_OK = re.compile(r"(^|/)(\.gitkeep|[^/]*\.example\.(csv|json))$")
_TYPE_EXT = (".xlsx", ".xls", ".xlsm", ".docx", ".doc", ".pdf", ".zip", ".tgz", ".gz", ".7z",
             ".tar", ".bundle", ".patch", ".sqlite", ".db", ".parquet", ".pkl", ".sav", ".dta",
             ".png", ".jpg", ".jpeg", ".tif", ".tiff", ".heic")
_NAME = re.compile(r"(patient|export|lab_?report|identifiable|private_|denylist|credential|"
                   r"secret|token|\.env$|\.pem$|\.key$)", re.I)

_EMAIL = re.compile(r"[A-Za-z0-9._%+\-]+@([A-Za-z0-9\-]+(?:\.[A-Za-z0-9\-]+)+)")
_EMAIL_OK_DOMAINS = ("example.com", "example.org", "example.net", "users.noreply.github.com",
                     "anthropic.com")


def _email_domain_ok(domain: str) -> bool:
    """The domain is an allowed one or a sub-domain of it (never a longer name that merely
    ends with the same letters: ``notexample.com`` is not ``example.com``)."""
    d = domain.lower().rstrip(".")
    return any(d == ok or d.endswith("." + ok) for ok in _EMAIL_OK_DOMAINS)


_PHONE = re.compile(r"(?<![\w+])\+\d{1,3}[ .\-]?\(?\d{1,4}\)?[ .\-]?\d{2,4}[ .\-]?\d{2,4}"
                    r"(?:[ .\-]?\d{2,4})?(?!\d)")
_LONGNUM = re.compile(r"(?<![A-Za-z0-9_.\-])\d{9,}(?![A-Za-z0-9_.\-])")
_DOB = re.compile(r"\b(?:dob|date of birth|birth ?date)\b\s*[:=]\s*[\"']?\w", re.I)
_HEADER = re.compile(r"(?i)(^|[;,\t])\s*\"?(patient_?id|patientid|patient_ID)\"?\s*([;,\t]|$)")
_GENO = re.compile(r"(^|[;,\t])(Het|Hom)([;,\t]|$)")
_CRED = re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----|AKIA[0-9A-Z]{16}|"
                   r"gh[pousr]_[A-Za-z0-9]{36,}|xox[baprs]-[A-Za-z0-9\-]{10,}|"
                   r"sk-[A-Za-z0-9]{32,}")
_TOKEN = re.compile(r"[A-Za-z0-9_.\-]{4,64}")


# ------------------------------------------------------------- denylist ------
class DenylistError(ValueError):
    """The private denylist is unusable. The message never contains an entry."""


@dataclass
class Denylist:
    literals: list[tuple[int, str]] = field(default_factory=list)
    regexes: list[tuple[int, re.Pattern]] = field(default_factory=list)
    digests: dict[str, int] = field(default_factory=dict)
    source: str = ""

    def __len__(self) -> int:
        return len(self.literals) + len(self.regexes) + len(self.digests)

    def scan(self, text: str) -> Counter:
        """Count matches per entry number; values are never returned."""
        hits: Counter = Counter()
        low = text.lower()
        for n, lit in self.literals:
            c = low.count(lit)
            if c:
                hits[n] += c
        for n, rx in self.regexes:
            c = len(rx.findall(text))
            if c:
                hits[n] += c
        if self.digests:
            for tok in _TOKEN.findall(text):
                n = self.digests.get(hashlib.sha256(tok.lower().encode()).hexdigest())
                if n:
                    hits[n] += 1
        return hits


def _inside(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def load_denylist(path: str | os.PathLike | None, repo_root: Path) -> Denylist | None:
    """Read the private denylist, or return None if none is configured or readable.

    Entry syntax, one per line (``#`` starts a comment):
        ``text``              case-insensitive literal (at least 4 characters);
        ``literal:text``      the same, explicitly;
        ``regex:pattern``     a Python regular expression;
        ``sha256:hex``        the SHA-256 of a lower-cased token (a whole word
                              made of letters, digits, ``_ . -``).
    The file must live outside the repository; a path inside it is an error,
    whatever its Git status, so the list cannot be committed by accident.
    """
    if not path:
        return None
    p = Path(path).expanduser()
    if not p.exists():
        return None
    if _inside(p, repo_root):
        raise DenylistError("the private denylist must be stored outside the repository")
    dl = Denylist(source=p.name)
    for n, raw in enumerate(p.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        kind, _, rest = line.partition(":")
        kind = kind.lower()
        if kind == "regex":
            try:
                dl.regexes.append((n, re.compile(rest)))
            except re.error:
                raise DenylistError(f"denylist entry {n}: invalid regular expression") from None
        elif kind == "sha256":
            h = rest.strip().lower()
            if not re.fullmatch(r"[0-9a-f]{64}", h):
                raise DenylistError(f"denylist entry {n}: not a SHA-256 digest")
            dl.digests[h] = n
        else:
            text = rest if kind == "literal" else line
            if len(text) < 4:
                raise DenylistError(f"denylist entry {n}: shorter than 4 characters")
            dl.literals.append((n, text.lower()))
    return dl


# -------------------------------------------------------------- scanning ----
@dataclass(frozen=True)
class Finding:
    path: str
    rule: str
    count: int
    entry: int | None = None   # denylist entry number only

    def label(self) -> str:
        return f"{self.rule}" + (f" (entry {self.entry})" if self.entry else "")


def _count_field(f: Finding) -> int | str:
    """Path rules flag a file (the count is always one); content counts are suppressed if small."""
    return f.count if f.rule.startswith("PATH-") else shown(f.count)


@dataclass
class ScanResult:
    findings: list[Finding] = field(default_factory=list)
    files_scanned: int = 0
    files_skipped_binary: int = 0
    files_skipped_large: int = 0
    denylist_status: str = "not configured"
    denylist_entries: int = 0

    @property
    def clean(self) -> bool:
        return not self.findings

    def to_dict(self) -> dict:
        return {
            "clean": self.clean,
            "files_scanned": self.files_scanned,
            "files_skipped_binary": self.files_skipped_binary,
            "files_skipped_large": self.files_skipped_large,
            "denylist": {"status": self.denylist_status, "entries": self.denylist_entries},
            "small_counts": "counts below 5 are shown as '<5'",
            "findings": [{"path": f.path, "rule": f.rule, "count": _count_field(f),
                          **({"denylist_entry": f.entry} if f.entry else {})}
                         for f in self.findings],
            "note": "Paths, rule identifiers and counts only. Pattern scanning is not a "
                    "guarantee that no sensitive data are present.",
        }

    def render(self) -> str:
        head = "no findings" if self.clean else f"{len(self.findings)} finding(s)"
        dl = self.denylist_status + (
            f" ({self.denylist_entries} entries)" if self.denylist_entries else "")
        lines = [f"privacy scan: {head} in {self.files_scanned} file(s); private denylist: {dl}"]
        if self.files_skipped_binary or self.files_skipped_large:
            lines.append(f"  not scanned: {self.files_skipped_binary} binary, "
                         f"{self.files_skipped_large} over {MAX_BYTES // 1_000_000} MB")
        for f in sorted(self.findings, key=lambda f: (f.path, f.rule)):
            lines.append(f"  {f.path}: {f.label()}"
                         + ("" if f.rule.startswith("PATH-") else f" x{_n(f.count)}"))
        lines.append("  (a clean result does not guarantee that no sensitive data are present)")
        return "\n".join(lines)


def path_findings(path: str) -> list[Finding]:
    """Rules about where a file is and what it is called, not what it contains."""
    out: list[Finding] = []
    norm = path.replace("\\", "/")
    name = norm.rsplit("/", 1)[-1]
    if norm.startswith(("data/", "results/", "logs/")) and norm not in _ALLOWED_DATA \
            and not _DATA_OK.search(norm):
        out.append(Finding(path, "PATH-DATA", 1))
    if norm.lower().endswith(_TYPE_EXT):
        out.append(Finding(path, "PATH-TYPE", 1))
    if _NAME.search(name) and not name.endswith(("privacy.py", "check_privacy.py")) \
            and not name.startswith(("test_privacy", "privacy_policy")) \
            and "example" not in name and norm != "release_allowlist.txt":
        out.append(Finding(path, "PATH-NAME", 1))
    return out


def content_findings(path: str, text: str, denylist: Denylist | None) -> list[Finding]:
    """Rules about file content. Only counts leave this function."""
    out: list[Finding] = []

    n = sum(1 for m in _EMAIL.finditer(text) if not _email_domain_ok(m.group(1)))
    if n:
        out.append(Finding(path, "CONTENT-EMAIL", n))
    n = len(_PHONE.findall(text))
    if n:
        out.append(Finding(path, "CONTENT-PHONE", n))
    # lines that name an epoch/timestamp are build constants, not record numbers
    n = sum(len(_LONGNUM.findall(ln)) for ln in text.splitlines()
            if "epoch" not in ln.lower() and "timestamp" not in ln.lower())
    if n:
        out.append(Finding(path, "CONTENT-LONGNUM", n))
    n = len(_DOB.findall(text))
    if n:
        out.append(Finding(path, "CONTENT-DOB", n))
    n = len(_CRED.findall(text))
    if n:
        out.append(Finding(path, "CONTENT-CREDENTIAL", n))

    # A "synthetic data" banner inside a file proves nothing: real data can carry the same
    # words. Table and genotype checks therefore always run. Only the files named in
    # SYNTHETIC_EXAMPLES (by path, and only if the path is of a documented example) are exempt.
    if not _synthetic_example(path):
        lines = text.splitlines()
        for i, line in enumerate(lines):
            if _HEADER.search(line) and any(d in line for d in ";\t,"):
                delim = max(";\t,", key=line.count)
                rows = [x for x in lines[i + 1:i + 6] if x.count(delim) >= 2]
                if len(rows) >= 3:
                    out.append(Finding(path, "CONTENT-TABLE", 1))
                    break
        geno = sum(1 for line in lines
                   if _GENO.search(line) and sum(map(line.count, ";\t,")) >= 3)
        if geno >= 3:
            out.append(Finding(path, "CONTENT-GENOTYPE", geno))

    if denylist is not None:
        for entry, count in sorted(denylist.scan(text).items()):
            out.append(Finding(path, "DENYLIST", count, entry))
    return out


def _synthetic_example(path: str) -> bool:
    return fnmatch_any(path.replace("\\", "/"), SYNTHETIC_EXAMPLES)


def _is_text(data: bytes) -> bool:
    return b"\x00" not in data[:8192]


class _Symlink:
    """Marks an item that is a symbolic link: it is reported and never followed."""

    def __repr__(self) -> str:
        return "SYMLINK"


SYMLINK = _Symlink()
Item = tuple[str, "bytes | None | _Symlink"]


def scan_items(items: Iterable[Item], denylist: Denylist | None = None,
               path_rules: bool = True) -> ScanResult:
    """Scan (path, content) pairs. ``content`` None means "check the path only"; ``SYMLINK``
    means the item is a symbolic link, which is reported and not read."""
    res = ScanResult()
    if denylist is not None:
        res.denylist_status = "active"
        res.denylist_entries = len(denylist)
    for path, data in items:
        if path_rules:
            res.findings += path_findings(path)
        if data is SYMLINK:
            res.findings.append(Finding(path, "PATH-SYMLINK", 1))
            continue
        if data is None:
            continue
        if len(data) > MAX_BYTES:
            res.files_skipped_large += 1
            continue
        if not _is_text(data):
            res.files_skipped_binary += 1
            continue
        res.files_scanned += 1
        res.findings += content_findings(path, data.decode("utf-8", errors="replace"), denylist)
    return res


# ------------------------------------------------------------ git helpers ----
class ScanError(RuntimeError):
    """The scan could not be completed. A scan that cannot see all its files never reports
    a clean result. The message holds no path, file name or Git output."""


def _git(root: Path, *args: str) -> subprocess.CompletedProcess:
    """Run Git; a failure to run it or a non-zero status raises :class:`ScanError`."""
    try:
        r = subprocess.run(["git", "-C", str(root), *args], capture_output=True)
    except OSError:
        raise ScanError("Git could not be run") from None
    if r.returncode != 0:
        raise ScanError(f"a Git command failed (git {args[0]}, status {r.returncode})")
    return r


def _git_probe(root: Path, *args: str) -> subprocess.CompletedProcess:
    """Run Git where a non-zero status is an answer, not a failure (a query)."""
    try:
        return subprocess.run(["git", "-C", str(root), *args], capture_output=True)
    except OSError:
        raise ScanError("Git could not be run") from None


def repo_root(start: Path) -> Path:
    r = _git_probe(start, "rev-parse", "--show-toplevel")
    if r.returncode != 0:
        raise ScanError("not inside a Git repository")
    return Path(r.stdout.decode().strip())


def _read(f: Path) -> bytes | None | _Symlink:
    """Content of a working-tree file. A symlink is never followed; a path that is not a
    regular file has no content; an unreadable file is a scan failure."""
    if f.is_symlink():
        return SYMLINK
    if not f.is_file():
        return None
    try:
        return f.read_bytes()
    except OSError:
        raise ScanError("a file could not be read") from None


def _names(out: bytes) -> list[str]:
    return [n for n in out.decode("utf-8", errors="surrogateescape").split("\0") if n]


def tracked_items(root: Path) -> list[Item]:
    """Tracked files, content from the working tree (what the next commit would carry)."""
    return [(n, _read(root / n)) for n in _names(_git(root, "ls-files", "-z").stdout)]


def working_tree_items(root: Path) -> list[Item]:
    """Tracked files plus untracked files that Git does not ignore (what `git add .` would take)."""
    names = _names(_git(root, "ls-files", "-z", "--cached", "--others",
                        "--exclude-standard").stdout)
    return [(n, _read(root / n)) for n in sorted(set(names))]


def staged_items(root: Path) -> list[Item]:
    """Added/modified files in the index, with the staged content (not the working copy)."""
    raw = _git(root, "diff", "--cached", "--raw", "-z", "--no-abbrev",
               "--diff-filter=ACMR").stdout.decode("utf-8", errors="surrogateescape")
    toks = [t for t in raw.split("\0") if t]
    out: list[Item] = []
    i = 0
    while i < len(toks):
        meta = toks[i].split()          # ":old_mode new_mode old_sha new_sha status"
        status = meta[4][0]
        n_paths = 2 if status in "CR" else 1
        path = toks[i + n_paths]        # the destination path
        i += 1 + n_paths
        if meta[1] == "120000":
            out.append((path, SYMLINK))
        else:
            out.append((path, _git(root, "show", f":{path}").stdout))
    return out


def path_items(paths: Sequence[str | os.PathLike], root: Path) -> list[Item]:
    """Files under the given paths (generated outputs, a notebook folder, ...).

    Symbolic links, to files or folders, are reported and never followed, so a link cannot
    carry the scan outside the folder it was asked to look at.
    """
    def label(f: Path) -> str:
        try:
            return str(f.resolve().relative_to(root.resolve()))
        except ValueError:
            return f.name

    out: list[Item] = []
    for p in map(Path, paths):
        if p.is_symlink():
            out.append((p.name, SYMLINK))
        elif p.is_file():
            out.append((label(p), _read(p)))
        elif p.is_dir():
            for dirpath, dirnames, filenames in os.walk(p, followlinks=False):
                here = Path(dirpath)
                for d in sorted(dirnames):
                    if (here / d).is_symlink():
                        out.append((label(here / d), SYMLINK))
                dirnames[:] = sorted(d for d in dirnames if not (here / d).is_symlink())
                for fn in sorted(filenames):
                    out.append((label(here / fn), _read(here / fn)))
        else:
            raise ScanError("a path given to --path does not exist")
    return out


def tracked_denylist_files(root: Path, denylist_path: str | os.PathLike | None) -> bool:
    """True if the denylist file is tracked by Git (which must never be the case)."""
    if not denylist_path:
        return False
    p = Path(denylist_path).expanduser().resolve()
    if not _inside(p, root):
        return False
    rel = str(p.relative_to(root.resolve()))
    return _git_probe(root, "ls-files", "--error-unmatch", rel).returncode == 0


def fnmatch_any(name: str, patterns: Iterable[str]) -> bool:
    return any(fnmatch.fnmatchcase(name, p) for p in patterns)
