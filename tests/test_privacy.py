"""Tests for the privacy scanner.

Sentinel strings are assembled at run time from fragments so this file does not itself contain
anything the scanner would flag, and every test that checks output asserts the sentinel is absent
from stdout, stderr and the JSON report. All values are invented.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

from mefv import privacy

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "check_privacy.py"

EMAIL = "jane.doe" + "@" + "clinic.test"
PHONE = "+" + "44 20 7946 0958"
LONGNUM = "1234" + "56789012"
NAME = "Zorblax" + "Quenfell"          # invented surname for the denylist
HET = "He" + "t"                       # assembled: a literal genotype column would flag this file
SENTINELS = [EMAIL, PHONE, LONGNUM, NAME]


def git(cwd, *args):
    subprocess.run(["git", "-C", str(cwd), *args], check=True, capture_output=True)


@pytest.fixture
def repo(tmp_path):
    r = tmp_path / "repo"
    r.mkdir()
    git(r, "init", "-q")
    git(r, "config", "user.email", "t@example.org")
    git(r, "config", "user.name", "t")
    (r / "README.md").write_text("# demo\nnothing here\n")
    git(r, "add", ".")
    git(r, "commit", "-qm", "init")
    return r


def run(repo, *args, env=None):
    import os
    e = {k: v for k, v in os.environ.items() if k != privacy.DENYLIST_ENV}
    e.update(env or {})
    return subprocess.run([sys.executable, str(SCRIPT), "--root", str(repo), *args],
                          capture_output=True, text=True, env=e)


def assert_no_values(*outputs):
    for out in outputs:
        for s in SENTINELS:
            assert s not in out
            assert s.lower() not in out.lower()


# --- clean / baseline ---------------------------------------------------------------------
def test_clean_repo_passes_and_says_no_denylist(repo):
    r = run(repo)
    assert r.returncode == 0
    assert "no findings" in r.stdout and "not configured" in r.stdout
    assert "not a guarantee" in r.stdout or "does not guarantee" in r.stdout


def test_this_repository_is_clean():
    if subprocess.run(["git", "-C", str(ROOT), "rev-parse"], capture_output=True).returncode:
        pytest.skip("not a Git checkout")
    r = subprocess.run([sys.executable, str(SCRIPT), "--root", str(ROOT)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout


# --- general content rules ----------------------------------------------------------------
@pytest.mark.parametrize("text,rule", [
    (f"contact {EMAIL}\n", "CONTENT-EMAIL"),
    (f"tel {PHONE}\n", "CONTENT-PHONE"),
    (f"mrn {LONGNUM}\n", "CONTENT-LONGNUM"),
    ("date of " + "birth: 1970-01-01\n", "CONTENT-DOB"),
    ("-----BEGIN " + "RSA PRIVATE KEY-----\n", "CONTENT-CREDENTIAL"),
    ("patient_id,variant,zyg\nA1,x,Het\nA2,y,Het\nA3,z,Hom\n", "CONTENT-TABLE"),
    ("a;b;c;Het\nd;e;f;Het\ng;h;i;Hom\n", "CONTENT-GENOTYPE"),
])
def test_content_rule_detects_and_never_echoes(repo, text, rule):
    (repo / "notes.txt").write_text(text)
    r = run(repo, "--json", str(repo.parent / "rep.json"))
    assert r.returncode == 1
    assert f"notes.txt: {rule}" in r.stdout
    report = (repo.parent / "rep.json").read_text()
    assert rule in report
    assert_no_values(r.stdout, r.stderr, report)


def test_example_domain_email_and_epoch_constant_are_not_flagged(repo):
    (repo / "a.txt").write_text("write to someone@example.org\nSOURCE_DATE_EPOCH=1767225600\n")
    assert run(repo).returncode == 0


def test_a_synthetic_banner_does_not_exempt_a_file_from_the_table_rules(repo):
    """Real data can carry the same words in a header line: the banner proves nothing."""
    rows = "".join(f"{k},x,{HET},1\n" for k in "ABC")
    (repo / "s.csv").write_text("# SYNTHETIC DATA\npatient_id,v,z,w\n" + rows)
    r = run(repo)
    assert r.returncode == 1
    assert "s.csv: CONTENT-TABLE" in r.stdout and "s.csv: CONTENT-GENOTYPE" in r.stdout
    (repo / "s2.csv").write_text(f"# SYNTHETIC DATA\nmail {EMAIL}\n")
    r = run(repo)
    assert r.returncode == 1 and "s2.csv: CONTENT-EMAIL" in r.stdout


def test_only_documented_example_paths_are_exempt_from_the_table_rules(repo):
    table = "patient_id,v,z,w\n" + "".join(f"{k},x,{HET},1\n" for k in "ABC")
    ok = repo / "data" / "external" / "population_frequencies.example.csv"
    ok.parent.mkdir(parents=True)
    ok.write_text(table)
    assert run(repo).returncode == 0
    elsewhere = repo / "data" / "external" / "population_frequencies.csv"
    elsewhere.write_text(table)                       # same content, not a documented example
    r = run(repo)
    assert r.returncode == 1 and "population_frequencies.csv: CONTENT-TABLE" in r.stdout


@pytest.mark.parametrize("domain,flagged", [
    ("example.org", False), ("mail.example.org", False), ("users.noreply.github.com", False),
    ("notexample.org", True),                          # merely ends with the allowed letters
    ("evil-example.com", True), ("example.org.evil.test", True),
    ("xusers.noreply.github.com", True),
])
def test_email_allowlist_matches_whole_domains_only(repo, domain, flagged):
    addr = "someone" + "@" + domain                    # assembled here: not an address in the file
    (repo / "a.txt").write_text(f"write to {addr}\n")
    r = run(repo)
    assert (r.returncode == 1) is flagged
    if flagged:
        assert "CONTENT-EMAIL" in r.stdout and addr not in r.stdout


# --- path rules ---------------------------------------------------------------------------
@pytest.mark.parametrize("rel,rule", [
    ("data/raw/cohort.csv", "PATH-DATA"),
    ("results/table1.csv", "PATH-DATA"),
    ("export.xlsx", "PATH-TYPE"),
    ("backup.zip", "PATH-TYPE"),
    ("patient_list.txt", "PATH-NAME"),
    ("private_notes.txt", "PATH-NAME"),
])
def test_path_rules(repo, rel, rule):
    f = repo / rel
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_bytes(b"x\n")
    r = run(repo)
    assert r.returncode == 1 and f"{rel}: {rule}" in r.stdout


def test_documented_placeholders_are_allowed(repo):
    for rel in ("data/raw/.gitkeep", "results/.gitkeep", "data/external/freq.example.csv"):
        f = repo / rel
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text("")
    assert run(repo).returncode == 0


def test_binary_files_are_counted_not_read(repo):
    (repo / "blob.dat").write_bytes(b"\x00\x01" + EMAIL.encode())
    r = run(repo)
    assert r.returncode == 0 and "1 binary" in r.stdout


# --- private denylist ---------------------------------------------------------------------
def make_denylist(tmp_path, *lines):
    p = tmp_path / "private" / "denylist.txt"
    p.parent.mkdir(exist_ok=True)
    p.write_text("\n".join(lines) + "\n")
    return p


def test_denylist_literal_is_detected_case_insensitively_without_echo(repo, tmp_path):
    dl = make_denylist(tmp_path, "# comment", NAME)
    (repo / "doc.md").write_text(f"the clinician {NAME.upper()} reviewed this\n")
    r = run(repo, "--denylist", str(dl), "--json", str(tmp_path / "rep.json"))
    assert r.returncode == 1
    assert "doc.md: DENYLIST (entry 2) x<5" in r.stdout
    assert "active" in r.stdout
    assert_no_values(r.stdout, r.stderr, (tmp_path / "rep.json").read_text())


def test_denylist_via_environment_variable(repo, tmp_path):
    dl = make_denylist(tmp_path, NAME)
    (repo / "doc.md").write_text(NAME + "\n")
    r = run(repo, env={privacy.DENYLIST_ENV: str(dl)})
    assert r.returncode == 1 and "DENYLIST" in r.stdout


def test_denylist_regex_and_hash_entries(repo, tmp_path):
    digest = hashlib.sha256(NAME.lower().encode()).hexdigest()
    dl = make_denylist(tmp_path, r"regex:ZZ-\d{4}", f"sha256:{digest}")
    (repo / "a.md").write_text("code ZZ-1234 appears\n")
    (repo / "b.md").write_text(f"person {NAME} appears\n")
    r = run(repo, "--denylist", str(dl))
    assert r.returncode == 1
    assert "a.md: DENYLIST (entry 1)" in r.stdout and "b.md: DENYLIST (entry 2)" in r.stdout
    assert_no_values(r.stdout, r.stderr)
    assert "ZZ-1234" not in r.stdout and digest not in r.stdout


def test_denylist_inside_repository_is_refused(repo):
    dl = repo / "denylist.txt"
    dl.write_text(NAME + "\n")
    r = run(repo, "--denylist", str(dl))
    assert r.returncode == 2
    assert "outside the repository" in r.stderr
    assert_no_values(r.stdout, r.stderr)


@pytest.mark.parametrize("line", ["abc", "regex:([", "sha256:nothex"])
def test_malformed_denylist_is_refused_without_echoing_entry(repo, tmp_path, line):
    dl = make_denylist(tmp_path, line)
    r = run(repo, "--denylist", str(dl))
    assert r.returncode == 2
    assert "entry 1" in r.stderr
    assert line not in r.stderr


def test_a_configured_but_missing_denylist_stops_the_scan(repo, tmp_path):
    """Asking for a denylist and not getting one must not become a quietly weaker scan."""
    r = run(repo, "--denylist", str(tmp_path / "absent.txt"))
    assert r.returncode == 2 and "not found" in r.stderr and "no findings" not in r.stdout
    r = run(repo, env={privacy.DENYLIST_ENV: str(tmp_path / "absent.txt")})
    assert r.returncode == 2


def test_no_denylist_at_all_degrades_and_says_so(repo):
    r = run(repo)
    assert r.returncode == 0 and "denylist: not configured" in r.stdout


def test_require_denylist_fails_closed(repo):
    r = run(repo, "--require-denylist")
    assert r.returncode == 2


# --- modes --------------------------------------------------------------------------------
def test_staged_mode_checks_index_content_not_working_copy(repo):
    f = repo / "n.txt"
    f.write_text(f"mail {EMAIL}\n")
    git(repo, "add", "n.txt")
    f.write_text("clean now\n")            # working copy cleaned after staging
    r = run(repo, "--staged")
    assert r.returncode == 1 and "n.txt: CONTENT-EMAIL" in r.stdout
    assert run(repo, "--tracked").returncode == 0     # n.txt not committed


def test_staged_mode_ignores_unstaged(repo):
    (repo / "n.txt").write_text(f"mail {EMAIL}\n")
    assert run(repo, "--staged").returncode == 0


def test_path_mode_scans_generated_outputs_content_only(repo):
    out = repo / "outputs"
    out.mkdir()
    (out / "table.csv").write_text(f"x\n{LONGNUM}\n")
    r = run(repo, "--path", str(out))
    assert r.returncode == 1 and "CONTENT-LONGNUM" in r.stdout
    assert "PATH-DATA" not in r.stdout
    assert_no_values(r.stdout, r.stderr)


def test_notebook_json_is_scanned(repo):
    nb = {"cells": [{"cell_type": "markdown", "source": [f"contact {EMAIL}"]}]}
    (repo / "analysis.ipynb").write_text(json.dumps(nb))
    r = run(repo)
    assert r.returncode == 1 and "analysis.ipynb: CONTENT-EMAIL" in r.stdout


def test_not_a_repository_exits_2(tmp_path):
    r = subprocess.run([sys.executable, str(SCRIPT), "--root", str(tmp_path)],
                       capture_output=True, text=True)
    assert r.returncode == 2


# --- library level ------------------------------------------------------------------------
def test_finding_objects_hold_no_values():
    res = privacy.scan_items([("f.txt", f"{EMAIL} {NAME}".encode())],
                             privacy.Denylist(literals=[(1, NAME.lower())]))
    blob = repr(res.findings) + json.dumps(res.to_dict()) + res.render()
    assert_no_values(blob)
    assert {f.rule for f in res.findings} == {"CONTENT-EMAIL", "DENYLIST"}


def test_every_rule_has_a_description():
    used = {"PATH-DATA", "PATH-TYPE", "PATH-NAME", "PATH-SYMLINK", "CONTENT-EMAIL", "CONTENT-PHONE",
            "CONTENT-LONGNUM", "CONTENT-DOB", "CONTENT-TABLE", "CONTENT-GENOTYPE",
            "CONTENT-CREDENTIAL", "DENYLIST"}
    assert used == set(privacy.RULES)


@pytest.mark.parametrize("name", ["denylist.txt", "private_denylist.txt", "x/my.denylist",
                                  "owner_disposition.json", "study.disposition.json"])
def test_private_owner_files_are_ignored_by_git(name):
    if subprocess.run(["git", "-C", str(ROOT), "rev-parse"], capture_output=True).returncode:
        pytest.skip("not a Git checkout")
    r = subprocess.run(["git", "-C", str(ROOT), "check-ignore", "-q", name])
    assert r.returncode == 0, f"{name} is not ignored by .gitignore"


# --- symbolic links ------------------------------------------------------------------------
@pytest.fixture
def outside(tmp_path):
    """A folder outside the repository holding something the scanner must never read."""
    d = tmp_path / "outside"
    d.mkdir()
    (d / "secret.txt").write_text(f"mail {EMAIL}\n")
    return d


def test_a_tracked_symlink_is_reported_and_its_target_is_not_read(repo, outside):
    (repo / "link.txt").symlink_to(outside / "secret.txt")
    git(repo, "add", "link.txt")
    for mode in ("--tracked", "--staged", None):
        r = run(repo, *([mode] if mode else []))
        assert r.returncode == 1 and "link.txt: PATH-SYMLINK" in r.stdout, mode
        assert "CONTENT-EMAIL" not in r.stdout           # the target was not followed
        assert_no_values(r.stdout, r.stderr)


def test_an_untracked_symlink_is_reported_in_the_working_tree_scan(repo, outside):
    (repo / "dir_link").symlink_to(outside, target_is_directory=True)
    r = run(repo)
    assert r.returncode == 1 and "dir_link: PATH-SYMLINK" in r.stdout
    assert_no_values(r.stdout, r.stderr)


def test_path_mode_never_follows_links_out_of_the_folder(repo, outside):
    out = repo / "gen"
    out.mkdir()
    (out / "ok.txt").write_text("fine\n")
    (out / "file_link").symlink_to(outside / "secret.txt")
    (out / "dir_link").symlink_to(outside, target_is_directory=True)
    r = run(repo, "--path", str(out))
    assert r.returncode == 1
    assert "file_link: PATH-SYMLINK" in r.stdout and "dir_link: PATH-SYMLINK" in r.stdout
    assert "CONTENT-EMAIL" not in r.stdout
    assert_no_values(r.stdout, r.stderr)
    r = run(repo, "--path", str(out / "dir_link"))      # a link given directly is not followed
    assert r.returncode == 1 and "PATH-SYMLINK" in r.stdout and "CONTENT-EMAIL" not in r.stdout


def test_path_mode_with_a_missing_path_fails(repo, tmp_path):
    r = run(repo, "--path", str(tmp_path / "nope"))
    assert r.returncode == 2 and "does not exist" in r.stderr


# --- a failing Git command is never a clean result ------------------------------------------
@pytest.fixture
def broken_git(tmp_path, monkeypatch):
    """A Git that fails every command except the one that locates the repository."""
    real = subprocess.run
    def fake(cmd, *a, **k):
        if cmd[:1] == ["git"] and "rev-parse" not in cmd:
            return subprocess.CompletedProcess(cmd, 128, b"", b"fatal: path /private/name\n")
        return real(cmd, *a, **k)
    monkeypatch.setattr(privacy.subprocess, "run", fake)


@pytest.mark.parametrize("fn", ["tracked_items", "working_tree_items", "staged_items"])
def test_a_failing_git_command_raises_instead_of_returning_no_files(repo, broken_git, fn):
    with pytest.raises(privacy.ScanError) as e:
        getattr(privacy, fn)(repo)
    assert "/private/name" not in str(e.value)           # no Git output is repeated


def test_the_cli_exits_2_when_git_fails(repo, tmp_path):
    fake = tmp_path / "bin"
    fake.mkdir()
    (fake / "git").write_text("#!/bin/sh\n"
                              "case \"$*\" in *rev-parse*) exec /usr/bin/env -u PATH "
                              "$REAL_GIT \"$@\" ;; *) echo fatal >&2; exit 128 ;; esac\n")
    (fake / "git").chmod(0o755)
    import os
    import shutil
    real = shutil.which("git")
    env = {"PATH": f"{fake}:{os.environ['PATH']}", "REAL_GIT": real}
    for mode in ("--tracked", "--staged"):
        r = run(repo, mode, env=env)
        assert r.returncode == 2 and "could not be completed" in r.stderr
        assert "no findings" not in r.stdout


def test_an_unreadable_file_fails_the_scan(repo):
    import os
    f = repo / "locked.txt"
    f.write_text("x\n")
    f.chmod(0)
    try:
        if os.access(f, os.R_OK):
            pytest.skip("permissions are not enforced for this user")
        r = run(repo)
        assert r.returncode == 2 and "could not be read" in r.stderr
    finally:
        f.chmod(0o644)


# --- small counts are never shown ----------------------------------------------------------
def test_small_counts_are_suppressed_in_text_and_json_but_large_ones_are_shown(repo, tmp_path):
    (repo / "few.txt").write_text(f"mail {EMAIL}\nmail {EMAIL}\n")                 # 2 matches
    (repo / "many.txt").write_text("".join(f"mail {i}" + EMAIL + "\n" for i in range(7)))
    out = tmp_path / "rep.json"
    r = run(repo, "--json", str(out))
    assert r.returncode == 1
    assert "few.txt: CONTENT-EMAIL x<5" in r.stdout and "many.txt: CONTENT-EMAIL x7" in r.stdout
    rep = json.loads(out.read_text())
    counts = {f["path"]: f["count"] for f in rep["findings"] if f["rule"] == "CONTENT-EMAIL"}
    assert counts == {"few.txt": "<5", "many.txt": 7}
    assert_no_values(r.stdout, r.stderr, out.read_text())
