# Reproducibility and environment

## Environment

* Python 3.11 or newer (`requires-python`). CI runs 3.11 and 3.13. Python 3.12
  is declared supported by the metadata but is not exercised by CI or by the
  remediation test runs.
* Runtime dependencies are declared with lower bounds in `pyproject.toml`
  (`numpy`, `pandas`, `scipy`, `matplotlib`). Fully pinned environments,
  including transitive dependencies and test tools, are in
  `requirements-lock-py3.11.txt` and `requirements-lock-py3.13.txt`. CI and
  `make clean-clone-check` install exactly those files, then `pip install
  --no-deps -e .`.
* Set up a local environment:

```bash
python3.13 -m venv .venv && . .venv/bin/activate
pip install -r requirements-lock-py3.13.txt
pip install --no-deps -e .
make check          # lint, privacy scan, tests
make demo-check     # synthetic demonstration, run twice, outputs compared
```

## What is deterministic

* Random seeds are fixed in the code and recorded in the manifest. The bootstrap
  replicate count is a setting (`MEFV_BOOTSTRAP_REPLICATES`) and is recorded.
* `SOURCE_DATE_EPOCH` is pinned so that PDF metadata is stable.
* Repeated runs on the same inputs give byte-identical tables, figures and
  intermediates. The run manifest's `deterministic` section is identical; its
  `execution` section is not, by design.
* Output depends on the exact versions of numpy, pandas, scipy and matplotlib.
  Byte identity is claimed per locked environment, not across versions.

## What the validation and manifest add

Stage 01 writes `validation_report.json` (counts only). It is an intermediate
recorded in the manifest, so a change in validation status shows up in
`scripts/compare_manifests.py`. A run whose blocking findings carry an owner disposition has the status
`owner_disposition_recorded`; one that only acknowledged them does not produce
outputs at all ([`validation.md`](validation.md)).

## Continuous integration

`.github/workflows/ci.yml` runs on every push and pull request: lint, the
privacy scan (built-in rules only; the private denylist is never available to
CI), the full tests and the synthetic demonstration, on Python 3.11 and 3.13,
with `contents: read` permission and no secrets.

**Not yet done: pinning the actions by commit SHA.** `actions/checkout@v4` and
`actions/setup-python@v5` are tags that can move. Pinning needs the SHAs
resolved from GitHub, which could not be done during the remediation; they must
not be guessed. See [`owner_decisions.md`](owner_decisions.md).

## Not claimed

Reproducibility is a statement about the software on given inputs. It says
nothing about whether the inputs, the assumptions or the interpretations are
correct.
