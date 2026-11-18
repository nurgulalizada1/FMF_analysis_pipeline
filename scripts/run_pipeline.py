#!/usr/bin/env python3
"""Run the pipeline stages in order and record what each one did.

Stages run as separate processes in the order given by their numeric prefix.
Each stage inherits the environment, so configuration is set once, outside.

Outputs
-------
* **Scientific artifacts** - every file under the tables, figures and
  classification directories. These are the results.
* **Intermediates** - the two processed analysis tables, the generated
  annotations and the merged variant dictionary.

Both are deterministic: the same inputs and settings give byte-identical
files.

The run manifest
----------------
The runner writes ``run_manifest.json`` (``MEFV_RUN_MANIFEST``; by default
next to the tables directory, never inside an output directory). It has two
sections:

* ``deterministic`` - digests of every artifact and intermediate, a combined
  digest of the artifacts, the analysis settings and each stage's exit code.
  Two runs on the same inputs must produce an identical ``deterministic``
  section; this is what the determinism tests compare.
* ``execution`` - intentionally non-deterministic or machine-specific
  metadata: the time of the run, stage durations, absolute paths, package
  versions, the git commit and the provenance of the reference data. It is
  kept because it is useful, and excluded from every reproducibility claim.

The manifest itself, any other file named ``run_manifest.json`` and hidden
files are never collected as outputs, so repeated runs cannot collect a
previous manifest.

Clean runs
----------
``--clean`` empties the output directories before a full run. To make that
safe, the runner only empties a directory that it created itself (marked by a
``.mefv-output`` file) or one with no visible contents; anything else is
refused. Hidden files such as ``.gitkeep`` are kept. Without
``--clean``, files in the output directories that the current full run did
not write are reported as stale and left out of the manifest's outputs.

Usage
-----
    python scripts/run_pipeline.py                 # every stage
    python scripts/run_pipeline.py --clean         # every stage, fresh output directories
    python scripts/run_pipeline.py --from 04       # from stage 04 onwards
    python scripts/run_pipeline.py --only 02 03    # selected stages
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

# Make the package importable when the script is run without installing it.
import pathlib as _pathlib
import sys as _sys

_sys.path.insert(0, str(_pathlib.Path(__file__).resolve().parents[1] / "src"))

from mefv import config  # noqa: E402

HERE = Path(__file__).resolve().parent
STAGES = sorted(p for p in HERE.glob("[0-9][0-9]_*.py"))
MARKER = ".mefv-output"
MANIFEST_NAME = "run_manifest.json"
#: Written into every output directory when a stage fails, so that files left from an
#: earlier run are not mistaken for the results of this one. Removed at the start of the
#: next run and never collected as an output.
FAILED_NAME = "RUN_FAILED_outputs_are_not_current.txt"
MANIFEST_VERSION = 2


def digest(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def artifact_dirs() -> list[Path]:
    return [config.TABLES, config.FIGURES, config.CLASSIFICATION]


def intermediate_files() -> list[Path]:
    return [config.PATIENTS, config.VARIANTS, config.VALIDATION_REPORT,
            config.ANNOTATIONS_GENERATED, config.DICTIONARY]


def _is_output(f: Path, manifest: Path) -> bool:
    return (f.is_file() and not f.name.startswith(".")
            and f.name not in (MANIFEST_NAME, FAILED_NAME) and f.resolve() != manifest.resolve())


def collect(directories: list[Path], manifest: Path, since: float | None
            ) -> tuple[dict[str, dict], list[str]]:
    """Digest the outputs; files not written since ``since`` are returned as stale."""
    out: dict[str, dict] = {}
    stale: list[str] = []
    for d in directories:
        if not d.exists():
            continue
        for f in sorted(d.rglob("*")):
            if not _is_output(f, manifest):
                continue
            key = str(f.relative_to(d.parent))
            if since is not None and f.stat().st_mtime < since:
                stale.append(key)
                continue
            out[key] = {"bytes": f.stat().st_size, "sha256": digest(f)}
    return out, stale


def combined(entries: dict[str, dict]) -> str:
    h = hashlib.sha256()
    for key in sorted(entries):
        h.update(f"{key}\0{entries[key]['sha256']}\n".encode())
    return h.hexdigest()


def clean_outputs(directories: list[Path]) -> None:
    """Remove the visible contents of each output directory.

    Hidden files (the marker, a tracked .gitkeep) are kept. A directory that
    has visible contents but no marker was not created by this runner, so it
    is left untouched and the run is refused.
    """
    for d in directories:
        if d.exists():
            visible = [e for e in d.iterdir() if not e.name.startswith(".")]
            if visible and not (d / MARKER).exists():
                raise SystemExit(
                    f"refusing to empty {d}: it was not created by run_pipeline.py "
                    f"(no {MARKER} marker) and is not empty. Remove it yourself or "
                    "choose another output directory."
                )
            for e in visible:
                shutil.rmtree(e) if e.is_dir() else e.unlink()
        else:
            d.mkdir(parents=True)
        (d / MARKER).write_text("Output directory managed by scripts/run_pipeline.py; "
                                "--clean empties it.\n")


def clear_failure_notes(directories: list[Path]) -> None:
    for d in directories:
        (d / FAILED_NAME).unlink(missing_ok=True)


def write_failure_notes(directories: list[Path], stage: str) -> None:
    """Say, next to any leftover files, that this run failed and they are not its results."""
    text = (f"The run failed at {stage}. Files in this directory were NOT produced by a complete\n"
            "run of the current inputs: any that are here are left over from an earlier run or\n"
            "from stages that finished before the failure. Do not use them as results.\n"
            "Fix the problem and run scripts/run_pipeline.py --clean.\n")
    for d in directories:
        d.mkdir(parents=True, exist_ok=True)
        (d / FAILED_NAME).write_text(text)


def mark(directories: list[Path]) -> None:
    for d in directories:
        if not d.exists():
            d.mkdir(parents=True)
            (d / MARKER).write_text("Output directory managed by scripts/run_pipeline.py; "
                                    "--clean empties it.\n")


def environment() -> dict[str, str]:
    import matplotlib
    import numpy
    import pandas
    import scipy

    try:
        commit = subprocess.run(
            ["git", "-C", str(config.ROOT), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:
        commit = "unknown"
    return {
        "python": sys.version.split()[0],
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "git_commit": commit,
    }


def reference_provenance() -> dict[str, str] | None:
    from mefv.reference import read_provenance

    path = config.POPULATION_FREQUENCIES
    if not path.exists():
        return None
    prov = read_provenance(path)
    return {k: prov.get(k, "") for k in ("source", "version", "retrieved")}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="start", default=None, help="first stage number to run")
    ap.add_argument("--only", nargs="*", default=None, help="run only these stage numbers")
    ap.add_argument("--manifest", default=None, help="where to write the run manifest")
    ap.add_argument("--clean", action="store_true",
                    help="empty the output directories first (full runs only)")
    ap.add_argument("--accept-reviewed-findings", default=None, metavar="RULES",
                    help="comma-separated blocking validation rules (V5, V6a, V6b, V7) the data "
                         "owner has seen; passed to stage 01 as MEFV_ACCEPT_REVIEWED_FINDINGS. "
                         "This does NOT release the run")
    ap.add_argument("--owner-disposition", default=None, metavar="FILE",
                    help="private JSON file recording the owner's explicit disposition of each "
                         "blocking finding, stored outside the repository; passed to stage 01 "
                         "as MEFV_OWNER_DISPOSITION")
    args = ap.parse_args()
    if args.accept_reviewed_findings is not None:
        os.environ["MEFV_ACCEPT_REVIEWED_FINDINGS"] = args.accept_reviewed_findings
    if args.owner_disposition is not None:
        os.environ["MEFV_OWNER_DISPOSITION"] = args.owner_disposition

    stages = STAGES
    partial = bool(args.only or args.start)
    if args.only:
        stages = [s for s in stages if s.name[:2] in set(args.only)]
    elif args.start:
        stages = [s for s in stages if s.name[:2] >= args.start]
    if args.clean and partial:
        ap.error("--clean empties every output directory, so it needs a full run")

    manifest_path = Path(args.manifest) if args.manifest else config.RUN_MANIFEST
    outputs = [*artifact_dirs(), config.PROCESSED_DIR]
    if any(manifest_path.resolve().is_relative_to(d.resolve()) for d in outputs):
        ap.error(f"the manifest {manifest_path} must not be inside an output directory")
    if args.clean:
        clean_outputs(outputs)
    else:
        mark(outputs)

    # Only a full run regenerates every output directory, so only it may withdraw an earlier
    # failure marker. A partial run does not show the other directories to be current.
    if not partial:
        clear_failure_notes(outputs)
    config.pin_build_timestamp()

    started = time.time()
    results = []
    for stage in stages:
        print(f"--- {stage.name}", flush=True)
        t0 = time.monotonic()
        proc = subprocess.run([sys.executable, str(stage)], env=os.environ.copy())
        results.append({"stage": stage.name, "returncode": proc.returncode,
                        "seconds": round(time.monotonic() - t0, 2)})
        if proc.returncode != 0:
            print(f"stage {stage.name} failed with code {proc.returncode}", file=sys.stderr)
            write_failure_notes(outputs, stage.name)
            break

    # Two seconds of slack absorb coarse file-system timestamps.
    since = None if partial else started - 2.0
    artifacts, stale = collect(artifact_dirs(), manifest_path, since)
    intermediates = {}
    for f in intermediate_files():
        if f.exists() and (since is None or f.stat().st_mtime >= since):
            intermediates[f.name] = {"bytes": f.stat().st_size, "sha256": digest(f)}

    manifest = {
        "manifest_version": MANIFEST_VERSION,
        "deterministic": {
            "artifacts": artifacts,
            "artifacts_sha256": combined(artifacts),
            "intermediates": intermediates,
            "settings": {
                "seed_bootstrap": config.SEED_BOOTSTRAP,
                "seed_figures": config.SEED_FIGURES,
                "bootstrap_replicates": config.BOOTSTRAP_REPLICATES,
                "source_date_epoch": os.environ.get("SOURCE_DATE_EPOCH"),
            },
            "stages": [{"stage": r["stage"], "returncode": r["returncode"]} for r in results],
            # true only for a full run in which every stage succeeded
            "run_complete": (not partial) and bool(results)
                            and all(r["returncode"] == 0 for r in results),
        },
        "execution": {
            "note": "Non-deterministic or machine-specific; excluded from "
                    "reproducibility comparisons.",
            "generated_at": datetime.now(UTC).isoformat(timespec="seconds"),
            "partial_run": partial,
            "stage_seconds": {r["stage"]: r["seconds"] for r in results},
            "environment": environment(),
            "paths": {
                "patients": str(config.PATIENTS),
                "variants": str(config.VARIANTS),
                "dictionary": str(config.DICTIONARY),
                "annotation_overlay": str(config.ANNOTATION_OVERLAY),
                "population_frequencies": str(config.POPULATION_FREQUENCIES),
                "tables": str(config.TABLES),
                "figures": str(config.FIGURES),
                "classification": str(config.CLASSIFICATION),
            },
            "reference_data": reference_provenance(),
            "stale_files_not_collected": stale,
        },
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"\nmanifest -> {manifest_path} ({len(artifacts)} artifacts, "
          f"{len(intermediates)} intermediates)")
    if stale:
        print(f"warning: {len(stale)} file(s) in the output directories were not written "
              "by this run and were not collected; use --clean for a fresh run.",
              file=sys.stderr)

    failed = [r for r in results if r["returncode"] != 0]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
