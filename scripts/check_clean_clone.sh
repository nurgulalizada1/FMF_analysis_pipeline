#!/usr/bin/env bash
# Verify that the committed state of a branch works from a fresh clone.
#
# Clones the given branch (default: the current one) into a temporary
# directory, so no untracked or private file from this checkout is present,
# installs the pinned dependencies for the given Python into a new virtual
# environment, and runs lint, the test suite and the demonstration twice.
#
# Usage: scripts/check_clean_clone.sh [python-executable] [branch]
set -euo pipefail
PY="${1:-python3}"
BRANCH="${2:-$(git rev-parse --abbrev-ref HEAD)}"
SRC="$(git rev-parse --show-toplevel)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

VERSION="$("$PY" -c 'import sys; print(f"{sys.version_info[0]}.{sys.version_info[1]}")')"
LOCK="requirements-lock-py${VERSION}.txt"

git clone --quiet --branch "$BRANCH" "$SRC" "$WORK/repo"
cd "$WORK/repo"
echo "clean clone of $BRANCH at $(git rev-parse HEAD) with Python $VERSION"
test -f "$LOCK" || { echo "no lock file $LOCK for Python $VERSION" >&2; exit 2; }
for private in study_parameters.json data/external/source_annotations.csv \
               data/reviewed/variant_annotations_reviewed.csv; do
    test ! -e "$private" || { echo "private file present in clone: $private" >&2; exit 3; }
done

"$PY" -m venv "$WORK/venv"
# shellcheck disable=SC1091
. "$WORK/venv/bin/activate"
python -m pip install --quiet --upgrade pip
python -m pip install --quiet -r "$LOCK"
python -m pip install --quiet --no-deps -e .
python -c "import mefv, mefv.config, mefv.reference; print('import ok:', mefv.__file__)"
ruff check src scripts tests
pytest -q
make demo-check
echo "clean-clone check passed"
