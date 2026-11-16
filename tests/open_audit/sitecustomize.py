"""Test helper: record the files a Python process opens.

Placed on PYTHONPATH by the integration tests only. When MEFV_OPEN_AUDIT_LOG
is set, every file opened under MEFV_OPEN_AUDIT_ROOT is appended to that log,
so a test can show that no private default location was read. It changes
nothing else about the process.
"""
import os
import sys

_LOG = os.environ.get("MEFV_OPEN_AUDIT_LOG")
_ROOT = os.environ.get("MEFV_OPEN_AUDIT_ROOT")

if _LOG and _ROOT:
    _root = os.path.realpath(_ROOT)
    _seen: set[str] = set()
    _busy = False

    def _hook(event, args):
        global _busy
        if event != "open" or _busy or not args:
            return
        path = args[0]
        if isinstance(path, int):
            return
        try:
            path = os.path.realpath(os.fsdecode(path))
        except Exception:
            return
        if path in _seen or not path.startswith(_root + os.sep):
            return
        _seen.add(path)
        _busy = True
        try:
            with open(_LOG, "a", encoding="utf-8") as fh:
                fh.write(path + "\n")
        finally:
            _busy = False

    sys.addaudithook(_hook)
