#!/usr/bin/env python3
"""
bin/pathview-cli.py
Direct-run launcher, kept for continuity.

The CLI itself lives in ``pathview.cli`` and is installed as the
``pathview-plus`` and ``pathview-cli`` commands.  This script lets you run it
straight from a clone without installing anything.
"""

import sys

def _ensure_pathview_importable() -> None:
    """
    Make ``lib/`` importable as ``pathview`` when the package is not installed.

    The source lives in ``lib/`` and is mapped to the ``pathview`` package name
    by setuptools at install time.  Running this file directly from a clone
    happens before that mapping exists, so fall back to loading it in place.
    """
    import importlib.util
    from pathlib import Path

    if importlib.util.find_spec("pathview") is not None:
        return

    lib = Path(__file__).resolve().parents[1] / "lib"
    init = lib / "__init__.py"
    if not init.exists():
        return
    spec = importlib.util.spec_from_file_location(
        "pathview", init, submodule_search_locations=[str(lib)]
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["pathview"] = module
    spec.loader.exec_module(module)


_ensure_pathview_importable()

from pathview.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
