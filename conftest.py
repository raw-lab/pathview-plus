"""
conftest.py (repository root)

Makes ``lib/`` importable as ``pathview`` when the package is not installed.

The source directory is ``lib/`` and setuptools maps it to the ``pathview``
package name at install time (``package-dir = { pathview = "lib" }``).  That
mapping only exists once the package is installed, so running ``pytest``
straight from a fresh clone would otherwise fail with ImportError before any
test ran.  Installing first (``pip install -e ".[dev]"``) is still the
recommended workflow; this only removes a confusing failure for anyone who
does not.

If ``pathview`` is already importable — an editable or regular install — this
does nothing at all, so the installed package always wins.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _install_lib_as_pathview() -> None:
    if importlib.util.find_spec("pathview") is not None:
        return                                  # already installed; leave it alone

    lib = Path(__file__).resolve().parent / "lib"
    init = lib / "__init__.py"
    if not init.exists():
        return

    spec = importlib.util.spec_from_file_location(
        "pathview", init, submodule_search_locations=[str(lib)]
    )
    if spec is None or spec.loader is None:      # pragma: no cover
        return
    module = importlib.util.module_from_spec(spec)
    sys.modules["pathview"] = module
    spec.loader.exec_module(module)


_install_lib_as_pathview()
