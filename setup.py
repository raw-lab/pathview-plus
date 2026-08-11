"""
setup.py
Compatibility shim for `python setup.py` and for tooling that expects it.

All packaging metadata — name, version, dependencies, entry points and the
lib/ -> pathview package mapping — lives in pyproject.toml. Nothing is
declared here, deliberately: this file carrying its own version string is
exactly how pathview-plus 2.x ended up with setup.py saying 2.0.2 while
pathview/__init__.py said 2.0.0 (bug 16). One source of truth avoids that.

Install with:

    pip install .            # or: pip install -e ".[dev]"
"""

from setuptools import setup

setup()
