"""Sphinx configuration for the pathview-plus documentation."""

from __future__ import annotations

import importlib.util
import os
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# The package source is lib/, mapped to the `pathview` name at install time.
# Sphinx may run before that mapping exists, so load it in place if needed.
if importlib.util.find_spec("pathview") is None:
    _lib = Path(__file__).resolve().parents[1] / "lib"
    _spec = importlib.util.spec_from_file_location(
        "pathview", _lib / "__init__.py", submodule_search_locations=[str(_lib)]
    )
    _mod = importlib.util.module_from_spec(_spec)
    sys.modules["pathview"] = _mod
    _spec.loader.exec_module(_mod)

# Docs must build without network access, exactly like the test suite.
os.environ.setdefault("PATHVIEW_OFFLINE", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import pathview  # noqa: E402

project = "pathview-plus"
author = "Richard Allen White III — RAW Lab, UNC Charlotte"
copyright = f"{date.today().year}, RAW Lab"
release = pathview.__version__
version = ".".join(release.split(".")[:2])

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx_copybutton",
    "sphinx_design",
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

# -- autodoc ---------------------------------------------------------------
autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_rtype = False

# -- theme -----------------------------------------------------------------
html_theme = "furo"
html_static_path = ["_static"]
html_title = f"pathview-plus {version}"
html_theme_options = {
    "source_repository": "https://github.com/raw-lab/pathview-plus/",
    "source_branch": "main",
    "source_directory": "docs/",
    "light_css_variables": {
        "color-brand-primary": "#1565C0",
        "color-brand-content": "#1565C0",
    },
    "dark_css_variables": {
        "color-brand-primary": "#64B5F6",
        "color-brand-content": "#64B5F6",
    },
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "polars": ("https://docs.pola.rs/api/python/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "networkx": ("https://networkx.org/documentation/stable/", None),
}

myst_enable_extensions = ["colon_fence", "deflist", "substitution"]
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True


# -- generated pages -------------------------------------------------------

def _write_parity_page(app) -> None:
    """
    Generate the parity table from ``pathview.parity`` at build time.

    Keeping the docs table generated rather than hand-written means it cannot
    disagree with the matrix the test suite enforces.
    """
    from pathview import parity_summary, print_parity
    from pathview.parity import FEATURE_MATRIX

    s = parity_summary()
    out = Path(app.srcdir) / "parity.md"

    lines = [
        "# Feature parity",
        "",
        "```{note}",
        "This page is generated from `pathview.parity.FEATURE_MATRIX` when the",
        "docs are built, and `tests/test_parity.py` asserts that every feature",
        "marked **full** resolves to a real, importable API. The table cannot",
        "drift from the code.",
        "```",
        "",
        "## Summary",
        "",
        f"- **{s['total_features']}** capabilities tracked",
        f"- **{s['full']}** full, **{s['partial']}** partial, **{s['none']}** not implemented",
        f"- Covers **{s['vs_pathview_R_pct']}%** of what pathview (R) does "
        f"({s['vs_pathview_R']})",
        f"- Covers **{s['vs_SBGNview_R_pct']}%** of what SBGNview (R) does "
        f"({s['vs_SBGNview_R']})",
        f"- **{s['beyond_both']}** capabilities present in neither R package",
        "",
        "## Known gaps",
        "",
    ]
    for f in FEATURE_MATRIX:
        if f.pathview_plus == "none":
            lines.append(f"- **{f.name}** — {f.note}")
    lines += ["", "## Full matrix", "", print_parity(markdown=True), ""]
    out.write_text("\n".join(lines), encoding="utf-8")


def setup(app):
    app.connect("builder-inited", _write_parity_page)
    return {"parallel_read_safe": True}
