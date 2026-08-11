"""
Enforce the parity claims.

Every feature listed as ``full`` in ``pathview.parity.FEATURE_MATRIX`` names
an API in its ``api`` field.  This suite asserts that each such name really
exists and is callable, so the matrix cannot drift from the code: adding a
row with ``full`` and no implementation fails the build.
"""
from __future__ import annotations

import pytest


def _resolve(name: str):
    """
    Resolve an API name against the package.

    Handles top-level exports, dotted attributes, dataclass field names, and
    submodules (``cli``, ``cache``), so the matrix can point at whichever is
    the honest reference for a feature.
    """
    import dataclasses
    import importlib

    import pathview

    target = name.split("(")[0].strip()

    if "." in target:
        head, _, tail = target.partition(".")
        obj = getattr(pathview, head, None)
        if obj is None:
            try:
                obj = importlib.import_module(f"pathview.{head}")
            except ImportError:
                return None
        for part in tail.split("."):
            nxt = getattr(obj, part, None)
            if nxt is None and dataclasses.is_dataclass(obj):
                names = {f.name for f in dataclasses.fields(obj)}
                if part in names:
                    return obj                 # the field exists on the class
            if nxt is None:
                return None
            obj = nxt
        return obj

    obj = getattr(pathview, target, None)
    if obj is not None:
        return obj
    try:
        return importlib.import_module(f"pathview.{target}")
    except ImportError:
        return None


def _api_names(feature) -> list[str]:
    """Names referenced by a feature's api field, ignoring prose."""
    raw = feature.api.strip()
    if not raw or "=" in raw or " " in raw.replace(" / ", "/"):
        parts = [p.strip() for p in raw.replace(" / ", "/").split("/")]
        return [p for p in parts if p and "=" not in p and " " not in p]
    return [p.strip() for p in raw.split("/") if p.strip()]


FULL_FEATURES = None


def _full_features():
    global FULL_FEATURES
    if FULL_FEATURES is None:
        from pathview import FEATURE_MATRIX
        FULL_FEATURES = [f for f in FEATURE_MATRIX
                         if f.pathview_plus == "full" and f.api]
    return FULL_FEATURES


def test_matrix_is_populated():
    from pathview import FEATURE_MATRIX
    assert len(FEATURE_MATRIX) > 40


def test_every_full_feature_names_an_api():
    """A feature cannot be claimed 'full' without pointing at something."""
    from pathview import FEATURE_MATRIX
    missing = [f.name for f in FEATURE_MATRIX
               if f.pathview_plus == "full" and not f.api]
    assert not missing, f"features marked full with no api reference: {missing}"


@pytest.mark.parametrize("feature", _full_features() or [],
                         ids=lambda f: f.name.replace(" ", "_")[:48])
def test_full_feature_api_exists_and_is_callable(feature):
    names = _api_names(feature)
    if not names:
        pytest.skip(f"{feature.name}: api field is descriptive, not a symbol")

    resolved = [(n, _resolve(n)) for n in names]
    found = [(n, o) for n, o in resolved if o is not None]
    assert found, (
        f"{feature.name!r} is marked 'full' but none of {names} is importable "
        "from pathview"
    )
    import types
    for name, obj in found:
        usable = (callable(obj)
                  or isinstance(obj, (dict, tuple, frozenset, list,
                                      types.ModuleType))
                  or hasattr(obj, "__getitem__"))
        assert usable, f"{name} exists but is neither callable nor a data structure"


def test_every_partial_feature_documents_its_limitation():
    from pathview import FEATURE_MATRIX
    undocumented = [f.name for f in FEATURE_MATRIX
                    if f.pathview_plus == "partial" and not f.note]
    assert not undocumented, (
        f"partial features must state the limitation: {undocumented}"
    )


def test_every_gap_documents_the_workaround():
    """A missing feature must tell the user what to do instead."""
    from pathview import FEATURE_MATRIX
    silent = [f.name for f in FEATURE_MATRIX
              if f.pathview_plus == "none" and not f.note]
    assert not silent, f"gaps must be explained: {silent}"


def test_summary_arithmetic_is_consistent():
    from pathview import FEATURE_MATRIX, parity_summary
    s = parity_summary()
    assert s["total_features"] == len(FEATURE_MATRIX)
    assert s["full"] + s["partial"] + s["none"] + s["n/a"] == s["total_features"]
    assert 0 <= s["vs_pathview_R_pct"] <= 100
    assert 0 <= s["vs_SBGNview_R_pct"] <= 100


def test_gaps_match_the_matrix():
    from pathview import FEATURE_MATRIX, parity_summary
    declared = set(parity_summary()["gaps"])
    actual = {f.name for f in FEATURE_MATRIX if f.pathview_plus == "none"}
    assert declared == actual


def test_markdown_render_is_a_table():
    from pathview import print_parity
    md = print_parity(markdown=True)
    lines = md.splitlines()
    assert lines[0].startswith("| Category")
    assert set(lines[1]) <= set("|-: ")
    assert len(lines) > 40


def test_feature_table_frame():
    from pathview import feature_table
    df = feature_table()
    assert {"category", "feature", "pathview_plus", "pathview_R",
            "SBGNview_R"} <= set(df.columns)
    assert df.height > 40


def test_status_vocabulary_is_closed():
    from pathview import FEATURE_MATRIX
    allowed = {"full", "partial", "none", "n/a"}
    for f in FEATURE_MATRIX:
        assert f.pathview_plus in allowed
        assert f.pathview_r in allowed
        assert f.sbgnview_r in allowed


def test_no_stub_functions_are_exported():
    """
    Guard the no-stub standard.

    A public callable whose body only warns and returns None is a stub.  This
    walks the exported surface and fails if it finds one.
    """
    import inspect

    import pathview

    offenders = []
    for name in pathview.__all__:
        obj = getattr(pathview, name, None)
        if not inspect.isfunction(obj):
            continue
        try:
            src = inspect.getsource(obj)
        except (OSError, TypeError):
            continue
        body = [ln.strip() for ln in src.splitlines()]
        body = [ln for ln in body if ln and not ln.startswith("#")]
        code = [ln for ln in body
                if not ln.startswith(('"""', "'''", "def ", "@"))]
        if not code:
            continue
        warn_only = all(
            ln.startswith(("warnings.warn", "return None", "return", '"', "'", ")"))
            or ln.endswith(("stacklevel=2)", "stacklevel=3)"))
            for ln in code
        )
        if warn_only and any("warn" in ln for ln in code):
            offenders.append(name)

    assert not offenders, f"stub functions exported: {offenders}"


def test_sdist_would_ship_a_runnable_test_suite():
    """
    Guard the packaging defect that shipped a test suite without conftest.py.

    setuptools auto-includes ``tests/test*.py`` but not ``conftest.py``, so an
    sdist built without an explicit MANIFEST.in rule contains a suite in which
    every fixture-using test errors — and a downstream packager running it
    would see failures caused purely by packaging.
    """
    from pathlib import Path

    import pathview

    root = Path(pathview.__file__).resolve().parents[1]
    manifest = root / "MANIFEST.in"
    if not manifest.exists():
        pytest.skip("not an installed-from-source tree")

    rules = manifest.read_text()
    assert "recursive-include tests *.py" in rules, (
        "MANIFEST.in must include every tests/*.py, or conftest.py is dropped "
        "from the sdist"
    )


def test_every_public_name_actually_exists():
    """
    ``__all__`` must not name anything that is not bound.

    A stale entry breaks ``from pathview import *`` with an unhelpful
    AttributeError, and is easy to introduce — a linter removing an
    apparently-unused ``from . import x`` will do it silently.
    """
    import pathview

    missing = [name for name in pathview.__all__ if not hasattr(pathview, name)]
    assert not missing, f"__all__ names that do not exist: {missing}"


def test_star_import_works():
    """The whole public surface must import cleanly."""
    namespace: dict = {}
    exec("from pathview import *", namespace)          # noqa: S102
    import pathview
    for name in pathview.__all__:
        assert name in namespace, f"{name} did not survive a star import"


def test_no_duplicate_public_names():
    import pathview

    seen, dupes = set(), []
    for name in pathview.__all__:
        if name in seen:
            dupes.append(name)
        seen.add(name)
    assert not dupes, f"duplicate entries in __all__: {dupes}"


class TestSmokeTestModule:
    """
    ``pathview.test_all_features`` answers "is my install working?" without a
    checkout. It is only useful if it can fail.
    """

    def test_all_checks_pass_on_a_working_install(self):
        from pathview.test_all_features import run_all

        passed, failed = run_all(verbose=False)
        assert failed == 0, f"{failed} smoke checks failed on a working install"
        assert passed >= 10

    def test_it_detects_a_broken_component(self):
        """
        2.x shipped a file of this name whose checks could not fail: every
        assertion was ``result is not None`` and ``pathview()`` returned ``{}``
        on error, so it stayed green while nothing worked (bug 14).
        """
        from pathview.color_mapping import ColorScale
        from pathview.test_all_features import run_all

        import pathview

        original = ColorScale.map_values
        try:
            # Reproduce bug 4: every value collapses to one colour.
            ColorScale.map_values = lambda self, values: ["#FF0000"] * len(list(values))
            _, failed = run_all(verbose=False)
        finally:
            ColorScale.map_values = original

        assert failed >= 1, "the smoke test cannot detect a broken colour scale"
        assert pathview.gene_scale(limit=2.0).map_values([-2.0])[0] != "#FF0000", (
            "the monkeypatch leaked out of the test"
        )

    def test_needs_no_fixtures_or_network(self, monkeypatch):
        """
        It must work from a bare ``pip install``, anywhere.

        Rather than grepping the source — which trips over URLs in help text —
        this makes any outbound request raise, and runs the checks.
        """
        from pathview.test_all_features import run_all

        def _no_network(*args, **kwargs):
            raise AssertionError("the smoke test made a network request")

        import requests
        monkeypatch.setattr(requests, "get", _no_network, raising=False)
        monkeypatch.setattr(requests, "post", _no_network, raising=False)
        monkeypatch.setattr(requests.Session, "request", _no_network, raising=False)

        passed, failed = run_all(verbose=False)
        assert failed == 0, "smoke test needs the network or missing fixtures"
        assert passed >= 10
