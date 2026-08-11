"""
test_all_features.py
A fast smoke test of an installed pathview-plus.

This is *not* the test suite. The real suite lives in ``tests/`` at the
repository root: 320 tests covering every module, every fixed bug and every
parity claim. Run it with ``pytest`` from a checkout.

What this file is for is the question "is my installation working?", answerable
from anywhere without a checkout:

    python -m pathview.test_all_features

It exercises the offline paths only — species resolution, bundled tables,
identifier routing, parsing, colour mapping and vector rendering — so it works
behind a firewall and finishes in a couple of seconds.

A note on the name
------------------
pathview-plus 2.x shipped a file with this name whose checks could not fail:
``pathview()`` returned ``{}`` on every error path and the assertions were
``assert result is not None``, so the suite passed green while nothing worked
(bug 14). Each check below asserts a *specific value*, and reports which ones
failed rather than stopping at the first.
"""

from __future__ import annotations

import sys
import traceback
from collections.abc import Callable


class SmokeTestFailure(AssertionError):
    """Raised when a smoke check does not hold."""


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def check_import() -> str:
    import pathview

    assert pathview.__version__, "no version string"
    missing = [n for n in pathview.__all__ if not hasattr(pathview, n)]
    assert not missing, f"__all__ names that do not exist: {missing}"
    return f"pathview {pathview.__version__}, {len(pathview.__all__)} exports"


def check_bundled_data() -> str:
    from pathview import bundled_files

    files = bundled_files()
    expected = {"korg.tsv.gz", "cpd_names.tsv.gz", "cpd_xref.tsv.gz",
                "edge_subtypes.tsv.gz", "bods.tsv.gz", "demo_gse16873.tsv.gz",
                "sbgn_index.tsv.gz", "sbgn_xref.tsv.gz"}
    missing = expected - set(files)
    assert not missing, f"bundled data missing from the install: {sorted(missing)}"
    return f"{len(files)} data files, {sum(files.values()) // 1024:,} KB"


def check_species() -> str:
    from pathview import get_species_code, organism_count

    n = organism_count()
    assert n > 8000, f"organism table looks truncated ({n} rows)"
    for query, expected in (("human", "hsa"), ("mouse", "mmu"),
                            ("9606", "hsa"), ("E. coli", "eco")):
        got = get_species_code(query).kegg_code
        assert got == expected, f"{query!r} resolved to {got}, expected {expected}"
    return f"{n:,} organisms, 4/4 lookups correct"


def check_compounds() -> str:
    from pathview import compound_name, cpd_name_to_kegg

    assert compound_name("C00022") == "Pyruvate"
    assert cpd_name_to_kegg(["Citrate"])["KEGG"][0] == "C00158"
    return "names and conjugate-base lookup correct"


def check_colour_scale() -> str:
    from pathview import colorpanel2, gene_scale

    ramp = colorpanel2(10, "#00FF00", "#BEBEBE", "#FF0000")
    assert ramp[0] == "#00FF00" and ramp[-1] == "#FF0000"

    scale = gene_scale(limit=2.0, bins=10)
    colours = scale.map_values([-3.0, 0.0, 3.0])
    assert colours[0] != colours[2], "up and down must differ"
    assert colours[0].upper() == "#00FF00" and colours[2].upper() == "#FF0000"
    assert scale.map_values([float("nan")])[0] == scale.na_col
    return "R-parity binning and clamping correct"


def check_aggregation() -> str:
    import polars as pl

    from pathview import mol_sum

    data = pl.DataFrame({"id": ["1", "2", "3"], "v": [1.0, -5.0, 3.0]})
    idmap = pl.DataFrame({"id": ["1", "2", "3"], "t": ["A", "A", "B"]})
    methods = ["sum", "mean", "median", "max", "min", "max_abs", "random", "first"]
    for method in methods:
        out = mol_sum(data, idmap, sum_method=method, rand_seed=1)
        assert out.height == 2, f"{method} produced {out.height} rows"
    assert mol_sum(data, idmap, sum_method="max_abs")["v"][0] == -5.0
    return f"all {len(methods)} aggregation methods work"


def check_identifier_routing() -> str:
    from pathview import id_route, map_ids_to_sbgn

    route = id_route("ENTREZ", "SYMBOL")
    assert route is not None, "no route from ENTREZ to SYMBOL"
    got = map_ids_to_sbgn(["1017"], "ENTREZ", "SYMBOL")["SYMBOL"][0]
    assert got == "CDK2", f"1017 mapped to {got}, expected CDK2"
    return f"routing works: {' -> '.join(route)}"


def check_sbgn_collection() -> str:
    from pathview import find_sbgn_pathway, sbgn_collection_info

    info = sbgn_collection_info()
    assert info["total"] > 5000, f"collection index looks truncated ({info['total']})"
    for pid, source in (("P00001", "panther"), ("SMP00001", "smpdb"),
                        ("GLYCOLYSIS", "metacyc")):
        assert find_sbgn_pathway(pid)["source"] == source
    return f"{info['total']:,} pathways across {len(info['by_source'])} sources"


def check_databases() -> str:
    from pathview import DATABASE_INFO, detect_database

    for pid, expected in (("hsa04110", "kegg"), ("R-HSA-109688", "reactome"),
                          ("P00001", "panther"), ("SMP00001", "smpdb"),
                          ("GLYCOLYSIS", "metacyc")):
        got = detect_database(pid)
        assert got == expected, f"{pid} detected as {got}, expected {expected}"
    stubs = [k for k, v in DATABASE_INFO.items()
             if v["available"] and not callable(v["downloader"])]
    assert not stubs, f"sources claiming availability with no downloader: {stubs}"
    return f"{len(DATABASE_INFO)} sources, all with working downloaders"


def check_splines() -> str:
    import numpy as np

    from pathview import catmull_rom_spline, route_edge_spline

    for points in ([(0, 0), (1, 2), (3, 1), (4, 3)],
                   [(0, 0), (0, 0), (1, 1)],
                   [(0, 0), (1, 1)]):
        curve = catmull_rom_spline(points, n_points=8)
        assert curve.size and np.isfinite(curve).all(), f"NaN for {points}"
    assert np.isfinite(route_edge_spline((0, 0), (10, 10))).all()
    return "curves finite for all inputs"


def check_parsing_and_render() -> str:
    """Parse a synthesised KGML and render it, with no network and no fixtures."""
    import tempfile
    from pathlib import Path

    import polars as pl

    from pathview import gene_scale, node_color, node_info, parse_kgml, pathway_edges

    kgml = """<?xml version="1.0"?>
<pathway name="path:tst00001" org="tst" number="00001" title="Smoke test">
  <entry id="1" name="tst:1017" type="gene">
    <graphics name="CDK2" x="100" y="100" width="46" height="17" type="rectangle"/>
  </entry>
  <entry id="2" name="tst:7157" type="gene">
    <graphics name="TP53" x="220" y="100" width="46" height="17" type="rectangle"/>
  </entry>
  <entry id="3" name="cpd:C00022" type="compound">
    <graphics name="C00022" x="160" y="180" width="8" height="8" type="circle"/>
  </entry>
  <relation entry1="1" entry2="2" type="PPrel">
    <subtype name="activation" value="--&gt;"/>
  </relation>
</pathway>
"""
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "tst00001.xml"
        path.write_text(kgml)

        pathway = parse_kgml(path)
        assert len(pathway) == 3, f"parsed {len(pathway)} entries, expected 3"

        nodes = node_info(pathway)
        edges = pathway_edges(pathway)
        assert edges.height == 1, f"parsed {edges.height} edges, expected 1"
        assert edges["subtype"][0] == "activation"

        values = pl.DataFrame({"id": ["1", "2"], "logfc": [2.0, -2.0]})
        colours = node_color(values, gene_scale(limit=2.0), id_col="id")
        assert colours["logfc_col"][0] != colours["logfc_col"][1]

        import matplotlib
        matplotlib.use("Agg", force=True)
        from pathview import keggview_vector

        out = keggview_vector(
            nodes, edges,
            color_map={"1": [colours["logfc_col"][0]],
                       "2": [colours["logfc_col"][1]]},
            pathway_name="tst00001", out_dir=tmp, output_format="png",
            plot_col_key=False,
        )
        assert out.exists() and out.stat().st_size > 1000, "render produced nothing"
        size_kb = out.stat().st_size // 1024
    return f"parsed 3 nodes / 1 edge, rendered {size_kb} KB"


CHECKS: list[tuple[str, Callable[[], str]]] = [
    ("import", check_import),
    ("bundled data", check_bundled_data),
    ("species lookup", check_species),
    ("compound names", check_compounds),
    ("colour scales", check_colour_scale),
    ("aggregation", check_aggregation),
    ("identifier routing", check_identifier_routing),
    ("SBGN collection", check_sbgn_collection),
    ("database dispatch", check_databases),
    ("splines", check_splines),
    ("parse + render", check_parsing_and_render),
]


def run_all(verbose: bool = True) -> tuple[int, int]:
    """
    Run every check.

    Returns ``(passed, failed)``.  Every check runs even if an earlier one
    fails, so one broken area does not hide the state of the rest.
    """
    passed = failed = 0
    width = max(len(name) for name, _ in CHECKS)

    for name, fn in CHECKS:
        try:
            detail = fn()
            passed += 1
            if verbose:
                print(f"  PASS  {name:<{width}}  {detail}")
        except Exception as exc:
            failed += 1
            if verbose:
                print(f"  FAIL  {name:<{width}}  {type(exc).__name__}: {exc}")
                if not isinstance(exc, AssertionError):
                    traceback.print_exc(limit=3)
    return passed, failed


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Smoke-test a pathview-plus install.")
    parser.add_argument("--quiet", "-q", action="store_true")
    args = parser.parse_args(argv)

    if not args.quiet:
        print("pathview-plus smoke test (offline)\n")

    passed, failed = run_all(verbose=not args.quiet)

    print(f"\n{passed} passed, {failed} failed")
    if failed:
        print("\nThis is a smoke test, not the test suite. For the full suite:\n"
              "    git clone https://github.com/raw-lab/pathview-plus\n"
              "    cd pathview-plus && pip install -e '.[dev]' && pytest")
    return 1 if failed else 0


if __name__ == "__main__":                                    # pragma: no cover
    sys.exit(main())
