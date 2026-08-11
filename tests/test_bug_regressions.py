"""
One test per bug fixed in 3.0.0.

Each test fails against 2.x.  Test ids match the numbering in the README's
bug checklist so a reviewer can trace claim -> test -> fix.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

# --- Bug 1 & 2: species lookup ---------------------------------------------

def test_bug01_species_lookup_works_offline():
    """2.x: RuntimeError('Failed to fetch KEGG organism list: 403') with no network."""
    from pathview import get_species_code, is_offline
    assert is_offline()
    assert get_species_code("hsa").kegg_code == "hsa"


def test_bug02_species_lookup_accepts_names_not_just_codes():
    """2.x matched whole space-split fields, so only the bare code ever matched."""
    from pathview import get_species_code
    for query in ("human", "Homo sapiens", "9606", "T01001", "hsa"):
        assert get_species_code(query).kegg_code == "hsa", query


def test_bug02b_common_name_prefers_the_common_organism():
    """'mouse' must be Mus musculus, not Microcebus murinus (mouse lemur)."""
    from pathview import get_species_code
    assert get_species_code("mouse").kegg_code == "mmu"
    assert get_species_code("rat").kegg_code == "rno"


def test_bug02c_organism_table_is_complete():
    from pathview import organism_count
    assert organism_count() > 8000


# --- Bug 3 & 4: bgcolor treated as data ------------------------------------

def test_bug03_bgcolor_is_metadata_not_a_value_column():
    from pathview.constants import NODE_META_COLS
    assert "bgcolor" in NODE_META_COLS
    assert "fgcolor" in NODE_META_COLS


def test_bug04_every_node_does_not_render_red():
    """
    2.x: bgcolor '#FFFFFF' was treated as data, int('FFFFFF',16)=16777215
    clipped to the scale maximum, so every node came out solid red whatever
    the user's data said.
    """
    from pathview import gene_scale, node_color

    df = pl.DataFrame({
        "id": ["a", "b", "c"],
        "bgcolor": ["#FFFFFF", "#BFFFBF", "#FFFFFF"],
        "label": ["A", "B", "C"],
        "exp1": [-2.0, 0.0, 2.0],
    })
    out = node_color(df, gene_scale(limit=2.0, bins=10), id_col="id")

    assert out.columns == ["id", "exp1_col"], "only numeric columns may be mapped"
    colors = out["exp1_col"].to_list()
    assert len(set(colors)) == 3, f"data must produce distinct colours, got {colors}"
    assert colors[0] != colors[2], "down- and up-regulated must differ"
    assert colors[2].upper() == "#FF0000" and colors[0].upper() == "#00FF00"


# --- Bug 5: aggregation methods --------------------------------------------

@pytest.mark.parametrize("method", ["sum", "mean", "median", "max", "min",
                                    "max_abs", "random", "first"])
def test_bug05_all_documented_sum_methods_work(method):
    """2.x: max_abs and random raised AttributeError inside group_by.agg."""
    from pathview import mol_sum

    data = pl.DataFrame({"id": ["1", "2", "3"], "v": [1.0, -5.0, 3.0]})
    idmap = pl.DataFrame({"id": ["1", "2", "3"], "t": ["A", "A", "B"]})
    out = mol_sum(data, idmap, sum_method=method)
    assert out.height == 2
    assert out["v"].null_count() == 0


def test_bug05b_max_abs_keeps_the_sign():
    from pathview import mol_sum
    data = pl.DataFrame({"id": ["1", "2"], "v": [1.0, -5.0]})
    idmap = pl.DataFrame({"id": ["1", "2"], "t": ["A", "A"]})
    assert mol_sum(data, idmap, sum_method="max_abs")["v"][0] == -5.0


def test_bug05c_random_is_reproducible_with_a_seed():
    from pathview import mol_sum
    data = pl.DataFrame({"id": [str(i) for i in range(20)],
                         "v": [float(i) for i in range(20)]})
    idmap = pl.DataFrame({"id": [str(i) for i in range(20)], "t": ["A"] * 20})
    a = mol_sum(data, idmap, sum_method="random", rand_seed=7)["v"][0]
    b = mol_sum(data, idmap, sum_method="random", rand_seed=7)["v"][0]
    assert a == b


# --- Bug 6: SBGN namespaces -------------------------------------------------

def test_bug06_namespaced_sbgn_parses_identically(fixtures_dir):
    """2.x returned 0 glyphs and 0 arcs for every namespaced (i.e. real) file."""
    from pathview import parse_sbgn

    bare = parse_sbgn(fixtures_dir / "P00001.new.layout.sbgn")
    ns = parse_sbgn(fixtures_dir / "P00001.namespaced.sbgn")

    assert len(ns.glyphs) > 0 and len(ns.arcs) > 0
    assert len(ns.glyphs) == len(bare.glyphs)
    assert len(ns.arcs) == len(bare.arcs)


def test_bug06b_sbgn_ports_resolve_to_glyphs(fixtures_dir):
    """Arcs reference <port> ids; unresolved ports leave the map edgeless."""
    from pathview import arc_resolution_report, parse_sbgn

    p = parse_sbgn(fixtures_dir / "ports_pd.sbgn")
    rep = arc_resolution_report(p)
    assert rep["resolution_rate"] == 1.0
    assert rep["arcs_via_port"] >= 2


def test_bug06c_auxiliary_glyphs_are_not_nodes(fixtures_dir):
    """State variables must not be promoted to top-level nodes."""
    from pathview import parse_sbgn
    p = parse_sbgn(fixtures_dir / "ports_pd.sbgn")
    assert "hk1_p" not in p.glyphs
    assert p.glyphs["hk1"].state_variables[0]["value"] == "P"
    assert p.glyphs["g6p"].clone_marker is True


# --- Bug 7: spline NaN ------------------------------------------------------

@pytest.mark.parametrize("pts", [
    [(0, 0), (1, 2), (3, 1), (4, 3)],
    [(0, 0), (0, 0), (1, 1), (2, 0)],
    [(0, 0), (1, 1)],
    [(5, 5)],
    [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)],
])
def test_bug07_catmull_rom_is_always_finite(pts):
    """2.x padded with duplicate endpoints -> divide-by-zero -> NaN."""
    from pathview import catmull_rom_spline
    curve = catmull_rom_spline(pts, n_points=10)
    assert curve.size > 0
    assert np.isfinite(curve).all()


def test_bug07b_spline_interpolates_its_control_points():
    from pathview import catmull_rom_spline
    pts = [(0.0, 0.0), (1.0, 2.0), (3.0, 1.0), (4.0, 3.0)]
    curve = catmull_rom_spline(pts, n_points=12)
    assert np.allclose(curve[0], pts[0], atol=1e-6)
    assert np.allclose(curve[-1], pts[-1], atol=1e-6)


# --- Bug 8: named colours ---------------------------------------------------

def test_bug08_named_colours_parse():
    """2.x: int('re', 16) -> ValueError, crashing highlight_nodes' own default."""
    from pathview import to_rgb
    assert to_rgb("red") == (255, 0, 0)
    assert to_rgb("#F00") == (255, 0, 0)
    assert to_rgb("rebeccapurple") == (102, 51, 153)


def test_bug08b_highlight_nodes_default_colour_does_not_crash(workdir, gene_data):
    from pathview import highlight_nodes, pathview

    res = pathview("04110", gene_data=gene_data, species="hsa",
                   kegg_dir=workdir, out_dir=workdir, render_mode="native",
                   quiet=True)
    out = res + highlight_nodes(["1017"])          # default color="red"
    assert out.image_array is not None


# --- Bugs 9 & 10: node geometry --------------------------------------------

def test_bug09_no_coordinate_flipping():
    """Compounds were painted at the mirrored y; geometry is now one rule."""
    from pathview import node_boxes
    df = pl.DataFrame({"entry_id": ["1"], "x": [100.0], "y": [50.0],
                       "width": [46.0], "height": [17.0], "type": ["gene"],
                       "shape": ["rectangle"], "label": ["X"], "bgcolor": ["#FFF"]})
    box = node_boxes(df)[0]
    assert box.y == 50.0 and box.top == 41.5 and box.bottom == 58.5


def test_bug10_compound_radius_is_half_the_diameter():
    """2.x used the full width as the radius, doubling every compound circle."""
    from pathview import node_boxes
    df = pl.DataFrame({"entry_id": ["1"], "x": [10.0], "y": [10.0],
                       "width": [8.0], "height": [8.0], "type": ["compound"],
                       "shape": ["circle"], "label": ["C"], "bgcolor": ["#FFF"]})
    assert node_boxes(df)[0].radius == 4.0


def test_bug10b_roundrectangle_is_not_a_circle():
    """KEGG pathway-link nodes are rounded rectangles, not discs."""
    from pathview import node_boxes
    df = pl.DataFrame({"entry_id": ["1"], "x": [10.0], "y": [10.0],
                       "width": [90.0], "height": [25.0], "type": ["map"],
                       "shape": ["roundrectangle"], "label": ["Glycolysis"],
                       "bgcolor": ["#FFF"]})
    assert node_boxes(df)[0].is_round is False


# --- Bug 11: graph view had no edges ---------------------------------------

def test_bug11_graph_view_has_edges(nodes_tca, edges_tca):
    """2.x called add_node in a loop and never added a single edge."""
    from pathview import build_graph
    G = build_graph(nodes_tca, edges_tca)
    assert G.number_of_edges() > 0
    assert G.number_of_edges() == edges_tca.height


def test_bug11b_kgml_yields_relation_and_reaction_edges(edges_tca):
    kinds = set(edges_tca["source_kind"].to_list())
    assert {"relation", "reaction"} <= kinds


# --- Bug 12: negative unmapped count ---------------------------------------

def test_bug12_unmapped_count_is_never_negative():
    """2.x computed mol_data.height - merged.height, which counts join expansion."""
    from pathview import mol_sum
    data = pl.DataFrame({"id": ["1", "2"], "v": [1.0, 2.0]})
    idmap = pl.DataFrame({"id": ["1", "1", "1", "2"],
                          "t": ["A", "B", "C", "D"]})     # one-to-many
    res = mol_sum(data, idmap, detailed=True)
    assert res.n_unmapped == 0
    assert res.n_mapped == 2
    assert res.n_targets == 4


# --- Bug 13: cpd-only runs --------------------------------------------------

def test_bug13_compound_data_alone_does_not_crash(workdir, cpd_data):
    """2.x's CLI cast gene_data unconditionally, so --cpd-data alone crashed."""
    from pathview import pathview
    res = pathview("00020", cpd_data=cpd_data, species="hsa",
                   kegg_dir=workdir, out_dir=workdir,
                   render_mode="vector", quiet=True)
    assert res.plot_data_cpd is not None
    assert res.plot_data_gene is None, "no gene data was supplied"
    assert res.plot_data_cpd.height > 0
    assert res.output_path.exists()


# --- Bug 14: tests that passed when everything failed ----------------------

def test_bug14_failure_is_distinguishable_from_success(workdir, gene_data):
    """
    2.x returned {} on every failure and the suite asserted 'is not None',
    so the tests passed even when nothing worked.
    """
    from pathview.errors import PathwayNotFoundError

    from pathview import PathwayResult, pathview

    good = pathview("00020", gene_data=gene_data, species="hsa",
                    kegg_dir=workdir, out_dir=workdir,
                    render_mode="vector", quiet=True)
    assert isinstance(good, PathwayResult) and bool(good) is True

    with pytest.raises(PathwayNotFoundError):
        pathview("99999", gene_data=gene_data, species="hsa",
                 kegg_dir=workdir, out_dir=workdir,
                 render_mode="vector", quiet=True)


def test_bug14b_unknown_species_raises_with_suggestions():
    from pathview.errors import SpeciesNotFoundError

    from pathview import get_species_code
    with pytest.raises(SpeciesNotFoundError) as exc:
        get_species_code("Homo sapein")
    assert "Homo sapiens" in str(exc.value) or "hsa" in str(exc.value)


# --- Bug 15: silent stubs ---------------------------------------------------

def test_bug15_no_silent_stubs_in_the_database_layer():
    """
    2.x exported ``download_panther`` and ``download_smpdb`` as
    warn-and-return-None: they advertised support and delivered nothing, and
    DATABASE_INFO listed them beside working downloaders so a caller could not
    tell which was which.

    3.0.0 makes them real, backed by the pre-generated SBGN collection, so the
    property guarded here is stronger than "the stubs are gone": every
    declared downloader must be a genuine callable that either produces a file
    or raises. A function that only warns and returns None fails this test
    whatever its name.
    """
    import inspect

    import pathview
    from pathview import DATABASE_INFO

    for key, info in DATABASE_INFO.items():
        assert isinstance(info["available"], bool), key
        assert info["note"], f"{key} must document its status"
        assert info["source"], f"{key} must name where its data comes from"

        fn = info["downloader"]
        if not info["available"]:
            assert fn is None, f"{key} is unavailable but declares a downloader"
            continue

        assert callable(fn), f"{key} is available but has no downloader"
        code = [
            line.strip() for line in inspect.getsource(fn).splitlines()
            if line.strip()
            and not line.strip().startswith(("#", '"""', "'''", "def ", "@"))
        ]
        stub = bool(code) and all(
            line.startswith(("warnings.warn", "return None", "return"))
            for line in code
        )
        assert not stub, f"{key}'s downloader is a stub"

    for name in ("download_kegg", "download_reactome", "download_panther",
                 "download_metacyc", "download_smpdb", "download_metacrop"):
        assert hasattr(pathview, name), f"{name} is not exported"
        assert name in pathview.__all__, f"{name} is missing from __all__"


def test_bug15a_capability_table_stays_self_describing():
    """
    ``available_databases()`` must agree with the table it summarises.

    The failure mode being guarded is a table entry that claims availability
    it does not have — which is how 2.x's stubs hid in plain sight.
    """
    from pathview import DATABASE_INFO, available_databases

    assert set(available_databases()) == {
        k for k, v in DATABASE_INFO.items() if v["available"]
    }


def test_bug15b_change_labels_and_opacity_actually_do_something(workdir, gene_data):
    from pathview import change_labels, highlight_nodes, pathview

    res = pathview("04110", gene_data=gene_data, species="hsa",
                   kegg_dir=workdir, out_dir=workdir,
                   render_mode="native", quiet=True)
    before = res.image_array.copy()

    relabelled = res + change_labels({"1017": "CDK2*"})
    assert not np.array_equal(before, relabelled.image_array)
    assert relabelled.label_changes["1017"] == "CDK2*"

    faint = res + highlight_nodes(["1017"], color="red", opacity=0.25)
    solid = res + highlight_nodes(["1017"], color="red", opacity=1.0)
    assert not np.array_equal(faint.image_array, solid.image_array)


def test_bug15c_list_reactome_pathways_respects_species():
    """2.x hard-coded Homo sapiens into the URL regardless of the argument."""
    from pathview import reactome_top_url

    assert reactome_top_url("Mus musculus").endswith("Mus%20musculus")
    assert reactome_top_url("Danio rerio").endswith("Danio%20rerio")
    assert "Homo" not in reactome_top_url("Mus musculus")


# --- Bug 16: version consistency -------------------------------------------

def test_bug16_version_is_consistent():
    from pathlib import Path

    import tomllib

    import pathview

    root = Path(pathview.__file__).resolve().parents[1]
    pyproject = root / "pyproject.toml"
    if pyproject.exists():
        meta = tomllib.loads(pyproject.read_text())
        assert meta["project"]["version"] == pathview.__version__


# --- Bug 17 (found during the rewrite): sim data must be mappable ----------

def test_bug17_simulated_ids_are_real_identifiers():
    """2.x emitted gene1..gene1000, which map to nothing at all."""
    from pathview import sim_mol_data

    genes = sim_mol_data("gene", n_mol=20)["entrez"].to_list()
    assert all(g.isdigit() for g in genes)

    cpds = sim_mol_data("cpd", n_mol=20)["compound"].to_list()
    assert all(c.startswith(("C", "G", "D")) for c in cpds)


# --- Bug 18 (found while generating the release figures): highlight geometry

def test_bug18_highlights_land_on_the_node_they_mark(workdir, gene_data):
    """
    Highlights are drawn in KGML coordinates, so the raster they are drawn on
    must be in KGML coordinates too.  Carrying the *composed figure* instead —
    which has a title band, colour keys, padding and a dpi scale — displaces
    every highlight, silently and plausibly.
    """
    from pathview import highlight_nodes, node_boxes, pathview

    res = pathview("04110", gene_data=gene_data, species="hsa",
                   kegg_dir=workdir, out_dir=workdir,
                   render_mode="native", quiet=True)

    target = "1017"
    # One gene id can occur in several KGML entries, so collect every box the
    # highlight is expected to touch.
    entries = {row["entry_id"] for row in res.node_data.iter_rows(named=True)
               if target in (row.get("kegg_names") or [])}
    assert entries, "fixture no longer contains the target gene"
    boxes = [b for b in node_boxes(res.node_data) if b.entry_id in entries]

    before = res.frame.array.copy()
    after = (res + highlight_nodes([target], color="#7C3AED", width=4)).frame

    changed = np.argwhere(np.any(after.array != before, axis=2))
    assert changed.size > 0, "highlight drew nothing"

    pad = 10
    regions = []
    for b in boxes:
        px0, py0 = after.to_pixels(b.left, b.top)
        px1, py1 = after.to_pixels(b.right, b.bottom)
        regions.append((px0 - pad, py0 - pad, px1 + pad, py1 + pad))

    stray = [
        (int(x), int(y)) for y, x in changed
        if not any(x0 <= x <= x1 and y0 <= y <= y1 for x0, y0, x1, y1 in regions)
    ]
    assert not stray, (
        f"{len(stray)} highlighted pixels fall outside every target node, "
        f"e.g. {stray[:3]}; expected regions {[tuple(round(v) for v in r) for r in regions]}"
    )


def test_bug18b_raster_frame_transform_is_consistent():
    from pathview import RasterFrame

    frame = RasterFrame(np.zeros((100, 200, 3), np.uint8), x0=50.0, y0=20.0, scale=2.0)
    assert frame.to_pixels(50.0, 20.0) == (0.0, 0.0)
    assert frame.to_pixels(150.0, 70.0) == (200.0, 100.0)
    assert frame.length(10.0) == 20.0


@pytest.mark.parametrize("mode", ["native", "vector", "svg", "graph"])
def test_bug18c_highlighting_works_in_every_mode(workdir, gene_data, mode):
    """Every mode must carry a usable, correctly located raster."""
    from pathview import highlight_nodes, pathview

    res = pathview("04110", gene_data=gene_data, species="hsa",
                   kegg_dir=workdir, out_dir=workdir, render_mode=mode,
                   output_format="svg" if mode == "svg" else "png",
                   out_suffix=mode, quiet=True)
    assert res.frame is not None, f"{mode} carried no raster"

    out = res + highlight_nodes(["1017"], color="red", width=3)
    assert not np.array_equal(out.frame.array, res.frame.array)
    saved = out.save(workdir / f"hl_{mode}.png")
    assert saved.exists() and saved.stat().st_size > 5_000
