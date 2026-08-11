"""Rendering tests: every mode must produce a real, non-trivial file."""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")


class TestRenderModes:
    @pytest.mark.parametrize("mode,fmt", [
        ("vector", "png"), ("vector", "pdf"), ("vector", "svg"),
        ("svg", "svg"), ("graph", "png"),
    ])
    def test_mode_writes_a_file(self, workdir, gene_data, cpd_data, mode, fmt):
        from pathview import pathview
        res = pathview("00020", gene_data=gene_data, cpd_data=cpd_data,
                       species="hsa", kegg_dir=workdir, out_dir=workdir,
                       render_mode=mode, output_format=fmt, quiet=True)
        assert res.output_path.exists()
        assert res.output_path.stat().st_size > 5_000

    def test_native_overlays_the_kegg_png(self, workdir, gene_data):
        from pathview import pathview
        res = pathview("04110", gene_data=gene_data, species="hsa",
                       kegg_dir=workdir, out_dir=workdir,
                       render_mode="native", quiet=True)
        assert res.output_path.exists()
        assert res.image_array is not None and res.image_array.ndim == 3

    def test_auto_falls_back_to_vector_without_a_background(self, workdir, gene_data):
        from pathview import pathview
        (workdir / "hsa00020.png").unlink(missing_ok=True)
        res = pathview("00020", gene_data=gene_data, species="hsa",
                       kegg_dir=workdir, out_dir=workdir,
                       render_mode="auto", quiet=True)
        assert res.output_path.exists()

    def test_unknown_mode_raises(self, workdir, gene_data):
        from pathview import pathview
        with pytest.raises(ValueError):
            pathview("00020", gene_data=gene_data, species="hsa",
                     kegg_dir=workdir, out_dir=workdir,
                     render_mode="hologram", quiet=True)


class TestColourApplication:
    def test_data_changes_the_picture(self, workdir, nodes_tca):
        """Two opposite datasets must not produce identical rasters."""
        from PIL import Image

        from pathview import pathview

        up = pl.DataFrame({"entrez": ["1431", "3417", "3418"], "v": [2.0, 2.0, 2.0]})
        down = pl.DataFrame({"entrez": ["1431", "3417", "3418"], "v": [-2.0, -2.0, -2.0]})

        a = pathview("00020", gene_data=up, species="hsa", kegg_dir=workdir,
                     out_dir=workdir, render_mode="vector", out_suffix="up",
                     quiet=True, limit=2.0)
        b = pathview("00020", gene_data=down, species="hsa", kegg_dir=workdir,
                     out_dir=workdir, render_mode="vector", out_suffix="down",
                     quiet=True, limit=2.0)

        ia = np.array(Image.open(a.output_path).convert("RGB"))
        ib = np.array(Image.open(b.output_path).convert("RGB"))
        assert not np.array_equal(ia, ib)

    def test_both_colour_keys_appear_for_dual_omics(self, workdir, gene_data, cpd_data):
        from pathview import pathview
        res = pathview("00020", gene_data=gene_data, cpd_data=cpd_data,
                       species="hsa", kegg_dir=workdir, out_dir=workdir,
                       render_mode="svg", output_format="svg", quiet=True)
        svg = res.output_path.read_text()
        assert "RNA-seq log2FC" in svg
        assert "Metabolite log2FC" in svg

    def test_gene_only_run_draws_one_key(self, workdir, gene_data):
        from pathview import pathview
        res = pathview("00020", gene_data=gene_data, species="hsa",
                       kegg_dir=workdir, out_dir=workdir, render_mode="svg",
                       output_format="svg", quiet=True)
        svg = res.output_path.read_text()
        assert "RNA-seq log2FC" in svg
        assert "Metabolite log2FC" not in svg


class TestSVGValidity:
    def test_svg_is_well_formed(self, workdir, gene_data, cpd_data):
        from xml.etree import ElementTree as ET

        from pathview import pathview
        res = pathview("00020", gene_data=gene_data, cpd_data=cpd_data,
                       species="hsa", kegg_dir=workdir, out_dir=workdir,
                       render_mode="svg", output_format="svg", quiet=True)
        ET.parse(res.output_path)             # raises if malformed

    def test_svg_element_ids_are_unique(self, workdir, gene_data):
        """2.x emitted one <marker id='marker_arrow'> per edge."""
        import re

        from pathview import pathview
        res = pathview("00020", gene_data=gene_data, species="hsa",
                       kegg_dir=workdir, out_dir=workdir, render_mode="svg",
                       output_format="svg", quiet=True)
        ids = re.findall(r'\sid="([^"]+)"', res.output_path.read_text())
        assert len(ids) == len(set(ids)), "duplicate element ids in SVG"

    def test_svg_contains_edges(self, workdir, gene_data):
        from pathview import pathview
        res = pathview("00020", gene_data=gene_data, species="hsa",
                       kegg_dir=workdir, out_dir=workdir, render_mode="svg",
                       output_format="svg", quiet=True)
        assert 'class="pv-edge"' in res.output_path.read_text()


class TestLegends:
    def test_kegg_legend_saves_without_a_display(self, tmp_path):
        """2.x ended in plt.show() and produced nothing headless."""
        from pathview import kegg_legend
        out = tmp_path / "legend.png"
        fig = kegg_legend(out_path=out)
        assert fig is not None and out.exists() and out.stat().st_size > 5_000

    def test_sbgn_legend_saves(self, tmp_path):
        from pathview import sbgn_legend
        out = tmp_path / "sbgn.png"
        sbgn_legend(out_path=out)
        assert out.exists() and out.stat().st_size > 5_000

    @pytest.mark.parametrize("kind", ["both", "edge", "node"])
    def test_legend_variants(self, tmp_path, kind):
        from pathview import kegg_legend
        out = tmp_path / f"{kind}.png"
        kegg_legend(legend_type=kind, out_path=out)
        assert out.exists()


class TestGraphView:
    @pytest.mark.parametrize("layout", ["kgml", "spring", "kamada_kawai",
                                        "circular", "shell"])
    def test_layouts(self, workdir, nodes_tca, edges_tca, layout):
        from pathview import available_layouts, keggview_graph
        if layout not in available_layouts():
            pytest.skip(f"{layout} needs an optional dependency")
        out = keggview_graph(nodes_tca, edges_tca, pathway_name="t",
                             out_dir=workdir, layout=layout,
                             output_format="png", out_suffix=layout)
        assert out.exists() and out.stat().st_size > 5_000

    def test_missing_optional_layout_dependency_is_explained(self, workdir,
                                                             nodes_tca, edges_tca):
        """A missing optional dep must name itself, not surface as ImportError."""
        from pathview.errors import RenderError

        from pathview import available_layouts, keggview_graph

        if "kamada_kawai" in available_layouts():
            pytest.skip("SciPy is installed")
        with pytest.raises(RenderError, match="SciPy"):
            keggview_graph(nodes_tca, edges_tca, pathway_name="t",
                           out_dir=workdir, layout="kamada_kawai")

    def test_unknown_layout_lists_the_valid_ones(self, workdir, nodes_tca, edges_tca):
        from pathview import keggview_graph
        with pytest.raises(ValueError, match="kgml"):
            keggview_graph(nodes_tca, edges_tca, pathway_name="t",
                           out_dir=workdir, layout="hyperbolic")

    def test_metrics(self, nodes_tca, edges_tca):
        from pathview import build_graph, pathway_metrics
        m = pathway_metrics(build_graph(nodes_tca, edges_tca))
        assert m["nodes"] > 0 and m["edges"] > 0
        assert 0 <= m["density"] <= 1
        assert len(m["hubs"]) == 5


class TestThemes:
    @pytest.mark.parametrize("theme", ["publication", "slate", "dark"])
    def test_theme_renders(self, workdir, gene_data, theme):
        from pathview import pathview
        res = pathview("00020", gene_data=gene_data, species="hsa",
                       kegg_dir=workdir, out_dir=workdir, render_mode="vector",
                       theme=theme, out_suffix=theme, quiet=True)
        assert res.output_path.exists()


class TestPainting:
    def test_paint_preserves_dark_pixels(self):
        """Gene symbols and glyph outlines must survive the overlay."""
        from pathview import gene_scale, node_color, paint_nodes

        img = np.full((60, 60, 3), 255, dtype=np.uint8)
        img[28:32, 20:40] = 0                              # simulated label text
        nodes = pl.DataFrame({"entry_id": ["1"], "x": [30.0], "y": [30.0],
                              "width": [40.0], "height": [18.0], "type": ["gene"],
                              "shape": ["rectangle"], "label": ["X"],
                              "bgcolor": ["#FFFFFF"], "v": [2.0]})
        cols = node_color(nodes.rename({"entry_id": "id"}),
                          gene_scale(limit=2.0), id_col="id", value_cols=["v"])
        out = paint_nodes(img.copy(), nodes, cols, "gene")
        assert (out[28:32, 20:40] == 0).all(), "dark label pixels were painted over"
        assert not np.array_equal(out, img), "node was not painted at all"

    def test_circle_painting_is_antialiased(self):
        from pathview import compound_scale, node_color, paint_nodes

        img = np.full((60, 60, 3), 255, dtype=np.uint8)
        nodes = pl.DataFrame({"entry_id": ["1"], "x": [30.0], "y": [30.0],
                              "width": [16.0], "height": [16.0],
                              "type": ["compound"], "shape": ["circle"],
                              "label": ["C"], "bgcolor": ["#FFFFFF"], "v": [2.0]})
        cols = node_color(nodes.rename({"entry_id": "id"}),
                          compound_scale(limit=2.0), id_col="id", value_cols=["v"])
        out = paint_nodes(img.copy(), nodes, cols, "compound")
        assert len(np.unique(out.reshape(-1, 3), axis=0)) > 3, "edges are aliased"
