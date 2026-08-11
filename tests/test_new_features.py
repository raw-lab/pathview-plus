"""Batch rendering, compartment shading, group splitting and node expansion."""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")


class TestMultiPathwayBatch:
    def test_sequence_renders_every_pathway(self, workdir, gene_data):
        from pathview import PathwayResultSet, pathview
        rs = pathview(["00020", "04110"], gene_data=gene_data, species="hsa",
                      kegg_dir=workdir, out_dir=workdir,
                      render_mode="vector", quiet=True)
        assert isinstance(rs, PathwayResultSet)
        assert len(rs) == 2
        assert all(r.output_path.exists() for r in rs)

    def test_single_id_still_returns_a_single_result(self, workdir, gene_data):
        from pathview import PathwayResult, pathview
        res = pathview("00020", gene_data=gene_data, species="hsa",
                       kegg_dir=workdir, out_dir=workdir,
                       render_mode="vector", quiet=True)
        assert isinstance(res, PathwayResult)

    def test_one_failure_does_not_lose_the_rest(self, workdir, gene_data):
        from pathview import pathview
        rs = pathview(["00020", "99999", "04110"], gene_data=gene_data,
                      species="hsa", kegg_dir=workdir, out_dir=workdir,
                      render_mode="vector", quiet=True)
        assert len(rs) == 2
        assert "99999" in rs.failures
        assert bool(rs) is True

    def test_continue_on_error_false_raises(self, workdir, gene_data):
        from pathview.errors import PathwayNotFoundError

        from pathview import pathview
        with pytest.raises(PathwayNotFoundError):
            pathview(["00020", "99999"], gene_data=gene_data, species="hsa",
                     kegg_dir=workdir, out_dir=workdir, render_mode="vector",
                     quiet=True, continue_on_error=False)

    def test_set_indexing_and_iteration(self, workdir, gene_data):
        from pathview import pathview
        rs = pathview(["00020", "04110"], gene_data=gene_data, species="hsa",
                      kegg_dir=workdir, out_dir=workdir,
                      render_mode="vector", quiet=True)
        assert rs["00020"].pathway_name == "hsa00020"
        assert rs[1].pathway_name == "hsa04110"
        assert "00020" in rs
        assert len(list(rs)) == 2
        assert len(rs.output_paths) == 2

    def test_modifier_broadcasts_across_the_batch(self, workdir, gene_data):
        from pathview import highlight_nodes, pathview
        rs = pathview(["00020", "04110"], gene_data=gene_data, species="hsa",
                      kegg_dir=workdir, out_dir=workdir,
                      render_mode="vector", quiet=True)
        before = {k: v.frame.array.copy() for k, v in rs.items()}
        out = rs + highlight_nodes(["1017", "1431"], color="red", width=4)
        changed = [not np.array_equal(before[k], v.frame.array)
                   for k, v in out.items()]
        assert any(changed), "highlight reached no pathway in the batch"

    def test_frame_report(self, workdir, gene_data):
        from pathview import pathview
        rs = pathview(["00020", "99999"], gene_data=gene_data, species="hsa",
                      kegg_dir=workdir, out_dir=workdir,
                      render_mode="vector", quiet=True)
        df = rs.to_frame()
        assert df.height == 2
        assert set(df["status"].to_list()) == {"ok"} | {
            s for s in df["status"].to_list() if s != "ok"}

    def test_empty_sequence_is_rejected(self, workdir, gene_data):
        from pathview import pathview
        with pytest.raises(ValueError, match="empty sequence"):
            pathview([], gene_data=gene_data, species="hsa",
                     kegg_dir=workdir, out_dir=workdir, quiet=True)


class TestGroupSplitting:
    def test_groups_are_replaced_by_components(self, nodes_cellcycle):
        from pathview import split_groups
        before = nodes_cellcycle.filter(pl.col("type") == "group").height
        assert before > 0
        out = split_groups(nodes_cellcycle)
        assert out.filter(pl.col("type") == "group").height == 0
        assert "parent_group" in out.columns
        assert out.filter(pl.col("parent_group").is_not_null()).height > 0

    def test_reports_what_it_did(self, nodes_cellcycle):
        from pathview import split_groups
        res = split_groups(nodes_cellcycle, detailed=True)
        assert res.n_groups_split > 0
        assert "complexes split" in res.summary()

    def test_is_a_no_op_without_groups(self, nodes_tca):
        from pathview import split_groups
        out = split_groups(nodes_tca)
        assert out.height == nodes_tca.height


class TestNodeExpansion:
    def test_multi_gene_nodes_become_several(self, nodes_cellcycle):
        from pathview import expand_nodes
        res = expand_nodes(nodes_cellcycle, detailed=True)
        assert res.n_after > res.n_before
        assert res.n_nodes_expanded > 0
        assert res.data.filter(pl.col("expanded")).height > 0

    def test_expanded_nodes_tile_the_original_box(self, nodes_cellcycle):
        """
        Every expanded node must cover exactly the area of the node it
        replaced — checked for all of them, not one, because group_by gives no
        ordering guarantee and a sampled check would pass or fail at random.
        """
        from pathview import expand_nodes

        out = expand_nodes(nodes_cellcycle)
        originals = {r["entry_id"]: r for r in nodes_cellcycle.iter_rows(named=True)}
        expanded = out.filter(pl.col("expanded"))
        assert expanded.height > 0

        areas: dict[str, float] = {}
        for row in expanded.iter_rows(named=True):
            areas[row["parent_entry"]] = areas.get(row["parent_entry"], 0.0) + (
                row["width"] * row["height"])

        for parent, area in areas.items():
            original = originals[parent]
            expected = (original["width"] or 0) * (original["height"] or 0)
            assert abs(area - expected) < 1e-6, (
                f"node {parent} covered {area:.2f} after expansion, "
                f"was {expected:.2f}")

    def test_each_expanded_node_has_one_identifier(self, nodes_cellcycle):
        from pathview import expand_nodes
        out = expand_nodes(nodes_cellcycle)
        expanded = out.filter(pl.col("expanded"))
        assert all(len(n) == 1 for n in expanded["kegg_names"].to_list())

    def test_huge_nodes_are_left_alone(self):
        """A 40-way subdivision would be unreadable; not expanding is better."""
        from pathview import expand_nodes
        df = pl.DataFrame({
            "entry_id": ["1"], "name": ["x"],
            "kegg_names": [[str(i) for i in range(40)]],
            "type": ["gene"], "x": [10.0], "y": [10.0], "width": [46.0],
            "height": [17.0], "bgcolor": ["#FFF"], "fgcolor": ["#000"],
            "label": ["x"], "shape": ["rectangle"], "reaction": [""],
            "component": [""], "size": [1], "link": [""],
        }, schema_overrides={"kegg_names": pl.List(pl.String)})
        assert expand_nodes(df).height == 1

    def test_pathview_wires_both_options(self, workdir, gene_data):
        from pathview import pathview
        plain = pathview("04110", gene_data=gene_data, species="hsa",
                         kegg_dir=workdir, out_dir=workdir,
                         render_mode="vector", out_suffix="plain", quiet=True)
        expanded = pathview("04110", gene_data=gene_data, species="hsa",
                            kegg_dir=workdir, out_dir=workdir,
                            render_mode="vector", out_suffix="exp",
                            split_group=True, expand_node=True, quiet=True)
        assert expanded.node_data.height > plain.node_data.height
        assert "expansion" in expanded.diagnostics

    def test_edges_survive_expansion(self, workdir, gene_data):
        """Unremapped edges silently vanish after expansion."""
        from pathview import pathview
        res = pathview("04110", gene_data=gene_data, species="hsa",
                       kegg_dir=workdir, out_dir=workdir, render_mode="vector",
                       split_group=True, expand_node=True, quiet=True)
        assert res.edge_data.height > 0
        ids = set(res.node_data["entry_id"].to_list())
        assert set(res.edge_data["source"].to_list()) <= ids
        assert set(res.edge_data["target"].to_list()) <= ids


class TestCompartments:
    def test_compartments_are_extracted(self, fixtures_dir):
        from pathview import parse_sbgn, sbgn_compartments
        comps = sbgn_compartments(parse_sbgn(fixtures_dir / "ports_pd.sbgn"))
        assert comps.height == 1
        assert comps["label"][0] == "cytosol"

    def test_ordered_largest_first(self, fixtures_dir):
        from pathview import parse_sbgn, sbgn_compartments
        comps = sbgn_compartments(parse_sbgn(fixtures_dir / "P00001.new.layout.sbgn"))
        areas = comps["area"].to_list()
        assert areas == sorted(areas, reverse=True)

    def test_shading_changes_the_render(self, fixtures_dir, tmp_path):
        from pathview import (
            keggview_vector,
            parse_sbgn,
            sbgn_compartments,
            sbgn_edges,
            sbgn_to_df,
        )
        pw = parse_sbgn(fixtures_dir / "P00001.new.layout.sbgn")
        nodes, edges = sbgn_to_df(pw), sbgn_edges(pw)
        comps = sbgn_compartments(pw)

        plain = keggview_vector(nodes, edges, pathway_name="p", out_dir=tmp_path,
                                out_suffix="plain", output_format="png",
                                plot_col_key=False)
        shaded = keggview_vector(nodes, edges, pathway_name="p", out_dir=tmp_path,
                                 out_suffix="shaded", output_format="png",
                                 plot_col_key=False, compartments=comps)
        assert plain.read_bytes() != shaded.read_bytes()

    def test_no_compartments_is_harmless(self, nodes_tca, edges_tca, tmp_path):
        from pathview import keggview_vector
        out = keggview_vector(nodes_tca, edges_tca, pathway_name="t",
                              out_dir=tmp_path, output_format="png",
                              compartments=None, plot_col_key=False)
        assert out.exists()


class TestSbgnview:
    def test_renders_a_local_file(self, fixtures_dir, tmp_path):
        from pathview import sbgnview
        res = sbgnview(fixtures_dir / "P00001.new.layout.sbgn",
                       out_dir=tmp_path, quiet=True)
        assert res.output_path.exists()
        assert res.node_data.height > 50
        assert res.diagnostics["compartments"] == 2

    def test_maps_gene_data_onto_glyphs(self, fixtures_dir, tmp_path):
        from pathview import sbgnview
        df = pl.DataFrame({"symbol": ["CDK2", "TP53", "CCND1"],
                           "log2fc": [1.5, -1.2, 0.8]})
        res = sbgnview(fixtures_dir / "P00001.new.layout.sbgn", gene_data=df,
                       gene_idtype="SYMBOL", out_dir=tmp_path, quiet=True)
        assert res.output_path.exists()

    def test_compartments_can_be_turned_off(self, fixtures_dir, tmp_path):
        from pathview import sbgnview
        a = sbgnview(fixtures_dir / "P00001.new.layout.sbgn", out_dir=tmp_path,
                     out_suffix="on", show_compartments=True, quiet=True)
        b = sbgnview(fixtures_dir / "P00001.new.layout.sbgn", out_dir=tmp_path,
                     out_suffix="off", show_compartments=False, quiet=True)
        assert a.output_path.read_bytes() != b.output_path.read_bytes()

    def test_batch(self, fixtures_dir, tmp_path):
        from pathview import PathwayResultSet, sbgnview_batch
        rs = sbgnview_batch(
            [fixtures_dir / "P00001.new.layout.sbgn", fixtures_dir / "ports_pd.sbgn"],
            out_dir=tmp_path, quiet=True)
        assert isinstance(rs, PathwayResultSet)
        assert len(rs) == 2

    def test_unknown_pathway_is_recorded_not_raised(self, tmp_path):
        from pathview import sbgnview
        rs = sbgnview(["NOT_A_PATHWAY_XYZ"], out_dir=tmp_path,
                      sbgn_dir=tmp_path, quiet=True)
        assert len(rs) == 0 and rs.failures
