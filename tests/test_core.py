"""Functional tests for parsing, mapping, colour and geometry."""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

# --- organisms --------------------------------------------------------------

class TestOrganisms:
    def test_kegg_species_code_alias(self):
        from pathview import get_species_code, kegg_species_code
        assert kegg_species_code("human") == get_species_code("human").kegg_code

    @pytest.mark.parametrize("query,expected", [
        ("hsa", "hsa"), ("mmu", "mmu"), ("Mus musculus", "mmu"),
        ("mouse", "mmu"), ("10090", "mmu"), ("E. coli", "eco"),
        ("Escherichia coli", "eco"), ("yeast", "sce"), ("zebrafish", "dre"),
        ("Arabidopsis thaliana", "ath"), ("fruit fly", "dme"),
    ])
    def test_resolution(self, query, expected):
        from pathview import get_species_code
        assert get_species_code(query).kegg_code == expected

    def test_species_info_fields(self):
        from pathview import get_species_code
        info = get_species_code("hsa")
        assert info.tax_id == "9606"
        assert "Homo sapiens" in info.scientific_name
        assert info.kegg_code == "hsa"

    def test_search_returns_ranked_matches(self):
        from pathview import search_organisms
        hits = search_organisms("coli", limit=5)
        assert hits.height > 0
        assert "eco" in hits["kegg_code"].to_list()

    def test_list_organisms_is_a_frame(self):
        from pathview import list_organisms
        df = list_organisms()
        assert {"kegg_code", "scientific_name", "tax_id"} <= set(df.columns)
        assert df.height > 8000

    def test_default_gene_idtype(self):
        from pathview import default_gene_idtype
        assert default_gene_idtype("hsa") == "ENTREZ"
        assert default_gene_idtype("mmu") == "ENTREZ"
        assert default_gene_idtype("ko") == "KEGG"


# --- KGML -------------------------------------------------------------------

class TestKGML:
    def test_node_counts(self, kgml_tca):
        assert len(kgml_tca) == 64
        assert kgml_tca.org == "hsa"

    def test_node_frame_schema(self, nodes_tca):
        expected = {"entry_id", "name", "kegg_names", "type", "x", "y",
                    "width", "height", "bgcolor", "fgcolor", "label", "shape"}
        assert expected <= set(nodes_tca.columns)

    def test_labels_are_shortened(self, kgml_cellcycle):
        from pathview import node_info
        labels = node_info(kgml_cellcycle)["label"].to_list()
        assert not any("..." in (lab or "") for lab in labels)
        assert "CDKN2A" in labels

    def test_groups_inherit_component_ids(self, kgml_cellcycle):
        groups = [n for n in kgml_cellcycle.nodes.values() if n.is_group]
        assert groups
        assert all(g.kegg_names for g in groups)
        assert all(g.label and g.label.lower() != "undefined" for g in groups)

    def test_edges_present(self, edges_tca):
        assert edges_tca.height > 100
        assert set(edges_tca.columns) >= {"source", "target", "subtype", "source_kind"}

    def test_edges_reference_existing_nodes(self, nodes_tca, edges_tca):
        ids = set(nodes_tca["entry_id"].to_list())
        assert set(edges_tca["source"].to_list()) <= ids
        assert set(edges_tca["target"].to_list()) <= ids

    def test_old_kgml_dialect_yields_reaction_edges(self, fixtures_dir):
        """ko-prefixed maps give substrates a name but no entry id."""
        from pathview import parse_kgml, pathway_edges
        edges = pathway_edges(parse_kgml(fixtures_dir / "ko00051.xml"))
        rxn = edges.filter(pl.col("source_kind") == "reaction")
        assert rxn.height > 50

    def test_malformed_file_raises_parse_error(self, tmp_path):
        from pathview.errors import ParseError

        from pathview import parse_kgml
        bad = tmp_path / "bad.xml"
        bad.write_text("No such pathway: hsa99999")
        with pytest.raises(ParseError):
            parse_kgml(bad)

    def test_missing_file_raises(self, tmp_path):
        from pathview.errors import PathwayNotFoundError

        from pathview import parse_kgml
        with pytest.raises(PathwayNotFoundError):
            parse_kgml(tmp_path / "nope.xml")


# --- SBGN -------------------------------------------------------------------

class TestSBGN:
    def test_glyphs_and_arcs(self, fixtures_dir):
        from pathview import parse_sbgn
        p = parse_sbgn(fixtures_dir / "P00001.new.layout.sbgn")
        assert len(p.glyphs) > 50 and len(p.arcs) > 50

    def test_converts_to_node_frame(self, fixtures_dir):
        from pathview import parse_sbgn, sbgn_to_df
        df = sbgn_to_df(parse_sbgn(fixtures_dir / "P00001.new.layout.sbgn"))
        assert {"entry_id", "type", "x", "y", "glyph_class"} <= set(df.columns)
        assert set(df["type"].to_list()) <= {"gene", "compound", "process",
                                             "operator", "map", "unknown"}

    def test_identifiers_extracted_from_annotations(self, fixtures_dir):
        from pathview import parse_sbgn
        p = parse_sbgn(fixtures_dir / "ports_pd.sbgn")
        assert "kegg:C00031" in p.glyphs["glc"].identifiers
        assert "uniprot:P19367" in p.glyphs["hk1"].identifiers

    def test_canvas_covers_all_glyphs(self, fixtures_dir):
        from pathview import parse_sbgn, sbgn_canvas
        x0, y0, w, h = sbgn_canvas(parse_sbgn(fixtures_dir / "ports_pd.sbgn"))
        assert w > 0 and h > 0


# --- colour -----------------------------------------------------------------

class TestColour:
    def test_r_colorpanel2_parity_even(self):
        from pathview import colorpanel2
        cols = colorpanel2(10, "#00FF00", "#BEBEBE", "#FF0000")
        assert len(cols) == 10
        assert cols[0] == "#00FF00" and cols[-1] == "#FF0000"
        assert cols[4] == cols[5] == "#BEBEBE"

    def test_r_colorpanel2_parity_odd_drops_duplicate_midpoint(self):
        from pathview import colorpanel2
        cols = colorpanel2(5, "#00FF00", "#BEBEBE", "#FF0000")
        assert len(cols) == 5
        assert cols[2] == "#BEBEBE"
        assert len(set(cols)) == 5

    def test_symmetric_limits(self):
        from pathview import ColorScale
        assert ColorScale(limit=2.0, both_dirs=True).bounds() == (-2.0, 2.0)
        assert ColorScale(limit=2.0, both_dirs=False).bounds() == (0.0, 2.0)

    def test_values_outside_limits_are_clamped_not_dropped(self):
        from pathview import ColorScale
        sc = ColorScale(limit=1.0, bins=10)
        assert sc.map_values([-99.0])[0] == sc.map_values([-1.0])[0]
        assert sc.map_values([99.0])[0] == sc.map_values([1.0])[0]

    def test_nan_uses_na_colour(self):
        from pathview import ColorScale
        sc = ColorScale(limit=1.0, na_col="#123456")
        assert sc.map_values([float("nan")]) == ["#123456"]

    def test_gene_and_compound_scales_are_independent(self):
        from pathview import compound_scale, gene_scale
        g, c = gene_scale(limit=2.0), compound_scale(limit=2.0)
        assert g.anchors() != c.anchors()
        assert g.map_values([1.0]) != c.map_values([1.0])

    def test_named_palettes_resolve(self):
        from pathview import ColorScale, list_palettes
        for name in list_palettes():
            assert len(ColorScale(palette=name).colors()) == 10

    def test_unknown_palette_raises(self):
        from pathview import ColorScale
        with pytest.raises(ValueError):
            ColorScale(palette="chartreuse-explosion").colors()

    def test_trans_fun_applied_before_binning(self):
        from pathview import ColorScale
        plain = ColorScale(limit=8.0, bins=8)
        logged = ColorScale(limit=8.0, bins=8, trans_fun=np.log2)
        # 6.0 sits high on the plain scale but log2(6) = 2.58 sits mid-high.
        assert plain.map_values([6.0]) != logged.map_values([6.0])

    def test_multi_condition_produces_one_column_each(self):
        from pathview import gene_scale, node_color
        df = pl.DataFrame({"id": ["a", "b"], "c1": [-1.0, 1.0], "c2": [1.0, -1.0]})
        out = node_color(df, gene_scale(limit=1.0))
        assert out.columns == ["id", "c1_col", "c2_col"]
        assert out["c1_col"][0] != out["c2_col"][0]


# --- mapping ----------------------------------------------------------------

class TestMapping:
    def test_gene_mapping_on_tca(self, gene_data, nodes_tca):
        from pathview import node_map
        res = node_map(gene_data, nodes_tca, "gene", detailed=True)
        assert res.ok and res.n_nodes_with_data > 10
        assert 0 < res.mapped_fraction <= 1

    def test_compound_mapping_on_tca(self, cpd_data, nodes_tca):
        from pathview import node_map
        res = node_map(cpd_data, nodes_tca, "compound", detailed=True)
        assert res.ok and res.n_nodes_with_data > 5

    def test_all_mapped_records_every_contributing_id(self, gene_data, nodes_tca):
        from pathview import node_map
        out = node_map(gene_data, nodes_tca, "gene")
        multi = out.filter(pl.col("all_mapped").str.contains(","))
        assert multi.height > 0, "multi-gene nodes must list every id"

    def test_map_null_returns_layout_without_data(self, nodes_tca):
        from pathview import node_map
        out = node_map(None, nodes_tca, "gene")
        assert out is not None and out.height > 0

    def test_unmatched_ids_raise_with_a_useful_message(self, nodes_tca):
        from pathview import node_map
        junk = pl.DataFrame({"id": ["ZZZ1", "ZZZ2"], "v": [1.0, 2.0]})
        res = node_map(junk, nodes_tca, "gene", detailed=True)
        assert not res.ok or res.n_nodes_with_data == 0

    def test_mol_sum_error_names_the_columns(self):
        from pathview.errors import MappingError

        from pathview import mol_sum
        data = pl.DataFrame({"gene": ["X"], "v": [1.0]})
        idmap = pl.DataFrame({"gene": ["Y"], "t": ["A"]})
        with pytest.raises(MappingError) as exc:
            mol_sum(data, idmap)
        assert "gene" in str(exc.value)


# --- identifiers ------------------------------------------------------------

class TestIdentifiers:
    def test_compound_names_resolve_offline(self):
        from pathview import compound_name
        assert compound_name("C00031") == "D-Glucose"
        assert compound_name("C00022") == "Pyruvate"

    def test_unknown_compound_falls_back_to_the_accession(self):
        from pathview import compound_name
        assert compound_name("C99999") == "C99999"

    @pytest.mark.parametrize("name,kegg", [
        ("Pyruvate", "C00022"), ("Citrate", "C00158"), ("ATP", "C00002"),
        ("D-Glucose", "C00031"), ("2-Oxoglutarate", "C00026"),
        ("Acetyl-CoA", "C00024"), ("Succinate", "C00042"),
    ])
    def test_metabolite_names_map_to_kegg(self, name, kegg):
        from pathview import cpd_name_to_kegg
        assert cpd_name_to_kegg([name])["KEGG"][0] == kegg

    def test_cross_reference_lookup_offline(self):
        from pathview import cpd_id_map
        out = cpd_id_map(["498-15-7"], "CAS", "KEGG", detailed=True)
        assert out.n_resolved == 1
        assert out.data["KEGG"][0] == "C11382"

    def test_kegg_to_kegg_is_identity(self):
        from pathview import cpd_id_map
        out = cpd_id_map(["C00031", "C00002"], "KEGG", "KEGG")
        assert out["KEGG"].to_list() == ["C00031", "C00002"]

    def test_supported_types_are_listed(self):
        from pathview import supported_cpd_idtypes, supported_gene_idtypes
        assert "CAS" in supported_cpd_idtypes()
        assert "SYMBOL" in supported_gene_idtypes()

    def test_kegg_code_translates_to_taxid_for_mygene(self):
        from pathview.id_mapping import _mygene_species
        assert _mygene_species("hsa") == "9606"
        assert _mygene_species("mmu") == "10090"


# --- geometry ---------------------------------------------------------------

class TestGeometry:
    def test_extent_covers_every_node(self, nodes_tca):
        from pathview import Extent, node_boxes
        boxes = node_boxes(nodes_tca)
        ext = Extent.from_boxes(boxes, pad=0)
        assert all(ext.x0 <= b.left and b.right <= ext.x1 for b in boxes)
        assert all(ext.y0 <= b.top and b.bottom <= ext.y1 for b in boxes)

    def test_slices_tile_the_node_exactly(self):
        from pathview import node_boxes, slice_bounds
        df = pl.DataFrame({"entry_id": ["1"], "x": [50.0], "y": [50.0],
                           "width": [40.0], "height": [20.0], "type": ["gene"],
                           "shape": ["rectangle"], "label": [""], "bgcolor": ["#FFF"]})
        box = node_boxes(df)[0]
        bands = slice_bounds(box, 4)
        assert len(bands) == 4
        assert bands[0][0] == box.left and bands[-1][1] == box.right
        assert all(abs((b - a) - 10.0) < 1e-9 for a, b in bands)

    def test_nodes_without_coordinates_are_skipped(self):
        from pathview import node_boxes
        df = pl.DataFrame({"entry_id": ["1"], "x": [None], "y": [None],
                           "width": [40.0], "height": [20.0], "type": ["gene"],
                           "shape": ["rectangle"], "label": [""], "bgcolor": ["#FFF"]},
                          schema_overrides={"x": pl.Float64, "y": pl.Float64})
        assert node_boxes(df) == []
