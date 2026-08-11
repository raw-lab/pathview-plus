"""
The pre-generated SBGN collection, its crosswalks, and the database
downloaders that were stubs in 2.x.

Index and crosswalk tests run offline (both ship in the wheel).  Tests that
need the network are marked and skipped when it is unavailable, so the suite
still passes behind a firewall.
"""
from __future__ import annotations

import pytest


def _online() -> bool:
    from pathview.cache import is_offline
    return not is_offline()


network = pytest.mark.skipif(True, reason="set PATHVIEW_TEST_NETWORK=1 to enable")

import os  # noqa: E402

if os.environ.get("PATHVIEW_TEST_NETWORK") == "1":
    network = pytest.mark.network


class TestCollectionIndex:
    def test_index_ships_in_the_wheel(self):
        from pathview import sbgn_index
        idx = sbgn_index()
        assert idx.height > 5000
        assert {"pathway_id", "source", "subdir", "filename"} <= set(idx.columns)

    def test_all_five_sources_present(self):
        from pathview import sbgn_collection_info
        info = sbgn_collection_info()
        assert set(info["by_source"]) == {"reactome", "smpdb", "panther",
                                          "metacyc", "metacrop"}
        assert info["by_source"]["panther"] > 100
        assert info["by_source"]["smpdb"] > 500
        assert info["by_source"]["metacyc"] > 2000
        assert info["total"] == sum(info["by_source"].values())

    @pytest.mark.parametrize("source,example", [
        ("panther", "P00001"), ("smpdb", "SMP00001"),
        ("metacyc", "GLYCOLYSIS"), ("reactome", "R-HSA-109688"),
    ])
    def test_known_pathways_are_indexed(self, source, example):
        from pathview import find_sbgn_pathway
        entry = find_sbgn_pathway(example)
        assert entry["source"] == source

    def test_browsing_is_filtered_and_offline(self):
        from pathview import list_sbgn_pathways
        panther = list_sbgn_pathways("panther")
        assert panther.height > 100
        assert set(panther["source"].to_list()) == {"panther"}

    def test_search_within_a_source(self):
        from pathview import list_sbgn_pathways
        hits = list_sbgn_pathways("metacyc", query="glycoly")
        assert hits.height >= 1
        assert all("GLYCOLY" in p.upper() for p in hits["pathway_id"].to_list())

    def test_unknown_source_is_rejected(self):
        from pathview import list_sbgn_pathways
        with pytest.raises(ValueError, match="Unknown source"):
            list_sbgn_pathways("wikipathways")

    def test_missing_pathway_suggests_alternatives(self):
        from pathview.errors import PathwayNotFoundError

        from pathview import find_sbgn_pathway
        with pytest.raises(PathwayNotFoundError) as exc:
            find_sbgn_pathway("P0000X")
        assert "P00001" in str(exc.value) or "collection" in str(exc.value)

    def test_url_is_well_formed(self):
        from pathview import sbgn_url
        url = sbgn_url("P00001")
        assert url.startswith("https://")
        assert url.endswith(".sbgn")


class TestCrosswalks:
    def test_crosswalk_ships_in_the_wheel(self):
        from pathview import sbgn_xref
        assert sbgn_xref().height > 500_000

    def test_routes_are_declared(self):
        from pathview import crosswalk_routes
        routes = crosswalk_routes()
        assert routes.height >= 5
        assert routes["pairs"].sum() > 500_000

    @pytest.mark.parametrize("src,dst", [
        ("SYMBOL", "pathwayCommons"), ("ENTREZ", "pathwayCommons"),
        ("KEGG", "chebi"), ("NAME", "kegg"), ("KO", "pathwayCommons"),
        ("ENTREZ", "SYMBOL"),
    ])
    def test_route_exists(self, src, dst):
        from pathview import id_route
        assert id_route(src, dst) is not None

    def test_multi_hop_routing(self):
        """Entrez reaches gene symbol only through KO and Pathway Commons."""
        from pathview import id_route
        route = id_route("ENTREZ", "SYMBOL")
        assert len(route) > 2
        assert route[0] == "entrez" and route[-1] == "symbol"

    def test_known_conversions(self):
        from pathview import map_ids_to_sbgn
        out = map_ids_to_sbgn(["CDK2"], "SYMBOL", "pathwayCommons")
        assert out["PATHWAYCOMMONS"][0] is not None
        assert map_ids_to_sbgn(["C00022"], "KEGG", "chebi")["CHEBI"][0] == "15361"

    def test_entrez_to_symbol_roundtrip(self):
        from pathview import map_ids_to_sbgn
        out = map_ids_to_sbgn(["1017", "7157"], "ENTREZ", "SYMBOL")
        assert out["SYMBOL"].to_list() == ["CDK2", "TP53"]

    def test_identity_returns_two_columns(self):
        """A same-type conversion must not collapse to a single column."""
        from pathview import map_ids_to_sbgn
        out = map_ids_to_sbgn(["C00022"], "KEGG", "KEGG")
        assert out.width == 2
        assert out[out.columns[1]][0] == "C00022"

    def test_impossible_route_names_the_options(self):
        from pathview import map_ids_to_sbgn
        with pytest.raises(ValueError, match="Available types"):
            map_ids_to_sbgn(["x"], "PDB", "chebi")

    def test_detailed_reports_coverage(self):
        from pathview import map_ids_to_sbgn
        res = map_ids_to_sbgn(["CDK2", "TP53", "NOTAGENE"], "SYMBOL",
                              "pathwayCommons", detailed=True)
        assert res.n_input == 3
        assert 0 < res.n_resolved < 3
        assert "crosswalk" in res.source


class TestDatabaseCapabilities:
    def test_every_database_is_now_available(self):
        """2.x listed PANTHER, MetaCyc and SMPDB as unsupported stubs."""
        from pathview import DATABASE_INFO, available_databases
        assert set(available_databases()) == set(DATABASE_INFO)
        for key, info in DATABASE_INFO.items():
            assert info["available"] is True, key
            assert callable(info["downloader"]), key
            assert info["note"] and info["source"], key

    @pytest.mark.parametrize("pid,expected", [
        ("hsa04110", "kegg"), ("04110", "kegg"),
        ("R-HSA-109688", "reactome"), ("P00001", "panther"),
        ("SMP00001", "smpdb"), ("GLYCOLYSIS", "metacyc"),
        ("1CMET2-PWY", "metacyc"), ("Alanine Degradation", "metacrop"),
    ])
    def test_dispatch_by_id_shape(self, pid, expected):
        from pathview import detect_database
        assert detect_database(pid) == expected

    def test_unrecognised_id_is_reported(self):
        from pathview import detect_database
        assert detect_database("!!!nonsense!!!") is None

    def test_downloader_rejects_mismatched_source(self):
        from pathview.errors import PathwayNotFoundError

        from pathview import download_panther
        with pytest.raises(PathwayNotFoundError, match="not panther"):
            download_panther("SMP00001", ".")

    def test_offline_download_explains_itself(self, tmp_path):
        from pathview.errors import NetworkError

        from pathview import download_sbgn
        with pytest.raises(NetworkError, match="[Oo]ffline"):
            download_sbgn("P00001", output_dir=tmp_path)


@network
class TestLiveCollectionDownload:
    @pytest.mark.parametrize("pid,source", [
        ("P00001", "panther"), ("SMP00001", "smpdb"), ("GLYCOLYSIS", "metacyc"),
    ])
    def test_download_and_parse(self, tmp_path, pid, source):
        from pathview.databases import DATABASE_INFO

        from pathview import arc_resolution_report, parse_sbgn
        path = DATABASE_INFO[source]["downloader"](pid, tmp_path)
        assert path.exists() and path.stat().st_size > 1000
        pw = parse_sbgn(path)
        assert len(pw.glyphs) > 10
        assert arc_resolution_report(pw)["resolution_rate"] > 0.9


class TestLocalSbgnFilesParseIdentically:
    def test_local_file_needs_no_download(self, fixtures_dir, tmp_path):
        """
        A file exported by hand from any SBGN source must work exactly like
        one from the collection — there is no second code path.
        """
        from pathview import parse_sbgn, sbgn_edges, sbgn_to_df

        src = fixtures_dir / "P00001.new.layout.sbgn"
        copy = tmp_path / "hand_exported.sbgn"
        copy.write_text(src.read_text())

        a, b = parse_sbgn(src), parse_sbgn(copy)
        assert len(a.glyphs) == len(b.glyphs)
        assert sbgn_to_df(a).height == sbgn_to_df(b).height
        assert sbgn_edges(a).height == sbgn_edges(b).height

    def test_sbgnview_accepts_a_path(self, fixtures_dir, tmp_path):
        from pathview import sbgnview
        res = sbgnview(fixtures_dir / "P00001.new.layout.sbgn",
                       out_dir=tmp_path, quiet=True)
        assert res.output_path.exists()
        assert res.diagnostics["glyphs"] > 50


class TestRDataReader:
    """The reader that made the bundled SBGN data possible."""

    def test_reads_a_named_vector_table(self):
        import pathlib

        from pathview import read_rdata
        src = pathlib.Path("/tmp/svd/SBGNview.data-master/data/chebi_kegg.RData")
        if not src.exists():
            pytest.skip("source RData not present in this environment")
        obj = read_rdata(src)
        table = list(obj.values())[0]
        assert set(table) == {"chebi", "kegg"}
        assert len(table["chebi"]) == len(table["kegg"]) > 10_000

    def test_unreadable_file_raises_rather_than_guessing(self, tmp_path):
        """
        A parser that guesses at an unknown type misaligns the byte stream and
        corrupts everything after it, so an unsupported type must raise.
        """
        from pathview import read_rdata
        bad = tmp_path / "not.rda"
        bad.write_bytes(b"this is not an RData file at all")
        with pytest.raises(ValueError, match="RDX2/RDX3"):
            read_rdata(bad)

    def test_object_listing(self):
        import pathlib

        from pathview import rdata_objects
        src = pathlib.Path("/tmp/svd/SBGNview.data-master/data/chebi_kegg.RData")
        if not src.exists():
            pytest.skip("source RData not present in this environment")
        assert "chebi_kegg" in rdata_objects(src)
