"""CLI tests, including the compound-only invocation that crashed in 2.x."""
from __future__ import annotations

import pytest


@pytest.fixture
def gene_csv(tmp_path, gene_data):
    """The full demo set: a head() slice contains no TCA-cycle genes."""
    p = tmp_path / "genes.csv"
    gene_data.write_csv(p)
    return p


@pytest.fixture
def cpd_csv(tmp_path, cpd_data):
    p = tmp_path / "cpds.csv"
    cpd_data.write_csv(p)
    return p


def _run(argv) -> int:
    from pathview.cli import main
    return main([str(a) for a in argv])


class TestRenderCommand:
    def test_gene_only(self, workdir, gene_csv, capsys):
        rc = _run(["render", "00020", "--species", "hsa",
                   "--gene-data", gene_csv, "--kegg-dir", workdir,
                   "--out-dir", workdir, "--render-mode", "vector",
                   "--offline", "--quiet"])
        assert rc == 0
        assert (workdir / "hsa00020.pathview.png").exists()

    def test_compound_only_does_not_crash(self, workdir, cpd_csv):
        """2.x: AttributeError on gene_data.cast(...) with --cpd-data alone."""
        rc = _run(["render", "00020", "--species", "hsa",
                   "--cpd-data", cpd_csv, "--kegg-dir", workdir,
                   "--out-dir", workdir, "--render-mode", "vector",
                   "--offline", "--quiet"])
        assert rc == 0
        assert (workdir / "hsa00020.pathview.png").exists()

    def test_both_data_types(self, workdir, gene_csv, cpd_csv):
        rc = _run(["render", "00020", "--species", "human",
                   "--gene-data", gene_csv, "--cpd-data", cpd_csv,
                   "--kegg-dir", workdir, "--out-dir", workdir,
                   "--render-mode", "svg", "--output-format", "svg",
                   "--offline", "--quiet"])
        assert rc == 0
        svg = (workdir / "hsa00020.pathview.svg").read_text()
        assert "RNA-seq log2FC" in svg and "Metabolite log2FC" in svg

    def test_no_data_and_no_map_null_is_an_error(self, workdir):
        with pytest.raises(SystemExit):
            _run(["render", "00020", "--kegg-dir", workdir,
                  "--out-dir", workdir, "--offline"])

    def test_map_null_renders_without_data(self, workdir):
        rc = _run(["render", "00020", "--map-null", "--kegg-dir", workdir,
                   "--out-dir", workdir, "--render-mode", "vector",
                   "--offline", "--quiet"])
        assert rc == 0

    def test_multiple_pathways(self, workdir, gene_csv):
        rc = _run(["render", "00020", "04110", "--gene-data", gene_csv,
                   "--kegg-dir", workdir, "--out-dir", workdir,
                   "--render-mode", "vector", "--offline", "--quiet"])
        assert rc == 0
        assert (workdir / "hsa00020.pathview.png").exists()
        assert (workdir / "hsa04110.pathview.png").exists()

    def test_split_limits(self, workdir, gene_csv, cpd_csv):
        rc = _run(["render", "00020", "--gene-data", gene_csv,
                   "--cpd-data", cpd_csv, "--limit", "gene=2,cpd=1",
                   "--kegg-dir", workdir, "--out-dir", workdir,
                   "--render-mode", "vector", "--offline", "--quiet"])
        assert rc == 0

    def test_missing_file_is_a_clean_error(self, workdir):
        with pytest.raises(SystemExit) as exc:
            _run(["render", "00020", "--gene-data", workdir / "nope.csv",
                  "--kegg-dir", workdir, "--out-dir", workdir, "--offline"])
        assert "not found" in str(exc.value)

    def test_single_column_file_is_rejected(self, tmp_path, workdir):
        bad = tmp_path / "bad.csv"
        bad.write_text("entrez\n1017\n1019\n")
        with pytest.raises(SystemExit) as exc:
            _run(["render", "00020", "--gene-data", bad,
                  "--kegg-dir", workdir, "--out-dir", workdir, "--offline"])
        assert "value column" in str(exc.value)

    def test_non_numeric_columns_are_dropped(self, tmp_path, workdir):
        f = tmp_path / "mixed.csv"
        f.write_text("entrez,note,logfc\n1431,ok,1.5\n3417,ok,-1.2\n")
        rc = _run(["render", "00020", "--gene-data", f, "--kegg-dir", workdir,
                   "--out-dir", workdir, "--render-mode", "vector",
                   "--offline", "--quiet"])
        assert rc == 0


class TestOtherCommands:
    def test_species(self, capsys):
        assert _run(["species", "human"]) == 0
        assert "hsa" in capsys.readouterr().out

    def test_species_unknown_returns_nonzero(self, capsys):
        assert _run(["species", "Notareal organism"]) == 1

    def test_search(self, capsys):
        assert _run(["search", "coli", "--limit", "3"]) == 0
        assert "eco" in capsys.readouterr().out

    def test_parity(self, capsys):
        assert _run(["parity"]) == 0
        out = capsys.readouterr().out
        assert "pathview" in out and "features tracked" in out

    def test_parity_markdown(self, capsys):
        assert _run(["parity", "--markdown"]) == 0
        assert capsys.readouterr().out.startswith("| Category")

    def test_info(self, capsys):
        assert _run(["info"]) == 0
        out = capsys.readouterr().out
        assert "organisms bundled" in out and "KEGG" in out

    def test_legend(self, tmp_path, capsys):
        out = tmp_path / "leg.png"
        assert _run(["legend", "--out", out]) == 0
        assert out.exists()

    def test_version(self, capsys):
        with pytest.raises(SystemExit) as exc:
            _run(["--version"])
        assert exc.value.code == 0


class TestCliRobustness:
    def test_piping_output_does_not_traceback(self, tmp_path):
        """`pathview-plus info | head` closes the pipe early; that is normal."""
        import subprocess
        import sys

        proc = subprocess.run(
            f"{sys.executable} -m pathview.cli info | head -2",
            shell=True, capture_output=True, text=True,
            cwd=str(tmp_path.parent), env={"PATH": "/usr/bin:/bin",
                                           "PYTHONPATH": _pkg_root()},
        )
        assert "BrokenPipeError" not in proc.stderr
        assert "Traceback" not in proc.stderr


def _pkg_root() -> str:
    from pathlib import Path

    import pathview
    return str(Path(pathview.__file__).resolve().parents[1])


class TestSbgnCommand:
    def test_lists_the_collection(self, capsys):
        assert _run(["sbgn", "--list"]) == 0
        out = capsys.readouterr().out
        assert "pre-generated SBGN pathways" in out
        assert "PANTHER" in out and "MetaCyc" in out

    def test_filters_by_source(self, capsys):
        assert _run(["sbgn", "--list", "--source", "panther", "-n", "5"]) == 0
        lines = [x for x in capsys.readouterr().out.splitlines() if x.strip()]
        assert lines and all(x.startswith("panther") for x in lines)

    def test_search(self, capsys):
        assert _run(["sbgn", "--list", "--source", "metacyc",
                     "--query", "glycoly", "-n", "5"]) == 0
        assert "GLYCOLY" in capsys.readouterr().out.upper()

    def test_renders_a_local_file(self, tmp_path, fixtures_dir):
        rc = _run(["sbgn", str(fixtures_dir / "P00001.new.layout.sbgn"),
                   "--out-dir", tmp_path, "--offline", "--quiet"])
        assert rc == 0
        assert list(tmp_path.glob("*.png"))

    def test_info_reports_the_collection(self, capsys):
        assert _run(["info"]) == 0
        out = capsys.readouterr().out
        assert "SBGN collection" in out
        assert "manual only" not in out, "every source should now be available"


class TestExpansionFlags:
    def test_split_group_and_expand_node(self, workdir, gene_csv):
        rc = _run(["render", "04110", "--gene-data", gene_csv,
                   "--split-group", "--expand-node",
                   "--kegg-dir", workdir, "--out-dir", workdir,
                   "--render-mode", "vector", "--offline", "--quiet"])
        assert rc == 0
