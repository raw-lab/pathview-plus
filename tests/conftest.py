"""Shared fixtures.  Every test runs offline: no test may touch the network."""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session", autouse=True)
def _offline():
    """Force offline mode for the whole session."""
    from pathview.cache import set_offline
    set_offline(True)
    yield
    set_offline(False)


@pytest.fixture(scope="session")
def fixtures_dir() -> Path:
    return FIXTURES


@pytest.fixture(scope="session")
def kgml_tca():
    from pathview import parse_kgml
    return parse_kgml(FIXTURES / "hsa00020.xml")


@pytest.fixture(scope="session")
def kgml_cellcycle():
    from pathview import parse_kgml
    return parse_kgml(FIXTURES / "hsa04110.xml")


@pytest.fixture(scope="session")
def nodes_tca(kgml_tca):
    from pathview import node_info
    return node_info(kgml_tca)


@pytest.fixture(scope="session")
def edges_tca(kgml_tca):
    from pathview import pathway_edges
    return pathway_edges(kgml_tca)


@pytest.fixture(scope="session")
def gene_data():
    from pathview import demo_gene_data
    return demo_gene_data(2)


@pytest.fixture(scope="session")
def cpd_data(nodes_tca):
    from pathview import demo_cpd_data
    ids = [c for row in nodes_tca.filter(nodes_tca["type"] == "compound")["kegg_names"].to_list()
           for c in row]
    return demo_cpd_data(ids, n_mol=18)


@pytest.fixture
def workdir(tmp_path) -> Path:
    """A scratch directory pre-loaded with the KGML/PNG fixtures."""
    for name in ("hsa00020.xml", "hsa04110.xml", "hsa04110.png"):
        src = FIXTURES / name
        if src.exists():
            shutil.copy(src, tmp_path / name)
    return tmp_path


@pytest.fixture(scope="session")
def nodes_cellcycle(kgml_cellcycle):
    from pathview import node_info
    return node_info(kgml_cellcycle)
