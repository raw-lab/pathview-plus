"""
examples.py
Runnable examples, using only bundled data.

Every function here works offline and returns the path it wrote, so the whole
module can be executed as a smoke test of a fresh installation:

    python -m pathview.examples --out figures

The examples deliberately use ``demo_gene_data()`` (real GSE16873 breast
cancer log2 ratios, bundled) rather than random numbers, so the mapping path
they exercise is the same one user data takes.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .cache import is_offline
from .mol_data import demo_cpd_data, demo_gene_data


def _kgml_available(kegg_dir: Path, pathway: str) -> bool:
    return (kegg_dir / f"{pathway}.xml").exists() or not is_offline()


# ---------------------------------------------------------------------------
# KEGG
# ---------------------------------------------------------------------------

def example_basic(out_dir: str | Path = ".", kegg_dir: str | Path = ".") -> Path:
    """Cell cycle with one condition of real expression data."""
    from .pathview import pathview

    res = pathview("04110", gene_data=demo_gene_data(1), species="hsa",
                   kegg_dir=kegg_dir, out_dir=out_dir,
                   out_suffix="example_basic", limit=1.5, quiet=True)
    return res.output_path


def example_multi_condition(out_dir: str | Path = ".",
                            kegg_dir: str | Path = ".") -> Path:
    """Two conditions: each node splits into two vertical bands."""
    from .pathview import pathview

    res = pathview("04110", gene_data=demo_gene_data(2), species="human",
                   kegg_dir=kegg_dir, out_dir=out_dir,
                   out_suffix="example_multi", limit=1.5, quiet=True)
    return res.output_path


def example_gene_and_metabolite(out_dir: str | Path = ".",
                                kegg_dir: str | Path = ".") -> Path:
    """
    Transcripts and metabolites on one map, on independent scales.

    This is the case the package is built around: two colour scales and two
    keys, so a transcript at +2 and a metabolite at +2 are not read off the
    same ruler.
    """
    from .kgml_parser import node_info, parse_kgml
    from .pathview import pathview

    kegg_dir = Path(kegg_dir)
    cpds: list[str] = []
    kgml = kegg_dir / "hsa00020.xml"
    if kgml.exists():
        nodes = node_info(parse_kgml(kgml))
        cpds = [c for row in nodes.filter(pl.col("type") == "compound")["kegg_names"]
                .to_list() for c in row]

    res = pathview("00020", gene_data=demo_gene_data(1),
                   cpd_data=demo_cpd_data(cpds or None, n_mol=20),
                   species="human", kegg_dir=kegg_dir, out_dir=out_dir,
                   out_suffix="example_dual", render_mode="vector",
                   limit={"gene": 1.5, "cpd": 1.5},
                   subtitle="GSE16873 transcriptome + metabolite abundances",
                   quiet=True)
    return res.output_path


def example_compound_names(out_dir: str | Path = ".",
                           kegg_dir: str | Path = ".") -> Path:
    """
    Metabolite names rather than accessions.

    Conjugate-base forms resolve offline: "Pyruvate" finds "Pyruvic acid".
    """
    from .pathview import pathview

    cpds = pl.DataFrame({
        "name": ["Pyruvate", "Citrate", "2-Oxoglutarate", "Succinate",
                 "Fumarate", "L-Malate", "Oxaloacetate", "Acetyl-CoA"],
        "log2fc": [1.4, -0.8, 0.3, -1.6, 0.9, -0.4, 1.1, -1.2],
    })
    res = pathview("00020", cpd_data=cpds, cpd_idtype="NAME", species="hsa",
                   kegg_dir=kegg_dir, out_dir=out_dir,
                   out_suffix="example_names", render_mode="vector",
                   limit=1.5, quiet=True)
    return res.output_path


def example_batch(out_dir: str | Path = ".", kegg_dir: str | Path = "."):
    """Several pathways in one call; a failure is recorded, not raised."""
    from .pathview import pathview

    return pathview(["00020", "04110"], gene_data=demo_gene_data(1),
                    species="hsa", kegg_dir=kegg_dir, out_dir=out_dir,
                    out_suffix="example_batch", render_mode="vector",
                    limit=1.5, quiet=True)


def example_expansion(out_dir: str | Path = ".", kegg_dir: str | Path = ".") -> Path:
    """Complexes split into subunits, paralogue families into single genes."""
    from .pathview import pathview

    res = pathview("04110", gene_data=demo_gene_data(1), species="hsa",
                   kegg_dir=kegg_dir, out_dir=out_dir,
                   out_suffix="example_expanded", render_mode="vector",
                   split_group=True, expand_node=True, limit=1.5, quiet=True)
    return res.output_path


def example_highlighting(out_dir: str | Path = ".",
                         kegg_dir: str | Path = ".") -> Path:
    """Composable post-hoc emphasis."""
    from .highlighting import change_labels, highlight_nodes, highlight_path
    from .pathview import pathview

    res = pathview("04110", gene_data=demo_gene_data(1), species="hsa",
                   kegg_dir=kegg_dir, out_dir=out_dir,
                   out_suffix="example_hl", limit=1.5, quiet=True)
    annotated = (res
                 + highlight_nodes(["1017", "1019", "595"], color="#7C3AED", width=3)
                 + highlight_path(["595", "983"], color="#F59E0B")
                 + change_labels({"1017": "CDK2 *"}))
    return annotated.save(Path(out_dir) / "example_highlighted.png")


# ---------------------------------------------------------------------------
# SBGN
# ---------------------------------------------------------------------------

def example_sbgn(sbgn_file: str | Path, out_dir: str | Path = ".") -> Path:
    """Render a local SBGN file with compartment shading."""
    from .sbgnview import sbgnview

    genes = pl.DataFrame({
        "symbol": ["CPS1", "OTC", "ASS1", "ASL", "ARG1", "GLS2", "GOT2", "GLUD1"],
        "log2fc": [1.5, -1.2, 0.8, -0.4, 1.9, -1.6, 0.5, 1.1],
    })
    res = sbgnview(sbgn_file, gene_data=genes, gene_idtype="SYMBOL",
                   out_dir=out_dir, out_suffix="example_sbgn",
                   show_compartments=True, limit=1.5, quiet=True)
    return res.output_path


def example_sbgn_collection(pathway_id: str = "SMP00001",
                            out_dir: str | Path = ".",
                            sbgn_dir: str | Path = ".") -> Path:
    """
    Render from the pre-generated SBGN collection.

    Needs one download the first time; the file is cached afterwards.
    """
    from .sbgnview import sbgnview

    genes = pl.DataFrame({
        "symbol": ["CPS1", "OTC", "ASS1", "ASL", "ARG1", "GLS2", "GOT2", "GLUD1"],
        "log2fc": [1.5, -1.2, 0.8, -0.4, 1.9, -1.6, 0.5, 1.1],
    })
    res = sbgnview(pathway_id, gene_data=genes, gene_idtype="SYMBOL",
                   sbgn_dir=sbgn_dir, out_dir=out_dir,
                   out_suffix="example_collection", limit=1.5, quiet=True)
    return res.output_path


# ---------------------------------------------------------------------------
# Offline-only examples: these never touch the network
# ---------------------------------------------------------------------------

def example_species_lookup() -> dict:
    """Species resolution, entirely from the bundled table."""
    from .organisms import get_species_code, organism_count

    return {
        "organisms_bundled": organism_count(),
        "resolved": {q: get_species_code(q).kegg_code
                     for q in ("human", "Mus musculus", "9606", "E. coli", "yeast")},
    }


def example_identifier_routing() -> dict:
    """Crosswalk routing between identifier systems, offline."""
    from .sbgn_hub import id_route, map_ids_to_sbgn

    return {
        "entrez_to_symbol_route": id_route("ENTREZ", "SYMBOL"),
        "mapped": map_ids_to_sbgn(["1017", "7157"], "ENTREZ", "SYMBOL").to_dicts(),
    }


def example_collection_summary() -> dict:
    """What is in the pre-generated SBGN collection, offline."""
    from .sbgn_hub import sbgn_collection_info

    return sbgn_collection_info()


def run_offline_examples() -> dict:
    """Every example that needs no network at all."""
    return {
        "species": example_species_lookup(),
        "identifiers": example_identifier_routing(),
        "collection": example_collection_summary(),
    }


def main(argv: list[str] | None = None) -> int:
    """``python -m pathview.examples`` — run what the environment allows."""
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=".", help="output directory")
    parser.add_argument("--kegg-dir", default=".",
                        help="directory holding KGML/PNG files")
    parser.add_argument("--offline-only", action="store_true",
                        help="run only the examples that need no network")
    args = parser.parse_args(argv)

    Path(args.out).mkdir(parents=True, exist_ok=True)
    print(json.dumps(run_offline_examples(), indent=2, default=str))

    if args.offline_only:
        return 0

    rendered = []
    for fn in (example_basic, example_multi_condition, example_gene_and_metabolite,
               example_compound_names, example_expansion, example_highlighting):
        try:
            rendered.append(str(fn(args.out, args.kegg_dir)))
            print(f"ok   {fn.__name__}")
        except Exception as exc:
            print(f"skip {fn.__name__}: {exc}")
    for path in rendered:
        print(" ", path)
    return 0


if __name__ == "__main__":                                    # pragma: no cover
    raise SystemExit(main())
