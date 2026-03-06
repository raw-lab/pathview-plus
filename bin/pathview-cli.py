#!/usr/bin/env python3
"""
pathview.py  –  master CLI entry point
=======================================
Visualise gene / compound expression data on KEGG pathway diagrams.

Usage
-----
    python pathview.py --pathway-id 04110 --gene-data gene_expr.tsv
    python pathview.py --pathway-id 04110 --species mmu --gene-data gene_expr.tsv
    python pathview.py --pathway-id 04110 --gene-data gd.tsv --cpd-data cpd.tsv \\
                       --gene-idtype SYMBOL --cpd-idtype KEGG
    python pathview.py --legend

Module layout
-------------
    pathview/               ← importable package
      __init__.py           ← public API re-exports
      constants.py          ← shared types and literals
      utils.py              ← string helpers, numeric aggregators
      id_mapping.py         ← gene / compound ID conversion
      mol_data.py           ← mol_sum, sim_mol_data
      kegg_api.py           ← species lookup, file download
      kgml_parser.py        ← KGML XML → dataclasses + DataFrame
      color_mapping.py      ← colormaps, node_color, draw_color_key
      node_mapping.py       ← node_map
      rendering.py          ← keggview_native, keggview_graph, kegg_legend
      pathview.py           ← core orchestrator function

    pathview.py             ← this file  (CLI front-end)

Dependencies
------------
    pip install polars requests matplotlib seaborn numpy Pillow networkx
"""

from __future__ import annotations

import argparse
import sys

import polars as pl

from pathview.rendering import kegg_legend
from pathview.orchestrator import pathview
from pathview.mol_data import sim_mol_data


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pathview",
        description="Overlay gene/compound data on KEGG pathway diagrams.",
        formatter_class=argparse.RawDescriptionHelpFormatter,#argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python pathview.py --pathway-id 04110 --gene-data expr.tsv\n"
            "  python pathview.py --pathway-id hsa04110 --species hsa "
            "--gene-idtype SYMBOL --gene-data expr.tsv\n"
            "  python pathview.py --legend\n"
            "  python pathview.py --simulate --pathway-id 04110"
        ),
    )

    # ---- Pathway -----------------------------------------------------------
    p.add_argument(
        "--pathway-id",
        help="KEGG pathway number, e.g. '04110' or 'hsa04110'.",
    )

    # ---- Input data --------------------------------------------------------
    data = p.add_argument_group("Input data")
    data.add_argument(
        "--gene-data", metavar="TSV",
        help="TSV file: first column = gene IDs, remaining = expression values.",
    )
    data.add_argument(
        "--cpd-data", metavar="TSV",
        help="TSV file: first column = compound IDs, remaining = abundance.",
    )
    data.add_argument(
        "--gene-idtype", default="ENTREZ",
        metavar="TYPE",
        help="Input gene ID type: ENTREZ, SYMBOL, UNIPROT, ENSEMBL, KEGG.",
    )
    data.add_argument(
        "--cpd-idtype", default="KEGG",
        metavar="TYPE",
        help="Input compound ID type: KEGG, PUBCHEM, CHEBI.",
    )

    # ---- Species & paths ---------------------------------------------------
    run = p.add_argument_group("Species and paths")
    run.add_argument("--species",  default="hsa",       help="KEGG species code.")
    run.add_argument("--kegg-dir", default=".", metavar="DIR",
                     help="Directory for downloaded KEGG files and output images.")
    run.add_argument("--out-suffix", default="pathview",
                     help="Suffix appended to each output filename.")

    # ---- Rendering ---------------------------------------------------------
    rend = p.add_argument_group("Rendering")
    rend.add_argument(
        "--kegg-native", action=argparse.BooleanOptionalAction, default=True,
        help="Use the KEGG PNG background (native) or a NetworkX graph layout.",
    )
    rend.add_argument(
        "--output-format", default="png",
        choices=["png", "pdf", "svg"],
        help="Output format: png (pixel-based), pdf (vector graph), or svg (vector native).",
    )
    rend.add_argument(
        "--map-symbol", action=argparse.BooleanOptionalAction, default=True,
        help="Replace Entrez IDs with gene symbols in node labels.",
    )
    rend.add_argument(
        "--node-sum", default="sum",
        choices=["sum", "mean", "median", "max", "max_abs", "random"],
        help="Aggregation method for multiple probes mapping to one node.",
    )
    rend.add_argument(
        "--min-nnodes", type=int, default=3,
        help="Skip pathways with fewer than this many mappable nodes.",
    )
    rend.add_argument(
        "--no-signature", action="store_true",
        help="Suppress the 'Rendered by pathview.py' watermark.",
    )
    rend.add_argument(
        "--no-col-key", action="store_true",
        help="Suppress the colour-scale legend bar.",
    )

    # ---- Colour scale ------------------------------------------------------
    col = p.add_argument_group("Colour scale")
    col.add_argument("--limit-gene", type=float, default=1.0,
                     help="Symmetric colour-scale limit for gene data (±value).")
    col.add_argument("--limit-cpd",  type=float, default=1.0,
                     help="Symmetric colour-scale limit for compound data.")
    col.add_argument("--bins-gene",  type=int,   default=10,
                     help="Colour bins for gene data.")
    col.add_argument("--bins-cpd",   type=int,   default=10,
                     help="Colour bins for compound data.")
    col.add_argument("--low-gene",  default="green",  help="Low-end gene colour.")
    col.add_argument("--mid-gene",  default="gray",   help="Mid-point gene colour.")
    col.add_argument("--high-gene", default="red",    help="High-end gene colour.")
    col.add_argument("--low-cpd",   default="blue",   help="Low-end compound colour.")
    col.add_argument("--mid-cpd",   default="gray",   help="Mid-point compound colour.")
    col.add_argument("--high-cpd",  default="yellow", help="High-end compound colour.")

    # ---- Utilities ---------------------------------------------------------
    util = p.add_argument_group("Utilities")
    util.add_argument(
        "--legend",
        action="store_true",
        help="Display the KEGG element legend and exit.",
    )
    util.add_argument(
        "--simulate",
        action="store_true",
        help="Generate and use simulated gene data (requires --pathway-id).",
    )
    util.add_argument(
        "--n-sim", type=int, default=200,
        help="Number of molecules in simulated data (used with --simulate).",
    )

    return p


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)

    # -- Legend only ---------------------------------------------------------
    if args.legend:
        kegg_legend()
        return

    # -- Require pathway-id for everything else ------------------------------
    if not args.pathway_id:
        _build_parser().error("--pathway-id is required (unless using --legend).")

    # -- Load or simulate data -----------------------------------------------
    gene_data = cpd_data = None

    if args.simulate:
        print(f"Info: Generating simulated gene data (n={args.n_sim}) …")
        gene_data = sim_mol_data(
            mol_type="gene",
            species=args.species,
            n_mol=args.n_sim,
        )
    else:
        if args.gene_data:
            gene_data = pl.read_csv(args.gene_data, separator="\t")
            print(f"Info: Loaded gene data — {gene_data.height} rows, "
                  f"{gene_data.width - 1} experiment column(s).")
        if args.cpd_data:
            cpd_data = pl.read_csv(args.cpd_data, separator="\t")
            print(f"Info: Loaded compound data — {cpd_data.height} rows, "
                  f"{cpd_data.width - 1} experiment column(s).")

    if gene_data is None and cpd_data is None:
        _build_parser().error(
            "Provide at least one of --gene-data, --cpd-data, or --simulate."
        )

    # -- Run pathview --------------------------------------------------------
    result = pathview(
        pathway_id    = args.pathway_id,
        gene_data     = gene_data,
        cpd_data      = cpd_data,
        species       = args.species,
        kegg_dir      = args.kegg_dir,
        kegg_native   = args.kegg_native,
        output_format = args.output_format,
        gene_idtype   = args.gene_idtype,
        cpd_idtype    = args.cpd_idtype,
        out_suffix    = args.out_suffix,
        node_sum      = args.node_sum,
        map_symbol    = args.map_symbol,
        min_nnodes    = args.min_nnodes,
        new_signature = not args.no_signature,
        plot_col_key  = not args.no_col_key,
        limit         = {"gene": args.limit_gene, "cpd": args.limit_cpd},
        bins          = {"gene": args.bins_gene,  "cpd": args.bins_cpd},
        both_dirs     = {"gene": True,            "cpd": True},
        low           = {"gene": args.low_gene,   "cpd": args.low_cpd},
        mid           = {"gene": args.mid_gene,   "cpd": args.mid_cpd},
        high          = {"gene": args.high_gene,  "cpd": args.high_cpd},
    )

    if result:
        gdf = result.get("plot_data_gene")
        cdf = result.get("plot_data_cpd")
        if gdf is not None:
            print(f"Info: Gene plot data — {gdf.height} nodes mapped.")
        if cdf is not None:
            print(f"Info: Compound plot data — {cdf.height} nodes mapped.")
    else:
        print("Warning: Pathway was skipped or no data could be mapped.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
