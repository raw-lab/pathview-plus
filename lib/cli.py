"""
cli.py
Command-line interface for pathview-plus.

Fixes over v2.x
---------------
The old ``bin/pathview-cli.py`` cast ``gene_data`` unconditionally, so
``--cpd-data`` on its own raised ``AttributeError: 'NoneType' object has no
attribute 'cast'`` — a metabolomics-only run, the obvious use case for a
metabolic map, could not be done from the command line at all.  It also
exposed neither render modes, themes, palettes, nor offline operation.

Entry point
-----------
Installed as ``pathview-plus`` (and ``pathview-cli`` for continuity).

    pathview-plus render 00020 --species human \\
        --gene-data rnaseq.csv --cpd-data metabolites.csv \\
        --render-mode vector --output-format pdf

    pathview-plus species mouse
    pathview-plus search coli
    pathview-plus parity --markdown
    pathview-plus legend --out legend.png
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

import polars as pl

# ---------------------------------------------------------------------------
# Input helpers
# ---------------------------------------------------------------------------

def _read_table(path: str | Path, label: str) -> pl.DataFrame:
    """
    Read a delimited data file, inferring the separator from the extension.

    The first column is treated as identifiers; every other column is coerced
    to Float64 so a stray text column cannot silently become "data".
    """
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"error: {label} file not found: {p}")

    sep = {".tsv": "\t", ".txt": "\t", ".tab": "\t"}.get(p.suffix.lower(), ",")
    try:
        # infer_schema_length=0 reads every column as text.  Type inference on
        # the identifier column is actively harmful: an expression matrix whose
        # first 2,000 gene ids are numeric but which later contains a probe id
        # such as "AFFX-BioB-3" makes the whole read fail.  Value columns are
        # cast explicitly below instead.
        df = pl.read_csv(p, separator=sep, infer_schema_length=0)
    except Exception as exc:
        raise SystemExit(f"error: could not read {label} file {p}: {exc}") from exc

    if df.width < 2:
        raise SystemExit(
            f"error: {label} file {p.name} has {df.width} column(s); an "
            "identifier column plus at least one value column is required."
        )

    id_col = df.columns[0]
    value_cols = df.columns[1:]
    df = df.with_columns(
        [pl.col(id_col).cast(pl.String)]
        + [pl.col(c).cast(pl.Float64, strict=False) for c in value_cols]
    )

    usable = [c for c in value_cols if df[c].null_count() < df.height]
    if not usable:
        raise SystemExit(
            f"error: no numeric value column in {p.name}. Columns after the "
            f"first must be numeric; found {value_cols}."
        )
    if len(usable) < len(value_cols):
        dropped = set(value_cols) - set(usable)
        print(f"[pathview] ignoring non-numeric column(s) in {p.name}: "
              f"{', '.join(sorted(dropped))}", file=sys.stderr)

    return df.select([id_col] + usable)


def _parse_limit(text: str | None) -> float | dict | None:
    """``--limit 1.5`` or ``--limit gene=2,cpd=1``."""
    if not text:
        return None
    if "=" not in text:
        return float(text)
    out: dict[str, float] = {}
    for part in text.split(","):
        key, _, val = part.partition("=")
        key = key.strip().lower()
        if key not in ("gene", "cpd"):
            raise SystemExit(
                f"error: --limit key must be gene or cpd, got {key!r}")
        out[key] = float(val)
    return out


# ---------------------------------------------------------------------------
# Sub-commands
# ---------------------------------------------------------------------------

def cmd_render(args: argparse.Namespace) -> int:
    from pathview.cache import set_offline
    from pathview.errors import PathviewError

    from pathview import pathview

    if args.offline:
        set_offline(True)

    gene_data = _read_table(args.gene_data, "--gene-data") if args.gene_data else None
    # Guarded independently: v2.x cast gene_data unconditionally here.
    cpd_data = _read_table(args.cpd_data, "--cpd-data") if args.cpd_data else None

    if gene_data is None and cpd_data is None and not args.map_null:
        raise SystemExit(
            "error: supply --gene-data, --cpd-data, or --map-null to draw an "
            "unmapped diagram."
        )

    failures = 0
    for pathway_id in args.pathway_id:
        try:
            res = pathview(
                pathway_id,
                gene_data=gene_data, cpd_data=cpd_data, species=args.species,
                gene_idtype=args.gene_idtype, cpd_idtype=args.cpd_idtype,
                kegg_dir=args.kegg_dir, out_dir=args.out_dir,
                out_suffix=args.out_suffix,
                render_mode=args.render_mode, output_format=args.output_format,
                theme=args.theme, title=args.title, subtitle=args.subtitle,
                figure_width=args.figure_width, dpi=args.dpi,
                gene_color=args.gene_palette, cpd_color=args.cpd_palette,
                limit=_parse_limit(args.limit), bins=args.bins,
                both_dirs=not args.one_sided, discrete=args.discrete,
                node_sum=args.node_sum, map_null=args.map_null,
                split_group=args.split_group, expand_node=args.expand_node,
                min_nnodes=args.min_nodes, rand_seed=args.seed,
                draw_edges=not args.no_edges,
                show_link_edges=args.link_edges,
                plot_col_key=not args.no_key,
                quiet=args.quiet,
            )
            if not args.quiet:
                print(res.summary())
            else:
                print(res.output_path)
        except PathviewError as exc:
            print(f"error [{pathway_id}]: {exc}", file=sys.stderr)
            failures += 1
    return 1 if failures else 0


def cmd_species(args: argparse.Namespace) -> int:
    from pathview.errors import SpeciesNotFoundError

    from pathview import get_species_code

    try:
        info = get_species_code(args.query)
    except SpeciesNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    rows = [("KEGG code", info.kegg_code),
            ("Scientific name", info.scientific_name),
            ("Common name", info.common_name or "-"),
            ("NCBI taxonomy", info.tax_id or "-"),
            ("KEGG T-number", info.ktax_id or "-"),
            ("Gene id type", "ENTREZ" if info.entrez_gnodes else "KEGG")]
    width = max(len(k) for k, _ in rows)
    for key, val in rows:
        print(f"{key:<{width}}  {val}")
    return 0


def cmd_search(args: argparse.Namespace) -> int:
    from pathview import search_organisms

    hits = search_organisms(args.query, limit=args.limit)
    if hits.is_empty():
        print(f"No organism matched {args.query!r}.", file=sys.stderr)
        return 1
    for row in hits.iter_rows(named=True):
        common = f"  ({row['common_name']})" if row.get("common_name") else ""
        print(f"{row['kegg_code']:<8} {row['scientific_name']}{common}")
    return 0


def cmd_parity(args: argparse.Namespace) -> int:
    from pathview import parity_summary, print_parity

    print(print_parity(category=args.category, markdown=args.markdown))
    if not args.category:
        s = parity_summary()
        print(f"\n{s['total_features']} features tracked: "
              f"{s['full']} full, {s['partial']} partial, {s['none']} missing.")
        print(f"vs pathview (R): {s['vs_pathview_R']} ({s['vs_pathview_R_pct']}%)")
        print(f"vs SBGNview (R): {s['vs_SBGNview_R']} ({s['vs_SBGNview_R_pct']}%)")
    return 0


def cmd_legend(args: argparse.Namespace) -> int:
    from pathview import kegg_legend, sbgn_legend

    out = Path(args.out)
    if args.standard == "sbgn":
        sbgn_legend(out_path=out, theme=args.theme, dpi=args.dpi)
    else:
        kegg_legend(legend_type=args.legend_type, out_path=out,
                    theme=args.theme, dpi=args.dpi)
    print(out)
    return 0


def cmd_download(args: argparse.Namespace) -> int:
    from pathview.errors import PathviewError

    from pathview import download_pathway

    failures = 0
    for pid in args.pathway_id:
        try:
            result = download_pathway(pid, output_dir=args.out_dir,
                                      species=args.species)
            print(f"{pid}: {result}")
        except PathviewError as exc:
            print(f"error [{pid}]: {exc}", file=sys.stderr)
            failures += 1
    return 1 if failures else 0


def cmd_sbgn(args: argparse.Namespace) -> int:
    """Browse the pre-generated SBGN collection, or render from it."""
    from pathview.cache import set_offline
    from pathview.errors import PathviewError

    from pathview import list_sbgn_pathways, sbgn_collection_info, sbgnview

    if args.list or not args.pathway_id:
        info = sbgn_collection_info()
        if not args.source and not args.query:
            print(f"{info['total']:,} pre-generated SBGN pathways:")
            for src, n in sorted(info["by_source"].items(),
                                 key=lambda kv: -kv[1]):
                print(f"  {info['sources'][src]:<42} {n:>6,}")
            return 0
        hits = list_sbgn_pathways(args.source, args.query, args.max_results)
        if hits.is_empty():
            print("No pathway matched.", file=sys.stderr)
            return 1
        for row in hits.iter_rows(named=True):
            print(f"{row['source']:<10} {row['pathway_id']}")
        return 0

    if args.offline:
        set_offline(True)
    gene_data = _read_table(args.gene_data, "--gene-data") if args.gene_data else None
    cpd_data = _read_table(args.cpd_data, "--cpd-data") if args.cpd_data else None

    try:
        res = sbgnview(
            args.pathway_id if len(args.pathway_id) > 1 else args.pathway_id[0],
            gene_data=gene_data, cpd_data=cpd_data,
            gene_idtype=args.gene_idtype, cpd_idtype=args.cpd_idtype,
            sbgn_dir=args.sbgn_dir, out_dir=args.out_dir,
            output_format=args.output_format, theme=args.theme,
            show_compartments=not args.no_compartments,
            limit=_parse_limit(args.limit), quiet=args.quiet,
        )
    except PathviewError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(res.summary())
    return 0


def cmd_info(args: argparse.Namespace) -> int:
    from pathview.cache import cache_dir

    import pathview
    from pathview import (
        DATABASE_INFO,
        list_palettes,
        organism_count,
        supported_cpd_idtypes,
        supported_gene_idtypes,
    )

    print(f"pathview-plus {pathview.__version__}")
    print(f"organisms bundled : {organism_count():,}")
    print(f"cache directory   : {cache_dir()}")
    print(f"palettes          : {', '.join(sorted(list_palettes()))}")
    print(f"gene id types     : {', '.join(supported_gene_idtypes()[:12])} ...")
    print(f"compound id types : {', '.join(supported_cpd_idtypes()[:12])} ...")
    from pathview import sbgn_collection_info
    ci = sbgn_collection_info()
    print(f"SBGN collection   : {ci['total']:,} pathways indexed offline")
    print("\nPathway sources:")
    for meta in DATABASE_INFO.values():
        mark = "available" if meta["available"] else "manual only"
        print(f"  {meta['name']:<18} {mark:<12} e.g. {meta['example']:<14} "
              f"({meta['source']})")
    return 0


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    import pathview

    parser = argparse.ArgumentParser(
        prog="pathview-plus",
        description="Overlay omics data on KEGG and SBGN pathway diagrams.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  pathview-plus render 00020 --species human --gene-data rna.csv\n"
            "  pathview-plus render 00020 --cpd-data metabolites.csv "
            "--cpd-idtype NAME\n"
            "  pathview-plus render 04110 04010 --gene-data rna.csv "
            "--render-mode vector --output-format pdf\n"
            "  pathview-plus species 'Mus musculus'\n"
            "  pathview-plus parity --markdown > PARITY.md\n"
        ),
    )
    parser.add_argument("--version", action="version",
                        version=f"pathview-plus {pathview.__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    # -- render -------------------------------------------------------------
    r = sub.add_parser("render", help="render one or more pathways")
    r.add_argument("pathway_id", nargs="+",
                   help="KEGG id(s), e.g. 00020 or hsa00020")
    r.add_argument("--species", "-s", default="hsa",
                   help="code, name or taxid (hsa, human, 'Homo sapiens', 9606)")
    r.add_argument("--gene-data", help="CSV/TSV: ids in column 1, values after")
    r.add_argument("--cpd-data", help="CSV/TSV: ids in column 1, values after")
    r.add_argument("--gene-idtype", default="ENTREZ",
                   help="ENTREZ, SYMBOL, ENSEMBL, UNIPROT, REFSEQ, KEGG")
    r.add_argument("--cpd-idtype", default="KEGG",
                   help="KEGG, NAME, CAS, CHEBI, HMDB, PUBCHEM, DRUGBANK")
    r.add_argument("--kegg-dir", default=".", help="KGML/PNG cache directory")
    r.add_argument("--out-dir", "-o", default=".", help="output directory")
    r.add_argument("--out-suffix", default="pathview")
    r.add_argument("--render-mode", default="auto",
                   choices=["auto", "native", "vector", "graph", "svg"])
    r.add_argument("--output-format", "-f", default="png",
                   choices=["png", "pdf", "svg"])
    r.add_argument("--theme", default="publication",
                   choices=["publication", "slate", "dark"])
    r.add_argument("--gene-palette", default=None,
                   help="palette name for transcript data")
    r.add_argument("--cpd-palette", default=None,
                   help="palette name for metabolite data")
    r.add_argument("--limit", default=None,
                   help="scale limit: '1.5' or 'gene=2,cpd=1'")
    r.add_argument("--bins", type=int, default=10)
    r.add_argument("--one-sided", action="store_true",
                   help="map 0..limit instead of -limit..+limit")
    r.add_argument("--discrete", action="store_true")
    r.add_argument("--node-sum", default="sum",
                   choices=["sum", "mean", "median", "max", "min",
                            "max_abs", "random", "first"])
    r.add_argument("--split-group", action="store_true",
                   help="replace complexes with their subunits")
    r.add_argument("--expand-node", action="store_true",
                   help="give each gene on a multi-gene node its own sub-node")
    r.add_argument("--map-null", action="store_true",
                   help="draw the diagram even with no data")
    r.add_argument("--min-nodes", type=int, default=3)
    r.add_argument("--seed", type=int, default=None,
                   help="seed for --node-sum random")
    r.add_argument("--title", default=None)
    r.add_argument("--subtitle", default=None)
    r.add_argument("--figure-width", type=float, default=14.0)
    r.add_argument("--dpi", type=int, default=220)
    r.add_argument("--no-edges", action="store_true")
    r.add_argument("--link-edges", action="store_true",
                   help="also draw edges to pathway-link nodes")
    r.add_argument("--no-key", action="store_true", help="omit colour keys")
    r.add_argument("--offline", action="store_true",
                   help="never attempt a network request")
    r.add_argument("--quiet", "-q", action="store_true")
    r.set_defaults(func=cmd_render)

    # -- species / search ---------------------------------------------------
    sp = sub.add_parser("species", help="resolve a species to its KEGG code")
    sp.add_argument("query")
    sp.set_defaults(func=cmd_species)

    se = sub.add_parser("search", help="search the bundled organism table")
    se.add_argument("query")
    se.add_argument("--limit", "-n", type=int, default=10)
    se.set_defaults(func=cmd_search)

    # -- parity -------------------------------------------------------------
    pa = sub.add_parser("parity", help="feature matrix vs the R packages")
    pa.add_argument("--category", default=None)
    pa.add_argument("--markdown", action="store_true")
    pa.set_defaults(func=cmd_parity)

    # -- legend -------------------------------------------------------------
    lg = sub.add_parser("legend", help="write a diagram legend")
    lg.add_argument("--out", default="legend.png")
    lg.add_argument("--standard", default="kegg", choices=["kegg", "sbgn"])
    lg.add_argument("--legend-type", default="both",
                    choices=["both", "edge", "node"])
    lg.add_argument("--theme", default="publication")
    lg.add_argument("--dpi", type=int, default=200)
    lg.set_defaults(func=cmd_legend)

    # -- download -----------------------------------------------------------
    dl = sub.add_parser("download", help="fetch pathway files without rendering")
    dl.add_argument("pathway_id", nargs="+")
    dl.add_argument("--species", "-s", default="hsa")
    dl.add_argument("--out-dir", "-o", default=".")
    dl.set_defaults(func=cmd_download)

    # -- sbgn ---------------------------------------------------------------
    sb = sub.add_parser("sbgn", help="browse or render the SBGN collection")
    sb.add_argument("pathway_id", nargs="*",
                    help="collection id(s) or path(s) to local .sbgn files")
    sb.add_argument("--list", action="store_true", help="list rather than render")
    sb.add_argument("--source", default=None,
                    choices=["reactome", "smpdb", "panther", "metacyc", "metacrop"])
    sb.add_argument("--query", default=None, help="substring filter on the id")
    # -n caps the listing; --limit is the colour-scale limit, matching
    # `render`.  Sharing one flag for both made the scale limit receive the
    # listing count.
    sb.add_argument("--max", "-n", type=int, default=40, dest="max_results",
                    help="maximum rows when listing")
    sb.add_argument("--limit", default=None,
                    help="colour scale limit: '1.5' or 'gene=2,cpd=1'")
    sb.add_argument("--gene-data")
    sb.add_argument("--cpd-data")
    sb.add_argument("--gene-idtype", default="SYMBOL",
                    help="SYMBOL gives the densest offline SBGN coverage")
    sb.add_argument("--cpd-idtype", default="KEGG")
    sb.add_argument("--sbgn-dir", default=".", help="SBGN cache directory")
    sb.add_argument("--out-dir", "-o", default=".")
    sb.add_argument("--output-format", "-f", default="png",
                    choices=["png", "pdf", "svg"])
    sb.add_argument("--theme", default="publication",
                    choices=["publication", "slate", "dark"])
    sb.add_argument("--no-compartments", action="store_true")
    sb.add_argument("--offline", action="store_true")
    sb.add_argument("--quiet", "-q", action="store_true")
    sb.set_defaults(func=cmd_sbgn)

    # -- info ---------------------------------------------------------------
    nf = sub.add_parser("info", help="show installation and capability summary")
    nf.set_defaults(func=cmd_info)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
