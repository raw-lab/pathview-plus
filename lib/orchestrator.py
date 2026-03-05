"""
orchestrator.py  (module inside the pathview package)
Core orchestrator: resolves IDs, downloads KEGG files, maps data to nodes,
and dispatches to the appropriate renderer.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional

import polars as pl

from .color_mapping import node_color
from .constants     import SumMethod, VALID_NODE_TYPES, NODE_META_COLS
from .id_mapping    import cpd_id_map, eg2id, id2eg
from .kegg_api      import download_kegg, kegg_species_code
from .kgml_parser   import node_info, parse_kgml
from .mol_data      import mol_sum
from .node_mapping  import node_map
from .rendering     import keggview_graph, keggview_native
from .svg_rendering import keggview_svg


# ---------------------------------------------------------------------------
# Defaults factory  (avoids mutable default arguments)
# ---------------------------------------------------------------------------

def _defaults() -> dict:
    return dict(
        limit     = {"gene": 1.0,    "cpd": 1.0},
        bins      = {"gene": 10,     "cpd": 10},
        both_dirs = {"gene": True,   "cpd": True},
        discrete  = {"gene": False,  "cpd": False},
        low       = {"gene": "green","cpd": "blue"},
        mid       = {"gene": "gray", "cpd": "gray"},
        high      = {"gene": "red",  "cpd": "yellow"},
        trans_fun = {"gene": None,   "cpd": None},
    )


# ---------------------------------------------------------------------------
# pathview
# ---------------------------------------------------------------------------

def pathview(
    pathway_id: str,
    gene_data: Optional[pl.DataFrame] = None,
    cpd_data: Optional[pl.DataFrame] = None,
    species: str = "hsa",
    kegg_dir: str | Path = ".",
    kegg_native: bool = True,
    output_format: str = "png",  # NEW: png, pdf, or svg
    gene_idtype: str = "ENTREZ",
    cpd_idtype: str = "KEGG",
    out_suffix: str = "pathview",
    node_sum: SumMethod = "sum",
    map_symbol: bool = True,
    map_null: bool = True,
    min_nnodes: int = 3,
    new_signature: bool = True,
    plot_col_key: bool = True,
    # Colour-scale parameters (all accept {"gene": …, "cpd": …} dicts)
    limit: dict | None = None,
    bins: dict | None = None,
    both_dirs: dict | None = None,
    discrete: dict | None = None,
    low: dict | None = None,
    mid: dict | None = None,
    high: dict | None = None,
    na_col: str = "transparent",
    trans_fun: dict | None = None,
    **kwargs,
) -> dict:
    """
    Overlay molecular data onto a KEGG pathway diagram.

    Parameters
    ----------
    pathway_id:    KEGG pathway number, e.g. ``"04110"`` or ``"hsa04110"``.
    gene_data:     DataFrame — first column = gene IDs, rest = numeric values.
    cpd_data:      DataFrame — first column = compound IDs, rest = numeric.
    species:       KEGG species code (default ``"hsa"``).
    kegg_dir:      Working directory for downloaded and output files.
    kegg_native:   True → overlay on KEGG PNG; False → NetworkX graph layout.
    gene_idtype:   Input gene ID type (``"ENTREZ"``, ``"SYMBOL"``, ``"KEGG"``…).
    cpd_idtype:    Input compound ID type (``"KEGG"``, ``"PUBCHEM"``…).
    out_suffix:    Suffix for output filenames.
    node_sum:      Aggregation method for multi-probe nodes.
    map_symbol:    Replace Entrez IDs with gene symbols in node labels.
    map_null:      Render nodes even when no data is provided.
    min_nnodes:    Skip pathway if fewer than this many mappable nodes exist.
    new_signature: Add a "Rendered by pathview.py" watermark.
    plot_col_key:  Draw the colour-scale legend.
    limit/bins/both_dirs/discrete/low/mid/high/trans_fun:
                   Colour-scale parameters, each a dict with "gene" and "cpd"
                   keys.
    na_col:        Colour for unmapped nodes (default ``"transparent"``).

    Returns
    -------
    dict with keys ``"plot_data_gene"`` and ``"plot_data_cpd"`` (Polars
    DataFrames), or an empty dict when the pathway could not be processed.

    Examples
    --------
    >>> import polars as pl
    >>> from pathview import pathview
    >>> gene_df = pl.read_csv("gene_expr.tsv", separator="\\t")
    >>> result  = pathview("04110", gene_data=gene_df, species="hsa")
    """
    if gene_data is None and cpd_data is None:
        raise ValueError("At least one of gene_data or cpd_data must be provided.")

    # Merge caller-supplied dicts over defaults
    cfg = _defaults()
    for key, val in dict(
        limit=limit, bins=bins, both_dirs=both_dirs, discrete=discrete,
        low=low, mid=mid, high=high, trans_fun=trans_fun,
    ).items():
        if val is not None:
            cfg[key] = val

    kegg_dir = Path(kegg_dir)

    # ---- Species resolution ------------------------------------------------
    species_info = kegg_species_code(species)
    species_code = species_info.kegg_code
    if species_code == "ko":
        gene_idtype = "KEGG"

    # ---- Normalise pathway ID ----------------------------------------------
    pathway_name = (
        pathway_id if pathway_id.startswith(species_code)
        else f"{species_code}{pathway_id}"
    )
    numeric_id = pathway_name.replace(species_code, "")

    # ---- Gene ID conversion ------------------------------------------------
    if gene_data is not None:
        gene_data = _maybe_convert_gene_ids(
            gene_data, gene_idtype, species_code, node_sum
        )

    # ---- Compound ID conversion --------------------------------------------
    if cpd_data is not None and "kegg" not in cpd_idtype.lower():
        cpd_data = _maybe_convert_cpd_ids(cpd_data, cpd_idtype, node_sum)

    # ---- Download missing files --------------------------------------------
    needed = ["xml", "png"] if kegg_native else ["xml"]
    existing = {f.name for f in kegg_dir.iterdir()} if kegg_dir.exists() else set()
    missing = [t for t in needed if f"{pathway_name}.{t}" not in existing]

    if missing:
        status = download_kegg(numeric_id, species=species_code,
                               kegg_dir=kegg_dir, file_type=missing)
        if status.get(pathway_name) == "failed":
            warnings.warn(f"Failed to download files for {pathway_name}; skipping.")
            return {}

    # ---- Parse KGML --------------------------------------------------------
    pathway  = parse_kgml(kegg_dir / f"{pathway_name}.xml")
    node_data = (
        node_info(pathway)
        .filter(
            pl.col("type").is_in(VALID_NODE_TYPES)
            & pl.col("x").is_not_null()
            & pl.col("y").is_not_null()
        )
    )

    if node_data.height < min_nnodes:
        warnings.warn(
            f"Only {node_data.height} mappable nodes for {pathway_name} "
            f"(minimum {min_nnodes}); skipping."
        )
        return {}

    # ---- Map gene data onto nodes ------------------------------------------
    gene_node_type = "ortholog" if species_code == "ko" else "gene"
    plot_data_gene, cols_gene = _map_and_color(
        mol_data=gene_data,
        node_data=node_data,
        node_types=gene_node_type,
        node_sum=node_sum,
        map_null=map_null,
        color_cfg={k: cfg[k]["gene"] for k in ("limit","bins","both_dirs","discrete","low","mid","high","trans_fun")},
        na_col=na_col,
    )

    # Optionally replace Entrez IDs with gene symbols in labels
    if plot_data_gene is not None and map_symbol and gene_data is not None:
        plot_data_gene = _add_symbol_labels(plot_data_gene, species_code)

    # ---- Map compound data onto nodes --------------------------------------
    plot_data_cpd, cols_cpd = _map_and_color(
        mol_data=cpd_data,
        node_data=node_data,
        node_types="compound",
        node_sum=node_sum,
        map_null=map_null,
        color_cfg={k: cfg[k]["cpd"] for k in ("limit","bins","both_dirs","discrete","low","mid","high","trans_fun")},
        na_col=na_col,
    )

    # ---- Render ------------------------------------------------------------
    render_kwargs = dict(
        plot_data_gene=plot_data_gene, cols_gene=cols_gene,
        plot_data_cpd=plot_data_cpd,  cols_cpd=cols_cpd,
        node_data=node_data,
        pathway_name=pathway_name,
        kegg_dir=kegg_dir,
        out_suffix=out_suffix,
        new_signature=new_signature,
        plot_col_key=plot_col_key,
        **{k: cfg[k] for k in ("limit","bins","both_dirs","discrete","low","mid","high")},
    )

    if output_format == "svg":
        # SVG vector output (works for both KEGG and SBGN)
        keggview_svg(**{k: v for k, v in render_kwargs.items()
                        if k not in ("discrete", "plot_col_key")})
    elif kegg_native and output_format == "png":
        # PNG with KEGG background (only for KEGG pathways)
        keggview_native(**render_kwargs)
    else:
        # PDF graph layout (works for both KEGG and SBGN)
        keggview_graph(**{k: v for k, v in render_kwargs.items()
                          if k not in ("discrete",)}, **kwargs)

    return {"plot_data_gene": plot_data_gene, "plot_data_cpd": plot_data_cpd}


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _maybe_convert_gene_ids(
    gene_data: pl.DataFrame,
    gene_idtype: str,
    species_code: str,
    node_sum: SumMethod,
) -> pl.DataFrame:
    """Convert non-Entrez gene IDs to Entrez before pathway mapping."""
    if gene_idtype.upper() in ("ENTREZ", "ENTREZID", "KEGG"):
        return gene_data
    id_col   = gene_data.columns[0]
    id_map   = id2eg(gene_data[id_col].to_list(), category=gene_idtype, org=species_code)
    return mol_sum(gene_data, id_map, sum_method=node_sum)


def _maybe_convert_cpd_ids(
    cpd_data: pl.DataFrame,
    cpd_idtype: str,
    node_sum: SumMethod,
) -> pl.DataFrame:
    """Convert non-KEGG compound IDs to KEGG before pathway mapping."""
    id_col   = cpd_data.columns[0]
    id_map   = cpd_id_map(cpd_data[id_col].to_list(), in_type=cpd_idtype, out_type="KEGG")
    return mol_sum(cpd_data, id_map, sum_method=node_sum)


def _map_and_color(
    mol_data: Optional[pl.DataFrame],
    node_data: pl.DataFrame,
    node_types: str,
    node_sum: SumMethod,
    map_null: bool,
    color_cfg: dict,
    na_col: str,
) -> tuple[Optional[pl.DataFrame], Optional[pl.DataFrame]]:
    """
    Map molecule data to nodes then compute per-node colours.

    Returns (plot_data, cols) where cols is a DataFrame of hex colour strings,
    or (None, None) when no nodes of the requested type exist.
    """
    if mol_data is None and not map_null:
        return None, None

    plot_data = node_map(mol_data, node_data, node_types=node_types, node_sum=node_sum)
    if plot_data is None:
        return None, None

    val_cols = [c for c in plot_data.columns if c not in NODE_META_COLS]
    if not val_cols:
        return plot_data, None

    cols = node_color(
        plot_data.select(["entry_id"] + val_cols).rename({"entry_id": "id"}),
        limit    = color_cfg["limit"],
        bins     = color_cfg["bins"],
        both_dirs= color_cfg["both_dirs"],
        discrete = color_cfg["discrete"],
        low      = color_cfg["low"],
        mid      = color_cfg["mid"],
        high     = color_cfg["high"],
        na_col   = na_col,
        trans_fun= color_cfg["trans_fun"],
    )
    return plot_data, cols


def _add_symbol_labels(
    plot_data: pl.DataFrame,
    species_code: str,
) -> pl.DataFrame:
    """Attempt to replace Entrez-based labels with gene symbols."""
    try:
        gene_ids = plot_data["kegg_names"].drop_nulls().to_list()
        sym_map  = eg2id(gene_ids, category="SYMBOL", org=species_code)
        return plot_data.join(sym_map, left_on="kegg_names", right_on="ENTREZID", how="left")
    except Exception as exc:
        warnings.warn(f"Symbol label mapping failed: {exc}")
        return plot_data
