"""
sbgnview.py
The SBGN counterpart to :func:`pathview.pathview`.

``sbgnview()`` is to SBGN maps what ``pathview()`` is to KEGG maps: give it a
pathway identifier (or a local SBGN file) and some omics data, and it
resolves the map, maps the data onto the glyphs, and renders it.

The hard part of SBGN is not drawing — it is that glyph identifiers are
opaque.  A Pathway Commons macromolecule is called
``Protein_03a9d8039cd87a8c55de7405670d4682``; no user has that in their count
matrix.  :mod:`pathview.sbgn_hub` ships the crosswalks that bridge Entrez,
gene symbol, KO, KEGG compound, ChEBI and compound name to those glyph ids,
so this works offline once the map itself is cached.

Compared with R SBGNview
------------------------
Covered: the pre-generated collection, glyph-id crosswalks, compartment
shading, multi-condition slicing, arc styling by class, highlighting, and
rendering to PNG/PDF/SVG.

Public API
----------
  sbgnview          : render one SBGN map with data overlaid
  sbgnview_batch    : render several
  sbgn_node_map     : map data onto SBGN glyphs (used internally, useful alone)
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import polars as pl

from .cache import is_offline
from .color_mapping import ColorScale
from .constants import DEFAULT_CPD_PALETTE, DEFAULT_GENE_PALETTE, SumMethod
from .errors import MappingError, PathviewError, PathwayNotFoundError
from .highlighting import PathwayResult, PathwayResultSet
from .mol_data import compound_name, mol_sum
from .sbgn_hub import _canon, download_sbgn, map_ids_to_sbgn
from .sbgn_parser import parse_sbgn, sbgn_compartments, sbgn_edges, sbgn_to_df

#: Glyph types that carry gene/protein data, and those that carry metabolites.
_GENE_TYPES = ("gene",)
_CPD_TYPES = ("compound",)


def _resolve_source(
    pathway_id: str | Path,
    sbgn_dir: str | Path,
    overwrite: bool = False,
) -> Path:
    """
    Turn a pathway id or file path into a local SBGN file.

    A path that exists is used as-is — an SBGN file exported by hand from
    PANTHER, MetaCyc or SMPDB parses exactly like one from the collection, so
    there is no second code path for "unsupported" sources.
    """
    candidate = Path(pathway_id)
    if candidate.exists() and candidate.is_file():
        return candidate

    local = Path(sbgn_dir) / f"{pathway_id}.sbgn"
    if local.exists() and local.stat().st_size > 0 and not overwrite:
        return local

    return download_sbgn(str(pathway_id), output_dir=sbgn_dir, overwrite=overwrite)


def sbgn_node_map(
    mol_data: pl.DataFrame | None,
    node_data: pl.DataFrame,
    id_type: str = "ENTREZ",
    node_types: Sequence[str] = _GENE_TYPES,
    node_sum: SumMethod = "sum",
    rand_seed: int | None = None,
    detailed: bool = False,
):
    """
    Map molecular data onto SBGN glyphs via the bundled crosswalks.

    User identifiers are translated to the identifier system the glyphs
    actually use, then aggregated onto glyph ids.
    """
    from .node_mapping import NodeMapResult

    targets = node_data.filter(pl.col("type").is_in(list(node_types)))
    if targets.is_empty():
        return NodeMapResult(None) if detailed else None
    if mol_data is None:
        res = NodeMapResult(data=targets, n_nodes=targets.height)
        return res if detailed else targets

    # Which identifier systems do this map's glyphs actually use?
    glyph_ids: set[str] = set()
    for names in targets["kegg_names"].to_list():
        glyph_ids.update(str(n) for n in (names or []))

    id_col = mol_data.columns[0]
    user_ids = [str(v) for v in mol_data[id_col].to_list()]

    # Try every identifier system the glyphs might use, and keep whichever
    # actually lands on this map.  Pathway Commons hashes, ChEBI ids and bare
    # gene symbols all appear as glyph ids depending on the source database,
    # so guessing one would fail silently on the others.
    candidates: list[tuple[str, list[str]]] = [(id_type, user_ids)]

    # The offline Entrez route runs through KO, which is sparse; gene symbol
    # has a direct 294k-pair route to Pathway Commons.  When a conversion
    # service is reachable, going Entrez -> symbol first is far denser.
    if _canon(id_type) not in ("symbol",) and not is_offline():
        try:
            import warnings

            from .id_mapping import eg2id
            with warnings.catch_warnings():
                # This detour is an optimisation, not a requirement: if the
                # service is unreachable the offline route still works, so a
                # warning here would be noise rather than information.
                warnings.simplefilter("ignore")
                conv = eg2id(user_ids, category="SYMBOL", org="hsa")
            syms = conv["SYMBOL"].to_list()
            if any(syms):
                candidates.append(("SYMBOL", [s or "" for s in syms]))
        except Exception:
            pass

    best: pl.DataFrame | None = None
    best_hits = 0
    for src_type, src_ids in candidates:
        for target_type in ("pathwaycommons", "chebi", "kegg", "symbol",
                            "entrez", "ko"):
            try:
                conv = map_ids_to_sbgn(src_ids, src_type, target_type)
            except ValueError:
                continue
            if conv.width < 2:
                continue
            col = conv.columns[1]
            hits = sum(1 for v in conv[col].to_list() if v and v in glyph_ids)
            if hits > best_hits:
                # Re-key onto the user's own identifiers so the join below
                # still works after a symbol detour.
                conv = conv.with_columns(pl.Series(conv.columns[0], user_ids))
                best, best_hits = conv, hits

    if best is None or best_hits == 0:
        res = NodeMapResult(None, n_nodes=targets.height, n_ids_input=len(user_ids))
        return res if detailed else None

    # glyph identifier -> entry_id
    exploded = targets.explode("kegg_names").filter(pl.col("kegg_names").is_not_null())
    conv_col = best.columns[1]
    bridge = (best.rename({best.columns[0]: id_col, conv_col: "kegg_names"})
                  .drop_nulls()
                  .join(exploded.select(["kegg_names", "entry_id"]), on="kegg_names")
                  .select([pl.col(id_col).cast(pl.String), pl.col("entry_id")])
                  .unique())
    if bridge.is_empty():
        res = NodeMapResult(None, n_nodes=targets.height, n_ids_input=len(user_ids))
        return res if detailed else None

    summed = mol_sum(mol_data, bridge, sum_method=node_sum,
                     rand_seed=rand_seed, detailed=True)
    values = summed.data.rename({summed.data.columns[0]: "entry_id"})
    plot_data = targets.join(values, on="entry_id", how="left")

    value_cols = [c for c in values.columns if c != "entry_id"]
    with_data = plot_data.filter(
        pl.any_horizontal([pl.col(c).is_not_null() for c in value_cols])
    ).height if value_cols else 0

    res = NodeMapResult(
        data=plot_data, n_nodes=targets.height, n_nodes_with_data=with_data,
        n_ids_input=len(user_ids), n_ids_mapped=summed.n_mapped,
        unmapped_ids=summed.unmapped_ids, value_columns=value_cols,
    )
    return res if detailed else plot_data


def sbgnview(
    pathway_id: str | Path | Sequence[str],
    gene_data: pl.DataFrame | None = None,
    cpd_data: pl.DataFrame | None = None,
    *,
    gene_idtype: str = "ENTREZ",
    cpd_idtype: str = "KEGG",
    sbgn_dir: str | Path = ".",
    out_dir: str | Path | None = None,
    out_suffix: str = "sbgnview",
    output_format: str = "png",
    theme: str = "publication",
    title: str | None = None,
    subtitle: str | None = None,
    show_compartments: bool = True,
    show_processes: bool = True,
    draw_edges: bool = True,
    gene_color: ColorScale | str | dict | None = None,
    cpd_color: ColorScale | str | dict | None = None,
    limit: float | dict | None = None,
    bins: int = 10,
    node_sum: SumMethod = "sum",
    map_cpd_name: bool = True,
    figure_width: float = 15.0,
    dpi: int = 220,
    plot_col_key: bool = True,
    new_signature: bool = True,
    rand_seed: int | None = None,
    overwrite: bool = False,
    quiet: bool = False,
    continue_on_error: bool = True,
) -> PathwayResult | PathwayResultSet:
    """
    Overlay omics data on an SBGN pathway map.

    Parameters
    ----------
    pathway_id:
        A collection id (``"P00001"``, ``"SMP00001"``, ``"GLYCOLYSIS"``,
        ``"R-HSA-109688"``), a path to a local ``.sbgn`` file, or a sequence
        of either.
    gene_data, cpd_data:
        First column identifiers, remaining numeric columns one per condition.
    gene_idtype, cpd_idtype:
        Identifier systems of the input; translated to glyph identifiers via
        the bundled crosswalks.
    show_compartments:
        Shade and label compartments behind the map.

    Returns
    -------
    PathwayResult, or PathwayResultSet when *pathway_id* is a sequence.

    Examples
    --------
    >>> from pathview import sbgnview                      # doctest: +SKIP
    >>> res = sbgnview("P00001", gene_data=df, gene_idtype="SYMBOL")
    """
    from .vector_rendering import keggview_vector, render_vector_array

    if not isinstance(pathway_id, (str, Path)) and isinstance(pathway_id, (list, tuple, set)):
        results, failures = {}, {}
        for pid in pathway_id:
            try:
                res = sbgnview(
                    pid, gene_data=gene_data, cpd_data=cpd_data,
                    gene_idtype=gene_idtype, cpd_idtype=cpd_idtype,
                    sbgn_dir=sbgn_dir, out_dir=out_dir, out_suffix=out_suffix,
                    output_format=output_format, theme=theme,
                    show_compartments=show_compartments,
                    show_processes=show_processes, draw_edges=draw_edges,
                    gene_color=gene_color, cpd_color=cpd_color, limit=limit,
                    bins=bins, node_sum=node_sum, map_cpd_name=map_cpd_name,
                    figure_width=figure_width, dpi=dpi,
                    plot_col_key=plot_col_key, new_signature=new_signature,
                    rand_seed=rand_seed, overwrite=overwrite, quiet=quiet,
                )
                results[str(pid)] = res
            except PathviewError as exc:
                if not continue_on_error:
                    raise
                failures[str(pid)] = f"{type(exc).__name__}: {exc}"
                if not quiet:
                    print(f"[sbgnview] {pid}: {exc}")
        return PathwayResultSet(results, failures)

    log = (lambda *a: None) if quiet else (lambda *a: print("[sbgnview]", *a))
    out_dir = Path(out_dir) if out_dir is not None else Path(sbgn_dir)

    path = _resolve_source(pathway_id, sbgn_dir, overwrite)
    pathway = parse_sbgn(path)
    if not pathway.glyphs:
        raise PathwayNotFoundError(f"{path.name} contains no glyphs.")

    node_data = sbgn_to_df(pathway, include_processes=show_processes)
    edge_data = sbgn_edges(pathway)
    comps = sbgn_compartments(pathway) if show_compartments else None

    def _pick(v, key):
        return v.get(key) if isinstance(v, dict) else v

    from .pathview import _as_scale
    gsc = _as_scale(gene_color, DEFAULT_GENE_PALETTE, "Gene / protein log2FC",
                    limit=_pick(limit, "gene") if limit is not None else 1.0,
                    bins=bins)
    csc = _as_scale(cpd_color, DEFAULT_CPD_PALETTE, "Metabolite log2FC",
                    limit=_pick(limit, "cpd") if limit is not None else 1.0,
                    bins=bins)

    gres = sbgn_node_map(gene_data, node_data, gene_idtype, _GENE_TYPES,
                         node_sum, rand_seed, detailed=True)
    cres = sbgn_node_map(cpd_data, node_data, cpd_idtype, _CPD_TYPES,
                         node_sum, rand_seed, detailed=True)

    plot_gene = gres.data if (gres.ok and gene_data is not None) else None
    plot_cpd = cres.data if (cres.ok and cpd_data is not None) else None

    diagnostics = {"glyphs": len(pathway.glyphs), "arcs": len(pathway.arcs),
                   "compartments": len(pathway.compartments),
                   "language": pathway.language, "source_file": str(path)}
    if gene_data is not None:
        diagnostics["gene"] = gres.summary()
        log(f"genes — {gres.summary()}")
    if cpd_data is not None:
        diagnostics["cpd"] = cres.summary()
        log(f"metabolites — {cres.summary()}")
    if gene_data is not None and cpd_data is not None \
            and gres.n_nodes_with_data == 0 and cres.n_nodes_with_data == 0:
        raise MappingError(
            f"Nothing mapped onto {path.stem}. SBGN glyph ids are database "
            f"specific; check gene_idtype ({gene_idtype}) and cpd_idtype "
            f"({cpd_idtype}), and see crosswalk_routes() for the conversions "
            "the bundled table supports."
        )

    from .pathview import _colors
    cols_gene, cols_cpd = _colors(plot_gene, gsc), _colors(plot_cpd, csc)

    color_map: dict[str, list[str]] = {}
    for cdf in (cols_gene, cols_cpd):
        if cdf is None:
            continue
        ccols = [c for c in cdf.columns if c.endswith("_col")]
        for row in cdf.iter_rows(named=True):
            color_map[str(row["id"])] = [row[c] for c in ccols]

    label_map: dict[str, str] = {}
    if map_cpd_name:
        for row in node_data.iter_rows(named=True):
            if row["type"] != "compound" or row.get("label"):
                continue
            for nm in (row.get("kegg_names") or []):
                resolved = compound_name(str(nm))
                if resolved != str(nm):
                    label_map[str(row["entry_id"])] = resolved
                    break

    name = Path(str(pathway_id)).stem
    plot_title = title if title is not None else (pathway.pathway_name or name)

    out_path = keggview_vector(
        node_data=node_data, edge_data=edge_data, color_map=color_map,
        label_map=label_map, pathway_name=name, title=plot_title,
        subtitle=subtitle, out_dir=out_dir, out_suffix=out_suffix,
        output_format=output_format,
        gene_scale=gsc if plot_gene is not None else None,
        cpd_scale=csc if plot_cpd is not None else None,
        theme=theme, dpi=dpi, figure_width=figure_width,
        new_signature=new_signature, plot_col_key=plot_col_key,
        draw_edges=draw_edges, compartments=comps,
    )
    log(f"wrote {out_path}")

    result = PathwayResult(
        pathway_id=name, pathway_name=plot_title, species="",
        plot_data_gene=plot_gene, plot_data_cpd=plot_cpd,
        cols_gene=cols_gene, cols_cpd=cols_cpd,
        node_data=node_data, edge_data=edge_data, output_path=out_path,
        gene_scale=gsc, cpd_scale=csc, diagnostics=diagnostics,
    )
    try:
        result.frame = render_vector_array(
            node_data, edge_data, color_map, label_map, theme=theme,
            draw_edges=draw_edges, compartments=comps,
        )
    except Exception:
        result.frame = None
    return result


def sbgnview_batch(pathway_ids: Sequence[str], **kw) -> PathwayResultSet:
    """Render several SBGN maps; equivalent to passing a list to :func:`sbgnview`."""
    return sbgnview(list(pathway_ids), **kw)
