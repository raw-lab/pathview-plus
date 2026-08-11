"""
pathview.py
The orchestrator: resolve species, obtain the pathway, map data, render.

Fixes over v2.x
---------------
* Returned ``{}`` on every failure path, so callers could not distinguish
  "network down" from "pathway has too few nodes" — and the test suite's
  ``assert result is not None`` passed for all of them.  Now returns a
  :class:`PathwayResult` (falsy when empty) and raises typed errors.
* ``numeric_id = pathway_name.replace(species_code, "")`` used ``str.replace``,
  which strips *every* occurrence; ``removeprefix`` is correct.
* Colour configuration was eight parallel ``{"gene": ..., "cpd": ...}`` dicts
  threaded through five call sites.  Replaced by two
  :class:`~pathview.color_mapping.ColorScale` objects — which is also what
  makes independent RNA-seq and metabolite scales possible.
* ``keggview_svg`` was called with a filtered kwargs dict whose keys did not
  match its signature.
* Compound labels showed raw accessions; they now resolve to names.

Public API
----------
  pathview  : the main entry point
"""

from __future__ import annotations

import warnings
from collections.abc import Callable, Sequence
from pathlib import Path

import polars as pl

from .color_mapping import ColorScale, node_color, value_columns
from .constants import DEFAULT_CPD_PALETTE, DEFAULT_GENE_PALETTE, NA_COLOR, SumMethod
from .databases import download_kegg
from .errors import MappingError, PathviewError, PathwayNotFoundError
from .expansion import expand_nodes, split_groups
from .highlighting import PathwayResult, PathwayResultSet
from .id_mapping import cpd_id_map, eg2id, id2eg
from .kgml_parser import node_info, parse_kgml, pathway_edges
from .mol_data import compound_name, mol_sum
from .node_mapping import node_map
from .organisms import get_species_code


def _as_scale(value, fallback: str, label: str, **kw) -> ColorScale:
    """Accept a ColorScale, a palette name, or a kwargs dict."""
    if isinstance(value, ColorScale):
        return value
    if isinstance(value, str):
        return ColorScale(palette=value, label=label, _fallback=fallback, **kw)
    if isinstance(value, dict):
        merged = {**kw, **value}
        merged.setdefault("label", label)
        merged.setdefault("_fallback", fallback)
        return ColorScale(**merged)
    return ColorScale(label=label, _fallback=fallback, **kw)


def pathview(
    pathway_id: str | Sequence[str],
    gene_data: pl.DataFrame | None = None,
    cpd_data: pl.DataFrame | None = None,
    species: str = "hsa",
    *,
    # --- input identifiers -------------------------------------------------
    gene_idtype: str = "ENTREZ",
    cpd_idtype: str = "KEGG",
    # --- paths -------------------------------------------------------------
    kegg_dir: str | Path = ".",
    out_dir: str | Path | None = None,
    out_suffix: str = "pathview",
    # --- rendering ---------------------------------------------------------
    render_mode: str = "auto",          # auto | native | vector | graph | svg
    output_format: str = "png",         # png | pdf | svg
    theme: str = "publication",
    title: str | None = None,
    subtitle: str | None = None,
    figure_width: float = 14.0,
    dpi: int = 220,
    draw_edges: bool = True,
    show_link_edges: bool = False,
    plot_col_key: bool = True,
    new_signature: bool = True,
    # --- colour scales -----------------------------------------------------
    gene_color: ColorScale | str | dict | None = None,
    cpd_color: ColorScale | str | dict | None = None,
    limit: float | dict | None = None,
    bins: int | dict = 10,
    both_dirs: bool | dict = True,
    discrete: bool | dict = False,
    na_col: str = NA_COLOR,
    trans_fun: Callable | None = None,
    # --- mapping -----------------------------------------------------------
    node_sum: SumMethod = "sum",
    map_symbol: bool = True,
    map_cpd_name: bool = True,
    map_null: bool = True,
    min_nnodes: int = 3,
    split_group: bool = False,
    expand_node: bool = False,
    rand_seed: int | None = None,
    # --- misc --------------------------------------------------------------
    quiet: bool = False,
    continue_on_error: bool = True,
    **kwargs,
) -> PathwayResult | PathwayResultSet:
    """
    Overlay transcript and/or metabolite data on a pathway diagram.

    Parameters
    ----------
    pathway_id:
        KEGG number (``"04110"``) or full id (``"hsa04110"``).
    gene_data, cpd_data:
        DataFrames whose first column holds identifiers and whose remaining
        numeric columns hold one value per condition.  Multiple value columns
        split each node into vertical slices.
    species:
        Anything :func:`~pathview.organisms.get_species_code` accepts —
        ``"hsa"``, ``"human"``, ``"Homo sapiens"``, ``"9606"``.
    render_mode:
        ``native`` overlays KEGG's PNG; ``vector`` draws the map from the KGML
        (works offline, scales losslessly); ``graph`` draws a NetworkX
        diagram; ``svg`` writes a standalone SVG.  ``auto`` picks ``native``
        when a background image is present and ``vector`` otherwise.
    gene_color, cpd_color:
        Independent colour scales.  Pass a palette name
        (``gene_color="rdbu"``), a dict of ColorScale kwargs, or a
        :class:`ColorScale`.  These are separate so a transcript log2FC and a
        metabolite log2FC are never read off the same key.

    Returns
    -------
    PathwayResult
        Composable with ``+ highlight_nodes(...)``, indexable like the v2.x
        dict, and falsy when nothing could be mapped.

    Examples
    --------
    >>> from pathview import pathview, demo_gene_data
    >>> res = pathview("00020", gene_data=demo_gene_data(2), species="human",
    ...                render_mode="vector", output_format="pdf")
    >>> res.output_path.name
    'hsa00020.pathview.pdf'
    """
    # Batch: R pathview accepts a vector of pathway ids, and looping in user
    # code loses the shared setup (species resolution, colour scales, one
    # download session).  A sequence renders them all and returns a
    # PathwayResultSet; a single id still returns a single PathwayResult, so
    # existing calls are unaffected.
    if not isinstance(pathway_id, str) and isinstance(pathway_id, (list, tuple, set)):
        ids = [str(p) for p in pathway_id]
        if not ids:
            raise ValueError("pathway_id is an empty sequence.")
        results: dict[str, PathwayResult] = {}
        failures: dict[str, str] = {}
        for pid in ids:
            try:
                res = pathview(
                    pid, gene_data=gene_data, cpd_data=cpd_data, species=species,
                    gene_idtype=gene_idtype, cpd_idtype=cpd_idtype,
                    kegg_dir=kegg_dir, out_dir=out_dir, out_suffix=out_suffix,
                    render_mode=render_mode, output_format=output_format,
                    theme=theme, title=title, subtitle=subtitle,
                    figure_width=figure_width, dpi=dpi, draw_edges=draw_edges,
                    show_link_edges=show_link_edges, plot_col_key=plot_col_key,
                    new_signature=new_signature, gene_color=gene_color,
                    cpd_color=cpd_color, limit=limit, bins=bins,
                    both_dirs=both_dirs, discrete=discrete, na_col=na_col,
                    trans_fun=trans_fun, node_sum=node_sum,
                    map_symbol=map_symbol, map_cpd_name=map_cpd_name,
                    map_null=map_null, min_nnodes=min_nnodes,
                    rand_seed=rand_seed, quiet=quiet,
                    split_group=split_group, expand_node=expand_node,
                    **kwargs,
                )
                if res:
                    results[pid] = res
                else:
                    failures[pid] = "no mappable nodes"
            except PathviewError as exc:
                if not continue_on_error:
                    raise
                failures[pid] = f"{type(exc).__name__}: {exc}"
                if not quiet:
                    print(f"[pathview] {pid}: {exc}")
        return PathwayResultSet(results, failures)

    if gene_data is None and cpd_data is None and not map_null:
        raise ValueError("Provide gene_data, cpd_data, or set map_null=True.")

    kegg_dir = Path(kegg_dir)
    out_dir = Path(out_dir) if out_dir is not None else kegg_dir
    log = (lambda *a: None) if quiet else (lambda *a: print("[pathview]", *a))

    # ---- species -----------------------------------------------------------
    info = get_species_code(species)
    code = info.kegg_code
    if code in ("ko", "ec", "rn"):
        gene_idtype = "KEGG"

    # ---- pathway id --------------------------------------------------------
    pid = str(pathway_id).strip()
    pathway_name = pid if pid.startswith(code) else f"{code}{pid}"
    numeric_id = pathway_name.removeprefix(code)      # not str.replace

    # ---- colour scales -----------------------------------------------------
    def _pick(v, key):
        return v.get(key) if isinstance(v, dict) else v

    gsc = _as_scale(gene_color, DEFAULT_GENE_PALETTE, "RNA-seq log2FC",
                    limit=_pick(limit, "gene") if limit is not None else 1.0,
                    bins=_pick(bins, "gene") or 10,
                    both_dirs=bool(_pick(both_dirs, "gene")),
                    discrete=bool(_pick(discrete, "gene")),
                    na_col=na_col, trans_fun=trans_fun)
    csc = _as_scale(cpd_color, DEFAULT_CPD_PALETTE, "Metabolite log2FC",
                    limit=_pick(limit, "cpd") if limit is not None else 1.0,
                    bins=_pick(bins, "cpd") or 10,
                    both_dirs=bool(_pick(both_dirs, "cpd")),
                    discrete=bool(_pick(discrete, "cpd")),
                    na_col=na_col, trans_fun=trans_fun)

    # ---- identifier conversion --------------------------------------------
    if gene_data is not None and gene_idtype.upper() not in ("ENTREZ", "ENTREZID", "KEGG"):
        gene_data = _convert(gene_data, id2eg, gene_idtype, code, node_sum, log)
    if cpd_data is not None and cpd_idtype.upper() != "KEGG":
        cpd_data = _convert_cpd(cpd_data, cpd_idtype, node_sum, log)

    # ---- obtain the pathway ------------------------------------------------
    xml_path = kegg_dir / f"{pathway_name}.xml"
    png_path = kegg_dir / f"{pathway_name}.png"
    want_png = render_mode in ("native", "auto")

    needed = [t for t, p in (("xml", xml_path), ("png", png_path))
              if not p.exists() and (t == "xml" or want_png)]
    if needed:
        status = download_kegg(numeric_id, species=code, kegg_dir=kegg_dir,
                               file_type=needed)
        log(f"{pathway_name}: download {status.get(pathway_name)}")

    if not xml_path.exists():
        raise PathwayNotFoundError(
            f"No KGML for {pathway_name} at {xml_path}. Place the file there "
            "manually if you are offline: pathview-plus renders it without "
            "any further network access in render_mode='vector'."
        )

    pathway = parse_kgml(xml_path)
    node_data = node_info(pathway)
    edge_data = pathway_edges(pathway)

    # R pathview's split.group / expand.node, with the expansion recorded.
    expansion_notes: list[str] = []
    if split_group:
        res = split_groups(node_data, detailed=True)
        node_data, edge_data = res.data, _remap_edges(edge_data, res.data)
        expansion_notes.append(res.summary())
    if expand_node:
        res = expand_nodes(node_data, detailed=True)
        node_data, edge_data = res.data, _remap_edges(edge_data, res.data)
        expansion_notes.append(res.summary())

    mappable = node_data.filter(
        pl.col("x").is_not_null() & pl.col("y").is_not_null()
    )
    if mappable.height < min_nnodes:
        warnings.warn(
            f"{pathway_name}: only {mappable.height} positioned nodes "
            f"(min_nnodes={min_nnodes}); skipping.", stacklevel=2)
        return PathwayResult(pathway_id=numeric_id, pathway_name=pathway_name,
                             species=code)
    node_data = mappable

    # ---- map data onto nodes ----------------------------------------------
    diagnostics: dict = {"species": info.display_name, "nodes": node_data.height,
                         "edges": edge_data.height}
    if expansion_notes:
        diagnostics["expansion"] = "; ".join(expansion_notes)

    gene_types = "ortholog" if code in ("ko", "ec") else "gene"
    gres = node_map(gene_data, node_data, gene_types, node_sum,
                    rand_seed=rand_seed, detailed=True)
    cres = node_map(cpd_data, node_data, "compound", node_sum,
                    rand_seed=rand_seed, detailed=True)

    # Only report a data frame for a class the caller actually supplied.
    # node_map() returns the bare layout when data is None (R's map.null
    # behaviour) and that layout is still used for rendering, but reporting
    # it as ``plot_data_gene`` would imply gene data had been mapped.
    plot_gene = gres.data if (gres.ok and gene_data is not None) else None
    plot_cpd = cres.data if (cres.ok and cpd_data is not None) else None
    if gene_data is not None:
        diagnostics["gene"] = gres.summary()
        log(f"genes  — {gres.summary()}")
    if cpd_data is not None:
        diagnostics["cpd"] = cres.summary()
        log(f"metabolites — {cres.summary()}")

    if gene_data is not None and cpd_data is not None \
            and gres.n_nodes_with_data == 0 and cres.n_nodes_with_data == 0:
        raise MappingError(
            f"Nothing mapped onto {pathway_name}. Check gene_idtype "
            f"({gene_idtype}), cpd_idtype ({cpd_idtype}) and species ({code})."
        )

    # ---- colours -----------------------------------------------------------
    cols_gene = _colors(plot_gene, gsc)
    cols_cpd = _colors(plot_cpd, csc)

    color_map: dict[str, list[str]] = {}
    for cdf in (cols_gene, cols_cpd):
        if cdf is None:
            continue
        ccols = [c for c in cdf.columns if c.endswith("_col")]
        for row in cdf.iter_rows(named=True):
            color_map[str(row["id"])] = [row[c] for c in ccols]

    # ---- labels ------------------------------------------------------------
    label_map: dict[str, str] = {}
    if map_cpd_name:
        for row in node_data.iter_rows(named=True):
            if row["type"] == "compound":
                names = row.get("kegg_names") or []
                if names:
                    label_map[str(row["entry_id"])] = compound_name(names[0])
    if map_symbol and gene_data is not None and code not in ("ko", "ec"):
        label_map.update(_symbol_labels(node_data, code, log))
    label_map.update(kwargs.pop("label_map", {}) or {})

    # ---- render ------------------------------------------------------------
    mode = render_mode
    if mode == "auto":
        mode = "native" if png_path.exists() else "vector"
    if mode == "native" and not png_path.exists():
        log("no KEGG background image available — falling back to vector mode")
        mode = "vector"

    plot_title = title if title is not None else (
        f"{pathway.title or pathway_name} — {info.display_name}")

    keys_gene = gsc if (plot_gene is not None and gene_data is not None) else None
    keys_cpd = csc if (plot_cpd is not None and cpd_data is not None) else None

    out_path = _render(
        mode=mode, output_format=output_format,
        node_data=node_data, edge_data=edge_data,
        color_map=color_map, label_map=label_map,
        plot_gene=plot_gene, cols_gene=cols_gene,
        plot_cpd=plot_cpd, cols_cpd=cols_cpd,
        pathway_name=pathway_name, kegg_dir=kegg_dir, out_dir=out_dir,
        out_suffix=out_suffix, title=plot_title, subtitle=subtitle,
        gene_scale=keys_gene, cpd_scale=keys_cpd, theme=theme,
        dpi=dpi, figure_width=figure_width, draw_edges=draw_edges,
        show_link_edges=show_link_edges, plot_col_key=plot_col_key,
        new_signature=new_signature,
    )
    log(f"wrote {out_path}")

    result = PathwayResult(
        pathway_id=numeric_id, pathway_name=pathway_name, species=code,
        plot_data_gene=plot_gene, plot_data_cpd=plot_cpd,
        cols_gene=cols_gene, cols_cpd=cols_cpd,
        node_data=node_data, edge_data=edge_data,
        output_path=out_path, gene_scale=gsc, cpd_scale=csc,
        diagnostics=diagnostics,
    )
    # The raster carried on the result is the *map* in KGML coordinates, not
    # the composed figure: highlighting draws in KGML space, and a figure has
    # a title band, colour keys, padding and a dpi scale that would silently
    # displace every highlight.
    result.frame = _map_raster(
        mode=mode, node_data=node_data, edge_data=edge_data,
        color_map=color_map, label_map=label_map, theme=theme,
        plot_gene=plot_gene, cols_gene=cols_gene,
        plot_cpd=plot_cpd, cols_cpd=cols_cpd,
        png_path=png_path, draw_edges=draw_edges,
        show_link_edges=show_link_edges,
    )
    return result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _remap_edges(edge_data: pl.DataFrame | None,
                 node_data: pl.DataFrame) -> pl.DataFrame | None:
    """
    Point edges at the nodes that replaced the ones they referenced.

    After splitting or expansion the original entry ids no longer exist, so an
    unremapped edge table silently drops every edge — the map would render
    with the right nodes and no connectivity at all.
    """
    if edge_data is None or edge_data.is_empty():
        return edge_data

    replacement: dict[str, list[str]] = {}
    for col in ("parent_entry", "parent_group"):
        if col not in node_data.columns:
            continue
        for row in node_data.iter_rows(named=True):
            parent = row.get(col)
            if parent and parent != row["entry_id"]:
                replacement.setdefault(str(parent), []).append(str(row["entry_id"]))

    present = set(node_data["entry_id"].to_list())

    def resolve(entry: str) -> list[str]:
        if entry in present:
            return [entry]
        return replacement.get(entry, [])

    rows = []
    for row in edge_data.iter_rows(named=True):
        for s_id in resolve(str(row["source"])):
            for t_id in resolve(str(row["target"])):
                if s_id != t_id:
                    new = dict(row)
                    new["source"], new["target"] = s_id, t_id
                    rows.append(new)
    if not rows:
        return edge_data.head(0)
    return pl.DataFrame(rows, schema=edge_data.schema).unique(
        subset=["source", "target", "subtype"], keep="first")


def _convert(data: pl.DataFrame, fn, idtype: str, code: str,
             node_sum: SumMethod, log) -> pl.DataFrame:
    res = fn(data[data.columns[0]].cast(pl.String).to_list(),
             category=idtype, org=code, detailed=True)
    log(f"gene ids — {res.summary()}")
    if res.n_resolved == 0:
        raise MappingError(
            f"No {idtype} identifier could be converted to Entrez. "
            "If the input is already Entrez, pass gene_idtype='ENTREZ'."
        )
    return mol_sum(data, res.data, sum_method=node_sum)


def _convert_cpd(data: pl.DataFrame, idtype: str, node_sum: SumMethod,
                 log) -> pl.DataFrame:
    res = cpd_id_map(data[data.columns[0]].cast(pl.String).to_list(),
                     in_type=idtype, out_type="KEGG", detailed=True)
    log(f"compound ids — {res.summary()}")
    if res.n_resolved == 0:
        raise MappingError(
            f"No {idtype} compound identifier could be converted to KEGG. "
            "Supported types: see supported_cpd_idtypes()."
        )
    return mol_sum(data, res.data, sum_method=node_sum)


def _colors(plot_data: pl.DataFrame | None, scale: ColorScale):
    if plot_data is None:
        return None
    vcols = value_columns(plot_data)
    if not vcols:
        return None
    return node_color(plot_data.rename({"entry_id": "id"}), scale,
                      id_col="id", value_cols=vcols)


def _symbol_labels(node_data: pl.DataFrame, code: str, log) -> dict[str, str]:
    """Replace Entrez-based gene labels with symbols, if a service is reachable."""
    rows = node_data.filter(pl.col("type").is_in(["gene", "ortholog"]))
    if rows.is_empty():
        return {}
    first: dict[str, str] = {}
    for row in rows.iter_rows(named=True):
        names = row.get("kegg_names") or []
        if names and not (row.get("label") or "").strip():
            first[str(row["entry_id"])] = str(names[0])
    if not first:
        return {}
    try:
        res = eg2id(list(first.values()), category="SYMBOL", org=code, detailed=True)
    except Exception as exc:
        log(f"symbol lookup skipped: {exc}")
        return {}
    table = dict(zip(res.data["ENTREZID"].to_list(), res.data["SYMBOL"].to_list()))
    return {eid: table[gid] for eid, gid in first.items() if table.get(gid)}


def _render(*, mode: str, output_format: str, **kw) -> Path:
    """Dispatch to the requested renderer with only the arguments it accepts."""
    if mode == "native":
        from .rendering import keggview_native
        return keggview_native(
            plot_data_gene=kw["plot_gene"], cols_gene=kw["cols_gene"],
            plot_data_cpd=kw["plot_cpd"], cols_cpd=kw["cols_cpd"],
            node_data=kw["node_data"], pathway_name=kw["pathway_name"],
            kegg_dir=kw["kegg_dir"], out_dir=kw["out_dir"],
            out_suffix=kw["out_suffix"], gene_scale=kw["gene_scale"],
            cpd_scale=kw["cpd_scale"], title=kw["title"],
            new_signature=kw["new_signature"], plot_col_key=kw["plot_col_key"],
            dpi=kw["dpi"], output_format=output_format,
        )
    if mode == "svg":
        from .svg_rendering import keggview_svg
        return keggview_svg(
            node_data=kw["node_data"], edge_data=kw["edge_data"],
            color_map=kw["color_map"], label_map=kw["label_map"],
            pathway_name=kw["pathway_name"], title=kw["title"],
            out_dir=kw["out_dir"], out_suffix=kw["out_suffix"],
            gene_scale=kw["gene_scale"], cpd_scale=kw["cpd_scale"],
            theme=kw["theme"], new_signature=kw["new_signature"],
            plot_col_key=kw["plot_col_key"], draw_edges=kw["draw_edges"],
            show_link_edges=kw["show_link_edges"],
        )
    if mode == "graph":
        from .graph_rendering import keggview_graph
        return keggview_graph(
            node_data=kw["node_data"], edge_data=kw["edge_data"],
            color_map=kw["color_map"], pathway_name=kw["pathway_name"],
            title=kw["title"], out_dir=kw["out_dir"],
            out_suffix=kw["out_suffix"],
            output_format="pdf" if output_format == "svg" else output_format,
            gene_scale=kw["gene_scale"], cpd_scale=kw["cpd_scale"],
            theme=kw["theme"], dpi=kw["dpi"],
            plot_col_key=kw["plot_col_key"], new_signature=kw["new_signature"],
        )
    if mode == "vector":
        from .vector_rendering import keggview_vector
        return keggview_vector(
            node_data=kw["node_data"], edge_data=kw["edge_data"],
            color_map=kw["color_map"], label_map=kw["label_map"],
            pathway_name=kw["pathway_name"], title=kw["title"],
            subtitle=kw["subtitle"], out_dir=kw["out_dir"],
            out_suffix=kw["out_suffix"], output_format=output_format,
            gene_scale=kw["gene_scale"], cpd_scale=kw["cpd_scale"],
            theme=kw["theme"], dpi=kw["dpi"], figure_width=kw["figure_width"],
            new_signature=kw["new_signature"], plot_col_key=kw["plot_col_key"],
            draw_edges=kw["draw_edges"], show_link_edges=kw["show_link_edges"],
        )
    raise ValueError(
        f"Unknown render_mode {mode!r}. Choose auto, native, vector, graph or svg."
    )


def _map_raster(*, mode: str, png_path: Path, **kw):
    """
    Produce the map raster in KGML coordinates for post-hoc modification.

    ``native`` reuses KEGG's own image, whose pixels *are* KGML units.
    Everything else rasterises the pathway axes alone at a known scale.  SVG
    output has no raster of its own, so it borrows the vector one, which keeps
    ``result + highlight_nodes(...)`` working in every mode.
    """
    try:
        if mode == "native" and png_path.exists():
            from .rendering import render_native_array
            return render_native_array(
                kw["plot_gene"], kw["cols_gene"],
                kw["plot_cpd"], kw["cols_cpd"], background=png_path,
            )
        from .vector_rendering import render_vector_array
        return render_vector_array(
            kw["node_data"], kw["edge_data"], kw["color_map"], kw["label_map"],
            theme=kw["theme"], draw_edges=kw["draw_edges"],
            show_link_edges=kw["show_link_edges"],
        )
    except Exception:
        return None
