"""
graph_rendering.py
NetworkX graph view.

Fixes over v2.x
---------------
v2.x's ``keggview_graph`` called ``G.add_node`` in a loop and **never added a
single edge**, then passed ``arrows=True`` to ``nx.draw_networkx``.  The
result was a scatter of disconnected boxes labelled "graph view".  It also
had no way to reach the relation data, because only the node DataFrame was
passed in.  This module takes the edge table and builds a real graph, with
subtype-coloured edges, optional layout algorithms, and graph metrics.

Public API
----------
  build_graph      : node + edge DataFrames -> nx.DiGraph
  keggview_graph   : render the graph and save
  pathway_metrics  : degree / centrality summary
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .color_mapping import ColorScale, draw_dual_key
from .errors import RenderError
from .layout import node_boxes
from .utils import is_transparent, to_hex
from .vector_rendering import _DEFAULT_EDGE, EDGE_STYLE, THEMES

#: Layouts that delegate to SciPy inside NetworkX.  SciPy is deliberately not
#: a hard dependency — it is a large install for one optional layout, and
#: ``layout="kgml"`` (the default) uses KEGG's own coordinates instead.
_SCIPY_LAYOUTS = frozenset({"kamada_kawai"})


def _has_scipy() -> bool:
    """True when SciPy is importable."""
    import importlib.util
    return importlib.util.find_spec("scipy") is not None


def available_layouts() -> list[str]:
    """
    Graph layouts usable in this installation.

    ``kamada_kawai`` appears only when SciPy is present, so callers can offer
    the choice honestly instead of failing at render time.
    """
    base = ["kgml", "spring", "circular", "shell"]
    return base + (["kamada_kawai"] if _has_scipy() else [])


def build_graph(node_data: pl.DataFrame, edge_data: pl.DataFrame | None = None):
    """Build a ``networkx.DiGraph`` carrying node attributes and edge subtypes."""
    try:
        import networkx as nx
    except ImportError as exc:                                # pragma: no cover
        raise RenderError("networkx is required for the graph view: "
                          "pip install networkx") from exc

    G = nx.DiGraph()
    for row in node_data.iter_rows(named=True):
        attrs = {k: v for k, v in row.items() if k != "kegg_names"}
        G.add_node(str(row["entry_id"]), **attrs)

    if edge_data is not None and not edge_data.is_empty():
        for row in edge_data.iter_rows(named=True):
            s, t = str(row["source"]), str(row["target"])
            if s in G and t in G:
                G.add_edge(s, t,
                           subtype=row.get("subtype", ""),
                           edge_type=row.get("edge_type", ""),
                           kind=row.get("source_kind", ""))
    return G


def pathway_metrics(G) -> dict:
    """Degree and centrality summary for a pathway graph."""
    import networkx as nx

    if G.number_of_nodes() == 0:
        return {"nodes": 0, "edges": 0}
    deg = dict(G.degree())
    top = sorted(deg.items(), key=lambda kv: kv[1], reverse=True)[:5]
    und = G.to_undirected()
    return {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "density": round(nx.density(G), 5),
        "components": nx.number_connected_components(und),
        "mean_degree": round(sum(deg.values()) / len(deg), 3),
        "hubs": [(G.nodes[n].get("label") or n, d) for n, d in top],
    }


def keggview_graph(
    node_data: pl.DataFrame,
    edge_data: pl.DataFrame | None = None,
    color_map: dict[str, list[str]] | None = None,
    pathway_name: str = "pathway",
    title: str | None = None,
    out_dir: str | Path = ".",
    out_suffix: str = "pathview",
    output_format: str = "pdf",
    layout: str = "kgml",
    gene_scale: ColorScale | None = None,
    cpd_scale: ColorScale | None = None,
    theme: str = "publication",
    node_size: float = 420.0,
    font_size: float = 5.5,
    dpi: int = 200,
    plot_col_key: bool = True,
    new_signature: bool = True,
) -> Path:
    """
    Render the pathway as a graph diagram.

    ``layout`` may be ``kgml`` (KEGG's own coordinates), ``spring``,
    ``kamada_kawai``, ``circular`` or ``shell``.
    """
    import matplotlib
    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt
    import networkx as nx

    th = THEMES.get(theme, THEMES["publication"])
    color_map = color_map or {}
    G = build_graph(node_data, edge_data)
    if G.number_of_nodes() == 0:
        raise RenderError(f"{pathway_name}: no nodes to draw.")

    if layout == "kgml":
        pos = {b.entry_id: (b.x, -b.y) for b in node_boxes(node_data)}
        missing = [n for n in G.nodes if n not in pos]
        if missing:
            G.remove_nodes_from(missing)
    else:
        fn = {"spring": nx.spring_layout, "kamada_kawai": nx.kamada_kawai_layout,
              "circular": nx.circular_layout, "shell": nx.shell_layout}.get(layout)
        if fn is None:
            raise ValueError(
                f"Unknown layout {layout!r}. Choose kgml, spring, "
                "kamada_kawai, circular or shell."
            )
        if layout in _SCIPY_LAYOUTS and not _has_scipy():
            raise RenderError(
                f"layout={layout!r} needs SciPy, which is an optional "
                "dependency: pip install 'pathview-plus[layouts]' (or "
                "pip install scipy). layout='kgml' uses KEGG's own "
                "coordinates and needs nothing extra."
            )
        pos = fn(G, seed=42) if layout == "spring" else fn(G)

    fills, labels = [], {}
    for n in G.nodes:
        cols = color_map.get(n)
        c = cols[0] if cols else None
        fills.append(th["unmapped"] if is_transparent(c) else to_hex(c, False))
        labels[n] = (G.nodes[n].get("label") or "")[:14]

    keys = [s for s in (gene_scale, cpd_scale) if s is not None] if plot_col_key else []
    fig, ax = plt.subplots(figsize=(15, 11), facecolor=th["bg"])
    ax.set_facecolor(th["map_bg"])

    by_style: dict[tuple, list] = {}
    for u, v, d in G.edges(data=True):
        st = EDGE_STYLE.get((d.get("subtype") or "").lower(),
                            EDGE_STYLE.get((d.get("edge_type") or "").lower(),
                                           _DEFAULT_EDGE))
        by_style.setdefault((st["color"], st["style"]), []).append((u, v))

    mpl_style = {"-": "solid", "--": "dashed", ":": "dotted"}
    for (color, style), edges in by_style.items():
        nx.draw_networkx_edges(
            G, pos, edgelist=edges, ax=ax, edge_color=color,
            style=mpl_style.get(style, "solid"), width=0.8, alpha=0.6,
            arrows=True, arrowsize=7, arrowstyle="-|>",
            node_size=node_size, connectionstyle="arc3,rad=0.08",
        )

    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=fills, node_size=node_size,
                           edgecolors=th["border"], linewidths=0.5, node_shape="o")
    nx.draw_networkx_labels(G, pos, labels=labels, ax=ax, font_size=font_size,
                            font_color=th["text"])

    ax.set_title(title or pathway_name, fontsize=13, fontweight="semibold",
                 color=th["title"], pad=12)
    ax.set_axis_off()

    if keys:
        draw_dual_key(fig, gene_scale, cpd_scale,
                      rect=(0.15, 0.035, 0.70, 0.022), label_size=7.5)
    if new_signature:
        fig.text(0.99, 0.01, "pathview-plus", ha="right", va="bottom",
                 fontsize=6.4, color=th["muted"], style="italic")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{pathway_name}.{out_suffix}.{output_format}"
    fig.savefig(out_path, dpi=dpi, facecolor=th["bg"], bbox_inches="tight",
                pad_inches=0.15)
    plt.close(fig)
    return out_path
