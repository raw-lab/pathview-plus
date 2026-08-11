"""
vector_rendering.py
Draw the pathway from its KGML/SBGN geometry as true vector output.

Why this exists
---------------
v2.x had exactly one good-looking mode, and it required KEGG's raster PNG.
Behind a firewall — or for any organism whose map image KEGG does not serve —
the only fallback was ``keggview_graph``, a NetworkX spring-ish scatter with
no edges.  This renderer draws the map itself from the coordinates already
present in the KGML, so:

  * it needs no background image and works entirely offline;
  * output is real vector geometry (PDF/SVG scale losslessly for figures);
  * labels are laid out with knowledge of node extents, so they fit;
  * edges are drawn as curves with proper KEGG arrowheads and inhibition bars;
  * gene and metabolite scales are drawn as two separate keys.

Public API
----------
  keggview_vector : render a pathway to png/pdf/svg
  draw_pathway    : draw onto a caller-supplied Axes
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import polars as pl

from .color_mapping import ColorScale, draw_dual_key
from .layout import (
    Extent,
    NodeBox,
    figure_size,
    fit_label,
    label_font_size,
    node_boxes,
    slice_bounds,
)
from .splines import offset_endpoints, route_edge_spline
from .utils import contrast_text_color, is_transparent, to_hex


def _halo(theme: dict, width: float = 1.6):
    """White halo so labels stay legible where they overlap edges."""
    from matplotlib import patheffects
    return [patheffects.withStroke(linewidth=width, foreground=theme["map_bg"])]

# ---------------------------------------------------------------------------
# Edge styling (KEGG subtype conventions)
# ---------------------------------------------------------------------------

EDGE_STYLE: dict[str, dict] = {
    "activation":          {"color": "#2E7D32", "style": "-",  "head": "arrow"},
    "expression":          {"color": "#388E3C", "style": "--", "head": "arrow"},
    "inhibition":          {"color": "#C62828", "style": "-",  "head": "bar"},
    "repression":          {"color": "#AD1457", "style": "--", "head": "bar"},
    "indirect effect":     {"color": "#78909C", "style": ":",  "head": "arrow"},
    "indirect":            {"color": "#78909C", "style": ":",  "head": "arrow"},
    "state change":        {"color": "#8E24AA", "style": ":",  "head": "none"},
    "binding/association": {"color": "#1565C0", "style": "-",  "head": "none"},
    "dissociation":        {"color": "#00838F", "style": "-",  "head": "none"},
    "phosphorylation":     {"color": "#EF6C00", "style": "-",  "head": "arrow"},
    "dephosphorylation":   {"color": "#F9A825", "style": "-",  "head": "arrow"},
    "glycosylation":       {"color": "#6D4C41", "style": "-",  "head": "arrow"},
    "ubiquitination":      {"color": "#AB47BC", "style": "-",  "head": "arrow"},
    "methylation":         {"color": "#0097A7", "style": "-",  "head": "arrow"},
    "compound":            {"color": "#9E9E9E", "style": "-",  "head": "none"},
    "substrate":           {"color": "#546E7A", "style": "-",  "head": "arrow"},
    "product":             {"color": "#546E7A", "style": "-",  "head": "arrow"},
    # SBGN arc classes
    "production":          {"color": "#546E7A", "style": "-",  "head": "arrow"},
    "consumption":         {"color": "#546E7A", "style": "-",  "head": "none"},
    "catalysis":           {"color": "#2E7D32", "style": "-",  "head": "circle"},
    "stimulation":         {"color": "#2E7D32", "style": "-",  "head": "arrow"},
    "modulation":          {"color": "#7E57C2", "style": "-",  "head": "diamond"},
    "necessary stimulation": {"color": "#00695C", "style": "-", "head": "arrow"},
    "logic arc":           {"color": "#90A4AE", "style": "-",  "head": "none"},
}
_DEFAULT_EDGE = {"color": "#B0BEC5", "style": "-", "head": "arrow"}

THEMES: dict[str, dict] = {
    "publication": {
        "bg": "#FFFFFF", "map_bg": "#FAFAFA", "border": "#37474F",
        "text": "#212121", "muted": "#78909C", "grid": "#ECEFF1",
        "unmapped": "#F5F5F5", "title": "#111111",
        "compartment": "#E3F2FD", "compartment_edge": "#90CAF9",
    },
    "slate": {
        "bg": "#FFFFFF", "map_bg": "#F7F9FB", "border": "#455A64",
        "text": "#263238", "muted": "#90A4AE", "grid": "#E3E9ED",
        "unmapped": "#EEF2F5", "title": "#1B262C",
        "compartment": "#E8EDF1", "compartment_edge": "#B0BEC5",
    },
    "dark": {
        "bg": "#12161C", "map_bg": "#1A1F27", "border": "#4A5568",
        "text": "#E6EDF3", "muted": "#8B99A8", "grid": "#232A34",
        "unmapped": "#2A313B", "title": "#F0F6FC",
        "compartment": "#232B36", "compartment_edge": "#3C4757",
    },
}

#: Fill opacities for successively nested compartments, so an inner
#: compartment reads as inside an outer one rather than replacing it.
_COMPARTMENT_ALPHA = (0.55, 0.42, 0.32, 0.24, 0.18)


def _cover_compartments(extent: Extent, compartments: pl.DataFrame) -> Extent:
    """Widen an extent so shaded compartments are not clipped at the edge."""
    xs0 = (compartments["x"] - compartments["width"] / 2).min()
    xs1 = (compartments["x"] + compartments["width"] / 2).max()
    ys0 = (compartments["y"] - compartments["height"] / 2).min()
    ys1 = (compartments["y"] + compartments["height"] / 2).max()
    pad = 20.0
    return Extent(min(extent.x0, float(xs0) - pad), min(extent.y0, float(ys0) - pad),
                  max(extent.x1, float(xs1) + pad), max(extent.y1, float(ys1) + pad))


def _style_for(subtype: str, edge_type: str = "") -> dict:
    return EDGE_STYLE.get((subtype or "").lower(),
                          EDGE_STYLE.get((edge_type or "").lower(), _DEFAULT_EDGE))


# ---------------------------------------------------------------------------
# Node drawing
# ---------------------------------------------------------------------------

def _rounded(ax, box: NodeBox, x0: float, x1: float, face: str,
             edge: str, lw: float, radius: float, zorder: float):
    from matplotlib.patches import FancyBboxPatch, Rectangle
    w = x1 - x0
    if w <= 0:
        return
    # Only the outermost slices get rounded corners, so a sliced node still
    # reads as one shape rather than a row of separate pills.
    if radius > 0 and w > 2.5 * radius:
        patch = FancyBboxPatch(
            (x0 + radius, box.top + radius),
            max(0.1, w - 2 * radius), max(0.1, box.height - 2 * radius),
            boxstyle=f"round,pad={radius}", linewidth=lw,
            facecolor=face, edgecolor=edge, zorder=zorder,
        )
    else:
        patch = Rectangle((x0, box.top), w, box.height, linewidth=lw,
                          facecolor=face, edgecolor=edge, zorder=zorder)
    ax.add_patch(patch)


def _wedge(ax, box: NodeBox, x0: float, x1: float, face: str, zorder: float):
    """A vertical band of a circle, clipped to the circle."""
    from matplotlib.patches import Circle, Rectangle
    clip = Circle((box.x, box.y), box.radius, transform=ax.transData)
    band = Rectangle((x0, box.top), max(0.1, x1 - x0), box.height,
                     facecolor=face, edgecolor="none", zorder=zorder)
    ax.add_patch(band)
    band.set_clip_path(clip)


def draw_node(
    ax,
    box: NodeBox,
    colors: Sequence[str],
    theme: dict,
    show_label: bool = True,
    label_override: str | None = None,
    border_width: float = 0.6,
    label_scale: float = 1.0,
) -> None:
    """Draw one node, sliced into as many bands as it has colours."""
    from matplotlib.patches import Circle

    fills = [c for c in colors] or [theme["unmapped"]]
    fills = [theme["unmapped"] if is_transparent(c) else to_hex(c, False) for c in fills]
    bands = slice_bounds(box, len(fills))

    if box.is_round:
        ax.add_patch(Circle((box.x, box.y), box.radius,
                            facecolor=fills[0], edgecolor="none", zorder=2.0))
        for (x0, x1), face in zip(bands, fills):
            _wedge(ax, box, x0, x1, face, zorder=2.1)
        ax.add_patch(Circle((box.x, box.y), box.radius, facecolor="none",
                            edgecolor=theme["border"], linewidth=border_width,
                            zorder=2.6))
    else:
        r = box.corner_radius
        for (x0, x1), face in zip(bands, fills):
            rad = r if len(fills) == 1 else 0.0
            _rounded(ax, box, x0, x1, face, "none", 0.0, rad, 2.1)
        _rounded(ax, box, box.left, box.right, "none",
                 theme["border"], border_width, r, 2.6)

    if not show_label:
        return

    label = label_override if label_override is not None else box.label
    if not label:
        return

    text = fit_label(box, label)
    size = label_font_size(box) * label_scale
    fg = contrast_text_color(fills[len(fills) // 2], light="#FFFFFF", dark=theme["text"])

    if box.is_round:
        # Metabolite labels sit below the circle: a KEGG compound node is
        # ~8 px across and cannot contain readable text.
        ax.text(box.x, box.bottom + 2.2, text, ha="center", va="top",
                fontsize=max(3.9, min(5.4, size * 1.15)), color=theme["text"],
                zorder=3.0, linespacing=0.95,
                path_effects=_halo(theme))
    else:
        ax.text(box.x, box.y, text, ha="center", va="center",
                fontsize=size, color=fg, zorder=3.0, linespacing=0.92)


# ---------------------------------------------------------------------------
# Edge drawing
# ---------------------------------------------------------------------------

def draw_edge(
    ax,
    src: NodeBox,
    tgt: NodeBox,
    subtype: str = "",
    edge_type: str = "",
    theme: dict | None = None,
    curvature: float = 0.12,
    width: float = 0.9,
    alpha: float = 0.75,
) -> None:
    """Draw one edge as a curve with the KEGG head for its subtype."""
    from matplotlib.patches import FancyArrowPatch

    style = _style_for(subtype, edge_type)
    s, t = offset_endpoints(
        (src.x, src.y), (tgt.x, tgt.y),
        source_radius=src.radius if src.is_round else src.height / 2 + 1,
        target_radius=(tgt.radius + 2.5) if tgt.is_round else tgt.height / 2 + 3.0,
    )
    if np.hypot(t[0] - s[0], t[1] - s[1]) < 1.0:
        return

    curve = route_edge_spline(s, t, routing_mode="curved",
                              curvature=curvature, n_points=24)
    ax.plot(curve[:, 0], curve[:, 1], color=style["color"],
            linestyle=style["style"], linewidth=width, alpha=alpha,
            solid_capstyle="round", zorder=1.4)

    head = style["head"]
    if head == "none" or len(curve) < 2:
        return

    p_end, p_prev = curve[-1], curve[-3 if len(curve) > 3 else -2]
    arrowstyle = {
        "arrow": "-|>", "bar": "|-|", "circle": "-|>", "diamond": "-|>",
    }.get(head, "-|>")
    ax.add_patch(FancyArrowPatch(
        tuple(p_prev), tuple(p_end), arrowstyle=arrowstyle,
        mutation_scale=6.5 if head == "arrow" else 5.0,
        color=style["color"], linewidth=width, alpha=alpha,
        shrinkA=0, shrinkB=0, zorder=1.5,
    ))


# ---------------------------------------------------------------------------
# Whole-pathway drawing
# ---------------------------------------------------------------------------

def draw_compartments(ax, compartments: pl.DataFrame, theme: dict,
                      label_size: float = 7.5) -> None:
    """
    Shade compartments behind the pathway.

    Drawn largest-first with decreasing opacity so nesting reads correctly,
    and labelled at the top-left where a label is least likely to collide
    with the biology.
    """
    from matplotlib.patches import FancyBboxPatch

    if compartments is None or compartments.is_empty():
        return

    for depth, row in enumerate(compartments.iter_rows(named=True)):
        w, h = float(row["width"]), float(row["height"])
        if w <= 0 or h <= 0:
            continue
        alpha = _COMPARTMENT_ALPHA[min(depth, len(_COMPARTMENT_ALPHA) - 1)]
        pad = min(6.0, min(w, h) / 8.0)
        ax.add_patch(FancyBboxPatch(
            (row["x"] - w / 2 + pad, row["y"] - h / 2 + pad),
            max(1.0, w - 2 * pad), max(1.0, h - 2 * pad),
            boxstyle=f"round,pad={pad}",
            facecolor=theme.get("compartment", "#E3F2FD"), alpha=alpha,
            edgecolor=theme.get("compartment_edge", "#90CAF9"),
            linewidth=0.9, linestyle="--", zorder=0.5,
        ))
        if row["label"]:
            ax.text(row["x"] - w / 2 + pad * 1.6, row["y"] - h / 2 + pad * 1.6,
                    str(row["label"]), ha="left", va="top",
                    fontsize=label_size, color=theme["muted"],
                    style="italic", zorder=0.6)


def draw_pathway(
    ax,
    node_data: pl.DataFrame,
    edge_data: pl.DataFrame | None = None,
    color_map: dict[str, list[str]] | None = None,
    label_map: dict[str, str] | None = None,
    theme: str | dict = "publication",
    draw_edges: bool = True,
    draw_map_nodes: bool = True,
    show_link_edges: bool = False,
    edge_alpha: float = 0.78,
    label_scale: float = 1.0,
    compartments: pl.DataFrame | None = None,
) -> Extent:
    """
    Draw a whole pathway onto *ax* and return the extent used.

    *color_map* maps entry_id -> list of fill colours (one per experiment).
    """
    th = THEMES.get(theme, THEMES["publication"]) if isinstance(theme, str) else theme
    color_map = color_map or {}
    label_map = label_map or {}

    boxes = [b for b in node_boxes(node_data) if not b.is_title]
    by_id = {b.entry_id: b for b in boxes}
    extent = Extent.from_boxes(boxes, pad=34.0)
    if compartments is not None and not compartments.is_empty():
        extent = _cover_compartments(extent, compartments)

    ax.set_xlim(extent.x0, extent.x1)
    # Inverted y: KGML coordinates point downwards.  The *axis* is flipped,
    # never the data — one place, one rule.
    ax.set_ylim(extent.y1, extent.y0)
    ax.set_facecolor(th["map_bg"])
    ax.set_axis_off()

    if compartments is not None:
        draw_compartments(ax, compartments, th)

    if draw_edges and edge_data is not None and not edge_data.is_empty():
        span = max(extent.width, extent.height)
        for row in edge_data.iter_rows(named=True):
            s = by_id.get(str(row["source"]))
            t = by_id.get(str(row["target"]))
            if s is None or t is None:
                continue
            # Edges into pathway-link nodes are cross-references, not
            # reactions; drawn at full weight they dominate the figure.
            link_edge = s.node_type == "map" or t.node_type == "map"
            if link_edge and not show_link_edges:
                continue
            reach = float(np.hypot(t.x - s.x, t.y - s.y))
            fade = 0.30 if link_edge else edge_alpha
            if reach > 0.55 * span:
                fade *= 0.55
            draw_edge(ax, s, t, row.get("subtype", ""), row.get("edge_type", ""),
                      th, alpha=fade,
                      width=0.65 if link_edge else 0.95,
                      curvature=0.06 if reach > 0.4 * span else 0.13)

    # Pathway-link nodes first, so they sit behind the biology.
    for box in boxes:
        if box.node_type == "map":
            if not draw_map_nodes:
                continue
            draw_node(ax, box, ["#ECEFF1"], th,
                      label_override=label_map.get(box.entry_id),
                      border_width=0.5, label_scale=label_scale * 0.95)

    for box in boxes:
        if box.node_type == "map":
            continue
        cols = color_map.get(box.entry_id) or [th["unmapped"]]
        draw_node(ax, box, cols, th,
                  label_override=label_map.get(box.entry_id),
                  label_scale=label_scale)

    return extent


def render_vector_array(
    node_data: pl.DataFrame,
    edge_data: pl.DataFrame | None = None,
    color_map: dict[str, list[str]] | None = None,
    label_map: dict[str, str] | None = None,
    theme: str = "publication",
    px_per_unit: float = 1.6,
    draw_edges: bool = True,
    show_link_edges: bool = False,
    compartments: pl.DataFrame | None = None,
):
    """
    Render only the pathway (no title, no keys, no padding) to a
    :class:`~pathview.layout.RasterFrame`.

    The axes fill the canvas exactly, so raster pixels map to KGML coordinates
    by a single known scale.  That is what allows highlighting to be applied
    to vector-mode output at the right place.
    """
    import matplotlib
    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    from .layout import RasterFrame

    th = THEMES.get(theme, THEMES["publication"])
    boxes = [b for b in node_boxes(node_data) if not b.is_title]
    extent = Extent.from_boxes(boxes, pad=34.0)
    if compartments is not None and not compartments.is_empty():
        extent = _cover_compartments(extent, compartments)

    dpi = 100.0
    fig_w = max(1.0, extent.width * px_per_unit / dpi)
    fig_h = max(1.0, extent.height * px_per_unit / dpi)

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi, facecolor=th["map_bg"])
    ax = fig.add_axes([0, 0, 1, 1])          # axes fill the canvas exactly
    draw_pathway(ax, node_data, edge_data, color_map, label_map, theme=th,
                 draw_edges=draw_edges, show_link_edges=show_link_edges,
                 compartments=compartments)
    ax.set_xlim(extent.x0, extent.x1)
    ax.set_ylim(extent.y1, extent.y0)

    fig.canvas.draw()
    arr = np.asarray(fig.canvas.buffer_rgba(), dtype=np.uint8)[..., :3].copy()
    plt.close(fig)

    scale = arr.shape[1] / extent.width
    return RasterFrame(arr, x0=extent.x0, y0=extent.y0, scale=scale)


def keggview_vector(
    node_data: pl.DataFrame,
    edge_data: pl.DataFrame | None = None,
    color_map: dict[str, list[str]] | None = None,
    label_map: dict[str, str] | None = None,
    pathway_name: str = "pathway",
    title: str | None = None,
    subtitle: str | None = None,
    out_dir: str | Path = ".",
    out_suffix: str = "pathview",
    output_format: str = "pdf",
    gene_scale: ColorScale | None = None,
    cpd_scale: ColorScale | None = None,
    theme: str = "publication",
    dpi: int = 220,
    figure_width: float = 14.0,
    new_signature: bool = True,
    plot_col_key: bool = True,
    draw_edges: bool = True,
    show_link_edges: bool = False,
    compartments: pl.DataFrame | None = None,
) -> Path:
    """
    Render a pathway as publication-quality vector output.

    Returns the path written.
    """
    import matplotlib
    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    th = THEMES.get(theme, THEMES["publication"])
    boxes = [b for b in node_boxes(node_data) if not b.is_title]
    extent = Extent.from_boxes(boxes, pad=34.0)
    if compartments is not None and not compartments.is_empty():
        extent = _cover_compartments(extent, compartments)
    fig_w, fig_h = figure_size(extent, target_width_in=figure_width)

    keys = [s for s in (gene_scale, cpd_scale) if s is not None] if plot_col_key else []
    key_h = 0.62 if keys else 0.0
    head_h = 0.62 if title else 0.16

    total_h = fig_h + key_h + head_h
    fig = plt.figure(figsize=(fig_w, total_h), facecolor=th["bg"])

    ax = fig.add_axes([0.012, (key_h + 0.10) / total_h, 0.976,
                       fig_h / total_h])
    draw_pathway(ax, node_data, edge_data, color_map, label_map,
                 theme=th, draw_edges=draw_edges,
                 show_link_edges=show_link_edges, compartments=compartments)

    if title:
        fig.text(0.5, 1 - 0.18 / total_h, title, ha="center", va="top",
                 fontsize=14, fontweight="semibold", color=th["title"])
    if subtitle:
        fig.text(0.5, 1 - 0.42 / total_h, subtitle, ha="center", va="top",
                 fontsize=8.6, color=th["muted"])

    if keys:
        draw_dual_key(fig, gene_scale, cpd_scale,
                      rect=(0.14, 0.30 / total_h, 0.72, 0.20 / total_h),
                      label_size=7.6)

    if new_signature:
        fig.text(0.992, 0.09 / total_h, "pathview-plus", ha="right", va="bottom",
                 fontsize=6.4, color=th["muted"], style="italic")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{pathway_name}.{out_suffix}.{output_format}"
    fig.savefig(out_path, dpi=dpi, facecolor=th["bg"],
                bbox_inches="tight", pad_inches=0.14,
                transparent=False)
    plt.close(fig)
    return out_path
