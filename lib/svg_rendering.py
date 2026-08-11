"""
svg_rendering.py
Standalone SVG output with edges, markers and an embedded colour key.

Fixes over v2.x
---------------
* Edges were never drawn.  ``render_edge_svg`` existed but ``keggview_svg``
  never called it, so the SVG was a field of disconnected boxes.
* ``render_edge_svg`` emitted a fresh ``<defs><marker id="marker_arrow">``
  every call, producing hundreds of duplicate element ids in one document —
  invalid SVG that renderers resolve unpredictably.  Markers are now defined
  once, in the document header.
* Node colours were looked up with ``cols.filter(pl.col("id") == node_id)``
  *inside* a loop over nodes, and the filter ran twice per colour column:
  O(nodes x columns) full scans.  Lookups are now a single dict build.
* The canvas ignored node extents, clipping the rightmost and bottom nodes.
* No colour key was emitted at all.

Public API
----------
  keggview_svg    : write a complete SVG document
  render_node_svg : SVG for a single node
  render_edge_svg : SVG for a single edge
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import polars as pl

from .color_mapping import ColorScale
from .layout import Extent, NodeBox, fit_label, node_boxes, slice_bounds
from .splines import offset_endpoints, points_to_bezier_path, route_edge_spline
from .utils import contrast_text_color, escape_xml, is_transparent, to_hex
from .vector_rendering import _DEFAULT_EDGE, EDGE_STYLE, THEMES

_MARKERS = ("arrow", "bar", "circle", "diamond")


def _defs(theme: dict) -> str:
    """Marker and filter definitions, emitted once per document."""
    out = ["  <defs>"]
    for name in _MARKERS:
        shapes = {
            "arrow": '<path d="M0,0 L10,5 L0,10 z"/>',
            "bar": '<path d="M8,0 L10,0 L10,10 L8,10 z"/>',
            "circle": '<circle cx="5" cy="5" r="4"/>',
            "diamond": '<path d="M0,5 L5,0 L10,5 L5,10 z"/>',
        }[name]
        out.append(
            f'    <marker id="pv-{name}" viewBox="0 0 10 10" refX="9" refY="5" '
            f'markerWidth="5" markerHeight="5" orient="auto-start-reverse" '
            f'markerUnits="strokeWidth">{shapes}</marker>'
        )
    out.append(
        '    <filter id="pv-halo" x="-20%" y="-20%" width="140%" height="140%">'
        '<feDropShadow dx="0" dy="0.4" stdDeviation="0.5" '
        'flood-color="#000" flood-opacity="0.18"/></filter>'
    )
    out.append("  </defs>")
    return "\n".join(out)


def _header(extent: Extent, title: str, theme: dict) -> str:
    w, h = extent.width, extent.height
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{w:.0f}" height="{h:.0f}" '
        f'viewBox="{extent.x0:.1f} {extent.y0:.1f} {w:.1f} {h:.1f}">\n'
        f'  <title>{escape_xml(title)}</title>\n'
        '  <style type="text/css"><![CDATA[\n'
        '    .pv-node { stroke-width: 0.7; }\n'
        f'    .pv-label {{ font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; '
        f'text-anchor: middle; fill: {theme["text"]}; }}\n'
        '    .pv-edge { fill: none; stroke-linecap: round; }\n'
        '  ]]></style>\n'
        f'  <rect x="{extent.x0:.1f}" y="{extent.y0:.1f}" width="{w:.1f}" '
        f'height="{h:.1f}" fill="{theme["map_bg"]}"/>\n'
        + _defs(theme) + "\n"
    )


def render_node_svg(box: NodeBox, fill_colors: Sequence[str], theme: dict,
                    label_override: str | None = None) -> str:
    """SVG fragment for one node, sliced by experiment."""
    fills = [theme["unmapped"] if is_transparent(c) else to_hex(c, False)
             for c in (fill_colors or [])] or [theme["unmapped"]]
    bands = slice_bounds(box, len(fills))
    nid = escape_xml(box.entry_id)
    parts: list[str] = [f'  <g id="node-{nid}">']

    if box.is_round:
        cid = f"clip-{nid}"
        parts.append(
            f'    <clipPath id="{cid}"><circle cx="{box.x:.1f}" cy="{box.y:.1f}" '
            f'r="{box.radius:.1f}"/></clipPath>'
        )
        parts.append(f'    <g clip-path="url(#{cid})">')
        for (x0, x1), face in zip(bands, fills):
            parts.append(
                f'      <rect x="{x0:.1f}" y="{box.top:.1f}" '
                f'width="{max(0.1, x1 - x0):.1f}" height="{box.height:.1f}" fill="{face}"/>'
            )
        parts.append("    </g>")
        parts.append(
            f'    <circle cx="{box.x:.1f}" cy="{box.y:.1f}" r="{box.radius:.1f}" '
            f'fill="none" stroke="{theme["border"]}" class="pv-node"/>'
        )
    else:
        r = box.corner_radius
        for (x0, x1), face in zip(bands, fills):
            rr = f' rx="{r:.1f}"' if len(fills) == 1 else ""
            parts.append(
                f'    <rect x="{x0:.1f}" y="{box.top:.1f}" '
                f'width="{max(0.1, x1 - x0):.1f}" height="{box.height:.1f}"{rr} fill="{face}"/>'
            )
        parts.append(
            f'    <rect x="{box.left:.1f}" y="{box.top:.1f}" width="{box.width:.1f}" '
            f'height="{box.height:.1f}" rx="{r:.1f}" fill="none" '
            f'stroke="{theme["border"]}" class="pv-node"/>'
        )

    label = label_override if label_override is not None else box.label
    if label:
        text = fit_label(box, label)
        lines = text.split("\n")
        if box.is_round:
            size = 4.6
            y = box.bottom + size + 0.8
            fill = theme["text"]
        else:
            size = max(3.4, min(6.5, (box.width * 1.5) / max(1, max(len(t) for t in lines))))
            y = box.y - (len(lines) - 1) * size * 0.55 + size * 0.35
            fill = contrast_text_color(fills[len(fills) // 2], dark=theme["text"])
        for line in lines:
            parts.append(
                f'    <text x="{box.x:.1f}" y="{y:.1f}" class="pv-label" '
                f'font-size="{size:.1f}" fill="{fill}">{escape_xml(line)}</text>'
            )
            y += size * 1.1

    parts.append("  </g>")
    return "\n".join(parts)


def render_edge_svg(src: NodeBox, tgt: NodeBox, subtype: str = "",
                    edge_type: str = "", width: float = 0.9,
                    opacity: float = 0.75) -> str:
    """SVG fragment for one edge; markers come from the shared ``<defs>``."""
    style = EDGE_STYLE.get((subtype or "").lower(),
                           EDGE_STYLE.get((edge_type or "").lower(), _DEFAULT_EDGE))
    s, t = offset_endpoints(
        (src.x, src.y), (tgt.x, tgt.y),
        src.radius if src.is_round else src.height / 2 + 1,
        (tgt.radius + 3.0) if tgt.is_round else tgt.height / 2 + 3.0,
    )
    curve = route_edge_spline(s, t, routing_mode="curved", curvature=0.12, n_points=14)
    path = points_to_bezier_path([tuple(p) for p in curve])
    if not path:
        return ""

    dash = {"-": "", "--": ' stroke-dasharray="4,2.5"', ":": ' stroke-dasharray="1.5,2"'}
    marker = ""
    if style["head"] != "none":
        marker = f' marker-end="url(#pv-{style["head"]})"'
    return (
        f'  <path d="{path}" class="pv-edge" stroke="{style["color"]}" '
        f'stroke-width="{width:.2f}" opacity="{opacity:.2f}"'
        f'{dash.get(style["style"], "")}{marker}/>'
    )


def _color_key_svg(scale: ColorScale, x: float, y: float,
                   width: float, height: float, theme: dict) -> str:
    """An inline colour key so the SVG is self-describing."""
    colors = scale.colors()
    n = len(colors)
    lo, hi = scale.bounds()
    parts = ['  <g class="pv-key">']
    bw = width / n
    for i, c in enumerate(colors):
        parts.append(f'    <rect x="{x + i * bw:.1f}" y="{y:.1f}" width="{bw:.2f}" '
                     f'height="{height:.1f}" fill="{c}"/>')
    parts.append(f'    <rect x="{x:.1f}" y="{y:.1f}" width="{width:.1f}" '
                 f'height="{height:.1f}" fill="none" stroke="{theme["border"]}" '
                 f'stroke-width="0.6"/>')
    for frac, val in ((0.0, lo), (0.5, (lo + hi) / 2), (1.0, hi)):
        parts.append(
            f'    <text x="{x + frac * width:.1f}" y="{y + height + 7:.1f}" '
            f'class="pv-label" font-size="6">{val:g}</text>'
        )
    if scale.label:
        parts.append(f'    <text x="{x + width / 2:.1f}" y="{y - 3:.1f}" '
                     f'class="pv-label" font-size="7" '
                     f'font-weight="600">{escape_xml(scale.label)}</text>')
    parts.append("  </g>")
    return "\n".join(parts)


def keggview_svg(
    node_data: pl.DataFrame,
    edge_data: pl.DataFrame | None = None,
    color_map: dict[str, list[str]] | None = None,
    label_map: dict[str, str] | None = None,
    pathway_name: str = "pathway",
    title: str | None = None,
    out_dir: str | Path = ".",
    out_suffix: str = "pathview",
    gene_scale: ColorScale | None = None,
    cpd_scale: ColorScale | None = None,
    theme: str = "publication",
    new_signature: bool = True,
    plot_col_key: bool = True,
    draw_edges: bool = True,
    show_link_edges: bool = False,
) -> Path:
    """Write a complete, self-contained SVG document."""
    th = THEMES.get(theme, THEMES["publication"])
    color_map = color_map or {}
    label_map = label_map or {}

    boxes = [b for b in node_boxes(node_data) if not b.is_title]
    by_id = {b.entry_id: b for b in boxes}

    keys = [s for s in (gene_scale, cpd_scale) if s is not None] if plot_col_key else []
    extent = Extent.from_boxes(boxes, pad=34.0)
    if keys:
        extent = Extent(extent.x0, extent.y0, extent.x1, extent.y1 + 46.0)

    parts = [_header(extent, title or pathway_name, th)]

    if draw_edges and edge_data is not None and not edge_data.is_empty():
        parts.append("  <!-- edges -->")
        for row in edge_data.iter_rows(named=True):
            s, t = by_id.get(str(row["source"])), by_id.get(str(row["target"]))
            if s is None or t is None:
                continue
            if (s.node_type == "map" or t.node_type == "map") and not show_link_edges:
                continue
            frag = render_edge_svg(s, t, row.get("subtype", ""), row.get("edge_type", ""))
            if frag:
                parts.append(frag)

    parts.append("  <!-- nodes -->")
    for box in sorted(boxes, key=lambda b: 0 if b.node_type == "map" else 1):
        parts.append(render_node_svg(box, color_map.get(box.entry_id, []), th,
                                     label_map.get(box.entry_id)))

    if keys:
        kw = min(240.0, extent.width * 0.32)
        y = extent.y1 - 30.0
        x = extent.x0 + extent.width * 0.10
        for i, sc in enumerate(keys):
            parts.append(_color_key_svg(sc, x + i * (kw + 40.0), y, kw, 8.0, th))

    if new_signature:
        parts.append(
            f'  <text x="{extent.x1 - 6:.1f}" y="{extent.y1 - 5:.1f}" '
            f'class="pv-label" text-anchor="end" font-size="6" '
            f'fill="{th["muted"]}" font-style="italic">pathview-plus</text>'
        )

    parts.append("</svg>")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{pathway_name}.{out_suffix}.svg"
    out_path.write_text("\n".join(parts), encoding="utf-8")
    return out_path
