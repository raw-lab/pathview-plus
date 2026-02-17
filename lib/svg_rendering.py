"""
svg_rendering.py
Generate SVG (Scalable Vector Graphics) output for KEGG pathways.

This complements the existing PNG (pixel-based) and PDF (graph-based) modes
with native vector graphics that are web-friendly, scalable, and editable.

Public API
----------
  keggview_svg    : Render pathway as SVG with data overlay
  render_node_svg : Generate SVG code for a single node
  render_edge_svg : Generate SVG code for a single edge
  
SVG advantages over PNG:
  - Scalable without quality loss
  - Smaller file size for simple diagrams
  - Editable in vector graphics software
  - Web-native format (no conversion needed)
  - Supports CSS styling and JavaScript interaction
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional
from xml.etree import ElementTree as ET

import polars as pl

from .color_mapping import draw_color_key, make_colormap
from .utils import wordwrap


# ---------------------------------------------------------------------------
# SVG header and footer
# ---------------------------------------------------------------------------

def _svg_header(width: int, height: int, title: str = "") -> str:
    """Generate SVG document header."""
    return f'''<?xml version="1.0" encoding="UTF-8"?>
<svg width="{width}" height="{height}" 
     xmlns="http://www.w3.org/2000/svg" 
     xmlns:xlink="http://www.w3.org/1999/xlink"
     viewBox="0 0 {width} {height}">
  <title>{title}</title>
  <defs>
    <style type="text/css">
      .node {{ stroke: #333; stroke-width: 1; }}
      .edge {{ stroke: #666; stroke-width: 1; fill: none; }}
      .label {{ font-family: Arial, sans-serif; font-size: 11px; fill: #000; text-anchor: middle; }}
    </style>
  </defs>
'''

def _svg_footer() -> str:
    """Generate SVG document footer."""
    return "</svg>"


# ---------------------------------------------------------------------------
# Node rendering (rectangles and ellipses)
# ---------------------------------------------------------------------------

def render_node_svg(
    node_id: str,
    x: float,
    y: float,
    width: float,
    height: float,
    shape: str,
    label: str,
    fill_colors: list[str],
    opacity: float = 1.0,
) -> str:
    """
    Render a single node as SVG.
    
    Parameters
    ----------
    node_id:     Unique identifier for the node
    x, y:        Center coordinates
    width, height: Node dimensions
    shape:       "rectangle", "ellipse", "roundedrectangle"
    label:       Text to display on node
    fill_colors: List of hex colors (one per data column/state)
    opacity:     Fill opacity (0-1)
    
    Returns SVG code string for this node.
    """
    svg_parts = []
    n_states = len(fill_colors)
    
    # Calculate bounding box
    x1 = x - width / 2
    y1 = y - height / 2
    
    if shape == "ellipse":
        # Slice ellipse vertically for multi-state
        rx, ry = width / 2, height / 2
        for i, color in enumerate(fill_colors):
            # Create clipped ellipse slices
            clip_x = x1 + (i * width / n_states)
            clip_width = width / n_states
            svg_parts.append(
                f'<g clip-path="url(#clip_{node_id}_{i})">'
                f'  <ellipse cx="{x}" cy="{y}" rx="{rx}" ry="{ry}" '
                f'    fill="{color}" opacity="{opacity}" class="node"/>'
                f'</g>'
            )
            # Define clip path
            svg_parts.insert(
                0,
                f'<clipPath id="clip_{node_id}_{i}">'
                f'  <rect x="{clip_x}" y="{y1}" width="{clip_width}" height="{height}"/>'
                f'</clipPath>'
            )
    else:
        # Rectangle or rounded rectangle
        rx_round = 5 if shape == "roundedrectangle" else 0
        for i, color in enumerate(fill_colors):
            slice_x = x1 + (i * width / n_states)
            slice_width = width / n_states
            svg_parts.append(
                f'<rect x="{slice_x}" y="{y1}" width="{slice_width}" height="{height}" '
                f'rx="{rx_round}" fill="{color}" opacity="{opacity}" class="node"/>'
            )
    
    # Add label
    wrapped = wordwrap(label, width=max(8, int(width / 8)))
    lines = wrapped.split("\n")
    y_text = y - (len(lines) - 1) * 5
    for line in lines:
        svg_parts.append(
            f'<text x="{x}" y="{y_text}" class="label">{_escape_xml(line)}</text>'
        )
        y_text += 12
    
    return "\n".join(svg_parts)


# ---------------------------------------------------------------------------
# Edge rendering
# ---------------------------------------------------------------------------

def render_edge_svg(
    source_x: float,
    source_y: float,
    target_x: float,
    target_y: float,
    edge_type: str = "arrow",
    color: str = "#666",
    width: float = 1.5,
) -> str:
    """
    Render a single edge as SVG.
    
    Parameters
    ----------
    source_x, source_y: Start coordinates
    target_x, target_y: End coordinates
    edge_type:          "arrow", "inhibition", "dotted"
    color:              Stroke color
    width:              Line width
    
    Returns SVG code string for this edge.
    """
    marker_id = f"marker_{edge_type}"
    path_style = f'stroke="{color}" stroke-width="{width}" fill="none" class="edge"'
    
    if edge_type == "dotted":
        path_style += ' stroke-dasharray="3,3"'
    
    # Define arrow markers
    markers = f'''
    <defs>
      <marker id="{marker_id}" viewBox="0 0 10 10" refX="8" refY="5"
              markerWidth="6" markerHeight="6" orient="auto">
        <path d="M 0 0 L 10 5 L 0 10 z" fill="{color}"/>
      </marker>
    </defs>
    '''
    
    # Draw line
    line = f'<line x1="{source_x}" y1="{source_y}" x2="{target_x}" y2="{target_y}" '\
           f'{path_style} marker-end="url(#{marker_id})"/>'
    
    return markers + line


# ---------------------------------------------------------------------------
# Main SVG rendering function
# ---------------------------------------------------------------------------

def keggview_svg(
    plot_data_gene: Optional[pl.DataFrame],
    cols_gene: Optional[pl.DataFrame],
    plot_data_cpd: Optional[pl.DataFrame],
    cols_cpd: Optional[pl.DataFrame],
    node_data: pl.DataFrame,
    pathway_name: str,
    kegg_dir: Path = Path("."),
    out_suffix: str = "pathview",
    new_signature: bool = True,
    **kwargs,
) -> None:
    """
    Render pathway as SVG with data overlay.
    
    This is an alternative to keggview_native (PNG) and keggview_graph (PDF).
    Generates a standalone SVG file with nodes colored by expression data.
    
    Parameters
    ----------
    plot_data_gene:  Gene node positions + data
    cols_gene:       Gene node color assignments
    plot_data_cpd:   Compound node positions + data
    cols_cpd:        Compound node color assignments
    node_data:       All pathway nodes
    pathway_name:    Pathway ID
    kegg_dir:        Output directory
    out_suffix:      Output filename suffix
    new_signature:   Add "Rendered by pathview.py" watermark
    """
    # Determine canvas size from node positions
    max_x = node_data["x"].max() or 1000
    max_y = node_data["y"].max() or 800
    canvas_width = int(max_x + 100)
    canvas_height = int(max_y + 100)
    
    svg_code = [_svg_header(canvas_width, canvas_height, pathway_name)]
    
    # Render gene nodes
    if plot_data_gene is not None and cols_gene is not None:
        svg_code.append("<!-- Gene nodes -->")
        color_cols = [c for c in cols_gene.columns if c.endswith("_col")]
        for row in plot_data_gene.iter_rows(named=True):
            node_id = row["entry_id"]
            colors = [cols_gene.filter(pl.col("id") == node_id)[c].item() 
                     for c in color_cols if not cols_gene.filter(pl.col("id") == node_id).is_empty()]
            if not colors or all(c == "transparent" for c in colors):
                colors = ["#CCCCCC"]
            svg_code.append(
                render_node_svg(
                    node_id=node_id,
                    x=row["x"],
                    y=row["y"],
                    width=row["width"],
                    height=row["height"],
                    shape=row.get("shape", "rectangle"),
                    label=row.get("label", ""),
                    fill_colors=colors,
                )
            )
    
    # Render compound nodes
    if plot_data_cpd is not None and cols_cpd is not None:
        svg_code.append("<!-- Compound nodes -->")
        color_cols = [c for c in cols_cpd.columns if c.endswith("_col")]
        for row in plot_data_cpd.iter_rows(named=True):
            node_id = row["entry_id"]
            colors = [cols_cpd.filter(pl.col("id") == node_id)[c].item() 
                     for c in color_cols if not cols_cpd.filter(pl.col("id") == node_id).is_empty()]
            if not colors or all(c == "transparent" for c in colors):
                colors = ["#DDDDFF"]
            svg_code.append(
                render_node_svg(
                    node_id=node_id,
                    x=row["x"],
                    y=row["y"],
                    width=row["width"],
                    height=row["height"],
                    shape=row.get("shape", "ellipse"),
                    label=row.get("label", ""),
                    fill_colors=colors,
                )
            )
    
    # Add signature
    if new_signature:
        svg_code.append(
            f'<text x="10" y="{canvas_height - 10}" '
            f'style="font-size: 9px; fill: #666;">Rendered by pathview.py</text>'
        )
    
    svg_code.append(_svg_footer())
    
    # Write to file
    out_path = Path(kegg_dir) / f"{pathway_name}.{out_suffix}.svg"
    out_path.write_text("\n".join(svg_code), encoding="utf-8")
    print(f"Info: Written → {out_path}")


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _escape_xml(text: str) -> str:
    """Escape XML special characters."""
    return (text
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace('"', "&quot;")
            .replace("'", "&apos;"))
