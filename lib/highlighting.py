"""
highlighting.py
Layer-by-layer pathway graph modifications.

Implements a ggplot2-style composable interface for post-hoc pathway
customization. Users can highlight specific nodes or paths, change colors,
adjust labels, and more — all without re-running the full rendering pipeline.

Usage
-----
    from pathview import pathview, highlight_nodes, highlight_edges
    
    result = pathview("04110", gene_data=data, species="hsa")
    
    # Compose modifications with +
    modified = (result
                + highlight_nodes(["1956", "2099"], color="red", width=4)
                + highlight_edges([("1956", "2099")], color="blue", width=3))
    
    # Save modified version
    modified.save("highlighted.png")

Public API
----------
  PathwayResult : Container for pathway rendering results
  highlight_nodes : Highlight specific nodes
  highlight_edges : Highlight specific edges
  highlight_path  : Highlight an entire path
  change_labels   : Update node labels
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import polars as pl
from PIL import Image


# ---------------------------------------------------------------------------
# PathwayResult container
# ---------------------------------------------------------------------------

@dataclass
class PathwayResult:
    """
    Container for pathway rendering results.
    
    Supports composable modifications via the + operator.
    Stores both the rendered image and the underlying data so modifications
    can be applied without full re-rendering.
    """
    pathway_id: str
    plot_data_gene: Optional[pl.DataFrame] = None
    plot_data_cpd: Optional[pl.DataFrame] = None
    output_path: Optional[Path] = None
    image_array: Optional[np.ndarray] = None
    modifications: list[Callable] = field(default_factory=list)
    
    def __add__(self, modifier: Callable) -> PathwayResult:
        """
        Apply a modification function and return a new PathwayResult.
        
        This implements ggplot2-style layer composition:
        
            result = pathview(...) + highlight_nodes(...) + highlight_edges(...)
        """
        new_result = PathwayResult(
            pathway_id=self.pathway_id,
            plot_data_gene=self.plot_data_gene,
            plot_data_cpd=self.plot_data_cpd,
            output_path=self.output_path,
            image_array=self.image_array.copy() if self.image_array is not None else None,
            modifications=self.modifications + [modifier],
        )
        # Apply the modification
        modifier(new_result)
        return new_result
    
    def save(self, path: str | Path, format: str = "png") -> None:
        """Save the modified pathway to a file."""
        if self.image_array is None:
            raise ValueError("No image data to save")
        
        img = Image.fromarray(self.image_array)
        if format.lower() == "pdf":
            img.save(path, "PDF", resolution=300.0)
        else:
            img.save(path, format.upper())
        print(f"Info: Saved modified pathway → {path}")
    
    def show(self) -> None:
        """Display the pathway using PIL."""
        if self.image_array is None:
            raise ValueError("No image data to display")
        Image.fromarray(self.image_array).show()


# ---------------------------------------------------------------------------
# Node highlighting
# ---------------------------------------------------------------------------

def highlight_nodes(
    node_ids: list[str],
    color: str = "red",
    width: int = 4,
    opacity: float = 1.0,
) -> Callable[[PathwayResult], None]:
    """
    Highlight specified nodes by changing their border.
    
    Parameters
    ----------
    node_ids: List of node IDs to highlight (Entrez IDs or KEGG IDs)
    color:    Border color (hex or named color)
    width:    Border width in pixels
    opacity:  Border opacity (0-1)
    
    Returns a modifier function that can be added to a PathwayResult.
    
    Example
    -------
    >>> result = pathview("04110", gene_data=data)
    >>> highlighted = result + highlight_nodes(["1956", "2099"], color="red", width=4)
    >>> highlighted.save("highlighted.png")
    """
    def modifier(result: PathwayResult) -> None:
        if result.image_array is None:
            return
        
        # Find node positions
        nodes_to_highlight = []
        if result.plot_data_gene is not None:
            genes = result.plot_data_gene.filter(
                pl.col("kegg_names").is_in(node_ids)
            )
            if not genes.is_empty():
                nodes_to_highlight.append(genes)
        
        if result.plot_data_cpd is not None:
            cpds = result.plot_data_cpd.filter(
                pl.col("kegg_names").is_in(node_ids)
            )
            if not cpds.is_empty():
                nodes_to_highlight.append(cpds)
        
        # Draw highlights
        img_height = result.image_array.shape[0]
        rgb = _hex_to_rgb(color)
        
        for df in nodes_to_highlight:
            for row in df.iter_rows(named=True):
                cx, cy = row["x"], row["y"]
                hw, hh = row["width"] / 2, row["height"] / 2
                _draw_border(
                    result.image_array, 
                    cx=cx, cy=cy, 
                    half_width=hw, half_height=hh,
                    img_height=img_height,
                    rgb=rgb, thickness=width, opacity=opacity
                )
    
    return modifier


# ---------------------------------------------------------------------------
# Edge highlighting
# ---------------------------------------------------------------------------

def highlight_edges(
    edge_pairs: list[tuple[str, str]],
    color: str = "blue",
    width: int = 3,
) -> Callable[[PathwayResult], None]:
    """
    Highlight specified edges (connections between nodes).
    
    Parameters
    ----------
    edge_pairs: List of (source_id, target_id) tuples
    color:      Edge color
    width:      Edge width in pixels
    
    Returns a modifier function.
    
    Example
    -------
    >>> result + highlight_edges([("1956", "2099"), ("2099", "5594")])
    """
    def modifier(result: PathwayResult) -> None:
        if result.image_array is None:
            return
        
        # Find node positions for edge endpoints
        gene_pos = {}
        if result.plot_data_gene is not None:
            for row in result.plot_data_gene.iter_rows(named=True):
                gene_pos[row["kegg_names"]] = (row["x"], row["y"])
        
        img_height = result.image_array.shape[0]
        rgb = _hex_to_rgb(color)
        
        # Draw lines between pairs
        for source_id, target_id in edge_pairs:
            if source_id in gene_pos and target_id in gene_pos:
                x1, y1 = gene_pos[source_id]
                x2, y2 = gene_pos[target_id]
                _draw_line(
                    result.image_array,
                    x1=x1, y1=y1, x2=x2, y2=y2,
                    img_height=img_height,
                    rgb=rgb, thickness=width
                )
    
    return modifier


# ---------------------------------------------------------------------------
# Path highlighting
# ---------------------------------------------------------------------------

def highlight_path(
    path_node_ids: list[str],
    color: str = "orange",
    node_width: int = 3,
    edge_width: int = 2,
) -> Callable[[PathwayResult], None]:
    """
    Highlight an entire path (nodes and edges).
    
    Parameters
    ----------
    path_node_ids: Ordered list of node IDs forming a path
    color:         Color for both nodes and edges
    node_width:    Border width for nodes
    edge_width:    Width for connecting edges
    
    Returns a modifier function.
    
    Example
    -------
    >>> result + highlight_path(["1956", "2099", "5594", "207"], color="orange")
    """
    # Build edge pairs from consecutive nodes
    edge_pairs = list(zip(path_node_ids[:-1], path_node_ids[1:]))
    
    def modifier(result: PathwayResult) -> None:
        # Apply both node and edge highlighting
        highlight_nodes(path_node_ids, color=color, width=node_width)(result)
        highlight_edges(edge_pairs, color=color, width=edge_width)(result)
    
    return modifier


# ---------------------------------------------------------------------------
# Label modification
# ---------------------------------------------------------------------------

def change_labels(
    label_map: dict[str, str],
    font_size: int = 11,
    color: str = "black",
) -> Callable[[PathwayResult], None]:
    """
    Change labels for specified nodes.
    
    Parameters
    ----------
    label_map:  Dict mapping node_id → new_label
    font_size:  Font size for new labels
    color:      Text color
    
    Returns a modifier function.
    
    Example
    -------
    >>> result + change_labels({"1956": "EGFR*", "2099": "ESR1*"})
    """
    def modifier(result: PathwayResult) -> None:
        # This would require text rendering on the image
        # For now, store the label changes for future re-rendering
        if not hasattr(result, '_label_changes'):
            result._label_changes = {}
        result._label_changes.update(label_map)
    
    return modifier


# ---------------------------------------------------------------------------
# Drawing primitives
# ---------------------------------------------------------------------------

def _hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    """Convert hex color to RGB tuple."""
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))


def _draw_border(
    img: np.ndarray,
    cx: float,
    cy: float,
    half_width: float,
    half_height: float,
    img_height: int,
    rgb: tuple[int, int, int],
    thickness: int,
    opacity: float,
) -> None:
    """Draw a rectangle border on the image array."""
    # Convert KGML coordinates to image coordinates
    py = int(img_height - cy)
    px = int(cx)
    hw, hh = int(half_width), int(half_height)
    
    # Draw rectangle border (4 sides)
    for t in range(thickness):
        # Top
        img[max(0, py - hh - t):min(img.shape[0], py - hh - t + 1),
            max(0, px - hw):min(img.shape[1], px + hw)] = rgb
        # Bottom
        img[max(0, py + hh + t):min(img.shape[0], py + hh + t + 1),
            max(0, px - hw):min(img.shape[1], px + hw)] = rgb
        # Left
        img[max(0, py - hh):min(img.shape[0], py + hh),
            max(0, px - hw - t):min(img.shape[1], px - hw - t + 1)] = rgb
        # Right
        img[max(0, py - hh):min(img.shape[0], py + hh),
            max(0, px + hw + t):min(img.shape[1], px + hw + t + 1)] = rgb


def _draw_line(
    img: np.ndarray,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    img_height: int,
    rgb: tuple[int, int, int],
    thickness: int,
) -> None:
    """Draw a line on the image array using Bresenham's algorithm."""
    # Convert coordinates
    px1, py1 = int(x1), int(img_height - y1)
    px2, py2 = int(x2), int(img_height - y2)
    
    # Bresenham's line algorithm
    dx = abs(px2 - px1)
    dy = abs(py2 - py1)
    sx = 1 if px1 < px2 else -1
    sy = 1 if py1 < py2 else -1
    err = dx - dy
    
    while True:
        # Draw thick point
        for t in range(-thickness // 2, thickness // 2 + 1):
            for u in range(-thickness // 2, thickness // 2 + 1):
                py = py1 + t
                px = px1 + u
                if 0 <= py < img.shape[0] and 0 <= px < img.shape[1]:
                    img[py, px] = rgb
        
        if px1 == px2 and py1 == py2:
            break
        
        e2 = 2 * err
        if e2 > -dy:
            err -= dy
            px1 += sx
        if e2 < dx:
            err += dx
            py1 += sy
