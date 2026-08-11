"""
layout.py
Coordinate handling and node geometry.

The v2.x y-flip bug
-------------------
KGML ``(x, y)`` are already **image pixel coordinates with the origin at the
top-left**, identical to the coordinate system of the KEGG PNG they annotate.
No flip is ever required.

v2.x flipped inconsistently, in three different ways in three places:

  * ``_paint_gene_nodes``  flipped twice (``cy = h - cy`` then ``h - cy``),
    which cancels out — accidentally correct, and therefore never noticed.
  * ``_paint_cpd_nodes``   flipped once, so every compound was painted at the
    vertically mirrored position.  On a metabolic map that puts every
    metabolite in the wrong place.
  * ``highlighting._draw_border`` / ``_draw_line`` also flipped once, so
    highlights landed nowhere near the nodes they were highlighting.

This module centralises the geometry so there is exactly one definition, and
the renderers never do arithmetic on coordinates themselves.

Public API
----------
  NodeBox        : resolved geometry for one node
  node_boxes     : DataFrame -> list[NodeBox]
  slice_bounds   : split a node into per-experiment slices
  Extent         : canvas extent helper
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import polars as pl

from .constants import DEFAULT_CPD_RADIUS, DEFAULT_GENE_HEIGHT, DEFAULT_GENE_WIDTH

#: KGML coordinates are top-left origin, y increasing downwards, same as the
#: KEGG PNG.  Renderers invert the *axis*, never the data.
Y_AXIS_POINTS_DOWN = True


@dataclass(frozen=True)
class NodeBox:
    """Resolved geometry and styling for a single node."""

    entry_id: str
    x: float                     # centre x, KGML pixel space
    y: float                     # centre y, KGML pixel space (origin top-left)
    width: float
    height: float
    shape: str = "rectangle"
    node_type: str = "gene"
    label: str = ""
    bgcolor: str = "#FFFFFF"

    @property
    def left(self) -> float:
        return self.x - self.width / 2.0

    @property
    def right(self) -> float:
        return self.x + self.width / 2.0

    @property
    def top(self) -> float:
        return self.y - self.height / 2.0

    @property
    def bottom(self) -> float:
        return self.y + self.height / 2.0

    @property
    def radius(self) -> float:
        """
        Radius for round nodes.

        Half the smaller dimension.  v2.x used the full ``width`` as the
        radius, drawing every compound circle at twice its true size so
        neighbouring metabolites overlapped and swallowed their labels.
        """
        return min(self.width, self.height) / 2.0

    @property
    def is_round(self) -> bool:
        """
        True only for genuinely circular glyphs.

        ``roundrectangle`` is a *rectangle with rounded corners* in KGML — it
        is what KEGG uses for pathway-link nodes.  Treating it as a circle
        (as the shape list originally did) drew every pathway link as a large
        disc and destroyed the map's visual grammar.
        """
        return self.shape in ("circle", "ellipse", "round") \
            or self.node_type == "compound"

    @property
    def corner_radius(self) -> float:
        """Corner radius for rectangular glyphs, 0 for hard corners."""
        if self.shape == "roundrectangle" or self.node_type == "map":
            return min(6.0, self.height / 2.4)
        return min(2.8, self.height / 5.0)

    @property
    def is_title(self) -> bool:
        """KEGG emits a pseudo-entry carrying the map title; it is not a node."""
        return self.label.upper().startswith("TITLE:")

    def contains(self, px: float, py: float) -> bool:
        if self.is_round:
            r = self.radius
            return ((px - self.x) ** 2 + (py - self.y) ** 2) <= r * r
        return self.left <= px <= self.right and self.top <= py <= self.bottom


def node_boxes(
    df: pl.DataFrame | None,
    default_type: str = "gene",
) -> list[NodeBox]:
    """Convert a node DataFrame into geometry objects, skipping unplaced nodes."""
    if df is None or df.is_empty():
        return []

    boxes: list[NodeBox] = []
    for row in df.iter_rows(named=True):
        x, y = row.get("x"), row.get("y")
        if x is None or y is None:
            continue
        ntype = row.get("type") or default_type
        is_cpd = ntype == "compound"
        w = row.get("width") or (DEFAULT_CPD_RADIUS * 2 if is_cpd else DEFAULT_GENE_WIDTH)
        h = row.get("height") or (DEFAULT_CPD_RADIUS * 2 if is_cpd else DEFAULT_GENE_HEIGHT)
        boxes.append(NodeBox(
            entry_id=str(row.get("entry_id", "")),
            x=float(x), y=float(y), width=float(w), height=float(h),
            shape=str(row.get("shape") or ("circle" if is_cpd else "rectangle")),
            node_type=str(ntype),
            label=str(row.get("label") or ""),
            bgcolor=str(row.get("bgcolor") or "#FFFFFF"),
        ))
    return boxes


def slice_bounds(box: NodeBox, n_slices: int) -> list[tuple[float, float]]:
    """
    Split a node horizontally into *n_slices* equal (x0, x1) bands.

    One band per experiment / condition, which is how a multi-condition
    comparison is shown on a single node.
    """
    n = max(1, int(n_slices))
    edges = np.linspace(box.left, box.right, n + 1)
    return [(float(edges[i]), float(edges[i + 1])) for i in range(n)]


def slice_angles(n_slices: int, start: float = 90.0) -> list[tuple[float, float]]:
    """Split a circle into *n_slices* pie wedges, in degrees."""
    n = max(1, int(n_slices))
    step = 360.0 / n
    return [(start + i * step, start + (i + 1) * step) for i in range(n)]


@dataclass
class RasterFrame:
    """
    A rendered raster together with its mapping to KGML coordinates.

    Highlighting draws in KGML pixel space, but a saved *figure* has a title
    band, colour keys and padding, and is scaled by dpi — so KGML (x, y) is
    not figure (px, py).  Storing the affine relationship makes the two
    unambiguous instead of leaving them to coincidentally agree, which is how
    highlights end up drawn nowhere near the nodes they mark.
    """

    array: np.ndarray
    x0: float = 0.0          # KGML x at raster column 0
    y0: float = 0.0          # KGML y at raster row 0
    scale: float = 1.0       # raster pixels per KGML unit

    def to_pixels(self, x: float, y: float) -> tuple[float, float]:
        """KGML coordinates -> raster pixel coordinates."""
        return ((x - self.x0) * self.scale, (y - self.y0) * self.scale)

    def length(self, value: float) -> float:
        """Scale a KGML distance into raster pixels."""
        return value * self.scale

    @property
    def shape(self) -> tuple[int, ...]:
        return self.array.shape

    def copy(self) -> RasterFrame:
        return RasterFrame(self.array.copy(), self.x0, self.y0, self.scale)


@dataclass
class Extent:
    """Canvas extent in KGML pixel space."""

    x0: float = 0.0
    y0: float = 0.0
    x1: float = 1200.0
    y1: float = 900.0

    @property
    def width(self) -> float:
        return max(1.0, self.x1 - self.x0)

    @property
    def height(self) -> float:
        return max(1.0, self.y1 - self.y0)

    @property
    def aspect(self) -> float:
        return self.width / self.height

    def padded(self, pad: float) -> Extent:
        return Extent(self.x0 - pad, self.y0 - pad, self.x1 + pad, self.y1 + pad)

    @classmethod
    def from_boxes(cls, boxes: Sequence[NodeBox], pad: float = 30.0) -> Extent:
        if not boxes:
            return cls()
        return cls(
            min(b.left for b in boxes) - pad,
            min(b.top for b in boxes) - pad,
            max(b.right for b in boxes) + pad,
            max(b.bottom for b in boxes) + pad,
        )

    @classmethod
    def from_image(cls, width: int, height: int) -> Extent:
        return cls(0.0, 0.0, float(width), float(height))


def figure_size(extent: Extent, target_width_in: float = 14.0,
                max_height_in: float = 18.0) -> tuple[float, float]:
    """Figure size in inches preserving the pathway's aspect ratio."""
    w = float(target_width_in)
    h = w / max(0.05, extent.aspect)
    if h > max_height_in:
        h = max_height_in
        w = h * extent.aspect
    return (round(w, 2), round(h, 2))


def label_font_size(box: NodeBox, base: float = 6.5, min_size: float = 3.2) -> float:
    """Scale a label to the node it sits in."""
    if not box.label:
        return base
    longest = max((len(line) for line in box.label.split("\n")), default=1)
    fitted = (box.width * 1.55) / max(1, longest)
    return float(max(min_size, min(base, fitted)))


def fit_label(box: NodeBox, label: str, max_lines: int = 2,
              width_scale: float = 1.0) -> str:
    """
    Wrap a label into the node's width, truncating with an ellipsis.

    Round nodes get a generous allowance because their label is drawn
    *outside* the glyph: a KEGG compound circle is 8 px across, so fitting
    text to the node width truncated every metabolite name to "C159…".
    """
    from .utils import wordwrap

    if not label:
        return ""
    effective = box.width * width_scale
    if box.is_round:
        effective = max(effective, 62.0)
    chars = max(4, int(effective / 4.2))
    wrapped = wordwrap(label, width=chars).split("\n")
    if len(wrapped) > max_lines:
        wrapped = wrapped[:max_lines]
        wrapped[-1] = wrapped[-1][:max(1, chars - 1)] + "…"
    return "\n".join(wrapped)
