"""
rendering.py
Native rendering: overlay data on the KEGG background PNG.

Fixes over v2.x
---------------
* Compound nodes were painted at the vertically mirrored position (a single
  y-flip in ``_paint_cpd_nodes`` where the gene painter flipped twice and
  cancelled out).  All geometry now comes from :mod:`layout`.
* Compound radius used the full node width, drawing every metabolite at twice
  its real size.
* Circles were rasterised with a hard boolean mask over a full-image meshgrid:
  aliased edges, and an O(image area) allocation *per compound node*.  Painting
  is now supersampled within the node's own bounding box only — smooth edges
  and orders of magnitude less work.
* The colour key was drawn with ``plt.colorbar(ax=ax_img)``, which steals space
  from the image axes; and it sampled a 256-step continuous ramp while nodes
  used 10 discrete bins, so the key did not describe the figure.  Keys are now
  drawn in their own axes from the same bins as the nodes, one per data class.
* The signature was placed in data coordinates with the ``transform`` argument
  commented out, so it landed in the top-left corner of the image.

Public API
----------
  keggview_native : overlay onto the KEGG PNG and save
  paint_nodes     : paint node colours into an RGB array
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

from .color_mapping import ColorScale, draw_dual_key
from .errors import RenderError
from .layout import NodeBox, node_boxes, slice_bounds
from .utils import is_transparent, to_rgb

_SUPERSAMPLE = 4


# ---------------------------------------------------------------------------
# Painting primitives
# ---------------------------------------------------------------------------

def _color_columns(df: pl.DataFrame | None) -> list[str]:
    return [c for c in df.columns if c.endswith("_col")] if df is not None else []


def _color_lookup(col_data: pl.DataFrame | None) -> dict[str, dict]:
    if col_data is None or col_data.is_empty():
        return {}
    key = "id" if "id" in col_data.columns else col_data.columns[0]
    return {str(row[key]): row for row in col_data.iter_rows(named=True)}


def _paint_rect(
    img: np.ndarray,
    x0: float, x1: float, y0: float, y1: float,
    rgb: np.ndarray,
    preserve_dark: bool = True,
    dark_threshold: int = 160,
) -> None:
    """
    Fill an axis-aligned band, leaving dark pixels (text, borders) intact.

    KEGG glyph outlines and gene symbols are near-black; painting over them
    erases the labels the figure depends on.
    """
    h, w = img.shape[:2]
    px0, px1 = max(0, int(round(x0))), min(w, int(round(x1)))
    py0, py1 = max(0, int(round(y0))), min(h, int(round(y1)))
    if px1 <= px0 or py1 <= py0:
        return

    region = img[py0:py1, px0:px1, :3]
    if preserve_dark:
        mask = region.max(axis=2) > dark_threshold
        region[mask] = rgb
    else:
        region[:] = rgb
    img[py0:py1, px0:px1, :3] = region


def _circle_coverage(box: NodeBox, x0: float, x1: float) -> tuple[np.ndarray, tuple[int, int, int, int]]:
    """
    Anti-aliased coverage for the slice of a circle between x0 and x1.

    Supersampled inside the node's bounding box only, rather than allocating a
    full-image meshgrid per node as v2.x did.
    """
    px0, py0 = int(np.floor(box.left)), int(np.floor(box.top))
    px1, py1 = int(np.ceil(box.right)), int(np.ceil(box.bottom))
    nw, nh = max(1, px1 - px0), max(1, py1 - py0)
    s = _SUPERSAMPLE

    xs = px0 + (np.arange(nw * s) + 0.5) / s
    ys = py0 + (np.arange(nh * s) + 0.5) / s
    gx, gy = np.meshgrid(xs, ys)

    r = box.radius
    inside = ((gx - box.x) ** 2 + (gy - box.y) ** 2) <= r * r
    inside &= (gx >= x0) & (gx < x1)

    cov = inside.reshape(nh, s, nw, s).mean(axis=(1, 3))
    return cov, (px0, py0, px1, py1)


def _paint_circle_slice(img: np.ndarray, box: NodeBox,
                        x0: float, x1: float, rgb: np.ndarray) -> None:
    """Alpha-composite one wedge of a compound circle onto the image."""
    cov, (px0, py0, px1, py1) = _circle_coverage(box, x0, x1)
    h, w = img.shape[:2]
    cx0, cy0 = max(0, px0), max(0, py0)
    cx1, cy1 = min(w, px1), min(h, py1)
    if cx1 <= cx0 or cy1 <= cy0:
        return

    sub = cov[cy0 - py0:cy1 - py0, cx0 - px0:cx1 - px0][..., None]
    if sub.size == 0:
        return

    region = img[cy0:cy1, cx0:cx1, :3].astype(np.float32)
    blended = region * (1.0 - sub) + rgb.astype(np.float32) * sub
    img[cy0:cy1, cx0:cx1, :3] = np.clip(blended, 0, 255).astype(np.uint8)


def paint_nodes(
    img: np.ndarray,
    plot_data: pl.DataFrame | None,
    col_data: pl.DataFrame | None,
    node_kind: str = "gene",
) -> np.ndarray:
    """
    Paint node colours onto an H x W x 3 uint8 image.

    Multi-column data becomes vertical slices across the node, one per
    experiment, in column order.
    """
    if plot_data is None or col_data is None or plot_data.is_empty():
        return img

    ccols = _color_columns(col_data)
    if not ccols:
        return img

    lookup = _color_lookup(col_data)
    boxes = node_boxes(plot_data, default_type=node_kind)

    for box in boxes:
        row = lookup.get(box.entry_id)
        if not row:
            continue
        bands = slice_bounds(box, len(ccols))
        for (x0, x1), ccol in zip(bands, ccols):
            color = row.get(ccol)
            if is_transparent(color):
                continue
            rgb = np.array(to_rgb(color), dtype=np.uint8)
            if box.is_round:
                _paint_circle_slice(img, box, x0, x1, rgb)
            else:
                _paint_rect(img, x0, x1, box.top, box.bottom, rgb)

    return img


# ---------------------------------------------------------------------------
# keggview_native
# ---------------------------------------------------------------------------

def render_native_array(
    plot_data_gene=None, cols_gene=None,
    plot_data_cpd=None, cols_cpd=None,
    background: str | Path = "",
):
    """
    Paint the data onto the KEGG PNG and return a :class:`RasterFrame`.

    The frame is in KGML pixel space at scale 1, because KEGG's map image and
    its KGML coordinates share an origin and a unit.
    """
    from PIL import Image

    from .layout import RasterFrame

    img = np.array(Image.open(background).convert("RGB"), dtype=np.uint8)
    img = paint_nodes(img, plot_data_gene, cols_gene, "gene")
    img = paint_nodes(img, plot_data_cpd, cols_cpd, "compound")
    return RasterFrame(img, x0=0.0, y0=0.0, scale=1.0)


def keggview_native(
    plot_data_gene: pl.DataFrame | None = None,
    cols_gene: pl.DataFrame | None = None,
    plot_data_cpd: pl.DataFrame | None = None,
    cols_cpd: pl.DataFrame | None = None,
    node_data: pl.DataFrame | None = None,
    pathway_name: str = "pathway",
    kegg_dir: str | Path = ".",
    out_dir: str | Path | None = None,
    out_suffix: str = "pathview",
    gene_scale: ColorScale | None = None,
    cpd_scale: ColorScale | None = None,
    title: str | None = None,
    new_signature: bool = True,
    plot_col_key: bool = True,
    dpi: int = 200,
    background: str | Path | None = None,
    output_format: str = "png",
) -> Path:
    """
    Overlay data onto the KEGG background PNG and save the figure.

    Returns the path written.
    """
    import matplotlib
    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    kegg_dir = Path(kegg_dir)
    out_dir = Path(out_dir) if out_dir is not None else kegg_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    png_path = Path(background) if background else kegg_dir / f"{pathway_name}.png"
    if not png_path.exists():
        raise RenderError(
            f"Background image not found: {png_path}. "
            "KEGG's PNG is required for kegg_native rendering; use "
            "render_mode='vector' to draw the map from the KGML instead, "
            "which needs no background and works offline."
        )

    frame = render_native_array(plot_data_gene, cols_gene,
                                plot_data_cpd, cols_cpd, background=png_path)
    img = frame.array

    h, w = img.shape[:2]
    fig_w = min(16.0, max(8.0, w / 110))
    fig_h = fig_w * h / w

    keys = [s for s in (gene_scale, cpd_scale) if s is not None] if plot_col_key else []
    key_band = 0.085 if keys else 0.0
    title_band = 0.05 if title else 0.0

    fig = plt.figure(figsize=(fig_w, fig_h * (1 + key_band + title_band)),
                     facecolor="white")
    ax = fig.add_axes([0.01, key_band + 0.02, 0.98,
                       0.96 - key_band - title_band])
    ax.imshow(img, interpolation="antialiased", aspect="equal")
    ax.set_axis_off()

    if title:
        fig.text(0.5, 0.985, title, ha="center", va="top",
                 fontsize=12.5, fontweight="semibold", color="#1A1A1A")

    if keys:
        draw_dual_key(fig, gene_scale, cpd_scale,
                      rect=(0.12, 0.035, 0.76, 0.028), label_size=7.5)

    if new_signature:
        fig.text(0.995, 0.006, "pathview-plus", ha="right", va="bottom",
                 fontsize=6.5, color="#9AA0A6", style="italic")

    out_path = out_dir / f"{pathway_name}.{out_suffix}.{output_format}"
    fig.savefig(out_path, dpi=dpi, facecolor="white",
                bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    return out_path
