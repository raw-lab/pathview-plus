"""
rendering.py
Pathway diagram rendering:
  - keggview_native : overlay data on a KEGG background PNG (pixel painting)
  - keggview_graph  : draw a NetworkX graph diagram styled with Seaborn
  - kegg_legend     : display a standalone KEGG diagram element legend
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns
from PIL import Image

from .color_mapping import draw_color_key, make_colormap
from .utils import wordwrap


# ---------------------------------------------------------------------------
# KEGG edge subtype reference table
# ---------------------------------------------------------------------------

_EDGE_SUBTYPES = [
    # (name,              colour,    label, style,  arrowhead)
    ("activation",        "#00CC00", "-->",  "solid",  "normal"),
    ("inhibition",        "#CC0000", "--|",  "solid",  "tee"),
    ("expression",        "#00AA00", "-->",  "dashed", "normal"),
    ("repression",        "#AA0000", "--|",  "dashed", "tee"),
    ("indirect",          "#888888", "..>",  "dotted", "normal"),
    ("binding",           "#0000CC", "---",  "solid",  "none"),
    ("compound",          "#8800AA", "---",  "solid",  "none"),
    ("phosphorylation",   "#FF6600", "+p",   "solid",  "normal"),
    ("dephosphorylation", "#FF6600", "-p",   "solid",  "normal"),
    ("ubiquitination",    "#FF00FF", "+u",   "solid",  "normal"),
    ("methylation",       "#00AAFF", "+m",   "solid",  "normal"),
    ("others/unknown",    "#888888", "?",    "solid",  "normal"),
]


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _hex_to_rgb255(hex_col: str) -> Optional[np.ndarray]:
    """
    Convert a hex colour string to a uint8 [R, G, B] array.
    Returns None for transparent / empty strings.
    """
    hex_col = hex_col.lstrip("#")
    if not hex_col or hex_col.lower() in ("transparent", "none"):
        return None
    return np.array([int(hex_col[i:i+2], 16) for i in (0, 2, 4)], dtype=np.uint8)


def _color_cols(df: pl.DataFrame) -> list[str]:
    """Return column names that end with '_col'."""
    return [c for c in df.columns if c.endswith("_col")]


# ---------------------------------------------------------------------------
# Native view  (overlay on KEGG PNG)
# ---------------------------------------------------------------------------

def _paint_gene_nodes(
    img: np.ndarray,
    plot_data: pl.DataFrame,
    col_data: pl.DataFrame,
) -> np.ndarray:
    """
    Paint gene-node rectangles onto a H×W×3 uint8 image array.

    The node width is divided evenly across multi-state colour columns so
    that each experiment gets a horizontal slice of the node box.
    """
    h = img.shape[0]
    ccols = _color_cols(col_data)
    n_states = len(ccols)
    col_lookup = {row["id"]: row for row in col_data.iter_rows(named=True)}

    for row in plot_data.iter_rows(named=True):
        cx, cy   = row["x"], row["y"]
        cy = h - cy # Flip y to convert from math to image coordinate space
        hw, hh   = row["width"] / 2, row["height"] / 2
        px_l     = max(0, int(cx - hw))
        px_r     = min(img.shape[1], int(cx + hw))
        py_t     = max(0, int(h - cy - hh))
        py_b     = min(h, int(h - cy + hh))
        x_breaks = np.linspace(px_l, px_r, n_states + 1, dtype=int)

        node_cols = col_lookup.get(row["entry_id"], {})
        for k, ccol in enumerate(ccols):
            rgb = _hex_to_rgb255(node_cols.get(ccol, ""))
            if rgb is None:
                continue
            sl_l, sl_r = int(x_breaks[k]), int(x_breaks[k + 1])
            region = img[py_t:py_b, sl_l:sl_r, :3]
            # Keep black pixels (borders / text)
            mask = region.sum(axis=2) > 0
            region[mask] = rgb
            img[py_t:py_b, sl_l:sl_r, :3] = region

    return img


def _paint_cpd_nodes(
    img: np.ndarray,
    plot_data: pl.DataFrame,
    col_data: pl.DataFrame,
) -> np.ndarray:
    """
    Paint compound-node ellipses onto a H×W×3 uint8 image array.

    Multi-state colours are applied as vertical slices through the circle.
    """
    h, w_img = img.shape[:2]
    ccols = _color_cols(col_data)
    n_states = len(ccols)
    yy, xx = np.mgrid[0:h, 0:w_img]
    col_lookup = {row["id"]: row for row in col_data.iter_rows(named=True)}

    for row in plot_data.iter_rows(named=True):
        cx, cy, r = row["x"], row["y"], row["width"]
        cy = h - cy # Flip y to convert from math to image coordinate space
        dist_sq   = (xx - cx) ** 2 + (yy - cy) ** 2
        inside    = dist_sq < r ** 2
        border    = (dist_sq >= (r - 2) ** 2) & inside
        x_breaks  = np.linspace(cx - r, cx + r, n_states + 1)

        node_cols = col_lookup.get(row["entry_id"], {})
        for k, ccol in enumerate(ccols):
            rgb = _hex_to_rgb255(node_cols.get(ccol, ""))
            if rgb is None:
                continue
            mask = inside & (xx >= x_breaks[k]) & (xx < x_breaks[k + 1])
            img[mask, :3] = rgb

        img[border, :3] = 0   # restore black border

    return img


def keggview_native(
    plot_data_gene: Optional[pl.DataFrame],
    cols_gene: Optional[pl.DataFrame],
    plot_data_cpd: Optional[pl.DataFrame],
    cols_cpd: Optional[pl.DataFrame],
    node_data: pl.DataFrame,
    pathway_name: str,
    kegg_dir: Path = Path("."),
    out_suffix: str = "pathview",
    limit: dict | None = None,
    bins: dict | None = None,
    both_dirs: dict | None = None,
    discrete: dict | None = None,
    low: dict | None = None,
    mid: dict | None = None,
    high: dict | None = None,
    new_signature: bool = True,
    plot_col_key: bool = True,
    dpi: int = 150,
) -> None:
    """
    Render expression data overlaid on the KEGG pathway PNG background.

    Reads ``<kegg_dir>/<pathway_name>.png``, paints gene and compound nodes
    with the supplied colour data, and writes
    ``<kegg_dir>/<pathway_name>.<out_suffix>.png``.
    """
    if limit     is None: limit     = {"gene": 1, "cpd": 1}
    if bins      is None: bins      = {"gene": 10, "cpd": 10}
    if both_dirs is None: both_dirs = {"gene": True, "cpd": True}
    if discrete  is None: discrete  = {"gene": False, "cpd": False}
    if low       is None: low       = {"gene": "green", "cpd": "blue"}
    if mid       is None: mid       = {"gene": "gray",  "cpd": "gray"}
    if high      is None: high      = {"gene": "red",   "cpd": "yellow"}

    png_path = Path(kegg_dir) / f"{pathway_name}.png"
    if not png_path.exists():
        raise FileNotFoundError(f"Background PNG not found: {png_path}")

    img = np.array(Image.open(png_path).convert("RGB"), dtype=np.uint8)

    if plot_data_gene is not None and cols_gene is not None:
        img = _paint_gene_nodes(img, plot_data_gene, cols_gene)
    if plot_data_cpd is not None and cols_cpd is not None:
        img = _paint_cpd_nodes(img, plot_data_cpd, cols_cpd)

    h, w = img.shape[:2]
    key_height = 0.6 if plot_col_key else 0.0
    fig, axes = plt.subplots(
        nrows=2 if plot_col_key else 1,
        figsize=(w / dpi, h / dpi + key_height),
        gridspec_kw={"height_ratios": [h, int(dpi * key_height)]} if plot_col_key else None,
    )
    ax_img = axes[0] if plot_col_key else axes

    ax_img.imshow(img, aspect="auto")
    ax_img.axis("off")

    if new_signature:
        ax_img.text(
            0.02, 0.02, "Rendered by pathview.py",
            #transform=ax_img.transAxes, #TODO: This lable looks better on top
            fontsize=6, color="black", fontweight="bold", va="bottom",
        )

    if plot_col_key and plot_data_gene is not None:
        draw_color_key(
            ax_img,
            limit=limit["gene"], bins=bins["gene"],
            both_dirs=both_dirs["gene"], discrete=discrete["gene"],
            low=low["gene"], mid=mid["gene"], high=high["gene"],
        )
        if plot_col_key:
            axes[1].set_visible(False)

    out_path = Path(kegg_dir) / f"{pathway_name}.{out_suffix}.png"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Info: Written → {out_path}")


# ---------------------------------------------------------------------------
# Graph view  (NetworkX / Seaborn)
# ---------------------------------------------------------------------------

def keggview_graph(
    plot_data_gene: Optional[pl.DataFrame],
    cols_gene: Optional[pl.DataFrame],
    plot_data_cpd: Optional[pl.DataFrame],
    cols_cpd: Optional[pl.DataFrame],
    node_data: pl.DataFrame,
    pathway_name: str,
    out_suffix: str = "pathview",
    kegg_dir: Path = Path("."),
    cex: float = 0.7,
    limit: dict | None = None,
    bins: dict | None = None,
    both_dirs: dict | None = None,
    low: dict | None = None,
    mid: dict | None = None,
    high: dict | None = None,
    new_signature: bool = True,
    plot_col_key: bool = True,
) -> None:
    """
    Render pathway as a NetworkX directed graph with Seaborn styling.

    Nodes are positioned using the KGML (x, y) coordinates.  Saves a PDF to
    ``<kegg_dir>/<pathway_name>.<out_suffix>.pdf``.
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx is required for graph view: pip install networkx")

    if limit     is None: limit     = {"gene": 1, "cpd": 1}
    if bins      is None: bins      = {"gene": 10, "cpd": 10}
    if both_dirs is None: both_dirs = {"gene": True, "cpd": True}
    if low       is None: low       = {"gene": "green", "cpd": "blue"}
    if mid       is None: mid       = {"gene": "gray",  "cpd": "gray"}
    if high      is None: high      = {"gene": "red",   "cpd": "yellow"}

    # Build colour lookup from both gene and compound colour DataFrames
    color_lookup: dict[str, str] = {}
    for col_df in (cols_gene, cols_cpd):
        if col_df is not None:
            first_col = next((c for c in col_df.columns if c.endswith("_col")), None)
            if first_col:
                for row in col_df.iter_rows(named=True):
                    color_lookup.setdefault(row["id"], row[first_col])

    # Build directed graph from node_data
    G = nx.DiGraph()
    for row in node_data.iter_rows(named=True):
        G.add_node(row["entry_id"], **row)

    pos = {
        row["entry_id"]: (
            row["x"] if row["x"] is not None else 0.0,
            -(row["y"] if row["y"] is not None else 0.0),
        )
        for row in node_data.iter_rows(named=True)
    }
    node_colors = [color_lookup.get(n, "#CCCCCC") for n in G.nodes]
    node_labels = {
        row["entry_id"]: wordwrap(row.get("label", ""), width=12)
        for row in node_data.iter_rows(named=True)
    }

    with sns.axes_style("white"):
        fig, ax = plt.subplots(figsize=(14, 10))
        ax.set_title(pathway_name, fontsize=12, fontweight="bold")

        #TODO: Temporary fix, update prior steps that use transparent instead of none
        node_colors = ['none' if x=='transparent' else x for x in node_colors]
        nx.draw_networkx(
            G,
            pos=pos,
            ax=ax,
            labels=node_labels,
            node_color=node_colors,
            node_size=800,
            font_size=cex * 10,
            arrows=True,
            arrowsize=12,
            edge_color="#555555",
        )

        if plot_col_key:
            draw_color_key(
                ax,
                limit=limit["gene"], bins=bins["gene"],
                both_dirs=both_dirs["gene"],
                low=low["gene"], mid=mid["gene"], high=high["gene"],
            )

        if new_signature:
            ax.text(
                0.01, 0.01, "Rendered by pathview.py",
                transform=ax.transAxes, fontsize=7, va="bottom",
            )
        ax.axis("off")

    out_path = Path(kegg_dir) / f"{pathway_name}.{out_suffix}.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Info: Written → {out_path}")


# ---------------------------------------------------------------------------
# KEGG legend
# ---------------------------------------------------------------------------

def kegg_legend(
    legend_type: str = "both",
) -> None:
    """
    Display a standalone reference legend for KEGG pathway elements.

    Parameters
    ----------
    legend_type: One of "both", "edge", or "node".
    """
    if legend_type not in ("both", "edge", "node"):
        warnings.warn(f"legend_type must be 'both', 'edge', or 'node'; got '{legend_type}'.")
        return

    n = len(_EDGE_SUBTYPES)
    with sns.axes_style("white"):
        fig, ax = plt.subplots(figsize=(9, 7))
        ax.set_xlim(-0.2, 4.5)
        ax.set_ylim(-0.5, n + 1.5)
        ax.axis("off")
        ax.set_title("KEGG Diagram Legend", fontweight="bold", fontsize=12)

        _line_styles = {"solid": "-", "dashed": "--", "dotted": ":"}

        if legend_type in ("both", "edge"):
            ax.text(0.9, n + 1.0, "Edge Types", fontsize=10, fontweight="bold", ha="right")
            for i, (name, col, label, style, arrow) in enumerate(_EDGE_SUBTYPES):
                y = n - i - 0.5
                ax.text(0.85, y, name, ha="right", va="center", fontsize=8)
                ax.annotate(
                    "",
                    xy=(1.8, y), xytext=(1.0, y),
                    arrowprops=dict(
                        arrowstyle="->" if arrow == "normal" else "-|>",
                        color=col,
                        linestyle=_line_styles.get(style, "-"),
                        lw=1.5,
                    ),
                )
                ax.text(1.4, y + 0.22, label, color=col, fontsize=7, ha="center")

        if legend_type in ("both", "node"):
            x_off = 2.5 if legend_type == "both" else 0.5
            ax.text(x_off + 1.2, n + 1.0, "Node Types", fontsize=10, fontweight="bold", ha="right")
            node_specs = [
                ("gene / protein / enzyme", "rectangle"),
                ("compound / metabolite",   "ellipse"),
                ("pathway link",            "text"),
            ]
            for i, (label, shape) in enumerate(node_specs):
                y = n - i * 3.5 - 0.5
                ax.text(x_off + 1.1, y, label, ha="right", va="center", fontsize=8)
                xc = x_off + 1.5
                if shape == "ellipse":
                    ax.add_patch(mpatches.Ellipse(
                        (xc, y), 0.45, 0.28, color="#DDDDDD", ec="black", lw=1,
                    ))
                elif shape == "rectangle":
                    ax.add_patch(mpatches.FancyBboxPatch(
                        (xc - 0.22, y - 0.14), 0.44, 0.28,
                        boxstyle="square", color="#DDDDDD", ec="black", lw=1,
                    ))
                else:
                    ax.text(xc, y, "Pathway Name", ha="center", va="center",
                            fontsize=8, style="italic")

    plt.tight_layout()
    plt.show()
