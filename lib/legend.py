"""
legend.py
Standalone legends for KEGG and SBGN diagram elements.

Fixes over v2.x
---------------
``kegg_legend`` ended in ``plt.show()`` with no save path, so on any headless
machine — every cluster, every CI job — it drew nothing and returned None.
The node-type block also placed rows at ``y = n - i*3.5 - 0.5``, which falls
below the axis limits after the second entry, so entries were clipped away.
Both legends now return the figure and can save to a file.

Public API
----------
  kegg_legend  : KEGG edge subtypes and node shapes
  sbgn_legend  : SBGN glyph and arc classes
  edge_subtypes: the bundled KEGG edge-subtype table
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import polars as pl

from .bundled import read_bundled_tsv
from .sbgn_parser import SBGN_ARC_CLASSES, SBGN_GLYPH_CLASSES
from .vector_rendering import EDGE_STYLE, THEMES

_DATA = Path(__file__).parent / "data"


@lru_cache(maxsize=1)
def edge_subtypes() -> pl.DataFrame:
    """KEGG's own edge-subtype table (name, glyph, explanation, styling)."""
    return read_bundled_tsv("edge_subtypes.tsv.gz")


def kegg_legend(
    legend_type: str = "both",
    out_path: str | Path | None = None,
    theme: str = "publication",
    dpi: int = 200,
    show: bool = False,
):
    """
    Draw the KEGG element legend.

    Returns the Matplotlib figure; writes a file when *out_path* is given.
    """
    import matplotlib
    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, FancyBboxPatch

    if legend_type not in ("both", "edge", "node"):
        raise ValueError("legend_type must be 'both', 'edge' or 'node'.")

    th = THEMES.get(theme, THEMES["publication"])
    tbl = edge_subtypes().filter(pl.col("value") != "attribute")
    rows = tbl.to_dicts()
    n = len(rows)

    two = legend_type == "both"
    fig, axes = plt.subplots(
        1, 2 if two else 1,
        figsize=(11, 0.36 * n + 1.6) if two else (6, 0.36 * n + 1.6),
        facecolor=th["bg"], gridspec_kw={"width_ratios": [1.35, 1]} if two else None,
    )
    ax_e, ax_n = (axes[0], axes[1]) if two else (axes, axes)

    if legend_type in ("both", "edge"):
        ax_e.set_xlim(0, 10)
        ax_e.set_ylim(-1, n + 1)
        ax_e.axis("off")
        ax_e.set_title("Edge subtypes", fontsize=11, fontweight="semibold",
                       color=th["title"], loc="left")
        dash = {"-": (None, None), "--": (4, 2), ":": (1, 2)}
        for i, r in enumerate(rows):
            y = n - i - 0.5
            name = r["name"]
            st = EDGE_STYLE.get(name.lower(), {"color": r.get("color") or "#888",
                                               "style": "-", "head": "arrow"})
            ax_e.text(3.7, y, name, ha="right", va="center", fontsize=8,
                      color=th["text"])
            ln, = ax_e.plot([4.1, 6.1], [y, y], color=st["color"], lw=1.5,
                            solid_capstyle="round")
            d = dash.get(st["style"], (None, None))
            if d[0]:
                ln.set_dashes(d)
            if st["head"] == "arrow":
                ax_e.annotate("", xy=(6.5, y), xytext=(6.1, y),
                              arrowprops=dict(arrowstyle="-|>", color=st["color"], lw=1.4))
            elif st["head"] == "bar":
                ax_e.plot([6.4, 6.4], [y - 0.22, y + 0.22], color=st["color"], lw=1.6)
            ax_e.text(6.9, y, r.get("value") or "", fontsize=7.5,
                      color=st["color"], va="center", family="monospace")
            ax_e.text(7.6, y, (r.get("Explanation") or "")[:44], fontsize=6.4,
                      color=th["muted"], va="center")

    if legend_type in ("both", "node"):
        specs = [
            ("gene / protein", "rect", "#BFFFBF"),
            ("enzyme (ortholog)", "rect", "#CCFFCC"),
            ("compound / metabolite", "circle", "#FFFFFF"),
            ("group / complex", "rect", "#E8EAF6"),
            ("pathway link", "round", "#ECEFF1"),
            ("no data mapped", "rect", th["unmapped"]),
        ]
        m = len(specs)
        ax_n.set_xlim(0, 10)
        ax_n.set_ylim(-1, m + 1)
        ax_n.axis("off")
        ax_n.set_title("Node types", fontsize=11, fontweight="semibold",
                       color=th["title"], loc="left")
        for i, (label, shape, fill) in enumerate(specs):
            y = m - i - 0.5                       # fits inside ylim, unlike v2.x
            ax_n.text(4.6, y, label, ha="right", va="center", fontsize=8,
                      color=th["text"])
            if shape == "circle":
                ax_n.add_patch(Circle((5.6, y), 0.26, facecolor=fill,
                                      edgecolor=th["border"], lw=0.8))
            else:
                style = "round,pad=0.06" if shape == "round" else "square,pad=0.02"
                ax_n.add_patch(FancyBboxPatch((5.15, y - 0.2), 0.95, 0.4,
                                              boxstyle=style, facecolor=fill,
                                              edgecolor=th["border"], lw=0.8))
        ax_n.text(0.1, -0.6, "Multi-condition data splits a node into vertical "
                             "slices, one per condition.",
                  fontsize=6.6, color=th["muted"])

    fig.suptitle("KEGG diagram legend", fontsize=12.5, fontweight="semibold",
                 color=th["title"])
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=dpi, facecolor=th["bg"], bbox_inches="tight")
    if show:                                                 # pragma: no cover
        plt.show()
    return fig


def sbgn_legend(out_path: str | Path | None = None, theme: str = "publication",
                dpi: int = 200, show: bool = False):
    """Draw the SBGN glyph/arc-class legend.  Returns the figure."""
    import matplotlib
    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    th = THEMES.get(theme, THEMES["publication"])
    glyphs = list(SBGN_GLYPH_CLASSES.items())
    arcs = list(SBGN_ARC_CLASSES.items())
    n = max(len(glyphs), len(arcs))

    fig, axes = plt.subplots(1, 2, figsize=(12, 0.32 * n + 1.4), facecolor=th["bg"])
    for ax, items, title in ((axes[0], glyphs, "SBGN glyph classes"),
                             (axes[1], arcs, "SBGN arc classes")):
        ax.set_xlim(0, 10)
        ax.set_ylim(-1, n + 1)
        ax.axis("off")
        ax.set_title(title, fontsize=11, fontweight="semibold",
                     color=th["title"], loc="left")
        for i, (k, v) in enumerate(items):
            y = n - i - 0.5
            ax.text(0.2, y, k, fontsize=8, color=th["text"], va="center",
                    fontweight="medium")
            ax.text(4.3, y, v, fontsize=6.8, color=th["muted"], va="center")

    fig.suptitle("SBGN reference", fontsize=12.5, fontweight="semibold",
                 color=th["title"])
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=dpi, facecolor=th["bg"], bbox_inches="tight")
    if show:                                                 # pragma: no cover
        plt.show()
    return fig
