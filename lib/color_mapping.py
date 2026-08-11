"""
color_mapping.py
Value -> colour mapping, with independent scales for transcript and
metabolite data.

Design
------
RNA-seq log2 fold-changes are signed and symmetric about zero: a node must
read as *down* or *up* at a glance, and zero must be visually neutral.
Metabolite abundances live on a different scale entirely and must not be
confounded with the transcript scale — a metabolite at +2 and a transcript at
+2 are not the same statement.  pathview-plus therefore keeps two fully
independent :class:`ColorScale` objects on one map and draws a colour key for
each.

Parity with R pathview
----------------------
``node_color`` reproduces R's binning exactly:

  * scalar ``limit`` with ``both_dirs=True``  -> ``(-|limit|, +|limit|)``
  * scalar ``limit`` with ``both_dirs=False`` -> ``(0, limit)``
  * ``bins`` discrete colours from ``colorpanel2(bins, low, mid, high)``
  * cut points ``linspace(lo, hi, bins + 1)``, left-closed / right-open,
    with the top bin closed (R's ``cut(..., right=FALSE,
    include.lowest=TRUE)``)
  * values outside the limits are clamped, not dropped
  * ``discrete=True`` is honoured only when the limits are integral and
    ``(hi - lo) %% bins == 0`` — otherwise a note is issued and the data is
    treated continuously, exactly as R does.

v2.x accepted ``discrete`` and silently ignored it, used a continuous
Matplotlib ``Normalize`` instead of binning (so colours never matched the
key), and reused the same scale for genes and compounds.

Public API
----------
  ColorScale       : a configured, reusable scale
  make_colormap    : R colorpanel2-compatible discrete colormap
  node_color       : DataFrame of values -> DataFrame of hex colours
  draw_color_key   : render one colour key
  draw_dual_key    : render gene + compound keys side by side
  list_palettes    : available named palettes
"""

from __future__ import annotations

import warnings
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field

import numpy as np
import polars as pl

from .constants import (
    DEFAULT_CPD_PALETTE,
    DEFAULT_GENE_PALETTE,
    NA_COLOR,
    NODE_META_COLS,
    PALETTES,
)
from .utils import to_hex, to_rgb

# ---------------------------------------------------------------------------
# Palettes
# ---------------------------------------------------------------------------

def list_palettes() -> dict[str, tuple[str, str, str]]:
    """Named (low, mid, high) presets available to ``palette=``."""
    return dict(PALETTES)


def resolve_palette(
    palette: str | None,
    low: str | None,
    mid: str | None,
    high: str | None,
    fallback: str = DEFAULT_GENE_PALETTE,
) -> tuple[str, str, str]:
    """
    Resolve a palette name plus explicit overrides into (low, mid, high).

    Explicit ``low``/``mid``/``high`` always win over the named palette, so a
    user can nudge one end without redefining the whole scale.
    """
    if palette:
        key = str(palette).lower()
        if key not in PALETTES:
            raise ValueError(
                f"Unknown palette {palette!r}. Available: {', '.join(sorted(PALETTES))}"
            )
        base = PALETTES[key]
    else:
        base = PALETTES[fallback]
    return (low or base[0], mid or base[1], high or base[2])


# ---------------------------------------------------------------------------
# colorpanel2 — R-compatible discrete ramp
# ---------------------------------------------------------------------------

def colorpanel2(n: int, low: object, mid: object = None, high: object = None) -> list[str]:
    """
    Discrete colour ramp matching R pathview's ``colorpanel2``.

    With ``mid`` supplied the ramp is two-legged (low -> mid -> high) and, for
    odd *n*, R drops the duplicated midpoint sample rather than interpolating
    across it; that behaviour is reproduced here so bin colours are identical
    between the two implementations.
    """
    n = max(1, int(n))
    if mid is None and high is None:
        return [to_hex(low)] * n

    if high is None:                       # two-point ramp
        lo, hi = np.array(to_rgb(low), float), np.array(to_rgb(mid), float)
        if n == 1:
            return [to_hex(tuple(lo.astype(int)))]
        ramp = np.linspace(lo, hi, n)
        return [f"#{int(round(r)):02X}{int(round(g)):02X}{int(round(b)):02X}" for r, g, b in ramp]

    isodd = n % 2 == 1
    n_even = n + 1 if isodd else n
    lo = np.array(to_rgb(low), float)
    md = np.array(to_rgb(mid), float)
    hi = np.array(to_rgb(high), float)

    lower = n_even // 2
    upper = n_even - lower
    first = np.linspace(lo, md, lower) if lower > 1 else np.array([lo])
    second = np.linspace(md, hi, upper) if upper > 1 else np.array([hi])
    ramp = np.vstack([first, second])

    if isodd:                              # R deletes element (lower + 1), 1-based
        ramp = np.delete(ramp, lower, axis=0)

    return [
        f"#{int(round(r)):02X}{int(round(g)):02X}{int(round(b)):02X}"
        for r, g, b in ramp[:n]
    ]


def make_colormap(low="green", mid="gray", high="red", n: int = 256):
    """Build a Matplotlib ``LinearSegmentedColormap`` from three anchors."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(
        "pathview_plus", [to_hex(low, False), to_hex(mid, False), to_hex(high, False)],
        N=max(2, int(n)),
    )


# ---------------------------------------------------------------------------
# ColorScale
# ---------------------------------------------------------------------------

@dataclass
class ColorScale:
    """
    A configured colour scale.

    One of these exists per data class on a map (transcripts, metabolites),
    which is what keeps the two from being conflated.

    Parameters
    ----------
    limit:     Scalar (symmetric) or (vmin, vmax).
    bins:      Number of discrete colour bins.
    both_dirs: Scalar limits become +/-limit when True, (0, limit) when False.
    discrete:  Treat values as integer levels when the limits permit it.
    low/mid/high, palette: Colour anchors.
    na_col:    Colour for unmapped / NaN nodes.
    trans_fun: Optional transform applied to values before binning.
    label:     Key title, e.g. "RNA-seq log2FC".
    """

    limit: float | tuple[float, float] = 1.0
    bins: int = 10
    both_dirs: bool = True
    discrete: bool = False
    low: str | None = None
    mid: str | None = None
    high: str | None = None
    palette: str | None = None
    na_col: str = NA_COLOR
    trans_fun: Callable[[np.ndarray], np.ndarray] | None = None
    label: str = ""
    _fallback: str = field(default=DEFAULT_GENE_PALETTE, repr=False)

    # -- derived ----------------------------------------------------------
    def anchors(self) -> tuple[str, str, str]:
        return resolve_palette(self.palette, self.low, self.mid, self.high, self._fallback)

    def bounds(self) -> tuple[float, float]:
        """Resolve ``limit`` into (vmin, vmax), honouring ``both_dirs``."""
        lim = self.limit
        if isinstance(lim, (int, float, np.floating, np.integer)):
            v = float(lim)
            if not np.isfinite(v) or v == 0:
                v = 1.0
            return (-abs(v), abs(v)) if self.both_dirs else (0.0, abs(v))
        lo, hi = float(lim[0]), float(lim[1])
        if hi < lo:
            lo, hi = hi, lo
        if lo == hi:
            hi = lo + 1.0
        return lo, hi

    def effective_bins(self) -> tuple[int, tuple[float, float], bool]:
        """
        Resolve (bins, bounds, discrete_applied).

        R only honours ``discrete`` when the limits are integral and the range
        divides evenly into the bin count; otherwise it emits a note and falls
        back to continuous.  Same rule here, same note.
        """
        lo, hi = self.bounds()
        bins = max(1, int(self.bins))
        if not self.discrete:
            return bins, (lo, hi), False

        integral = float(lo).is_integer() and float(hi).is_integer()
        divides = ((hi - lo) % bins) == 0
        if integral and divides:
            return bins + 1, (lo, hi + 1.0), True

        warnings.warn(
            "discrete=True ignored: limits must be integers and (high - low) "
            f"must divide bins evenly (got limit=({lo}, {hi}), bins={bins}). "
            "Data treated as continuous.",
            stacklevel=3,
        )
        return bins, (lo, hi), False

    def colors(self) -> list[str]:
        """The ordered list of bin colours."""
        low, mid, high = self.anchors()
        bins, _, _ = self.effective_bins()
        if self.both_dirs:
            return colorpanel2(bins, low, mid, high)
        return colorpanel2(bins, mid, high)          # one-sided: mid -> high

    def cut_points(self) -> np.ndarray:
        bins, (lo, hi), _ = self.effective_bins()
        return np.linspace(lo, hi, bins + 1)

    # -- mapping ----------------------------------------------------------
    def map_values(self, values: Iterable[float]) -> list[str]:
        """Map an iterable of numbers to hex colours (NaN -> ``na_col``)."""
        arr = np.asarray(list(values), dtype=float)
        if arr.size == 0:
            return []
        if self.trans_fun is not None:
            arr = np.asarray(self.trans_fun(arr), dtype=float)

        bins, (lo, hi), is_disc = self.effective_bins()
        cols = self.colors()
        na = to_hex(self.na_col)

        if is_disc:
            arr = np.where(np.isnan(arr), np.nan, np.floor(arr))

        finite = ~np.isnan(arr)
        clamped = np.clip(arr, lo, hi)

        # Left-closed, right-open bins; the top bin is closed so the maximum
        # lands in the last bin instead of overflowing (R's include.lowest).
        edges = np.linspace(lo, hi, bins + 1)
        idx = np.digitize(clamped, edges[1:-1], right=False)
        idx = np.clip(idx, 0, bins - 1)

        return [cols[i] if ok else na for i, ok in zip(idx, finite)]

    def sample(self, n: int = 256) -> list[str]:
        """A finer sampling of the same ramp, for smooth colour bars."""
        low, mid, high = self.anchors()
        return colorpanel2(n, low, mid, high) if self.both_dirs else colorpanel2(n, mid, high)

    def to_dict(self) -> dict:
        lo, hi = self.bounds()
        low, mid, high = self.anchors()
        return {
            "limit": [lo, hi], "bins": int(self.bins), "both_dirs": self.both_dirs,
            "discrete": self.discrete, "low": low, "mid": mid, "high": high,
            "na_col": self.na_col, "label": self.label,
        }


def gene_scale(**kw) -> ColorScale:
    """A scale preconfigured for transcript / protein fold-changes."""
    kw.setdefault("label", "RNA-seq log2FC")
    kw.setdefault("_fallback", DEFAULT_GENE_PALETTE)
    return ColorScale(**kw)


def compound_scale(**kw) -> ColorScale:
    """A scale preconfigured for metabolite abundances."""
    kw.setdefault("label", "Metabolite log2FC")
    kw.setdefault("_fallback", DEFAULT_CPD_PALETTE)
    return ColorScale(**kw)


# ---------------------------------------------------------------------------
# node_color
# ---------------------------------------------------------------------------

def node_color(
    plot_data: pl.DataFrame | None,
    scale: ColorScale | None = None,
    *,
    id_col: str = "id",
    value_cols: Sequence[str] | None = None,
    limit: float | tuple[float, float] = 1.0,
    bins: int = 10,
    both_dirs: bool = True,
    discrete: bool = False,
    low: str | None = None,
    mid: str | None = None,
    high: str | None = None,
    palette: str | None = None,
    na_col: str = NA_COLOR,
    trans_fun: Callable[[np.ndarray], np.ndarray] | None = None,
) -> pl.DataFrame | None:
    """
    Map every numeric column of *plot_data* to a matching ``*_col`` column.

    Only genuinely numeric columns are mapped.  v2.x selected value columns by
    subtracting a metadata set that was missing ``bgcolor``, then coerced the
    resulting ``"#FFFFFF"`` strings with ``int(hex, 16)`` — 16777215, which
    clamps to the top of every scale.  That single omission painted every node
    on every map solid red regardless of the user's data.  Selection is now
    dtype-driven, so a stray string column cannot be mistaken for data.

    Returns a DataFrame with *id_col* plus one ``<value>_col`` per value
    column, or None when *plot_data* is None.
    """
    if plot_data is None or plot_data.is_empty():
        return None

    sc = scale or ColorScale(
        limit=limit, bins=bins, both_dirs=both_dirs, discrete=discrete,
        low=low, mid=mid, high=high, palette=palette, na_col=na_col,
        trans_fun=trans_fun,
    )

    if id_col not in plot_data.columns:
        raise ValueError(f"node_color: id column {id_col!r} not in {plot_data.columns}")

    if value_cols is None:
        value_cols = [
            c for c in plot_data.columns
            if c != id_col
            and c not in NODE_META_COLS
            and plot_data.schema[c].is_numeric()
        ]
    value_cols = list(value_cols)

    out: dict[str, list] = {id_col: plot_data[id_col].to_list()}
    for col in value_cols:
        out[f"{col}_col"] = sc.map_values(plot_data[col].to_list())

    return pl.DataFrame(out)


def value_columns(df: pl.DataFrame | None, id_col: str = "entry_id") -> list[str]:
    """Numeric, non-metadata columns of *df* — the columns that carry data."""
    if df is None:
        return []
    return [
        c for c in df.columns
        if c != id_col and c not in NODE_META_COLS and df.schema[c].is_numeric()
    ]


# ---------------------------------------------------------------------------
# Colour keys
# ---------------------------------------------------------------------------

def draw_color_key(
    ax,
    scale: ColorScale,
    *,
    title: str | None = None,
    n_ticks: int = 5,
    label_size: float = 8.0,
    show_bins: bool = True,
    border: str = "#333333",
) -> None:
    """
    Draw a colour key onto *ax*, which the caller has already positioned.

    v2.x called ``plt.colorbar(..., ax=ax_img)``, which steals space from the
    pathway axes and rescales the image; it also built the key from a 256-step
    continuous ramp while the nodes were painted from a 10-bin ramp, so the
    key did not describe the figure it sat next to.  Here the key is drawn
    from ``scale.colors()`` — the same list the nodes used.
    """
    lo, hi = scale.bounds()
    colors = scale.colors() if show_bins else scale.sample(256)
    n = len(colors)

    for i, c in enumerate(colors):
        ax.add_patch(_rect(i / n, 0.0, 1.0 / n, 1.0, c))

    ax.add_patch(_rect(0, 0, 1, 1, "none", edge=border, lw=0.8))

    ticks = np.linspace(lo, hi, max(2, int(n_ticks)))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks(np.linspace(0, 1, len(ticks)))
    ax.set_xticklabels([_fmt(t) for t in ticks], fontsize=label_size)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=2, pad=1.5, colors="#333333")
    for spine in ax.spines.values():
        spine.set_visible(False)

    label = title if title is not None else scale.label
    if label:
        ax.set_title(label, fontsize=label_size + 0.5, pad=3.5,
                     color="#222222", fontweight="medium")


def draw_dual_key(
    fig,
    gene_scale_: ColorScale | None,
    cpd_scale_: ColorScale | None,
    *,
    rect: tuple[float, float, float, float] = (0.10, 0.03, 0.80, 0.045),
    label_size: float = 8.0,
) -> list:
    """
    Draw the transcript and metabolite keys side by side beneath a figure.

    This is what makes a combined RNA-seq + metabolomics map readable: two
    keys, each labelled, each describing only its own data class.  Returns the
    axes created.
    """
    scales = [s for s in (gene_scale_, cpd_scale_) if s is not None]
    if not scales:
        return []

    x0, y0, w, h = rect
    axes = []
    gap = 0.06
    each = (w - gap * (len(scales) - 1)) / len(scales)
    for i, sc in enumerate(scales):
        ax = fig.add_axes([x0 + i * (each + gap), y0, each, h])
        draw_color_key(ax, sc, label_size=label_size)
        axes.append(ax)
    return axes


def _rect(x, y, w, h, face, edge="none", lw=0.0):
    from matplotlib.patches import Rectangle
    return Rectangle((x, y), w, h, facecolor=face, edgecolor=edge,
                     linewidth=lw, zorder=2 if face == "none" else 1)


def _fmt(v: float) -> str:
    if abs(v) >= 1000 or (v != 0 and abs(v) < 0.01):
        return f"{v:.1e}"
    if float(v).is_integer():
        return str(int(v))
    return f"{v:.2f}".rstrip("0").rstrip(".")
