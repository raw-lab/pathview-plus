"""
color_mapping.py
Colour-scale utilities:
  - make_colormap  : build a three-point diverging LinearSegmentedColormap
  - node_color     : map numeric node values → hex colour strings
  - draw_color_key : render a colour-bar legend onto a Matplotlib Axes
"""

from __future__ import annotations

from typing import Callable, Optional

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.colors import LinearSegmentedColormap, Normalize

from .constants import SumMethod


# ---------------------------------------------------------------------------
# Colormap construction
# ---------------------------------------------------------------------------

def make_colormap(
    low: str = "green",
    mid: str = "gray",
    high: str = "red",
    n: int = 256,
) -> LinearSegmentedColormap:
    """
    Build a three-point diverging colour map: low → mid → high.

    Parameters
    ----------
    low, mid, high: Matplotlib colour strings or hex codes.
    n:              Number of discrete colour levels.
    """
    return LinearSegmentedColormap.from_list("pv_cmap", [low, mid, high], N=n)


# ---------------------------------------------------------------------------
# Node colour mapping
# ---------------------------------------------------------------------------

def node_color(
    plot_data: pl.DataFrame,
    limit: float | tuple[float, float] = 1.0,
    bins: int = 10,
    both_dirs: bool = True,
    discrete: bool = False,
    low: str = "green",
    mid: str = "gray",
    high: str = "red",
    na_col: str = "transparent",
    trans_fun: Optional[Callable[[np.ndarray], np.ndarray]] = None,
) -> pl.DataFrame:
    """
    Convert numeric node values to hex colour strings.

    Parameters
    ----------
    plot_data: DataFrame with an 'id' column and one or more numeric value
               columns.  Each numeric column produces a paired '*_col' column.
    limit:     Scalar (symmetric ±limit) or (vmin, vmax) tuple.
    bins:      Number of colour bins.
    both_dirs: When True and *limit* is scalar, use ±limit range.
    discrete:  Reserved for future discrete-colour support.
    low/mid/high: Colour endpoints.
    na_col:    Colour string for NaN values (default "transparent").
    trans_fun: Optional transformation applied to values before colouring.

    Returns a DataFrame with 'id' and one '*_col' column per input value column.
    """
    vmin, vmax = _resolve_limits(limit, both_dirs)
    cmap = make_colormap(low, mid, high, n=bins)
    norm = Normalize(vmin=vmin, vmax=vmax, clip=True)

    value_cols = [c for c in plot_data.columns if c != "id"]
    result: dict[str, list] = {"id": plot_data["id"].to_list()}

    #TODO: Quick fix, please verify for accuracy
    vec = np.vectorize(lambda x: int(x.lstrip("#"), 16) if isinstance(x, str) else x)
    for col in value_cols:
        vals = vec(plot_data[col].to_numpy()).astype(float)
        if trans_fun is not None:
            vals = trans_fun(vals)
        result[f"{col}_col"] = [_value_to_hex(v, cmap, norm, na_col) for v in vals]

    return pl.DataFrame(result)


def _resolve_limits(
    limit: float | tuple[float, float],
    both_dirs: bool,
) -> tuple[float, float]:
    """Convert a scalar or tuple limit into (vmin, vmax)."""
    if isinstance(limit, (int, float)):
        return (-abs(limit), abs(limit)) if both_dirs else (0.0, float(limit))
    return float(limit[0]), float(limit[1])


def _value_to_hex(
    v: float,
    cmap: LinearSegmentedColormap,
    norm: Normalize,
    na_col: str,
) -> str:
    """Map a single float to a hex colour string, returning *na_col* for NaN."""
    if np.isnan(v):
        return na_col
    r, g, b, _ = cmap(norm(v))
    return "#{:02X}{:02X}{:02X}".format(int(r * 255), int(g * 255), int(b * 255))


# ---------------------------------------------------------------------------
# Colour-key legend
# ---------------------------------------------------------------------------

def draw_color_key(
    ax: plt.Axes,
    limit: float | tuple[float, float] = 1.0,
    bins: int = 10,
    both_dirs: bool = True,
    discrete: bool = False,
    low: str = "green",
    mid: str = "gray",
    high: str = "red",
    label_size: float = 8,
) -> None:
    """
    Draw a horizontal colour-bar legend as a Matplotlib colorbar.

    Parameters
    ----------
    ax:         Axes on which to anchor the colorbar.
    limit:      Colour scale limits (scalar or tuple).
    bins:       Number of colour bins.
    both_dirs:  Whether to show negative values.
    label_size: Font size for tick labels.
    """
    vmin, vmax = _resolve_limits(limit, both_dirs)
    cmap = make_colormap(low, mid, high)
    norm = Normalize(vmin=vmin, vmax=vmax)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = plt.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.046, pad=0.04)
    ticks = [vmin, (vmin + vmax) / 2.0, vmax]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{t:.2g}" for t in ticks])
    cbar.ax.tick_params(labelsize=label_size)
