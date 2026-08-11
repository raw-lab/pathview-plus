"""
utils.py
String, colour and numeric helpers shared across the package.
"""

from __future__ import annotations

import textwrap
from collections.abc import Iterable, Sequence

import numpy as np

# ---------------------------------------------------------------------------
# Colour parsing
#
# v2.x used ``hex.lstrip("#")`` then ``int(..., 16)``, which raised
# ValueError("invalid literal for int() with base 16: 're'") on any named
# colour — including the default ``color="red"`` of highlight_nodes(), so the
# highlighting API crashed on its own defaults.  Everything now goes through
# to_rgb / to_hex, which accept names, #rgb, #rrggbb, #rrggbbaa and tuples.
# ---------------------------------------------------------------------------

TRANSPARENT_TOKENS = frozenset({"transparent", "none", "na", "", "#00000000"})


def is_transparent(color: object) -> bool:
    """True when *color* denotes 'do not paint'."""
    if color is None:
        return True
    if isinstance(color, str):
        return color.strip().lower() in TRANSPARENT_TOKENS
    return False


def to_rgb(color: object, default: tuple[int, int, int] = (255, 255, 255)) -> tuple[int, int, int]:
    """
    Convert any Matplotlib-compatible colour to an (R, G, B) 0-255 tuple.

    Accepts named colours ("red"), hex ("#F00", "#FF0000", "#FF0000CC"),
    float tuples, and returns *default* for transparent/unparseable input.
    """
    if is_transparent(color):
        return default
    if isinstance(color, (tuple, list, np.ndarray)):
        vals = list(color)[:3]
        if vals and max(float(v) for v in vals) <= 1.0:
            return tuple(int(round(float(v) * 255)) for v in vals)      # type: ignore[return-value]
        return tuple(int(v) for v in vals)                              # type: ignore[return-value]

    s = str(color).strip()
    if s.startswith("#"):
        h = s[1:]
        if len(h) == 3:
            h = "".join(c * 2 for c in h)
        if len(h) >= 6:
            try:
                return (int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16))
            except ValueError:
                return default
    try:
        from matplotlib.colors import to_rgb as mpl_to_rgb
        r, g, b = mpl_to_rgb(s)
        return (int(round(r * 255)), int(round(g * 255)), int(round(b * 255)))
    except Exception:
        return default


def to_hex(color: object, keep_transparent: bool = True) -> str:
    """Normalise any colour spec to '#RRGGBB' (or 'transparent')."""
    if is_transparent(color):
        return "transparent" if keep_transparent else "#FFFFFF"
    r, g, b = to_rgb(color)
    return f"#{r:02X}{g:02X}{b:02X}"


def contrast_text_color(background: object, light: str = "#FFFFFF", dark: str = "#111111") -> str:
    """
    Pick a legible text colour for *background* using WCAG relative luminance.

    Node labels in v2.x were always black, which is unreadable on a saturated
    down-regulated node.  This keeps labels legible at both ends of the scale.
    """
    r, g, b = to_rgb(background, default=(255, 255, 255))

    def chan(c: float) -> float:
        c = c / 255.0
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

    lum = 0.2126 * chan(r) + 0.7152 * chan(g) + 0.0722 * chan(b)
    return dark if lum > 0.45 else light


def blend(fg: object, bg: object, alpha: float) -> str:
    """Alpha-composite *fg* over *bg* and return a hex string."""
    fr, fg_, fb = to_rgb(fg)
    br, bg_, bb = to_rgb(bg)
    a = min(1.0, max(0.0, float(alpha)))
    return f"#{int(round(fr * a + br * (1 - a))):02X}{int(round(fg_ * a + bg_ * (1 - a))):02X}{int(round(fb * a + bb * (1 - a))):02X}"


# ---------------------------------------------------------------------------
# String utilities
# ---------------------------------------------------------------------------

def wordwrap(text: object, width: int = 20, break_word: bool = False) -> str:
    """Wrap *text* to *width* columns, breaking inside words only if asked."""
    s = "" if text is None else str(text)
    if not s:
        return ""
    if not break_word:
        wrapped = textwrap.wrap(s, width)
        return "\n".join(wrapped) if wrapped else s
    return strfit(s, width)


def strfit(s: object, width: int = 20) -> str:
    """
    Hard-wrap *s* to *width* characters, preferring nearby whitespace.

    Mirrors R pathview's ``strfit``.  Guarded against the zero-width and
    empty-input infinite loops the v2.x version could enter.
    """
    s = " ".join(str(s or "").split())
    if not s:
        return ""
    width = max(1, int(width))
    chars = list(s)
    lines: list[str] = []

    while chars:
        if len(chars) <= width + 3:
            lines.append("".join(chars))
            break
        for delta in (0, 1, 2, -1, -2):
            pos = width + delta
            if 0 < pos < len(chars) and chars[pos] == " ":
                lines.append("".join(chars[:pos]))
                chars = chars[pos + 1:]
                break
        else:
            lines.append("".join(chars[:width]) + "\\")
            chars = chars[width:]
    return "\n".join(lines)


def short_label(label: object, node_type: str = "gene") -> str:
    """
    Reduce a KEGG multi-name label to its first synonym.

    KGML stores labels as "CDK4, CMM3, PSK-J3..."; R pathview's
    ``node.info(short.name=TRUE)`` keeps only the leading synonym and strips
    the ellipsis.  v2.x kept the whole string, so nodes rendered with
    overflowing multi-gene labels.
    """
    s = str(label or "").strip()
    if not s:
        return ""
    if node_type == "map":
        return s.replace("...", "")
    first = s.split(", ")[0]
    return first.replace("...", "").strip()


def escape_xml(text: object) -> str:
    """Escape the five XML predefined entities."""
    return (str(text or "")
            .replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;").replace("'", "&apos;"))


def strip_kegg_prefix(name: object) -> str:
    """'hsa:1017' -> '1017'; 'cpd:C00022' -> 'C00022'; passthrough otherwise."""
    s = str(name or "")
    return s.split(":", 1)[1] if ":" in s else s


# ---------------------------------------------------------------------------
# Numeric aggregators
# ---------------------------------------------------------------------------

def _clean(values: Iterable[float]) -> np.ndarray:
    arr = np.asarray(list(values), dtype=float)
    return arr[~np.isnan(arr)]


def max_abs(values: Iterable[float]) -> float:
    """Element with the largest absolute value, NaN-safe (signed result)."""
    clean = _clean(values)
    if clean.size == 0:
        return float("nan")
    return float(clean[np.argmax(np.abs(clean))])


def random_pick(values: Iterable[float], rng: np.random.Generator | None = None) -> float:
    """
    Randomly choose one element, NaN-safe and reproducible.

    Values are sorted before sampling.  A seeded generator alone is not
    enough: ``group_by`` gives no guarantee about within-group row order, so
    the same seed applied to differently ordered rows picked different
    elements.  Sorting makes the choice depend only on the *set* of values.
    v2.x used the global ``np.random.choice`` and could not be seeded at all.
    """
    clean = np.sort(_clean(values))
    if clean.size == 0:
        return float("nan")
    gen = rng if rng is not None else np.random.default_rng()
    return float(clean[int(gen.integers(0, clean.size))])


def nan_safe(fn, values: Sequence[float]) -> float:
    """Apply *fn* to the non-NaN subset of *values*; NaN when nothing remains."""
    clean = _clean(values)
    return float("nan") if clean.size == 0 else float(fn(clean))
