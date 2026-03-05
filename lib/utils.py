"""
utils.py
General-purpose utility functions:
  - String wrapping / fitting  (wordwrap, _strfit)
  - Numeric aggregators        (max_abs, random_pick)
"""

from __future__ import annotations

import textwrap

import numpy as np


# ---------------------------------------------------------------------------
# String utilities
# ---------------------------------------------------------------------------

def wordwrap(text: str, width: int = 20, break_word: bool = False) -> str:
    """
    Wrap *text* to *width* columns.

    When *break_word* is False (default) wrapping only occurs at whitespace.
    When True, long words are split and a backslash continuation marker is
    inserted at hard break points.
    """
    if not break_word:
        return "\n".join(textwrap.wrap(text, width))
    return _strfit(text, width)


def _strfit(s: str, width: int = 20) -> str:
    """
    Hard-wrap *s* to *width* characters per line.

    Prefers whitespace break points within ±2 characters of *width*.
    Falls back to a forced break with a trailing '\\' when none is found.
    """
    s = " ".join(s.split())   # collapse all whitespace to single spaces
    chars = list(s)
    lines: list[str] = []

    while chars:
        if len(chars) <= width + 3:
            lines.append("".join(chars))
            break

        # Prefer a whitespace break closest to the target width
        for delta in (0, 1, 2, -1, -2):
            pos = width + delta
            if 0 < pos < len(chars) and chars[pos] == " ":
                lines.append("".join(chars[:pos]))
                chars = chars[pos + 1:]
                break
        else:
            # No nearby whitespace — force a mid-word break
            lines.append("".join(chars[:width]) + "\\")
            chars = chars[width:]

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Numeric aggregators
# ---------------------------------------------------------------------------

def max_abs(values: np.ndarray) -> float:
    """Return the element with the largest absolute value, ignoring NaNs."""
    clean = values[~np.isnan(values)]
    if clean.size == 0:
        return float("nan")
    return float(clean[np.argmax(np.abs(clean))])


def random_pick(values: np.ndarray) -> float:
    """Return a randomly chosen element, ignoring NaNs.  NaN if empty."""
    clean = values[~np.isnan(values)]
    if clean.size == 0:
        return float("nan")
    return float(np.random.choice(clean))
