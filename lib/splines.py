"""
splines.py
Curve construction for pathway edges.

Fixes over v2.x
---------------
* ``catmull_rom_spline`` produced NaN.  It padded the control points with
  duplicates (``[p[0]] + points + [p[-1]]``), which makes the first and last
  knot intervals zero-length; every basis term then divided by zero.  On a
  four-point input, 20 of 30 output coordinates were NaN.  Padding is now
  done by *reflection*, and degenerate intervals are guarded, so the result
  is finite for any input including duplicated points.
* The evaluation loop ran in Python per sample point; it is now vectorised
  over each segment.
* ``smooth_path_svg`` emitted ``Q x y x y`` (a control point identical to the
  endpoint, i.e. a straight line) and then chained ``T`` commands off it, so
  "smoothing" produced a polyline.  It now computes real Catmull-Rom-derived
  cubic control points.

Public API
----------
  cubic_bezier, quadratic_bezier, catmull_rom_spline
  route_edge_spline, bezier_to_svg_path, smooth_path_svg
  points_to_bezier_path, offset_endpoints
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

_EPS = 1e-9
Point = tuple[float, float]


# ---------------------------------------------------------------------------
# Bezier curves
# ---------------------------------------------------------------------------

def cubic_bezier(p0: Point, p1: Point, p2: Point, p3: Point,
                 n_points: int = 50) -> np.ndarray:
    """Sample a cubic Bezier defined by four control points."""
    t = np.linspace(0.0, 1.0, max(2, int(n_points)))[:, None]
    a, b, c, d = (np.asarray(p, float) for p in (p0, p1, p2, p3))
    return ((1 - t) ** 3 * a + 3 * (1 - t) ** 2 * t * b
            + 3 * (1 - t) * t ** 2 * c + t ** 3 * d)


def quadratic_bezier(p0: Point, p1: Point, p2: Point,
                     n_points: int = 50) -> np.ndarray:
    """Sample a quadratic Bezier defined by three control points."""
    t = np.linspace(0.0, 1.0, max(2, int(n_points)))[:, None]
    a, b, c = (np.asarray(p, float) for p in (p0, p1, p2))
    return (1 - t) ** 2 * a + 2 * (1 - t) * t * b + t ** 2 * c


# ---------------------------------------------------------------------------
# Catmull-Rom
# ---------------------------------------------------------------------------

def _reflect_pad(pts: np.ndarray) -> np.ndarray:
    """
    Pad with reflected phantom points.

    Duplicating the endpoints (what v2.x did) yields zero-length knot
    intervals and a division by zero in the Barry-Goldman recurrence.
    Reflecting keeps every interval strictly positive.
    """
    first = pts[0] + (pts[0] - pts[1])
    last = pts[-1] + (pts[-1] - pts[-2])
    return np.vstack([first, pts, last])


def catmull_rom_spline(points: Sequence[Point], n_points: int = 20,
                       alpha: float = 0.5) -> np.ndarray:
    """
    Interpolating Catmull-Rom spline through *points*.

    alpha=0 uniform, 0.5 centripetal (no cusps or self-intersections),
    1.0 chordal.  Returns an (N, 2) array that is always finite.
    """
    pts = np.asarray([(float(x), float(y)) for x, y in points], dtype=float)
    if pts.shape[0] == 0:
        return np.empty((0, 2))
    if pts.shape[0] == 1:
        return pts.copy()
    if pts.shape[0] == 2:
        t = np.linspace(0, 1, max(2, int(n_points)))[:, None]
        return pts[0] * (1 - t) + pts[1] * t

    # Collapse consecutive duplicates: they carry no shape information and
    # would create zero-length intervals.
    keep = [0] + [i for i in range(1, len(pts))
                  if np.linalg.norm(pts[i] - pts[i - 1]) > _EPS]
    pts = pts[keep]
    if pts.shape[0] < 3:
        t = np.linspace(0, 1, max(2, int(n_points)))[:, None]
        return pts[0] * (1 - t) + pts[-1] * t

    p = _reflect_pad(pts)
    n = max(2, int(n_points))
    out: list[np.ndarray] = []

    for i in range(len(p) - 3):
        p0, p1, p2, p3 = p[i], p[i + 1], p[i + 2], p[i + 3]

        def knot(prev: float, a_: np.ndarray, b_: np.ndarray) -> float:
            d = float(np.linalg.norm(b_ - a_))
            return prev + max(d, _EPS) ** alpha        # never zero-length

        t0 = 0.0
        t1 = knot(t0, p0, p1)
        t2 = knot(t1, p1, p2)
        t3 = knot(t2, p2, p3)

        t = np.linspace(t1, t2, n)[:, None]
        a1 = (t1 - t) / (t1 - t0) * p0 + (t - t0) / (t1 - t0) * p1
        a2 = (t2 - t) / (t2 - t1) * p1 + (t - t1) / (t2 - t1) * p2
        a3 = (t3 - t) / (t3 - t2) * p2 + (t - t2) / (t3 - t2) * p3
        b1 = (t2 - t) / (t2 - t0) * a1 + (t - t0) / (t2 - t0) * a2
        b2 = (t3 - t) / (t3 - t1) * a2 + (t - t1) / (t3 - t1) * a3
        seg = (t2 - t) / (t2 - t1) * b1 + (t - t1) / (t2 - t1) * b2

        out.append(seg if not out else seg[1:])        # avoid duplicate joints

    curve = np.vstack(out)
    return curve[np.isfinite(curve).all(axis=1)]


# ---------------------------------------------------------------------------
# Edge routing
# ---------------------------------------------------------------------------

def offset_endpoints(
    source: Point,
    target: Point,
    source_radius: float = 0.0,
    target_radius: float = 0.0,
) -> tuple[Point, Point]:
    """
    Shorten a segment so it starts and ends at node boundaries.

    Without this an arrowhead is hidden underneath the target node, which is
    why v2.x's graph view looked like it had no arrowheads.
    """
    s = np.asarray(source, float)
    t = np.asarray(target, float)
    d = t - s
    dist = float(np.linalg.norm(d))
    if dist < _EPS:
        return tuple(s), tuple(t)                      # type: ignore[return-value]
    u = d / dist
    s2 = s + u * min(source_radius, dist / 2)
    t2 = t - u * min(target_radius, dist / 2)
    return tuple(s2), tuple(t2)                        # type: ignore[return-value]


def route_edge_spline(
    source: Point,
    target: Point,
    obstacles: Sequence[tuple[float, float, float, float]] | None = None,
    routing_mode: str = "curved",
    curvature: float = 0.15,
    n_points: int = 30,
) -> np.ndarray:
    """
    Route an edge between two points.

    Modes: ``straight``, ``orthogonal`` (Manhattan with rounded corners),
    ``curved`` (arc bowed perpendicular to the chord), and ``avoid``
    (bows around the first intersecting obstacle).

    *obstacles* are (cx, cy, width, height) rectangles.
    """
    s = np.asarray(source, float)
    t = np.asarray(target, float)

    if routing_mode == "straight":
        return np.vstack([s, t])

    if routing_mode == "orthogonal":
        midx = (s[0] + t[0]) / 2.0
        return catmull_rom_spline(
            [tuple(s), (midx, s[1]), (midx, t[1]), tuple(t)], n_points=n_points
        )

    chord = t - s
    dist = float(np.linalg.norm(chord))
    if dist < _EPS:
        return np.vstack([s, t])
    normal = np.array([-chord[1], chord[0]]) / dist

    bow = curvature * dist
    if routing_mode == "avoid" and obstacles:
        for (cx, cy, w, h) in obstacles:
            if _segment_hits_rect(s, t, cx, cy, w, h):
                bow = max(bow, (max(w, h) / 2.0) + 12.0)
                break

    mid = (s + t) / 2.0 + normal * bow
    c1 = s + (mid - s) * 0.7
    c2 = t + (mid - t) * 0.7
    return cubic_bezier(tuple(s), tuple(c1), tuple(c2), tuple(t), n_points=n_points)


def _segment_hits_rect(s: np.ndarray, t: np.ndarray,
                       cx: float, cy: float, w: float, h: float,
                       samples: int = 24) -> bool:
    """Cheap sampled test for segment/rectangle intersection."""
    xs = np.linspace(s[0], t[0], samples)
    ys = np.linspace(s[1], t[1], samples)
    return bool(np.any(
        (np.abs(xs - cx) <= w / 2.0) & (np.abs(ys - cy) <= h / 2.0)
    ))


# ---------------------------------------------------------------------------
# SVG path emission
# ---------------------------------------------------------------------------

def bezier_to_svg_path(curve: np.ndarray, close: bool = False,
                       precision: int = 2) -> str:
    """Convert sampled points to an SVG path (``M`` then ``L`` commands)."""
    arr = np.asarray(curve, float)
    if arr.size == 0:
        return ""
    arr = arr[np.isfinite(arr).all(axis=1)]
    if arr.size == 0:
        return ""
    fmt = f"{{:.{precision}f}}"
    parts = ["M " + fmt.format(arr[0, 0]) + " " + fmt.format(arr[0, 1])]
    parts += ["L " + fmt.format(x) + " " + fmt.format(y) for x, y in arr[1:]]
    if close:
        parts.append("Z")
    return " ".join(parts)


def points_to_bezier_path(points: Sequence[Point], tension: float = 0.5,
                          precision: int = 2) -> str:
    """
    Emit a smooth cubic SVG path through *points*.

    Control points are derived from neighbouring knots (a Catmull-Rom to
    Bezier conversion), so the curve genuinely passes through every waypoint.
    """
    pts = [(float(x), float(y)) for x, y in points]
    if len(pts) < 2:
        return ""
    if len(pts) == 2:
        f = f"{{:.{precision}f}}"
        return (f"M {f.format(pts[0][0])} {f.format(pts[0][1])} "
                f"L {f.format(pts[1][0])} {f.format(pts[1][1])}")

    p = np.asarray(pts, float)
    ext = np.vstack([p[0] + (p[0] - p[1]), p, p[-1] + (p[-1] - p[-2])])
    k = max(0.0, min(1.0, float(tension))) / 3.0
    f = f"{{:.{precision}f}}"
    out = [f"M {f.format(p[0][0])} {f.format(p[0][1])}"]

    for i in range(1, len(ext) - 2):
        p0, p1, p2, p3 = ext[i - 1], ext[i], ext[i + 1], ext[i + 2]
        c1 = p1 + (p2 - p0) * k
        c2 = p2 - (p3 - p1) * k
        out.append(
            f"C {f.format(c1[0])} {f.format(c1[1])} "
            f"{f.format(c2[0])} {f.format(c2[1])} "
            f"{f.format(p2[0])} {f.format(p2[1])}"
        )
    return " ".join(out)


def smooth_path_svg(points: Sequence[Point], tension: float = 0.5) -> str:
    """Backwards-compatible alias for :func:`points_to_bezier_path`."""
    return points_to_bezier_path(points, tension=tension)
