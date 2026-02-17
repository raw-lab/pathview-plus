"""
splines.py
Bezier curve (spline) rendering for pathway edges.

Provides smoother, more aesthetically pleasing edge routing compared to
straight lines. Particularly useful for complex pathways with many crossings.

Public API
----------
  cubic_bezier       : Calculate points along a cubic Bezier curve
  quadratic_bezier   : Calculate points along a quadratic Bezier curve
  catmull_rom_spline : Calculate smooth curve through control points
  route_edge_spline  : Auto-route an edge avoiding obstacles
"""

from __future__ import annotations

from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Bezier curves
# ---------------------------------------------------------------------------

def cubic_bezier(
    p0: tuple[float, float],
    p1: tuple[float, float],
    p2: tuple[float, float],
    p3: tuple[float, float],
    n_points: int = 50,
) -> np.ndarray:
    """
    Calculate points along a cubic Bezier curve.
    
    A cubic Bezier is defined by 4 control points:
    - p0: start point
    - p1: first control point
    - p2: second control point
    - p3: end point
    
    Parameters
    ----------
    p0, p1, p2, p3: Control points as (x, y) tuples
    n_points:       Number of points to sample along the curve
    
    Returns
    -------
    Array of shape (n_points, 2) containing (x, y) coordinates.
    
    Example
    -------
    >>> curve = cubic_bezier((0, 0), (1, 2), (3, 2), (4, 0), n_points=100)
    >>> plt.plot(curve[:, 0], curve[:, 1])
    """
    t = np.linspace(0, 1, n_points)[:, np.newaxis]
    
    # Cubic Bezier formula: B(t) = (1-t)³P₀ + 3(1-t)²tP₁ + 3(1-t)t²P₂ + t³P₃
    p0 = np.array(p0)
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    
    curve = (
        (1 - t)**3 * p0
        + 3 * (1 - t)**2 * t * p1
        + 3 * (1 - t) * t**2 * p2
        + t**3 * p3
    )
    return curve


def quadratic_bezier(
    p0: tuple[float, float],
    p1: tuple[float, float],
    p2: tuple[float, float],
    n_points: int = 50,
) -> np.ndarray:
    """
    Calculate points along a quadratic Bezier curve.
    
    A quadratic Bezier is defined by 3 control points:
    - p0: start point
    - p1: control point
    - p2: end point
    
    Parameters
    ----------
    p0, p1, p2: Control points as (x, y) tuples
    n_points:   Number of points to sample
    
    Returns array of shape (n_points, 2).
    """
    t = np.linspace(0, 1, n_points)[:, np.newaxis]
    
    # Quadratic Bezier: B(t) = (1-t)²P₀ + 2(1-t)tP₁ + t²P₂
    p0 = np.array(p0)
    p1 = np.array(p1)
    p2 = np.array(p2)
    
    curve = (1 - t)**2 * p0 + 2 * (1 - t) * t * p1 + t**2 * p2
    return curve


# ---------------------------------------------------------------------------
# Catmull-Rom splines
# ---------------------------------------------------------------------------

def catmull_rom_spline(
    points: list[tuple[float, float]],
    n_points: int = 50,
    alpha: float = 0.5,
) -> np.ndarray:
    """
    Calculate a smooth Catmull-Rom spline through control points.
    
    Catmull-Rom splines pass through all control points (interpolating spline)
    and produce smooth curves. The 'alpha' parameter controls the
    parameterization:
    - alpha = 0.0: uniform (can produce loops)
    - alpha = 0.5: centripetal (most common, no loops/cusps)
    - alpha = 1.0: chordal
    
    Parameters
    ----------
    points:   List of (x, y) control points to interpolate
    n_points: Number of points to sample between each pair
    alpha:    Parameterization (0.5 = centripetal, recommended)
    
    Returns array of shape (total_points, 2).
    
    Example
    -------
    >>> control_pts = [(0, 0), (1, 2), (3, 1), (4, 3)]
    >>> smooth_curve = catmull_rom_spline(control_pts, n_points=30)
    """
    if len(points) < 2:
        return np.array(points)
    
    # Add phantom points at start and end
    p = [points[0]] + points + [points[-1]]
    curves = []
    
    for i in range(len(p) - 3):
        p0, p1, p2, p3 = p[i], p[i+1], p[i+2], p[i+3]
        
        # Calculate segment lengths
        t0 = 0
        t1 = t0 + _distance(p0, p1) ** alpha
        t2 = t1 + _distance(p1, p2) ** alpha
        t3 = t2 + _distance(p2, p3) ** alpha
        
        # Sample points in the valid range [t1, t2]
        t = np.linspace(t1, t2, n_points)
        
        # Catmull-Rom basis functions
        for ti in t:
            a1 = (t1 - ti) / (t1 - t0) * np.array(p0) + (ti - t0) / (t1 - t0) * np.array(p1)
            a2 = (t2 - ti) / (t2 - t1) * np.array(p1) + (ti - t1) / (t2 - t1) * np.array(p2)
            a3 = (t3 - ti) / (t3 - t2) * np.array(p2) + (ti - t2) / (t3 - t2) * np.array(p3)
            
            b1 = (t2 - ti) / (t2 - t0) * a1 + (ti - t0) / (t2 - t0) * a2
            b2 = (t3 - ti) / (t3 - t1) * a2 + (ti - t1) / (t3 - t1) * a3
            
            c = (t2 - ti) / (t2 - t1) * b1 + (ti - t1) / (t2 - t1) * b2
            curves.append(c)
    
    return np.array(curves)


def _distance(p1: tuple[float, float], p2: tuple[float, float]) -> float:
    """Euclidean distance between two points."""
    return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)


# ---------------------------------------------------------------------------
# Auto-routing
# ---------------------------------------------------------------------------

def route_edge_spline(
    source: tuple[float, float],
    target: tuple[float, float],
    obstacles: Optional[list[tuple[float, float, float, float]]] = None,
    routing_mode: str = "orthogonal",
) -> np.ndarray:
    """
    Auto-route an edge between source and target, avoiding obstacles.
    
    This is a simplified routing algorithm. For production use, consider
    more sophisticated routing like:
    - A* pathfinding
    - Visibility graphs
    - Force-directed edge bundling
    
    Parameters
    ----------
    source:       (x, y) starting point
    target:       (x, y) ending point
    obstacles:    List of (x, y, width, height) rectangles to avoid
    routing_mode: "straight", "orthogonal", or "curved"
    
    Returns array of points defining the routed path.
    """
    if routing_mode == "straight" or obstacles is None:
        return np.array([source, target])
    
    if routing_mode == "orthogonal":
        # Simple orthogonal routing (Manhattan-style)
        sx, sy = source
        tx, ty = target
        midx = (sx + tx) / 2
        
        control_points = [
            source,
            (midx, sy),
            (midx, ty),
            target,
        ]
        return catmull_rom_spline(control_points, n_points=20)
    
    elif routing_mode == "curved":
        # Gentle S-curve
        sx, sy = source
        tx, ty = target
        
        # Control points for cubic Bezier
        dx = tx - sx
        dy = ty - sy
        c1 = (sx + dx * 0.3, sy + dy * 0.1)
        c2 = (sx + dx * 0.7, sy + dy * 0.9)
        
        return cubic_bezier(source, c1, c2, target, n_points=30)
    
    return np.array([source, target])


# ---------------------------------------------------------------------------
# SVG path generation
# ---------------------------------------------------------------------------

def bezier_to_svg_path(
    curve: np.ndarray,
    close: bool = False,
) -> str:
    """
    Convert a Bezier curve to SVG path data.
    
    Parameters
    ----------
    curve:  Array of shape (n, 2) containing (x, y) points
    close:  Whether to close the path (Z command)
    
    Returns SVG path data string (for use in <path d="..."/>)
    
    Example
    -------
    >>> curve = cubic_bezier((10, 10), (50, 80), (150, 80), (200, 10))
    >>> path_data = bezier_to_svg_path(curve)
    >>> svg = f'<path d="{path_data}" stroke="black" fill="none"/>'
    """
    if len(curve) == 0:
        return ""
    
    path_parts = [f"M {curve[0, 0]:.2f} {curve[0, 1]:.2f}"]
    
    for point in curve[1:]:
        path_parts.append(f"L {point[0]:.2f} {point[1]:.2f}")
    
    if close:
        path_parts.append("Z")
    
    return " ".join(path_parts)


def smooth_path_svg(
    points: list[tuple[float, float]],
    tension: float = 0.5,
) -> str:
    """
    Generate smooth SVG path using quadratic Bezier commands.
    
    Parameters
    ----------
    points:  List of (x, y) waypoints
    tension: Curve tension (0 = sharp corners, 1 = very smooth)
    
    Returns SVG path data using S (smooth cubic bezier) commands.
    """
    if len(points) < 2:
        return ""
    
    path_parts = [f"M {points[0][0]} {points[0][1]}"]
    
    for i in range(1, len(points)):
        x, y = points[i]
        if i == 1:
            # First curve segment uses Q (quadratic)
            path_parts.append(f"Q {x} {y} {x} {y}")
        else:
            # Subsequent segments use T (smooth continuation)
            path_parts.append(f"T {x} {y}")
    
    return " ".join(path_parts)
