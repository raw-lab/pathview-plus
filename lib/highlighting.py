"""
highlighting.py
Composable post-hoc pathway modification.

Fixes over v2.x
---------------
* ``pathview()`` returned a plain ``dict``, so the entire documented
  ``result + highlight_nodes(...)`` API could never be used — ``dict`` has no
  ``__add__``.  The orchestrator now returns a :class:`PathwayResult`.
* ``image_array`` was never populated by anything, so every modifier hit its
  ``if result.image_array is None: return`` guard and silently did nothing.
  The result now carries the rendered raster.
* ``_hex_to_rgb("red")`` raised ``ValueError: invalid literal for int() with
  base 16: 're'`` — the function crashed on ``highlight_nodes``' own default
  argument.  Colour parsing goes through :func:`utils.to_rgb`.
* Highlights were drawn at ``img_height - y``, mirroring them vertically away
  from the nodes they were meant to mark.
* ``opacity`` was accepted and ignored; ``change_labels`` stored a dict and
  never rendered anything.  Both now do what they say.

Usage
-----
    result = pathview("00020", gene_data=df, species="hsa")
    (result
     + highlight_nodes(["1431", "3417"], color="#7C3AED", width=3)
     + highlight_path(["1431", "3417", "3418"], color="orange")
     + change_labels({"1431": "CS *"})
    ).save("annotated.png")
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import polars as pl

from .layout import NodeBox, node_boxes
from .utils import to_rgb

Modifier = Callable[["PathwayResult"], None]


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class PathwayResult:
    """
    Everything a pathview run produced, plus a composable modification API.

    Also behaves like the v2.x return dict (``result["plot_data_gene"]``) so
    existing scripts keep working.
    """

    pathway_id: str = ""
    pathway_name: str = ""
    species: str = ""
    plot_data_gene: pl.DataFrame | None = None
    plot_data_cpd: pl.DataFrame | None = None
    cols_gene: pl.DataFrame | None = None
    cols_cpd: pl.DataFrame | None = None
    node_data: pl.DataFrame | None = None
    edge_data: pl.DataFrame | None = None
    output_path: Path | None = None
    frame: object = None          # RasterFrame: raster + KGML transform
    gene_scale: object = None
    cpd_scale: object = None
    diagnostics: dict = field(default_factory=dict)
    modifications: list[Modifier] = field(default_factory=list)
    label_changes: dict[str, str] = field(default_factory=dict)

    # -- dict compatibility ------------------------------------------------
    def __getitem__(self, key: str):
        return getattr(self, key)

    def get(self, key: str, default=None):
        return getattr(self, key, default)

    def keys(self) -> list[str]:
        return ["plot_data_gene", "plot_data_cpd", "node_data", "edge_data",
                "output_path", "diagnostics"]

    def __bool__(self) -> bool:
        return self.node_data is not None and not self.node_data.is_empty()

    # -- composition -------------------------------------------------------
    def __add__(self, modifier: Modifier) -> PathwayResult:
        if not callable(modifier):
            raise TypeError(
                f"Expected a modifier function, got {type(modifier).__name__}. "
                "Use highlight_nodes(...) rather than highlight_nodes."
            )
        new = PathwayResult(
            pathway_id=self.pathway_id, pathway_name=self.pathway_name,
            species=self.species,
            plot_data_gene=self.plot_data_gene, plot_data_cpd=self.plot_data_cpd,
            cols_gene=self.cols_gene, cols_cpd=self.cols_cpd,
            node_data=self.node_data, edge_data=self.edge_data,
            output_path=self.output_path,
            frame=None if self.frame is None else self.frame.copy(),
            gene_scale=self.gene_scale, cpd_scale=self.cpd_scale,
            diagnostics=dict(self.diagnostics),
            modifications=self.modifications + [modifier],
            label_changes=dict(self.label_changes),
        )
        modifier(new)
        return new

    @property
    def image_array(self) -> np.ndarray | None:
        """The rendered raster, or None when the mode produced only vectors."""
        return None if self.frame is None else self.frame.array

    @image_array.setter
    def image_array(self, value) -> None:
        from .layout import RasterFrame
        self.frame = None if value is None else RasterFrame(value)

    # -- output ------------------------------------------------------------
    def save(self, path: str | Path, dpi: int = 200) -> Path:
        """Write the (possibly modified) raster to *path*."""
        from PIL import Image

        if self.frame is None:
            raise ValueError(
                "No raster to save. Highlighting operates on a rendered "
                "image; run pathview(..., render_mode='native') or "
                "render_mode='vector' first."
            )
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        img = Image.fromarray(self.frame.array.astype(np.uint8))
        fmt = out.suffix.lstrip(".").upper() or "PNG"
        if fmt == "PDF":
            img.convert("RGB").save(out, "PDF", resolution=float(dpi))
        else:
            img.save(out, "JPEG" if fmt in ("JPG", "JPEG") else fmt)
        return out

    def summary(self) -> str:
        d = self.diagnostics
        parts = [f"{self.pathway_name or self.pathway_id}"]
        if "gene" in d:
            parts.append(f"genes: {d['gene']}")
        if "cpd" in d:
            parts.append(f"compounds: {d['cpd']}")
        if self.output_path:
            parts.append(f"-> {self.output_path.name}")
        return " | ".join(parts)

    def _boxes(self) -> list[NodeBox]:
        return node_boxes(self.node_data) if self.node_data is not None else []


class PathwayResultSet:
    """
    The results of rendering several pathways in one call.

    Behaves like an ordered mapping of pathway id to
    :class:`PathwayResult`, iterates over the results, and broadcasts ``+``
    to every member so a highlight can be applied across a whole batch.

    Failures are kept rather than raised: one unavailable pathway in a batch
    of twenty should not discard the nineteen that worked.  ``failures``
    records what went wrong, and the set is falsy only when nothing succeeded.
    """

    def __init__(self, results: dict[str, PathwayResult] | None = None,
                 failures: dict[str, str] | None = None):
        self._results: dict[str, PathwayResult] = dict(results or {})
        self.failures: dict[str, str] = dict(failures or {})

    # -- container protocol ------------------------------------------------
    def __len__(self) -> int:
        return len(self._results)

    def __iter__(self):
        return iter(self._results.values())

    def __getitem__(self, key):
        if isinstance(key, int):
            return list(self._results.values())[key]
        return self._results[key]

    def __contains__(self, key) -> bool:
        return key in self._results

    def __bool__(self) -> bool:
        return bool(self._results)

    def __repr__(self) -> str:
        return (f"PathwayResultSet({len(self._results)} rendered, "
                f"{len(self.failures)} failed)")

    def keys(self):
        return self._results.keys()

    def values(self):
        return self._results.values()

    def items(self):
        return self._results.items()

    def get(self, key, default=None):
        return self._results.get(key, default)

    # -- convenience -------------------------------------------------------
    @property
    def output_paths(self) -> list[Path]:
        return [r.output_path for r in self._results.values() if r.output_path]

    def __add__(self, modifier: Modifier) -> PathwayResultSet:
        """Apply a modifier to every result, returning a new set."""
        return PathwayResultSet(
            {k: (v + modifier) for k, v in self._results.items()},
            self.failures,
        )

    def summary(self) -> str:
        lines = [r.summary() for r in self._results.values()]
        lines += [f"{k}: FAILED — {v}" for k, v in self.failures.items()]
        head = f"{len(self._results)} pathway(s) rendered"
        if self.failures:
            head += f", {len(self.failures)} failed"
        return "\n".join([head, *lines])

    def to_frame(self):
        """One row per pathway: id, name, output path, status."""
        import polars as pl

        rows = [{"pathway": k, "name": r.pathway_name,
                 "output": str(r.output_path) if r.output_path else "",
                 "status": "ok"} for k, r in self._results.items()]
        rows += [{"pathway": k, "name": "", "output": "", "status": v}
                 for k, v in self.failures.items()]
        return pl.DataFrame(rows) if rows else pl.DataFrame(
            schema={"pathway": pl.String, "name": pl.String,
                    "output": pl.String, "status": pl.String})


# ---------------------------------------------------------------------------
# Node resolution
# ---------------------------------------------------------------------------

def _resolve(result: PathwayResult, ids: Sequence[str]) -> list[NodeBox]:
    """
    Find node boxes for *ids*, which may be entry ids or biological ids.

    Accepting both matters: users think in Entrez/KEGG accessions, while the
    layout is keyed on KGML entry ids.
    """
    wanted = {str(i) for i in ids}
    boxes = result._boxes()
    by_entry = {b.entry_id: b for b in boxes}

    hits: list[NodeBox] = [by_entry[i] for i in wanted if i in by_entry]

    nd = result.node_data
    if nd is not None and "kegg_names" in nd.columns:
        for row in nd.iter_rows(named=True):
            names = row.get("kegg_names") or []
            if any(str(n) in wanted for n in names):
                b = by_entry.get(str(row["entry_id"]))
                if b is not None and b not in hits:
                    hits.append(b)
    return hits


# ---------------------------------------------------------------------------
# Drawing primitives (top-left origin throughout — no flipping)
# ---------------------------------------------------------------------------

def _blend(img: np.ndarray, ys: slice, xs: slice,
           rgb: tuple[int, int, int], opacity: float) -> None:
    region = img[ys, xs, :3]
    if region.size == 0:
        return
    a = float(np.clip(opacity, 0.0, 1.0))
    if a >= 1.0:
        region[:] = rgb
    else:
        region[:] = np.clip(region.astype(np.float32) * (1 - a)
                            + np.asarray(rgb, np.float32) * a, 0, 255).astype(np.uint8)
    img[ys, xs, :3] = region


def _draw_border(frame, box: NodeBox, rgb: tuple[int, int, int],
                 thickness: int, opacity: float) -> None:
    """
    Outline a node.

    Coordinates go through the frame's transform: KGML y already points down,
    so there is no flip, but the raster may be offset and scaled relative to
    KGML space and the highlight must follow.
    """
    img = frame.array
    h, w = img.shape[:2]
    t = max(1, int(round(thickness * max(1.0, frame.scale))))
    px0, py0 = frame.to_pixels(box.left, box.top)
    px1, py1 = frame.to_pixels(box.right, box.bottom)
    x0, x1 = int(round(px0)) - t, int(round(px1)) + t
    y0, y1 = int(round(py0)) - t, int(round(py1)) + t
    x0c, x1c = max(0, x0), min(w, x1)
    y0c, y1c = max(0, y0), min(h, y1)
    if x1c <= x0c or y1c <= y0c:
        return
    _blend(img, slice(y0c, min(h, y0c + t)), slice(x0c, x1c), rgb, opacity)
    _blend(img, slice(max(0, y1c - t), y1c), slice(x0c, x1c), rgb, opacity)
    _blend(img, slice(y0c, y1c), slice(x0c, min(w, x0c + t)), rgb, opacity)
    _blend(img, slice(y0c, y1c), slice(max(0, x1c - t), x1c), rgb, opacity)


def _draw_line(frame, p0: tuple[float, float], p1: tuple[float, float],
               rgb: tuple[int, int, int], thickness: int, opacity: float) -> None:
    """Draw a thick line by sampling along the segment (no coordinate flip)."""
    img = frame.array
    h, w = img.shape[:2]
    x0, y0 = frame.to_pixels(*p0)
    x1, y1 = frame.to_pixels(*p1)
    thickness = max(1, int(round(thickness * max(1.0, frame.scale))))
    steps = int(max(abs(x1 - x0), abs(y1 - y0))) + 1
    t = max(1, int(thickness)) // 2
    for xf, yf in zip(np.linspace(x0, x1, steps), np.linspace(y0, y1, steps)):
        px, py = int(round(xf)), int(round(yf))
        _blend(img, slice(max(0, py - t), min(h, py + t + 1)),
               slice(max(0, px - t), min(w, px + t + 1)), rgb, opacity)


# ---------------------------------------------------------------------------
# Modifiers
# ---------------------------------------------------------------------------

def highlight_nodes(node_ids: Sequence[str], color: str = "red", width: int = 3,
                    opacity: float = 1.0) -> Modifier:
    """Outline the given nodes.  Accepts named colours, hex, or RGB tuples."""
    rgb = to_rgb(color, default=(220, 38, 38))

    def modifier(result: PathwayResult) -> None:
        if result.frame is None:
            return
        for box in _resolve(result, node_ids):
            _draw_border(result.frame, box, rgb, width, opacity)

    return modifier


def highlight_edges(edge_pairs: Sequence[tuple[str, str]], color: str = "blue",
                    width: int = 3, opacity: float = 1.0) -> Modifier:
    """Draw a line between each pair of nodes."""
    rgb = to_rgb(color, default=(37, 99, 235))

    def modifier(result: PathwayResult) -> None:
        if result.frame is None:
            return
        for src, tgt in edge_pairs:
            s = _resolve(result, [src])
            t = _resolve(result, [tgt])
            if not s or not t:
                continue
            _draw_line(result.frame, (s[0].x, s[0].y), (t[0].x, t[0].y),
                       rgb, width, opacity)

    return modifier


def highlight_path(path_node_ids: Sequence[str], color: str = "orange",
                   node_width: int = 3, edge_width: int = 3,
                   opacity: float = 1.0) -> Modifier:
    """Outline an ordered chain of nodes and connect them."""
    pairs = list(zip(path_node_ids[:-1], path_node_ids[1:]))

    def modifier(result: PathwayResult) -> None:
        highlight_edges(pairs, color=color, width=edge_width, opacity=opacity)(result)
        highlight_nodes(path_node_ids, color=color, width=node_width,
                        opacity=opacity)(result)

    return modifier


def change_labels(label_map: dict[str, str], font_size: int = 9,
                  color: str = "#111111") -> Modifier:
    """
    Replace node labels on the rendered image.

    v2.x stored the mapping on an undeclared attribute and never drew
    anything.  This actually repaints the node and writes the new text.
    """
    rgb = to_rgb(color, default=(17, 17, 17))

    def modifier(result: PathwayResult) -> None:
        result.label_changes.update({str(k): str(v) for k, v in label_map.items()})
        if result.frame is None:
            return
        frame = result.frame
        try:
            from PIL import Image, ImageDraw, ImageFont
        except ImportError:                                   # pragma: no cover
            return

        img = Image.fromarray(frame.array.astype(np.uint8))
        draw = ImageDraw.Draw(img)
        size = max(6, int(round(font_size * max(1.0, frame.scale))))
        try:
            font = ImageFont.truetype(
                "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", size)
        except OSError:
            font = ImageFont.load_default()

        for key, text in label_map.items():
            for box in _resolve(result, [key]):
                left, top = frame.to_pixels(box.left, box.top)
                right, bottom = frame.to_pixels(box.right, box.bottom)
                cx, cy = frame.to_pixels(box.x, box.y)
                draw.rectangle([left + 1, top + 1, right - 1, bottom - 1],
                               fill=(255, 255, 255))
                draw.rectangle([left, top, right, bottom],
                               outline=(70, 70, 70), width=1)
                try:
                    bb = draw.textbbox((0, 0), text, font=font)
                    tw, th_ = bb[2] - bb[0], bb[3] - bb[1]
                except Exception:                             # pragma: no cover
                    tw, th_ = len(text) * size * 0.55, size
                draw.text((cx - tw / 2, cy - th_ / 2), text,
                          fill=tuple(rgb), font=font)

        frame.array = np.array(img, dtype=np.uint8)

    return modifier


def annotate(text: str, xy: tuple[float, float], color: str = "#111111",
             font_size: int = 11) -> Modifier:
    """Write free text at a pathway coordinate."""
    rgb = to_rgb(color, default=(17, 17, 17))

    def modifier(result: PathwayResult) -> None:
        if result.frame is None:
            return
        from PIL import Image, ImageDraw, ImageFont

        frame = result.frame
        img = Image.fromarray(frame.array.astype(np.uint8))
        draw = ImageDraw.Draw(img)
        size = max(6, int(round(font_size * max(1.0, frame.scale))))
        try:
            font = ImageFont.truetype(
                "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", size)
        except OSError:
            font = ImageFont.load_default()
        draw.text(frame.to_pixels(*xy), text, fill=tuple(rgb), font=font)
        frame.array = np.array(img, dtype=np.uint8)

    return modifier
