"""
expansion.py
Splitting complexes and expanding multi-identifier nodes.

Why this matters
----------------
A KEGG signalling map represents a complex as a *group* entry whose
components are separate entries stacked inside one box, and represents a
family of paralogues as a *single* entry carrying several gene ids
("CDK4, CDK6").  Painting one colour across either loses information: a
complex where one subunit is up and another down reads as a single
aggregated value, and a paralogue family reads as whichever member the
aggregation happened to pick.

R pathview solves this with ``split.group`` and ``expand.node``.  This module
implements both, and adds the piece R does not: the expansion is *recorded*,
so a caller can tell which node a row came from and audit the figure.

  * ``split_groups``  — replace each group entry with its component entries,
    positioned inside the group's box.
  * ``expand_nodes``  — replace each multi-identifier entry with one entry per
    identifier, subdividing the original box.

Both are lossless with respect to geometry: the expanded nodes tile exactly
the area the original occupied, so the map's layout is unchanged.

Public API
----------
  split_groups   : group entries -> their components
  expand_nodes   : multi-id entries -> one entry per identifier
  expansion_report : what changed
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import polars as pl

from .constants import GENE_NODE_TYPES


@dataclass
class ExpansionResult:
    """What an expansion did, so a figure's provenance can be checked."""

    data: pl.DataFrame
    n_before: int
    n_after: int
    n_groups_split: int = 0
    n_nodes_expanded: int = 0

    def summary(self) -> str:
        parts = [f"{self.n_before} -> {self.n_after} nodes"]
        if self.n_groups_split:
            parts.append(f"{self.n_groups_split} complexes split")
        if self.n_nodes_expanded:
            parts.append(f"{self.n_nodes_expanded} multi-gene nodes expanded")
        return "; ".join(parts)


def _tile(x: float, y: float, w: float, h: float, n: int) -> list[tuple[float, float, float, float]]:
    """
    Divide a box into *n* sub-boxes that tile it exactly.

    Splits along the longer axis when the box is clearly oblong, otherwise
    uses a near-square grid, so an expanded node keeps the aspect ratio of
    the glyph it replaced instead of becoming a row of slivers.
    """
    n = max(1, int(n))
    if n == 1:
        return [(x, y, w, h)]

    if h >= w * 1.4:
        rows, cols = n, 1
    elif w >= h * 1.4:
        rows, cols = 1, n
    else:
        cols = int(np.ceil(np.sqrt(n)))
        rows = int(np.ceil(n / cols))

    ch = h / rows
    left, top = x - w / 2.0, y - h / 2.0
    out: list[tuple[float, float, float, float]] = []

    for r in range(rows):
        # A grid sized for `rows * cols` leaves the last row short whenever n
        # is not a multiple of cols — three items in a 2x2 grid would cover
        # only three quarters of the node.  Widening the final row's cells to
        # span the full width keeps the tiling exact for every n, which is
        # what makes an expanded node occupy precisely the area of the node it
        # replaced.
        in_row = min(cols, n - r * cols)
        if in_row <= 0:
            break
        cw = w / in_row
        for c in range(in_row):
            out.append((left + (c + 0.5) * cw, top + (r + 0.5) * ch, cw, ch))
    return out


def split_groups(
    node_data: pl.DataFrame,
    detailed: bool = False,
) -> pl.DataFrame | ExpansionResult:
    """
    Replace each group (complex) entry with its component entries.

    KGML group entries list their members in the ``component`` column and the
    members are themselves entries, so splitting means dropping the group and
    keeping the components — which KEGG already positions inside the group
    box.  Components not present as their own entries are synthesised from the
    group's identifiers so nothing is lost.

    Adds a ``parent_group`` column recording where each row came from.
    """
    if node_data is None or node_data.is_empty():
        return ExpansionResult(node_data, 0, 0) if detailed else node_data

    if "component" not in node_data.columns:
        out = node_data.with_columns(pl.lit(None, pl.String).alias("parent_group"))
        return ExpansionResult(out, out.height, out.height) if detailed else out

    groups = node_data.filter(
        (pl.col("type") == "group") & (pl.col("component") != "")
    )
    if groups.is_empty():
        out = node_data.with_columns(pl.lit(None, pl.String).alias("parent_group"))
        return ExpansionResult(out, out.height, out.height) if detailed else out

    by_id = {row["entry_id"]: row for row in node_data.iter_rows(named=True)}
    parent_of: dict[str, str] = {}
    synthesised: list[dict] = []

    for grp in groups.iter_rows(named=True):
        members = [m for m in str(grp["component"]).split(";") if m]
        present = [m for m in members if m in by_id]
        for m in present:
            parent_of[m] = grp["entry_id"]

        missing = [m for m in members if m not in by_id]
        if missing or not present:
            # A group whose components are not separate entries still has to
            # yield something: tile the group box over its own identifiers.
            names = grp.get("kegg_names") or []
            boxes = _tile(grp["x"] or 0.0, grp["y"] or 0.0,
                          grp["width"] or 46.0, grp["height"] or 17.0,
                          max(1, len(names)))
            for i, (nm, (bx, by, bw, bh)) in enumerate(zip(names, boxes)):
                row = dict(grp)
                row.update({
                    "entry_id": f"{grp['entry_id']}.g{i}",
                    "type": "gene", "kegg_names": [nm],
                    "x": bx, "y": by, "width": bw, "height": bh,
                    "component": "", "size": 1,
                    "label": nm,
                })
                synthesised.append(row)
                parent_of[row["entry_id"]] = grp["entry_id"]

    kept = node_data.filter(
        ~((pl.col("type") == "group") & (pl.col("component") != ""))
    )
    if synthesised:
        kept = pl.concat(
            [kept, pl.DataFrame(synthesised, schema=node_data.schema)],
            how="vertical_relaxed",
        )

    out = kept.with_columns(
        pl.col("entry_id")
        .map_elements(lambda e: parent_of.get(e), return_dtype=pl.String)
        .alias("parent_group")
    )
    res = ExpansionResult(out, node_data.height, out.height,
                          n_groups_split=groups.height)
    return res if detailed else out


def expand_nodes(
    node_data: pl.DataFrame,
    node_types: list[str] | None = None,
    max_per_node: int = 12,
    detailed: bool = False,
) -> pl.DataFrame | ExpansionResult:
    """
    Give every identifier on a multi-identifier node its own sub-node.

    A KEGG entry such as ``CDK4, CDK6`` becomes two nodes tiling the original
    box, so paralogues that disagree are visible instead of being collapsed by
    ``node_sum``.  Entries with more than *max_per_node* identifiers are left
    alone — a 40-way subdivision is unreadable, and silently producing one
    would be worse than not expanding.

    Adds ``parent_entry`` (the original entry id) and ``expanded`` columns.
    """
    if node_data is None or node_data.is_empty():
        return ExpansionResult(node_data, 0, 0) if detailed else node_data
    if "kegg_names" not in node_data.columns:
        raise ValueError("expand_nodes needs the kegg_names column from node_info().")

    wanted = set(node_types) if node_types else set(GENE_NODE_TYPES) | {"compound"}

    rows: list[dict] = []
    n_expanded = 0
    for row in node_data.iter_rows(named=True):
        names = list(row.get("kegg_names") or [])
        eligible = (row["type"] in wanted and 1 < len(names) <= max_per_node
                    and row.get("x") is not None)
        if not eligible:
            out = dict(row)
            out["parent_entry"] = row["entry_id"]
            out["expanded"] = False
            rows.append(out)
            continue

        boxes = _tile(row["x"], row["y"], row["width"] or 46.0,
                      row["height"] or 17.0, len(names))
        for i, (nm, (bx, by, bw, bh)) in enumerate(zip(names, boxes)):
            out = dict(row)
            out.update({
                "entry_id": f"{row['entry_id']}.{i}",
                "kegg_names": [nm],
                "x": bx, "y": by, "width": bw, "height": bh,
                "size": 1, "label": nm,
                "parent_entry": row["entry_id"],
                "expanded": True,
            })
            rows.append(out)
        n_expanded += 1

    schema = dict(node_data.schema)
    schema["parent_entry"] = pl.String
    schema["expanded"] = pl.Boolean
    out = pl.DataFrame(rows, schema=schema)

    res = ExpansionResult(out, node_data.height, out.height,
                          n_nodes_expanded=n_expanded)
    return res if detailed else out


def expansion_report(before: pl.DataFrame, after: pl.DataFrame) -> dict:
    """Compare two node frames — useful for asserting an expansion did something."""
    return {
        "nodes_before": before.height,
        "nodes_after": after.height,
        "added": after.height - before.height,
        "groups_before": before.filter(pl.col("type") == "group").height,
        "groups_after": after.filter(pl.col("type") == "group").height,
    }
