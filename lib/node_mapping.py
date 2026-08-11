"""
node_mapping.py
Join molecular data onto pathway nodes.

Fixes over v2.x
---------------
* The prefix-stripping regex was ``^[a-z]+:``, which misses organism codes
  containing digits (``taes:``, ``hsa4:``) and uppercase namespaces.  KEGG
  identifiers are now split at parse time and matched exactly.
* v2.x re-split the space-separated ``name`` field on every call and applied a
  no-op double ``rename``; ``kegg_names`` is now a real list column produced
  once by the parser.
* Multi-ID nodes reported only the first identifier.  ``all_mapped`` now
  records every input ID that landed on the node, which is what R pathview
  shows and what users need in order to audit a figure.
* Returns node-level diagnostics so a caller can tell "nothing mapped" from
  "no nodes of this type exist".

Public API
----------
  node_map      : mol_data + node_data -> per-node values
  NodeMapResult : mapping diagnostics
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import polars as pl

from .constants import CPD_NODE_TYPES, GENE_NODE_TYPES, SumMethod
from .errors import MappingError
from .mol_data import mol_sum


@dataclass
class NodeMapResult:
    """Outcome of mapping one data class onto one set of node types."""

    data: pl.DataFrame | None
    n_nodes: int = 0
    n_nodes_with_data: int = 0
    n_ids_input: int = 0
    n_ids_mapped: int = 0
    unmapped_ids: list[str] = None            # type: ignore[assignment]
    value_columns: list[str] = None           # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.unmapped_ids is None:
            self.unmapped_ids = []
        if self.value_columns is None:
            self.value_columns = []

    @property
    def mapped_fraction(self) -> float:
        """Fraction of input identifiers that landed on a node."""
        return (self.n_ids_mapped / self.n_ids_input) if self.n_ids_input else 0.0

    @property
    def ok(self) -> bool:
        return self.data is not None and self.n_nodes > 0

    def summary(self) -> str:
        if not self.ok:
            return "no mappable nodes"
        return (f"{self.n_nodes_with_data}/{self.n_nodes} nodes carry data; "
                f"{self.n_ids_mapped}/{self.n_ids_input} input IDs used")


def _explode_names(node_data: pl.DataFrame) -> pl.DataFrame:
    """One row per (entry_id, kegg identifier)."""
    if "kegg_names" not in node_data.columns:
        return (node_data
                .with_columns(pl.col("name").str.split(" ").alias("kegg_names"))
                .explode("kegg_names")
                .with_columns(
                    pl.col("kegg_names").str.replace(r"^[A-Za-z][A-Za-z0-9_]*:", "")
                ))
    return node_data.explode("kegg_names").filter(pl.col("kegg_names").is_not_null())


def node_map(
    mol_data: pl.DataFrame | None,
    node_data: pl.DataFrame,
    node_types: str | Sequence[str] = "gene",
    node_sum: SumMethod = "sum",
    rand_seed: int | None = None,
    detailed: bool = False,
) -> pl.DataFrame | None | NodeMapResult:
    """
    Map *mol_data* onto the nodes of *node_data* whose type is in *node_types*.

    Returns the node layout plus one column per experiment.  When *mol_data*
    is None the layout is returned with no value columns, so unmapped maps
    still render (R pathview's ``map.null`` behaviour).
    """
    if isinstance(node_types, str):
        wanted = {node_types}
    else:
        wanted = set(node_types)
    # "gene" implies the KEGG gene-like types; "compound" implies compounds.
    if wanted == {"gene"}:
        wanted = set(GENE_NODE_TYPES)
    elif wanted == {"compound"}:
        wanted = set(CPD_NODE_TYPES)

    if node_data is None or node_data.is_empty():
        return NodeMapResult(None) if detailed else None

    targets = node_data.filter(pl.col("type").is_in(list(wanted)))
    if targets.is_empty():
        return NodeMapResult(None) if detailed else None

    if mol_data is None:
        res = NodeMapResult(data=targets, n_nodes=targets.height)
        return res if detailed else targets

    exploded = _explode_names(targets)
    if exploded.is_empty():
        return NodeMapResult(None, n_nodes=targets.height) if detailed else None

    id_col = mol_data.columns[0]
    id_map = exploded.select([
        pl.col("kegg_names").cast(pl.String).alias(id_col),
        pl.col("entry_id").alias("__target"),
    ]).drop_nulls().unique()

    try:
        summed = mol_sum(mol_data, id_map, sum_method=node_sum,
                         rand_seed=rand_seed, detailed=True)
    except MappingError:
        res = NodeMapResult(None, n_nodes=targets.height,
                            n_ids_input=mol_data.height)
        return res if detailed else None

    values = summed.data.rename({summed.data.columns[0]: "entry_id"})

    # Record which input identifiers actually landed on each node.
    input_ids = set(mol_data[id_col].cast(pl.String).to_list())
    all_mapped = (
        exploded
        .filter(pl.col("kegg_names").is_in(list(input_ids)))
        .group_by("entry_id")
        .agg(pl.col("kegg_names").unique().sort().str.join(",").alias("all_mapped"))
    )

    plot_data = (targets
                 .join(values, on="entry_id", how="left")
                 .join(all_mapped, on="entry_id", how="left"))

    value_cols = [c for c in values.columns if c != "entry_id"]
    with_data = plot_data.filter(
        pl.any_horizontal([pl.col(c).is_not_null() for c in value_cols])
    ).height if value_cols else 0

    res = NodeMapResult(
        data=plot_data,
        n_nodes=targets.height,
        n_nodes_with_data=with_data,
        n_ids_input=summed.n_input,
        n_ids_mapped=summed.n_mapped,
        unmapped_ids=summed.unmapped_ids,
        value_columns=value_cols,
    )
    return res if detailed else plot_data
