"""
node_mapping.py
Map molecular expression / abundance data onto KEGG pathway nodes:
  - node_map : join mol_data to node_data via KEGG gene/compound IDs
"""

from __future__ import annotations

from typing import Optional

import polars as pl

from .constants import SumMethod
from .mol_data import mol_sum
from .utils import wordwrap


# ---------------------------------------------------------------------------
# Node mapping
# ---------------------------------------------------------------------------

def node_map(
    mol_data: Optional[pl.DataFrame],
    node_data: pl.DataFrame,
    node_types: str | list[str] = "gene",
    node_sum: SumMethod = "sum",
    entrez_gnodes: bool = True,
) -> Optional[pl.DataFrame]:
    """
    Map *mol_data* onto pathway nodes of the specified *node_types*.

    Parameters
    ----------
    mol_data:      DataFrame whose first column contains molecule IDs and
                   remaining columns contain numeric values.  Pass None to
                   produce a position-only result with NaN values.
    node_data:     Tidy node DataFrame produced by kgml_parser.node_info().
    node_types:    Node type string(s) to include (e.g. "gene", "compound").
    node_sum:      Aggregation method when multiple probes map to one node.
    entrez_gnodes: True when gene nodes use Entrez IDs (vs KEGG gene IDs).

    Returns a merged DataFrame of node positions and molecular values, or
    None when no nodes of the requested type exist.
    """
    if isinstance(node_types, str):
        node_types = [node_types]

    target_nodes = node_data.filter(pl.col("type").is_in(node_types))
    if target_nodes.is_empty():
        return None

    # Expand the space-separated "name" field into individual KEGG IDs
    exploded = (
        target_nodes
        .with_columns(pl.col("name").str.split(" ").alias("kegg_names"))
        .explode("kegg_names")
        .with_columns(
            # Strip species prefix (e.g. "hsa:1234" → "1234")
            pl.col("kegg_names").str.replace(r"^[a-z]+:", "", literal=False)
        )
    )

    if mol_data is None:
        # Return node layout only, with a placeholder NaN value column
        return (
            exploded
            .group_by("entry_id")
            .agg([
                pl.col("kegg_names").first(),
                pl.col("x").first(),
                pl.col("y").first(),
                pl.col("width").first(),
                pl.col("height").first(),
                pl.col("label").first(),
                pl.col("type").first(),
                pl.col("size").first(),
            ])
            .with_columns(pl.lit(float("nan")).alias("mol_val"))
        )

    id_col = mol_data.columns[0]
    id_map = (
        exploded
        .select(["kegg_names", "entry_id"])
        .rename({"kegg_names": id_col, "entry_id": "__target"})
    )

    try:
        summed = mol_sum(mol_data, id_map.rename({"__target": "target_id"}).rename({"target_id": "__target"}), sum_method=node_sum)
    except ValueError:
        return None

    # Re-join aggregated values back to the full node layout
    plot_data = target_nodes.join(
        summed.rename({id_col: "entry_id"}),
        on="entry_id",
        how="left",
    )
    return plot_data
