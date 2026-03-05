"""
mol_data.py
Molecular data handling:
  - mol_sum      : aggregate multi-probe data to target IDs (Polars-based)
  - sim_mol_data : generate simulated expression / abundance data for testing
"""

from __future__ import annotations

import warnings
from typing import Callable, Optional

import numpy as np
import polars as pl
import requests

from .constants import KEGG_BASE, SumMethod
from .utils import max_abs, random_pick


# ---------------------------------------------------------------------------
# Aggregation dispatch
# ---------------------------------------------------------------------------

def _make_agg_expr(col: str, method: SumMethod):
    """Return a Polars aggregation expression for a single column."""
    match method:
        case "sum":    return pl.col(col).sum()
        case "mean":   return pl.col(col).mean()
        case "median": return pl.col(col).median()
        case "max":    return pl.col(col).max()
        case "max_abs":
            return pl.col(col).map_elements(
                lambda s: max_abs(s.to_numpy()), return_dtype=pl.Float64
            )
        case "random":
            return pl.col(col).map_elements(
                lambda s: random_pick(s.to_numpy()), return_dtype=pl.Float64
            )
        case _:
            raise ValueError(
                f"Unknown sum_method '{method}'. "
                "Choose from: sum, mean, median, max, max_abs, random."
            )


# ---------------------------------------------------------------------------
# mol_sum
# ---------------------------------------------------------------------------

def mol_sum(
    mol_data: pl.DataFrame,
    id_map: pl.DataFrame,
    sum_method: SumMethod = "sum",
) -> pl.DataFrame:
    """
    Aggregate *mol_data* from source IDs to target IDs defined by *id_map*.

    Parameters
    ----------
    mol_data:   DataFrame whose **first column** contains source IDs; all
                remaining columns are treated as numeric expression values.
    id_map:     Two-column DataFrame [source_id, target_id].
    sum_method: How to combine multiple source rows mapping to one target.

    Returns a DataFrame keyed by target IDs with the same numeric columns
    as *mol_data*.  Raises ValueError when no IDs can be mapped.
    """
    id_col     = mol_data.columns[0]
    src_col, tgt_col = id_map.columns[:2]

    # Rename id_map columns to neutral names for the join
    mapping = id_map.rename({src_col: id_col, tgt_col: "__target"})

    merged = mol_data.join(mapping, on=id_col, how="inner")
    if merged.is_empty():
        raise ValueError(
            f"No IDs from '{id_col}' could be mapped using the provided id_map."
        )

    n_unmapped = mol_data.height - merged.height
    if n_unmapped > 0:
        print(f"Note: {n_unmapped} of {mol_data.height} input IDs unmapped.")

    numeric_cols = [c for c in merged.columns if c not in (id_col, "__target")]
    aggregated = (
        merged
        .drop(id_col)
        .group_by("__target")
        .agg([_make_agg_expr(c, sum_method).alias(c) for c in numeric_cols])
        .rename({"__target": id_col})
    )
    return aggregated


# ---------------------------------------------------------------------------
# sim_mol_data
# ---------------------------------------------------------------------------

def sim_mol_data(
    mol_type: str = "gene",
    species: str = "hsa",
    n_mol: int = 100,
    n_exp: int = 1,
    rand_seed: int = 100,
    discrete: bool = False,
) -> pl.DataFrame:
    """
    Generate simulated molecular abundance data for testing and demos.

    Parameters
    ----------
    mol_type:  "gene" (fetches real KEGG gene IDs) or "cpd" (fake KEGG IDs).
    species:   KEGG species code used when *mol_type* is "gene".
    n_mol:     Number of molecules to sample.
    n_exp:     Number of simulated experiment columns.
    rand_seed: NumPy RNG seed for reproducibility.
    discrete:  When True, return only the sampled IDs (no numeric values).

    Returns a DataFrame with an 'id' column and *n_exp* numeric columns named
    'exp1', 'exp2', … (or just 'id' when *discrete* is True).
    """
    rng = np.random.default_rng(rand_seed)

    if mol_type == "gene":
        ids = _fetch_kegg_gene_ids(species)
    elif mol_type == "cpd":
        ids = [f"C{i:05d}" for i in range(1, 5001)]
    else:
        raise ValueError(f"mol_type must be 'gene' or 'cpd', got '{mol_type}'.")

    n_available = len(ids)
    if n_mol > n_available:
        warnings.warn(
            f"Requested {n_mol} molecules but only {n_available} available; "
            "using all available IDs."
        )
        n_mol = n_available

    sampled = list(rng.choice(ids, size=n_mol, replace=False))

    if discrete:
        return pl.DataFrame({"id": sampled})

    data: dict[str, list] = {"id": sampled}
    for i in range(1, n_exp + 1):
        data[f"exp{i}"] = rng.standard_normal(n_mol).tolist()

    return pl.DataFrame(data)


def _fetch_kegg_gene_ids(species: str) -> list[str]:
    """Fetch all gene IDs for *species* from KEGG; fall back to dummy IDs."""
    url = f"{KEGG_BASE}/list/{species}"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return [
            line.split("\t")[0].split(":")[1]
            for line in resp.text.strip().splitlines()
            if "\t" in line
        ]
    except Exception as exc:
        warnings.warn(f"Failed to fetch KEGG gene list for '{species}': {exc}. Using dummy IDs.")
        return [f"gene{i}" for i in range(1, 1001)]
