"""
mol_data.py
Molecular data aggregation and test-data generation.

Fixes over v2.x
---------------
* ``max_abs`` and ``random`` crashed.  Inside a Polars ``group_by.agg``
  context, ``pl.col(c).map_elements(f)`` hands *f* one scalar at a time, so
  the v2.x lambda calling ``s.to_numpy()`` raised
  ``AttributeError: 'float' object has no attribute 'to_numpy'``.  Two of the
  six documented ``node_sum`` methods were therefore dead on arrival.  They
  now use list-aggregation semantics and are covered by tests.
* The unmapped-ID count went negative whenever one source ID mapped to
  several targets (``mol_data.height - merged.height`` counts join expansion,
  not misses).  Unmapped IDs are now counted as a set difference and returned
  to the caller rather than printed.
* ``sim_mol_data`` needed a live KEGG call to invent gene IDs; it now draws
  from bundled tables and never touches the network.
* Real demo data (GSE16873) ships with the package, so examples exercise the
  same code path as user data.

Public API
----------
  mol_sum        : aggregate values from source IDs onto target IDs
  sim_mol_data   : simulated expression / abundance data
  demo_gene_data : real GSE16873 breast-cancer expression
  demo_cpd_data  : simulated metabolite abundances on real KEGG compound IDs
  compound_names : KEGG compound id -> name
"""

from __future__ import annotations

import gzip
import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import polars as pl

from .bundled import read_bundled_tsv
from .constants import SUM_METHODS, SumMethod
from .errors import MappingError

_DATA = Path(__file__).parent / "data"


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

@dataclass
class MolSumResult:
    """Aggregation outcome plus the diagnostics needed to trust it."""

    data: pl.DataFrame
    n_input: int
    n_mapped: int
    n_unmapped: int
    unmapped_ids: list[str]
    n_targets: int

    @property
    def mapped_fraction(self) -> float:
        return (self.n_mapped / self.n_input) if self.n_input else 0.0

    def summary(self) -> str:
        return (f"{self.n_mapped}/{self.n_input} IDs mapped "
                f"({self.mapped_fraction:.0%}) onto {self.n_targets} targets")


def _agg_expr(col: str, method: SumMethod, rng: np.random.Generator | None):
    """Polars aggregation expression for one column under *method*."""
    from .utils import max_abs, random_pick

    c = pl.col(col)
    if method == "sum":
        return c.sum()
    if method == "mean":
        return c.mean()
    if method == "median":
        return c.median()
    if method == "max":
        return c.max()
    if method == "min":
        return c.min()
    if method == "first":
        return c.first()
    if method == "max_abs":
        # Aggregate the group to a list first, then reduce it: map_elements in
        # an agg context is applied per-element, which is what broke v2.x.
        return c.map_batches(
            lambda s: pl.Series([max_abs(s.to_list())], dtype=pl.Float64),
            return_dtype=pl.Float64, returns_scalar=True,
        )
    if method == "random":
        return c.map_batches(
            lambda s: pl.Series([random_pick(s.to_list(), rng)], dtype=pl.Float64),
            return_dtype=pl.Float64, returns_scalar=True,
        )
    raise ValueError(
        f"Unknown sum_method {method!r}. Choose from: {', '.join(SUM_METHODS)}."
    )


def mol_sum(
    mol_data: pl.DataFrame,
    id_map: pl.DataFrame,
    sum_method: SumMethod = "sum",
    rand_seed: int | None = None,
    detailed: bool = False,
) -> pl.DataFrame | MolSumResult:
    """
    Aggregate *mol_data* from source identifiers onto the targets in *id_map*.

    Parameters
    ----------
    mol_data:   First column = source IDs, remaining columns = numeric values.
    id_map:     Two columns [source_id, target_id]; may be many-to-many.
    sum_method: One of sum, mean, median, max, min, max_abs, random, first.
    rand_seed:  Seed for the ``random`` method, so runs are reproducible.
    detailed:   Return a :class:`MolSumResult` instead of a bare DataFrame.

    Raises MappingError when nothing maps — with the diagnostics needed to see
    why, rather than v2.x's bare "No IDs could be mapped".
    """
    if sum_method not in SUM_METHODS:
        raise ValueError(
            f"Unknown sum_method {sum_method!r}. Choose from: {', '.join(SUM_METHODS)}."
        )
    if mol_data is None or mol_data.is_empty():
        raise MappingError("mol_sum: mol_data is empty.")
    if id_map is None or id_map.is_empty():
        raise MappingError("mol_sum: id_map is empty; nothing to map onto.")
    if mol_data.width < 2:
        raise MappingError(
            "mol_sum: mol_data needs an ID column plus at least one value column."
        )

    id_col = mol_data.columns[0]
    src_col, tgt_col = id_map.columns[:2]

    mol = mol_data.with_columns(pl.col(id_col).cast(pl.String).str.strip_chars())
    mapping = (
        id_map.select([pl.col(src_col).cast(pl.String).str.strip_chars().alias(id_col),
                       pl.col(tgt_col).cast(pl.String).alias("__target")])
        .drop_nulls()
        .unique()
    )

    numeric_cols = [c for c in mol.columns
                    if c != id_col and mol.schema[c].is_numeric()]
    if not numeric_cols:
        mol = mol.with_columns([
            pl.col(c).cast(pl.Float64, strict=False) for c in mol.columns if c != id_col
        ])
        numeric_cols = [c for c in mol.columns if c != id_col]

    merged = mol.join(mapping, on=id_col, how="inner")

    src_ids = set(mol[id_col].to_list())
    mapped_ids = set(merged[id_col].to_list()) if not merged.is_empty() else set()
    unmapped = sorted(src_ids - mapped_ids)

    if merged.is_empty():
        sample_in = list(src_ids)[:3]
        sample_map = mapping[id_col].to_list()[:3]
        raise MappingError(
            f"No identifier in column {id_col!r} matched the id_map. "
            f"Input looks like {sample_in}; the map expects {sample_map}. "
            "Check gene_idtype/cpd_idtype and the species code."
        )

    rng = np.random.default_rng(rand_seed) if rand_seed is not None else None
    aggregated = (
        merged.drop(id_col)
        .group_by("__target")
        .agg([_agg_expr(c, sum_method, rng).alias(c) for c in numeric_cols])
        .rename({"__target": id_col})
        .sort(id_col)
    )

    result = MolSumResult(
        data=aggregated,
        n_input=len(src_ids),
        n_mapped=len(mapped_ids),
        n_unmapped=len(unmapped),
        unmapped_ids=unmapped,
        n_targets=aggregated.height,
    )
    return result if detailed else aggregated


# ---------------------------------------------------------------------------
# Bundled reference tables
# ---------------------------------------------------------------------------

def _read_gz_tsv(path: Path) -> pl.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return pl.read_csv(fh.read().encode(), separator="\t", infer_schema_length=0)


@lru_cache(maxsize=1)
def compound_names() -> pl.DataFrame:
    """KEGG compound accession -> primary name (7,312 compounds, offline)."""
    return read_bundled_tsv("cpd_names.tsv.gz")


@lru_cache(maxsize=1)
def _compound_name_map() -> dict[str, str]:
    """Accession -> preferred display name (primary name wins over synonyms)."""
    df = compound_names().filter(pl.col("primary") == "1")
    out = dict(zip(df["compound_id"].to_list(), df["name"].to_list()))
    for cid, nm in zip(compound_names()["compound_id"].to_list(),
                       compound_names()["name"].to_list()):
        out.setdefault(cid, nm)
    return out


def compound_name(cid: str) -> str:
    """Human-readable name for a KEGG compound id, or the id itself."""
    return _compound_name_map().get(str(cid).strip(), str(cid))


@lru_cache(maxsize=1)
def _known_compound_ids() -> list[str]:
    return compound_names()["compound_id"].unique().sort().to_list()


# ---------------------------------------------------------------------------
# Demo / simulated data
# ---------------------------------------------------------------------------

def demo_gene_data(n_samples: int = 1, as_log_ratio: bool = True) -> pl.DataFrame:
    """
    Real expression data: GSE16873 breast-cancer DCIS samples.

    This is the dataset R pathview's vignette uses, shipped with the package
    so the examples exercise the real mapping path rather than random noise.

    Returns a DataFrame with an ``entrez`` column plus *n_samples* value
    columns of log2 ratios.
    """
    df = read_bundled_tsv("demo_gse16873.tsv.gz")
    value_cols = [c for c in df.columns if c != "entrez"][:max(1, int(n_samples))]
    out = df.select(
        [pl.col("entrez").cast(pl.String)]
        + [pl.col(c).cast(pl.Float64) for c in value_cols]
    )
    if not as_log_ratio:
        out = out.with_columns([pl.col(c) + 8.0 for c in value_cols])
    return out


def demo_cpd_data(
    pathway_compounds: Sequence[str] | None = None,
    n_mol: int = 30,
    n_exp: int = 1,
    rand_seed: int = 42,
) -> pl.DataFrame:
    """
    Metabolite abundances on real KEGG compound accessions.

    When *pathway_compounds* is given the IDs are drawn from that pathway, so
    the data actually lands on the map instead of scattering across KEGG.
    """
    rng = np.random.default_rng(rand_seed)
    pool = [c for c in (pathway_compounds or _known_compound_ids()) if c]
    if not pool:
        pool = _known_compound_ids()
    n = min(int(n_mol), len(pool))
    ids = list(rng.choice(np.array(pool, dtype=object), size=n, replace=False))

    data: dict[str, list] = {"compound": [str(i) for i in ids]}
    for i in range(1, max(1, int(n_exp)) + 1):
        data[f"exp{i}"] = (rng.standard_normal(n) * 1.2).round(4).tolist()
    return pl.DataFrame(data)


def sim_mol_data(
    mol_type: str = "gene",
    species: str = "hsa",
    n_mol: int = 100,
    n_exp: int = 1,
    rand_seed: int = 100,
    discrete: bool = False,
    id_pool: Sequence[str] | None = None,
) -> pl.DataFrame:
    """
    Generate simulated molecular data for demos and tests.

    Works fully offline: gene IDs are drawn from the bundled demo expression
    set (real Entrez IDs) and compound IDs from the bundled KEGG compound
    table (real accessions).  v2.x called KEGG for a gene list and, on
    failure, emitted IDs like ``gene1``…``gene1000`` that map to nothing, so
    "simulated" runs silently produced empty maps.
    """
    rng = np.random.default_rng(rand_seed)

    if id_pool is not None:
        pool = [str(x) for x in id_pool if str(x)]
    elif mol_type == "gene":
        pool = demo_gene_data()["entrez"].to_list()
    elif mol_type in ("cpd", "compound", "metabolite"):
        pool = _known_compound_ids()
    else:
        raise ValueError(f"mol_type must be 'gene' or 'cpd', got {mol_type!r}.")

    if not pool:
        raise MappingError(f"No identifier pool available for mol_type={mol_type!r}.")

    n = int(n_mol)
    if n > len(pool):
        warnings.warn(
            f"Requested {n} molecules but only {len(pool)} are available; using all.",
            stacklevel=2,
        )
        n = len(pool)

    ids = [str(x) for x in rng.choice(np.array(pool, dtype=object), size=n, replace=False)]
    id_name = "entrez" if mol_type == "gene" else "compound"

    if discrete:
        return pl.DataFrame({id_name: ids})

    data: dict[str, list] = {id_name: ids}
    for i in range(1, max(1, int(n_exp)) + 1):
        data[f"exp{i}"] = rng.standard_normal(n).round(4).tolist()
    return pl.DataFrame(data)
