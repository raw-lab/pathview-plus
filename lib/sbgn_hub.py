"""
sbgn_hub.py
The pre-generated SBGN collection and its cross-database identifier
crosswalks.

What this closes
----------------
SBGNview's advantage over a from-scratch reimplementation was never its
renderer — it was ``SBGNview.data``: thousands of pre-built SBGN-ML files for
databases with no usable per-pathway API, plus the identifier crosswalks
needed to join omics data to their opaque glyph ids.  Both are published on
GitHub, so both are usable from Python.

  * **Index** — 5,206 pathways across Reactome, SMPDB, PANTHER, MetaCyc and
    MetaCrop ship inside the wheel as a 32 KB table.  Listing and searching
    the collection therefore works offline.
  * **Files** — the SBGN-ML itself (about 690 MB in total) is fetched per
    pathway on first use and cached on disk, so nothing is downloaded that
    is not asked for.
  * **Crosswalks** — 770k identifier pairs (ChEBI, KEGG, Entrez, KO, gene
    symbol, compound name, Pathway Commons) ship in the wheel, so mapping
    data onto an SBGN map needs no network at all.

This is what makes ``download_panther``, ``download_metacyc`` and
``download_smpdb`` real functions rather than the warn-and-return-None stubs
they were in 2.x: those databases publish no per-pathway SBGN endpoint, but
the pre-generated collection does.

Public API
----------
  list_sbgn_pathways   : browse/search the collection (offline)
  find_sbgn_pathway    : resolve one id, with suggestions on a miss
  download_sbgn        : fetch and cache one pathway's SBGN-ML
  sbgn_collection_info : counts per source database
  sbgn_xref            : the bundled crosswalk table
  map_ids_to_sbgn      : user identifiers -> SBGN glyph identifiers
  crosswalk_routes     : which identifier conversions are available offline
"""

from __future__ import annotations

import gzip
from collections.abc import Iterable, Sequence
from functools import lru_cache
from pathlib import Path

import polars as pl

from .bundled import read_bundled_tsv
from .cache import cache_dir, http_get, is_offline
from .errors import NetworkError, PathwayNotFoundError

_DATA = Path(__file__).parent / "data"

#: Where the pre-generated collection lives.
SBGNHUB_RAW = ("https://raw.githubusercontent.com/datapplab/SBGNhub/master/"
               "data/SBGN.with.stamp")

#: Source database -> human-readable name.
SBGN_SOURCES: dict[str, str] = {
    "reactome": "Reactome",
    "smpdb": "SMPDB (Small Molecule Pathway Database)",
    "panther": "PANTHER",
    "metacyc": "MetaCyc",
    "metacrop": "MetaCrop",
}


# ---------------------------------------------------------------------------
# Index
# ---------------------------------------------------------------------------

def _read_gz_tsv(path: Path) -> pl.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return pl.read_csv(fh.read().encode(), separator="\t", infer_schema_length=0)


@lru_cache(maxsize=1)
def sbgn_index() -> pl.DataFrame:
    """
    The bundled catalogue of pre-generated SBGN pathways.

    Columns: pathway_id, source, subdir, filename.  Ships in the wheel, so
    this works with no network access.
    """
    return read_bundled_tsv("sbgn_index.tsv.gz")


def sbgn_collection_info() -> dict:
    """Counts per source database, plus the total."""
    idx = sbgn_index()
    counts = {
        row["source"]: row["len"]
        for row in idx.group_by("source").len().sort("source").to_dicts()
    }
    return {"total": idx.height, "by_source": counts, "sources": dict(SBGN_SOURCES)}


def list_sbgn_pathways(
    source: str | None = None,
    query: str | None = None,
    limit: int | None = None,
) -> pl.DataFrame:
    """
    Browse the collection, optionally filtered by *source* and/or *query*.

    >>> list_sbgn_pathways(source="panther", limit=3).height
    3
    """
    idx = sbgn_index()
    if source:
        key = str(source).lower()
        if key not in SBGN_SOURCES:
            raise ValueError(
                f"Unknown source {source!r}. Available: {', '.join(sorted(SBGN_SOURCES))}."
            )
        idx = idx.filter(pl.col("source") == key)
    if query:
        q = str(query).lower()
        idx = idx.filter(pl.col("pathway_id").str.to_lowercase().str.contains(q, literal=True))
    out = idx.select(["pathway_id", "source", "filename"])
    return out.head(limit) if limit else out


def find_sbgn_pathway(pathway_id: str) -> dict:
    """
    Resolve a pathway id to its catalogue entry.

    Raises PathwayNotFoundError with near-miss suggestions rather than
    returning None, so a typo is diagnosable.
    """
    pid = str(pathway_id).strip()
    hit = sbgn_index().filter(pl.col("pathway_id") == pid)
    if not hit.is_empty():
        return hit.row(0, named=True)

    ci = sbgn_index().filter(pl.col("pathway_id").str.to_lowercase() == pid.lower())
    if not ci.is_empty():
        return ci.row(0, named=True)

    import difflib
    close = difflib.get_close_matches(pid, sbgn_index()["pathway_id"].to_list(),
                                      n=5, cutoff=0.6)
    hint = f" Did you mean: {', '.join(close)}?" if close else ""
    raise PathwayNotFoundError(
        f"{pid!r} is not in the pre-generated SBGN collection "
        f"({sbgn_index().height} pathways).{hint} "
        "Use list_sbgn_pathways(source=...) to browse, or parse a local file "
        "with parse_sbgn()."
    )


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def sbgn_url(pathway_id: str) -> str:
    """Public URL of a pathway's SBGN-ML in the collection."""
    from urllib.parse import quote

    entry = find_sbgn_pathway(pathway_id)
    return f"{SBGNHUB_RAW}/{entry['subdir']}/{quote(entry['filename'])}"


def download_sbgn(
    pathway_id: str,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
) -> Path:
    """
    Fetch one pathway's SBGN-ML from the collection and cache it.

    Files land in *output_dir* when given, otherwise in the package cache, so
    a second call for the same pathway does no network I/O.

    Raises PathwayNotFoundError when the id is not in the catalogue, and
    NetworkError when it is but cannot be fetched.
    """
    entry = find_sbgn_pathway(pathway_id)
    pid = entry["pathway_id"]

    target_dir = Path(output_dir) if output_dir else (cache_dir() / "sbgn")
    target_dir.mkdir(parents=True, exist_ok=True)
    safe = "".join(c if c.isalnum() or c in "-._" else "_" for c in pid)
    out = target_dir / f"{safe}.sbgn"

    if out.exists() and out.stat().st_size > 0 and not overwrite:
        return out

    if is_offline():
        raise NetworkError(
            f"Offline mode and no cached copy of {pid}. The SBGN collection "
            "index is bundled, but the files themselves are fetched on "
            "demand; pre-fetch with download_sbgn() while online."
        )

    body = http_get(sbgn_url(pid), suffix=".sbgn")
    if "<sbgn" not in body[:4000]:
        raise NetworkError(
            f"The collection returned no SBGN-ML for {pid}; the response did "
            "not contain an <sbgn> element."
        )
    out.write_text(body, encoding="utf-8")
    return out


def download_sbgn_batch(
    pathway_ids: Sequence[str],
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    quiet: bool = True,
) -> dict[str, Path | str]:
    """
    Fetch several pathways, reporting per-pathway outcome.

    One failure does not abort the rest — the return value records which
    succeeded, so a partial download is usable and diagnosable.
    """
    out: dict[str, Path | str] = {}
    for pid in pathway_ids:
        try:
            out[str(pid)] = download_sbgn(pid, output_dir, overwrite)
        except Exception as exc:
            out[str(pid)] = f"failed: {exc}"
            if not quiet:
                print(f"[sbgn] {pid}: {exc}")
    return out


# ---------------------------------------------------------------------------
# Crosswalks
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def sbgn_xref() -> pl.DataFrame:
    """
    The bundled cross-database identifier table.

    Columns: from_type, from_id, to_type, to_id.  Around 770k pairs covering
    ChEBI, KEGG, Entrez, KO, gene symbol, compound name and Pathway Commons.
    """
    return read_bundled_tsv("sbgn_xref.tsv.gz")


def crosswalk_routes() -> pl.DataFrame:
    """Which identifier conversions the bundled crosswalk supports, and how many pairs."""
    return (sbgn_xref()
            .group_by(["from_type", "to_type"]).len()
            .rename({"len": "pairs"})
            .sort("pairs", descending=True))


@lru_cache(maxsize=1)
def _adjacency() -> dict[str, set[str]]:
    """Which identifier types the bundled crosswalk connects, both ways."""
    xr = sbgn_xref()
    adj: dict[str, set[str]] = {}
    for a, b in zip(xr["from_type"].to_list(), xr["to_type"].to_list()):
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    return adj


def id_route(src: str, dst: str) -> list[str] | None:
    """
    Shortest chain of identifier types joining *src* to *dst*.

    A fixed two-hop rule is not enough: Entrez reaches gene symbol only as
    Entrez -> KO -> Pathway Commons -> symbol.  Breadth-first search finds
    whatever chain the bundled crosswalk actually supports, and returns None
    rather than guessing when there is none.
    """
    src, dst = _canon(src), _canon(dst)
    if src == dst:
        return [src]
    adj = _adjacency()
    if src not in adj or dst not in adj:
        return None

    from collections import deque
    queue = deque([[src]])
    seen = {src}
    while queue:
        path = queue.popleft()
        for nxt in sorted(adj.get(path[-1], ())):
            if nxt == dst:
                return path + [nxt]
            if nxt not in seen:
                seen.add(nxt)
                queue.append(path + [nxt])
    return None


@lru_cache(maxsize=64)
def _route(from_type: str, to_type: str) -> dict[str, list[str]]:
    """Build a lookup for one conversion direction, using reverse pairs too."""
    xr = sbgn_xref()
    fwd = xr.filter((pl.col("from_type") == from_type) & (pl.col("to_type") == to_type))
    rev = xr.filter((pl.col("from_type") == to_type) & (pl.col("to_type") == from_type))

    table: dict[str, list[str]] = {}
    for src, dst in zip(fwd["from_id"].to_list(), fwd["to_id"].to_list()):
        table.setdefault(src, []).append(dst)
    for dst, src in zip(rev["from_id"].to_list(), rev["to_id"].to_list()):
        table.setdefault(src, []).append(dst)
    return table


#: Identifier type aliases accepted by :func:`map_ids_to_sbgn`.
_ID_ALIASES = {
    "entrez": "entrez", "entrezid": "entrez", "eg": "entrez",
    "symbol": "symbol", "genesymbol": "symbol", "hgnc": "symbol",
    "ko": "ko", "keggortholog": "ko",
    "kegg": "kegg", "keggcompound": "kegg", "cpd": "kegg",
    "chebi": "chebi",
    "name": "name", "compoundname": "name",
    "pathwaycommons": "pathwaycommons", "pc": "pathwaycommons",
    "ensembl": "ensembl",
}


def _canon(id_type: str) -> str:
    key = str(id_type or "").strip().lower().replace("_", "").replace(".", "")
    return _ID_ALIASES.get(key, key)


def supported_sbgn_idtypes() -> list[str]:
    """Identifier types the bundled crosswalk understands."""
    types = set(sbgn_xref()["from_type"].to_list()) | set(sbgn_xref()["to_type"].to_list())
    return sorted(types)


def map_ids_to_sbgn(
    ids: Iterable[str],
    id_type: str = "ENTREZ",
    target: str = "pathwayCommons",
    detailed: bool = False,
):
    """
    Map user identifiers onto the identifiers used inside SBGN glyphs.

    Pathway Commons SBGN glyph ids are opaque hashes, so without a crosswalk
    an SBGN map cannot carry omics data at all.  Direct routes are used where
    the bundled table has them, and a two-hop route through Pathway Commons or
    KEGG is tried otherwise.

    Returns a DataFrame [id_type, target]; with ``detailed=True`` an
    :class:`~pathview.id_mapping.IdMapResult` carrying the resolution counts.
    """
    from .id_mapping import IdMapResult

    src, dst = _canon(id_type), _canon(target)
    ids = [str(i).strip() for i in ids]

    if src == dst:
        # Identity. Column names must still differ, or the DataFrame collapses
        # to one column and every caller indexing columns[1] fails.
        df = pl.DataFrame({f"{src.upper()}_IN": ids, f"{dst.upper()}": ids},
                          schema={f"{src.upper()}_IN": pl.String,
                                  f"{dst.upper()}": pl.String})
        if not detailed:
            return df
        return IdMapResult(df, len(ids), len(ids), "identity")

    direct = _route(src, dst)
    values: list[str | None] = []
    route_used = "direct"

    if direct:
        for i in ids:
            hits = direct.get(i)
            values.append(hits[0] if hits else None)
    else:
        path = id_route(src, dst)
        if path is None:
            raise ValueError(
                f"No bundled route from {id_type!r} to {target!r}. "
                f"Available types: {', '.join(supported_sbgn_idtypes())}. "
                "See crosswalk_routes() for the pairs that exist."
            )
        route_used = " -> ".join(path)
        steps = [_route(a, b) for a, b in zip(path, path[1:])]
        for i in ids:
            current = [i]
            for step in steps:
                nxt: list[str] = []
                for c in current:
                    nxt.extend(step.get(c, ()))
                    if len(nxt) > 64:            # keep fan-out bounded
                        break
                current = list(dict.fromkeys(nxt))
                if not current:
                    break
            values.append(current[0] if current else None)

    df = pl.DataFrame({src.upper(): ids, dst.upper(): values},
                      schema={src.upper(): pl.String, dst.upper(): pl.String})
    if not detailed:
        return df
    return IdMapResult(df, len(ids), sum(v is not None for v in values),
                       f"bundled SBGN crosswalk ({route_used})")
