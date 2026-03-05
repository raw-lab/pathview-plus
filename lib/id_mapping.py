"""
id_mapping.py
Gene and compound identifier mapping:
  - id2eg        : arbitrary gene ID  → Entrez Gene ID  (MyGene.info)
  - eg2id        : Entrez Gene ID     → any gene ID     (MyGene.info)
  - cpd_id_map   : compound ID        → KEGG compound   (KEGG REST conv)
"""

from __future__ import annotations

import warnings
from typing import Optional

import polars as pl
import requests

from .constants import KEGG_BASE


# ---------------------------------------------------------------------------
# Gene ID mapping  (MyGene.info REST API)
# ---------------------------------------------------------------------------

_SCOPE_MAP = {
    "symbol":  "symbol",
    "alias":   "alias",
    "uniprot": "uniprot",
    "ensembl": "ensembl.gene",
    "refseq":  "refseq",
}

_MYGENE_URL = "https://mygene.info/v3/querymany"


def _query_mygene(
    ids: list[str],
    scopes: str,
    field: str,
    species: str,
) -> dict[str, Optional[str]]:
    """
    POST a batch query to MyGene.info and return a {query_id: field_value} dict.
    Returns an empty-value dict on network failure.
    """
    payload = {
        "q":         ",".join(ids),
        "scopes":    scopes,
        "species":   species,
        "fields":    field,
        "returnall": "true",
    }
    try:
        resp = requests.post(_MYGENE_URL, data=payload, timeout=30)
        resp.raise_for_status()
        hits = resp.json()
    except Exception as exc:
        warnings.warn(f"MyGene.info query failed: {exc}")
        return {i: None for i in ids}

    lookup: dict[str, Optional[str]] = {}
    for hit in hits:
        qid = hit.get("query", "")
        val = hit.get(field)
        if isinstance(val, list):
            val = val[0] if val else None
        if qid and qid not in lookup:
            lookup[qid] = str(val) if val is not None else None

    return lookup


def id2eg(ids: list[str], category: str, org: str = "Hs") -> pl.DataFrame:
    """
    Map arbitrary gene IDs to Entrez Gene IDs.

    Parameters
    ----------
    ids:      Input gene identifiers.
    category: ID type of *ids* (e.g. "SYMBOL", "ENSEMBL", "UNIPROT").
    org:      Species for MyGene.info (e.g. "Hs", "Mm", "hsa").

    Returns a two-column DataFrame: [category, "ENTREZID"].

    Raises ValueError if *category* is already an Entrez type.
    """
    if category.lower() in ("entrez", "eg", "entrezid"):
        raise ValueError("Input IDs are already Entrez Gene IDs.")

    scope = _SCOPE_MAP.get(category.lower(), category.lower())
    lookup = _query_mygene(ids, scopes=scope, field="entrezgene", species=org)
    return pl.DataFrame({category: ids, "ENTREZID": [lookup.get(i) for i in ids]})


def eg2id(
    eg_ids: list[str],
    category: str = "SYMBOL",
    org: str = "Hs",
) -> pl.DataFrame:
    """
    Map Entrez Gene IDs to another identifier type.

    Parameters
    ----------
    eg_ids:   Entrez Gene IDs to convert.
    category: Target ID type (e.g. "SYMBOL", "UNIPROT", "ENSEMBL").
    org:      Species for MyGene.info.

    Returns a two-column DataFrame: ["ENTREZID", category].

    Raises ValueError if *category* is an Entrez type.
    """
    if category.lower() in ("entrez", "eg", "entrezid"):
        raise ValueError("Output category cannot be Entrez Gene ID.")

    field_map = {
        "symbol":  "symbol",
        "name":    "name",
        "uniprot": "uniprot",
        "ensembl": "ensembl.gene",
        "alias":   "alias",
    }
    field = field_map.get(category.lower(), category.lower())
    lookup = _query_mygene(eg_ids, scopes="entrezgene", field=field, species=org)
    return pl.DataFrame({"ENTREZID": eg_ids, category: [lookup.get(i) for i in eg_ids]})


# ---------------------------------------------------------------------------
# Compound ID mapping  (KEGG REST conv endpoint)
# ---------------------------------------------------------------------------

_CPD_TYPE_MAP = {
    "pubchem": "pubchem",
    "chebi":   "chebi",
    "kegg":    "cpd",
}


def cpd_id_map(
    in_ids: list[str],
    in_type: str,
    out_type: str = "KEGG",
) -> pl.DataFrame:
    """
    Map compound IDs between identifier systems using KEGG REST.

    Parameters
    ----------
    in_ids:   Input compound identifiers.
    in_type:  Source ID type (e.g. "PUBCHEM", "CHEBI", "KEGG").
    out_type: Target ID type (default "KEGG").

    Returns a two-column DataFrame: [in_type, out_type].
    """
    src = _CPD_TYPE_MAP.get(in_type.lower(),  in_type.lower())
    dst = _CPD_TYPE_MAP.get(out_type.lower(), out_type.lower())

    out_ids: list[Optional[str]] = []
    for cid in in_ids:
        url = f"{KEGG_BASE}/conv/{dst}/{src}:{cid}"
        try:
            resp = requests.get(url, timeout=15)
            if resp.ok and resp.text.strip():
                parts = resp.text.strip().split("\t")
                out_ids.append(parts[1].split(":")[1] if len(parts) > 1 else None)
            else:
                out_ids.append(None)
        except Exception:
            out_ids.append(None)

    return pl.DataFrame({in_type: in_ids, out_type: out_ids})
