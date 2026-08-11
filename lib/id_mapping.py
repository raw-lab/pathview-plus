"""
id_mapping.py
Gene and compound identifier conversion.

Fixes over v2.x
---------------
* ``cpd_id_map`` issued **one HTTP request per compound**.  A 500-metabolite
  study meant 500 sequential round-trips; KEGG throttles long before that
  finishes.  Requests are now batched (KEGG's ``conv`` endpoint accepts up to
  100 ids per call) and the bundled cross-reference table answers the common
  conversions with no network at all.
* The species argument was passed straight through to MyGene.info, so
  ``org="hsa"`` — the KEGG code the orchestrator supplies — was not a species
  MyGene recognises and every lookup silently returned nulls.  KEGG codes are
  now translated to taxonomy ids via the bundled organism table.
* MyGene caps ``querymany`` at 1,000 ids per POST; v2.x sent the whole list in
  one request, so large inputs failed wholesale.  Requests are chunked.
* Failures were warnings that produced an all-null map indistinguishable from
  "these genes genuinely have no Entrez id".  Conversions now report how many
  IDs resolved.

Public API
----------
  id2eg       : any gene ID   -> Entrez
  eg2id       : Entrez        -> any gene ID
  cpd_id_map  : compound ID   -> KEGG compound (offline where possible)
  cpd_name_to_kegg : compound name -> KEGG accession (offline)
  supported_gene_idtypes / supported_cpd_idtypes
"""

from __future__ import annotations

import re
import warnings
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import polars as pl

from .bundled import read_bundled_tsv
from .cache import http_get, http_post_json, is_offline
from .constants import KEGG_BASE, MYGENE_URL
from .errors import NetworkError

_DATA = Path(__file__).parent / "data"

_MYGENE_CHUNK = 900
_KEGG_CHUNK = 100

#: Input id type -> MyGene.info scope.
_SCOPE_MAP = {
    "symbol": "symbol", "alias": "alias,symbol", "uniprot": "uniprot",
    "ensembl": "ensembl.gene", "ensemblprot": "ensembl.protein",
    "ensembltrans": "ensembl.transcript", "refseq": "refseq",
    "accnum": "accession", "unigene": "unigene", "entrez": "entrezgene",
    "entrezid": "entrezgene", "eg": "entrezgene", "genename": "name",
    "hgnc": "HGNC", "mgi": "MGI", "tair": "symbol", "orf": "symbol",
    "enzyme": "ec",
}

#: Output id type -> MyGene.info field.
_FIELD_MAP = {
    "symbol": "symbol", "genename": "name", "name": "name",
    "uniprot": "uniprot.Swiss-Prot", "ensembl": "ensembl.gene",
    "ensemblprot": "ensembl.protein", "ensembltrans": "ensembl.transcript",
    "refseq": "refseq.rna", "alias": "alias", "unigene": "unigene",
    "enzyme": "ec", "hgnc": "HGNC",
}

#: Compound id type -> the label used in the bundled cross-reference table.
_CPD_XREF_TYPES = {
    "cas": "CAS Registry Number",
    "chebi": "ChEBI accession",
    "chembl": "ChEMBL COMPOUND",
    "drugbank": "DrugBank accession",
    "hmdb": "HMDB accession",
    "lipidmaps": "LIPID MAPS accession",
    "pubchem": "PubChem accession",
    "knapsack": "KNApSAcK accession",
    "3dmet": "3DMET accession",
    "nikkaji": "NIKKAJI accession",
    "beilstein": "Beilstein Registry Number",
    "gmelin": "Gmelin Registry Number",
    "glycan": "KEGG GLYCAN accession",
    "drug": "KEGG DRUG accession",
}


def supported_gene_idtypes() -> list[str]:
    """Gene identifier types accepted by ``gene_idtype``."""
    return sorted({"ENTREZ", "KEGG"} | {k.upper() for k in _SCOPE_MAP}
                  | {k.upper() for k in _FIELD_MAP})


def supported_cpd_idtypes() -> list[str]:
    """Compound identifier types accepted by ``cpd_idtype``."""
    return sorted({"KEGG"} | {k.upper() for k in _CPD_XREF_TYPES})


@dataclass
class IdMapResult:
    """Conversion outcome with the counts needed to trust it."""

    data: pl.DataFrame
    n_input: int
    n_resolved: int
    source: str = "network"

    @property
    def resolved_fraction(self) -> float:
        return (self.n_resolved / self.n_input) if self.n_input else 0.0

    def summary(self) -> str:
        return (f"{self.n_resolved}/{self.n_input} identifiers resolved "
                f"({self.resolved_fraction:.0%}) via {self.source}")


# ---------------------------------------------------------------------------
# Species translation
# ---------------------------------------------------------------------------

def _mygene_species(org: str) -> str:
    """
    Translate a KEGG organism code to something MyGene.info understands.

    MyGene accepts common names and NCBI taxonomy ids, not KEGG codes; the
    bundled organism table carries the taxonomy id for exactly this purpose.
    """
    from .errors import SpeciesNotFoundError
    from .organisms import get_species_code

    raw = str(org or "").strip()
    if raw.isdigit():
        return raw
    try:
        info = get_species_code(raw)
    except SpeciesNotFoundError:
        return raw
    if info.tax_id:
        return info.tax_id
    return info.common_name or info.scientific_name or raw


def _chunks(items: Sequence[str], size: int) -> Iterable[list[str]]:
    for i in range(0, len(items), size):
        yield list(items[i:i + size])


# ---------------------------------------------------------------------------
# MyGene.info
# ---------------------------------------------------------------------------

def _query_mygene(
    ids: Sequence[str],
    scopes: str,
    field: str,
    species: str,
) -> dict[str, str | None]:
    """Batched MyGene.info querymany, chunked below the service's 1,000 cap."""
    lookup: dict[str, str | None] = {str(i): None for i in ids}
    unique = sorted({str(i).strip() for i in ids if str(i).strip()})
    if not unique:
        return lookup

    root_field = field.split(".")[0]

    for chunk in _chunks(unique, _MYGENE_CHUNK):
        payload = {
            "q": ",".join(chunk),
            "scopes": scopes,
            "species": species,
            "fields": field,
            "returnall": "false",
        }
        try:
            hits = http_post_json(MYGENE_URL, payload)
        except NetworkError as exc:
            warnings.warn(
                f"MyGene.info unavailable ({exc}); {len(chunk)} identifiers left "
                "unmapped. Supply gene_idtype='ENTREZ' to skip conversion.",
                stacklevel=3,
            )
            continue

        if isinstance(hits, dict):
            hits = hits.get("out", []) or []
        for hit in hits:
            if not isinstance(hit, dict) or hit.get("notfound"):
                continue
            qid = str(hit.get("query", ""))
            val = hit.get(field)
            if val is None:
                val = hit.get(root_field)
            if isinstance(val, dict):
                val = next((v for v in val.values() if v), None)
            if isinstance(val, list):
                val = next((v for v in val if v), None)
                if isinstance(val, dict):
                    val = next((v for v in val.values() if v), None)
            if qid and val is not None and lookup.get(qid) is None:
                lookup[qid] = str(val)

    return lookup


def id2eg(
    ids: Sequence[str],
    category: str = "SYMBOL",
    org: str = "hsa",
    detailed: bool = False,
) -> pl.DataFrame | IdMapResult:
    """
    Map arbitrary gene identifiers to Entrez Gene IDs.

    Returns a two-column DataFrame [category, "ENTREZID"].
    """
    cat = str(category).upper()
    if cat in ("ENTREZ", "EG", "ENTREZID"):
        raise ValueError(
            "id2eg: input identifiers are already Entrez; pass "
            "gene_idtype='ENTREZ' to skip conversion entirely."
        )

    ids = [str(i) for i in ids]
    scope = _SCOPE_MAP.get(cat.lower(), cat.lower())
    lookup = _query_mygene(ids, scope, "entrezgene", _mygene_species(org))
    values = [lookup.get(i) for i in ids]

    df = pl.DataFrame({cat: ids, "ENTREZID": values},
                      schema={cat: pl.String, "ENTREZID": pl.String})
    res = IdMapResult(df, len(ids), sum(v is not None for v in values),
                      "cache" if is_offline() else "MyGene.info")
    return res if detailed else df


def eg2id(
    eg_ids: Sequence[str],
    category: str = "SYMBOL",
    org: str = "hsa",
    detailed: bool = False,
) -> pl.DataFrame | IdMapResult:
    """
    Map Entrez Gene IDs to another identifier type.

    Returns a two-column DataFrame ["ENTREZID", category].
    """
    cat = str(category).upper()
    if cat in ("ENTREZ", "EG", "ENTREZID"):
        raise ValueError("eg2id: output category cannot be Entrez.")

    eg_ids = [str(i) for i in eg_ids]
    field = _FIELD_MAP.get(cat.lower(), cat.lower())
    lookup = _query_mygene(eg_ids, "entrezgene", field, _mygene_species(org))
    values = [lookup.get(i) for i in eg_ids]

    df = pl.DataFrame({"ENTREZID": eg_ids, cat: values},
                      schema={"ENTREZID": pl.String, cat: pl.String})
    res = IdMapResult(df, len(eg_ids), sum(v is not None for v in values),
                      "cache" if is_offline() else "MyGene.info")
    return res if detailed else df


# ---------------------------------------------------------------------------
# Compound identifiers
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def compound_xref() -> pl.DataFrame:
    """Bundled KEGG compound cross-reference table (kegg_id, id_type, accession)."""
    return read_bundled_tsv("cpd_xref.tsv.gz")


@lru_cache(maxsize=1)
def _name_to_kegg() -> dict[str, str]:
    from .mol_data import compound_names
    df = compound_names()
    # Synonyms and derived conjugate-base forms are both included, because
    # metabolomics platforms export "Pyruvate" where the reference table
    # holds "Pyruvic acid".  Keys are normalised so punctuation, Greek
    # letters and spacing do not defeat the lookup.
    out: dict[str, str] = {}
    for cid, nm in zip(df["compound_id"].to_list(), df["name"].to_list()):
        out.setdefault(normalize_compound_name(nm), cid)
        out.setdefault(str(nm).strip().lower(), cid)
    return out


_GREEK = {"\u03b1": "alpha", "\u03b2": "beta", "\u03b3": "gamma",
          "\u03b4": "delta", "\u03c9": "omega", "\u03b5": "epsilon"}


def normalize_compound_name(name: object) -> str:
    """
    Normalise a compound name for lookup.

    Lower-cases, expands Greek letters, and drops spaces, hyphens, commas and
    primes so that "alpha-D-Glucose", "a-D-glucose" and "alpha D glucose" all
    collide on one key.
    """
    s = str(name or "").strip().lower()
    for g, w in _GREEK.items():
        s = s.replace(g, w)
    return re.sub(r"[\s\-,'\u2019()\[\]]+", "", s)


def cpd_name_to_kegg(names: Sequence[str]) -> pl.DataFrame:
    """
    Map compound *names* to KEGG accessions using the bundled name table.

    Metabolomics platforms usually export names, not accessions, and this
    resolves them with no network access.
    """
    table = _name_to_kegg()
    names = [str(n) for n in names]
    resolved = [table.get(normalize_compound_name(n)) or table.get(n.strip().lower())
                for n in names]
    return pl.DataFrame({"NAME": names, "KEGG": resolved},
                        schema={"NAME": pl.String, "KEGG": pl.String})


def _offline_cpd_map(in_ids: Sequence[str], in_type: str) -> dict[str, str | None]:
    """Resolve compound ids through the bundled cross-reference table."""
    label = _CPD_XREF_TYPES.get(in_type.lower())
    if label is None:
        return {}
    xref = compound_xref().filter(pl.col("id_type") == label)
    if xref.is_empty():
        return {}
    table = dict(zip(xref["accession"].to_list(), xref["kegg_id"].to_list()))
    return {i: table.get(str(i).strip()) for i in in_ids}


def cpd_id_map(
    in_ids: Sequence[str],
    in_type: str = "KEGG",
    out_type: str = "KEGG",
    detailed: bool = False,
) -> pl.DataFrame | IdMapResult:
    """
    Map compound identifiers between systems.

    Resolution order: identity (KEGG -> KEGG), the bundled cross-reference
    table, then batched KEGG ``conv`` calls for anything still unresolved.
    Returns a two-column DataFrame [in_type, out_type].
    """
    in_ids = [str(i) for i in in_ids]
    src, dst = in_type.upper(), out_type.upper()

    if src == dst:
        df = pl.DataFrame({src: in_ids, f"{dst}_OUT": in_ids})
        res = IdMapResult(df.rename({f"{dst}_OUT": dst}) if src != dst else
                          pl.DataFrame({src: in_ids, dst: in_ids}),
                          len(in_ids), len(in_ids), "identity")
        return res if detailed else res.data

    if src == "NAME":
        df = cpd_name_to_kegg(in_ids).rename({"NAME": src, "KEGG": dst})
        n = df[dst].drop_nulls().len()
        res = IdMapResult(df, len(in_ids), n, "bundled name table")
        return res if detailed else df

    resolved = _offline_cpd_map(in_ids, src) if dst == "KEGG" else {}
    source = "bundled cross-reference"

    missing = [i for i in in_ids if not resolved.get(i)]
    if missing and not is_offline():
        source = "bundled cross-reference + KEGG conv"
        kegg_src = {"pubchem": "pubchem", "chebi": "chebi", "kegg": "cpd"}.get(
            src.lower(), src.lower())
        kegg_dst = {"kegg": "cpd", "pubchem": "pubchem", "chebi": "chebi"}.get(
            dst.lower(), dst.lower())
        for chunk in _chunks(missing, _KEGG_CHUNK):
            query = "+".join(f"{kegg_src}:{c}" for c in chunk)
            try:
                text = http_get(f"{KEGG_BASE}/conv/{kegg_dst}/{query}")
            except NetworkError as exc:
                warnings.warn(
                    f"KEGG conv unavailable ({exc}); {len(chunk)} compound ids "
                    "left unmapped.", stacklevel=3,
                )
                break
            for line in text.strip().splitlines():
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                a = parts[0].split(":", 1)[-1]
                b = parts[1].split(":", 1)[-1]
                if a in resolved or a in missing:
                    resolved[a] = b

    values = [resolved.get(i) for i in in_ids]
    df = pl.DataFrame({src: in_ids, dst: values},
                      schema={src: pl.String, dst: pl.String})
    res = IdMapResult(df, len(in_ids), sum(v is not None for v in values), source)
    return res if detailed else df
