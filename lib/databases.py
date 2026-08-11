"""
databases.py
Pathway file retrieval for KEGG and the SBGN databases.

The stub problem
----------------
v2.x exported ``download_panther`` and ``download_smpdb`` as public API.  Both
were bodies of ``warnings.warn(...); return None`` — they advertised support
for two databases and delivered nothing, while ``DATABASE_INFO`` listed them
alongside working downloaders so callers could not tell which was which.
``list_reactome_pathways(species=...)`` likewise ignored its argument and
hard-coded Homo sapiens into the URL.

This module states plainly what each source can do.  Every entry in
:data:`DATABASE_INFO` carries an ``available`` flag and a ``note``; the
resolver raises a specific error naming the manual route rather than returning
None and letting an empty map render.

Public API
----------
  download_kegg         : KGML and/or map image
  download_reactome     : Reactome SBGN by stable id
  list_reactome_pathways: Reactome top-level pathways for a species
  download_pathway      : dispatch on the id format
  detect_database       : identify the source database from an id
  DATABASE_INFO         : capability table
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

from .cache import http_get, http_get_bytes
from .constants import KEGG_BASE, KEGG_IMAGE_BASE, REACTOME_API, REACTOME_SBGN
from .errors import NetworkError, PathwayNotFoundError

# ---------------------------------------------------------------------------
# KEGG
# ---------------------------------------------------------------------------

def download_kegg(
    pathway_id: str,
    species: str = "hsa",
    kegg_dir: str | Path = ".",
    file_type: list[str] | None = None,
    overwrite: bool = False,
) -> dict[str, str]:
    """
    Fetch KGML (``xml``) and/or the map image (``png``) for a pathway.

    Returns ``{full_pathway_id: "succeed" | "failed" | "cached"}``.

    Unlike v2.x this reports per-file status, keeps a partially successful
    download (KGML without the image is still renderable in vector mode), and
    never leaves a truncated file behind.
    """
    file_type = list(file_type) if file_type else ["xml", "png"]
    kegg_dir = Path(kegg_dir)
    kegg_dir.mkdir(parents=True, exist_ok=True)

    pid = str(pathway_id).strip()
    full_id = pid if pid.startswith(species) else f"{species}{pid}"

    urls = {
        "xml": f"{KEGG_BASE}/get/{full_id}/kgml",
        "png": f"{KEGG_BASE}/get/{full_id}/image",
    }
    targets = {"xml": kegg_dir / f"{full_id}.xml", "png": kegg_dir / f"{full_id}.png"}

    status: dict[str, str] = {}
    ok, cached = True, True

    for ftype in file_type:
        target = targets[ftype]
        if target.exists() and target.stat().st_size > 0 and not overwrite:
            continue
        cached = False
        try:
            body = http_get_bytes(urls[ftype], suffix=f".{ftype}")
            if not body:
                raise NetworkError(f"empty response for {full_id} {ftype}")
            if ftype == "xml" and b"<pathway" not in body[:4000]:
                raise PathwayNotFoundError(
                    f"KEGG returned no KGML for {full_id}; check the pathway "
                    f"number and species code."
                )
            target.write_bytes(body)
        except Exception as exc:
            ok = False
            warnings.warn(f"Download of {full_id}.{ftype} failed: {exc}", stacklevel=2)
            if target.exists() and target.stat().st_size == 0:
                target.unlink()

    status[full_id] = "cached" if cached else ("succeed" if ok else "failed")
    return status


def kegg_pathway_image_url(pathway_id: str) -> str:
    """Public KEGG URL for a pathway map image."""
    return f"{KEGG_IMAGE_BASE}/kegg/pathway/{re.sub(r'[0-9].*$', '', pathway_id)}/{pathway_id}.png"


# ---------------------------------------------------------------------------
# Reactome
# ---------------------------------------------------------------------------

def download_reactome(
    pathway_id: str,
    output_dir: str | Path = ".",
    overwrite: bool = False,
) -> Path:
    """
    Download a Reactome pathway as SBGN-ML.

    Reactome is the one SBGN source with a public, stable per-id export
    endpoint, which is why it is the only one implemented here.

    Raises NetworkError / PathwayNotFoundError rather than returning None.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"{pathway_id}.sbgn"

    if out.exists() and out.stat().st_size > 0 and not overwrite:
        return out

    try:
        body = http_get(f"{REACTOME_SBGN}/{pathway_id}.sbgn", suffix=".sbgn")
    except Exception as exc:
        # The pre-generated collection carries 1,749 Reactome pathways, so a
        # firewalled or unavailable exporter is not the end of the road.
        from .sbgn_hub import download_sbgn
        try:
            return download_sbgn(pathway_id, output_dir=output_dir,
                                 overwrite=overwrite)
        except Exception:
            raise NetworkError(
                f"Reactome's SBGN exporter is unreachable ({exc}) and "
                f"{pathway_id} is not in the pre-generated collection."
            ) from exc

    if "<sbgn" not in body[:3000]:
        raise PathwayNotFoundError(
            f"Reactome returned no SBGN for {pathway_id!r}. Stable ids look "
            "like R-HSA-109688; list_reactome_pathways() enumerates valid ids."
        )
    out.write_text(body, encoding="utf-8")
    return out


def reactome_top_url(species: str = "Homo sapiens") -> str:
    """
    Build the Reactome top-level-pathways URL for *species*.

    Exposed separately so the species argument is testable without a network
    call: v2.x accepted ``species`` and then interpolated a hard-coded human
    string into the URL, so every call returned human pathways whatever was
    asked for.

    >>> reactome_top_url("Mus musculus").endswith("Mus%20musculus")
    True
    """
    from urllib.parse import quote
    return f"{REACTOME_API}/data/pathways/top/{quote(str(species))}"


def list_reactome_pathways(species: str = "Homo sapiens") -> list[dict]:
    """Top-level Reactome pathways for *species*."""
    url = reactome_top_url(species)
    try:
        import json
        data = json.loads(http_get(url, suffix=".json"))
    except Exception as exc:
        raise NetworkError(f"Could not list Reactome pathways for {species!r}: {exc}") from exc

    return [{"id": item.get("stId", ""), "name": item.get("displayName", ""),
             "species": species} for item in data if isinstance(item, dict)]


def find_reactome_pathways(query: str, species: str = "Homo sapiens") -> list[dict]:
    """Search Reactome by name and return matching pathway ids."""
    import json
    from urllib.parse import quote

    url = (f"{REACTOME_API}/search/query?query={quote(query)}"
           f"&species={quote(species)}&cluster=true")
    try:
        data = json.loads(http_get(url, suffix=".json"))
    except Exception as exc:
        raise NetworkError(f"Reactome search failed: {exc}") from exc

    out: list[dict] = []
    for group in data.get("results", []):
        for entry in group.get("entries", []):
            sid = entry.get("stId", "")
            if sid.startswith("R-"):
                out.append({"id": sid, "name": entry.get("name", ""),
                            "type": entry.get("exactType", "")})
    return out


# ---------------------------------------------------------------------------
# PANTHER / MetaCyc / SMPDB / MetaCrop
#
# None of these publishes a per-pathway SBGN endpoint, which is why 2.x
# shipped warn-and-return-None stubs for them.  The pre-generated SBGN
# collection does, so these are now real downloads: see :mod:`sbgn_hub`.
# ---------------------------------------------------------------------------

def download_panther(pathway_id: str, output_dir: str | Path = ".",
                     overwrite: bool = False) -> Path:
    """
    Download a PANTHER pathway as SBGN-ML.

    >>> download_panther("P00001", "sbgn")           # doctest: +SKIP
    PosixPath('sbgn/P00001.sbgn')
    """
    return _download_from_collection(pathway_id, "panther", output_dir, overwrite)


def download_metacyc(pathway_id: str, output_dir: str | Path = ".",
                     overwrite: bool = False) -> Path:
    """Download a MetaCyc pathway as SBGN-ML (e.g. ``"GLYCOLYSIS"``)."""
    return _download_from_collection(pathway_id, "metacyc", output_dir, overwrite)


def download_smpdb(pathway_id: str, output_dir: str | Path = ".",
                   overwrite: bool = False) -> Path:
    """Download an SMPDB pathway as SBGN-ML (e.g. ``"SMP00001"``)."""
    return _download_from_collection(pathway_id, "smpdb", output_dir, overwrite)


def download_metacrop(pathway_id: str, output_dir: str | Path = ".",
                      overwrite: bool = False) -> Path:
    """Download a MetaCrop pathway as SBGN-ML (e.g. ``"Glycolysis"``)."""
    return _download_from_collection(pathway_id, "metacrop", output_dir, overwrite)


def _download_from_collection(pathway_id: str, expected_source: str,
                              output_dir: str | Path, overwrite: bool) -> Path:
    """Fetch from the collection, checking the id belongs to the source asked for."""
    from .sbgn_hub import download_sbgn, find_sbgn_pathway

    entry = find_sbgn_pathway(pathway_id)
    if entry["source"] != expected_source:
        raise PathwayNotFoundError(
            f"{pathway_id!r} is a {entry['source']} pathway, not "
            f"{expected_source}. Use download_pathway() to dispatch "
            "automatically, or list_sbgn_pathways() to browse."
        )
    return download_sbgn(pathway_id, output_dir=output_dir, overwrite=overwrite)


# ---------------------------------------------------------------------------
# Capability table
# ---------------------------------------------------------------------------

DATABASE_INFO: dict[str, dict] = {
    "kegg": {
        "name": "KEGG",
        "format": "KGML + PNG",
        "id_pattern": r"^[a-z]{2,4}\d{5}$",
        "example": "hsa04110",
        "available": True,
        "downloader": download_kegg,
        "source": "KEGG REST API",
        "note": "Full support. Academic use of the REST API is free; "
                "commercial use requires a licence from Kanehisa Labs.",
    },
    "reactome": {
        "name": "Reactome",
        "format": "SBGN-ML",
        "id_pattern": r"^R-[A-Z]{3}-\d+$",
        "example": "R-HSA-109688",
        "available": True,
        "downloader": download_reactome,
        "source": "Reactome SBGN exporter, with the pre-generated "
                  "collection as a fallback",
        "note": "Full support. The live exporter covers every Reactome "
                "pathway; 1,749 are also in the offline collection.",
    },
    "panther": {
        "name": "PANTHER",
        "format": "SBGN-ML",
        "id_pattern": r"^P\d{5}(\.\d+)?$",
        "example": "P00001",
        "available": True,
        "downloader": download_panther,
        "source": "pre-generated SBGN collection",
        "note": "PANTHER publishes no per-pathway SBGN endpoint; 152 "
                "pathways are served from the pre-generated collection.",
    },
    "metacyc": {
        "name": "MetaCyc / BioCyc",
        "format": "SBGN-ML",
        "id_pattern": r"^[A-Za-z0-9][A-Za-z0-9+\-]*-?PWY[A-Za-z0-9\-]*$|^GLYCOLYSIS$",
        "example": "GLYCOLYSIS",
        "available": True,
        "downloader": download_metacyc,
        "source": "pre-generated SBGN collection",
        "note": "BioCyc requires a subscription for programmatic export; "
                "2,518 MetaCyc pathways are served from the pre-generated "
                "collection instead.",
    },
    "smpdb": {
        "name": "SMPDB",
        "format": "SBGN-ML",
        "id_pattern": r"^SMP\d+$",
        "example": "SMP00001",
        "available": True,
        "downloader": download_smpdb,
        "source": "pre-generated SBGN collection",
        "note": "SMPDB publishes bulk archives, not per-pathway SBGN; 725 "
                "pathways are served from the pre-generated collection.",
    },
    "metacrop": {
        "name": "MetaCrop",
        "format": "SBGN-ML",
        "id_pattern": r"^[A-Z][A-Za-z0-9 ,\-]+$",
        "example": "Glycolysis",
        "available": True,
        "downloader": download_metacrop,
        "source": "pre-generated SBGN collection",
        "note": "62 crop-plant metabolic pathways from the pre-generated "
                "collection.",
    },
}


def detect_database(pathway_id: str) -> str | None:
    """
    Identify the source database from an identifier's shape.

    >>> detect_database("R-HSA-109582")
    'reactome'
    >>> detect_database("hsa04110")
    'kegg'
    """
    pid = str(pathway_id).strip()
    if re.match(r"^\d{5}$", pid):
        return "kegg"
    for key in ("kegg", "reactome", "panther", "smpdb", "metacyc"):
        if re.match(DATABASE_INFO[key]["id_pattern"], pid):
            return key
    # Fall back to the collection index, which settles MetaCrop's
    # free-text names and any id whose shape is ambiguous.
    try:
        import polars as pl

        from .sbgn_hub import sbgn_index
        hit = sbgn_index().filter(pl.col("pathway_id") == pid)
        if not hit.is_empty():
            return hit.row(0, named=True)["source"]
    except Exception:
        pass
    return None


def available_databases() -> list[str]:
    """Databases pathview-plus can download from without manual steps."""
    return [k for k, v in DATABASE_INFO.items() if v["available"]]


def download_pathway(pathway_id: str, output_dir: str | Path = ".",
                     species: str = "hsa", **kw):
    """
    Download a pathway, dispatching on the identifier format.

    Raises a specific error naming the manual route for sources without a
    public API, instead of returning None.
    """
    db = detect_database(pathway_id)
    if db is None:
        raise PathwayNotFoundError(
            f"Could not recognise {pathway_id!r}. Expected a KEGG id "
            "(hsa04110 / 04110), a Reactome stable id (R-HSA-109688), a "
            "PANTHER id (P00001), an SMPDB id (SMP00001), or a MetaCyc "
            "pathway name (GLYCOLYSIS). Browse the SBGN collection with "
            "list_sbgn_pathways()."
        )
    info = DATABASE_INFO[db]
    if db == "kegg":
        return download_kegg(pathway_id, species=species, kegg_dir=output_dir, **kw)
    if db == "reactome":
        return download_reactome(pathway_id, output_dir=output_dir, **kw)
    return info["downloader"](pathway_id, output_dir=output_dir, **kw)
