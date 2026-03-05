"""
kegg_api.py
KEGG REST API interactions:
  - SpeciesInfo      : dataclass holding per-species KEGG metadata
  - kegg_species_code: resolve a species name / abbreviation to SpeciesInfo
  - download_kegg    : fetch KGML (xml) and/or pathway image (png) files
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import requests

from .constants import KEGG_BASE


# ---------------------------------------------------------------------------
# Species resolution
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SpeciesInfo:
    """Immutable container for KEGG species metadata."""
    kegg_code:      str
    entrez_gnodes:  bool
    kegg_geneid:    Optional[str]
    ncbi_geneid:    Optional[str]
    ncbi_proteinid: Optional[str]
    uniprot:        Optional[str]


_KO_SPECIES = SpeciesInfo(
    kegg_code="ko",
    entrez_gnodes=False,
    kegg_geneid="K01488",
    ncbi_geneid=None,
    ncbi_proteinid=None,
    uniprot=None,
)


def kegg_species_code(species: str = "hsa") -> SpeciesInfo:
    """
    Resolve *species* (KEGG code, common name, or taxon) to a SpeciesInfo.

    Queries ``rest.kegg.jp/list/organism`` and matches any column.
    Raises ValueError for unknown species.
    """
    if species == "ko":
        return _KO_SPECIES

    url = f"{KEGG_BASE}/list/organism"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
    except Exception as exc:
        raise RuntimeError(f"Failed to fetch KEGG organism list: {exc}") from exc

    query = species.lower()
    for line in resp.text.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        if query in (p.lower() for p in parts):
            return SpeciesInfo(
                kegg_code=parts[1],
                entrez_gnodes=True,
                kegg_geneid=None,
                ncbi_geneid=None,
                ncbi_proteinid=None,
                uniprot=None,
            )

    raise ValueError(
        f"Unknown species '{species}'. "
        "Check https://rest.kegg.jp/list/organism for valid codes."
    )


# ---------------------------------------------------------------------------
# File download
# ---------------------------------------------------------------------------

def download_kegg(
    pathway_id: str,
    species: str = "hsa",
    kegg_dir: Path = Path("."),
    file_type: list[str] | None = None,
) -> dict[str, str]:
    """
    Download KEGG KGML and/or pathway image for *pathway_id*.

    Parameters
    ----------
    pathway_id: Numeric pathway ID, e.g. "04110" (species prefix is added
                automatically if absent).
    species:    KEGG species code used to build the full pathway ID.
    kegg_dir:   Directory where files are saved.
    file_type:  Subset of ["xml", "png"] to download (default: both).

    Returns a dict mapping the full pathway ID to "succeed" or "failed".
    """
    if file_type is None:
        file_type = ["xml", "png"]

    kegg_dir = Path(kegg_dir)
    kegg_dir.mkdir(parents=True, exist_ok=True)

    full_id = pathway_id if pathway_id.startswith(species) else f"{species}{pathway_id}"

    _url_templates = {
        "xml": f"{KEGG_BASE}/get/{full_id}/kgml",
        "png": f"{KEGG_BASE}/get/{full_id}/image",
    }
    _targets = {
        "xml": kegg_dir / f"{full_id}.xml",
        "png": kegg_dir / f"{full_id}.png",
    }

    status = {full_id: "succeed"}

    for ftype in file_type:
        url    = _url_templates[ftype]
        target = _targets[ftype]
        print(f"Info: Downloading {ftype} for {full_id} …")
        try:
            resp = requests.get(url, timeout=60)
            resp.raise_for_status()
            if ftype == "png":
                target.write_bytes(resp.content)
            else:
                target.write_text(resp.text, encoding="utf-8")
        except Exception as exc:
            warnings.warn(f"Download of {full_id} {ftype} failed: {exc}")
            status[full_id] = "failed"
            if target.exists():
                target.unlink()

    return status
