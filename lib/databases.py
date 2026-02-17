"""
databases.py
Download SBGN-ML files from multiple pathway databases:
  - Reactome (human pathways)
  - MetaCyc (metabolic pathways)
  - PANTHER (protein pathways)
  - SMPDB (small molecule pathways)

Public API
----------
  download_reactome  : Download Reactome pathway
  download_metacyc   : Download MetaCyc pathway
  list_pathways      : List available pathways from a database
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional

import requests


# ---------------------------------------------------------------------------
# Reactome downloader
# ---------------------------------------------------------------------------

_REACTOME_BASE = "https://reactome.org/ContentService/exporter/sbgn"

def download_reactome(
    pathway_id: str,
    output_dir: Path = Path("."),
    species: str = "Homo sapiens",
) -> Optional[Path]:
    """
    Download a Reactome pathway in SBGN-ML format.
    
    Parameters
    ----------
    pathway_id:  Reactome stable ID (e.g., "R-HSA-109582" for Hemostasis)
    output_dir:  Directory to save the .sbgn file
    species:     Species name (default: "Homo sapiens")
    
    Returns
    -------
    Path to downloaded file, or None if download failed
    
    Example
    -------
    >>> path = download_reactome("R-HSA-109582", output_dir=Path("./pathways"))
    >>> print(f"Downloaded to {path}")
    
    Note
    ----
    Reactome pathway IDs follow the format: R-[species code]-[number]
    - R-HSA-* : Homo sapiens
    - R-MMU-* : Mus musculus
    - R-RNO-* : Rattus norvegicus
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Reactome API endpoint
    url = f"{_REACTOME_BASE}/{pathway_id}.sbgn"
    
    output_path = output_dir / f"{pathway_id}.sbgn"
    
    print(f"Info: Downloading Reactome pathway {pathway_id}...")
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        output_path.write_text(resp.text, encoding="utf-8")
        print(f"Info: Downloaded → {output_path}")
        return output_path
    except Exception as exc:
        warnings.warn(f"Failed to download {pathway_id}: {exc}")
        return None


# ---------------------------------------------------------------------------
# MetaCyc downloader
# ---------------------------------------------------------------------------

_METACYC_BASE = "https://biocyc.org/META/pathway"

def download_metacyc(
    pathway_id: str,
    output_dir: Path = Path("."),
) -> Optional[Path]:
    """
    Download a MetaCyc pathway in SBGN-ML format.
    
    Parameters
    ----------
    pathway_id:  MetaCyc pathway ID (e.g., "PWY-7210" for pyrimidine deoxyribonucleotides biosynthesis)
    output_dir:  Directory to save the .sbgn file
    
    Returns
    -------
    Path to downloaded file, or None if download failed
    
    Example
    -------
    >>> path = download_metacyc("PWY-7210", output_dir=Path("./pathways"))
    
    Note
    ----
    MetaCyc requires registration for API access. This function uses the
    public web interface and may not work for all pathways.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Try BioCyc SBGN export (may require authentication)
    url = f"{_METACYC_BASE}?id={pathway_id}&export=sbgn"
    output_path = output_dir / f"{pathway_id}.sbgn"
    
    print(f"Info: Downloading MetaCyc pathway {pathway_id}...")
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        output_path.write_text(resp.text, encoding="utf-8")
        print(f"Info: Downloaded → {output_path}")
        return output_path
    except Exception as exc:
        warnings.warn(
            f"Failed to download {pathway_id}: {exc}\n"
            "Note: MetaCyc may require authentication for some pathways."
        )
        return None


# ---------------------------------------------------------------------------
# PANTHER downloader
# ---------------------------------------------------------------------------

def download_panther(
    pathway_id: str,
    output_dir: Path = Path("."),
) -> Optional[Path]:
    """
    Download a PANTHER pathway in SBGN-ML format.
    
    Parameters
    ----------
    pathway_id:  PANTHER pathway ID (e.g., "P00001" for p53 pathway)
    output_dir:  Directory to save the .sbgn file
    
    Returns
    -------
    Path to downloaded file, or None if download failed
    
    Note
    ----
    PANTHER pathways use pre-generated SBGN-ML files.
    This function expects them to be hosted or provided locally.
    """
    warnings.warn(
        "PANTHER SBGN downloads not yet implemented. "
        "Please download SBGN-ML files manually from PANTHER website."
    )
    return None


# ---------------------------------------------------------------------------
# SMPDB downloader
# ---------------------------------------------------------------------------

_SMPDB_BASE = "https://smpdb.ca/pathways"

def download_smpdb(
    pathway_id: str,
    output_dir: Path = Path("."),
) -> Optional[Path]:
    """
    Download an SMPDB (Small Molecule Pathway Database) pathway.
    
    Parameters
    ----------
    pathway_id:  SMPDB pathway ID (e.g., "SMP0000001" for Glycolysis)
    output_dir:  Directory to save the .sbgn file
    
    Returns
    -------
    Path to downloaded file, or None if download failed
    
    Note
    ----
    SMPDB provides downloadable pathway files. This function may need
    adjustment based on current SMPDB API availability.
    """
    warnings.warn(
        "SMPDB SBGN downloads not yet fully implemented. "
        "Check https://smpdb.ca for pathway files."
    )
    return None


# ---------------------------------------------------------------------------
# Pathway listing
# ---------------------------------------------------------------------------

def list_reactome_pathways(species: str = "Homo sapiens") -> list[dict]:
    """
    List available Reactome pathways for a species.
    
    Parameters
    ----------
    species: Species name (e.g., "Homo sapiens", "Mus musculus")
    
    Returns
    -------
    List of dicts with keys: id, name, species
    
    Example
    -------
    >>> pathways = list_reactome_pathways("Homo sapiens")
    >>> for pw in pathways[:5]:
    ...     print(f"{pw['id']}: {pw['name']}")
    """
    url = "https://reactome.org/ContentService/data/pathways/top/Homo%20sapiens"
    
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        
        pathways = []
        for item in data:
            pathways.append({
                "id": item.get("stId", ""),
                "name": item.get("displayName", ""),
                "species": species,
            })
        return pathways
    except Exception as exc:
        warnings.warn(f"Failed to list Reactome pathways: {exc}")
        return []


# ---------------------------------------------------------------------------
# Database information
# ---------------------------------------------------------------------------

DATABASE_INFO = {
    "reactome": {
        "name": "Reactome",
        "description": "Curated human pathway database",
        "url": "https://reactome.org",
        "id_pattern": "R-[SPECIES]-[NUMBER]",
        "example": "R-HSA-109582",
        "downloader": download_reactome,
    },
    "metacyc": {
        "name": "MetaCyc",
        "description": "Metabolic pathway database",
        "url": "https://metacyc.org",
        "id_pattern": "PWY-[NUMBER]",
        "example": "PWY-7210",
        "downloader": download_metacyc,
    },
    "panther": {
        "name": "PANTHER",
        "description": "Protein analysis through evolutionary relationships",
        "url": "http://www.pantherdb.org",
        "id_pattern": "P[NUMBER]",
        "example": "P00001",
        "downloader": download_panther,
    },
    "smpdb": {
        "name": "SMPDB",
        "description": "Small Molecule Pathway Database",
        "url": "https://smpdb.ca",
        "id_pattern": "SMP[NUMBER]",
        "example": "SMP0000001",
        "downloader": download_smpdb,
    },
}


def detect_database(pathway_id: str) -> Optional[str]:
    """
    Detect which database a pathway ID belongs to.
    
    Parameters
    ----------
    pathway_id: Pathway identifier
    
    Returns
    -------
    Database name ("reactome", "metacyc", etc.) or None
    
    Example
    -------
    >>> detect_database("R-HSA-109582")
    'reactome'
    >>> detect_database("PWY-7210")
    'metacyc'
    """
    if pathway_id.startswith("R-") and "-" in pathway_id[2:]:
        return "reactome"
    elif pathway_id.startswith("PWY-"):
        return "metacyc"
    elif pathway_id.startswith("P") and pathway_id[1:].isdigit():
        return "panther"
    elif pathway_id.startswith("SMP"):
        return "smpdb"
    return None
