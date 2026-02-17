"""
sbgn_parser.py
Parse SBGN-ML (Systems Biology Graphical Notation ML) files.

SBGN is used by Reactome, MetaCyc, MetaCrop, PANTHER, and SMPDB.
Supports Process Description (PD), Entity Relationship (ER), and Activity Flow (AF) languages.

Public API
----------
  parse_sbgn     : Path → SBGNPathway
  sbgn_to_df     : SBGNPathway → pl.DataFrame (unified with KGML format)
  
SBGN vs KGML differences:
  - Glyphs (nodes) instead of entries
  - Arcs (edges) with Bezier splines
  - Compartments (cellular locations)
  - Clone markers for repeated entities
  - Process nodes (reactions, associations)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from xml.etree import ElementTree as ET

import polars as pl


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class SBGNGlyph:
    """One <glyph> element from SBGN-ML."""
    glyph_id: str
    glyph_class: str  # macromolecule, simple chemical, process, etc.
    label: str = ""
    x: Optional[float] = None
    y: Optional[float] = None
    width: Optional[float] = None
    height: Optional[float] = None
    compartment: Optional[str] = None
    clone_marker: bool = False
    state_variables: list[dict] = field(default_factory=list)
    unit_of_information: list[dict] = field(default_factory=list)


@dataclass
class SBGNArc:
    """One <arc> element from SBGN-ML."""
    arc_id: str
    arc_class: str  # production, consumption, catalysis, inhibition, etc.
    source: str
    target: str
    spline_points: list[tuple[float, float]] = field(default_factory=list)


@dataclass
class SBGNPathway:
    """Container for all parsed elements of an SBGN-ML file."""
    pathway_id: str
    pathway_name: str
    glyphs: dict[str, SBGNGlyph] = field(default_factory=dict)
    arcs: list[SBGNArc] = field(default_factory=list)
    compartments: dict[str, SBGNGlyph] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# SBGN glyph classes mapping
# ---------------------------------------------------------------------------

# Map SBGN glyph classes to simplified types for unified interface
_GLYPH_TYPE_MAP = {
    # Entity pool nodes (EPN)
    "macromolecule": "gene",
    "simple chemical": "compound",
    "nucleic acid feature": "gene",
    "complex": "gene",
    "multimer": "gene",
    "unspecified entity": "gene",
    
    # Process nodes (PN)
    "process": "process",
    "omitted process": "process",
    "uncertain process": "process",
    "association": "process",
    "dissociation": "process",
    "phenotype": "process",
    
    # Container nodes
    "compartment": "compartment",
    "submap": "map",
    
    # Logical operators
    "and": "operator",
    "or": "operator",
    "not": "operator",
}


# ---------------------------------------------------------------------------
# XML parsing helpers
# ---------------------------------------------------------------------------

def _parse_bbox(elem: ET.Element) -> dict:
    """Extract bounding box from <bbox> child element."""
    bbox = elem.find(".//bbox", namespaces={"": "http://sbgn.org/libsbgn/0.2"})
    if bbox is None:
        bbox = elem.find("bbox")  # Try without namespace
    
    if bbox is not None:
        return {
            "x": float(bbox.get("x", 0)) + float(bbox.get("w", 46)) / 2,
            "y": float(bbox.get("y", 0)) + float(bbox.get("h", 17)) / 2,
            "width": float(bbox.get("w", 46)),
            "height": float(bbox.get("h", 17)),
        }
    return {}


def _parse_label(elem: ET.Element) -> str:
    """Extract label text from <label> child element."""
    label = elem.find(".//label", namespaces={"": "http://sbgn.org/libsbgn/0.2"})
    if label is None:
        label = elem.find("label")
    
    if label is not None:
        return label.get("text", "")
    return ""


def _parse_state_variables(elem: ET.Element) -> list[dict]:
    """Parse state variable glyphs (child glyphs with class='state variable')."""
    states = []
    for child in elem.findall(".//glyph[@class='state variable']"):
        states.append({
            "variable": child.get("variable", ""),
            "value": child.get("value", ""),
            "label": _parse_label(child),
        })
    return states


def _parse_spline(arc_elem: ET.Element) -> list[tuple[float, float]]:
    """
    Parse spline curve from arc.
    
    SBGN-ML can have:
    1. Straight lines: <start> and <end> points
    2. Bezier curves: <start>, multiple <next>, <end> with control points
    """
    points = []
    
    # Start point
    start = arc_elem.find("start")
    if start is not None:
        points.append((float(start.get("x", 0)), float(start.get("y", 0))))
    
    # Intermediate points (could be Bezier control points)
    for pt in arc_elem.findall("next"):
        points.append((float(pt.get("x", 0)), float(pt.get("y", 0))))
    
    # End point
    end = arc_elem.find("end")
    if end is not None:
        points.append((float(end.get("x", 0)), float(end.get("y", 0))))
    
    return points


# ---------------------------------------------------------------------------
# Main parsing functions
# ---------------------------------------------------------------------------

def parse_sbgn(filepath: str | Path) -> SBGNPathway:
    """
    Parse an SBGN-ML file and return a populated SBGNPathway.
    
    Parameters
    ----------
    filepath: Path to the .sbgn or .xml SBGN-ML file
    
    Returns
    -------
    SBGNPathway object with all glyphs, arcs, and compartments
    
    Example
    -------
    >>> pathway = parse_sbgn("R-HSA-109582.sbgn")
    >>> print(f"Found {len(pathway.glyphs)} glyphs")
    >>> df = sbgn_to_df(pathway)
    """
    tree = ET.parse(filepath)
    root = tree.getroot()
    
    # Handle namespace (SBGN files often use xmlns)
    ns = {"sbgn": "http://sbgn.org/libsbgn/0.2"}
    map_elem = root.find(".//sbgn:map", ns)
    if map_elem is None:
        map_elem = root.find("map")
    if map_elem is None:
        map_elem = root  # Fall back to root
    
    pathway = SBGNPathway(
        pathway_id=map_elem.get("id", Path(filepath).stem),
        pathway_name=map_elem.get("language", "SBGN-PD"),
    )
    
    # Parse all glyphs
    for glyph_elem in map_elem.findall(".//glyph"):
        glyph_id = glyph_elem.get("id", "")
        glyph_class = glyph_elem.get("class", "")
        
        if not glyph_id:
            continue
        
        bbox = _parse_bbox(glyph_elem)
        label = _parse_label(glyph_elem)
        
        glyph = SBGNGlyph(
            glyph_id=glyph_id,
            glyph_class=glyph_class,
            label=label,
            x=bbox.get("x"),
            y=bbox.get("y"),
            width=bbox.get("width"),
            height=bbox.get("height"),
            compartment=glyph_elem.get("compartmentRef"),
            clone_marker="clone" in glyph_elem.attrib.get("clone", "").lower(),
            state_variables=_parse_state_variables(glyph_elem),
        )
        
        if glyph_class == "compartment":
            pathway.compartments[glyph_id] = glyph
        else:
            pathway.glyphs[glyph_id] = glyph
    
    # Parse all arcs
    for arc_elem in map_elem.findall(".//arc"):
        arc_id = arc_elem.get("id", "")
        arc_class = arc_elem.get("class", "")
        source = arc_elem.get("source", "")
        target = arc_elem.get("target", "")
        
        if not all([arc_id, source, target]):
            continue
        
        arc = SBGNArc(
            arc_id=arc_id,
            arc_class=arc_class,
            source=source,
            target=target,
            spline_points=_parse_spline(arc_elem),
        )
        pathway.arcs.append(arc)
    
    return pathway


def sbgn_to_df(pathway: SBGNPathway) -> pl.DataFrame:
    """
    Convert SBGN pathway to a unified Polars DataFrame (compatible with KGML format).
    
    Parameters
    ----------
    pathway: Parsed SBGNPathway object
    
    Returns
    -------
    DataFrame with columns matching KGML node_info format:
    entry_id, name, type, x, y, width, height, label, shape, etc.
    
    This allows SBGN pathways to use the same rendering pipeline as KEGG.
    """
    records = []
    
    for glyph_id, glyph in pathway.glyphs.items():
        # Map SBGN class to simplified type
        node_type = _GLYPH_TYPE_MAP.get(glyph.glyph_class, "unknown")
        
        # Determine shape
        shape_map = {
            "macromolecule": "roundedrectangle",
            "simple chemical": "ellipse",
            "complex": "octagon",
            "process": "square",
        }
        shape = shape_map.get(glyph.glyph_class, "rectangle")
        
        records.append({
            "entry_id": glyph_id,
            "name": glyph_id,  # SBGN IDs are typically database IDs
            "type": node_type,
            "x": glyph.x,
            "y": glyph.y,
            "width": glyph.width,
            "height": glyph.height,
            "bgcolor": "#FFFFFF",
            "label": glyph.label or glyph_id,
            "shape": shape,
            "reaction": "",
            "component": "",
            "size": 1,
            "kegg_names": glyph_id,  # For ID mapping
        })
    
    return pl.DataFrame(records) if records else pl.DataFrame()


# ---------------------------------------------------------------------------
# SBGN glyph class reference
# ---------------------------------------------------------------------------

SBGN_GLYPH_CLASSES = {
    # Entity Pool Nodes (EPN)
    "macromolecule": "Protein, gene product",
    "simple chemical": "Small molecule, metabolite",
    "nucleic acid feature": "DNA, RNA fragment",
    "complex": "Molecular complex",
    "multimer": "Homogeneous multimer",
    "unspecified entity": "Unknown entity type",
    
    # Process Nodes (PN)
    "process": "Biochemical process",
    "omitted process": "Process details omitted",
    "uncertain process": "Uncertain process",
    "association": "Complex formation",
    "dissociation": "Complex dissociation",
    "phenotype": "Observable phenotype",
    
    # Containers
    "compartment": "Cellular compartment",
    "submap": "Link to another map",
    
    # Logical operators
    "and": "Logical AND",
    "or": "Logical OR",
    "not": "Logical NOT",
}

SBGN_ARC_CLASSES = {
    "production": "Product of process",
    "consumption": "Consumed by process",
    "catalysis": "Catalyzes process",
    "modulation": "Modulates process",
    "stimulation": "Stimulates process",
    "inhibition": "Inhibits process",
    "necessary stimulation": "Required stimulator",
    "logic arc": "Logical operator connection",
}
