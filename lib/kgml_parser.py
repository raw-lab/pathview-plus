"""
kgml_parser.py
Parse KEGG KGML (XML) pathway files into Python dataclasses and a tidy
Polars DataFrame suitable for downstream rendering.

Public API
----------
  parse_kgml  : Path → KGMLPathway
  node_info   : KGMLPathway → pl.DataFrame
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
class KGMLNode:
    """One <entry> element from a KGML file."""
    entry_id:  str
    name:      str
    node_type: str
    link:      str
    reaction:  str
    x:         Optional[float] = None
    y:         Optional[float] = None
    width:     Optional[float] = None
    height:    Optional[float] = None
    bgcolor:   str = "#FFFFFF"
    label:     str = ""
    shape:     str = "rectangle"
    component: list[str] = field(default_factory=list)


@dataclass
class KGMLEdge:
    """One <relation> element from a KGML file."""
    entry1:    str
    entry2:    str
    edge_type: str
    subtypes:  list[tuple[str, str]] = field(default_factory=list)


@dataclass
class KGMLReaction:
    """One <reaction> element from a KGML file."""
    name:       str
    rxn_type:   str                     # "reversible" | "irreversible"
    substrates: list[str] = field(default_factory=list)
    products:   list[str] = field(default_factory=list)


@dataclass
class KGMLPathway:
    """Container for all parsed elements of a KGML pathway file."""
    pathway_id:   str
    pathway_name: str
    nodes:        dict[str, KGMLNode]  = field(default_factory=dict)
    edges:        list[KGMLEdge]       = field(default_factory=list)
    reactions:    list[KGMLReaction]   = field(default_factory=list)


# ---------------------------------------------------------------------------
# Element parsers  (private helpers)
# ---------------------------------------------------------------------------

def _parse_graphics(elem: ET.Element) -> dict:
    """Extract display attributes from a <graphics> child element."""
    a = elem.attrib
    return {
        "x":      float(a.get("x", 0)),
        "y":      float(a.get("y", 0)),
        "width":  float(a.get("width", 46)),
        "height": float(a.get("height", 17)),
        "bgcolor": a.get("bgcolor", "#FFFFFF"),
        "shape":   a.get("type", "rectangle"),
        "label":   a.get("name", ""),
    }


def _parse_entry(elem: ET.Element) -> KGMLNode:
    """Parse a single <entry> element."""
    gfx_elem = elem.find("graphics")
    gfx = _parse_graphics(gfx_elem) if gfx_elem is not None else {}

    return KGMLNode(
        entry_id  = elem.attrib["id"],
        name      = elem.attrib.get("name", ""),
        node_type = elem.attrib.get("type", "gene"),
        link      = elem.attrib.get("link", ""),
        reaction  = elem.attrib.get("reaction", ""),
        x       = gfx.get("x"),
        y       = gfx.get("y"),
        width   = gfx.get("width"),
        height  = gfx.get("height"),
        bgcolor = gfx.get("bgcolor", "#FFFFFF"),
        label   = gfx.get("label", elem.attrib.get("name", "")),
        shape   = gfx.get("shape", "rectangle"),
        component = [c.attrib["id"] for c in elem.findall("component")],
    )


def _parse_relation(elem: ET.Element) -> KGMLEdge:
    """Parse a single <relation> element."""
    return KGMLEdge(
        entry1    = elem.attrib["entry1"],
        entry2    = elem.attrib["entry2"],
        edge_type = elem.attrib.get("type", ""),
        subtypes  = [
            (s.attrib.get("name", ""), s.attrib.get("value", ""))
            for s in elem.findall("subtype")
        ],
    )


def _parse_reaction(elem: ET.Element) -> KGMLReaction:
    """Parse a single <reaction> element."""
    return KGMLReaction(
        name       = elem.attrib.get("name", ""),
        rxn_type   = elem.attrib.get("type", "irreversible"),
        substrates = [s.attrib["id"] for s in elem.findall("substrate")],
        products   = [p.attrib["id"] for p in elem.findall("product")],
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_kgml(filepath: str | Path) -> KGMLPathway:
    """
    Parse a KEGG KGML file and return a populated KGMLPathway.

    Parameters
    ----------
    filepath: Path to the .xml KGML file.
    """
    root = ET.parse(filepath).getroot()
    pathway = KGMLPathway(
        pathway_id   = root.attrib.get("number", ""),
        pathway_name = root.attrib.get("name", ""),
    )
    _dispatch = {
        "entry":    lambda e: pathway.nodes.update({(n := _parse_entry(e)).entry_id: n}),
        "relation": lambda e: pathway.edges.append(_parse_relation(e)),
        "reaction": lambda e: pathway.reactions.append(_parse_reaction(e)),
    }
    for child in root:
        if child.tag in _dispatch:
            _dispatch[child.tag](child)

    return pathway


def node_info(pathway: KGMLPathway) -> pl.DataFrame:
    """
    Flatten KGMLPathway nodes into a tidy Polars DataFrame.

    Columns: entry_id, name, type, x, y, width, height, bgcolor,
             label, shape, reaction, component, size.
    """
    records = [
        {
            "entry_id":  n.entry_id,
            "name":      n.name,
            "type":      n.node_type,
            "x":         n.x,
            "y":         n.y,
            "width":     n.width,
            "height":    n.height,
            "bgcolor":   n.bgcolor,
            "label":     n.label,
            "shape":     n.shape,
            "reaction":  n.reaction,
            "component": ";".join(n.component),
            "size":      max(1, len(n.component)),
        }
        for n in pathway.nodes.values()
    ]
    return pl.DataFrame(records)
