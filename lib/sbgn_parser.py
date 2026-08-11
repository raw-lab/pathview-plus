"""
sbgn_parser.py
Parse SBGN-ML (Process Description, Entity Relationship, Activity Flow).

Fixes over v2.x
---------------
* **Namespaces.** Real Reactome/SBGNview exports declare a default namespace
  (``xmlns="http://sbgn.org/libsbgn/pd/0.1"``), so every tag is
  ``{ns}glyph``.  v2.x searched for the bare tag ``.//glyph`` and returned
  **zero glyphs and zero arcs** for every namespaced file — i.e. for every
  real Reactome download.  It also passed ``namespaces={"": ...}``, which
  ElementTree rejects outright ("empty namespace prefix is not supported in
  ElementPath").  Matching is now done on local names, so namespaced and
  bare documents parse identically.
* **Port resolution.** SBGN arcs usually reference ``<port>`` elements rather
  than glyphs.  v2.x stored the raw port id as source/target, so almost no
  arc resolved to a real node.  Ports are now indexed to their owning glyph
  and arcs resolve through them.
* **Nested glyphs.** ``.//glyph`` also matches state variables and units of
  information, which v2.x promoted to top-level nodes, injecting dozens of
  junk "nodes" per map.  Auxiliary glyph classes are now classified, not
  emitted as entities.
* **Clone markers.** ``clone`` is a child *element*, not an attribute; the
  v2.x attribute test was always False.
* **Identifiers.** Database cross-references are pulled from ``<annotation>``
  RDF and from the id itself, so SBGN nodes can carry omics data instead of
  being keyed on an opaque hash.

Public API
----------
  parse_sbgn      : path -> SBGNPathway
  sbgn_to_df      : SBGNPathway -> node DataFrame (KGML-compatible schema)
  sbgn_edges      : SBGNPathway -> edge DataFrame
  sbgn_canvas     : SBGNPathway -> (x0, y0, width, height)
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET

import polars as pl

from .errors import ParseError, PathwayNotFoundError

# ---------------------------------------------------------------------------
# Class vocabularies
# ---------------------------------------------------------------------------

#: Glyph classes that are auxiliary decorations, not entities in their own
#: right.  These must never become nodes.
AUXILIARY_CLASSES = frozenset({
    "state variable", "unit of information", "entity", "terminal",
    "outcome", "existence", "location", "cardinality",
})

ENTITY_CLASSES = frozenset({
    "macromolecule", "simple chemical", "nucleic acid feature", "complex",
    "multimer", "unspecified entity", "perturbing agent", "source and sink",
    "macromolecule multimer", "complex multimer", "simple chemical multimer",
    "nucleic acid feature multimer", "biological activity",
})

PROCESS_CLASSES = frozenset({
    "process", "omitted process", "uncertain process", "association",
    "dissociation", "phenotype",
})

OPERATOR_CLASSES = frozenset({"and", "or", "not", "delay"})

#: SBGN glyph class -> the simplified node type used by the render pipeline.
GLYPH_TYPE_MAP: dict[str, str] = {
    **{c: "gene" for c in (
        "macromolecule", "macromolecule multimer", "nucleic acid feature",
        "nucleic acid feature multimer", "complex", "complex multimer",
        "multimer", "unspecified entity", "perturbing agent",
        "biological activity",
    )},
    **{c: "compound" for c in (
        "simple chemical", "simple chemical multimer", "source and sink",
    )},
    **{c: "process" for c in PROCESS_CLASSES},
    **{c: "operator" for c in OPERATOR_CLASSES},
    "compartment": "compartment",
    "submap": "map",
    "tag": "map",
}

SBGN_GLYPH_CLASSES: dict[str, str] = {
    "macromolecule": "Protein or gene product",
    "simple chemical": "Small molecule / metabolite",
    "nucleic acid feature": "DNA or RNA fragment",
    "complex": "Molecular complex",
    "multimer": "Homogeneous multimer",
    "unspecified entity": "Entity of unspecified type",
    "perturbing agent": "External perturbation",
    "source and sink": "Empty set / boundary",
    "process": "Biochemical process",
    "omitted process": "Process with omitted detail",
    "uncertain process": "Process of uncertain nature",
    "association": "Complex formation",
    "dissociation": "Complex dissociation",
    "phenotype": "Observable phenotype",
    "compartment": "Cellular compartment",
    "submap": "Link to another map",
    "and": "Logical AND", "or": "Logical OR", "not": "Logical NOT",
    "state variable": "Post-translational state",
    "unit of information": "Auxiliary annotation",
}

SBGN_ARC_CLASSES: dict[str, str] = {
    "production": "Product of a process",
    "consumption": "Consumed by a process",
    "catalysis": "Catalyses a process",
    "modulation": "Modulates a process",
    "stimulation": "Stimulates a process",
    "inhibition": "Inhibits a process",
    "necessary stimulation": "Required stimulator",
    "logic arc": "Connection to a logical operator",
    "equivalence arc": "Equivalence between entities",
    "interaction": "Undirected interaction",
    "absolute inhibition": "Complete inhibition",
    "absolute stimulation": "Complete stimulation",
}

#: Arc classes that should be drawn with an inhibitory (bar) head.
INHIBITORY_ARCS = frozenset({"inhibition", "absolute inhibition"})

_ID_RX = re.compile(
    r"(?:^|[_/:])(?P<db>uniprot|chebi|kegg|reactome|ensembl|hgnc|ncbigene|entrez)"
    r"[:_](?P<acc>[A-Za-z0-9\-\.]+)", re.I
)


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class SBGNGlyph:
    """One ``<glyph>`` element."""

    glyph_id: str
    glyph_class: str
    label: str = ""
    x: float | None = None
    y: float | None = None
    width: float | None = None
    height: float | None = None
    compartment: str | None = None
    clone_marker: bool = False
    state_variables: list[dict] = field(default_factory=list)
    units_of_information: list[dict] = field(default_factory=list)
    identifiers: list[str] = field(default_factory=list)
    parent: str | None = None

    @property
    def node_type(self) -> str:
        return GLYPH_TYPE_MAP.get(self.glyph_class, "unknown")

    @property
    def is_entity(self) -> bool:
        return self.glyph_class in ENTITY_CLASSES

    @property
    def is_process(self) -> bool:
        return self.glyph_class in PROCESS_CLASSES


@dataclass
class SBGNArc:
    """One ``<arc>`` element, with source/target resolved through ports."""

    arc_id: str
    arc_class: str
    source: str
    target: str
    spline_points: list[tuple[float, float]] = field(default_factory=list)
    raw_source: str = ""
    raw_target: str = ""

    @property
    def inhibitory(self) -> bool:
        return self.arc_class in INHIBITORY_ARCS


@dataclass
class SBGNPathway:
    """Everything parsed from one SBGN-ML document."""

    pathway_id: str = ""
    pathway_name: str = ""
    language: str = "process description"
    glyphs: dict[str, SBGNGlyph] = field(default_factory=dict)
    arcs: list[SBGNArc] = field(default_factory=list)
    compartments: dict[str, SBGNGlyph] = field(default_factory=dict)
    ports: dict[str, str] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.glyphs)

    @property
    def resolved_arcs(self) -> list[SBGNArc]:
        """Arcs whose endpoints both resolve to a known glyph."""
        known = set(self.glyphs) | set(self.compartments)
        return [a for a in self.arcs if a.source in known and a.target in known]


# ---------------------------------------------------------------------------
# Namespace-insensitive helpers
# ---------------------------------------------------------------------------

def _localname(tag: object) -> str:
    t = str(tag)
    return t.rsplit("}", 1)[-1] if "}" in t else t


def _children(elem: ET.Element, name: str) -> Iterator[ET.Element]:
    """Direct children with local name *name*."""
    for c in elem:
        if _localname(c.tag) == name:
            yield c


def _descendants(elem: ET.Element, name: str) -> Iterator[ET.Element]:
    """All descendants with local name *name* (excluding *elem* itself)."""
    for c in elem.iter():
        if c is not elem and _localname(c.tag) == name:
            yield c


def _first(elem: ET.Element, name: str) -> ET.Element | None:
    return next(_children(elem, name), None)


def _fnum(v: object, default: float | None = None) -> float | None:
    try:
        return float(str(v))
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# Element parsing
# ---------------------------------------------------------------------------

def _parse_bbox(elem: ET.Element) -> dict:
    """
    Read the glyph's own ``<bbox>``.

    Uses a *direct child* lookup: a descendant search picks up the bbox of a
    nested state-variable glyph and mislocates the parent node.
    """
    bbox = _first(elem, "bbox")
    if bbox is None:
        return {}
    x = _fnum(bbox.get("x"), 0.0) or 0.0
    y = _fnum(bbox.get("y"), 0.0) or 0.0
    w = _fnum(bbox.get("w"), 0.0) or 0.0
    h = _fnum(bbox.get("h"), 0.0) or 0.0
    # SBGN bboxes are top-left anchored; the pipeline works in centres.
    return {"x": x + w / 2.0, "y": y + h / 2.0, "width": w, "height": h}


def _parse_label(elem: ET.Element) -> str:
    lab = _first(elem, "label")
    if lab is None:
        return ""
    text = lab.get("text", "")
    if not text and lab.text:
        text = lab.text
    return (text or "").strip()


def _parse_aux(elem: ET.Element) -> tuple[list[dict], list[dict]]:
    """Collect state variables and units of information from direct children."""
    states, units = [], []
    for child in _children(elem, "glyph"):
        cls = (child.get("class") or "").strip().lower()
        if cls == "state variable":
            sv = _first(child, "state")
            states.append({
                "variable": (sv.get("variable", "") if sv is not None
                             else child.get("variable", "")),
                "value": (sv.get("value", "") if sv is not None
                          else child.get("value", "")),
                "label": _parse_label(child),
            })
        elif cls == "unit of information":
            units.append({"label": _parse_label(child)})
    return states, units


def _parse_identifiers(elem: ET.Element, glyph_id: str) -> list[str]:
    """
    Pull database accessions from RDF annotations and from the glyph id.

    SBGN ids from Reactome/SBGNview are opaque hashes; without cross-references
    an SBGN node cannot be joined to omics data at all.
    """
    found: list[str] = []
    for anno in _descendants(elem, "annotation"):
        for node in anno.iter():
            for val in list(node.attrib.values()) + [node.text or ""]:
                for m in _ID_RX.finditer(str(val)):
                    found.append(f"{m.group('db').lower()}:{m.group('acc')}")
    for m in _ID_RX.finditer(glyph_id):
        found.append(f"{m.group('db').lower()}:{m.group('acc')}")
    return list(dict.fromkeys(found))


def _parse_spline(arc: ET.Element) -> list[tuple[float, float]]:
    """Start, intermediate ``<next>`` and end points of an arc."""
    pts: list[tuple[float, float]] = []
    start = _first(arc, "start")
    if start is not None:
        pts.append((_fnum(start.get("x"), 0.0) or 0.0, _fnum(start.get("y"), 0.0) or 0.0))
    for nxt in _children(arc, "next"):
        pts.append((_fnum(nxt.get("x"), 0.0) or 0.0, _fnum(nxt.get("y"), 0.0) or 0.0))
    end = _first(arc, "end")
    if end is not None:
        pts.append((_fnum(end.get("x"), 0.0) or 0.0, _fnum(end.get("y"), 0.0) or 0.0))
    return pts


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_sbgn(filepath: str | Path) -> SBGNPathway:
    """
    Parse an SBGN-ML document.

    Works with or without XML namespaces, resolves arc endpoints through
    ``<port>`` elements, and keeps auxiliary glyphs attached to their parent
    rather than promoting them to nodes.
    """
    path = Path(filepath)
    if not path.exists():
        raise PathwayNotFoundError(f"SBGN file not found: {path}")

    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        raise ParseError(f"{path.name} is not well-formed XML: {exc}") from exc

    map_elem = next(_descendants(root, "map"), None)
    if map_elem is None:
        map_elem = root if _localname(root.tag) == "map" else root

    pathway = SBGNPathway(
        pathway_id=map_elem.get("id") or path.stem,
        pathway_name=map_elem.get("name") or path.stem,
        language=map_elem.get("language", "process description"),
    )

    # -- Pass 1: ports.  Must precede arcs so endpoints can be resolved. ----
    for glyph_elem in _descendants(map_elem, "glyph"):
        gid = glyph_elem.get("id")
        if not gid:
            continue
        for port in _children(glyph_elem, "port"):
            pid = port.get("id")
            if pid:
                pathway.ports[pid] = gid

    # -- Pass 2: glyphs -----------------------------------------------------
    aux_ids: set[str] = set()
    for glyph_elem in _descendants(map_elem, "glyph"):
        gid = glyph_elem.get("id", "")
        cls = (glyph_elem.get("class") or "").strip().lower()
        if not gid:
            continue
        if cls in AUXILIARY_CLASSES:
            aux_ids.add(gid)
            continue

        bbox = _parse_bbox(glyph_elem)
        states, units = _parse_aux(glyph_elem)
        glyph = SBGNGlyph(
            glyph_id=gid,
            glyph_class=cls,
            label=_parse_label(glyph_elem),
            x=bbox.get("x"), y=bbox.get("y"),
            width=bbox.get("width"), height=bbox.get("height"),
            compartment=glyph_elem.get("compartmentRef"),
            clone_marker=_first(glyph_elem, "clone") is not None,
            state_variables=states,
            units_of_information=units,
            identifiers=_parse_identifiers(glyph_elem, gid),
        )
        if cls == "compartment":
            pathway.compartments[gid] = glyph
        else:
            pathway.glyphs[gid] = glyph

    # Record parentage for nested entities (subunits of a complex).
    for glyph_elem in _descendants(map_elem, "glyph"):
        parent_id = glyph_elem.get("id")
        if not parent_id or parent_id not in pathway.glyphs:
            continue
        for child in _children(glyph_elem, "glyph"):
            cid = child.get("id")
            if cid and cid in pathway.glyphs:
                pathway.glyphs[cid].parent = parent_id

    # -- Pass 3: arcs, resolved through ports -------------------------------
    for arc_elem in _descendants(map_elem, "arc"):
        aid = arc_elem.get("id", "")
        raw_src = arc_elem.get("source", "")
        raw_tgt = arc_elem.get("target", "")
        if not (raw_src and raw_tgt):
            continue
        src = pathway.ports.get(raw_src, raw_src)
        tgt = pathway.ports.get(raw_tgt, raw_tgt)
        if src in aux_ids or tgt in aux_ids:
            continue
        pathway.arcs.append(SBGNArc(
            arc_id=aid or f"arc_{len(pathway.arcs)}",
            arc_class=(arc_elem.get("class") or "").strip().lower(),
            source=src, target=tgt,
            spline_points=_parse_spline(arc_elem),
            raw_source=raw_src, raw_target=raw_tgt,
        ))

    return pathway


def sbgn_to_df(pathway: SBGNPathway, include_processes: bool = True) -> pl.DataFrame:
    """
    Convert an SBGN pathway to the same node schema as ``node_info``.

    This is what lets SBGN maps reuse the whole KEGG colour/render pipeline.
    """
    shape_map = {
        "macromolecule": "roundrectangle",
        "macromolecule multimer": "roundrectangle",
        "simple chemical": "circle",
        "simple chemical multimer": "circle",
        "complex": "octagon",
        "complex multimer": "octagon",
        "nucleic acid feature": "roundrectangle",
        "process": "square", "omitted process": "square",
        "uncertain process": "square",
        "association": "circle", "dissociation": "circle",
        "phenotype": "hexagon",
        "and": "circle", "or": "circle", "not": "circle",
        "source and sink": "circle",
        "perturbing agent": "hexagon",
    }

    records = []
    for gid, g in pathway.glyphs.items():
        if not include_processes and g.is_process:
            continue
        names = g.identifiers or ([g.label] if g.label else [gid])
        records.append({
            "entry_id": gid,
            "name": " ".join(names),
            "kegg_names": [n.split(":")[-1] for n in names],
            "type": g.node_type,
            "x": g.x, "y": g.y,
            "width": g.width, "height": g.height,
            "bgcolor": "#FFFFFF", "fgcolor": "#000000",
            "label": g.label or "",
            "shape": shape_map.get(g.glyph_class, "rectangle"),
            "reaction": "", "component": g.parent or "",
            "size": 1, "link": "",
            "glyph_class": g.glyph_class,
            "clone_marker": g.clone_marker,
            "compartment": g.compartment or "",
        })

    if not records:
        return pl.DataFrame(schema={
            "entry_id": pl.String, "name": pl.String, "kegg_names": pl.List(pl.String),
            "type": pl.String, "x": pl.Float64, "y": pl.Float64,
            "width": pl.Float64, "height": pl.Float64, "bgcolor": pl.String,
            "fgcolor": pl.String, "label": pl.String, "shape": pl.String,
            "reaction": pl.String, "component": pl.String, "size": pl.Int64,
            "link": pl.String, "glyph_class": pl.String,
            "clone_marker": pl.Boolean, "compartment": pl.String,
        })
    return pl.DataFrame(records, schema_overrides={"kegg_names": pl.List(pl.String)})


def sbgn_compartments(pathway: SBGNPathway) -> pl.DataFrame:
    """
    Compartment glyphs as a drawable frame.

    Compartments carry the biology's spatial organisation — a reaction in the
    mitochondrial matrix is a different statement from the same reaction in
    the cytosol — so they are returned as first-class geometry rather than
    being discarded with the rest of the non-entity glyphs.

    Columns: entry_id, label, x, y, width, height, area (largest first, so a
    renderer can paint nested compartments in the right order).
    """
    rows = []
    for gid, g in pathway.compartments.items():
        if g.x is None or g.y is None:
            continue
        w, h = float(g.width or 0), float(g.height or 0)
        rows.append({"entry_id": gid, "label": g.label or "",
                     "x": float(g.x), "y": float(g.y),
                     "width": w, "height": h, "area": w * h})
    if not rows:
        return pl.DataFrame(schema={
            "entry_id": pl.String, "label": pl.String, "x": pl.Float64,
            "y": pl.Float64, "width": pl.Float64, "height": pl.Float64,
            "area": pl.Float64,
        })
    return pl.DataFrame(rows).sort("area", descending=True)


def sbgn_edges(pathway: SBGNPathway) -> pl.DataFrame:
    """Edge table for an SBGN pathway, using resolved (port-aware) endpoints."""
    rows = [{
        "source": a.source, "target": a.target,
        "edge_type": a.arc_class, "subtype": a.arc_class,
        "value": a.arc_id, "source_kind": "arc",
        "reversible": False,
        "n_points": len(a.spline_points),
    } for a in pathway.resolved_arcs]

    if not rows:
        return pl.DataFrame(schema={
            "source": pl.String, "target": pl.String, "edge_type": pl.String,
            "subtype": pl.String, "value": pl.String, "source_kind": pl.String,
            "reversible": pl.Boolean, "n_points": pl.Int64,
        })
    return pl.DataFrame(rows)


def sbgn_canvas(pathway: SBGNPathway, pad: float = 40.0) -> tuple[float, float, float, float]:
    """Return (x0, y0, width, height) covering every glyph and compartment."""
    xs0, ys0, xs1, ys1 = [], [], [], []
    for g in list(pathway.glyphs.values()) + list(pathway.compartments.values()):
        if g.x is None or g.y is None:
            continue
        w, h = (g.width or 0) / 2, (g.height or 0) / 2
        xs0.append(g.x - w)
        xs1.append(g.x + w)
        ys0.append(g.y - h)
        ys1.append(g.y + h)
    if not xs0:
        return (0.0, 0.0, 1200.0, 900.0)
    x0, y0 = min(xs0) - pad, min(ys0) - pad
    return (x0, y0, max(xs1) + pad - x0, max(ys1) + pad - y0)


def arc_resolution_report(pathway: SBGNPathway) -> dict:
    """
    Diagnostic: how many arcs resolved to real glyphs.

    Port resolution is the difference between a handful of edges and a
    complete map, so it is worth being able to assert on.
    """
    total = len(pathway.arcs)
    resolved = len(pathway.resolved_arcs)
    via_port = sum(1 for a in pathway.arcs
                   if a.raw_source in pathway.ports or a.raw_target in pathway.ports)
    return {
        "arcs_total": total,
        "arcs_resolved": resolved,
        "arcs_via_port": via_port,
        "ports_indexed": len(pathway.ports),
        "resolution_rate": (resolved / total) if total else 0.0,
    }
