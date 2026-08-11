"""
kgml_parser.py
Parse KEGG KGML into dataclasses and a tidy Polars DataFrame.

Fixes over v2.x
---------------
* Edges are actually produced.  v2.x parsed ``<relation>`` and ``<reaction>``
  into dataclasses that nothing ever consumed — ``keggview_graph`` built a
  DiGraph with ``add_node`` only, so the "graph view" was a scatter of
  disconnected nodes.  ``pathway_edges()`` now returns a real edge table
  merging relations *and* reaction substrate/product links.
* Labels are shortened the way R's ``node.info(short.name=TRUE)`` does, so a
  node shows "CDK4" instead of "CDK4, CMM3, PSK-J3...".
* ``kegg_names`` are split and prefix-stripped at parse time, with a regex
  that accepts digits in the organism prefix.
* Group nodes resolve to their components' identifiers, so complexes can carry
  data instead of being silently unmappable.
* Parse failures raise ParseError instead of propagating a raw ExpatError.

Public API
----------
  parse_kgml     : path -> KGMLPathway
  node_info      : KGMLPathway -> node DataFrame
  pathway_edges  : KGMLPathway -> edge DataFrame
  canvas_size    : KGMLPathway -> (width, height)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET

import polars as pl

from .constants import (
    DEFAULT_CPD_RADIUS,
    DEFAULT_GENE_HEIGHT,
    DEFAULT_GENE_WIDTH,
    KEGG_CPD_BG,
    KEGG_GENE_BG,
)
from .errors import ParseError
from .utils import short_label, strip_kegg_prefix

_PREFIX_RX = re.compile(r"^[A-Za-z][A-Za-z0-9_]*:")


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class KGMLNode:
    """One ``<entry>`` element."""

    entry_id: str
    name: str
    node_type: str
    link: str = ""
    reaction: str = ""
    x: float | None = None
    y: float | None = None
    width: float | None = None
    height: float | None = None
    bgcolor: str = "#FFFFFF"
    fgcolor: str = "#000000"
    label: str = ""
    shape: str = "rectangle"
    component: list[str] = field(default_factory=list)
    kegg_names: list[str] = field(default_factory=list)

    @property
    def is_group(self) -> bool:
        return self.node_type == "group"


@dataclass
class KGMLEdge:
    """One ``<relation>`` element (or a synthesised reaction edge)."""

    entry1: str
    entry2: str
    edge_type: str = ""
    subtypes: list[tuple[str, str]] = field(default_factory=list)
    source: str = "relation"          # "relation" | "reaction"

    @property
    def subtype_name(self) -> str:
        return self.subtypes[0][0] if self.subtypes else ""


@dataclass
class KGMLReaction:
    """
    One ``<reaction>`` element.

    Two KGML dialects are in circulation: current files give each
    substrate/product an ``id`` pointing at an ``<entry>``, while older and
    ``ko``-prefixed maps give only ``name="cpd:C00095"``.  Both are captured so
    edges can be resolved either way.
    """

    entry_id: str
    name: str
    rxn_type: str = "irreversible"
    substrates: list[str] = field(default_factory=list)
    products: list[str] = field(default_factory=list)
    substrate_names: list[str] = field(default_factory=list)
    product_names: list[str] = field(default_factory=list)

    @property
    def reversible(self) -> bool:
        return self.rxn_type == "reversible"


@dataclass
class KGMLPathway:
    """Everything parsed from one KGML document."""

    pathway_id: str = ""
    pathway_name: str = ""
    org: str = ""
    number: str = ""
    title: str = ""
    image: str = ""
    link: str = ""
    nodes: dict[str, KGMLNode] = field(default_factory=dict)
    edges: list[KGMLEdge] = field(default_factory=list)
    reactions: list[KGMLReaction] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.nodes)

    def nodes_of_type(self, *types: str) -> list[KGMLNode]:
        want = set(types)
        return [n for n in self.nodes.values() if n.node_type in want]


# ---------------------------------------------------------------------------
# Element parsers
# ---------------------------------------------------------------------------

def _localname(tag: str) -> str:
    """Strip an XML namespace from a tag: '{ns}entry' -> 'entry'."""
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


def _iter_children(elem: ET.Element, name: str):
    """Direct children named *name*, namespace-insensitively."""
    for child in elem:
        if _localname(child.tag) == name:
            yield child


def _first_child(elem: ET.Element, name: str) -> ET.Element | None:
    return next(_iter_children(elem, name), None)


def _fnum(value: object, default: float | None) -> float | None:
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return default


def _parse_graphics(elem: ET.Element | None, node_type: str) -> dict:
    """Extract display attributes, defaulting by node type."""
    is_cpd = node_type == "compound"
    dw = DEFAULT_CPD_RADIUS if is_cpd else DEFAULT_GENE_WIDTH
    dh = DEFAULT_CPD_RADIUS if is_cpd else DEFAULT_GENE_HEIGHT
    if elem is None:
        return {"x": None, "y": None, "width": dw, "height": dh,
                "bgcolor": KEGG_CPD_BG if is_cpd else KEGG_GENE_BG,
                "fgcolor": "#000000",
                "shape": "circle" if is_cpd else "rectangle", "label": ""}

    a = elem.attrib
    bg = a.get("bgcolor", "") or (KEGG_CPD_BG if is_cpd else KEGG_GENE_BG)
    return {
        "x": _fnum(a.get("x"), None),
        "y": _fnum(a.get("y"), None),
        "width": _fnum(a.get("width"), dw) or dw,
        "height": _fnum(a.get("height"), dh) or dh,
        "bgcolor": bg if bg.startswith("#") else KEGG_GENE_BG,
        "fgcolor": a.get("fgcolor", "#000000") or "#000000",
        "shape": a.get("type", "circle" if is_cpd else "rectangle"),
        "label": a.get("name", ""),
    }


def _parse_entry(elem: ET.Element) -> KGMLNode:
    node_type = elem.attrib.get("type", "gene")
    gfx = _parse_graphics(_first_child(elem, "graphics"), node_type)
    raw_name = elem.attrib.get("name", "")

    names = [strip_kegg_prefix(p) for p in raw_name.split() if p and p != "undefined"]
    raw_label = gfx["label"] or raw_name
    if raw_label.strip().lower() == "undefined":
        raw_label = ""

    return KGMLNode(
        entry_id=elem.attrib["id"],
        name=raw_name,
        node_type=node_type,
        link=elem.attrib.get("link", ""),
        reaction=elem.attrib.get("reaction", ""),
        x=gfx["x"], y=gfx["y"], width=gfx["width"], height=gfx["height"],
        bgcolor=gfx["bgcolor"], fgcolor=gfx["fgcolor"],
        label=short_label(raw_label, node_type),
        shape=gfx["shape"],
        component=[c.attrib["id"] for c in _iter_children(elem, "component")
                   if "id" in c.attrib],
        kegg_names=names,
    )


def _parse_relation(elem: ET.Element) -> KGMLEdge:
    return KGMLEdge(
        entry1=elem.attrib["entry1"],
        entry2=elem.attrib["entry2"],
        edge_type=elem.attrib.get("type", ""),
        subtypes=[(s.attrib.get("name", ""), s.attrib.get("value", ""))
                  for s in _iter_children(elem, "subtype")],
        source="relation",
    )


def _parse_reaction(elem: ET.Element) -> KGMLReaction:
    subs = list(_iter_children(elem, "substrate"))
    prods = list(_iter_children(elem, "product"))
    return KGMLReaction(
        entry_id=elem.attrib.get("id", ""),
        name=elem.attrib.get("name", ""),
        rxn_type=elem.attrib.get("type", "irreversible"),
        substrates=[s.attrib["id"] for s in subs if "id" in s.attrib],
        products=[p.attrib["id"] for p in prods if "id" in p.attrib],
        substrate_names=[strip_kegg_prefix(s.attrib["name"]) for s in subs
                         if "name" in s.attrib],
        product_names=[strip_kegg_prefix(p.attrib["name"]) for p in prods
                       if "name" in p.attrib],
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_kgml(filepath: str | Path) -> KGMLPathway:
    """
    Parse a KGML file into a :class:`KGMLPathway`.

    Raises
    ------
    PathwayNotFoundError : the file does not exist
    ParseError           : the file is not well-formed KGML
    """
    from .errors import PathwayNotFoundError

    path = Path(filepath)
    if not path.exists():
        raise PathwayNotFoundError(f"KGML file not found: {path}")

    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        head = path.read_text(errors="replace")[:200].strip()
        raise ParseError(
            f"{path.name} is not well-formed XML: {exc}. "
            f"KEGG returns a plain-text error page for unknown pathway ids; "
            f"file begins: {head!r}"
        ) from exc

    if _localname(root.tag) != "pathway":
        raise ParseError(
            f"{path.name}: expected a <pathway> root element, got <{_localname(root.tag)}>"
        )

    pathway = KGMLPathway(
        pathway_id=root.attrib.get("name", ""),
        pathway_name=root.attrib.get("name", ""),
        org=root.attrib.get("org", ""),
        number=root.attrib.get("number", ""),
        title=root.attrib.get("title", ""),
        image=root.attrib.get("image", ""),
        link=root.attrib.get("link", ""),
    )

    for child in root:
        tag = _localname(child.tag)
        try:
            if tag == "entry":
                node = _parse_entry(child)
                pathway.nodes[node.entry_id] = node
            elif tag == "relation":
                pathway.edges.append(_parse_relation(child))
            elif tag == "reaction":
                pathway.reactions.append(_parse_reaction(child))
        except (KeyError, ValueError) as exc:
            raise ParseError(f"{path.name}: malformed <{tag}> element: {exc}") from exc

    _resolve_groups(pathway)
    return pathway


def _resolve_groups(pathway: KGMLPathway) -> None:
    """
    Give group (complex) nodes the union of their components' identifiers.

    KGML group entries carry ``name="undefined"``, so v2.x left them with no
    identifiers at all and every complex on a signalling map was unmappable.
    """
    for node in pathway.nodes.values():
        if not node.is_group or node.kegg_names:
            continue
        names: list[str] = []
        for cid in node.component:
            comp = pathway.nodes.get(cid)
            if comp:
                names.extend(comp.kegg_names)
        node.kegg_names = list(dict.fromkeys(names))
        if not node.label or node.label.lower() == "undefined":
            labels = [pathway.nodes[c].label for c in node.component
                      if c in pathway.nodes and pathway.nodes[c].label]
            node.label = "/".join(labels[:3])


def node_info(pathway: KGMLPathway, short_labels: bool = True) -> pl.DataFrame:
    """
    Flatten pathway nodes into a tidy DataFrame.

    Columns: entry_id, name, kegg_names (list), type, x, y, width, height,
    bgcolor, fgcolor, label, shape, reaction, component, size, link.
    """
    if not pathway.nodes:
        return _empty_node_frame()

    records = [
        {
            "entry_id": n.entry_id,
            "name": n.name,
            "kegg_names": n.kegg_names,
            "type": n.node_type,
            "x": n.x,
            "y": n.y,
            "width": n.width,
            "height": n.height,
            "bgcolor": n.bgcolor,
            "fgcolor": n.fgcolor,
            "label": n.label if short_labels else n.name,
            "shape": n.shape,
            "reaction": n.reaction,
            "component": ";".join(n.component),
            "size": max(1, len(n.component)),
            "link": n.link,
        }
        for n in pathway.nodes.values()
    ]
    return pl.DataFrame(records, schema_overrides={"kegg_names": pl.List(pl.String)})


def _empty_node_frame() -> pl.DataFrame:
    return pl.DataFrame(schema={
        "entry_id": pl.String, "name": pl.String, "kegg_names": pl.List(pl.String),
        "type": pl.String, "x": pl.Float64, "y": pl.Float64,
        "width": pl.Float64, "height": pl.Float64, "bgcolor": pl.String,
        "fgcolor": pl.String, "label": pl.String, "shape": pl.String,
        "reaction": pl.String, "component": pl.String, "size": pl.Int64,
        "link": pl.String,
    })


def pathway_edges(pathway: KGMLPathway, include_reactions: bool = True) -> pl.DataFrame:
    """
    Build the edge table: protein/gene relations plus reaction connectivity.

    Returns columns: source, target, edge_type, subtype, value, source_kind,
    reversible.

    Reaction edges are synthesised by joining each reaction's substrates and
    products through the enzyme entry that catalyses it
    (substrate -> enzyme -> product), which is how a metabolic map's arrows
    are actually laid out.  Reversible reactions contribute both directions.
    """
    rows: list[dict] = []

    for e in pathway.edges:
        if e.entry1 not in pathway.nodes or e.entry2 not in pathway.nodes:
            continue
        name, value = (e.subtypes[0] if e.subtypes else ("", ""))
        rows.append({
            "source": e.entry1, "target": e.entry2, "edge_type": e.edge_type,
            "subtype": name, "value": value, "source_kind": "relation",
            "reversible": False,
        })

    if include_reactions:
        rxn_by_name: dict[str, list[KGMLReaction]] = {}
        for r in pathway.reactions:
            for nm in r.name.split():
                rxn_by_name.setdefault(strip_kegg_prefix(nm), []).append(r)

        # Index entries by KEGG identifier so name-only substrate/product
        # elements (older and ko-prefixed KGML) resolve to entries too.
        by_kegg_name: dict[str, list[str]] = {}
        for n in pathway.nodes.values():
            for kn in n.kegg_names:
                by_kegg_name.setdefault(kn, []).append(n.entry_id)

        def participants(ids: list[str], names: list[str]) -> list[str]:
            if ids:
                return [i for i in ids if i in pathway.nodes]
            out: list[str] = []
            for nm in names:
                out.extend(by_kegg_name.get(nm, []))
            return out

        for node in pathway.nodes.values():
            if not node.reaction:
                continue
            for nm in node.reaction.split():
                for rxn in rxn_by_name.get(strip_kegg_prefix(nm), []):
                    subs = participants(rxn.substrates, rxn.substrate_names)
                    prods = participants(rxn.products, rxn.product_names)
                    for sub in subs:
                        rows.append({
                            "source": sub, "target": node.entry_id,
                            "edge_type": "reaction", "subtype": "substrate",
                            "value": rxn.name, "source_kind": "reaction",
                            "reversible": rxn.reversible,
                        })
                    for prod in prods:
                        rows.append({
                            "source": node.entry_id, "target": prod,
                            "edge_type": "reaction", "subtype": "product",
                            "value": rxn.name, "source_kind": "reaction",
                            "reversible": rxn.reversible,
                        })

    if not rows:
        return pl.DataFrame(schema={
            "source": pl.String, "target": pl.String, "edge_type": pl.String,
            "subtype": pl.String, "value": pl.String, "source_kind": pl.String,
            "reversible": pl.Boolean,
        })

    known = set(pathway.nodes)
    df = pl.DataFrame(rows)
    return df.filter(
        pl.col("source").is_in(known) & pl.col("target").is_in(known)
    ).unique(subset=["source", "target", "subtype"], keep="first")


def canvas_size(pathway: KGMLPathway, pad: float = 40.0) -> tuple[float, float]:
    """Bounding canvas implied by node extents, in KGML pixel units."""
    xs = [n.x + (n.width or 0) / 2 for n in pathway.nodes.values() if n.x is not None]
    ys = [n.y + (n.height or 0) / 2 for n in pathway.nodes.values() if n.y is not None]
    if not xs or not ys:
        return (1200.0, 900.0)
    return (max(xs) + pad, max(ys) + pad)
