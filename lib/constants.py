"""
constants.py
Shared type aliases, literals, and package-wide constants.
"""

from __future__ import annotations

from typing import Literal

# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

KEGG_BASE = "https://rest.kegg.jp"
KEGG_IMAGE_BASE = "https://www.kegg.jp"
MYGENE_URL = "https://mygene.info/v3/querymany"
REACTOME_SBGN = "https://reactome.org/ContentService/exporter/sbgn"
REACTOME_API = "https://reactome.org/ContentService"

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------

SumMethod = Literal["sum", "mean", "median", "max", "min", "max_abs", "random", "first"]
SUM_METHODS: tuple[str, ...] = (
    "sum", "mean", "median", "max", "min", "max_abs", "random", "first",
)

NodeType = Literal["gene", "compound", "map", "ortholog", "group", "enzyme", "reaction"]
OutputFormat = Literal["png", "pdf", "svg", "html"]
RenderMode = Literal["native", "vector", "graph", "svg"]

# ---------------------------------------------------------------------------
# Node schema
#
# BUGFIX (v3): "bgcolor" and "fgcolor" were absent from NODE_META_COLS in
# v2.x, so they were treated as *data* columns.  Every node was therefore
# coloured by int("FFFFFF", 16) = 16777215, clipped to the scale maximum,
# turning every node solid red regardless of the user's data.
# ---------------------------------------------------------------------------

NODE_META_COLS = frozenset({
    "entry_id", "name", "type", "x", "y", "width", "height",
    "bgcolor", "fgcolor", "label", "shape", "reaction", "component",
    "size", "kegg_names", "all_mapped", "link", "compartment",
    "glyph_class", "clone_marker", "mol_val",
})

# Node types that can carry molecular data
VALID_NODE_TYPES = frozenset({"gene", "enzyme", "compound", "ortholog", "group"})
GENE_NODE_TYPES = frozenset({"gene", "enzyme", "ortholog", "group"})
CPD_NODE_TYPES = frozenset({"compound"})

# ---------------------------------------------------------------------------
# Rendering defaults
# ---------------------------------------------------------------------------

DEFAULT_GENE_WIDTH = 46.0
DEFAULT_GENE_HEIGHT = 17.0
DEFAULT_CPD_RADIUS = 8.0

# KEGG's own palette, used when no data maps to a node
KEGG_GENE_BG = "#BFFFBF"
KEGG_CPD_BG = "#FFFFFF"
KEGG_MAP_BG = "#FFFFFF"
KEGG_BORDER = "#000000"

NA_COLOR = "#F2F2F2"

# ---------------------------------------------------------------------------
# Colour presets — named, colour-blind-aware diverging scales
# ---------------------------------------------------------------------------

PALETTES: dict[str, tuple[str, str, str]] = {
    # (low, mid, high)
    # Exact R pathview defaults: R's col2rgb("green")/("gray")/("red").
    # This is the package default so out-of-the-box output is directly
    # comparable to R pathview; the softer scales below are opt-in.
    "pathview":   ("#00FF00", "#BEBEBE", "#FF0000"),
    "pathview_soft": ("#00CC00", "#DCDCDC", "#FF0000"),
    "rdbu":       ("#2166AC", "#F7F7F7", "#B2182B"),   # ColorBrewer RdBu
    "rdylbu":     ("#4575B4", "#FFFFBF", "#D73027"),
    "viridis":    ("#440154", "#21908C", "#FDE725"),
    "cividis":    ("#00204D", "#7C7B78", "#FFEA46"),   # colour-blind safe
    "purpleorange": ("#5E3C99", "#F7F7F7", "#E66101"), # deuteranopia safe
    "tealrose":   ("#009392", "#F1EAC8", "#CF597E"),
    "bluered":    ("#0571B0", "#F7F7F7", "#CA0020"),
    "metabolite": ("#3B82F6", "#F5F5F5", "#F59E0B"),   # blue -> amber
    "rnaseq":     ("#1B7837", "#F7F7F7", "#762A83"),   # green -> purple
}

DEFAULT_GENE_PALETTE = "pathview"
DEFAULT_CPD_PALETTE = "metabolite"

# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

CACHE_TTL_SECONDS = 60 * 60 * 24 * 30      # 30 days
HTTP_TIMEOUT = 30
HTTP_RETRIES = 3
HTTP_BACKOFF = 0.75
