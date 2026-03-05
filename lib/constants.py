"""
constants.py
Shared type aliases, literals, and package-wide constants.
"""

from typing import Callable, Literal, Optional

import numpy as np

# KEGG REST base URL
KEGG_BASE = "https://rest.kegg.jp"

# Supported aggregation methods for multi-probe summarisation
SumMethod = Literal["sum", "mean", "median", "max", "max_abs", "random"]

# KEGG node types recognised by the renderer
NodeType = Literal["gene", "compound", "map", "ortholog", "group"]

# Non-data columns present on every node DataFrame
NODE_META_COLS = frozenset({
    "entry_id", "name", "type", "x", "y",
    "width", "height", "label", "shape",
    "reaction", "component", "size", "kegg_names",
})

# Valid biological node types (used for filtering)
VALID_NODE_TYPES = {"gene", "enzyme", "compound", "ortholog"}
