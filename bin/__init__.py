"""
pathview — Python implementation of R pathview + SBGNview features.

Complete pathway visualization system supporting:
- KEGG pathways (KGML format)
- SBGN pathways (Reactome, MetaCyc, PANTHER, SMPDB)
- Multiple rendering modes (PNG overlay, SVG vector, PDF graph)
- Highlighting and post-processing
- Spline curve rendering

Public API
----------
    # Core visualization
    from pathview import pathview
    
    # Data utilities
    from pathview import sim_mol_data, mol_sum, node_color
    
    # ID mapping
    from pathview import id2eg, eg2id, cpd_id_map
    
    # Parsing (KEGG)
    from pathview import parse_kgml, node_info
    
    # Parsing (SBGN)
    from pathview import parse_sbgn, sbgn_to_df
    
    # Database downloaders
    from pathview import (
        download_kegg, download_reactome, download_metacyc,
        list_reactome_pathways, detect_database
    )
    
    # Highlighting & post-processing
    from pathview import (
        PathwayResult, highlight_nodes, highlight_edges,
        highlight_path, change_labels
    )
    
    # Rendering modes
    from pathview import keggview_native, keggview_graph, keggview_svg
    
    # Spline curves
    from pathview import (
        cubic_bezier, quadratic_bezier, catmull_rom_spline,
        route_edge_spline, bezier_to_svg_path
    )
"""

from .color_mapping import draw_color_key, make_colormap, node_color
from .databases import (DATABASE_INFO, detect_database, download_metacyc,
                         download_panther, download_reactome, download_smpdb,
                         list_reactome_pathways)
from .highlighting import (PathwayResult, change_labels, highlight_edges,
                            highlight_nodes, highlight_path)
from .id_mapping import cpd_id_map, eg2id, id2eg
from .kegg_api import SpeciesInfo, download_kegg, kegg_species_code
from .kgml_parser import (KGMLEdge, KGMLNode, KGMLPathway, KGMLReaction,
                           node_info, parse_kgml)
from .mol_data import mol_sum, sim_mol_data
from .node_mapping import node_map
from .rendering import kegg_legend, keggview_graph, keggview_native
from .sbgn_parser import (SBGN_ARC_CLASSES, SBGN_GLYPH_CLASSES, SBGNArc,
                           SBGNGlyph, SBGNPathway, parse_sbgn, sbgn_to_df)
from .splines import (bezier_to_svg_path, catmull_rom_spline, cubic_bezier,
                       quadratic_bezier, route_edge_spline, smooth_path_svg)
from .svg_rendering import keggview_svg, render_edge_svg, render_node_svg
from .utils import max_abs, random_pick, wordwrap

__all__ = [
    # Core pipeline
    "pathview",
    
    # Data simulation & aggregation
    "sim_mol_data", "mol_sum",
    
    # ID mapping
    "id2eg", "eg2id", "cpd_id_map",
    
    # KEGG API
    "kegg_species_code", "download_kegg", "SpeciesInfo",
    
    # Database downloads (SBGN)
    "download_reactome", "download_metacyc", "download_panther", "download_smpdb",
    "list_reactome_pathways", "detect_database", "DATABASE_INFO",
    
    # Parsing (KGML)
    "parse_kgml", "node_info",
    "KGMLPathway", "KGMLNode", "KGMLEdge", "KGMLReaction",
    
    # Parsing (SBGN)
    "parse_sbgn", "sbgn_to_df",
    "SBGNPathway", "SBGNGlyph", "SBGNArc",
    "SBGN_GLYPH_CLASSES", "SBGN_ARC_CLASSES",
    
    # Node mapping
    "node_map",
    
    # Colors
    "node_color", "make_colormap", "draw_color_key",
    
    # Rendering (PNG, PDF, SVG)
    "keggview_native", "keggview_graph", "keggview_svg", "kegg_legend",
    "render_node_svg", "render_edge_svg",
    
    # Highlighting & post-processing
    "PathwayResult", "highlight_nodes", "highlight_edges",
    "highlight_path", "change_labels",
    
    # Spline curves
    "cubic_bezier", "quadratic_bezier", "catmull_rom_spline",
    "route_edge_spline", "bezier_to_svg_path", "smooth_path_svg",
    
    # Utilities
    "wordwrap", "max_abs", "random_pick",
]

# Import pathview last to avoid circular imports
from .pathview import pathview  # noqa: E402

__version__ = "2.0.0"
__author__ = "pathview.py contributors"
__description__ = "KEGG + SBGN pathway visualization with Python"
