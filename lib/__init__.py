"""
pathview-plus — KEGG and SBGN pathway visualisation for Python.

A Python reimplementation of R's ``pathview`` and ``SBGNview``, with
independent colour scales for transcript and metabolite data, true vector
output, and full offline operation.

Quick start
-----------
    from pathview import pathview, demo_gene_data, demo_cpd_data

    res = pathview(
        "00020",                       # TCA cycle
        gene_data=demo_gene_data(2),   # RNA-seq log2FC, 2 conditions
        cpd_data=demo_cpd_data(),      # metabolite abundances
        species="human",               # names, codes and taxids all work
        render_mode="vector",          # no KEGG PNG needed
        output_format="pdf",
    )
    print(res.summary())

Offline use
-----------
Species resolution, compound naming, identifier cross-referencing and vector
rendering all work with no network access.  Set ``PATHVIEW_OFFLINE=1`` (or
call ``set_offline()``) to guarantee no call is attempted.
"""

from __future__ import annotations

__version__ = "3.1.0"
__author__ = "Richard Allen White III"
__license__ = "CC-BY-NC-4.0"
__url__ = "https://github.com/raw-lab/pathview-plus"

# -- errors ------------------------------------------------------------------
# -- cache / offline ---------------------------------------------------------
from .bundled import bundled_files, read_bundled_tsv
from .cache import cache_dir, clear_cache, is_offline, set_offline

# -- colour ------------------------------------------------------------------
from .color_mapping import (
                     ColorScale,
                     colorpanel2,
                     compound_scale,
                     draw_color_key,
                     draw_dual_key,
                     gene_scale,
                     list_palettes,
                     make_colormap,
                     node_color,
                     resolve_palette,
                     value_columns,
)

# -- databases ---------------------------------------------------------------
from .databases import (
                     DATABASE_INFO,
                     available_databases,
                     detect_database,
                     download_kegg,
                     download_metacrop,
                     download_metacyc,
                     download_panther,
                     download_pathway,
                     download_reactome,
                     download_smpdb,
                     find_reactome_pathways,
                     list_reactome_pathways,
                     reactome_top_url,
)
from .errors import (
                     MappingError,
                     NetworkError,
                     ParseError,
                     PathviewError,
                     PathwayNotFoundError,
                     RenderError,
                     SpeciesNotFoundError,
)
from .examples import run_offline_examples
from .expansion import ExpansionResult, expand_nodes, expansion_report, split_groups
from .graph_rendering import available_layouts, build_graph, keggview_graph, pathway_metrics

# -- highlighting ------------------------------------------------------------
from .highlighting import (
                     PathwayResult,
                     PathwayResultSet,
                     annotate,
                     change_labels,
                     highlight_edges,
                     highlight_nodes,
                     highlight_path,
)
from .id_mapping import (
                     IdMapResult,
                     compound_xref,
                     cpd_id_map,
                     cpd_name_to_kegg,
                     eg2id,
                     id2eg,
                     supported_cpd_idtypes,
                     supported_gene_idtypes,
)

# -- parsing -----------------------------------------------------------------
from .kgml_parser import (
                     KGMLEdge,
                     KGMLNode,
                     KGMLPathway,
                     KGMLReaction,
                     canvas_size,
                     node_info,
                     parse_kgml,
                     pathway_edges,
)

# -- geometry ----------------------------------------------------------------
from .layout import Extent, NodeBox, RasterFrame, node_boxes, slice_bounds
from .legend import edge_subtypes, kegg_legend, sbgn_legend

# -- data --------------------------------------------------------------------
from .mol_data import (
                     MolSumResult,
                     compound_name,
                     compound_names,
                     demo_cpd_data,
                     demo_gene_data,
                     mol_sum,
                     sim_mol_data,
)
from .node_mapping import NodeMapResult, node_map

# -- organisms ---------------------------------------------------------------
from .organisms import (
                     SpeciesInfo,
                     bioconductor_orgdb_table,
                     default_gene_idtype,
                     get_species_code,
                     kegg_species_code,
                     list_organisms,
                     organism_count,
                     refresh_organism_table,
                     search_organisms,
)

# -- parity ------------------------------------------------------------------
from .parity import FEATURE_MATRIX, feature_table, parity_summary, print_parity

# -- orchestrator (last: it imports from nearly everything above) ------------
from .pathview import pathview  # noqa: E402
from .rdata import rdata_objects, read_rdata

# -- rendering ---------------------------------------------------------------
from .rendering import keggview_native, paint_nodes, render_native_array

# -- SBGN collection ---------------------------------------------------------
from .sbgn_hub import (
                     SBGN_SOURCES,
                     crosswalk_routes,
                     download_sbgn,
                     download_sbgn_batch,
                     find_sbgn_pathway,
                     id_route,
                     list_sbgn_pathways,
                     map_ids_to_sbgn,
                     sbgn_collection_info,
                     sbgn_index,
                     sbgn_url,
                     sbgn_xref,
                     supported_sbgn_idtypes,
)
from .sbgn_parser import (
                     SBGN_ARC_CLASSES,
                     SBGN_GLYPH_CLASSES,
                     SBGNArc,
                     SBGNGlyph,
                     SBGNPathway,
                     arc_resolution_report,
                     parse_sbgn,
                     sbgn_canvas,
                     sbgn_compartments,
                     sbgn_edges,
                     sbgn_to_df,
)
from .sbgnview import sbgn_node_map, sbgnview, sbgnview_batch  # noqa: E402

# -- splines -----------------------------------------------------------------
from .splines import (
                     bezier_to_svg_path,
                     catmull_rom_spline,
                     cubic_bezier,
                     offset_endpoints,
                     points_to_bezier_path,
                     quadratic_bezier,
                     route_edge_spline,
                     smooth_path_svg,
)
from .svg_rendering import keggview_svg, render_edge_svg, render_node_svg

# -- utilities ---------------------------------------------------------------
from .utils import (
                     contrast_text_color,
                     max_abs,
                     random_pick,
                     short_label,
                     strfit,
                     to_hex,
                     to_rgb,
                     wordwrap,
)
from .vector_rendering import (
                     EDGE_STYLE,
                     THEMES,
                     draw_compartments,
                     draw_pathway,
                     keggview_vector,
                     render_vector_array,
)

__all__ = [
    "__version__",
    # orchestrator
    "pathview", "PathwayResult", "PathwayResultSet",
    "sbgnview", "sbgnview_batch", "sbgn_node_map",
    # organisms
    "get_species_code", "kegg_species_code", "SpeciesInfo", "list_organisms",
    "search_organisms", "organism_count", "refresh_organism_table",
    "bioconductor_orgdb_table", "default_gene_idtype",
    # parsing
    "parse_kgml", "node_info", "pathway_edges", "canvas_size",
    "KGMLPathway", "KGMLNode", "KGMLEdge", "KGMLReaction",
    "parse_sbgn", "sbgn_to_df", "sbgn_edges", "sbgn_canvas",
    "arc_resolution_report", "sbgn_compartments",
    "SBGNPathway", "SBGNGlyph", "SBGNArc",
    "SBGN_GLYPH_CLASSES", "SBGN_ARC_CLASSES",
    # data
    "mol_sum", "sim_mol_data", "demo_gene_data", "demo_cpd_data",
    "compound_name", "compound_names", "MolSumResult",
    "node_map", "NodeMapResult",
    "id2eg", "eg2id", "cpd_id_map", "cpd_name_to_kegg", "compound_xref",
    "supported_gene_idtypes", "supported_cpd_idtypes", "IdMapResult",
    # colour
    "ColorScale", "gene_scale", "compound_scale", "node_color", "colorpanel2",
    "make_colormap", "draw_color_key", "draw_dual_key", "list_palettes",
    "resolve_palette", "value_columns",
    # geometry
    "NodeBox", "Extent", "RasterFrame", "node_boxes", "slice_bounds",
    # rendering
    "keggview_native", "keggview_vector", "keggview_graph", "keggview_svg",
    "render_native_array", "render_vector_array",
    "paint_nodes", "draw_pathway", "draw_compartments", "build_graph",
    "pathway_metrics",
    "available_layouts",
    "render_node_svg", "render_edge_svg", "THEMES", "EDGE_STYLE",
    "kegg_legend", "sbgn_legend", "edge_subtypes",
    # splines
    "cubic_bezier", "quadratic_bezier", "catmull_rom_spline",
    "route_edge_spline", "bezier_to_svg_path", "points_to_bezier_path",
    "smooth_path_svg", "offset_endpoints",
    # databases
    "download_kegg", "download_reactome", "download_pathway",
    "download_panther", "download_metacyc", "download_smpdb",
    "download_metacrop",
    # SBGN collection
    "sbgn_index", "list_sbgn_pathways", "find_sbgn_pathway", "download_sbgn",
    "download_sbgn_batch", "sbgn_collection_info", "sbgn_url", "sbgn_xref",
    "crosswalk_routes", "map_ids_to_sbgn", "id_route",
    "supported_sbgn_idtypes", "SBGN_SOURCES",
    # expansion
    "split_groups", "expand_nodes", "expansion_report", "ExpansionResult",
    # bundled data
    "read_bundled_tsv", "bundled_files", "read_rdata", "rdata_objects",
    "run_offline_examples",
    "list_reactome_pathways", "find_reactome_pathways", "reactome_top_url",
    "detect_database",
    "available_databases", "DATABASE_INFO",
    # highlighting
    "highlight_nodes", "highlight_edges", "highlight_path", "change_labels",
    "annotate",
    # utilities
    "wordwrap", "strfit", "short_label", "max_abs", "random_pick",
    "to_rgb", "to_hex", "contrast_text_color",
    # cache
    "set_offline", "is_offline", "cache_dir", "clear_cache",
    # parity
    "FEATURE_MATRIX", "feature_table", "parity_summary", "print_parity",
    # errors
    "PathviewError", "SpeciesNotFoundError", "PathwayNotFoundError",
    "NetworkError", "ParseError", "MappingError", "RenderError",
]
