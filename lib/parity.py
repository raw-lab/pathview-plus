"""
parity.py
Feature parity matrix against R ``pathview`` and R ``SBGNview``.

This is a first-class module rather than a paragraph in the README so the
claims are testable: ``tests/test_parity.py`` asserts that every feature
marked ``full`` in :data:`FEATURE_MATRIX` is actually importable and callable,
and the docs build renders the table from this same source.  A feature cannot
be marked supported here without a corresponding test.

Status vocabulary
-----------------
``full``     implemented and covered by a test
``partial``  implemented with a documented limitation
``none``     not implemented
``n/a``      not applicable to that package

The R columns record what the R packages do, taken from their published
sources (pathview 1.46, SBGNview 1.20): ``pathview::pathview`` and
``SBGNview::SBGNview`` argument lists and exported functions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

Status = Literal["full", "partial", "none", "n/a"]


@dataclass(frozen=True)
class Feature:
    """One comparable capability."""

    category: str
    name: str
    pathview_plus: Status
    pathview_r: Status
    sbgnview_r: Status
    api: str = ""
    note: str = ""


# ---------------------------------------------------------------------------
# The matrix
# ---------------------------------------------------------------------------

FEATURE_MATRIX: tuple[Feature, ...] = (
    # -- species / identifiers ---------------------------------------------
    Feature("Species", "KEGG organism code lookup", "full", "full", "full",
            "get_species_code", "10,718 organisms bundled"),
    Feature("Species", "Lookup works offline", "full", "partial", "none",
            "get_species_code",
            "R pathview bundles korg but falls back to a network fetch for misses"),
    Feature("Species", "Common/scientific name and taxid input", "full", "full", "partial",
            "get_species_code"),
    Feature("Species", "Fuzzy match with suggestions", "full", "none", "none",
            "search_organisms"),
    Feature("Species", "Refresh organism table from KEGG", "full", "full", "n/a",
            "refresh_organism_table"),

    # -- input formats -------------------------------------------------------
    Feature("Input", "KEGG KGML parsing", "full", "full", "partial",
            "parse_kgml"),
    Feature("Input", "SBGN-ML parsing", "full", "none", "full", "parse_sbgn"),
    Feature("Input", "SBGN namespaced documents", "full", "n/a", "full",
            "parse_sbgn", "Reactome exports declare a default namespace"),
    Feature("Input", "SBGN arc port resolution", "full", "n/a", "full",
            "arc_resolution_report"),
    Feature("Input", "SBGN state variables / clone markers", "full", "n/a", "full",
            "SBGNGlyph"),
    Feature("Input", "KGML reaction -> edge synthesis", "full", "full", "n/a",
            "pathway_edges"),
    Feature("Input", "Group/complex node splitting", "full", "full", "n/a",
            "split_groups", "R's split.group; expansion is recorded so a "
            "figure's provenance can be audited"),
    Feature("Input", "Multi-identifier node expansion", "full", "full", "n/a",
            "expand_nodes", "R's expand.node; sub-nodes tile the original box "
            "exactly and edges are remapped onto them"),

    # -- identifier mapping --------------------------------------------------
    Feature("ID mapping", "Entrez / KEGG gene ids", "full", "full", "full", "node_map"),
    Feature("ID mapping", "Symbol, Ensembl, UniProt, RefSeq", "full", "full", "full",
            "id2eg", "via MyGene.info instead of Bioconductor OrgDb"),
    Feature("ID mapping", "Symbol/Ensembl/UniProt conversion without R",
            "full", "n/a", "n/a", "id2eg",
            "MyGene.info covers what OrgDb provides, with no R dependency"),
    Feature("ID mapping", "Bioconductor OrgDb packages themselves", "n/a", "full", "full",
            "", "An R library cannot be imported from Python; the "
                "conversions it provides are covered by id2eg/eg2id and the "
                "bundled crosswalks"),
    Feature("ID mapping", "SBGN glyph-id crosswalks", "full", "n/a", "full",
            "map_ids_to_sbgn",
            "770k pairs bundled: ChEBI, KEGG, Entrez, KO, symbol, compound "
            "name, Pathway Commons"),
    Feature("ID mapping", "Multi-hop identifier routing", "full", "none", "none",
            "id_route", "Breadth-first search over the crosswalk graph"),
    Feature("ID mapping", "Compound cross-references (CAS, ChEBI, HMDB, ...)",
            "full", "full", "full", "cpd_id_map", "14k pairs bundled, offline"),
    Feature("ID mapping", "Compound name -> KEGG accession", "full", "full", "full",
            "cpd_name_to_kegg", "includes conjugate-base derivation"),
    Feature("ID mapping", "Batched / cached lookups", "full", "partial", "partial",
            "cpd_id_map"),
    Feature("ID mapping", "Multi-probe aggregation", "full", "full", "full",
            "mol_sum", "8 methods incl. max_abs and seeded random"),

    # -- colour --------------------------------------------------------------
    Feature("Colour", "Diverging scale with midpoint", "full", "full", "full",
            "ColorScale"),
    Feature("Colour", "R colorpanel2-identical binning", "full", "full", "n/a",
            "colorpanel2"),
    Feature("Colour", "Discrete levels", "full", "full", "partial", "ColorScale"),
    Feature("Colour", "Independent gene and compound scales", "full", "full", "full",
            "gene_scale"),
    Feature("Colour", "Two colour keys on one figure", "full", "partial", "partial",
            "draw_dual_key",
            "R pathview draws a single key; SBGNview draws separate legends per file"),
    Feature("Colour", "Named colour-blind-safe palettes", "full", "none", "none",
            "list_palettes"),
    Feature("Colour", "Value transform before binning", "full", "full", "full",
            "ColorScale"),
    Feature("Colour", "Multi-condition node slicing", "full", "full", "full",
            "node_color", "one column per condition"),

    # -- rendering -----------------------------------------------------------
    Feature("Render", "Overlay on KEGG PNG", "full", "full", "n/a",
            "keggview_native"),
    Feature("Render", "Vector map drawn from coordinates", "full", "partial", "full",
            "keggview_vector",
            "R pathview's graph mode uses Rgraphviz; this draws KEGG's own layout"),
    Feature("Render", "Standalone SVG output", "full", "none", "full",
            "keggview_svg"),
    Feature("Render", "PDF output", "full", "full", "full", "pathview"),
    Feature("Render", "PNG output", "full", "full", "full", "pathview"),
    Feature("Render", "Graph/network view with edges", "full", "full", "n/a",
            "keggview_graph"),
    Feature("Render", "KEGG edge subtype styling", "full", "full", "n/a",
            "EDGE_STYLE", "17 subtypes from KEGG's own table"),
    Feature("Render", "Spline / Bezier edge routing", "full", "none", "full",
            "route_edge_spline"),
    Feature("Render", "Compartment shading", "full", "none", "full",
            "draw_compartments",
            "Nested compartments shaded largest-first with decreasing opacity "
            "and labelled"),
    Feature("Render", "Renders with no network access", "full", "none", "none",
            "keggview_vector"),
    Feature("Render", "Themes", "full", "none", "partial", "THEMES"),
    Feature("Render", "Element legend", "full", "full", "partial", "kegg_legend"),

    # -- post-processing -----------------------------------------------------
    Feature("Post", "Composable modification (result + layer)", "full", "none", "full",
            "PathwayResult"),
    Feature("Post", "Highlight nodes", "full", "none", "full", "highlight_nodes"),
    Feature("Post", "Highlight edges/arcs", "full", "none", "full", "highlight_edges"),
    Feature("Post", "Highlight a path", "full", "none", "full", "highlight_path"),
    Feature("Post", "Change node labels", "full", "none", "full", "change_labels"),
    Feature("Post", "Free-text annotation", "full", "none", "partial", "annotate"),

    # -- data sources --------------------------------------------------------
    Feature("Sources", "KEGG download", "full", "full", "full", "download_kegg"),
    Feature("Sources", "Reactome SBGN download", "full", "none", "full",
            "download_reactome"),
    Feature("Sources", "Pre-generated SBGN collection", "full", "n/a", "full",
            "list_sbgn_pathways",
            "5,206 pathways indexed in the wheel; files fetched on demand and "
            "cached, so nothing unused is downloaded"),
    Feature("Sources", "Collection browsable offline", "full", "n/a", "partial",
            "sbgn_collection_info",
            "The index ships in the wheel; SBGNview.data must be installed in "
            "full (69 MB) to list its contents"),
    Feature("Sources", "PANTHER download", "full", "n/a", "full",
            "download_panther", "152 pathways via the pre-generated collection"),
    Feature("Sources", "MetaCyc download", "full", "n/a", "full",
            "download_metacyc", "2,518 pathways; BioCyc's own API needs a "
            "subscription"),
    Feature("Sources", "SMPDB download", "full", "n/a", "full",
            "download_smpdb", "725 pathways; SMPDB publishes only bulk archives"),
    Feature("Sources", "MetaCrop download", "full", "n/a", "full",
            "download_metacrop", "62 crop-plant metabolic pathways"),
    Feature("Sources", "Local SBGN files from any source", "full", "n/a", "full",
            "parse_sbgn", "A hand-exported file parses identically to a "
            "downloaded one; there is no second code path"),
    Feature("Sources", "Pathway search by name", "full", "none", "full",
            "find_reactome_pathways"),

    # -- SBGN rendering ------------------------------------------------------
    Feature("SBGN", "Top-level SBGN render entry point", "full", "n/a", "full",
            "sbgnview", "sbgnview() is to SBGN what pathview() is to KEGG"),
    Feature("SBGN", "Omics overlay on SBGN glyphs", "full", "none", "full",
            "sbgn_node_map",
            "Tries every identifier system the glyphs might use and keeps "
            "whichever lands on the map"),
    Feature("SBGN", "SBGN batch rendering", "full", "n/a", "full",
            "sbgnview_batch"),
    Feature("SBGN", "Arc-class styling", "full", "n/a", "full", "EDGE_STYLE"),
    Feature("SBGN", "Process / operator glyph rendering", "full", "n/a", "full",
            "sbgn_to_df"),
    Feature("SBGN", "Compartment-aware layout extent", "full", "n/a", "full",
            "sbgn_compartments",
            "The canvas is widened so shading is never clipped"),

    # -- engineering ---------------------------------------------------------
    Feature("Engineering", "Typed errors", "full", "none", "none", "PathviewError"),
    Feature("Engineering", "Mapping diagnostics returned", "full", "partial", "partial",
            "PathwayResult"),
    Feature("Engineering", "Disk cache with TTL", "full", "partial", "partial",
            "cache_dir", "R caches KGML files only"),
    Feature("Engineering", "Retry with backoff", "full", "none", "none", "clear_cache"),
    Feature("Engineering", "Enforced offline mode", "full", "none", "none", "set_offline"),
    Feature("Engineering", "Command-line interface", "full", "none", "none",
            "cli"),
    Feature("Engineering", "Graph metrics", "full", "none", "none", "pathway_metrics"),
    Feature("Engineering", "Multi-pathway batch in one call", "full", "full", "full",
            "PathwayResultSet",
            "A failed pathway is recorded rather than aborting the batch, and "
            "modifiers broadcast across the set"),
    Feature("Engineering", "Reads R .RData without R", "full", "n/a", "n/a",
            "read_rdata", "XDR reader covering the vectors and lists "
                          "reference data is published as; factors, closures "
                          "and S4 raise rather than guess"),
)


# ---------------------------------------------------------------------------
# Views
# ---------------------------------------------------------------------------

def feature_table(category: str | None = None, status: Status | None = None):
    """Return the matrix as a Polars DataFrame, optionally filtered."""
    import polars as pl

    rows = [
        {"category": f.category, "feature": f.name,
         "pathview_plus": f.pathview_plus, "pathview_R": f.pathview_r,
         "SBGNview_R": f.sbgnview_r, "api": f.api, "note": f.note}
        for f in FEATURE_MATRIX
        if (category is None or f.category == category)
        and (status is None or f.pathview_plus == status)
    ]
    return pl.DataFrame(rows)


def parity_summary() -> dict:
    """Counts per status, plus the two headline parity percentages."""
    n = len(FEATURE_MATRIX)
    counts = {s: sum(1 for f in FEATURE_MATRIX if f.pathview_plus == s)
              for s in ("full", "partial", "none", "n/a")}

    def covered(attr: str) -> tuple[int, int]:
        rel = [f for f in FEATURE_MATRIX if getattr(f, attr) in ("full", "partial")]
        have = [f for f in rel if f.pathview_plus in ("full", "partial")]
        return len(have), len(rel)

    pv_have, pv_total = covered("pathview_r")
    sv_have, sv_total = covered("sbgnview_r")

    return {
        "total_features": n,
        **counts,
        "vs_pathview_R": f"{pv_have}/{pv_total}",
        "vs_pathview_R_pct": round(100 * pv_have / pv_total, 1) if pv_total else 0.0,
        "vs_SBGNview_R": f"{sv_have}/{sv_total}",
        "vs_SBGNview_R_pct": round(100 * sv_have / sv_total, 1) if sv_total else 0.0,
        "beyond_both": sum(1 for f in FEATURE_MATRIX
                           if f.pathview_plus == "full"
                           and f.pathview_r in ("none", "n/a")
                           and f.sbgnview_r in ("none", "n/a")),
        "gaps": [f.name for f in FEATURE_MATRIX if f.pathview_plus == "none"],
    }


_MARK = {"full": "yes", "partial": "partial", "none": "no", "n/a": "n/a"}


def print_parity(category: str | None = None, markdown: bool = False) -> str:
    """Render the matrix as text (or GitHub-flavoured Markdown)."""
    rows = [f for f in FEATURE_MATRIX if category is None or f.category == category]
    out: list[str] = []

    if markdown:
        out.append("| Category | Feature | pathview-plus | pathview (R) | SBGNview (R) | Notes |")
        out.append("|---|---|:--:|:--:|:--:|---|")
        for f in rows:
            out.append(
                f"| {f.category} | {f.name} | {_MARK[f.pathview_plus]} | "
                f"{_MARK[f.pathview_r]} | {_MARK[f.sbgnview_r]} | {f.note} |"
            )
    else:
        w = max(len(f.name) for f in rows) + 2
        header = f"{'Feature':<{w}} {'plus':<9}{'pathview':<10}{'SBGNview':<10}"
        out.append(header)
        out.append("-" * len(header))
        current = None
        for f in rows:
            if f.category != current:
                current = f.category
                out.append(f"\n[{current}]")
            out.append(f"{f.name:<{w}} {_MARK[f.pathview_plus]:<9}"
                       f"{_MARK[f.pathview_r]:<10}{_MARK[f.sbgnview_r]:<10}")

    text = "\n".join(out)
    return text
