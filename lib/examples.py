"""
examples.py
Comprehensive examples demonstrating pathview.py features including the new
SBGNview-inspired additions: SVG rendering, highlighting, and spline curves.

Run these examples to explore all capabilities.
"""

import polars as pl
from pathview import (
    pathview,
    highlight_nodes,
    highlight_edges,
    highlight_path,
    cubic_bezier,
    kegg_legend,
    sim_mol_data,
)


# ===========================================================================
# Example 1: Basic pathway visualization (PNG output)
# ===========================================================================

def example_basic_png():
    """Most common use case: overlay expression data on KEGG pathway (PNG)."""
    print("\n" + "="*70)
    print("Example 1: Basic pathway visualization (PNG)")
    print("="*70)
    
    # Simulated data
    gene_df = sim_mol_data(mol_type="gene", species="hsa", n_mol=100, n_exp=1)
    
    result = pathview(
        pathway_id="04110",  # Cell cycle
        gene_data=gene_df,
        species="hsa",
        output_format="png",  # Default
        out_suffix="example1_basic",
    )
    print("✓ Generated: hsa04110.example1_basic.png")


# ===========================================================================
# Example 2: SVG output (vector graphics)
# ===========================================================================

def example_svg_output():
    """NEW: Generate scalable SVG instead of pixel-based PNG."""
    print("\n" + "="*70)
    print("Example 2: SVG vector output")
    print("="*70)
    
    gene_df = sim_mol_data(mol_type="gene", species="hsa", n_mol=100, n_exp=1)
    
    result = pathview(
        pathway_id="04110",
        gene_data=gene_df,
        species="hsa",
        output_format="svg",  # NEW: SVG mode
        out_suffix="example2_svg",
    )
    print("✓ Generated: hsa04110.example2_svg.svg")
    print("  → Open in browser or vector graphics editor")
    print("  → Scalable, smaller file size, web-native")


# ===========================================================================
# Example 3: Multi-condition comparison
# ===========================================================================

def example_multi_condition():
    """Visualize multiple experimental conditions side-by-side."""
    print("\n" + "="*70)
    print("Example 3: Multi-condition comparison")
    print("="*70)
    
    # Three experimental conditions
    gene_df = sim_mol_data(mol_type="gene", species="hsa", n_mol=150, n_exp=3)
    gene_df = gene_df.rename({
        "exp1": "Control",
        "exp2": "Treatment_A",
        "exp3": "Treatment_B",
    })
    
    result = pathview(
        pathway_id="04010",  # MAPK signaling
        gene_data=gene_df,
        species="hsa",
        out_suffix="example3_multistate",
        limit={"gene": 2.0, "cpd": 1.0},
    )
    print("✓ Each node is sliced into 3 color bands (one per condition)")


# ===========================================================================
# Example 4: Custom color scales
# ===========================================================================

def example_custom_colors():
    """Use custom color schemes (e.g., ColorBrewer palettes)."""
    print("\n" + "="*70)
    print("Example 4: Custom color scales")
    print("="*70)
    
    gene_df = sim_mol_data(mol_type="gene", species="hsa", n_mol=100, n_exp=1)
    
    result = pathview(
        pathway_id="04151",  # PI3K-Akt signaling
        gene_data=gene_df,
        species="hsa",
        out_suffix="example4_custom_colors",
        low={"gene": "#2166AC", "cpd": "#4575B4"},  # Blue scale
        mid={"gene": "#F7F7F7", "cpd": "#F7F7F7"},  # White
        high={"gene": "#D6604D", "cpd": "#B2182B"}, # Red scale
        limit={"gene": 2.5, "cpd": 1.5},
    )
    print("✓ ColorBrewer RdBu diverging palette applied")


# ===========================================================================
# Example 5: Highlighting (NEW feature)
# ===========================================================================

def example_highlighting():
    """
    NEW: Layer-by-layer modifications after rendering.
    Highlight specific genes, edges, or paths without re-rendering.
    """
    print("\n" + "="*70)
    print("Example 5: Highlighting with composable layers")
    print("="*70)
    
    gene_df = pl.DataFrame({
        "entrez": ["1956", "2099", "5594", "207", "4609"],
        "lfc":    [ 2.3,   -1.1,    1.8,  -0.5,   3.1],
    })
    
    # Base visualization
    result = pathview(
        pathway_id="04010",
        gene_data=gene_df,
        species="hsa",
        out_suffix="example5_base",
    )
    
    # Apply highlighting layers (ggplot2-style)
    print("\nApplying highlighting layers...")
    print("  1. Highlight EGFR and ESR1 in red")
    print("  2. Highlight edge between them in blue")
    print("  3. Highlight entire MAPK cascade path in orange")
    
    # Note: Highlighting requires PathwayResult object which we'll implement fully
    # This is a preview of the API design
    # highlighted = (
    #     result
    #     + highlight_nodes(["1956", "2099"], color="red", width=4)
    #     + highlight_edges([("1956", "2099")], color="blue", width=3)
    #     + highlight_path(["1956", "2099", "5594"], color="orange")
    # )
    # highlighted.save("hsa04010.example5_highlighted.png")
    
    print("✓ Highlighting API demonstrated (full implementation pending)")


# ===========================================================================
# Example 6: Gene symbol IDs
# ===========================================================================

def example_gene_symbols():
    """Use gene symbols instead of Entrez IDs."""
    print("\n" + "="*70)
    print("Example 6: Gene symbol IDs")
    print("="*70)
    
    gene_df = pl.DataFrame({
        "symbol": ["TP53", "EGFR", "KRAS", "PIK3CA", "AKT1"],
        "lfc":    [ -1.8,   2.4,    1.1,    1.5,     0.9],
    })
    
    result = pathview(
        pathway_id="04151",
        gene_data=gene_df,
        species="hsa",
        gene_idtype="SYMBOL",  # Automatic conversion to Entrez
        out_suffix="example6_symbols",
    )
    print("✓ Symbols automatically converted via MyGene.info")


# ===========================================================================
# Example 7: Combined gene + metabolite data
# ===========================================================================

def example_gene_plus_compound():
    """Overlay both gene expression and metabolite abundance."""
    print("\n" + "="*70)
    print("Example 7: Combined gene + metabolite data")
    print("="*70)
    
    gene_df = sim_mol_data(mol_type="gene", species="hsa", n_mol=80, n_exp=1)
    cpd_df = sim_mol_data(mol_type="cpd", n_mol=30, n_exp=1)
    
    result = pathview(
        pathway_id="00010",  # Glycolysis / Gluconeogenesis
        gene_data=gene_df,
        cpd_data=cpd_df,
        species="hsa",
        out_suffix="example7_gene_cpd",
        limit={"gene": 2.0, "cpd": 1.5},
        low={"gene": "green", "cpd": "blue"},
        high={"gene": "red", "cpd": "yellow"},
    )
    print("✓ Gene (rectangles) and compound (ellipses) overlays combined")


# ===========================================================================
# Example 8: Graph layout mode (no PNG background)
# ===========================================================================

def example_graph_layout():
    """Use NetworkX graph layout instead of KEGG PNG background."""
    print("\n" + "="*70)
    print("Example 8: Graph layout mode (PDF output)")
    print("="*70)
    
    gene_df = sim_mol_data(mol_type="gene", species="hsa", n_mol=100, n_exp=1)
    
    result = pathview(
        pathway_id="04010",
        gene_data=gene_df,
        species="hsa",
        kegg_native=False,  # Switch to graph mode
        output_format="pdf",
        out_suffix="example8_graph",
    )
    print("✓ Generated: hsa04010.example8_graph.pdf")
    print("  → NetworkX layout with Seaborn styling")


# ===========================================================================
# Example 9: Spline curves (NEW)
# ===========================================================================

def example_spline_curves():
    """
    NEW: Demonstrate Bezier curve generation.
    (Future: integrate into pathway edge rendering)
    """
    print("\n" + "="*70)
    print("Example 9: Spline curves (Bezier paths)")
    print("="*70)
    
    import matplotlib.pyplot as plt
    
    # Generate cubic Bezier curve
    curve = cubic_bezier(
        p0=(0, 0),      # Start
        p1=(1, 2),      # Control point 1
        p2=(3, 2),      # Control point 2
        p3=(4, 0),      # End
        n_points=100
    )
    
    plt.figure(figsize=(8, 4))
    plt.plot(curve[:, 0], curve[:, 1], 'b-', linewidth=2, label='Cubic Bezier')
    plt.plot([0, 1, 3, 4], [0, 2, 2, 0], 'ro--', alpha=0.5, label='Control points')
    plt.title("Bezier Curve Example")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig("example9_bezier_curve.png", dpi=150)
    plt.close()
    
    print("✓ Generated: example9_bezier_curve.png")
    print("  → Smooth curves for aesthetic edge routing")
    print("  → Future: automatic edge routing with obstacle avoidance")


# ===========================================================================
# Example 10: Display KEGG legend
# ===========================================================================

def example_legend():
    """Show the KEGG pathway element legend."""
    print("\n" + "="*70)
    print("Example 10: KEGG pathway legend")
    print("="*70)
    
    # This displays an interactive legend
    # kegg_legend(legend_type="both")  # Uncomment to show interactive plot
    
    print("✓ Run kegg_legend() to see node/edge reference diagram")


# ===========================================================================
# Run all examples
# ===========================================================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print("PATHVIEW.PY COMPREHENSIVE EXAMPLES")
    print("="*70)
    print("\nThis script demonstrates all features including new SBGNview additions:")
    print("  • SVG vector output")
    print("  • Highlighting layers")
    print("  • Spline curves")
    print("  • Multi-condition visualization")
    print("  • Custom color schemes")
    print("  • Gene + metabolite overlays")
    
    try:
        example_basic_png()
        example_svg_output()
        example_multi_condition()
        example_custom_colors()
        example_highlighting()
        example_gene_symbols()
        example_gene_plus_compound()
        example_graph_layout()
        example_spline_curves()
        example_legend()
        
        print("\n" + "="*70)
        print("✓ All examples completed successfully!")
        print("="*70)
        print("\nCheck the current directory for output files:")
        print("  • hsa04110.example1_basic.png")
        print("  • hsa04110.example2_svg.svg")
        print("  • hsa04010.example3_multistate.png")
        print("  • hsa04151.example4_custom_colors.png")
        print("  • hsa04010.example5_base.png")
        print("  • hsa04151.example6_symbols.png")
        print("  • hsa00010.example7_gene_cpd.png")
        print("  • hsa04010.example8_graph.pdf")
        print("  • example9_bezier_curve.png")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
