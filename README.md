# Pathview-plus — Complete Pathway Visualization

**Full-featured Python implementation of R pathview + SBGNview with support for KEGG, Reactome, MetaCyc, and more.**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---

## 🎯 Features

### Core Capabilities
- ✅ **KEGG Pathways** — Download and visualize any KEGG pathway
- ✅ **SBGN Pathways** — Support for Reactome, MetaCyc, PANTHER, SMPDB
- ✅ **Multiple Formats** — PNG (native overlay), SVG (vector), PDF (graph layout)
- ✅ **Gene & Metabolite Data** — Overlay expression and abundance data
- ✅ **Multi-Condition** — Visualize multiple experiments side-by-side
- ✅ **ID Conversion** — Automatic mapping: Entrez ↔ Symbol ↔ UniProt ↔ Ensembl
- ✅ **Highlighting** — Post-hoc emphasis of specific nodes/edges/paths
- ✅ **Spline Curves** — Smooth Bezier edge routing
- ✅ **Custom Colors** — Configurable diverging color scales

### New in v2.0
- 🆕 **Full SBGN-ML support** — Parse and render SBGN Process Description files
- 🆕 **Database integration** — Direct download from Reactome, MetaCyc
- 🆕 **SVG vector output** — Scalable graphics for web and publication
- 🆕 **Highlighting system** — ggplot2-style composable modifications
- 🆕 **Spline rendering** — Cubic Bezier and Catmull-Rom curves

---

## 📦 Installation

### Quick install

```bash
pip install pathview-plus
```

### Custom install

```bash
# Clone repository
git clone https://github.com/raw-lab/pathview-plus
cd pathview-plus

# Install dependencies
pip install -r requirements.txt
pip install .

# Or install specific packages
pip install polars numpy matplotlib seaborn Pillow networkx requests
```

**Dependencies:**
- Python ≥ 3.10
- polars ≥ 0.19.0
- matplotlib ≥ 3.7.0
- seaborn ≥ 0.12.0
- numpy ≥ 1.24.0
- Pillow ≥ 10.0.0
- networkx ≥ 3.1
- requests ≥ 2.31.0

---

## 🚀 Quick Start

### 1. Basic KEGG Pathway

```python
import polars as pl
from pathview import pathview

# Load your data
gene_data = pl.read_csv("gene_expr.tsv", separator="\t")

# Visualize on KEGG pathway
result = pathview(
    pathway_id="04110",      # Cell cycle
    gene_data=gene_data,
    species="hsa",
    output_format="png"
)
```

### 2. Reactome SBGN Pathway

```python
from pathview import download_reactome, parse_sbgn, sbgn_to_df, pathview

# Download Reactome pathway
path = download_reactome("R-HSA-109582")  # Hemostasis

# Parse and visualize
pathway = parse_sbgn(path)
node_df = sbgn_to_df(pathway)

# Overlay data
result = pathview(
    pathway_id="R-HSA-109582",
    gene_data=gene_data,
    output_format="svg"  # Vector graphics
)
```

### 3. Multi-Condition Comparison

```python
# Three experimental conditions
gene_data = pl.DataFrame({
    "entrez": ["1956", "2099", "5594", "207"],
    "Control": [0.5, -0.3, 1.2, -0.8],
    "Treatment_A": [2.1, -1.5, 0.4, 1.3],
    "Treatment_B": [1.8, -0.9, 2.3, 0.7],
})

result = pathview(
    pathway_id="04010",  # MAPK signaling
    gene_data=gene_data,
    species="hsa",
    limit={"gene": 2.5, "cpd": 1.5},
)
# Each node shows 3 color bands (one per condition)
```

### 4. Custom Color Schemes

```python
result = pathview(
    pathway_id="04151",
    gene_data=gene_data,
    species="hsa",
    low={"gene": "#2166AC", "cpd": "#4575B4"},   # Blue
    mid={"gene": "#F7F7F7", "cpd": "#F7F7F7"},   # White
    high={"gene": "#D6604D", "cpd": "#B2182B"},  # Red
)
```

---

## 📖 Complete Examples

### Example 1: Gene Symbol IDs

```python
gene_data = pl.DataFrame({
    "symbol": ["TP53", "EGFR", "KRAS", "PIK3CA", "AKT1"],
    "log2fc": [-1.8, 2.4, 1.1, 1.5, 0.9],
})

result = pathview(
    pathway_id="04151",
    gene_data=gene_data,
    species="hsa",
    gene_idtype="SYMBOL",  # Automatic conversion to Entrez
)
```

### Example 2: Combined Gene + Metabolite

```python
from pathview import sim_mol_data

gene_data = sim_mol_data(mol_type="gene", species="hsa", n_mol=80)
cpd_data = sim_mol_data(mol_type="cpd", n_mol=30)

result = pathview(
    pathway_id="00010",  # Glycolysis
    gene_data=gene_data,
    cpd_data=cpd_data,
    species="hsa",
    low={"gene": "green", "cpd": "blue"},
    high={"gene": "red", "cpd": "yellow"},
)
```

### Example 3: SVG Vector Output

```python
result = pathview(
    pathway_id="04110",
    gene_data=gene_data,
    species="hsa",
    output_format="svg",  # Scalable vector graphics
)
# Output: hsa04110.pathview.svg
# - Scalable without quality loss
# - Smaller file size
# - Editable in Inkscape/Illustrator
```

### Example 4: Graph Layout (No PNG Background)

```python
result = pathview(
    pathway_id="04010",
    gene_data=gene_data,
    species="hsa",
    kegg_native=False,     # Use NetworkX layout
    output_format="pdf",
)
# Output: hsa04010.pathview.pdf
```

### Example 5: Highlighting (API Preview)

```python
from pathview import highlight_nodes, highlight_path

result = pathview("04010", gene_data=data)

# Composable modifications (ggplot2-style)
highlighted = (result
               + highlight_nodes(["1956", "2099"], color="red", width=4)
               + highlight_path(["1956", "2099", "5594"], color="orange"))

highlighted.save("highlighted.png")
```

### Example 6: Spline Curves

```python
from pathview import cubic_bezier, catmull_rom_spline
import matplotlib.pyplot as plt

# Smooth Bezier curve
curve = cubic_bezier((0,0), (1,2), (3,2), (4,0), n_points=100)

plt.plot(curve[:, 0], curve[:, 1], linewidth=2)
plt.title("Bezier Curve Edge Routing")
plt.savefig("bezier_example.png")
```

### Example 7: Batch Processing

```python
pathways = ["04110", "04010", "04151", "00010"]

for pw_id in pathways:
    try:
        result = pathview(
            pathway_id=pw_id,
            gene_data=gene_data,
            species="hsa",
            out_suffix=f"batch_{pw_id}",
        )
        print(f"✓ Completed {pw_id}")
    except Exception as e:
        print(f"✗ Failed {pw_id}: {e}")
```

---

## 🖥️ Command Line Interface

```bash
# Basic usage
python pathview_cli.py --pathway-id 04110 --gene-data expr.tsv

# Specify species and ID type
python pathview_cli.py \
    --pathway-id 04110 \
    --species hsa \
    --gene-data expr.tsv \
    --gene-idtype SYMBOL

# Custom colors
python pathview_cli.py \
    --pathway-id 04010 \
    --gene-data expr.tsv \
    --low-gene '#2166AC' \
    --high-gene '#D6604D' \
    --output-format svg

# Simulate data (for testing)
python pathview_cli.py \
    --pathway-id 04110 \
    --simulate \
    --n-sim 200

# Display KEGG legend
python pathview_cli.py --legend
```

**CLI Arguments:**

```
Pathway:
  --pathway-id ID          KEGG pathway number (e.g., '04110')

Input data:
  --gene-data TSV          Gene expression file (TSV)
  --cpd-data TSV           Compound abundance file (TSV)
  --gene-idtype TYPE       Gene ID type: ENTREZ, SYMBOL, UNIPROT, ENSEMBL
  --cpd-idtype TYPE        Compound ID type: KEGG, PUBCHEM, CHEBI

Species & paths:
  --species CODE           KEGG species code (default: hsa)
  --kegg-dir DIR           Directory for files (default: .)
  --out-suffix SUFFIX      Output filename suffix (default: pathview)

Rendering:
  --kegg-native            Use KEGG PNG background (default: True)
  --output-format FORMAT   Output format: png, pdf, svg (default: png)
  --map-symbol             Replace Entrez with symbols (default: True)
  --node-sum METHOD        Aggregation: sum, mean, median, max
  --no-signature           Suppress watermark
  --no-col-key             Suppress color legend

Color scale:
  --limit-gene FLOAT       Color scale limit (default: 1.0)
  --bins-gene INT          Color bins (default: 10)
  --low-gene COLOR         Low-end color (default: green)
  --mid-gene COLOR         Mid-point color (default: gray)
  --high-gene COLOR        High-end color (default: red)
  --low-cpd COLOR          Low compound color (default: blue)
  --high-cpd COLOR         High compound color (default: yellow)

Utilities:
  --legend                 Display KEGG legend and exit
  --simulate               Generate simulated data
  --n-sim INT              Number of simulated molecules (default: 200)
```

---

## 📊 Input File Formats

### Gene Data (TSV)

First column = gene IDs, remaining columns = numeric expression values.

```tsv
entrez	Control	Treatment_A	Treatment_B
1956	2.31	0.45	1.82
2099	-1.14	-0.88	0.33
5594	0.72	1.33	-0.51
207	-0.88	1.21	0.94
```

### Gene Symbols

```tsv
gene_symbol	log2fc	p_value
TP53	-1.8	0.001
EGFR	2.4	0.0001
KRAS	1.1	0.01
```

### Compound Data (TSV)

```tsv
kegg	abundance
C00031	1.45
C00118	-0.83
C00022	2.11
```

---

## 🎨 Color Scale Configuration

### Three-Point Diverging Scale

```python
pathview(
    pathway_id="04110",
    gene_data=data,
    limit={"gene": 2.0, "cpd": 1.5},      # ±2.0 for genes, ±1.5 for compounds
    bins={"gene": 20, "cpd": 10},          # Color resolution
    low={"gene": "blue", "cpd": "green"},
    mid={"gene": "white", "cpd": "gray"},
    high={"gene": "red", "cpd": "yellow"},
)
```

The scale maps:
- `low value` → `low color` (default: green/blue)
- `0` → `mid color` (default: gray)
- `high value` → `high color` (default: red/yellow)

### One-Directional Scale

```python
both_dirs={"gene": False, "cpd": False}
# Maps: 0 (mid) → max (high)
```

---

## 🗂️ Supported ID Types

### Gene IDs

| Type | Value | Example |
|------|-------|---------|
| Entrez | `ENTREZ` | `1956` |
| Symbol | `SYMBOL` | `EGFR` |
| UniProt | `UNIPROT` | `P00533` |
| Ensembl | `ENSEMBL` | `ENSG00000146648` |
| KEGG | `KEGG` | `hsa:1956` |

### Compound IDs

| Type | Value | Example |
|------|-------|---------|
| KEGG | `KEGG` | `C00031` |
| PubChem | `PUBCHEM` | `5793` |
| ChEBI | `CHEBI` | `4167` |

---

## 🧬 Supported Databases

### KEGG
- **Format:** KGML (XML)
- **Species:** 500+ organisms
- **Download:** Automatic via KEGG REST API
- **Example:** `pathway_id="hsa04110"`

### Reactome
- **Format:** SBGN-ML
- **Species:** Human, mouse, rat, and more
- **Download:** `download_reactome("R-HSA-109582")`
- **Example:** Hemostasis, Immune System, Signaling

### MetaCyc
- **Format:** SBGN-ML
- **Coverage:** 2,800+ metabolic pathways
- **Download:** `download_metacyc("PWY-7210")`
- **Example:** Pyrimidine biosynthesis

### PANTHER
- **Format:** SBGN-ML
- **Coverage:** 177 signaling and metabolic pathways
- **Note:** Manual download required

### SMPDB
- **Format:** SBGN-ML
- **Coverage:** Small molecule pathways
- **Note:** Manual download from website

---

## 🏗️ Architecture

```
pathview/
├── __init__.py           # Public API exports
├── constants.py          # Type definitions
├── utils.py              # String/numeric utilities
│
├── id_mapping.py         # Gene/compound ID conversion
├── mol_data.py           # Data aggregation, simulation
│
├── kegg_api.py           # KEGG REST API
├── databases.py          # Reactome, MetaCyc downloaders
│
├── kgml_parser.py        # KEGG KGML (XML) parser
├── sbgn_parser.py        # SBGN-ML (XML) parser
│
├── color_mapping.py      # Colormaps, node coloring
├── node_mapping.py       # Map data onto nodes
│
├── rendering.py          # PNG/PDF renderers
├── svg_rendering.py      # SVG vector renderer
├── highlighting.py       # Post-hoc modifications
├── splines.py            # Bezier curve math
│
└── pathview.py           # Core orchestrator

pathview_cli.py           # Command-line interface
requirements.txt          # Dependencies
README.md                 # This file
```

**Module Statistics:**
- **15 modules** | **3,506 lines of code**
- Functional programming style
- Full type hints
- Comprehensive docstrings

---

## 🔧 API Reference

### Core Function

```python
pathview(
    pathway_id: str,
    gene_data: Optional[pl.DataFrame] = None,
    cpd_data: Optional[pl.DataFrame] = None,
    species: str = "hsa",
    kegg_dir: Path = ".",
    kegg_native: bool = True,
    output_format: str = "png",  # "png", "pdf", "svg"
    gene_idtype: str = "ENTREZ",
    cpd_idtype: str = "KEGG",
    out_suffix: str = "pathview",
    node_sum: str = "sum",
    map_symbol: bool = True,
    map_null: bool = True,
    min_nnodes: int = 3,
    new_signature: bool = True,
    plot_col_key: bool = True,
    # Color scale parameters
    limit: dict = {"gene": 1.0, "cpd": 1.0},
    bins: dict = {"gene": 10, "cpd": 10},
    both_dirs: dict = {"gene": True, "cpd": True},
    low: dict = {"gene": "green", "cpd": "blue"},
    mid: dict = {"gene": "gray", "cpd": "gray"},
    high: dict = {"gene": "red", "cpd": "yellow"},
    na_col: str = "transparent",
) -> dict
```

### Data Functions

```python
sim_mol_data(mol_type="gene", species="hsa", n_mol=100, n_exp=1) → pl.DataFrame
mol_sum(mol_data, id_map, sum_method="sum") → pl.DataFrame
```

### ID Mapping

```python
id2eg(ids, category, org="Hs") → pl.DataFrame
eg2id(eg_ids, category="SYMBOL", org="Hs") → pl.DataFrame
cpd_id_map(in_ids, in_type, out_type="KEGG") → pl.DataFrame
```

### Parsing

```python
# KEGG
parse_kgml(filepath) → KGMLPathway
node_info(pathway) → pl.DataFrame

# SBGN
parse_sbgn(filepath) → SBGNPathway
sbgn_to_df(pathway) → pl.DataFrame
```

### Database Downloads

```python
download_kegg(pathway_id, species="hsa", kegg_dir=".") → dict
download_reactome(pathway_id, output_dir=".") → Path
download_metacyc(pathway_id, output_dir=".") → Path
list_reactome_pathways(species="Homo sapiens") → list[dict]
detect_database(pathway_id) → str
```

### Highlighting

```python
# API design (full implementation in progress)
result = pathview(...)
highlighted = result + highlight_nodes(["1956", "2099"], color="red")
highlighted.save("output.png")
```

### Splines

```python
cubic_bezier(p0, p1, p2, p3, n_points=50) → np.ndarray
quadratic_bezier(p0, p1, p2, n_points=50) → np.ndarray
catmull_rom_spline(points, n_points=50, alpha=0.5) → np.ndarray
route_edge_spline(source, target, obstacles, mode="orthogonal") → np.ndarray
bezier_to_svg_path(curve, close=False) → str
```

---

## 📈 Performance

- **KEGG pathways:** ~2-5 seconds (download + render)
- **SBGN pathways:** ~3-8 seconds (more complex)
- **Multi-condition:** Linear scaling with # conditions
- **Batch processing:** Parallel processing possible

**Optimization tips:**
- Cache downloaded files (automatic)
- Use `output_format="svg"` for faster rendering
- Disable color key for batch jobs: `plot_col_key=False`

---

## 🤝 Contributing

Contributions welcome! Areas for improvement:

1. **SBGN rendering** — Improve glyph shape variety
2. **Edge routing** — Implement A* pathfinding for splines
3. **Database integration** — Add PANTHER, SMPDB auto-download
4. **Highlighting** — Wire up image modification backend
5. **Performance** — Parallel pathway processing

---

## 📄 License

GPL v3.0 — See LICENSE file

**Citations:**
- Original R pathview: Luo & Brouwer (2013). Bioinformatics
- SBGNview: Shashikant et al. (2022). Bioinformatics
- KEGG: Kanehisa et al. (2023). Nucleic Acids Research
- Reactome: Gillespie et al. (2022). Nucleic Acids Research

---

## 📞 Support

- **Issues:** https://github.com/raw-lab/pathview-plus/issues
- **Email:** rwhit101@charlotte.edu 
---

**Made with ❤️ for the pathway visualization community**
