# Pathview-plus — Complete Pathway Visualization

**Full-featured Python implementation of R pathview + SBGNview, with support for KEGG, Reactome, MetaCyc, PANTHER, SMPDB and MetaCrop.**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/pathview-plus?period=total&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=downloads)](https://pepy.tech/projects/pathview-plus)
[![CI](https://github.com/raw-lab/pathview-plus/actions/workflows/ci.yml/badge.svg)](https://github.com/raw-lab/pathview-plus/actions/workflows/ci.yml)
[![Docs](https://readthedocs.org/projects/pathview-plus/badge/?version=latest)](https://pathview-plus.readthedocs.io)
[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

---

## 🎯 Features

### Core capabilities

- ✅ **KEGG pathways** — download and visualise any KEGG pathway, 10,718 organisms
- ✅ **SBGN pathways** — Reactome, MetaCyc, PANTHER, SMPDB, MetaCrop
- ✅ **Five render modes** — native overlay, vector, SVG, graph, auto
- ✅ **Gene *and* metabolite data** — with **independent colour scales and two colour keys**
- ✅ **Multi-condition** — each node splits into one band per experiment
- ✅ **ID conversion** — Entrez ↔ Symbol ↔ UniProt ↔ Ensembl ↔ KEGG ↔ ChEBI ↔ Pathway Commons
- ✅ **Highlighting** — composable, ggplot2-style post-hoc emphasis
- ✅ **Batch rendering** — many pathways in one call, partial failure preserved
- ✅ **Compartment shading** — nested, labelled, opacity-ordered
- ✅ **Complex splitting** — `split_group` and `expand_node`, area-preserving
- ✅ **Works offline** — species, compound names, crosswalks and vector rendering need no network

### New in v3.1

- 🆕 **Pre-generated SBGN collection** — 5,206 pathways indexed *inside the wheel*
- 🆕 **PANTHER / MetaCyc / SMPDB / MetaCrop auto-download** — no manual step
- 🆕 **`sbgnview()`** — the SBGN counterpart to `pathview()`
- 🆕 **770k identifier crosswalk pairs** bundled, with breadth-first routing
- 🆕 **Multi-pathway batch** returning a `PathwayResultSet`
- 🆕 **`read_rdata()`** — read R `.RData` files without R

### New in v3.0

- 🆕 **21 confirmed bugs fixed**, each with a regression test — see [BUG_CHECKLIST.md](BUG_CHECKLIST.md)
- 🆕 **Vector renderer** — draws the map from KGML coordinates, no background image needed
- 🆕 **Typed errors and mapping diagnostics** instead of silent empty results
- 🆕 **Full CLI** and a **316-test offline suite**

---

## 📦 Installation

### Quick install

```bash
pip install pathview-plus
```

### From source

```bash
git clone https://github.com/raw-lab/pathview-plus
cd pathview-plus
pip install .
```

### With optional extras

```bash
pip install "pathview-plus[layouts]"   # SciPy, for the kamada_kawai graph layout
pip install "pathview-plus[dev]"       # pytest, ruff, build
pip install "pathview-plus[docs]"      # Sphinx toolchain
```

**Dependencies:** Python ≥ 3.10 · polars ≥ 1.0 · numpy ≥ 1.24 · matplotlib ≥ 3.7 · networkx ≥ 3.0 · Pillow ≥ 10.0 · requests ≥ 2.28

Verify the install and see what shipped with it:

```bash
pathview-plus info
pathview-plus sbgn --list
python -m pathview.test_all_features    # offline smoke test, ~2 seconds
```

---

![workflow](https://raw.githubusercontent.com/raw-lab/pathview-plus/main/pathview_plus_workflow.jpg)

---

## 🚀 Quick start

### 1. Basic KEGG pathway

```python
import polars as pl
from pathview import pathview

gene_data = pl.read_csv("gene_expr.tsv", separator="\t")

result = pathview(
    "04110",                 # Cell cycle
    gene_data=gene_data,
    species="hsa",           # or "human", "Homo sapiens", "9606"
    output_format="png",
)
print(result.summary())
```

### 2. Genes *and* metabolites, on separate scales

This is the case the package is built around. A transcript at +2 and a
metabolite at +2 are not the same statement, so they get **independent scales
and independent keys**.

```python
result = pathview(
    "00020",                             # TCA cycle
    gene_data=rna_df,
    cpd_data=metabolite_df,
    species="human",
    limit={"gene": 2.0, "cpd": 1.5},     # different scales, different limits
    gene_color="rnaseq",                 # green → grey → red
    cpd_color="metabolite",              # blue → white → amber
    render_mode="vector",
    output_format="pdf",
)
```

Metabolomics platforms export names, not accessions. `cpd_idtype="NAME"`
resolves them offline, including conjugate-base forms — *Pyruvate* finds
*Pyruvic acid*, *Citrate* finds *Citric acid*:

```python
cpds = pl.DataFrame({
    "name":   ["Pyruvate", "Citrate", "2-Oxoglutarate", "Succinate"],
    "log2fc": [ 1.4,       -0.8,       0.3,              -1.6],
})
result = pathview("00020", cpd_data=cpds, cpd_idtype="NAME", species="hsa")
```

### 3. SBGN pathways

`sbgnview()` is to SBGN what `pathview()` is to KEGG.

```python
from pathview import sbgnview, list_sbgn_pathways

list_sbgn_pathways("panther", limit=5)      # browse offline

result = sbgnview(
    "SMP00001",                # or P00001, GLYCOLYSIS, R-HSA-109688
    gene_data=gene_data,
    gene_idtype="SYMBOL",      # densest offline route — see the note below
    show_compartments=True,
)
```

A file exported by hand from any SBGN source parses identically — there is no
second code path:

```python
result = sbgnview("downloads/my_export.sbgn", gene_data=gene_data)
```

### 4. Multi-condition comparison

```python
gene_data = pl.DataFrame({
    "entrez":      ["1956", "2099", "5594", "207"],
    "Control":     [ 0.5,   -0.3,    1.2,   -0.8],
    "Treatment_A": [ 2.1,   -1.5,    0.4,    1.3],
    "Treatment_B": [ 1.8,   -0.9,    2.3,    0.7],
})

result = pathview("04010", gene_data=gene_data, species="hsa", limit=2.5)
# Each node shows three colour bands, one per condition
```

### 5. Batch rendering

Pass a sequence and get a `PathwayResultSet`. One unavailable pathway does not
discard the rest.

```python
rs = pathview(["04110", "04010", "04151", "00010"],
              gene_data=gene_data, species="hsa")

print(rs.summary())
rs.to_frame()                  # one row per pathway with status
rs["04110"].output_path        # index by id, by position, or iterate
```

### 6. Highlighting

```python
from pathview import highlight_nodes, highlight_path, change_labels

result = pathview("04010", gene_data=gene_data, species="hsa")

annotated = (result
             + highlight_nodes(["1956", "2099"], color="red", width=4)
             + highlight_path(["1956", "2099", "5594"], color="orange")
             + change_labels({"1956": "EGFR *"}))

annotated.save("highlighted.png")
```

Each `+` returns a new result; the original is untouched. Highlights are drawn
through an explicit KGML↔raster transform, so they land on the nodes they mark
in every render mode.

### 7. Complexes and paralogue families

A KEGG group entry stacks several subunits in one box; a single entry may carry
several gene ids (`CDK4, CDK6`). One colour across either hides disagreement
between members.

```python
result = pathview("04110", gene_data=gene_data, species="hsa",
                  split_group=True,     # complexes → their subunits
                  expand_node=True)     # multi-gene nodes → one node per gene

print(result.diagnostics["expansion"])
# 115 -> 96 nodes; 19 complexes split; 96 -> 170 nodes; 35 multi-gene nodes expanded
```

Sub-nodes tile the original box exactly, so the layout is unchanged, and edges
are remapped onto them.

### 8. Custom colour scale

```python
from pathview import ColorScale

scale = ColorScale(limit=(-3, 3), bins=12,
                   low="#2166AC", mid="#F7F7F7", high="#B2182B",
                   label="log2 FC (tumour / normal)")

result = pathview("04151", gene_data=gene_data, species="hsa",
                  gene_color=scale)
```

### 9. Working offline

```python
from pathview import pathview, set_offline

set_offline(True)                  # or export PATHVIEW_OFFLINE=1

result = pathview("00020", gene_data=gene_data, species="hsa",
                  kegg_dir="/data/kgml",     # contains hsa00020.xml
                  render_mode="vector")      # draws the map itself
```

Species resolution, compound naming, identifier crosswalks and vector
rendering all keep working — those tables ship inside the package.

---

## 🖥️ Command line

```bash
# KEGG
pathview-plus render 00020 --species human \
    --gene-data rnaseq.csv --cpd-data metabolites.csv \
    --render-mode vector --output-format pdf --limit gene=2,cpd=1

# Several pathways at once
pathview-plus render 04110 04010 00010 --gene-data rna.csv -o figures

# Complexes and paralogues
pathview-plus render 04110 --gene-data rna.csv --split-group --expand-node

# SBGN
pathview-plus sbgn --list                          # what's in the collection
pathview-plus sbgn --list --source panther -n 20
pathview-plus sbgn SMP00001 --gene-data rna.csv --gene-idtype SYMBOL

# Utilities
pathview-plus species 'Mus musculus'
pathview-plus search coli --limit 5
pathview-plus download hsa04110 -o kgml
pathview-plus legend --out legend.png
pathview-plus parity --markdown > PARITY.md
pathview-plus info
```

**Key `render` options**

| Option | Meaning |
|---|---|
| `--species`, `-s` | code, name or taxid (`hsa`, `human`, `9606`) |
| `--gene-data` / `--cpd-data` | CSV/TSV: ids in column 1, values after |
| `--gene-idtype` | `ENTREZ`, `SYMBOL`, `ENSEMBL`, `UNIPROT`, `REFSEQ`, `KEGG` |
| `--cpd-idtype` | `KEGG`, `NAME`, `CAS`, `CHEBI`, `HMDB`, `PUBCHEM` |
| `--render-mode` | `auto`, `native`, `vector`, `graph`, `svg` |
| `--output-format`, `-f` | `png`, `pdf`, `svg` |
| `--limit` | `1.5`, or `gene=2,cpd=1` for separate scales |
| `--node-sum` | `sum`, `mean`, `median`, `max`, `min`, `max_abs`, `random`, `first` |
| `--split-group` / `--expand-node` | split complexes / expand paralogue families |
| `--theme` | `publication`, `slate`, `dark` |
| `--offline` | never attempt a network request |

Either data file may be given alone — a metabolomics-only run is fully
supported.

---

## 📊 Input file formats

First column = identifiers, remaining numeric columns = one per condition.

**Gene data**

```tsv
entrez	Control	Treatment_A	Treatment_B
1956	2.31	0.45	1.82
2099	-1.14	-0.88	0.33
5594	0.72	1.33	-0.51
```

**Gene symbols** — pass `gene_idtype="SYMBOL"`

```tsv
symbol	log2fc
TP53	-1.8
EGFR	2.4
KRAS	1.1
```

**Compound data** — accessions or names

```tsv
kegg	abundance
C00031	1.45
C00118	-0.83
C00022	2.11
```

Non-numeric columns after the first are reported and skipped rather than
silently misread.

---

## 🖼️ Render modes

| Mode | Draws | Needs KEGG PNG | Vector |
|---|---|:--:|:--:|
| `native` | data painted onto KEGG's own map image | yes | no |
| `vector` | the map redrawn from KGML coordinates | no | yes |
| `svg` | a standalone SVG document | no | yes |
| `graph` | a NetworkX node-link diagram | no | yes |
| `auto` | `native` if the PNG is present, else `vector` | no | depends |

`native` preserves KEGG's own labels — dark pixels are left untouched, so gene
symbols survive the overlay. `vector` is the one for figures that get resized,
and the only one that works with no KEGG image at all.

---

## 🎨 Colour scales

Default anchors are **R pathview's exact defaults** (`#00FF00` / `#BEBEBE` /
`#FF0000`), so output is directly comparable to an R figure. Binning
reproduces R's `cut(..., right = FALSE, include.lowest = TRUE)` and
`colorpanel2()` exactly, including the odd-*n* midpoint rule.

```python
from pathview import list_palettes
list_palettes()
# pathview, pathview_soft, rdbu, rdylbu, viridis, cividis,
# rnaseq, metabolite, bluered, tealrose, purpleorange
```

Values beyond the limits are **clamped, not dropped**. `discrete=True` is
honoured only when the limits are integral and the range divides evenly into
the bins — the same gate R applies, with a warning when it does not.

---

## 🗂️ Supported ID types

**Gene** — `ENTREZ`, `SYMBOL`, `UNIPROT`, `ENSEMBL`, `REFSEQ`, `KEGG`, `KO`,
`ALIAS`, `HGNC`, `MGI`, and more via `supported_gene_idtypes()`

**Compound** — `KEGG`, `NAME`, `CAS`, `CHEBI`, `HMDB`, `PUBCHEM`, `DRUGBANK`,
`LIPIDMAPS`, `CHEMBL`, `KNAPSACK`, and more via `supported_cpd_idtypes()`

**SBGN glyph ids** — Pathway Commons glyph ids are opaque hashes. 770k
crosswalk pairs ship in the wheel, and routing between any two types is
breadth-first over the crosswalk graph:

```python
from pathview import id_route, map_ids_to_sbgn, crosswalk_routes

id_route("ENTREZ", "SYMBOL")
# ['entrez', 'ko', 'pathwaycommons', 'symbol']

map_ids_to_sbgn(["1017"], "ENTREZ", "SYMBOL")
crosswalk_routes()          # every conversion the bundled table supports
```

---

## 🧬 Supported databases

| Database | Format | Coverage | Download |
|---|---|--:|---|
| **KEGG** | KGML + PNG | 10,718 organisms | `download_kegg` — REST API |
| **Reactome** | SBGN-ML | all, +1,749 offline | `download_reactome` — live exporter, collection fallback |
| **MetaCyc** | SBGN-ML | 2,518 | `download_metacyc` |
| **SMPDB** | SBGN-ML | 725 | `download_smpdb` |
| **PANTHER** | SBGN-ML | 152 | `download_panther` |
| **MetaCrop** | SBGN-ML | 62 | `download_metacrop` |

**5,206 pre-generated SBGN pathways** are indexed inside the wheel, so browsing
and searching work offline. The SBGN-ML itself (~690 MB in total) is fetched
per pathway on first use and cached.

None of PANTHER, MetaCyc or SMPDB publishes a per-pathway SBGN endpoint. The
pre-generated collection does — which is why these are now real downloads
rather than a manual step.

---

## 🏗️ Repository layout

```
pathview-plus/
├── lib/                      # package source — installs as `pathview`
│   ├── __init__.py           #   public API (155 exports)
│   ├── pathview.py           #   KEGG orchestrator
│   ├── sbgnview.py           #   SBGN orchestrator
│   ├── cli.py                #   command-line interface
│   │
│   ├── organisms.py          #   10,718-organism offline table
│   ├── constants.py  errors.py  utils.py  cache.py  bundled.py
│   │
│   ├── kgml_parser.py        #   KEGG KGML
│   ├── sbgn_parser.py        #   SBGN-ML (namespace- and port-aware)
│   ├── expansion.py          #   split_group / expand_node
│   │
│   ├── id_mapping.py         #   gene & compound ID conversion
│   ├── sbgn_hub.py           #   SBGN collection + crosswalks
│   ├── mol_data.py           #   aggregation, demo data
│   ├── node_mapping.py       #   data → nodes
│   ├── rdata.py              #   read R .RData without R
│   │
│   ├── color_mapping.py      #   ColorScale, R-parity binning
│   ├── layout.py             #   geometry, RasterFrame
│   ├── splines.py            #   Bezier / Catmull-Rom
│   │
│   ├── rendering.py          #   native raster overlay
│   ├── vector_rendering.py   #   publication vector renderer
│   ├── svg_rendering.py      #   standalone SVG
│   ├── graph_rendering.py    #   NetworkX graph view
│   ├── legend.py  highlighting.py  databases.py  parity.py  examples.py
│   │
│   └── data/                 #   bundled tables: organisms, compounds,
│                             #     crosswalks, SBGN index, demo data
├── bin/
│   └── pathview-cli.py       # direct-run launcher
├── tests/                    # 317 tests, all offline
├── docs/                     # Sphinx source
├── recipe/                   # Bioconda recipe
├── .github/workflows/ci.yml  # 4 Python versions × 3 operating systems
├── setup.py                  # shim; metadata lives in pyproject.toml
├── pyproject.toml            # maps lib/ → the `pathview` package
├── requirements.txt
├── PARITY.md                 # feature matrix vs the R packages
├── BUG_CHECKLIST.md          # 21 fixed bugs, each with a test id
└── CHANGELOG.md
```

> **Note on `lib/`** — the directory is `lib/` but the importable package is
> `pathview`, via `package-dir = { pathview = "lib" }` in `pyproject.toml`.
> This is the layout the project has always used. You always
> `import pathview`, never `import lib`.

```

---

## 🔧 API reference

```python
pathview(
    pathway_id,                    # str or sequence of str
    gene_data=None, cpd_data=None,
    species="hsa",
    gene_idtype="ENTREZ", cpd_idtype="KEGG",
    kegg_dir=".", out_dir=None, out_suffix="pathview",
    render_mode="auto",            # auto | native | vector | graph | svg
    output_format="png",           # png | pdf | svg
    theme="publication",           # publication | slate | dark
    gene_color=None, cpd_color=None,   # ColorScale, palette name, or dict
    limit=None, bins=10, both_dirs=True, discrete=False,
    node_sum="sum", rand_seed=None,
    split_group=False, expand_node=False,
    map_symbol=True, map_cpd_name=True, map_null=True,
    min_nnodes=3, quiet=False,
) -> PathwayResult | PathwayResultSet

sbgnview(
    pathway_id,                    # collection id, local .sbgn path, or sequence
    gene_data=None, cpd_data=None,
    gene_idtype="ENTREZ", cpd_idtype="KEGG",
    sbgn_dir=".", out_dir=None,
    show_compartments=True, show_processes=True,
    ...
) -> PathwayResult | PathwayResultSet
```

Full reference: **[pathview-plus.readthedocs.io](https://pathview-plus.readthedocs.io)**

Selected functions:

```python
# Species
get_species_code("human")            # SpeciesInfo
search_organisms("coli", limit=5)

# Parsing
parse_kgml(path); node_info(pw); pathway_edges(pw)
parse_sbgn(path); sbgn_to_df(pw); sbgn_edges(pw); sbgn_compartments(pw)

# SBGN collection
list_sbgn_pathways(source=None, query=None)
download_sbgn(pathway_id); sbgn_collection_info()
download_panther / download_metacyc / download_smpdb / download_metacrop

# Identifiers
id2eg(ids, "SYMBOL", org="hsa"); eg2id(ids, "SYMBOL")
cpd_id_map(ids, "CAS", "KEGG"); cpd_name_to_kegg(names)
map_ids_to_sbgn(ids, "ENTREZ", "pathwayCommons"); id_route(a, b)

# Expansion & data
split_groups(node_df); expand_nodes(node_df)
mol_sum(mol_data, id_map, sum_method="max_abs")
demo_gene_data(2); demo_cpd_data(); sim_mol_data("cpd")

# R interop
read_rdata("cpd.names.rda")
```

---

## ✅ Feature parity

**74 capabilities tracked: 73 full, 0 partial, 0 missing.**
Covers **97.0%** of pathview (R) and **98.1%** of SBGNview (R); **11**
capabilities are in neither.

The matrix lives in [`lib/parity.py`](lib/parity.py) and
`tests/test_parity.py` **fails the build if a feature claims support without a
working implementation**. Full table: [PARITY.md](PARITY.md), or
`pathview-plus parity --markdown`.

### The one thing marked not-applicable

**Bioconductor OrgDb annotation packages.** A Python package cannot import an R
library. The *conversions* OrgDb provides are covered by `id2eg`/`eg2id` and by
the bundled crosswalks — but if you already work in R, pathview's OrgDb
integration is more convenient than either.

### Caveats on features that do work

- **ENSEMBL has no offline route to SBGN glyph ids.** Two SBGNview crosswalk
  files use an ALTREP variant the bundled reader does not decode. ENSEMBL still
  converts through `id2eg` online.
- **Offline Entrez → Pathway Commons is sparse**, routing via KO (21k pairs).
  `gene_idtype="SYMBOL"` uses a direct 294k-pair route and maps far more
  glyphs. **Prefer SYMBOL for offline SBGN work.**
- **SBGN files are fetched, not bundled.** The index ships in the wheel; the
  SBGN-ML is downloaded per pathway on first use, then cached.
- **The collection is a subset** — 1,749 Reactome pathways, not all of them.
  The live exporter covers the rest.

---

## 🧪 Development

```bash
git clone https://github.com/raw-lab/pathview-plus
cd pathview-plus
pip install -e ".[dev,docs]"

pytest -q                          # 316 tests, all offline
ruff check lib tests
cd docs && sphinx-build -b html . _build/html
```

The suite makes **no network calls**. CI runs with `PATHVIEW_OFFLINE=1` across
4 Python versions × 3 operating systems, so a regression that reintroduces a
hidden request fails the build.

Every bug fixed in 3.x has a regression test that fails against the version
that had it — see [BUG_CHECKLIST.md](BUG_CHECKLIST.md).

---

## 🤝 Contributing

Contributions welcome. Areas that would help most:

1. **ENSEMBL crosswalks** — decode the remaining ALTREP variants in
   `lib/rdata.py` so ENSEMBL reaches SBGN glyph ids offline
2. **Edge routing** — A* pathfinding for splines around obstacles
3. **SBGN glyph shapes** — a fuller SBGN-PD shape vocabulary
4. **Layout** — automatic relayout for maps without coordinates
5. **Performance** — parallel batch rendering

See [CONTRIBUTING.md](CONTRIBUTING.md). Please open an issue before large
changes. We welcome contributions from other experts expanding features in
Pathview-plus, including the R and Python versions.

---

## 📄 License

**Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)**
— see [LICENSE](LICENSE).

Bundled and retrieved third-party data carries its own terms, listed in the
LICENSE file. In particular: KEGG is © Kanehisa Laboratories — academic use of
the REST API is free, commercial use requires a licence.

---

## 📚 Citing

If you are publishing results obtained using Pathview-plus, please cite:

- **Pathview-plus (pre-print):** Figueroa III JL, Brouwer CR, White III RA.
  2026. *Pathview-plus: unlocking the metabolic pathways from cells to
  ecosystems.* bioRxiv.

If you use the R versions, please cite:

- **Pathview (R):** Luo W, Brouwer C. 2013. *Pathview: an R/Bioconductor
  package for pathway-based data integration and visualization.* Bioinformatics
  29(14):1830–1831. [doi:10.1093/bioinformatics/btt285](https://doi.org/10.1093/bioinformatics/btt285)
- **SBGNview (R):** Dong X, Vegesna K, Brouwer C, Luo W. 2022. *SBGNview: data
  analysis, integration and visualization on all pathways.* Bioinformatics
  38(5):1473–1476. [doi:10.1093/bioinformatics/btab793](https://doi.org/10.1093/bioinformatics/btab793)

The pre-generated SBGN collection and identifier crosswalks are derived from
the [SBGNview / SBGNhub](https://github.com/datapplab/SBGNhub) project.

---

## 📞 Support

- **Issues:** [open an issue](https://github.com/raw-lab/pathview-plus/issues)
- **Email:** [Dr. Richard Allen White III](mailto:rwhit101@uncc.edu)
- **Lab:** [RAW Lab](https://github.com/raw-lab), UNC Charlotte

---

**Made with ❤️ for the pathway visualization community**
