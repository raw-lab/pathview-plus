# Migrating

## From pathview-plus 2.x

3.0.0 is a rewrite. The main entry point keeps its name and its common
arguments, and the returned object still indexes like the old dict, so many
scripts run unchanged. Where behaviour differs, it differs because 2.x was
wrong.

### Output will look different — and that is the point

In 2.x, `bgcolor` was missing from the metadata column list, so the string
`"#FFFFFF"` was treated as data, converted with `int(hex, 16)` to 16,777,215,
and clipped to the top of the colour scale. **Every node rendered solid red
regardless of the data.** If your 2.x figures were uniformly red, that is why,
and 3.0.0 figures are the corrected ones.

Compound nodes were also drawn at the vertically mirrored position and at
twice their correct radius.

### Renamed and changed

| 2.x | 3.0.0 | Why |
|---|---|---|
| `kegg_native=True` | `render_mode="native"` | five modes, not a boolean |
| `low=`, `mid=`, `high=` dicts | `gene_color=`, `cpd_color=` | one object per data class |
| returns `dict` | returns `PathwayResult` | `dict` has no `__add__`, so the documented `result + highlight_nodes(...)` API could never work |
| returns `{}` on failure | raises a typed error | `{}` made every failure look like a different failure |
| `kegg_species_code()` → `SpeciesInfo` | → `str` | matches R's `kegg.species.code()` |
| `download_panther`, `download_smpdb` | removed | they only warned and returned `None` |

`PathwayResult` supports `result["plot_data_gene"]` and `.get(...)`, so
dict-style access continues to work.

### Removed

`download_panther()` and `download_smpdb()` advertised support they did not
have. PANTHER, MetaCyc and SMPDB have no public per-pathway SBGN endpoint;
{data}`~pathview.DATABASE_INFO` now records that plainly and
`download_pathway()` raises an error naming the manual route. Local SBGN files
from those sources parse and render exactly like Reactome's.

### New

- `render_mode="vector"` — draws the map from KGML coordinates, works offline,
  produces true vector PDF and SVG
- Independent gene and compound colour scales with two colour keys
- `set_offline()` / `PATHVIEW_OFFLINE=1`
- Typed exceptions and mapping diagnostics
- A working CLI, including metabolomics-only runs

## From R pathview

Argument names map almost directly; dots become underscores.

| R pathview | pathview-plus |
|---|---|
| `gene.data`, `cpd.data` | `gene_data`, `cpd_data` |
| `pathway.id` | first positional argument |
| `kegg.native = TRUE` | `render_mode="native"` |
| `kegg.dir` | `kegg_dir` |
| `node.sum` | `node_sum` |
| `limit = list(gene=1, cpd=1)` | `limit={"gene": 1, "cpd": 1}` |
| `both.dirs`, `trans.fun` | `both_dirs`, `trans_fun` |
| `low`/`mid`/`high` lists | `gene_color=` / `cpd_color=` |
| `map.symbol`, `map.cpdname` | `map_symbol`, `map_cpd_name` |
| `min.nnodes` | `min_nnodes` |

R takes a matrix with row names; pathview-plus takes a DataFrame whose first
column holds the identifiers:

```r
# R
pathview(gene.data = mat, pathway.id = "04110", species = "hsa")
```

```python
# Python — ids are a column, not an index
pathview("04110", gene_data=df, species="hsa")
```

Default colours are R's exact defaults, so figures are directly comparable.

## From R SBGNview

```r
SBGNview(gene.data = d, input.sbgn = "P00001", output.file = "out")
```

```python
from pathview import parse_sbgn, sbgn_to_df, sbgn_edges, keggview_vector

pw = parse_sbgn("P00001.sbgn")
keggview_vector(sbgn_to_df(pw), sbgn_edges(pw), pathway_name="P00001")
```

The one substantial gap is `SBGNview.data`: SBGNview ships thousands of
pre-generated SBGN files and cross-database identifier crosswalks. That bundle
is not replicated. Reactome's SBGN exporter is supported directly, and SBGN
files obtained any other way parse identically.
