# Concepts

## Colour scales

A {class}`~pathview.ColorScale` binds together everything that decides what
colour a value becomes: the limits, the number of bins, the anchor colours,
and how NaN is drawn.

```python
from pathview import ColorScale, gene_scale, compound_scale

sc = ColorScale(limit=2.0, bins=10, both_dirs=True, palette="rdbu")
sc.bounds()        # (-2.0, 2.0)
sc.map_values([-3, 0, 3])
```

Two constructors give sensible defaults per data class:

```python
g = gene_scale(limit=2.0)        # green -> grey -> red, "RNA-seq log2FC"
c = compound_scale(limit=1.0)    # blue -> white -> amber, "Metabolite log2FC"
```

The two are deliberately different so that a reader can tell at a glance which
key applies to which node.

### Binning matches R

Values are placed into `bins` discrete colours, cut at
`linspace(low, high, bins + 1)`, left-closed and right-open with the top bin
closed. Values beyond the limits are **clamped, not dropped**. This reproduces
R pathview's `cut(..., right = FALSE, include.lowest = TRUE)` exactly, and
`colorpanel2()` reproduces its ramp including the odd-`n` midpoint rule.

The default anchors are R's own — `#00FF00` / `#BEBEBE` / `#FF0000` — so
out-of-the-box output is directly comparable to an R figure. Softer and
colour-blind-safe alternatives are one keyword away:

```python
from pathview import list_palettes
list_palettes()
# pathview, pathview_soft, rdbu, rdylbu, viridis, cividis, rnaseq,
# metabolite, bluered, tealrose, purpleorange
```

### Discrete data

`discrete=True` is honoured only when the limits are integers and the range
divides evenly into the bins — the same gate R applies. Otherwise a warning
explains why, and the data is treated continuously.

## Render modes

| Mode | Draws | Needs KEGG PNG | Vector output |
|---|---|:--:|:--:|
| `native` | data painted onto KEGG's own map image | yes | no |
| `vector` | the map redrawn from KGML coordinates | no | yes |
| `svg` | a standalone SVG document | no | yes |
| `graph` | a NetworkX node-link diagram | no | yes |
| `auto` | `native` if the PNG is present, else `vector` | no | depends |

`native` gives the familiar KEGG look and preserves KEGG's own labels: dark
pixels are left untouched, so gene symbols survive the overlay. `vector` is
the one to use for figures that will be resized, and the only one that works
with no KEGG image at all.

## Identifier flow

```
your ids ──(id2eg / cpd_id_map)──> Entrez or KEGG accession
         ──(node_map)────────────> KGML entry ids
         ──(mol_sum)────────────-> one value per node
```

Each stage reports what it did. `detailed=True` returns a result object with
counts and the identifiers that failed to map, rather than leaving you to
guess why a map came out blank:

```python
from pathview import node_map
res = node_map(gene_df, node_data, "gene", detailed=True)
print(res.summary())        # 21/27 nodes carry data; 27/11911 input IDs used
res.unmapped_ids[:5]
```

## Aggregation

Several probes or transcripts often land on one node. `node_sum` chooses how
they combine:

`sum`, `mean`, `median`, `max`, `min`, `max_abs`, `random`, `first`.

`max_abs` keeps the most extreme value **with its sign**, which is usually what
you want for fold-changes. `random` takes `rand_seed` and is reproducible: it
sorts before sampling, so the result depends on the set of values rather than
on row order.

## Results are composable

`pathview()` returns a {class}`~pathview.PathwayResult`. It indexes like a
dict for backwards compatibility, and supports layered post-processing:

```python
from pathview import highlight_nodes, highlight_path, change_labels

annotated = (res
             + highlight_nodes(["1431", "3417"], color="#7C3AED", width=3)
             + highlight_path(["1431", "3417", "3418"], color="orange")
             + change_labels({"1431": "CS *"}))
annotated.save("figure_annotated.png")
```

Each `+` returns a new result; the original is untouched.
