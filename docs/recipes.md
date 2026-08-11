# Recipes

## Multi-condition time course

Every numeric column becomes one vertical slice of each node, in column order.

```python
import polars as pl
from pathview import pathview

df = pl.read_csv("timecourse.csv")     # entrez, t0, t6, t12, t24
res = pathview("04110", gene_data=df, species="hsa",
               limit=2.0, render_mode="vector", output_format="pdf",
               subtitle="0, 6, 12, 24 h after treatment")
```

## Metabolomics only

No gene data required — a metabolic map with only metabolite abundances is a
legitimate figure, and works from the CLI too.

```bash
pathview-plus render 00020 --cpd-data metabolites.csv --cpd-idtype NAME \
    --render-mode vector --cpd-palette metabolite
```

## Non-model organism

Anything in KEGG resolves, by code, name or taxonomy id:

```python
from pathview import pathview, search_organisms

search_organisms("tuberculosis", limit=3)
res = pathview("00020", gene_data=df, species="Mycobacterium tuberculosis")
```

## Reactome SBGN

```python
from pathview import download_reactome, parse_sbgn, sbgn_edges, sbgn_to_df

path = download_reactome("R-HSA-109582", output_dir="sbgn")
pw = parse_sbgn(path)
nodes, edges = sbgn_to_df(pw), sbgn_edges(pw)
```

SBGN files from any source parse the same way, so a PANTHER or MetaCyc export
downloaded by hand renders exactly like a Reactome one:

```python
pw = parse_sbgn("downloads/P00001.sbgn")
```

## Pathway topology

```python
from pathview import build_graph, pathway_metrics, node_info, parse_kgml, pathway_edges

pw = parse_kgml("hsa04110.xml")
G = build_graph(node_info(pw), pathway_edges(pw))
pathway_metrics(G)
# {'nodes': 115, 'edges': 79, 'density': ..., 'hubs': [('CDK2', 12), ...]}
```

## Batch rendering

```python
from pathview import pathview
from pathview.errors import PathviewError

for pid in ["00010", "00020", "00030", "04110"]:
    try:
        pathview(pid, gene_data=df, species="hsa", out_dir="figures",
                 render_mode="vector", quiet=True)
    except PathviewError as exc:
        print(f"{pid}: {exc}")
```

Or in one command:

```bash
pathview-plus render 00010 00020 00030 04110 --gene-data rna.csv -o figures
```

## Custom colour scale

```python
from pathview import ColorScale, pathview

sc = ColorScale(limit=(-3, 3), bins=12, low="#2166AC",
                mid="#F7F7F7", high="#B2182B", label="log2 FC (tumour/normal)")
res = pathview("04110", gene_data=df, species="hsa", gene_color=sc)
```

## Reproducible random aggregation

```python
res = pathview("04110", gene_data=df, species="hsa",
               node_sum="random", rand_seed=42)
```
