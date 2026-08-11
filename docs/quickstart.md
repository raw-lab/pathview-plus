# Quickstart

## Install

```bash
pip install pathview-plus
```

Verify the installation and see what is bundled:

```bash
pathview-plus info
```

## Your first figure

`pathview()` needs a pathway id, some data, and a species.

```python
from pathview import pathview, demo_gene_data

res = pathview(
    "04110",                      # Cell cycle
    gene_data=demo_gene_data(1),  # real GSE16873 breast-cancer log2 ratios
    species="hsa",
    out_dir="figures",
)
print(res.summary())
# hsa04110 | genes: 110/111 nodes carry data; ... | -> hsa04110.pathview.png
```

## Bring your own data

The first column holds identifiers; every other numeric column is one
condition. Extra conditions split each node into vertical slices.

```python
import polars as pl

df = pl.DataFrame({
    "entrez":   ["1017", "1019", "595", "983"],
    "control":  [ 0.10,  -0.40,  1.20,  0.05],
    "treated":  [ 1.80,  -1.90,  2.40, -0.30],
})

res = pathview("04110", gene_data=df, species="human", limit=2.0)
```

If your identifiers are not Entrez, say so and they are converted for you:

```python
res = pathview("04110", gene_data=df, species="human",
               gene_idtype="SYMBOL")     # or ENSEMBL, UNIPROT, REFSEQ, ...
```

## Transcripts and metabolites together

This is the case the package is built around. The two data classes get
separate scales and separate keys.

```python
from pathview import pathview, demo_cpd_data, demo_gene_data

res = pathview(
    "00020",
    gene_data=demo_gene_data(1),
    cpd_data=demo_cpd_data(),
    species="human",
    limit={"gene": 2.0, "cpd": 1.5},     # different scales, different limits
    gene_palette="rnaseq",
    cpd_palette="metabolite",
    render_mode="vector",
    output_format="pdf",
)
```

Metabolomics platforms usually export names rather than accessions. Pass
`cpd_idtype="NAME"` and they are resolved offline, including conjugate-base
forms such as *Pyruvate* for *Pyruvic acid*:

```python
cpds = pl.DataFrame({
    "name": ["Pyruvate", "Citrate", "2-Oxoglutarate", "Succinate", "L-Malate"],
    "log2fc": [1.4, -0.8, 0.3, -1.6, 0.9],
})
res = pathview("00020", cpd_data=cpds, cpd_idtype="NAME", species="hsa")
```

## Working offline

Behind a firewall, put the KGML where pathview-plus expects it and use the
vector renderer:

```python
from pathview import pathview, set_offline

set_offline(True)          # or export PATHVIEW_OFFLINE=1
res = pathview("00020", gene_data=df, species="hsa",
               kegg_dir="/data/kgml",    # contains hsa00020.xml
               render_mode="vector")     # draws the map itself
```

Species lookup, compound naming and ID cross-referencing keep working: those
tables ship inside the package.

## From the command line

```bash
pathview-plus render 00020 \
    --species human \
    --gene-data rnaseq.csv \
    --cpd-data metabolites.csv \
    --render-mode vector \
    --output-format pdf \
    --limit gene=2,cpd=1
```
