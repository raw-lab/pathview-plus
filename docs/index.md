# pathview-plus

**KEGG and SBGN pathway visualisation with multi-omics data overlay.**

A Python reimplementation of R's [pathview](https://bioconductor.org/packages/pathview/)
and [SBGNview](https://bioconductor.org/packages/SBGNview/), built for people who need
transcript *and* metabolite data on the same diagram, publication-quality vector output,
and a tool that keeps working behind a firewall.

```{code-block} python
from pathview import pathview, demo_gene_data, demo_cpd_data

res = pathview(
    "00020",                       # TCA cycle
    gene_data=demo_gene_data(2),   # RNA-seq log2FC, two conditions
    cpd_data=demo_cpd_data(),      # metabolite abundances
    species="human",               # codes, names and taxids all work
    render_mode="vector",          # no KEGG background image needed
    output_format="pdf",
)
print(res.summary())
```

## What makes it different

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} Two data classes, two scales
RNA-seq fold-changes and metabolite abundances get **independent colour scales
and independent keys**. A transcript at +2 and a metabolite at +2 are not the
same statement, and the figure says so.
:::

:::{grid-item-card} Works offline
Species resolution, compound naming, identifier cross-referencing and vector
rendering need **no network access**. 10,718 KEGG organisms ship in the wheel.
:::

:::{grid-item-card} Real vector output
The `vector` renderer draws the map from the KGML coordinates, so PDF and SVG
scale losslessly — no upscaled raster in your figure.
:::

:::{grid-item-card} Honest about gaps
The {doc}`parity` table is generated from the code, and the test suite fails
if a feature claims support without an implementation.
:::

::::

## Install

```bash
pip install pathview-plus
# or
conda install -c bioconda pathview-plus
```

```{toctree}
:maxdepth: 2
:caption: Guide

quickstart
concepts
recipes
cli
```

```{toctree}
:maxdepth: 2
:caption: Reference

api
parity
migration
```

## Citation

If you are publishing results obtained using Pathview-plus, please cite:

- **Pathview-plus (pre-print):** Figueroa III JL, Brouwer CR, White III RA.
  2026. *Pathview-plus: unlocking the metabolic pathways from cells to
  ecosystems.* bioRxiv.

If you use the R versions, please cite:

- Luo W, Brouwer C. *Pathview: an R/Bioconductor package for pathway-based data
  integration and visualization.* Bioinformatics 29(14):1830–1831, 2013.
- Dong X, Vegesna K, Brouwer C, Luo W. *SBGNview: data analysis, integration
  and visualization on all pathways.* Bioinformatics 38(5):1473–1476, 2022.

The pre-generated SBGN collection and identifier crosswalks are derived from
the [SBGNview / SBGNhub](https://github.com/datapplab/SBGNhub) project.

Licensed under CC BY-NC 4.0.
