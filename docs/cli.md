# Command line

```bash
pathview-plus --help
```

Sub-commands: `render`, `species`, `search`, `download`, `legend`, `parity`,
`info`.

## render

```bash
pathview-plus render PATHWAY_ID [PATHWAY_ID ...] [options]
```

| Option | Meaning |
|---|---|
| `--species`, `-s` | code, name or taxid (`hsa`, `human`, `9606`) |
| `--gene-data` | CSV/TSV: identifiers in column 1, values after |
| `--cpd-data` | CSV/TSV: identifiers in column 1, values after |
| `--gene-idtype` | `ENTREZ`, `SYMBOL`, `ENSEMBL`, `UNIPROT`, `REFSEQ`, `KEGG` |
| `--cpd-idtype` | `KEGG`, `NAME`, `CAS`, `CHEBI`, `HMDB`, `PUBCHEM` |
| `--render-mode` | `auto`, `native`, `vector`, `graph`, `svg` |
| `--output-format`, `-f` | `png`, `pdf`, `svg` |
| `--limit` | `1.5`, or `gene=2,cpd=1` for separate scales |
| `--gene-palette` / `--cpd-palette` | named palettes |
| `--node-sum` | `sum`, `mean`, `median`, `max`, `min`, `max_abs`, `random`, `first` |
| `--theme` | `publication`, `slate`, `dark` |
| `--offline` | never attempt a network request |

Either data file may be given alone; a metabolomics-only run is fully
supported.

## Other commands

```bash
pathview-plus species 'Mus musculus'      # resolve to a KEGG code
pathview-plus search coli --limit 5       # browse the bundled table
pathview-plus download hsa04110 -o kgml   # fetch without rendering
pathview-plus legend --out legend.png     # diagram element legend
pathview-plus parity --markdown           # feature matrix vs the R packages
pathview-plus info                        # what is installed and bundled
```
