# libraries here

This directory is the **package source**. It is named `lib/` but installs and
imports as `pathview`:

```toml
# pyproject.toml
[tool.setuptools]
packages = ["pathview"]
package-dir = { pathview = "lib" }
```

So you always write `import pathview`, never `import lib`.

## Two entry points you can run directly

```bash
python -m pathview.test_all_features    # is my install working? (offline, ~2 s)
python -m pathview.examples --out figs  # runnable examples
```

`test_all_features.py` is a **smoke test, not the test suite**. The real suite
is `tests/` at the repository root — 323 tests covering every module, every
fixed bug and every parity claim.

## Module map

| Module | Role |
|---|---|
| `pathview.py` | KEGG orchestrator — the `pathview()` function |
| `sbgnview.py` | SBGN orchestrator — the `sbgnview()` function |
| `cli.py` | command-line interface |
| `organisms.py` | offline organism table, species resolution |
| `kgml_parser.py` / `sbgn_parser.py` | pathway file parsing |
| `expansion.py` | `split_group` / `expand_node` |
| `id_mapping.py` / `sbgn_hub.py` | identifier conversion and crosswalks |
| `mol_data.py` / `node_mapping.py` | aggregation, mapping data onto nodes |
| `color_mapping.py` | `ColorScale`, R-parity binning |
| `layout.py` / `splines.py` | geometry and curves |
| `rendering.py` | native raster overlay on KEGG's PNG |
| `vector_rendering.py` | publication vector renderer |
| `svg_rendering.py` / `graph_rendering.py` | SVG and graph views |
| `highlighting.py` | composable post-hoc modification |
| `databases.py` | pathway downloads |
| `parity.py` | feature matrix, enforced by `tests/test_parity.py` |
| `rdata.py` | read R `.RData` without R |
| `bundled.py` | one reader for the tables in `data/` |
| `data/` | bundled reference tables (organisms, compounds, crosswalks, SBGN index) |
