# Changelog

## 3.1.0 — 2026-08-10

Closes every remaining feature gap against R pathview and R SBGNview.
Parity: **74 capabilities, 73 full, 0 missing** — 97.0% of pathview (R),
98.1% of SBGNview (R), 11 capabilities in neither.

### Added

- **Pre-generated SBGN collection.** 5,206 pathways indexed inside the
  wheel (Reactome 1,749, MetaCyc 2,518,
  SMPDB 725, PANTHER 152,
  MetaCrop 62). Browsing and searching work offline;
  files are fetched per pathway on first use and cached.
- **Real PANTHER, MetaCyc, SMPDB and MetaCrop downloaders**, replacing the
  2.x warn-and-return-None stubs. Every entry in `DATABASE_INFO` is now
  `available: True` with a working downloader.
- **`sbgnview()` / `sbgnview_batch()`** — the SBGN counterpart to `pathview()`.
- **770k identifier crosswalk pairs** bundled (ChEBI, KEGG, Entrez, KO, gene
  symbol, compound name, Pathway Commons), with breadth-first routing between
  any two via `map_ids_to_sbgn` and `id_route`.
- **Multi-pathway batch**: `pathview([...])` and `sbgnview([...])` return a
  `PathwayResultSet`. A failure is recorded, not raised; modifiers broadcast.
- **Compartment shading** — nested, labelled, opacity-ordered, with the canvas
  widened so shading is never clipped.
- **`split_group` / `expand_node`** — complexes split into subunits and
  multi-gene nodes expanded, area-preserving, with edges remapped.
- **`read_rdata`** — reads R `.RData`/`.rda` without R. Written to build the
  bundled data; exported because reading an `.rda` is a normal need.
- **`bundled.py`** — one reader for bundled tables, replacing four copies.

### Fixed

- Compound names containing an apostrophe broke the bundled-TSV load, because
  one of the four duplicated readers had CSV quoting enabled.
- Expanded sub-nodes did not tile their parent exactly (three items in a 2x2
  grid covered three quarters of the node).
- The RData reader dropped the `names` attribute on string vectors, which also
  misaligned the byte stream so every following object was silently mis-read.

### Changed

- `pathway_id` accepts a sequence as well as a string.
- Bioconductor OrgDb is now marked `n/a` rather than a gap: a Python package
  cannot import an R library, and the conversions it provides are covered.


## 3.0.0 — 2026-08-05

A rewrite. Every module was replaced. 18 confirmed bugs fixed, each with a
regression test that fails against 2.x.

### Fixed — correctness

- Species lookup no longer requires network access (10,718 organisms bundled)
- Species lookup accepts names, common names and taxonomy ids, not only codes
- `bgcolor` is metadata, not data — **this was making every node render solid
  red on every figure**
- `max_abs` and `random` aggregation work (both raised `AttributeError`)
- SBGN files with XML namespaces parse — previously 0 glyphs and 0 arcs for
  every real Reactome export
- SBGN arcs resolve through `<port>` elements
- Catmull-Rom splines no longer produce NaN
- Colour parsing accepts named colours; `highlight_nodes()` no longer crashes
  on its own default
- Compound nodes are drawn at the right position and the right size
- The graph view has edges
- Unmapped-identifier counts are never negative
- `--cpd-data` alone no longer crashes the CLI
- `list_reactome_pathways(species=...)` honours its argument
- Version numbers agree
- Highlights land on the nodes they mark (they were drawn in KGML
  coordinates onto a composed-figure raster, displacing every one)

### Removed

- `download_panther()` and `download_smpdb()` — warn-and-return-`None` stubs.
  `DATABASE_INFO` now records availability honestly and `download_pathway()`
  raises an error naming the manual route.

### Added

- `render_mode="vector"` — draws the map from KGML coordinates, works offline,
  true vector PDF/SVG
- Independent gene and compound colour scales with two colour keys
- `ColorScale` with R-`colorpanel2`-identical binning
- Standalone SVG renderer with real edges, shared marker defs and a colour key
- `PathwayResult` — composable with `+ highlight_nodes(...)`, still dict-like
- Typed exception hierarchy and mapping diagnostics
- `set_offline()` / `PATHVIEW_OFFLINE=1`
- Disk cache with TTL, retry and backoff; batched identifier lookups
- Full CLI: `render`, `species`, `search`, `download`, `legend`, `parity`, `info`
- `pathview.parity` — machine-readable feature matrix, enforced by tests
- Bundled offline data: organisms, compound names and synonyms, compound
  cross-references, KEGG edge subtypes, real GSE16873 demo expression data
- Themes, named colour-blind-safe palettes, graph metrics
- 220 tests, all offline; GitHub Actions CI across 4 Python versions and 3 OSes
- Sphinx documentation with a generated parity page

### Changed

- `kegg_native=True` → `render_mode="native"`
- `low`/`mid`/`high` dicts → `gene_color=` / `cpd_color=`
- `kegg_species_code()` returns a code string, as R's does
- Compound labels resolve to names rather than raw accessions
- Node labels shortened as R's `short.name=TRUE` does
- Default palette is R pathview's exact `#00FF00`/`#BEBEBE`/`#FF0000`
