# pathview-plus — feature parity

Generated from `pathview.parity.FEATURE_MATRIX`. Regenerate with
`pathview-plus parity --markdown`.

Every feature marked **yes** is asserted by `tests/test_parity.py` to resolve
to a real, importable API, so this table cannot drift from the code.

## Summary

| | |
|---|---|
| Capabilities tracked | 74 |
| Full | 73 |
| Partial | 0 |
| Not implemented | 0 |
| Not applicable | 1 |
| Coverage vs pathview (R) | 32/33 (97.0%) |
| Coverage vs SBGNview (R) | 53/54 (98.1%) |
| Present in neither R package | 11 |

## The one thing that is not applicable

**Bioconductor OrgDb annotation packages.** A Python package cannot import an
R library, so this is marked `n/a` rather than `no`. The *conversions* OrgDb
provides are covered by `id2eg`/`eg2id` (via MyGene.info) and by the bundled
crosswalks, which work offline. If you are already working in R and
Bioconductor, R pathview's OrgDb integration remains more convenient.

## Honest limitations that are not feature gaps

These do not appear in the matrix because they are caveats on features that
work, not missing features:

- **ENSEMBL has no offline route to SBGN glyph ids.** Two of SBGNview's
  crosswalk files use an ALTREP variant the bundled reader does not decode.
  ENSEMBL still converts through `id2eg` when a network is available.
- **Offline Entrez → Pathway Commons is sparse.** It routes via KO
  (21k pairs). `gene_idtype="SYMBOL"` uses a direct 294k-pair route and maps
  far more glyphs; prefer it for offline SBGN work.
- **SBGN files are fetched, not bundled.** The 5,206-pathway index ships in
  the wheel and is searchable offline, but the SBGN-ML itself totals ~690 MB
  and is downloaded per pathway on first use, then cached.
- **The pre-generated collection is a subset.** It carries 1749
  Reactome pathways, not all of them; the live Reactome exporter covers the rest.

## The SBGN collection

| Source | Pathways | Downloader |
|---|--:|---|
| Reactome | 1,749 | `download_reactome` (live exporter, collection fallback) |
| MetaCyc | 2,518 | `download_metacyc` |
| SMPDB | 725 | `download_smpdb` |
| PANTHER | 152 | `download_panther` |
| MetaCrop | 62 | `download_metacrop` |
| **Total** | **5,206** | |

## Full matrix

| Category | Feature | pathview-plus | pathview (R) | SBGNview (R) | Notes |
|---|---|:--:|:--:|:--:|---|
| Species | KEGG organism code lookup | yes | yes | yes | 10,718 organisms bundled |
| Species | Lookup works offline | yes | partial | no | R pathview bundles korg but falls back to a network fetch for misses |
| Species | Common/scientific name and taxid input | yes | yes | partial |  |
| Species | Fuzzy match with suggestions | yes | no | no |  |
| Species | Refresh organism table from KEGG | yes | yes | n/a |  |
| Input | KEGG KGML parsing | yes | yes | partial |  |
| Input | SBGN-ML parsing | yes | no | yes |  |
| Input | SBGN namespaced documents | yes | n/a | yes | Reactome exports declare a default namespace |
| Input | SBGN arc port resolution | yes | n/a | yes |  |
| Input | SBGN state variables / clone markers | yes | n/a | yes |  |
| Input | KGML reaction -> edge synthesis | yes | yes | n/a |  |
| Input | Group/complex node splitting | yes | yes | n/a | R's split.group; expansion is recorded so a figure's provenance can be audited |
| Input | Multi-identifier node expansion | yes | yes | n/a | R's expand.node; sub-nodes tile the original box exactly and edges are remapped onto them |
| ID mapping | Entrez / KEGG gene ids | yes | yes | yes |  |
| ID mapping | Symbol, Ensembl, UniProt, RefSeq | yes | yes | yes | via MyGene.info instead of Bioconductor OrgDb |
| ID mapping | Symbol/Ensembl/UniProt conversion without R | yes | n/a | n/a | MyGene.info covers what OrgDb provides, with no R dependency |
| ID mapping | Bioconductor OrgDb packages themselves | n/a | yes | yes | An R library cannot be imported from Python; the conversions it provides are covered by id2eg/eg2id and the bundled crosswalks |
| ID mapping | SBGN glyph-id crosswalks | yes | n/a | yes | 770k pairs bundled: ChEBI, KEGG, Entrez, KO, symbol, compound name, Pathway Commons |
| ID mapping | Multi-hop identifier routing | yes | no | no | Breadth-first search over the crosswalk graph |
| ID mapping | Compound cross-references (CAS, ChEBI, HMDB, ...) | yes | yes | yes | 14k pairs bundled, offline |
| ID mapping | Compound name -> KEGG accession | yes | yes | yes | includes conjugate-base derivation |
| ID mapping | Batched / cached lookups | yes | partial | partial |  |
| ID mapping | Multi-probe aggregation | yes | yes | yes | 8 methods incl. max_abs and seeded random |
| Colour | Diverging scale with midpoint | yes | yes | yes |  |
| Colour | R colorpanel2-identical binning | yes | yes | n/a |  |
| Colour | Discrete levels | yes | yes | partial |  |
| Colour | Independent gene and compound scales | yes | yes | yes |  |
| Colour | Two colour keys on one figure | yes | partial | partial | R pathview draws a single key; SBGNview draws separate legends per file |
| Colour | Named colour-blind-safe palettes | yes | no | no |  |
| Colour | Value transform before binning | yes | yes | yes |  |
| Colour | Multi-condition node slicing | yes | yes | yes | one column per condition |
| Render | Overlay on KEGG PNG | yes | yes | n/a |  |
| Render | Vector map drawn from coordinates | yes | partial | yes | R pathview's graph mode uses Rgraphviz; this draws KEGG's own layout |
| Render | Standalone SVG output | yes | no | yes |  |
| Render | PDF output | yes | yes | yes |  |
| Render | PNG output | yes | yes | yes |  |
| Render | Graph/network view with edges | yes | yes | n/a |  |
| Render | KEGG edge subtype styling | yes | yes | n/a | 17 subtypes from KEGG's own table |
| Render | Spline / Bezier edge routing | yes | no | yes |  |
| Render | Compartment shading | yes | no | yes | Nested compartments shaded largest-first with decreasing opacity and labelled |
| Render | Renders with no network access | yes | no | no |  |
| Render | Themes | yes | no | partial |  |
| Render | Element legend | yes | yes | partial |  |
| Post | Composable modification (result + layer) | yes | no | yes |  |
| Post | Highlight nodes | yes | no | yes |  |
| Post | Highlight edges/arcs | yes | no | yes |  |
| Post | Highlight a path | yes | no | yes |  |
| Post | Change node labels | yes | no | yes |  |
| Post | Free-text annotation | yes | no | partial |  |
| Sources | KEGG download | yes | yes | yes |  |
| Sources | Reactome SBGN download | yes | no | yes |  |
| Sources | Pre-generated SBGN collection | yes | n/a | yes | 5,206 pathways indexed in the wheel; files fetched on demand and cached, so nothing unused is downloaded |
| Sources | Collection browsable offline | yes | n/a | partial | The index ships in the wheel; SBGNview.data must be installed in full (69 MB) to list its contents |
| Sources | PANTHER download | yes | n/a | yes | 152 pathways via the pre-generated collection |
| Sources | MetaCyc download | yes | n/a | yes | 2,518 pathways; BioCyc's own API needs a subscription |
| Sources | SMPDB download | yes | n/a | yes | 725 pathways; SMPDB publishes only bulk archives |
| Sources | MetaCrop download | yes | n/a | yes | 62 crop-plant metabolic pathways |
| Sources | Local SBGN files from any source | yes | n/a | yes | A hand-exported file parses identically to a downloaded one; there is no second code path |
| Sources | Pathway search by name | yes | no | yes |  |
| SBGN | Top-level SBGN render entry point | yes | n/a | yes | sbgnview() is to SBGN what pathview() is to KEGG |
| SBGN | Omics overlay on SBGN glyphs | yes | no | yes | Tries every identifier system the glyphs might use and keeps whichever lands on the map |
| SBGN | SBGN batch rendering | yes | n/a | yes |  |
| SBGN | Arc-class styling | yes | n/a | yes |  |
| SBGN | Process / operator glyph rendering | yes | n/a | yes |  |
| SBGN | Compartment-aware layout extent | yes | n/a | yes | The canvas is widened so shading is never clipped |
| Engineering | Typed errors | yes | no | no |  |
| Engineering | Mapping diagnostics returned | yes | partial | partial |  |
| Engineering | Disk cache with TTL | yes | partial | partial | R caches KGML files only |
| Engineering | Retry with backoff | yes | no | no |  |
| Engineering | Enforced offline mode | yes | no | no |  |
| Engineering | Command-line interface | yes | no | no |  |
| Engineering | Graph metrics | yes | no | no |  |
| Engineering | Multi-pathway batch in one call | yes | yes | yes | A failed pathway is recorded rather than aborting the batch, and modifiers broadcast across the set |
| Engineering | Reads R .RData without R | yes | n/a | n/a | XDR reader covering the vectors and lists reference data is published as; factors, closures and S4 raise rather than guess |
