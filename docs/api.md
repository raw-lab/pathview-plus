# API reference

```{eval-rst}
.. currentmodule:: pathview
```

## Main entry point

```{eval-rst}
.. autofunction:: pathview
.. autoclass:: PathwayResult
   :members:
```

## Species

```{eval-rst}
.. autofunction:: get_species_code
.. autofunction:: kegg_species_code
.. autofunction:: search_organisms
.. autofunction:: list_organisms
.. autofunction:: organism_count
.. autofunction:: refresh_organism_table
.. autoclass:: SpeciesInfo
   :members:
```

## Parsing

```{eval-rst}
.. autofunction:: parse_kgml
.. autofunction:: node_info
.. autofunction:: pathway_edges
.. autofunction:: parse_sbgn
.. autofunction:: sbgn_to_df
.. autofunction:: sbgn_edges
.. autofunction:: arc_resolution_report
```

## Data and identifiers

```{eval-rst}
.. autofunction:: mol_sum
.. autofunction:: node_map
.. autofunction:: id2eg
.. autofunction:: eg2id
.. autofunction:: cpd_id_map
.. autofunction:: cpd_name_to_kegg
.. autofunction:: compound_name
.. autofunction:: demo_gene_data
.. autofunction:: demo_cpd_data
.. autofunction:: sim_mol_data
```

## Colour

```{eval-rst}
.. autoclass:: ColorScale
   :members:
.. autofunction:: gene_scale
.. autofunction:: compound_scale
.. autofunction:: node_color
.. autofunction:: colorpanel2
.. autofunction:: list_palettes
.. autofunction:: draw_dual_key
```

## Rendering

```{eval-rst}
.. autofunction:: keggview_native
.. autofunction:: keggview_vector
.. autofunction:: keggview_svg
.. autofunction:: keggview_graph
.. autofunction:: draw_pathway
.. autofunction:: build_graph
.. autofunction:: pathway_metrics
.. autofunction:: kegg_legend
.. autofunction:: sbgn_legend
```

## Post-processing

```{eval-rst}
.. autofunction:: highlight_nodes
.. autofunction:: highlight_edges
.. autofunction:: highlight_path
.. autofunction:: change_labels
.. autofunction:: annotate
```

## Databases

```{eval-rst}
.. autofunction:: download_kegg
.. autofunction:: download_reactome
.. autofunction:: download_pathway
.. autofunction:: list_reactome_pathways
.. autofunction:: find_reactome_pathways
.. autofunction:: detect_database
```

## Geometry and curves

```{eval-rst}
.. autoclass:: NodeBox
   :members:
.. autoclass:: Extent
   :members:
.. autofunction:: catmull_rom_spline
.. autofunction:: route_edge_spline
.. autofunction:: points_to_bezier_path
```

## Errors

```{eval-rst}
.. autoexception:: PathviewError
.. autoexception:: SpeciesNotFoundError
.. autoexception:: PathwayNotFoundError
.. autoexception:: NetworkError
.. autoexception:: ParseError
.. autoexception:: MappingError
.. autoexception:: RenderError
```

## SBGN maps

```{eval-rst}
.. autofunction:: sbgnview
.. autofunction:: sbgnview_batch
.. autofunction:: sbgn_node_map
.. autofunction:: sbgn_compartments
```

## The SBGN collection

```{eval-rst}
.. autofunction:: list_sbgn_pathways
.. autofunction:: find_sbgn_pathway
.. autofunction:: sbgn_collection_info
.. autofunction:: download_sbgn
.. autofunction:: download_sbgn_batch
.. autofunction:: sbgn_url
.. autofunction:: download_panther
.. autofunction:: download_metacyc
.. autofunction:: download_smpdb
.. autofunction:: download_metacrop
```

## Identifier crosswalks

```{eval-rst}
.. autofunction:: map_ids_to_sbgn
.. autofunction:: id_route
.. autofunction:: crosswalk_routes
.. autofunction:: supported_sbgn_idtypes
.. autofunction:: sbgn_xref
```

## Expansion

```{eval-rst}
.. autofunction:: split_groups
.. autofunction:: expand_nodes
.. autoclass:: ExpansionResult
   :members:
```

## Batch results

```{eval-rst}
.. autoclass:: PathwayResultSet
   :members:
```

## Reading R data files

```{eval-rst}
.. autofunction:: read_rdata
.. autofunction:: rdata_objects
.. autofunction:: read_bundled_tsv
```
