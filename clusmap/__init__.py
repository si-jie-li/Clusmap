"""clusmap: bulk RNA-seq module discovery, visualisation and annotation.

Typical use
-----------
    import clusmap as cm
    rna = cm.import_data("counts.tsv")
    rna = cm.preprocess(rna)
    state = cm.gen_mod(rna, deepSplit=1)          # -> ModuleState
    hm = cm.bulk_hm(rna, state, col_cat=...)      # clustermap
    cm.mod_GO("HM_ModGene.csv", organism="Mouse")

``state`` (a :class:`ModuleState`) is the single source of truth and supports
``split`` / ``merge`` / ``reassign`` / ``gene_module`` edits that propagate to
every downstream output.
"""
from .io import (set_working_directory, import_data, preprocess,
                 extract_color_cat, TableStructure)
from .state import ModuleState
from .cluster import gen_mod
from .plot import bulk_hm, cluster_sample_stats, module_color_map
from .spatial import (SpatialDataset, import_spatial, from_adata, select_hvgs,
                      spatial_module_scores, assign_spots_to_modules,
                      plot_spatial_modules, plot_spatial_expression)
from .annotate import (compute_pseudo_bulk, pseudo_bulk_hm, sc_marker_hm,
                       celltype_selection, mod_GO)
from .motif import module_motif
from .motif_run import motif_pipeline, prepare_module_fastas, run_ame
from .analysis import (module_eigengenes, hub_genes, top_hubs,
                       module_trait_correlation, project_modules, module_report,
                       module_preservation)
from .pipeline import run_pipeline_from_config
from .launch import launch_app
from . import config

__version__ = "0.2.0"

__all__ = [
    "set_working_directory", "import_data", "preprocess", "extract_color_cat",
    "TableStructure", "ModuleState", "gen_mod", "bulk_hm", "cluster_sample_stats",
    "module_color_map", "SpatialDataset", "import_spatial", "from_adata",
    "select_hvgs", "spatial_module_scores", "assign_spots_to_modules",
    "plot_spatial_modules", "plot_spatial_expression",
    "compute_pseudo_bulk", "pseudo_bulk_hm", "sc_marker_hm",
    "celltype_selection", "mod_GO", "module_motif", "run_pipeline_from_config",
    "motif_pipeline", "prepare_module_fastas", "run_ame", "config",
    "module_eigengenes", "hub_genes", "top_hubs", "module_trait_correlation",
    "project_modules", "module_report", "module_preservation", "launch_app",
]
