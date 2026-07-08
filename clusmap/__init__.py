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
                 extract_color_cat)
from .state import ModuleState
from .cluster import gen_mod
from .plot import bulk_hm, cluster_sample_stats
from .annotate import (compute_pseudo_bulk, pseudo_bulk_hm, sc_marker_hm,
                       celltype_selection, mod_GO)

__version__ = "0.2.0"

__all__ = [
    "set_working_directory", "import_data", "preprocess", "extract_color_cat",
    "ModuleState", "gen_mod", "bulk_hm", "cluster_sample_stats",
    "compute_pseudo_bulk", "pseudo_bulk_hm", "sc_marker_hm",
    "celltype_selection", "mod_GO"
]
