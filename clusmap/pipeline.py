"""One-call pipeline driver from a YAML/JSON config (kept for batch use).

Interactive / agent front-ends call the individual functions instead.
"""
from __future__ import annotations

import pandas as pd

from .io import import_data, preprocess, extract_color_cat, set_working_directory
from .cluster import gen_mod
from .plot import bulk_hm, cluster_sample_stats
from .annotate import (compute_pseudo_bulk, pseudo_bulk_hm, sc_marker_hm,
                       celltype_selection, mod_GO)
from .analysis import (module_eigengenes, hub_genes, module_trait_correlation,
                       module_preservation, project_modules, module_report)
from .motif_run import motif_pipeline


def run_pipeline_from_config(config_path: str) -> dict:
    """Run the pipeline described by a YAML (or JSON) config.

    Returns a dict of intermediate objects (``rna_df``, ``state``, ``hm`` ...)
    so a caller/agent can keep working with them.
    """
    import yaml
    with open(config_path) as fh:
        config = (__import__("json").load(fh) if config_path.endswith(".json")
                  else yaml.safe_load(fh))

    ctx: dict = {}
    if "set_working_directory" in config:
        set_working_directory(**config["set_working_directory"])

    if "import_data" not in config:
        raise ValueError("import_data section is required in config.")
    rna_df = ctx["rna_df"] = import_data(**config["import_data"])

    col_cat = None
    if "extract_color_cat" in config:
        col_cat = extract_color_cat(rna_df, **config["extract_color_cat"])

    if "preprocess" in config:
        rna_df = ctx["rna_df"] = preprocess(rna_df, **config["preprocess"])

    state = None
    if "gen_mod" in config:
        state = ctx["state"] = gen_mod(rna_df, **config["gen_mod"])

    hm = None
    if "bulk_hm" in config:
        args = dict(config["bulk_hm"])
        if col_cat and "col_cat" not in args:
            args["col_cat"] = col_cat
        hm = ctx["hm"] = bulk_hm(rna_df, state, **args)

    if "cluster_sample_stats" in config:
        args = dict(config["cluster_sample_stats"])
        if hm is not None:
            args["hm"] = hm
        ctx["stats"] = cluster_sample_stats(rna_df, **args)

    pb_df = None
    if "compute_pseudo_bulk" in config:
        pb_df = ctx["pseudo_bulk_df"] = compute_pseudo_bulk(**config["compute_pseudo_bulk"])
    if "pseudo_bulk_hm" in config:
        args = dict(config["pseudo_bulk_hm"])
        if pb_df is not None and "pseudo_bulk_df" not in args:
            args["pseudo_bulk_df"] = pb_df
        ctx["pb_hm"] = pseudo_bulk_hm(rna_df, hm, **args)

    selected = None
    if "celltype_selection" in config:
        selected = ctx["selected_celltypes"] = celltype_selection(**config["celltype_selection"])
    if "sc_marker_hm" in config:
        args = dict(config["sc_marker_hm"])
        if selected is not None and "celltype" not in args:
            args["celltype"] = selected
        ctx["swarm"] = sc_marker_hm(rna_df, hm, **args)

    if "mod_GO" in config:
        mod_GO(**config["mod_GO"])

    if "motif_pipeline" in config:
        ctx["motif"] = motif_pipeline(**config["motif_pipeline"])

    eigengenes = None
    if "module_eigengenes" in config:
        eigengenes = ctx["eigengenes"] = module_eigengenes(rna_df, state, **config["module_eigengenes"])

    hub_df = None
    if "hub_genes" in config:
        args = dict(config["hub_genes"])
        if eigengenes is not None and "eigengenes" not in args:
            args["eigengenes"] = eigengenes
        hub_df = ctx["hub_genes"] = hub_genes(rna_df, state, **args)

    if "module_trait_correlation" in config:
        args = dict(config["module_trait_correlation"])
        trait_file = args.pop("trait_file")
        trait_sep = args.pop("trait_sep", "\t")
        trait_index_col = args.pop("trait_index_col", 0)
        traits = pd.read_csv(trait_file, sep=trait_sep, index_col=trait_index_col)
        ctx["trait_cor"], ctx["trait_fdr"] = module_trait_correlation(eigengenes, traits, **args)

    if "module_preservation" in config:
        args = dict(config["module_preservation"])
        test_args = dict(args.pop("test_import_data"))
        test_rna = import_data(**test_args)
        ctx["preservation"] = module_preservation(rna_df, state, test_rna, **args)

    if "project_modules" in config:
        args = dict(config["project_modules"])
        new_args = dict(args.pop("new_import_data"))
        new_rna = import_data(**new_args)
        ctx["projection"] = project_modules(new_rna, eigengenes, **args)

    if "module_report" in config:
        args = dict(config["module_report"])
        if eigengenes is not None and "eigengenes" not in args:
            args["eigengenes"] = eigengenes
        if hub_df is not None and "hub_df" not in args:
            args["hub_df"] = hub_df
        ctx["report"] = module_report(args.pop("outdir", "."), rna_df, state, **args)

    print("Pipeline completed.")
    return ctx
