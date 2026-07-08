"""Pipeline step definitions for the no-code GUI — Streamlit-free and testable.

Each step is a :class:`StepSpec` carrying its tunable :class:`ParamSpec`s (with
defaults, help text and a ``advanced`` flag) and a *pure* ``run(session, params)``
function. ``session`` is a plain dict carrying intermediate objects (``rna``,
``state``, ``hm`` …); ``run`` mutates it and returns an artifact dict the UI
renders: ``{"log": str, "df": DataFrame, "image": path, "html": path,
"text": str, "level": "success"|"info"|"warning"|"error"}``.

Keeping this layer free of ``streamlit`` means every step (and the curation
helpers) can be unit-tested head-less, and ``gui.py`` stays a thin widget loop.
The parameter list of every step mirrors ``config_example.yaml`` one-to-one
(including the ``*_args`` pass-through dicts).
"""
from __future__ import annotations

import contextlib
import io
import os
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

import clusmap as cm
from . import config


# --------------------------------------------------------------------------- #
# spec dataclasses
# --------------------------------------------------------------------------- #
@dataclass
class ParamSpec:
    name: str
    kind: str          # int|float|str|bool|enum|list_str|json|intnone|floatnone|path
    default: Any
    help: str = ""
    choices: Optional[List[str]] = None
    advanced: bool = False                                       # fold by default
    show_if: Optional[Callable[[Dict[str, Any]], bool]] = None   # (params)->bool


@dataclass
class StepSpec:
    key: str
    label: str
    help: str
    run: Callable[[Dict[str, Any], Dict[str, Any]], Dict[str, Any]]
    params: List[ParamSpec] = field(default_factory=list)
    needs: List[str] = field(default_factory=list)   # session keys required
    group: str = "core"
    default_on: bool = True


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _capture(fn, *a, **k):
    """Run fn, returning (result, captured_stdout)."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        res = fn(*a, **k)
    return res, buf.getvalue()


def _missing(session, needs):
    return [n for n in needs if session.get(n) is None]


def _outdir(session):
    return session.get("outdir", "clusmap_out")


def _is_hm(obj) -> bool:
    """True if ``obj`` is a real seaborn ClusterGrid (not None / a stray widget value)."""
    return obj is not None and hasattr(obj, "fig") and hasattr(obj, "ax_heatmap")


def _ensure_hm(S):
    """Guarantee ``S['hm']`` is a valid ClusterGrid, (re)rendering if needed."""
    if not _is_hm(S.get("hm")):
        _run_bulk_hm(S, S.get("_hm_params") or _defaults(STEP_BY_KEY["bulk_hm"]))


def _tri(v):
    """Map an 'auto'/'true'/'false' enum to None/True/False."""
    return None if v in (None, "auto") else (v in (True, "true"))


# --------------------------------------------------------------------------- #
# run functions (pure)
# --------------------------------------------------------------------------- #
def _run_import(S, p):
    kw = dict(index_col=p["index_col"])
    if p.get("file_delimiter"):
        kw["file_delimiter"] = p["file_delimiter"]
    if p.get("header_path"):
        kw["header_path"] = p["header_path"]
        kw["header_delimiter"] = p.get("header_delimiter", "\t")
    for k in ("str_col_num", "float_col_num"):
        if p.get(k) is not None:
            kw[k] = p[k]
    if p.get("layer"):
        kw["layer"] = p["layer"]
    rna, log = _capture(cm.import_data, p["file_path"], **kw)
    S["rna_raw"] = rna
    S["rna"] = rna           # so downstream works even if preprocess skipped
    return {"log": log, "df": rna.iloc[:5, :min(8, rna.shape[1])].round(2),
            "text": f"Imported {rna.shape[0]} genes × {rna.shape[1]} samples"}


def _run_preprocess(S, p):
    if _missing(S, ["rna_raw"]):
        return {"text": "Run Import data first.", "level": "info"}
    rna, log = _capture(cm.preprocess, S["rna_raw"], start=p["start"], end=p["end"],
                        min_expr=p["min_expr"], log_base=p["log_base"],
                        plus_constant=p["plus_constant"], std_filter=p["std_filter"])
    S["rna"] = rna
    return {"log": log, "text": f"{rna.shape[0]} genes × {rna.shape[1]} samples after preprocessing"}


def _run_gen_mod(S, p):
    if _missing(S, ["rna"]):
        return {"text": "Run Import (and Preprocess) first.", "level": "info"}
    # advanced dynamicTreeCut knobs — only forward the ones the user actually set
    extra = {k: p[k] for k in ("cutHeight", "maxCoreScatter", "minGap",
                               "maxAbsCoreScatter", "minAbsGap", "minSplitHeight",
                               "minAbsSplitHeight", "minExternalSplit", "maxPamDist")
             if p.get(k) is not None}
    for k in ("pamRespectsDendro", "useMedoids", "respectSmallClusters", "verbose"):
        if k in p:
            extra[k] = p[k]
    state, log = _capture(cm.gen_mod, S["rna"], deepSplit=p["deepSplit"],
                          minClusterSize=p["minClusterSize"], metric=p["metric"],
                          method=p["method"], pamStage=p["pamStage"],
                          save_raw=p["save_raw"], outdir=_outdir(S), **extra)
    S["state"] = state
    S["hm"] = None           # force a fresh heatmap render
    S["eigengenes"] = None
    return {"log": log, "text": f"{state.n_modules} modules · "
            f"{int((state.raw_labels == 0).sum())} genes unassigned"}


def _run_bulk_hm(S, p):
    if _missing(S, ["rna", "state"]):
        return {"text": "Run Generate modules first.", "level": "info"}
    S["_hm_params"] = p      # remembered so curation can re-render with the same look
    col_cat = cm.extract_color_cat(S["rna"], p["col_keywords"]) if p.get("col_keywords") else None
    hm_args = dict(p.get("hm_args") or {})
    hm_args.setdefault("cmap", p["cmap"])
    if p.get("show_sample_names"):
        hm_args["xticklabels"] = True
    hm, log = _capture(cm.bulk_hm, S["rna"], S["state"], norm_method=p["norm_method"],
                       title=(p["title"] or None), mod_palette=p["mod_palette"],
                       goi_list=(p["goi_list"] or None), goi_size=p["goi_size"],
                       col_cat=col_cat, col_palette=p["col_palette"],
                       col_color_manual=(p.get("col_color_manual") or None),
                       col_legend=list(p["col_keywords"]) if p.get("col_keywords") else None,
                       row_band=p["row_band"], xtickangle=p["xtickangle"],
                       yticksize=p["yticksize"],
                       mod_num_left=p["mod_num_left"], mod_num_right=p["mod_num_right"],
                       mod_num_font=(p.get("mod_num_font") or None),
                       row_band_font=(p.get("row_band_font") or None),
                       hm_args=hm_args, save_format=p["save_format"], outdir=_outdir(S))
    S["hm"] = hm
    ext = p["save_format"].lower()
    img = os.path.join(_outdir(S), f"heatmap.{ext}")
    out = {"log": log, "text": f"Heatmap rendered ({S['state'].n_modules} modules)."}
    if ext == "png":
        out["image"] = img
    return out


def _run_cluster_stats(S, p):
    if _missing(S, ["rna", "state"]):
        return {"text": "Run Generate modules first.", "level": "info"}
    out = _outdir(S)
    modgene = os.path.join(out, "HM_ModGene.csv")
    res, log = _capture(cm.cluster_sample_stats, S["rna"], modgene,
                        hm=S.get("hm") if _is_hm(S.get("hm")) else None,
                        module_col=p["module_col"], gene_col=p["gene_col"],
                        stats=(p["stats"] or None), save_heatmap=p["save_heatmap"],
                        outdir=out)
    df0 = next(iter(res.values())) if isinstance(res, dict) and res else None
    return {"log": log, "df": df0.round(2) if df0 is not None else None,
            "text": f"Per-cluster sample stats written under {out}/"}


def _run_markers(S, p):
    if _missing(S, ["rna", "state"]):
        return {"text": "Run Generate modules first.", "level": "info"}
    _ensure_hm(S)
    out = _outdir(S)
    modgene = os.path.join(out, "HM_ModGene.csv")
    selected, log1 = _capture(cm.celltype_selection, p["marker_path"], modgene,
                              marker_sep=p["marker_sep"],
                              marker_col_celltype=p["marker_col_celltype"],
                              marker_col_gene=p["marker_col_gene"],
                              module_col_gene=p["module_col_gene"],
                              module_col_module=p["module_col_module"],
                              p_val_threshold=p["p_val_threshold"],
                              max_celltypes_per_module=p["max_celltypes_per_module"],
                              hm_cmap=p["hm_cmap"],
                              heatmap_top_n_celltypes=p["heatmap_top_n_celltypes"],
                              outdir=out)
    S["selected_celltypes"] = selected
    if not selected:
        return {"log": log1, "level": "warning",
                "text": "No cell types enriched (markers may not match your gene IDs). "
                "Check the marker file's gene-symbol column."}
    _, log2 = _capture(cm.sc_marker_hm, S["rna"], S["hm"], p["marker_path"],
                       sep=p["marker_sep"], celltype_col=p["marker_col_celltype"],
                       gene_col=p["marker_col_gene"], celltype=selected,
                       title=p["title"], marker_show=p["marker_show"],
                       axhspan_args=(p.get("axhspan_args") or None),
                       swarm_args=(p.get("swarm_args") or None), outdir=out)
    return {"log": log1 + log2, "text": f"Selected cell types: {selected}",
            "image": os.path.join(out, "heatmap_swarm.png")}


def _run_pseudobulk(S, p):
    if _missing(S, ["rna", "state"]):
        return {"text": "Run Generate modules first.", "level": "info"}
    h5ad = p.get("h5ad_file") or S.get("_pb_h5ad")
    if not h5ad:
        return {"text": "Provide a single-cell .h5ad path (here or in the top panel).",
                "level": "info"}
    celltype_key = p.get("celltype_key") or S.get("_pb_celltype") or "celltype"
    _ensure_hm(S)
    out = _outdir(S)
    pb, log1 = _capture(cm.compute_pseudo_bulk, h5ad, celltype_key,
                        layer=(p["layer"] or None), chunk_size=(p["chunk_size"] or None),
                        log_base=p["log_base"], plus_constant=p["plus_constant"], outdir=out)
    S["pseudo_bulk_df"] = pb
    _, log2 = _capture(cm.pseudo_bulk_hm, S["rna"], S["hm"], pseudo_bulk_df=pb,
                       norm_method=p["norm_method"], xtickangle=p["xtickangle"],
                       left_gap=p["left_gap"], col_width=p["col_width"],
                       title=p["title"], title_size=p["title_size"],
                       pseudo_hm_args=(p.get("pseudo_hm_args") or None), outdir=out)
    return {"log": log1 + log2,
            "text": f"Pseudo-bulk: {pb.shape[0]} genes × {pb.shape[1]} cell types",
            "image": os.path.join(out, "pb_heatmap.png")}


def _run_go(S, p):
    if _missing(S, ["state"]):
        return {"text": "Run Generate modules first.", "level": "info"}
    out = _outdir(S)
    modgene = os.path.join(out, "HM_ModGene.csv")
    mods = [int(x) for x in p["modules"]] if p["modules"] else None
    _, log = _capture(cm.mod_GO, modgene, organism=p["organism"],
                      module_col=p["module_col"], gene_col=p["gene_col"],
                      library=(p["library"] or None),
                      GO_category=tuple(p["categories"]) if p["categories"] else ("BP", "MF"),
                      mod=mods, barplot=_tri(p["barplot"]), dotplot=_tri(p["dotplot"]),
                      barplot_args=(p.get("barplot_args") or None),
                      dotplot_args=(p.get("dotplot_args") or None),
                      enrichr_args=(p.get("enrichr_args") or None),
                      bg=S["rna"] if p["use_background"] and S.get("rna") is not None else None,
                      outdir=os.path.join(out, "GO"))
    return {"log": log, "text": f"GO results written under {out}/GO/"}


def _run_motif(S, p):
    if _missing(S, ["state"]):
        return {"text": "Run Generate modules first.", "level": "info"}
    import shutil
    out = _outdir(S)
    modgene = os.path.join(out, "HM_ModGene.csv")
    org = p["organism"]
    if org == "other":
        promoter, dbs = p["promoter_fasta"], ([p["motif_db"]] if p["motif_db"] else [])
    else:
        promoter, dbs = config.builtin_motif(org)
    if not promoter or not dbs:
        return {"level": "warning",
                "text": f"No motif references for '{org}'. Bundled human/mouse refs "
                "ship at publish time; for now select 'other' and give a promoter "
                "FASTA + a .meme motif DB, or set them with clusmap-config."}
    ame_ok = shutil.which("ame") is not None
    if p["run_mode"] == "local" and not ame_ok and not p["dry_run"]:
        return {"level": "warning",
                "text": "MEME-suite `ame` not found on PATH. Install MEME / use the "
                "clusmap Docker image, or set run_mode to slurm/ssh, or enable dry_run."}
    _, log = _capture(cm.motif_pipeline, modgene, outdir=os.path.join(out, "motif"),
                      promoter_fasta=promoter, motif_dbs=dbs, run_mode=p["run_mode"],
                      evalue=p["evalue"], min_genes=p["min_genes"],
                      wrap_mode=p["wrap_mode"], dry_run=p["dry_run"])
    note = "" if ame_ok else "\n(note: `ame` not detected — install MEME-suite to actually run)"
    return {"log": log + note, "text": f"Motif pipeline done ({p['run_mode']}, dry_run={p['dry_run']})."}


def _run_hub(S, p):
    if _missing(S, ["rna", "state"]):
        return {"text": "Run Generate modules first.", "level": "info"}
    out = _outdir(S)
    eig, _ = _capture(cm.module_eigengenes, S["rna"], S["state"], outdir=out)
    hub, log = _capture(cm.hub_genes, S["rna"], S["state"], eigengenes=eig,
                        top_n=p["top_n"], outdir=out)
    S["eigengenes"] = eig
    S["hub_df"] = hub
    top = cm.top_hubs(hub, top_n=p["top_n"])
    return {"log": log, "df": hub[hub["is_hub"]][["module", "gene", "kME"]].head(60),
            "text": f"Hub genes per module: { {k: v[:5] for k, v in top.items()} }"}


def _run_trait(S, p):
    if _missing(S, ["rna", "state"]):
        return {"text": "Run Generate modules first.", "level": "info"}
    out = _outdir(S)
    eig = S.get("eigengenes")
    if eig is None:
        eig, _ = _capture(cm.module_eigengenes, S["rna"], S["state"], outdir=out)
        S["eigengenes"] = eig
    if p["trait_keywords"]:
        traits = pd.DataFrame(cm.extract_color_cat(S["rna"], p["trait_keywords"]),
                              index=S["rna"].columns)
    elif p["trait_file"]:
        traits = pd.read_csv(p["trait_file"], sep=None, engine="python", index_col=0)
    else:
        return {"text": "Provide trait_keywords (from sample names) or a trait_file.",
                "level": "info"}
    (cor, _fdr), log = _capture(cm.module_trait_correlation, eig, traits, outdir=out,
                                save_heatmap=p["save_heatmap"])
    return {"log": log, "df": cor.round(2),
            "image": os.path.join(out, "module_trait_heatmap.png")}


def _run_preservation(S, p):
    if _missing(S, ["rna", "state"]):
        return {"text": "Run Generate modules first.", "level": "info"}
    if not p["test_file"]:
        return {"text": "Provide an independent expression matrix to test against.",
                "level": "info"}
    out = _outdir(S)
    test = cm.import_data(p["test_file"])
    if p["preprocess_test"]:
        test = cm.preprocess(test, min_expr=p["min_expr"], log_base=p["log_base"])
    pres, log = _capture(cm.module_preservation, S["rna"], S["state"], test,
                         n_perm=p["n_perm"], min_genes=p["min_genes"], seed=p["seed"],
                         outdir=out)
    return {"log": log, "df": pres, "image": os.path.join(out, "module_preservation.png")}


def _run_project(S, p):
    if _missing(S, ["rna", "state"]):
        return {"text": "Run Generate modules first.", "level": "info"}
    if not p["new_file"]:
        return {"text": "Provide a new dataset to project the modules onto.", "level": "info"}
    out = _outdir(S)
    eig = S.get("eigengenes")
    if eig is None:
        eig, _ = _capture(cm.module_eigengenes, S["rna"], S["state"], outdir=out)
        S["eigengenes"] = eig
    new = cm.import_data(p["new_file"])
    proj, log = _capture(cm.project_modules, new, eig, min_cor=p["min_cor"], outdir=out)
    assigned = int((proj["assigned_module"] != 0).sum())
    return {"log": log, "df": proj.head(60),
            "text": f"Projected {new.shape[0]} genes — {assigned} assigned to a module."}


def _run_report(S, p):
    if _missing(S, ["rna", "state"]):
        return {"text": "Run Generate modules first.", "level": "info"}
    out = _outdir(S)
    path, log = _capture(cm.module_report, out, S["rna"], S["state"],
                         eigengenes=S.get("eigengenes"), hub_df=S.get("hub_df"),
                         report_name=p["report_name"])
    return {"log": log, "html": path, "text": f"Report written: {path}"}


# --------------------------------------------------------------------------- #
# module curation (text-based, integrated into the heatmap block)
# --------------------------------------------------------------------------- #
def ordered_module_blocks(state) -> List[int]:
    """Distinct module ids in heatmap (dendrogram) order, excluding 0."""
    ordered = np.asarray(state.hm_labels)[state.order]
    blocks: List[int] = []
    for m in ordered:
        m = int(m)
        if (not blocks or blocks[-1] != m):
            blocks.append(m)
    return [m for m in blocks if m != 0]


def _rerender_hm(S):
    _run_bulk_hm(S, S.get("_hm_params") or _defaults(STEP_BY_KEY["bulk_hm"]))


def _hm_image(S):
    ext = (S.get("_hm_params") or {}).get("save_format", "png").lower()
    return os.path.join(_outdir(S), f"heatmap.{ext}")


def curate_merge(S, a: int, b: int):
    from .interactive import are_neighbors
    state = S["state"]
    ordered = np.asarray(state.hm_labels)[state.order]
    if int(a) == int(b):
        return {"text": "Pick two different modules to merge.", "level": "info"}
    if not are_neighbors(ordered, int(a), int(b)):
        return {"level": "warning",
                "text": f"Modules {a} and {b} are not neighbours in the heatmap — "
                "only adjacent module blocks can be merged."}
    state.merge(int(a), int(b))
    _rerender_hm(S)
    return {"text": f"Merged modules {a} + {b} → {state.n_modules} modules now.",
            "image": _hm_image(S)}


def curate_split(S, m: int):
    state = S["state"]
    try:
        state.split(int(m))
    except Exception as e:                 # too-small / no-data
        return {"text": f"Cannot split module {m}: {e}", "level": "warning"}
    _rerender_hm(S)
    return {"text": f"Split module {m} → {state.n_modules} modules now.",
            "image": _hm_image(S)}


def curate_reassign(S, genes: List[str], target: int):
    state = S["state"]
    genes = [g for g in genes if g]
    if not genes:
        return {"text": "Enter at least one gene to reassign.", "level": "info"}
    known = [g for g in genes if state.gene_module(g) is not None]
    state.reassign(genes, int(target))
    _rerender_hm(S)
    miss = set(genes) - set(known)
    note = f" ({len(miss)} not found: {sorted(miss)})" if miss else ""
    return {"text": f"Moved {len(known)} gene(s) → module {target}{note}.",
            "image": _hm_image(S)}


def save_curation(S):
    state = S["state"]
    out = _outdir(S)
    state.save(out)
    _rerender_hm(S)
    return {"text": f"Saved module_state.pkl + module_state.json + HM_ModGene.csv + "
            f"heatmap to {out}/. Edits: {state.history or 'none'}.",
            "image": _hm_image(S)}


# --------------------------------------------------------------------------- #
# step registry  (params mirror config_example.yaml, incl. *_args dicts)
# --------------------------------------------------------------------------- #
_NORMS = ["z_score", "centralize", "min_max", "max", "none"]

STEPS: List[StepSpec] = [
    StepSpec("import", "1 · Import data",
             "Load the bulk expression matrix (genes × samples). Format is auto-"
             "detected (.h5ad/.tsv/.csv/.xlsx/.pkl); delimiter is sniffed and the leading "
             "annotation columns are detected. Override the column options only for "
             "non-standard files.", _run_import, group="core", params=[
        ParamSpec("file_path", "path", "", "Path to the expression file — or upload below."),
        ParamSpec("index_col", "int", 1, "1-based column to use as the gene index."),
        ParamSpec("file_delimiter", "str", "", "Column delimiter; blank = auto-sniff.", advanced=True),
        ParamSpec("header_path", "path", "", "Optional separate file with column (sample) names.", advanced=True),
        ParamSpec("header_delimiter", "str", "\t", "Delimiter used inside header_path.", advanced=True),
        ParamSpec("str_col_num", "intnone", None, "Number of leading non-numeric annotation columns; blank = auto.", advanced=True),
        ParamSpec("float_col_num", "intnone", None, "Number of trailing numeric sample columns; blank = auto.", advanced=True),
        ParamSpec("layer", "str", "", "For .h5ad only: which layer to read; blank = .X.", advanced=True),
    ]),

    StepSpec("preprocess", "2 · Preprocess",
             "Filter low-expression genes, log-transform, drop zero-variance genes. "
             "Set min_expr=0 / log_base blank to skip a step.", _run_preprocess,
             needs=["rna_raw"], group="core", params=[
        ParamSpec("min_expr", "float", 10.0, "Drop genes whose max expression < this (0 = keep all)."),
        ParamSpec("log_base", "intnone", 2, "Log base (e.g. 2); blank = no log transform."),
        ParamSpec("std_filter", "bool", True, "Drop genes with zero variance across samples."),
        ParamSpec("plus_constant", "float", 0.1, "Pseudocount added before log to avoid log(0).", advanced=True),
        ParamSpec("start", "intnone", None, "First sample column to keep; blank = all.", advanced=True),
        ParamSpec("end", "intnone", None, "Last sample column to keep; blank = all.", advanced=True),
    ]),

    StepSpec("gen_mod", "3 · Generate modules",
             "Hierarchical clustering + dynamic tree cut. deepSplit (0 coarse → 4 fine) "
             "and minClusterSize are the two knobs you'll usually touch.", _run_gen_mod,
             needs=["rna"], group="core", params=[
        ParamSpec("deepSplit", "int", 1, "Cut sensitivity: 0 = few large modules, 4 = many small."),
        ParamSpec("minClusterSize", "int", 30, "Minimum genes per module."),
        ParamSpec("metric", "enum", "correlation", "scipy pdist distance metric.",
                  choices=["correlation", "euclidean", "cosine"]),
        ParamSpec("method", "enum", "average", "scipy linkage method.",
                  choices=["average", "ward", "complete", "single"]),
        ParamSpec("pamStage", "bool", False, "Run the PAM stage of dynamic tree cut.", advanced=True),
        ParamSpec("save_raw", "bool", False, "Also save raw cutree dict + linkage (redundant; off by default).", advanced=True),
        ParamSpec("cutHeight", "floatnone", None, "Max dendrogram height to cut at; blank = auto.", advanced=True),
        ParamSpec("maxCoreScatter", "floatnone", None, "dynamicTreeCut maxCoreScatter; blank = default.", advanced=True),
        ParamSpec("minGap", "floatnone", None, "dynamicTreeCut minGap; blank = default.", advanced=True),
        ParamSpec("maxAbsCoreScatter", "floatnone", None, "dynamicTreeCut maxAbsCoreScatter; blank = default.", advanced=True),
        ParamSpec("minAbsGap", "floatnone", None, "dynamicTreeCut minAbsGap; blank = default.", advanced=True),
        ParamSpec("minSplitHeight", "floatnone", None, "dynamicTreeCut minSplitHeight; blank = default.", advanced=True),
        ParamSpec("minAbsSplitHeight", "floatnone", None, "dynamicTreeCut minAbsSplitHeight; blank = default.", advanced=True),
        ParamSpec("minExternalSplit", "floatnone", None, "dynamicTreeCut minExternalSplit; blank = default.", advanced=True),
        ParamSpec("pamRespectsDendro", "bool", True, "PAM assignments respect the dendrogram.", advanced=True),
        ParamSpec("useMedoids", "bool", False, "Use medoids instead of means in the PAM stage.", advanced=True),
        ParamSpec("maxPamDist", "floatnone", None, "Max distance for PAM assignment; blank = default.", advanced=True),
        ParamSpec("respectSmallClusters", "bool", True, "Keep small clusters from the initial cut.", advanced=True),
        ParamSpec("verbose", "int", 2, "dynamicTreeCut verbosity (0–4).", advanced=True),
    ]),

    StepSpec("bulk_hm", "4 · Clusterheatmap",
             "Draw the annotated module heatmap (robust colours + module band by default). "
             "After it renders you can curate modules just below.", _run_bulk_hm,
             needs=["rna", "state"], group="core", params=[
        ParamSpec("norm_method", "enum", "z_score", "Row normalisation for display.", choices=_NORMS),
        ParamSpec("title", "str", "", "Figure title; blank = auto."),
        ParamSpec("mod_palette", "enum", "hsv", "Module colour palette.",
                  choices=["hsv", "tab20", "Set3", "husl"]),
        ParamSpec("cmap", "enum", "turbo", "Heatmap colormap.",
                  choices=["turbo", "viridis", "RdBu_r", "magma"]),
        ParamSpec("show_sample_names", "bool", False, "Show every sample name on the x-axis (hm_args.xticklabels=True)."),
        ParamSpec("goi_list", "list_str", [], "Genes of interest to label on the right (comma-separated)."),
        ParamSpec("col_keywords", "json", {}, 'Column colour bands from sample-name keywords, '
                  'e.g. {"organ": ["thymus", "limb"]}.'),
        ParamSpec("row_band", "bool", True, "Draw the module colour band on the right."),
        ParamSpec("xtickangle", "int", 90, "x tick label rotation (degrees).", advanced=True),
        ParamSpec("yticksize", "int", 5, "Font size for gene (y) labels when shown.", advanced=True),
        ParamSpec("goi_size", "int", 10, "Font size for genes-of-interest labels.", advanced=True),
        ParamSpec("mod_num_left", "float", 0.01, "Left offset (axis fraction) for odd module numbers.", advanced=True),
        ParamSpec("mod_num_right", "float", 0.03, "Right offset (axis fraction) for even module numbers.", advanced=True),
        ParamSpec("col_palette", "enum", "tab20", "Palette for the column (sample) colour bands.",
                  choices=["tab20", "Set3", "hsv", "husl"], advanced=True),
        ParamSpec("save_format", "enum", "png", "Figure file format.", choices=["png", "pdf", "svg"], advanced=True),
        ParamSpec("hm_args", "json", {}, 'Passed to seaborn.clustermap, e.g. '
                  '{"vmin": -2, "vmax": 2, "col_cluster": false}.', advanced=True),
        ParamSpec("mod_num_font", "json", {}, 'Font kwargs for module numbers, e.g. {"fontsize": 5, "fontweight": "bold"}.', advanced=True),
        ParamSpec("row_band_font", "json", {}, 'Module-band appearance, e.g. {"band_width": 0.03, "alpha": 0.35}.', advanced=True),
        ParamSpec("col_color_manual", "json", {}, 'Override band colours: {band: {category: "#hex"}}.', advanced=True),
    ]),

    StepSpec("cluster_stats", "4b · Per-cluster sample stats",
             "Summarise each module's expression per sample (mean/median/std…) and draw "
             "a stat heatmap.", _run_cluster_stats, needs=["rna", "state"], group="core",
             default_on=False, params=[
        ParamSpec("stats", "list_str", [], "Stats to compute (comma-separated): mean, median, std, "
                  "iqr, mad, cv, min, max, n; blank = mean+std."),
        ParamSpec("module_col", "str", "hm_mod", "Module column in HM_ModGene.csv.", advanced=True),
        ParamSpec("gene_col", "str", "gene", "Gene column in HM_ModGene.csv.", advanced=True),
        ParamSpec("save_heatmap", "bool", True, "Also save a heatmap per statistic.", advanced=True),
    ]),

    StepSpec("markers", "5a · Cell-type markers (PanglaoDB)",
             "No single-cell data needed: enrich each module for PanglaoDB cell-type "
             "markers and draw a marker swarm plot beside the heatmap.", _run_markers,
             needs=["rna", "state"], group="annotate", default_on=False, params=[
        ParamSpec("marker_path", "path", "PanglaoDB_markers.tsv", "Path to the PanglaoDB marker TSV."),
        ParamSpec("marker_col_celltype", "str", "cell type", "Cell-type column in the marker file."),
        ParamSpec("marker_col_gene", "str", "official gene symbol", "Gene-symbol column in the marker file."),
        ParamSpec("max_celltypes_per_module", "int", 5, "Max significant cell types per module."),
        ParamSpec("p_val_threshold", "float", 0.05, "FDR threshold for selecting cell types.", advanced=True),
        ParamSpec("marker_sep", "str", "\t", "Delimiter of the marker file.", advanced=True),
        ParamSpec("module_col_gene", "str", "gene", "Gene column in HM_ModGene.csv.", advanced=True),
        ParamSpec("module_col_module", "str", "hm_mod", "Module column in HM_ModGene.csv.", advanced=True),
        ParamSpec("hm_cmap", "str", "Reds", "Colormap for the enrichment heatmap.", advanced=True),
        ParamSpec("heatmap_top_n_celltypes", "int", 30, "Cell types shown in the enrichment heatmap.", advanced=True),
        ParamSpec("title", "str", "Marker Gene Swarm Plot", "Swarm-plot title.", advanced=True),
        ParamSpec("marker_show", "bool", False, "Label individual marker genes on the y-axis.", advanced=True),
        ParamSpec("swarm_args", "json", {}, "Passed to seaborn.swarmplot, e.g. {\"size\": 4}.", advanced=True),
        ParamSpec("axhspan_args", "json", {}, "Passed to ax.axhspan (module background bands).", advanced=True),
    ]),

    StepSpec("pseudobulk", "5b · Single-cell pseudo-bulk",
             "With a single-cell .h5ad: compute per-cell-type pseudo-bulk and draw it "
             "beside the heatmap.", _run_pseudobulk, needs=["rna", "state"],
             group="annotate", default_on=False, params=[
        ParamSpec("h5ad_file", "path", "", "Path to the single-cell .h5ad (or set it in the top panel)."),
        ParamSpec("celltype_key", "str", "celltype", "obs column with cell-type labels."),
        ParamSpec("norm_method", "enum", "z_score", "Row normalisation for display.", choices=_NORMS),
        ParamSpec("log_base", "intnone", None, "Log base for pseudo-bulk; blank = leave as-is (set 2 for raw counts)."),
        ParamSpec("layer", "str", "", "Layer to use; blank = .X.", advanced=True),
        ParamSpec("chunk_size", "intnone", None, "Cells per chunk for big atlases; blank = load all.", advanced=True),
        ParamSpec("plus_constant", "float", 0.1, "Pseudocount before log.", advanced=True),
        ParamSpec("xtickangle", "int", 0, "Cell-type label rotation.", advanced=True),
        ParamSpec("left_gap", "float", 0.08, "Gap between the heatmap and this panel.", advanced=True),
        ParamSpec("col_width", "float", 0.1, "Per-cell-type column width.", advanced=True),
        ParamSpec("title", "str", "Pseudo Bulk Heatmap", "Panel title.", advanced=True),
        ParamSpec("title_size", "int", 10, "Panel title font size.", advanced=True),
        ParamSpec("pseudo_hm_args", "json", {}, "Passed to seaborn.heatmap.", advanced=True),
    ]),

    StepSpec("go", "6 · GO enrichment",
             "Over-representation GO analysis per module (Enrichr).", _run_go,
             needs=["state"], group="downstream", default_on=False, params=[
        ParamSpec("organism", "enum", "Mouse", "Organism for Enrichr libraries.",
                  choices=["Mouse", "Human", "Fly", "Yeast", "Fish", "Worm"]),
        ParamSpec("categories", "list_str", ["BP", "MF"], "GO categories (comma-separated): BP, MF, CC."),
        ParamSpec("modules", "list_str", [], "Restrict to these module ids; blank = all."),
        ParamSpec("use_background", "bool", True, "Use all heatmap genes as the statistical background."),
        ParamSpec("library", "list_str", [], "Explicit Enrichr libraries; blank = auto-pick latest GO_*.", advanced=True),
        ParamSpec("module_col", "str", "hm_mod", "Module column in HM_ModGene.csv.", advanced=True),
        ParamSpec("gene_col", "str", "gene", "Gene column in HM_ModGene.csv.", advanced=True),
        ParamSpec("barplot", "enum", "auto", "Draw bar plots (auto = only with background).",
                  choices=["auto", "true", "false"], advanced=True),
        ParamSpec("dotplot", "enum", "auto", "Draw dot plots (auto = only without background).",
                  choices=["auto", "true", "false"], advanced=True),
        ParamSpec("barplot_args", "json", {}, "Passed to gseapy.plot.barplot.", advanced=True),
        ParamSpec("dotplot_args", "json", {}, "Passed to gseapy.plot.dotplot.", advanced=True),
        ParamSpec("enrichr_args", "json", {}, "Passed to gseapy.enrichr.", advanced=True),
    ]),

    StepSpec("motif", "7 · Motif analysis (AME)",
             "TF-motif enrichment on module promoters via MEME-suite AME. Select an "
             "organism (human/mouse bundled at publish time) or 'other' for custom "
             "promoter FASTA + motif DB. Requires `ame` on PATH for local runs.",
             _run_motif, needs=["state"], group="downstream", default_on=False, params=[
        ParamSpec("organism", "enum", "mouse", "Built-in organism, or 'other' for custom paths.",
                  choices=["mouse", "human", "other"]),
        ParamSpec("promoter_fasta", "path", "", "Promoter FASTA (gene-named headers).",
                  show_if=lambda p: p.get("organism") == "other"),
        ParamSpec("motif_db", "path", "", "Motif database (.meme).",
                  show_if=lambda p: p.get("organism") == "other"),
        ParamSpec("run_mode", "enum", "local", "Where to run AME.", choices=["local", "slurm", "ssh"]),
        ParamSpec("evalue", "float", 1.0, "AME E-value report threshold."),
        ParamSpec("min_genes", "int", 50, "Skip modules with fewer promoters than this."),
        ParamSpec("dry_run", "bool", True, "Preview commands/jobs without executing AME."),
        ParamSpec("wrap_mode", "enum", "top_n", "How AME results are summarised.",
                  choices=["top_n", "p", "e", "pe"], advanced=True),
    ]),

    StepSpec("hub", "8 · Hub genes (eigengenes)",
             "Compute module eigengenes and rank each module's hub genes by kME "
             "(intramodular connectivity).", _run_hub, needs=["rna", "state"],
             group="downstream", default_on=False, params=[
        ParamSpec("top_n", "int", 20, "Top hub genes to flag per module."),
    ]),

    StepSpec("trait", "9 · Module–trait correlation",
             "Correlate module eigengenes with sample traits (FDR-corrected) and draw a "
             "heatmap.", _run_trait, needs=["rna", "state"], group="downstream",
             default_on=False, params=[
        ParamSpec("trait_keywords", "json", {}, 'Traits from sample-name keywords, '
                  'e.g. {"organ": ["thymus", "limb"]}.'),
        ParamSpec("trait_file", "path", "", "…or a samples×traits table (overridden by keywords)."),
        ParamSpec("save_heatmap", "bool", True, "Save the correlation heatmap.", advanced=True),
    ]),

    StepSpec("preservation", "10 · Module preservation",
             "Test whether the modules reproduce in an independent dataset (Zsummary: "
             "≥10 strong, ≥2 moderate, <2 none).", _run_preservation, needs=["rna", "state"],
             group="downstream", default_on=False, params=[
        ParamSpec("test_file", "path", "", "Independent expression matrix (replicate/cohort)."),
        ParamSpec("preprocess_test", "bool", True, "Preprocess the test matrix the same way."),
        ParamSpec("n_perm", "int", 200, "Permutations for the null distribution."),
        ParamSpec("min_expr", "float", 10.0, "min_expr for test preprocessing.", advanced=True),
        ParamSpec("log_base", "intnone", 2, "log_base for test preprocessing.", advanced=True),
        ParamSpec("min_genes", "int", 5, "Skip modules smaller than this.", advanced=True),
        ParamSpec("seed", "int", 0, "Random seed for the permutations.", advanced=True),
    ]),

    StepSpec("project", "10b · Project modules onto a new dataset",
             "Assign genes of a new dataset to the existing modules by eigengene "
             "correlation (needs shared sample columns).", _run_project,
             needs=["rna", "state"], group="downstream", default_on=False, params=[
        ParamSpec("new_file", "path", "", "New expression matrix to project."),
        ParamSpec("min_cor", "float", 0.5, "Minimum eigengene correlation to accept an assignment."),
    ]),

    StepSpec("report", "11 · HTML report",
             "Bundle module size, preservation, eigengene sparkline, hub genes and any "
             "GO/cell-type/motif results into one self-contained HTML.", _run_report,
             needs=["rna", "state"], group="downstream", default_on=False, params=[
        ParamSpec("report_name", "str", "module_report.html", "Output HTML file name.", advanced=True),
    ]),
]

STEP_BY_KEY = {s.key: s for s in STEPS}


def _defaults(step: StepSpec) -> Dict[str, Any]:
    return {pp.name: pp.default for pp in step.params}
