"""Provider-neutral tool layer for the clusmap agent.

`TOOLS` is a list of JSON-Schema tool definitions (the format every major LLM
tool-calling API accepts). `ToolSession` holds the live pipeline state
(rna_df, ModuleState, clustermap) across calls and dispatches a tool name +
arguments to the matching clusmap function. Nothing here imports an LLM SDK, so
the same tool layer works with any provider.
"""
from __future__ import annotations

import json
import os
from typing import Any, Dict

import matplotlib
matplotlib.use("Agg")  # agent runs head-less; figures are written to disk

import clusmap as cm


# --------------------------------------------------------------------------- #
# tool schemas
# --------------------------------------------------------------------------- #
TOOLS = [
    {
        "name": "load_data",
        "description": "Load a bulk RNA-seq expression matrix from a file "
                       "(.h5ad, .tsv, .csv, .txt, .xlsx auto-detected). Call this first. "
                       "Returns the gene/sample counts and a small preview so you can confirm orientation.",
        "input_schema": {
            "type": "object",
            "properties": {
                "file_path": {"type": "string", "description": "Path to the expression file on disk."},
            },
            "required": ["file_path"],
        },
    },
    {
        "name": "preprocess",
        "description": "Filter low-expression genes, optionally log-transform, drop zero-variance genes. "
                       "Use defaults unless the user asks otherwise.",
        "input_schema": {
            "type": "object",
            "properties": {
                "min_expr": {"type": "number", "description": "Drop genes whose max expression <= this. Default 10; use 0 to skip."},
                "log_base": {"type": ["integer", "null"], "description": "Log base (e.g. 2). null to skip log transform. Default 2."},
            },
        },
    },
    {
        "name": "cluster",
        "description": "Hierarchically cluster genes and cut into modules. "
                       "deepSplit 0 (few large modules) .. 4 (many small modules) is the main knob; "
                       "if the user wants 'about N modules', start at deepSplit 1-2 and adjust after seeing the count.",
        "input_schema": {
            "type": "object",
            "properties": {
                "deepSplit": {"type": "integer", "minimum": 0, "maximum": 4, "description": "Sensitivity 0-4. Default 1."},
                "minClusterSize": {"type": "integer", "description": "Minimum genes per module. Default 30."},
                "metric": {"type": "string", "description": "Distance metric. Default 'correlation'."},
                "method": {"type": "string", "description": "Linkage method. Default 'average'."},
            },
        },
    },
    {
        "name": "make_heatmap",
        "description": "Render the module clustermap to the output directory. Optional annotations: a title, "
                       "genes-of-interest labels, and column color bands derived from keywords in sample names.",
        "input_schema": {
            "type": "object",
            "properties": {
                "title": {"type": "string"},
                "norm_method": {"type": "string", "description": "z_score|centralize|min_max|max|none. Default z_score."},
                "goi_list": {"type": "array", "items": {"type": "string"}, "description": "Genes to label on the right."},
                "col_keywords": {
                    "type": "object",
                    "description": "Column color bands, e.g. {'organ': ['thymus','limb']} — matched against sample names.",
                    "additionalProperties": {"type": "array", "items": {"type": "string"}},
                },
                "show_sample_names": {"type": "boolean", "description": "Show every sample name on the x-axis. Default false."},
            },
        },
    },
    {
        "name": "module_overview",
        "description": "Summarise current modules: how many, sizes, and how many genes are unassigned.",
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "search_gene",
        "description": "Report which module a gene belongs to.",
        "input_schema": {
            "type": "object",
            "properties": {"gene": {"type": "string"}},
            "required": ["gene"],
        },
    },
    {
        "name": "list_module_genes",
        "description": "List the genes in a given module id.",
        "input_schema": {
            "type": "object",
            "properties": {"module": {"type": "integer"}},
            "required": ["module"],
        },
    },
    {
        "name": "edit_modules",
        "description": "Curate modules: split one module, merge two neighbours, or reassign genes. "
                       "Re-render the heatmap afterwards if the user wants to see the effect.",
        "input_schema": {
            "type": "object",
            "properties": {
                "op": {"type": "string", "enum": ["split", "merge", "reassign"]},
                "module": {"type": "integer", "description": "For split: the module to split."},
                "module_a": {"type": "integer", "description": "For merge: first module."},
                "module_b": {"type": "integer", "description": "For merge: second (neighbouring) module."},
                "genes": {"type": "array", "items": {"type": "string"}, "description": "For reassign: genes to move."},
                "target": {"type": "integer", "description": "For reassign: destination module (0 = unassign)."},
            },
            "required": ["op"],
        },
    },
    {
        "name": "celltype_swarm",
        "description": "Annotate modules with PanglaoDB cell-type markers (no single-cell data needed): "
                       "select the most correlated cell types per module and draw a marker swarm plot beside the heatmap.",
        "input_schema": {
            "type": "object",
            "properties": {
                "marker_path": {"type": "string", "description": "Path to PanglaoDB marker tsv."},
                "organism": {"type": "string", "description": "Used only to note context; selection is enrichment-based."},
                "max_celltypes_per_module": {"type": "integer", "description": "Default 5."},
            },
            "required": ["marker_path"],
        },
    },
    {
        "name": "go_enrichment",
        "description": "Run GO over-representation analysis per module (or a subset of modules).",
        "input_schema": {
            "type": "object",
            "properties": {
                "organism": {"type": "string", "description": "Human|Mouse|... Default Mouse."},
                "modules": {"type": "array", "items": {"type": "integer"}, "description": "Subset of modules; omit for all."},
                "categories": {"type": "array", "items": {"type": "string"}, "description": "BP|MF|CC. Default [BP, MF]."},
                "use_background": {"type": "boolean", "description": "Use all heatmap genes as background. Default true."},
            },
        },
    },
    {
        "name": "hub_genes",
        "description": "Compute module eigengenes and rank each module's hub genes by intramodular "
                       "connectivity (kME). Returns the top hubs per module.",
        "input_schema": {
            "type": "object",
            "properties": {"top_n": {"type": "integer", "description": "Top hubs per module. Default 10."}},
        },
    },
    {
        "name": "module_preservation",
        "description": "Test whether the current modules reproduce in an independent dataset "
                       "(replicate/cohort). Returns a Zsummary per module (>=10 strong, >=2 moderate, <2 none).",
        "input_schema": {
            "type": "object",
            "properties": {
                "test_file": {"type": "string", "description": "Path to the independent expression matrix."},
                "n_perm": {"type": "integer", "description": "Permutations. Default 200."},
            },
            "required": ["test_file"],
        },
    },
    {
        "name": "module_report",
        "description": "Write a self-contained HTML report bundling each module's size, preservation, "
                       "eigengene sparkline, hub genes, and any GO / cell-type / motif results found in the output dir.",
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "save_state",
        "description": "Persist the current ModuleState + HM_ModGene.csv so downstream/analysis can reuse the curation.",
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "set_output_dir",
        "description": "Set the directory where all outputs (figures, tables, state) are written.",
        "input_schema": {
            "type": "object",
            "properties": {"outdir": {"type": "string"}},
            "required": ["outdir"],
        },
    },
]


# --------------------------------------------------------------------------- #
# session executor
# --------------------------------------------------------------------------- #
class ToolSession:
    """Holds pipeline state and runs a single tool call. One per conversation."""

    def __init__(self, outdir: str = "clusmap_out"):
        self.outdir = outdir
        self.rna = None
        self.state = None
        self.hm = None
        self.marker_path = None

    # -- dispatch --------------------------------------------------------- #
    def run(self, name: str, args: Dict[str, Any]) -> str:
        """Execute a tool; always returns a string (JSON or message) for the LLM."""
        try:
            handler = getattr(self, f"_t_{name}", None)
            if handler is None:
                return f"ERROR: unknown tool {name!r}"
            return handler(**(args or {}))
        except Exception as e:  # surface errors back to the model, don't crash the loop
            import traceback
            return f"ERROR running {name}: {e}\n{traceback.format_exc(limit=2)}"

    def _need_state(self):
        if self.state is None:
            raise RuntimeError("No clustering yet — call cluster first.")
        if self.hm is None and self.rna is not None:
            pass

    # -- tools ------------------------------------------------------------ #
    def _t_set_output_dir(self, outdir):
        self.outdir = outdir
        os.makedirs(outdir, exist_ok=True)
        return f"Output directory set to {outdir}"

    def _t_load_data(self, file_path):
        self.rna = cm.import_data(file_path)
        prev = self.rna.iloc[:3, :min(4, self.rna.shape[1])]
        out = {
            "genes": int(self.rna.shape[0]), "samples": int(self.rna.shape[1]),
            "sample_names": list(map(str, self.rna.columns[:12])),
            "preview": prev.round(2).to_dict(),
        }
        detected = self.rna.attrs.get("clusmap_detected")
        if detected is not None:
            out["detected"] = {
                "header_row": detected.header_row,
                "gene_col": detected.gene_col,
                "transposed": detected.transposed,
                "confidence": detected.confidence,
                "notes": detected.notes,
            }
        return json.dumps(out)

    def _t_preprocess(self, min_expr=10, log_base=2):
        if self.rna is None:
            raise RuntimeError("Load data first.")
        self.rna = cm.preprocess(self.rna, min_expr=min_expr, log_base=log_base)
        return json.dumps({"genes": int(self.rna.shape[0]), "samples": int(self.rna.shape[1])})

    def _t_cluster(self, deepSplit=1, minClusterSize=30, metric="correlation", method="average"):
        if self.rna is None:
            raise RuntimeError("Load (and ideally preprocess) data first.")
        self.state = cm.gen_mod(self.rna, deepSplit=deepSplit, minClusterSize=minClusterSize,
                                metric=metric, method=method, outdir=self.outdir)
        self.hm = None  # invalidate previous figure
        return json.dumps({
            "n_modules": self.state.n_modules,
            "unassigned_genes": int((self.state.raw_labels == 0).sum()),
            "deepSplit": deepSplit, "minClusterSize": minClusterSize,
        })

    def _t_make_heatmap(self, title=None, norm_method="z_score", goi_list=None,
                        col_keywords=None, show_sample_names=False):
        self._need_state()
        col_cat = cm.extract_color_cat(self.rna, col_keywords) if col_keywords else None
        hm_args = {"xticklabels": True} if show_sample_names else None
        self.hm = cm.bulk_hm(self.rna, self.state, title=title, norm_method=norm_method,
                             goi_list=goi_list, col_cat=col_cat,
                             col_legend=list(col_keywords) if col_keywords else None,
                             hm_args=hm_args, outdir=self.outdir)
        return f"Heatmap written to {os.path.join(self.outdir, 'heatmap.png')} " \
               f"({self.state.n_modules} modules)."

    def _t_module_overview(self):
        self._need_state()
        sizes = {int(m): len(self.state.module_genes(int(m)))
                 for m in sorted(set(self.state.hm_labels)) if m != 0}
        return json.dumps({"n_modules": self.state.n_modules, "sizes": sizes,
                           "unassigned": int((self.state.hm_labels == 0).sum())})

    def _t_search_gene(self, gene):
        self._need_state()
        m = self.state.gene_module(gene)
        if m is None:
            return f"'{gene}' is not in the dataset."
        return f"'{gene}' is in module {m}." if m else f"'{gene}' is unassigned (module 0)."

    def _t_list_module_genes(self, module):
        self._need_state()
        genes = self.state.module_genes(int(module))
        return json.dumps({"module": int(module), "n": len(genes), "genes": genes[:200]})

    def _t_edit_modules(self, op, module=None, module_a=None, module_b=None,
                        genes=None, target=None):
        self._need_state()
        if op == "split":
            self.state.split(int(module)); msg = f"Split module {module}."
        elif op == "merge":
            self.state.merge(int(module_a), int(module_b)); msg = f"Merged {module_a}+{module_b}."
        elif op == "reassign":
            self.state.reassign(genes or [], int(target)); msg = f"Moved {len(genes or [])} genes -> {target}."
        else:
            return f"Unknown op {op}"
        self.state.save(self.outdir)
        return f"{msg} Now {self.state.n_modules} modules. (Re-run make_heatmap to refresh the figure.)"

    def _t_celltype_swarm(self, marker_path, organism=None, max_celltypes_per_module=5):
        self._need_state()
        if self.hm is None:
            self._t_make_heatmap()
        self.marker_path = marker_path
        modgene = os.path.join(self.outdir, "HM_ModGene.csv")
        selected = cm.celltype_selection(marker_path, modgene,
                                         max_celltypes_per_module=max_celltypes_per_module,
                                         outdir=self.outdir)
        cm.sc_marker_hm(self.rna, self.hm, marker_path, celltype_col="cell type",
                        gene_col="official gene symbol", celltype=selected, outdir=self.outdir)
        return json.dumps({"selected_celltypes": selected,
                           "figure": os.path.join(self.outdir, "heatmap_swarm.png")})

    def _t_go_enrichment(self, organism="Mouse", modules=None, categories=None, use_background=True):
        self._need_state()
        modgene = os.path.join(self.outdir, "HM_ModGene.csv")
        cm.mod_GO(modgene, organism=organism, mod=modules,
                  GO_category=tuple(categories) if categories else ("BP", "MF"),
                  bg=self.rna if use_background else None, outdir=self.outdir)
        return f"GO results written under {self.outdir}/ (per category, per module)."

    def _t_hub_genes(self, top_n=10):
        self._need_state()
        eig = cm.module_eigengenes(self.rna, self.state, outdir=self.outdir)
        hub_df = cm.hub_genes(self.rna, self.state, eigengenes=eig, top_n=top_n, outdir=self.outdir)
        return json.dumps({"top_hubs": cm.top_hubs(hub_df, top_n=top_n),
                           "table": os.path.join(self.outdir, "hub_genes.tsv")})

    def _t_module_preservation(self, test_file, n_perm=200):
        self._need_state()
        test = cm.import_data(test_file)
        pres = cm.module_preservation(self.rna, self.state, test, n_perm=n_perm, outdir=self.outdir)
        return pres.to_csv(sep="\t", index=False)

    def _t_module_report(self):
        self._need_state()
        path = cm.module_report(self.outdir, self.rna, self.state)
        return f"HTML report written to {path}"

    def _t_save_state(self):
        self._need_state()
        self.state.save(self.outdir)
        return f"State + HM_ModGene.csv saved to {self.outdir}/ (edits: {self.state.history})."
