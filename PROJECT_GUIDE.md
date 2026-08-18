# clusmap — Complete Project Guide

**Bulk RNA-seq module discovery → annotated clusterheatmaps → functional annotation**,
with an interactive editor and a conversational agent.

This guide is written for two audiences:

- **Part 2 — For users**: what the tool can do, with fully runnable examples of
  every analysis.
- **Part 3 — For developers**: a file-by-file map of the code so you know exactly
  where to revise or dive into a specific function.

There is also a short architecture overview in **Part 1** that both audiences
should read first.

> Everything in this repo can be run with the project conda env:
> `/opt/anaconda3/envs/clusmap/bin/python`
> (`import clusmap as cm` from the folder that contains `clusmap/`, or after
> `pip install -e ".[all]"`).

---

## Part 1 — What this project is

### 1.1 One-paragraph summary

clusmap clusters the genes of a **bulk RNA-seq expression matrix** into
co-expressed **modules** (hierarchical clustering + dynamic tree cut, the
WGCNA-style approach), draws an annotated **clusterheatmap**, and then helps you
interpret each module through four kinds of annotation:

1. **Cell-type annotation** — either PanglaoDB marker swarmplots (no
   single-cell data needed) or a **single-cell pseudo-bulk** panel.
2. **GO enrichment** — per-module over-representation analysis (Enrichr via
   gseapy).
3. **TF-motif enrichment** — MEME-suite `ame` on per-module promoter sequences
   (runs locally, on SLURM, or over SSH to a cluster).
4. **WGCNA-style module analyses** — eigengenes, hub genes (kME), module–trait
   correlation, module preservation across datasets, cross-dataset projection,
   and a self-contained HTML report.

Everything flows through **one `ModuleState` object**, so manual curation
(split / merge / reassign modules) propagates to every downstream output.

### 1.2 The pipeline flow

```
import_data (io.py) → preprocess (io.py) → gen_mod (cluster.py) → ModuleState
      → bulk_hm (plot.py) → heatmap + HM_ModGene.csv
      → annotate.py (pseudo_bulk_hm | sc_marker_hm + celltype_selection | mod_GO)
      → analysis.py (eigengenes, hub_genes, trait_correlation, preservation, projection, report)
      → motif_run.py (motif_pipeline → AME) → motif.py (module_motif wrap-up)
```

### 1.3 Three front-ends, one engine

| Front-end | Entry point | Use when |
|---|---|---|
| **Python API / batch** | `cm.import_data(...)` … or `cm.run_pipeline_from_config("config.yaml")` | Scripted, reproducible analyses |
| **No-code Streamlit app** | `/opt/anaconda3/envs/clusmap/bin/streamlit run clusmap/app.py` | Point-and-click; module curation in the browser |
| **Conversational agent** | `/opt/anaconda3/envs/clusmap/bin/python -m clusmap.agent` | "Just describe what you want" |

All three call the same engine functions; nothing is duplicated.

### 1.4 The core idea: `ModuleState` is the single source of truth

`clusmap/state.py::ModuleState` is a dataclass holding the gene list, the scipy
linkage matrix, the raw cluster labels, and an edit history. Every downstream
step reads from it, and the curation operations (`merge`/`split`/`reassign`)
mutate it **in place** — then `state.save(outdir)` rewrites `HM_ModGene.csv` +
`module_state.pkl` and every downstream step re-runs on the edited modules.

**Two label numberings coexist (critical to understand):**

- `raw_labels` — straight from `cutreeHybrid` (0 = unassigned).
- `hm_labels` — 1-based, renumbered in **heatmap (dendrogram) order**; this is
  what `HM_ModGene.csv`, the figures, and all annotation steps use.

When joining any table against another, confirm which numbering it uses.

---

## Part 2 — For users

### 2.1 Install & quick start

```bash
pip install clusmap                 # core: clustering + heatmaps + GO + motif wrap-up
pip install "clusmap[all]"          # + single-cell, interactive app, agent, xlsx
# local development from this repo:
/opt/anaconda3/envs/clusmap/bin/pip install -e ".[all]"
```

Five lines to go from a counts file to an annotated heatmap:

```python
import clusmap as cm
rna   = cm.import_data("counts.tsv")          # auto-detects .h5ad/.tsv/.csv/.xlsx/.pkl
rna   = cm.preprocess(rna)                    # min_expr filter + log2 + drop zero-variance
state = cm.gen_mod(rna, deepSplit=1)          # -> ModuleState (modules)
hm    = cm.bulk_hm(rna, state, outdir="out")  # clusterheatmap + HM_ModGene.csv
cm.mod_GO("out/HM_ModGene.csv", organism="Mouse")
```

### 2.2 Demo data included in this working copy

- `demo_out/demo_counts.tsv` — 286 genes × 24 samples (12 `thymus_*` + 12
  `limb_*` replicates). Handy for trying every example below.
- `demo_out/demo_replicate.tsv` — an independent replicate dataset for the
  **module preservation** test.
- `demo_out/` and `clusmap_out/` also contain **already-computed outputs** so you
  can see what each artifact looks like before generating your own.

> Note: `*.out`/`demo_out`/`clusmap_out` directories are git-ignored
> (`.gitignore` has `*_out/`), so these demo files live on this machine but will
> not appear in a fresh clone.

### 2.3 The complete catalog of analyses

| # | Analysis | Function | Output |
|---|---|---|---|
| 1 | Load any expression matrix | `cm.import_data` | genes×samples `DataFrame` |
| 2 | Filter / log-transform | `cm.preprocess` | cleaned `DataFrame` |
| 3 | Color bands from sample names | `cm.extract_color_cat` | `{band: [category…]}` |
| 4 | Cluster into modules | `cm.gen_mod` | `ModuleState` + `ModGene.csv` + `HM_ModGene.csv` + `module_state.pkl/.json` |
| 5 | Clusterheatmap | `cm.bulk_hm` | `heatmap.png/pdf/svg` + `HM_ModGene.csv` |
| 6 | Per-module × sample stats | `cm.cluster_sample_stats` | `cluster_sample_<stat>.tsv` + heatmap |
| 7 | Cell-type markers (no sc data) | `cm.celltype_selection` + `cm.sc_marker_hm` | `heatmap_swarm.*`, `celltype_selection/` |
| 8 | Single-cell pseudo-bulk | `cm.compute_pseudo_bulk` + `cm.pseudo_bulk_hm` | `pseudo_bulk.tsv`, `pb_heatmap.*` |
| 9 | GO enrichment | `cm.mod_GO` | `<GO cat>/module_<m>_GO.csv` (+ bar/dot plots) |
| 10 | TF-motif enrichment | `cm.motif_pipeline` + `cm.module_motif` | per-module FASTA, `AME_out/`, `motif_results/` |
| 11 | Eigengenes | `cm.module_eigengenes` | `module_eigengenes.tsv` |
| 12 | Hub genes (kME) | `cm.hub_genes`, `cm.top_hubs` | `hub_genes.tsv` |
| 13 | Module–trait correlation | `cm.module_trait_correlation` | `module_trait_cor.tsv` / `_fdr.tsv` + heatmap |
| 14 | Module preservation | `cm.module_preservation` | `module_preservation.tsv` + barplot |
| 15 | Project new dataset | `cm.project_modules` | `projected_modules.tsv` |
| 16 | HTML report | `cm.module_report` | `module_report.html` |
| 17 | Curate modules | `state.split/merge/reassign/save` | updated state + `HM_ModGene.csv` |
| 18 | Everything at once (YAML) | `cm.run_pipeline_from_config` | all of the above |

### 2.4 Concrete end-to-end examples

#### Example A — basic run: load → preprocess → cluster → heatmap → GO

```python
import clusmap as cm

rna = cm.import_data("demo_out/demo_counts.tsv")      # (1) auto-detect format
rna = cm.preprocess(rna)                              # (2) min_expr=10, log2, drop zero-var
state = cm.gen_mod(rna, deepSplit=1, minClusterSize=30, outdir="clusmap_out")  # (3) 5 modules on demo
print(state.n_modules, state.gene_module("prog0_g0")) # how many / where a gene lives

hm = cm.bulk_hm(rna, state, outdir="clusmap_out",     # (4) the clusterheatmap
                title="My experiment",
                norm_method="z_score",                # z_score|centralize|min_max|max|none
                goi_list=["prog0_g0"], goi_size=10,
                hm_args={"figsize": (14, 10), "col_cluster": False},
                row_band=True, save_format="png")

cm.mod_GO("clusmap_out/HM_ModGene.csv", organism="Mouse")   # (5) GO per module
```

Optional: color the top annotation bands by sample name keywords.

```python
col_cat = cm.extract_color_cat(rna, {"organ": ["thymus", "limb"]})  # auto from sample names
hm = cm.bulk_hm(rna, state, col_cat=col_cat, col_legend=["organ"],
                col_color_manual={"organ": {"thymus": "#e41a1c", "limb": "#377eb8"}})
```

#### Example B — cell-type annotation WITHOUT single-cell data (PanglaoDB markers)

```python
selected = cm.celltype_selection(
    "PanglaoDB_marker.tsv", "clusmap_out/HM_ModGene.csv",
    marker_col_celltype="cell type", marker_col_gene="official gene symbol",
    max_celltypes_per_module=5, outdir="clusmap_out")
# -> celltype_selection/p_val_table.tsv, module_top_celltypes.tsv,
#    selected_celltypes.txt, module_celltype_heatmap.pdf

cm.sc_marker_hm(rna, hm, "PanglaoDB_marker.tsv",
                celltype_col="cell type", gene_col="official gene symbol",
                celltype=selected, outdir="clusmap_out")  # -> heatmap_swarm.png
```

The bundled `PanglaoDB_marker.tsv` (8286 rows) matches these column names.

#### Example C — cell-type annotation WITH single-cell data (pseudo-bulk)

```python
pb = cm.compute_pseudo_bulk("atlas.h5ad", celltype_key="celltype",
                            log_base=2, chunk_size=20000)   # chunk_size for big atlases
cm.pseudo_bulk_hm(rna, hm, pseudo_bulk_df=pb, outdir="clusmap_out")  # -> pb_heatmap.png
```

`pseudo_bulk_hm` appends a cell-type heatmap to the **same figure** as the bulk
heatmap (`hm`), with the same row order, so you can compare module expression
against cell types.

#### Example D — WGCNA-style module analyses

```python
# (1) eigengenes  — first PC of each module (modules × samples)
eig = cm.module_eigengenes(rna, state, outdir="clusmap_out")

# (2) hub genes  — rank genes by kME = correlation with their module eigengene
hubs = cm.hub_genes(rna, state, eigengenes=eig, top_n=20, outdir="clusmap_out")
cm.top_hubs(hubs, top_n=10)          # {module: [top hub genes]}

# (3) module–trait correlation  (categorical traits are one-hot encoded, FDR-corrected)
traits = pd.DataFrame({"organ": ["thymus"]*12 + ["limb"]*12}, index=rna.columns)
cor, fdr = cm.module_trait_correlation(eig, traits, outdir="clusmap_out")  # + heatmap

# (4) preservation  — are the modules real in an independent dataset?
test = cm.import_data("demo_out/demo_replicate.tsv")
pres = cm.module_preservation(rna, state, test, n_perm=200, outdir="clusmap_out")
# Zsummary: >=10 strong, 2..10 moderate, <2 not preserved

# (5) project a NEW dataset's genes onto these modules (needs shared sample columns)
new = cm.import_data("new_cohort.tsv")
cm.project_modules(new, eig, min_cor=0.5, outdir="clusmap_out")

# (6) one self-contained HTML report bundling size + hubs + GO/celltype/motif/preservation
cm.module_report("clusmap_out", rna, state)
```

#### Example E — manual module curation (split / merge / reassign)

```python
state.merge(2, 3)                          # merge two modules (by heatmap id)
state.split(1)                             # split module 1 into two (re-clusters its genes)
state.reassign(["prog0_g0", "prog0_g1"], 4)  # move genes into module 4 (0 = unassign)
state.save("clusmap_out")                  # rewrite HM_ModGene.csv + module_state.pkl
state2 = cm.ModuleState.load("clusmap_out")  # reload later, e.g. in a new session
```

Any downstream step re-run after an edit automatically uses the edited modules.

#### Example F — TF-motif enrichment (genome & motif DBs live on your machine/HPC)

First configure paths once (the multi-GB reference data is **never** packaged):

```bash
clusmap-config set promoter_fasta /data/mm10_promoter_500bp.fa   # gene-named headers
clusmap-config set motif_db.JASPAR /db/JASPAR2026_vertebrates.meme
clusmap-config set run_mode slurm        # local | slurm | ssh
clusmap-config set ssh_host me@cluster   # only for run_mode=ssh
clusmap-config show
```

Then run the whole motif pipeline in one call:

```python
cm.motif_pipeline("clusmap_out/HM_ModGene.csv", outdir="motif_out",
                  run_mode="slurm", dry_run=True)      # dry_run previews the sbatch scripts
```

This splits `HM_ModGene.csv` → per-module promoter FASTAs (pure Python) → one
AME job per module (local `ame`, or `sbatch`, or SSH+`sbatch`). For
`slurm`/`ssh` the jobs run asynchronously; collect results afterwards:

```python
cm.module_motif("motif_out/AME_out", "motif_out/motif_results",
                mode="top_n", top_n=10)   # modes: top_n | p | e | pe
```

No cluster? The bundled `Dockerfile` installs MEME-suite + bedtools + samtools,
so the motif step runs locally:

```bash
docker build -t clusmap .
docker run --rm -v $PWD:/work -w /work -v /data:/data clusmap \
    python -c "import clusmap as cm; cm.motif_pipeline('HM_ModGene.csv', run_mode='local')"
```

One-time genome/promoter prep scripts live in [`motif_analysis/`](motif_analysis)
(see Part 3.6).

#### Example G — batch mode: the whole workflow from one YAML

```python
ctx = cm.run_pipeline_from_config("config_example.yaml")
# ctx holds rna_df, state, hm, eigengenes, hub_genes, trait_cor, preservation,
# projection, report, motif ... whatever sections you included.
```

[`config_example.yaml`](config_example.yaml) is the **complete parameter
reference** — every tunable of every function at its default, with comments. A
step only runs if its top-level key is present; delete a section to skip it.

#### Example H — no-code app (Streamlit), with notebook hand-off

```bash
/opt/anaconda3/envs/clusmap/bin/streamlit run clusmap/app.py
```

One page, one block per pipeline step; each block shows every parameter with its
default as a gray placeholder and has a ▶ Run button. The top panel runs the
whole pipeline in one click. **Module curation is folded into the heatmap
block** — merge / split / reassign / find gene, then 💾 Save.

To skip file re-parsing (important for non-standard matrices), build `rna` +
`state` once in a notebook and hand them straight to the app:

```python
rna   = cm.import_data("weird_layout.tsv", str_col_num=2, index_col=2,
                       header_path="header.txt")   # whatever your file needs
rna   = cm.preprocess(rna)
state = cm.gen_mod(rna)
cm.launch_app(rna, state, outdir="results", port=8502, headless=True)  # remote/HPC: headless
```

`launch_app` returns the Streamlit subprocess (`.terminate()` to stop); on a
remote machine tunnel the port (`ssh -L 8501:localhost:8501 …`).

#### Example I — conversational agent (no YAML, no notebook)

```bash
export ANTHROPIC_API_KEY=sk-ant-...
/opt/anaconda3/envs/clusmap/bin/python -m clusmap.agent --outdir results
# other providers:  --provider <name> --api-key <key>
```

Then just type, e.g.: *"Cluster `mouse_thymus.tsv` into about 8 modules, draw
the heatmap coloured by organ, then run GO on the modules."*

Programmatic entry:

```python
from clusmap.agent import chat
chat(provider="anthropic", outdir="results")
```

The agent exposes: `load_data`, `preprocess`, `cluster`, `make_heatmap`,
`module_overview`, `search_gene`, `list_module_genes`, `edit_modules`,
`celltype_swarm`, `go_enrichment`, `hub_genes`, `module_preservation`,
`module_report`, `save_state`, `set_output_dir`. To add another LLM provider,
implement the `LLMBackend` protocol (see `OpenAIBackend` in
[`clusmap/agent/backends.py`](clusmap/agent/backends.py)).

### 2.5 All output files, explained

| File | What it is |
|---|---|
| `module_state.pkl` | **The full `ModuleState`** (linkage + labels + gene order + edit history). Reload with `cm.ModuleState.load(outdir)`. **Primary artifact.** |
| `module_state.json` | Human-readable summary (n_modules, metric, method, history) — read without unpickling. |
| `ModGene.csv` | gene → **raw** `cutreeHybrid` module id (0 = unassigned). |
| `HM_ModGene.csv` | gene → **heatmap** module id (`hm_mod`) + colour. **Feeds GO / motif / celltype.** |
| `heatmap.<png/pdf/svg>` | The annotated clusterheatmap. |
| `heatmap_swarm.*` | Heatmap + PanglaoDB marker swarmplot panel. |
| `pb_heatmap.*` | Heatmap + single-cell pseudo-bulk panel. |
| `celltype_selection/` | `p_val_table.tsv` (module × celltype FDR), `module_top_celltypes.tsv`, `selected_celltypes.txt`, `module_celltype_heatmap.pdf`. |
| `<GO cat>/module_<m>_GO.csv` (+ plots) | Per-module Enrichr GO results (e.g. `BP/module_1_GO.csv`). |
| `mod_fasta/Module<m>.fa` | Per-module promoter FASTA for motif analysis. |
| `AME_out/Module<m>/ame.tsv` | Raw MEME-suite `ame` output per module. |
| `motif_results/module_motif_*.tsv` | Wrapped-up motif hits (long format, summary, p/E matrices). |
| `module_eigengenes.tsv` | modules × samples, first PC per module. |
| `hub_genes.tsv` | gene, module, kME, hub_rank, is_hub. |
| `module_trait_cor.tsv` / `module_trait_fdr.tsv` / `module_trait_heatmap.png` | Module–trait correlation (one-hot encoded) + BH FDR + annotated heatmap. |
| `module_preservation.tsv` + `.png` | Zsummary per module (+ green/amber/red barplot). |
| `projected_modules.tsv` | gene, assigned_module, best_cor for a new dataset. |
| `module_report.html` | Self-contained report: size, preservation, eigengene sparkline, hub genes, GO/celltype/motif. |
| `cluster_sample_<stat>.tsv` (+ heatmap) | Per-module × sample summary stats. |
| `pseudo_bulk.tsv` | genes × cell types mean expression. |

### 2.6 Config & knobs you'll actually touch

- **`deepSplit`** (0–4) and **`minClusterSize`** are the two clustering knobs —
  coarse vs fine modules.
- **`norm_method`** for display: `z_score` (default) / `centralize` / `min_max` /
  `max` / `none`.
- **`*_args` dicts** (`hm_args`, `swarm_args`, `enrichr_args`, …) are passed
  straight through as `**kwargs` to seaborn / matplotlib / gseapy — clusmap does
  not validate them, so anything the library accepts works.
- Reference paths (genome FASTA, motif DBs) are resolved via `clusmap-config`
  and degrade gracefully (warn, not crash) when unset.
- GO enrichment uses **Enrichr over the network** (gseapy) — needs internet and
  `organism` spelled as Enrichr expects (`Mouse`, `Human`, `Fly`, …).

---

## Part 3 — For developers

### 3.1 Repository layout

```
Clusmap/
├── clusmap/                 # the package (the only thing installed)
│   ├── __init__.py          # public API re-exports + __version__
│   ├── io.py                # import_data / preprocess / extract_color_cat / TableStructure
│   ├── state.py             # ModuleState — the single source of truth
│   ├── cluster.py           # gen_mod — linkage + cutreeHybrid -> ModuleState
│   ├── plot.py              # bulk_hm, cluster_sample_stats
│   ├── annotate.py          # pseudo-bulk, swarm, celltype_selection, mod_GO
│   ├── analysis.py          # eigengenes, hubs, traits, preservation, projection, report
│   ├── motif.py             # module_motif — AME result wrap-up
│   ├── motif_run.py         # motif_pipeline / prepare_module_fastas / run_ame
│   ├── config.py            # ~/.clusmap/config.yaml system + clusmap-config CLI
│   ├── pipeline.py          # run_pipeline_from_config (YAML driver)
│   ├── util.py              # canon_norm / normalize_rows
│   ├── interactive.py       # Plotly figure + pure selection helpers (no Streamlit)
│   ├── app.py               # Streamlit entry (thin)
│   ├── gui.py               # Streamlit widget loop
│   ├── gui_steps.py         # spec-driven step definitions (Streamlit-free)
│   ├── launch.py            # launch_app — notebook hand-off
│   ├── py.typed             # PEP-561 marker
│   └── agent/               # conversational agent
│       ├── __init__.py
│       ├── tools.py         # JSON-Schema tools + ToolSession executor
│       ├── backends.py      # LLMBackend protocol + AnthropicBackend (+ OpenAI stub)
│       ├── runner.py        # chat() + CLI main
│       └── __main__.py
├── tests/                   # pytest (import auto-detection + table detection)
├── motif_analysis/          # HPC/genome prep scripts (source of truth for motif prep)
├── PanglaoDB_marker.tsv     # bundled marker reference (8286 rows)
├── config_example.yaml      # full parameter reference for the batch pipeline
├── demo_out/, clusmap_out/  # example data + already-computed outputs (git-ignored)
├── pyproject.toml           # packaging, extras, console scripts
├── Dockerfile               # clusmap + MEME-suite + bedtools + samtools image
├── MANIFEST.in, LICENSE, README.md, USAGE.md, CLAUDE.md
└── .gitignore               # ignores *_out/, tests/, CLAUDE.md
```

### 3.2 The core engine — file by file

#### [`clusmap/__init__.py`](clusmap/__init__.py) — public API surface
Re-exports 26+ names (see the `__all__` list). **Add a new public function here**
if it should be importable as `cm.<name>`. Also holds `__version__ = "0.2.0"`.

#### [`clusmap/state.py`](clusmap/state.py) — `ModuleState`
The heart of the package. A `@dataclass` holding `genes`, `linkage`,
`raw_labels`, `metric`, `method`, `history`, `_data` (expression values, needed
for `split()`).

- **Build**: `ModuleState.from_cutree(link, mod, genes, ..., data=rna_df)` — the
  only factory; `gen_mod` calls it.
- **Queries**: `order` (= `leaves_list(linkage)`, the heatmap row order — fully
  reproducible, independent of any drawing step), `ordered_genes`, `hm_labels`
  (1-based, renumbered in dendrogram order via `_renumber()`), `to_modgene_df()`
  (emits `HM_ModGene.csv`), `gene_module()`, `module_genes()`, `n_modules`.
- **Mutations (in place)**: `merge(a,b)`, `split(m)` (re-clusters the module's
  genes — **requires `_data`**, else raises), `reassign(genes, hm_id)`,
  each appending to `history`.
- **IO**: `save(outdir)` writes `module_state.pkl` + `module_state.json` +
  `HM_ModGene.csv`; `load(path)` unpickles.

**Where to touch**: to change how module ids are numbered, edit `_renumber`/
`hm_labels`; to add a new curation op, add a method here (then wire it into
`interactive.py`/`gui_steps.py`/`agent/tools.py`).

#### [`clusmap/io.py`](clusmap/io.py) — loading & preprocessing
The "just pass a path" promise.

- `import_data()` — detects format by extension (`.h5ad` via anndata → transposed
  to genes×samples; `.pkl`; `.xlsx`; delimited text). For text it either uses
  the **legacy path** (when `str_col_num`/`index_col`/`header_path`/`float_col_num`
  are pinned) or the **auto-detect path**.
- Auto-detection pipeline (the interesting part): `_sniff_delimiter` →
  `_numeric_mask` (handles Excel-date gene names like `1-Mar`, `SEPT2` so they
  aren't mistaken for numbers) → `_find_dense_block` (largest contiguous numeric
  rectangle, O(R²·C)) → `_refine_block_edges` → `_infer_orientation` (is it
  samples×genes and needs transposing?) → `_auto_detect_structure`. Result is a
  [`TableStructure`](clusmap/io.py#L78) dataclass + a note stored in
  `rna_df.attrs["clusmap_detected"]`.
- `preprocess()` — column slice, `min_expr` max-expression filter, `log_base`
  log transform with pseudocount, zero-variance drop. `min_expr=0`/`log_base=None`
  skip the respective step.
- `extract_color_cat()` — keyword matching in sample names → `{band: [cat…]}`
  aligned to `rna_df.columns`; feeds `bulk_hm(col_cat=...)`.

**Where to touch**: to support a new file format, add a branch in `import_data`;
to make auto-detection smarter, edit the `_find_dense_block` /
`_infer_orientation` heuristics (they are unit-tested in `tests/test_table_detect.py`).

#### [`clusmap/cluster.py`](clusmap/cluster.py) — `gen_mod`
The only entry point that builds a `ModuleState`. Deduplicates gene names
(appends `_dup`), computes `pdist` + `linkage`, calls
`dynamicTreeCut.cutreeHybrid` (forwarding any extra kwargs that match its
signature — so `cutHeight`, `maxCoreScatter`, etc. work), and builds the state
**with the expression data attached** (`data=rna_df`) so `split()` works.

Writes to `outdir`: `ModGene.csv` (raw ids), then `state.save()` →
`module_state.pkl` + `.json` + `HM_ModGene.csv`. With `save_raw=True` it also
dumps `modules.pkl` (raw cutree dict) and `linkage.pkl`.

#### [`clusmap/plot.py`](clusmap/plot.py) — heatmap & per-cluster stats
- `bulk_hm(rna_df, state, ...)` — seaborn `clustermap` with `row_linkage` from the
  state. **Normalization is applied by clusmap itself** (`normalize_rows`) rather
  than via clustermap's `z_score` flag so every figure uses the same transform.
  Adds: module-number labels alternating left/right of the heatmap
  (`_draw_module_numbers`, numbers in **heatmap order**, 0 unassigned skipped),
  a right-hand module colour band (`_draw_row_band`, attached as `hm.ax_band`),
  optional `col_cat` colour bands + legend, genes-of-interest labels, auto
  shrinking of x labels (`_shrink_xticklabels`), robust colour limits.
  Writes `heatmap.<fmt>` and `HM_ModGene.csv` (with colours).
- `cluster_sample_stats()` — groups genes by `hm_mod`, computes the requested
  stats per (module × sample), writes `<stat>.tsv` + a row-z-scored stat heatmap.

#### [`clusmap/annotate.py`](clusmap/annotate.py) — annotation
- `compute_pseudo_bulk()` — mean expression per cell type from an `.h5ad`;
  supports `chunk_size` **backed reads** for huge atlases. Writes `pseudo_bulk.tsv`.
- `pseudo_bulk_hm()` — appends a cell-type heatmap **onto the same figure** as
  the bulk clustermap (`hm.fig`), reusing its row order. Missing genes 0-filled.
- `sc_marker_hm()` — appends a PanglaoDB-style **swarmplot** panel, with the
  module-colour background bands and optional per-marker gene labels.
- `celltype_selection()` — **Fisher-exact** enrichment of cell-type markers per
  module, BH-FDR corrected; returns the selected cell types and writes the
  `celltype_selection/` folder.
- `mod_GO()` — per-module Enrichr GO via gseapy. Auto-picks the latest
  `GO_<BP|MF|CC>_<year>` library per category (or explicit `library=` list),
  honours a `bg` background universe (which switches barplot/dotplot defaults),
  and writes `<cat>/module_<m>_GO.csv` + plots.

#### [`clusmap/analysis.py`](clusmap/analysis.py) — WGCNA-style analyses
All built on `ModuleState`, no WGCNA dependency:
- `module_eigengenes()` — per module: z-score genes across samples, take the
  first right singular vector (SVD), sign-align to the module mean (WGCNA
  convention). Returns modules × samples.
- `hub_genes()` — kME = Pearson correlation of each gene with its eigengene;
  tidy table with `hub_rank` + `is_hub`. `top_hubs()` → `{module: [genes]}`.
- `module_trait_correlation()` — numeric traits direct, categorical one-hot
  (`get_dummies`), Pearson + **BH FDR**, writes cor/fdr TSVs + annotated heatmap
  (`*` `<.05`, `**` `<.01`, `***` `<.001`).
- `module_preservation()` — permutation-based (Langfelder-inspired compact
  reimplementation): density + connectivity stats of each module in a test set,
  z-scored against `n_perm` random gene sets → `Zsummary` (≥10 strong, ≥2
  moderate, <2 none) + barplot.
- `project_modules()` — assign each gene of a new dataset to the module whose
  eigengene it correlates with best (≥ `min_cor`); requires ≥3 shared sample
  columns.
- `module_report()` — builds one self-contained HTML (inline CSS + SVG
  sparklines), pulling in whatever GO/celltype/motif/preservation files exist in
  `outdir`.

#### [`clusmap/motif.py`](clusmap/motif.py) — AME result wrap-up
`module_motif(root_dir, outdir, mode=...)` — globs `root_dir/Module*/ame.tsv`,
parses/normalizes columns, filters by `top_n` / `p` / `e` / `pe`, and writes
`module_motif_longformat.tsv`, `module_motif_summary.tsv`, and
`module_motif_adjpval_matrix.tsv` / `module_motif_eval_matrix.tsv`.

#### [`clusmap/motif_run.py`](clusmap/motif_run.py) — motif driver
The per-run pipeline wired to `HM_ModGene.csv`:
- `prepare_module_fastas()` — **pure-Python** subsetting of a gene-named promoter
  FASTA into `mod_fasta/Module<m>.fa` (indexes by the part before `::` in each
  header, skips modules with < `min_genes` promoters, skips module 0).
- `run_ame()` — runs `ame` with a fixed command
  (`--control --shuffle-- --kmer 2 --method fisher --scoring totalhits
  --evalue-report-threshold <e> --noseq`) either `local` (subprocess, needs `ame`
  on PATH), `slurm` (writes + submits one `.slurm` script per module), or `ssh`
  (scp FASTAs, submit remotely). `dry_run` prints without executing.
- `motif_pipeline()` — the end-to-end driver. For `local` it returns the wrapped
  tables; for `slurm`/`ssh` it prints the collect command to run later.

#### [`clusmap/config.py`](clusmap/config.py) — reference-path config
Holds paths that must **not** ship in the package (genome FASTA, motif DBs) and
HPC settings. Resolution order: explicit arg → env `CLUSMAP_<KEY>` →
`~/.clusmap/config.yaml` → built-in default. `get("motif_db.JASPAR")` supports
nested keys. `builtin_motif(organism)` looks in `clusmap/data/motif/<org>/`
(bundled at publish time) then falls back to config. `require(*keys)` raises a
clear "set it with clusmap-config" error. The `_cli()` is the `clusmap-config`
console script (`show | get KEY | set KEY VALUE`).

#### [`clusmap/pipeline.py`](clusmap/pipeline.py) — batch driver
`run_pipeline_from_config(path)` — one function per top-level YAML key, executed
in order, with automatic wiring between steps (e.g. `col_cat` from
`extract_color_cat` → `bulk_hm`; `pseudo_bulk_df` → `pseudo_bulk_hm`;
`eigengenes` → `hub_genes`/`module_trait_correlation`/`project_modules`).
Returns a dict `ctx` of intermediate objects. A step only runs if its key is
present.

#### [`clusmap/util.py`](clusmap/util.py) — shared helpers
`canon_norm()` (validates + aliases `norm_method` names) and `normalize_rows()`
(the single row-wise normalizer used by **every** figure). If you change how
normalization works, this is the one place.

### 3.3 Interactive heatmap & GUI

#### [`clusmap/interactive.py`](clusmap/interactive.py) — Plotly figure + pure logic
Deliberately **Streamlit-free** so it can be unit-tested headlessly.
- `module_color_lut`, `ordered_frame` (normalized values in heatmap row order),
  `build_figure` (imports Plotly **lazily**, builds a 2-panel figure: module
  colour band + expression heatmap with robust colour limits, hover shows gene).
- Pure selection helpers: `rows_from_selection`, `genes_at_rows`,
  `modules_at_rows`, `are_neighbors` (module adjacency in heatmap order).

> Note: the older drag-to-select Plotly editor was **removed** (it clobbered
> `session_state["hm"]` — a seaborn `ClusterGrid` — with a Streamlit value and
> broke the app). Curation is now text-based only.

#### [`clusmap/app.py`](clusmap/app.py) — Streamlit entry (thin)
Sets page config, handles the `CLUSMAP_SESSION` hand-off bundle from
`launch_app` (loads `rna`/`state`/`outdir` once without re-parsing the file),
then calls `gui.render_pipeline()`.

#### [`clusmap/gui.py`](clusmap/gui.py) — widget loop (thin)
Renders each `StepSpec` as a bordered block: parameter widgets (default shown as
a gray placeholder; bool = toggle, enum = selectbox, everything else = text
input coerced back to the right type), a ▶ Run button, in-block output rendering
(dataframe / image / HTML download / captured console log), status icons, and a
closable notification dialog (`_note`/`_flush_notes`). The top control panel
(`_top_controls`) runs the whole pipeline in one click and toggles the optional
steps. Module curation UI (`_render_curation`) is folded into the heatmap block.

#### [`clusmap/gui_steps.py`](clusmap/gui_steps.py) — the real logic (Streamlit-free)
Defines `ParamSpec`/`StepSpec` dataclasses and the **`STEPS` registry** (14
steps, mirrored one-to-one against `config_example.yaml`). Every `run(session,
params)` is a **pure function** taking a plain `session` dict and returning an
artifact dict `{"log", "df", "image", "html", "text", "level"}` — this is what
makes the GUI testable without Streamlit. Also holds the curation helpers
`ordered_module_blocks`, `curate_merge`, `curate_split`, `curate_reassign`,
`save_curation`.

**Where to touch**: to add a new GUI step, (1) add a `StepSpec` to `STEPS`, (2)
write a pure `_run_*` function, (3) add the corresponding YAML section if it
should be batchable — the params mirror each other.

#### [`clusmap/launch.py`](clusmap/launch.py) — notebook hand-off
`launch_app(rna_df, state, hm=None, outdir=..., port=..., headless=...)` —
pickles a bundle to a temp dir (optionally a reference heatmap PNG), sets
`CLUSMAP_SESSION`, and `Popen`s `streamlit run app.py`. Only stdlib imported at
module load (Streamlit resolved lazily).

### 3.4 The agent

#### [`clusmap/agent/tools.py`](clusmap/agent/tools.py) — tools + executor
`TOOLS` is a list of **JSON-Schema** tool definitions (the format every major
LLM tool-calling API accepts — no LLM SDK import). `ToolSession` holds the live
state (`rna`, `state`, `hm`, `outdir`, `marker_path`) across calls and
dispatches `name`/`args` to `_t_<name>` handlers. Errors are returned to the
model as strings (not raised) so the loop keeps going.

#### [`clusmap/agent/backends.py`](clusmap/agent/backends.py) — LLM providers
`LLMBackend` protocol (`send(user_text, on_tool) -> str`). `AnthropicBackend`
implements the full tool loop with the anthropic SDK (`claude-opus-4-8`,
adaptive thinking, `max_tokens=8000`). `OpenAIBackend` is a **documented stub**
showing the contract for adding other providers. `make_backend()` picks by name.

#### [`clusmap/agent/runner.py`](clusmap/agent/runner.py) — CLI
`chat()` is a REPL (`input()` loop) that wires a `ToolSession` + backend; `main()`
parses `--provider/--api-key/--model/--outdir`. `__main__.py` delegates to
`main` (that's `python -m clusmap.agent`; the `clusmap-agent` console script
points at `runner:main`).

### 3.5 Tests

- `tests/test_table_detect.py` — unit tests for the `io` auto-detection
  internals (`_is_number_cell`, `_find_dense_block`, `_auto_detect_structure`)
  against the messy fixtures.
- `tests/test_import_data_detect.py` — end-to-end `import_data` tests (standard,
  junk rows before header, gene in col 2, transposed, numeric gene names, numeric
  sample names, legacy override path, and a `gen_mod` smoke test).
- `tests/fixtures.py` — the 6 CSV fixtures as raw strings + a `grid(name)`
  helper.

Run: `/opt/anaconda3/envs/clusmap/bin/python -m pytest tests/ -v`
(single test: `pytest tests/test_foo.py::test_name -v`).

### 3.6 `motif_analysis/` — one-time genome/promoter prep (HPC)

These are the source of truth for preparing the reference FASTA once:

| Script | Job |
|---|---|
| `TSS_extract.sh` | Extract TSS positions from a Gencode GTF → `mm10_TSS.bed`. |
| `promoter_extract.sh` | `bedtools slop` around each TSS → `mm10_promoter_<bp>bp.bed`. |
| `promoter_seq_extract.slurm` | `bedtools getfasta -s` → gene-named promoter FASTA (headers `GENE::coords`). |
| `batch_AME.sh` | The original per-module AME batching (awk-split `HM_ModGene.csv`, `seqkit grep`, skip <50 genes, `sbatch` one `AME.slurm` per module). |
| `AME.slurm` | One AME job (JASPAR/HOCOMOCO/CIS-BP paths hard-coded to a cluster). |
| `AME_result_wrapup.py` / `_lite.py` | The original result collectors that `clusmap/motif.py` was refactored from. |

`motif_run.py` re-implements the FASTA-splitting in pure Python (no seqkit) and
generalizes the run modes, but the prep scripts here remain the reference for
generating `promoter_fasta`.

### 3.7 Packaging & infra

- **`pyproject.toml`** — setuptools build, extras `sc`/`app`/`agent`/`excel`/`all`
  /`dev`, console scripts `clusmap-config` (= `config:_cli`) and `clusmap-agent`
  (= `agent.runner:main`). New submodules must be added to `[tool.setuptools]
  packages`.
- **`Dockerfile`** — python:3.11-slim + MEME-suite (compiled) + bedtools +
  samtools + seqkit + `pip install .[all]`. Lets the motif step run without a
  cluster.
- **`MANIFEST.in`** — ships `LICENSE`, `README.md`, `USAGE.md`,
  `config_example.yaml`, `PanglaoDB_marker.tsv`, `motif_analysis/`, `py.typed`.
  (The heavy genome/motif reference data is intentionally **never** packaged.)
- **`.gitignore`** — ignores `*_out/` (all run outputs, incl. `demo_out/` and
  `clusmap_out/`), `.claudeignore`, `CLAUDE.md`, `tests/`. **Caveat**: this means
  the demo data in `demo_out/` is *not* committed.

### 3.8 Quick lookup: "where do I change X?"

| I want to… | File → function |
|---|---|
| Add a new file format | [clusmap/io.py](clusmap/io.py) → `import_data` (branch on `ext`) |
| Tune auto-detection of messy tables | [clusmap/io.py](clusmap/io.py) → `_find_dense_block`, `_infer_orientation` |
| Change clustering | [clusmap/cluster.py](clusmap/cluster.py) → `gen_mod` |
| Add a module edit operation | [clusmap/state.py](clusmap/state.py) → `ModuleState` (then wire into gui/agent) |
| Change heatmap looks | [clusmap/plot.py](clusmap/plot.py) → `bulk_hm` (+ `util.normalize_rows`) |
| Change row normalization | [clusmap/util.py](clusmap/util.py) → `normalize_rows` |
| Add a cell-type annotation method | [clusmap/annotate.py](clusmap/annotate.py) |
| Add a GO library/organism rule | [clusmap/annotate.py](clusmap/annotate.py) → `mod_GO` |
| Add a module-level analysis | [clusmap/analysis.py](clusmap/analysis.py) (mirror existing functions) |
| Change the AME command / run modes | [clusmap/motif_run.py](clusmap/motif_run.py) → `_ame_cmd`, `run_ame` |
| Change motif result filtering | [clusmap/motif.py](clusmap/motif.py) → `module_motif` |
| Add a config key / CLI flag | [clusmap/config.py](clusmap/config.py) → `DEFAULTS`, `_cli` |
| Add a batch-pipeline step | [clusmap/pipeline.py](clusmap/pipeline.py) + `config_example.yaml` |
| Add a GUI block | [clusmap/gui_steps.py](clusmap/gui_steps.py) → `STEPS` (+ pure `_run_*`) |
| Change the Plotly figure | [clusmap/interactive.py](clusmap/interactive.py) → `build_figure` |
| Add an agent tool | [clusmap/agent/tools.py](clusmap/agent/tools.py) → `TOOLS` + `_t_*` |
| Add an LLM provider | [clusmap/agent/backends.py](clusmap/agent/backends.py) (implement `LLMBackend`) |
| Add/extend tests | `tests/` (see `test_table_detect.py` for the pattern) |
| Change packaged files | `pyproject.toml` + `MANIFEST.in` |

### 3.9 Concepts & conventions you must know

1. **Two label numberings** — `raw_labels` (cutree) vs `hm_labels` (1-based,
   dendrogram order). Always confirm which a table/file uses before joining.
2. **Everything is derived from `ModuleState`** — don't store module identity in
   a parallel structure; mutate the state and re-run the step.
3. **`*_args` dicts are pass-through `**kwargs`** to the underlying library —
   unvalidated by design.
4. **`hm` is a seaborn `ClusterGrid`**, not a matplotlib Figure. Its `.fig` is
   the figure; `.ax_heatmap` the heatmap axis; `.dendrogram_row.reordered_ind`
   the row order. The module-colour band axis is attached as `hm.ax_band`.
   (The GUI checks `hasattr(obj, 'fig') and hasattr(obj, 'ax_heatmap')` to tell
   a real clustermap from a stray widget value — see `gui_steps._is_hm`.)
5. **Config resolution** is explicit arg → `CLUSMAP_<KEY>` env →
   `~/.clusmap/config.yaml` → default; heavy reference paths degrade gracefully.
6. **Normalization is done by clusmap, once** (`util.normalize_rows`); the
   seaborn `z_score`/`standard_scale` flags are stripped out of `hm_args`.
7. **Order of saved files**: `gen_mod` writes `ModGene.csv` (raw ids) then
   `HM_ModGene.csv` (heatmap ids) via `state.save()`; `bulk_hm` rewrites
   `HM_ModGene.csv` **with colours**. If you see two different `HM_ModGene.csv`
   timestamps, remember `bulk_hm` overwrites it.
8. **`state.split()` needs expression data** — build the state through `gen_mod`
   (which attaches `data=rna_df`), not by constructing `ModuleState` by hand
   without `data`.
9. **GO needs network** (Enrichr); motif needs reference files on disk; both
   degrade to clear warnings/errors, not crashes.

### 3.10 Known quirks / gaps (worth being aware of when building on this)

- **No `.streamlit/config.toml` in the repo**, although `USAGE.md` says the app
  ships one (it should raise the upload limit to 4 GB). Create
  `.streamlit/config.toml` with `[server] maxUploadSize = 4000` if you want the
  bigger in-browser upload limit.
- **`demo_out/` is git-ignored** (`*_out/`), so a fresh clone has no demo data or
  example outputs. Decide whether to commit a `data/` folder + reference outputs
  (and carve it out of `.gitignore`) so new contributors have runnable examples.
- **Two motif wrap-up implementations** now coexist: the original
  `motif_analysis/AME_result_wrapup*.py` (cluster-era, has a `p-value`/`adj_p-value`
  inconsistency) and the refactored `clusmap/motif.py::module_motif`. Prefer the
  package version; consider deleting the `_wrapup` scripts to avoid drift.
- **`gen_mod`'s function default `outdir` is `"."`** while `config_example.yaml`
  uses `"clusmap_out"`; batch runs write to `clusmap_out` but a bare
  `cm.gen_mod(rna)` writes into the current directory.
- `module_preservation` and `project_modules` require aligned gene/sample axes
  and raise/warn when they can't compute; the GUI and agent surface those as
  clear messages.
- Version is `0.2.0`; there is **no CI and no coverage config** — only the two
  test files under `tests/`. The `dev` extra pins `pytest`, `build`, `twine`.
