# clusmap — Usage Guide (v0.2)

A toolkit for bulk RNA-seq **module discovery → clusterheatmap → annotation**.
The whole package is driven by one object: **`ModuleState`**, returned by
`gen_mod`. Every plot and every annotation step reads from it, so module
identity stays consistent across the pipeline.

> Run with the project env: `/opt/anaconda3/envs/clusmap/bin/python`

---

## 0. Install / import

```python
import clusmap as cm        # run from the folder that contains clusmap/
```

The public API: `import_data`, `preprocess`, `extract_color_cat`, `gen_mod`,
`bulk_hm`, `cluster_sample_stats`, `module_color_map`, `compute_pseudo_bulk`,
`pseudo_bulk_hm`, `sc_marker_hm`, `celltype_selection`, `mod_GO`, `module_motif`,
`run_pipeline_from_config`, `ModuleState`, plus the spatial API (§14):
`import_spatial`, `from_adata`, `run_leiden`, `select_hvgs`,
`spatial_module_scores`, `assign_spots_to_modules`, `add_module_expression`,
`sort_columns_by_leiden`, `spatial_hm`, `plot_spatial_modules`,
`plot_spatial_expression`.

---

## 1. Load data — only the path is required

Format is auto-detected (`.h5ad`, `.tsv`, `.csv`, `.txt`, `.xlsx`, `.pkl`),
delimiter is sniffed, and leading annotation columns are detected automatically.

```python
rna = cm.import_data("counts.tsv")           # genes x samples DataFrame
rna = cm.import_data("atlas.h5ad")           # h5ad -> transposed to genes x samples
rna = cm.import_data("data.gene", header_path="header.gene", str_col_num=2, index_col=2)
```

Only reach for `str_col_num` / `index_col` / `header_path` when auto-detection
guesses wrong (e.g. unusual multi-column annotation files).

## 2. Preprocess

```python
rna = cm.preprocess(rna)                      # min_expr=10, log2, drop zero-variance
rna = cm.preprocess(rna, min_expr=0, log_base=None)   # skip filtering / log
```

## 3. Cluster into modules → `ModuleState`

```python
state = cm.gen_mod(rna, outdir="out")         # only rna required
state = cm.gen_mod(rna, deepSplit=2, minClusterSize=20)   # the two knobs you'll touch
```

`deepSplit` 0 (few big modules) → 4 (many small). `gen_mod` saves
`modules.pkl`, `linkage.pkl`, `ModGene.csv`, `HM_ModGene.csv`,
`module_state.pkl` to `outdir`.

Inspect the state:

```python
state.n_modules                 # number of modules
state.gene_module("Sox2")       # which module a gene is in (case-insensitive)
state.module_genes(3)           # genes in module 3 (heatmap order)
```

## 4. Heatmap — only `rna` + `state` required

```python
hm = cm.bulk_hm(rna, state, outdir="out")
```

Common options (all optional, sensible defaults):

```python
hm = cm.bulk_hm(
    rna, state,
    title="My experiment",                       # custom title
    norm_method="z_score",                       # z_score|centralize|min_max|max|none
    goi_list=["Sox2", "Pax6"], goi_size=10,      # genes of interest on the right
    hm_args={"figsize": (14, 10), "col_cluster": False,
             "dendrogram_ratio": (0.15, 0.1), "xticklabels": True},
    mod_palette="hsv", mod_num_font={"fontsize": 6, "fontweight": "bold"},
    row_band=True, row_band_font={"band_width": 0.03, "alpha": 0.35},
    save_format="pdf",                           # png|pdf|svg
)
```

### Column color matrices (annotation bands on top)

```python
col_cat = cm.extract_color_cat(rna, {"organ": ["thymus", "limb"]})   # auto from sample names
# or supply your own per-column lists / Series aligned to rna.columns
hm = cm.bulk_hm(
    rna, state,
    col_cat=col_cat,
    col_color_manual={"organ": {"thymus": "#e41a1c", "limb": "#377eb8"}},  # optional
    col_legend=["organ"],
)
```

Writes `heatmap.<fmt>` and `HM_ModGene.csv` (gene → heatmap module id + color).

## 5. Get more information about modules

**Without single cell — PanglaoDB marker swarmplot** (most-correlated cell types):

```python
selected = cm.celltype_selection("PanglaoDB_marker.tsv", "out/HM_ModGene.csv",
                                 max_celltypes_per_module=5)
cm.sc_marker_hm(rna, hm, "PanglaoDB_marker.tsv",
                celltype_col="cell type", gene_col="official gene symbol",
                celltype=selected, outdir="out")
```

**With single cell — pseudo-bulk heatmap beside the main one:**

```python
pb = cm.compute_pseudo_bulk("sc.h5ad", celltype_key="celltype",
                            chunk_size=20000)        # chunk_size only for big files
cm.pseudo_bulk_hm(rna, hm, pseudo_bulk_df=pb, log_base=2, outdir="out")
```

## 6. GO enrichment (per module or selected modules)

```python
cm.mod_GO("out/HM_ModGene.csv", organism="Mouse")          # module_col defaults to hm_mod
cm.mod_GO("out/HM_ModGene.csv", organism="Human",
          mod=[1, 3, 5], GO_category=["BP", "MF", "CC"], bg=rna)   # rna -> background genes
```

## 7. Motif analysis (genome on HPC/Docker, driven from Python)

The genome FASTA and motif DBs (JASPAR/HOCOMOCO/CIS-BP) are large and live on
your machine or cluster — clusmap stores only **paths**, via `clusmap-config`:

```bash
clusmap-config set promoter_fasta /data/mm10_promoter_500bp.fa   # gene-named headers
clusmap-config set motif_db.JASPAR /db/JASPAR2026_vertebrates.meme
clusmap-config set run_mode slurm        # local | slurm | ssh
clusmap-config set ssh_host me@cluster   # for run_mode=ssh
clusmap-config show
```

One call splits modules → per-module promoter FASTAs → AME jobs → result tables:

```python
cm.motif_pipeline("out/HM_ModGene.csv", outdir="motif_out")   # uses config defaults
cm.motif_pipeline("out/HM_ModGene.csv", run_mode="slurm", dry_run=True)  # preview jobs
```

- `run_mode="local"` needs `ame` on PATH (use the Docker image below).
- `run_mode="slurm"` writes + submits one `sbatch` script per module.
- `run_mode="ssh"` scp's FASTAs to `remote_workdir` and submits there.

For SLURM/SSH the jobs run asynchronously; once they finish, collect results:

```python
cm.module_motif("motif_out/AME_out", "motif_out/motif_results", mode="e", e_thresh=1e-2)
# modes: top_n (top_n=10) | p (p_thresh) | e (e_thresh) | pe
```

Genome/promoter prep (run once) is in [`motif_analysis/`](motif_analysis):
`TSS_extract.sh` → `promoter_extract.sh` → `promoter_seq_extract.slurm`.

**No cluster?** Build the Docker image (bundles MEME-suite + bedtools + seqkit):

```bash
docker build -t clusmap .
docker run --rm -v $PWD:/work -w /work -v /data:/data clusmap \
    python -c "import clusmap as cm; cm.motif_pipeline('HM_ModGene.csv', run_mode='local')"
```

## 8. Per-cluster × sample statistics

```python
cm.cluster_sample_stats(rna, "out/HM_ModGene.csv", hm=hm,
                        stats=["mean", "std", "cv"])
```

## 9. Editing modules (manual curation)

`ModuleState` edits propagate to every downstream output — just re-run the step
after editing and re-save.

```python
state.merge(2, 3)                 # merge module 2 and 3
state.split(1)                    # split module 1 into two (re-clusters its genes)
state.reassign(["Gene1", "Gene2"], 4)   # move genes into module 4 (0 = unassign)
state.save("out")                 # rewrites HM_ModGene.csv + module_state.pkl
cm.ModuleState.load("out")        # reload later
```

## 10. Batch mode (YAML, unchanged workflow)

```python
ctx = cm.run_pipeline_from_config("config_example.yaml")
# ctx holds rna_df, state, hm, eigengenes, hub_genes, trait_cor, preservation,
# projection, report, motif ... whichever sections you included, for further work
```

[`config_example.yaml`](config_example.yaml) at the repo root is a **complete
reference** — every tunable parameter of every pipeline function (including the
`*_args` pass-through dicts like `hm_args`, `swarm_args`, `enrichr_args`, and
the `dynamicTreeCut` extras under `gen_mod`) is listed at its actual default,
with an explanatory comment. A step only runs if its top-level key is present
in the file — delete or comment out sections you don't need rather than
editing every value. Covers the full pipeline: import → preprocess → cluster →
heatmap → stats → {pseudo-bulk | marker} annotation → GO → motif → eigengenes →
hub genes → trait correlation → preservation → projection → HTML report.

## 11. No-code app (Streamlit) — run the whole pipeline in the browser

```bash
/opt/anaconda3/envs/clusmap/bin/streamlit run clusmap/app.py
```

One page, **one block per pipeline function** — no sidebar, no tabs. Each block
shows **every tunable parameter** with its default as a **gray placeholder**
(leave a field blank to use that default; the commonly-tuned params are inline,
the rarely-used ones fold under **⚙ Advanced parameters**). True/False params are
click-to-toggle switches; **hover a parameter name or block title for its
explanation**. Each block has a **▶ Run** button and shows its result (table /
figure) in place; status **messages pop up in a closable window** instead of
piling up. Output files go to the "Output folder" set at the top.

**Run the whole thing in one click.** The top panel lets you:

- pick an **Annotation mode** — *Markers (PanglaoDB)* (no single-cell data) or
  *Single-cell pseudo-bulk*; for pseudo-bulk a **`.h5ad` path box appears right
  there** so you can point at the single-cell file without scrolling down;
- **toggle which downstream analyses** (GO, motif, hub genes, module–trait,
  preservation, projection, report, per-cluster stats) are included;
- press **⚡ Run pipeline** to run import → preprocess → cluster → heatmap → your
  chosen annotation → enabled downstream steps in order.

Core steps 1–4 always run on one click; everything else is governed by those
toggles. You can also run any block on its own with its **▶ Run** button.

Blocks: Import → Preprocess → Generate modules → Clusterheatmap (+ **module
curation**, see below) → Per-cluster stats → {Cell-type markers / Single-cell
pseudo-bulk} → GO → Motif (AME) → Hub genes → Module–trait → Preservation →
Project modules → HTML report. Every parameter here matches
[`config_example.yaml`](config_example.yaml) one-to-one (including the `*_args`
dicts — fill them as JSON, e.g. `{"vmin": -2, "vmax": 2}`).

### Curate modules — folded into the Clusterheatmap block

After the heatmap renders, an **✏️ Curate modules** expander appears under it.
Curation is **text-based and reliable** (the old drag-to-select Plotly figure was
fiddly and is gone):

- **Merge** two neighbouring modules (pick two from the list),
- **Split** a module into two sub-branches,
- **Reassign** named genes into a module (or `0` to unassign),
- **Find** which module a gene is in.

Each edit mutates the same `ModuleState` and **immediately re-renders the heatmap
above**. Press **💾 Save** to write `module_state.pkl`, `module_state.json`,
`HM_ModGene.csv` and a refreshed `heatmap.png` to the output folder — those feed
straight into the downstream blocks and back into your notebook
(`cm.ModuleState.load("out")`).

> **Motif block:** choose an **organism** — human/mouse references are bundled
> into the software at publish time, so you just select the organism; pick
> **other** to supply your own promoter FASTA + `.meme` motif DB. Local runs need
> MEME-suite `ame` on PATH (the app warns if it's missing; use the Docker image
> or `run_mode = slurm/ssh`). Keep `dry_run` on to preview the jobs first.

### Loading data: in the app, or handed off from a notebook

Load in the **Import block** — type a path or **upload** the file (uploads work
for any block that takes a file). For **non-standard matrices** (custom
header/index/annotation columns) either set the Import block's Advanced options,
or build `rna`/`state` once in a notebook and hand them straight to the app —
which reuses the live objects and never re-parses the file:

```python
import clusmap as cm
rna   = cm.import_data("weird_layout.tsv", str_col_num=2, index_col=2,
                       header_path="header.txt")   # whatever your file needs
rna   = cm.preprocess(rna)
state = cm.gen_mod(rna)
cm.launch_app(rna, state)            # opens the app, already populated
# cm.launch_app(rna, state, outdir="results", port=8502, headless=True)  # remote/HPC
```

It returns the Streamlit subprocess (`.terminate()` to stop). On a remote
machine, pass `headless=True` and tunnel the port
(`ssh -L 8501:localhost:8501 ...`). Then just click **▶ Run** on the
Clusterheatmap block to draw and curate.

### Upload size — is there a limit?

Yes. The in-browser uploader copies the whole file through the server, with a
default cap of **200 MB**. The app ships a [`.streamlit/config.toml`](.streamlit/config.toml)
that raises it to **4 GB** (`server.maxUploadSize = 4000`); change that number, or
override per run:

```bash
/opt/anaconda3/envs/clusmap/bin/streamlit run clusmap/app.py --server.maxUploadSize 8000
```

For **very large files** (big `.h5ad` atlases), don't upload — use the **path
box** in the block instead. It reads straight from disk with **no size limit**
and no copy. The notebook hand-off (`cm.launch_app`) and the plain Python API
have no limit either. (Tip: the pseudo-bulk step also has a `chunk_size` option
for memory-efficient backed reads of huge atlases.)

> Programmatic equivalent of the curation edits (no app needed):
> ```python
> state.split(1); state.merge(2, 3); state.reassign(["Sox2"], 4)
> state.save("out")     # rewrites HM_ModGene.csv for downstream steps
> ```

---

## 12. Conversational agent (no YAML, no notebook)

Talk to the pipeline instead of editing parameters. The agent asks for what it
needs ("about how many modules?", "what's the input file?") and calls the
clusmap functions for you. **Provider-agnostic** — pick the LLM by supplying a
provider + API key; Claude works out of the box. Runs on a laptop, **no GPU**
(the model runs on the provider's servers).

```bash
export ANTHROPIC_API_KEY=sk-ant-...
/opt/anaconda3/envs/clusmap/bin/python -m clusmap.agent --outdir results
# other providers:  python -m clusmap.agent --provider <name> --api-key <key>
```

Then just type, e.g.:

> Cluster `mouse_thymus.tsv` into about 8 modules, draw the heatmap coloured by
> organ, then run GO on the modules.

Programmatic entry point:

```python
from clusmap.agent import chat
chat(provider="anthropic", outdir="results")
```

The agent exposes these tools: `load_data`, `preprocess`, `cluster`,
`make_heatmap`, `module_overview`, `search_gene`, `list_module_genes`,
`edit_modules` (split/merge/reassign), `celltype_swarm`, `go_enrichment`,
`save_state`, `set_output_dir`. To add another LLM provider, implement one
adapter class — see `OpenAIBackend` in `clusmap/agent/backends.py`. Requires
`pip install anthropic` for the default Claude backend.

## 13. Module-level analyses (eigengenes, hubs, traits, report)

WGCNA-style analyses, all built on the `ModuleState`:

```python
eig  = cm.module_eigengenes(rna, state, outdir="out")     # modules x samples (first PC)
hubs = cm.hub_genes(rna, state, eigengenes=eig, top_n=20) # kME ranking; hub genes per module
cm.top_hubs(hubs, top_n=10)                               # {module: [top hub genes]}

# correlate modules with sample traits (numeric or categorical -> one-hot), FDR-corrected
traits = pd.DataFrame({"organ": [...], "age": [...]}, index=rna.columns)
cor, fdr = cm.module_trait_correlation(eig, traits, outdir="out")   # + heatmap

# assign a NEW dataset's genes to these modules (shared sample columns required)
cm.project_modules(new_rna, eig, min_cor=0.5, outdir="out")

# are the modules real? test reproducibility in an independent dataset
pres = cm.module_preservation(rna, state, test_rna, n_perm=200, outdir="out")
# -> Zsummary per module: >=10 strong, >=2 moderate, <2 not preserved (+ barplot)

# one self-contained HTML bundling size + eigengene + hubs + GO/celltype/motif (if present)
cm.module_report("out", rna, state)
```

Outputs: `module_eigengenes.tsv`, `hub_genes.tsv`, `module_trait_cor.tsv` +
`module_trait_fdr.tsv` + `module_trait_heatmap.png`, `projected_modules.tsv`,
`module_report.html`. The agent (§12) exposes `hub_genes` and `module_report` too.

## 14. Spatial transcriptomics (10x Visium)

A spot is a pseudo-bulk sample, so spatial data runs through the exact same
engine as bulk — `gene × sample` becomes `gene × spot`. Two kinds of module
coexist:

- **spot modules** — the *Leiden* clusters of `import_spatial`'s scanpy
  pipeline (a module **formed by spots**), used for the spatial cluster map;
- **gene modules** — the usual `gen_mod` dynamic-tree-cut modules (a module
  **formed by genes**), used for the heatmap and per-module expression plots.

Everything is scanpy-native (`sc.pl.spatial` over the H&E image). Requires
`pip install "clusmap[spatial]"` (adds `h5py`, `scanpy`, `leidenalg`).

### Import (folder or `.h5ad`)

`import_spatial` accepts **either** a 10x Visium output folder **or** a pre-made
`.h5ad` file, and — if the data carries no preprocessing — runs the standard
scanpy pipeline (normalize → log1p → HVG → scale → PCA → neighbors → Leiden) so
the returned object always has `adata.obs['leiden']` (spot clusters):

```python
sdata = cm.import_spatial("V1_Mouse_Brain_Sagittal_Posterior")  # folder
sdata = cm.import_spatial("brain.h5ad")                         # .h5ad file
sdata.adata                # scanpy AnnData: obs['leiden'], obsm['spatial'],
                           #   uns['spatial'] (H&E images + scale factors)
sdata.rna                  # genes x spots log1p-HVG matrix (for the clusmap engine)
sdata.leiden               # spot -> Leiden cluster (Series)
sdata.coords               # barcode -> x, y (pixels)
```

Already loading via scanpy? Hand the AnnData straight in:

```python
import scanpy as sc
adata = sc.datasets.visium_sge("V1_Mouse_Brain_Sagittal_Posterior")
sdata = cm.from_adata(adata, n_top_genes=2000)
```

`n_top_genes` (default 2000) controls the highly-variable-gene subset;
`resolution` (default 0.5) the Leiden resolution. A preprocessed `.h5ad` that
already has `leiden` is used as-is (no re-preprocessing).

### Gene modules (module formed by genes)

Cluster exactly as for bulk — `gen_mod` already receives the log1p HVG matrix
from `import_spatial`:

```python
state = cm.gen_mod(sdata.rna, deepSplit=1, minClusterSize=30, outdir="spatial_out")
```

### Two-version clusterheatmap (reuses `bulk_hm`)

The heatmap output is kept for spatial data, in **two versions**, with the
**Leiden** annotation as the default column colour band:

- **v1** — rows *and* columns both clustered (`col_cluster=True`);
- **v2** — rows clustered, columns grouped by Leiden cluster and clustered
  *within* each group (`col_cluster=False` with pre-sorted columns).

```python
hm_v1, hm_v2 = cm.spatial_hm(sdata.adata, sdata.rna, state, outdir="spatial_out")
# -> spatial_heatmap_v1.png, spatial_heatmap_v2.png
```

Add more column colour bars exactly like `bulk_hm`, e.g. a per-spot metadata
Series: `cm.spatial_hm(..., col_cat={'leiden': ..., 'region': region_series})`.

### Marker-celltype swarmplot (single-cell annotation)

The bulk single-cell annotation steps work unchanged on the gene modules:

```python
selected = cm.celltype_selection(marker_path, "spatial_out/HM_ModGene.csv",
                                 outdir="spatial_out")
cm.sc_marker_hm(sdata.rna, hm_v1, marker_path, celltype=selected,
                celltype_col="cell type", gene_col="official gene symbol",
                outdir="spatial_out")          # -> heatmap_swarm.png
```

### Spatial plots (sc.pl.spatial over the H&E image)

`spatial_module_scores` scores each spot × **gene module** (default
`method="mean"` = mean of the module's genes, or `method="eigengene"`), z-scored
per module so modules are comparable; `add_module_expression` drops those
columns into `adata.obs` as `module_<id>_expr`:

```python
scores = cm.add_module_expression(sdata.adata, sdata.rna, state)  # -> spots x modules
```

**`spatial_modules.png`** — the **Leiden** spot clusters (module formed by
spots):

```python
cm.plot_spatial_modules(sdata.adata, img_key="hires", outdir="spatial_out")
```

**`spatial_module_expression.png`** — one subplot per **gene module**, its mean
per-spot expression across the tissue:

```python
cm.plot_spatial_expression(sdata.adata, img_key="hires", outdir="spatial_out")
```

Both are plain `sc.pl.spatial` calls, so you can colour the tissue yourself in a
notebook:

```python
sc.pl.spatial(sdata.adata, color="leiden", img_key="hires")
sc.pl.spatial(sdata.adata, color=["module_1_expr", "module_2_expr"], img_key="hires", ncols=4)
```

### Inspect one module interactively

To zoom into a single gene module, export its genes and open them in
**cellxgene** (https://cellxgene.cziscience.com) against the same Visium
`.h5ad`:

- get the module's gene list from `spatial_out/HM_ModGene.csv` (rows with
  `hm_mod == N`), e.g. `awk -F, '$1==3 {print $2}' spatial_out/HM_ModGene.csv`;
- in cellxgene, add those genes to the colour-by panel — the spatial plot will
  highlight where that module's genes are expressed.

Full end-to-end example: run `python demo_spatial.py [DATA_DIR_OR_H5AD]`.

## 15. Packaging & publishing

clusmap is a proper installable package ([pyproject.toml](pyproject.toml)):

```bash
pip install -e ".[all]"          # local dev (core + sc + app + agent + excel)
python -m build                  # build sdist + wheel into dist/
twine upload dist/*              # publish to PyPI
```

Console scripts installed: `clusmap-config` (path/HPC config) and
`clusmap-agent` (chat agent). Optional extras keep the core lean:
`sc` (single-cell), `app` (Streamlit/Plotly), `agent` (anthropic), `excel`,
`spatial` (h5py + scanpy + leidenalg).

To publish the tool:
1. **PyPI** — `python -m build && twine upload dist/*`.
2. **bioconda** — add a recipe so users get `conda install -c bioconda clusmap`
   (pulls MEME-suite/bedtools as conda deps for the motif step).
3. **Docker** — `docker build -t clusmap .` then push to a registry; bundles the
   MEME toolchain so the motif step runs without a cluster.
4. **Zenodo** — link your GitHub repo and tag a release to mint a citable DOI.
5. **Docs** — host USAGE.md on Read the Docs / GitHub Pages.

The heavy reference data (genome, motif DBs) is **never** packaged — it stays on
the user's disk/cluster and is referenced through `clusmap-config`.

---

### Outputs cheat-sheet
`ModGene.csv` (original module ids) · `HM_ModGene.csv` (heatmap ids + colors) ·
`heatmap.*` · `pb_heatmap.*` · `heatmap_swarm.*` · `celltype_selection/` ·
`<GO category>/module_*_GO.csv` · `module_motif_*.tsv` · `module_state.pkl` ·
`module_eigengenes.tsv` · `hub_genes.tsv` · `module_trait_*.tsv` ·
`projected_modules.tsv` · `module_preservation.tsv` · `module_report.html` ·
`spatial_module_expression.png` (per-module expression grid) ·
`spatial_modules.png` (spot → module map) · `spatial_module_scores.tsv`.

### Which clustering files do I get, and why? (no redundancy)
`gen_mod` writes four things and **`module_state.pkl` is the one that matters**:

| File | What | Keep? |
|---|---|---|
| `module_state.pkl` | The full `ModuleState` — **contains the linkage matrix + labels + gene order + edit history**. Reload with `cm.ModuleState.load(outdir)`. | ✅ primary |
| `module_state.json` | Human-readable summary (n_modules, metric, method, edit history) — read it without unpickling. Content is derivable from the pkl; it's a tiny convenience. | ✅ (tiny) |
| `ModGene.csv` | gene → **original** cutree module id. | ✅ |
| `HM_ModGene.csv` | gene → **heatmap** module id (+ colours). Different numbering; used by GO/motif/celltype. | ✅ |

`linkage.pkl` and `modules.pkl` are **no longer written by default** — they were
redundant (the linkage lives in `module_state.pkl` as `state.linkage`; the raw
labels as `state.raw_labels`). Pass `cm.gen_mod(..., save_raw=True)` only if you
want the raw `cutreeHybrid` diagnostics dict or a standalone linkage file.
