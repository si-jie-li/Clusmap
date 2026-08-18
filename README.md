# clusmap

**Bulk RNA-seq module discovery → clusterheatmaps → annotation**, with an
interactive editor and a conversational agent.

clusmap clusters genes (hierarchical clustering + dynamic tree cut), draws
annotated clusterheatmaps, and helps you interpret each module via PanglaoDB
marker swarmplots, single-cell pseudo-bulk, GO enrichment, and TF-motif (AME)
analysis. Every result flows through one `ModuleState` object, so manual
curation (split / merge / reassign) propagates to all downstream outputs.

## Install

```bash
pip install clusmap                 # core: clustering + heatmaps + GO + motif wrap-up
pip install "clusmap[all]"          # + single-cell, interactive app, agent, xlsx
```

Optional extras: `sc` (single-cell pseudo-bulk), `app` (Streamlit/Plotly
interactive heatmap), `agent` (Claude agent), `excel` (.xlsx input), `spatial`
(10x Visium `.h5` reading).

## Quick start

```python
import clusmap as cm
rna   = cm.import_data("counts.tsv")          # auto-detects .h5ad/.tsv/.csv/.xlsx
rna   = cm.preprocess(rna)                    # filter + log2 + drop zero-variance
state = cm.gen_mod(rna, deepSplit=1)          # -> ModuleState (modules)
hm    = cm.bulk_hm(rna, state, outdir="out")  # clusterheatmap + HM_ModGene.csv
cm.mod_GO("out/HM_ModGene.csv", organism="Mouse")
```

Three front-ends, all on the same engine:

| Want… | Use |
|---|---|
| Scripted / reproducible | the Python API, or `cm.run_pipeline_from_config("config.yaml")` |
| Point-and-click curation | `streamlit run clusmap/app.py` (click a module, split/merge, zoom, search) |
| Just describe what you want | `python -m clusmap.agent` (provider-agnostic LLM agent) |

Full walkthrough: see [USAGE.md](USAGE.md).

## Motif analysis (HPC / Docker)

The genome FASTA and motif DBs (JASPAR/HOCOMOCO/CIS-BP) are large and live on
your machine or cluster — clusmap only stores paths:

```bash
clusmap-config set promoter_fasta /data/mm10_promoter_500bp.fa
clusmap-config set motif_db.JASPAR /db/JASPAR2026_vertebrates.meme
clusmap-config set run_mode slurm        # local | slurm | ssh
```

```python
cm.motif_pipeline("out/HM_ModGene.csv", outdir="motif_out")
# -> per-module promoter FASTAs, AME jobs (local/sbatch/ssh), result tables
```

No cluster? Build the Docker image (bundles MEME-suite + bedtools + seqkit) and
run `run_mode=local`:

```bash
docker build -t clusmap .
docker run --rm -v $PWD:/work -w /work -v /data:/data clusmap \
    python -c "import clusmap as cm; cm.motif_pipeline('HM_ModGene.csv', run_mode='local')"
```

Genome/promoter prep scripts are in [`motif_analysis/`](motif_analysis).

## Spatial transcriptomics (Visium)

Treat a spot as a pseudo-bulk sample: `import_spatial` reads a 10x Visium output
directory into a gene × spot matrix, then the same clustering engine finds gene
modules, and `plot_spatial_expression` / `plot_spatial_modules` render each
module across the tissue over the H&E image (colours match the heatmap).

```python
import clusmap as cm
sdata = cm.import_spatial("V1_Mouse_Brain_Sagittal_Posterior", image="lowres")
rna   = cm.select_hvgs(sdata.rna, n_top=2000)     # 32k genes -> 2k for tractability
rna   = cm.preprocess(rna)
state = cm.gen_mod(rna, deepSplit=1, minClusterSize=30, outdir="spatial_out")

scores = cm.spatial_module_scores(rna, state, method="mean")   # spots x modules
cm.plot_spatial_expression(scores, sdata.coords, image=sdata.image,
                           scale_factors=sdata.scale_factors, outdir="spatial_out")
cm.plot_spatial_modules(sdata.coords, cm.assign_spots_to_modules(scores),
                        state=state, image=sdata.image, outdir="spatial_out")
```

Run the full demo end-to-end: `python demo_spatial.py`. Requires
`pip install "clusmap[spatial]"` (adds `h5py`).

## Publishing checklist

`python -m build` → `twine upload dist/*` (PyPI); add a bioconda recipe; tag a
GitHub release and archive on Zenodo for a citable DOI; push the Docker image to
a registry. See [USAGE.md](USAGE.md) for details.

## License

MIT — see [LICENSE](LICENSE).
