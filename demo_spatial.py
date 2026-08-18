#!/usr/bin/env python
"""End-to-end spatial-transcriptomics demo for clusmap (scanpy-centric).

Treats a 10x Visium dataset as bulk data (gene x spot): the scanpy pipeline
finds **spot clusters** (Leiden), the clusmap engine finds **gene modules**
(hierarchical clustering + dynamic tree cut), and both are rendered over the
H&E image with ``sc.pl.spatial``.

Pipeline
--------
1. ``import_spatial`` — accepts a **Visium folder** OR a pre-made **.h5ad**
   file; runs normalize -> log1p -> HVG -> scale -> PCA -> neighbors -> Leiden
   when the data is raw, so ``adata.obs['leiden']`` always holds spot clusters.
2. ``gen_mod`` — gene modules (module formed *by genes*) on the log1p HVG
   gene x spot matrix.
3. Bulk analyses — two-version clusterheatmap (``spatial_hm``, reusing
   ``bulk_hm``) + optional marker swarmplot (``celltype_selection`` +
   ``sc_marker_hm``).
4. Spatial plots — ``sc.pl.spatial(adata, color=..., img_key='hires')``:
   ``spatial_modules.png`` (Leiden spot clusters) and
   ``spatial_module_expression.png`` (per-gene-module mean expression).

Usage
-----
    /opt/anaconda3/envs/clusmap/bin/python demo_spatial.py [DATA] [options]

``DATA`` is a Visium output folder (default
``spatial/data/V1_Mouse_Brain_Sagittal_Posterior``) or a path to a .h5ad file.
Download the default dataset with
``sc.datasets.visium_sge('V1_Mouse_Brain_Sagittal_Posterior')``.

Outputs (in ``spatial_out/`` by default):
    spatial_modules.png             Leiden spot clusters over the H&E image
    spatial_module_expression.png   per-gene-module mean expression (H&E)
    spatial_heatmap_v1.png          gene x spot clustermap, rows+cols clustered
    spatial_heatmap_v2.png          rows clustered, cols grouped by Leiden
    heatmap_swarm.png               (if --marker-file given) marker swarm plot
    HM_ModGene.csv / module_state.pkl / module_eigengenes.tsv  usual artifacts
"""
from __future__ import annotations

import argparse
import os
import sys

# Force a non-interactive backend so the script also runs headless (SSH/CI).
import matplotlib
matplotlib.use("Agg")

import clusmap as cm


DEFAULT_DATA = os.path.join("spatial", "data", "V1_Mouse_Brain_Sagittal_Posterior")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", nargs="?", default=DEFAULT_DATA,
                    help="Visium output folder or a .h5ad file "
                         "(default: %(default)s)")
    ap.add_argument("--outdir", default="spatial_out", help="where to write outputs")
    ap.add_argument("--n-hvgs", type=int, default=2000,
                    help="top highly-variable genes for clustering (default: %(default)s)")
    ap.add_argument("--resolution", type=float, default=0.5,
                    help="Leiden resolution for spot clusters (default: %(default)s)")
    ap.add_argument("--marker-file", default=None,
                    help="PanglaoDB-style marker table (tab-separated, 'cell type' + "
                         "'official gene symbol' columns) for the swarmplot annotation "
                         "(default: skip the swarmplot)")
    ap.add_argument("--img-key", default="hires", choices=["hires", "lowres"],
                    help="H&E backdrop resolution for sc.pl.spatial (default: %(default)s)")
    ap.add_argument("--no-spatial", action="store_true",
                    help="skip the sc.pl.spatial figures (heatmap only)")
    args = ap.parse_args(argv)

    if not os.path.exists(args.path):
        sys.exit(f"Data not found: {args.path}\n"
                 f"Download it first, e.g. "
                 f"sc.datasets.visium_sge('V1_Mouse_Brain_Sagittal_Posterior')")

    print("=" * 70)
    print("1. Import spatial data (folder or .h5ad) + scanpy preprocess/Leiden")
    print("=" * 70)
    sdata = cm.import_spatial(args.path, n_top_genes=args.n_hvgs,
                              resolution=args.resolution)
    print(f"   {sdata.rna.shape[0]} highly-variable genes x "
          f"{sdata.rna.shape[1]} spots")
    print(f"   {sdata.leiden.nunique()} Leiden spot clusters")
    print(f"   coords: {sdata.coords.shape[0]} spots")

    print("\n" + "=" * 70)
    print(f"2. Cluster genes into modules (module formed by genes)")
    print("=" * 70)
    state = cm.gen_mod(sdata.rna, deepSplit=1, minClusterSize=30,
                       outdir=args.outdir)
    print(f"   {state.n_modules} gene modules")

    print("\n" + "=" * 70)
    print("3. Two-version gene x spot clustermap (reuses bulk_hm)")
    print("=" * 70)
    print("   v1: rows and columns both clustered")
    print("   v2: rows clustered, columns grouped by Leiden + clustered within group")
    hm_v1, hm_v2 = cm.spatial_hm(sdata.adata, sdata.rna, state,
                                 outdir=args.outdir)

    if args.marker_file:
        print("\n" + "=" * 70)
        print("4. Marker-celltype annotation (swarmplot)")
        print("=" * 70)
        module_path = os.path.join(args.outdir, "HM_ModGene.csv")
        selected = cm.celltype_selection(args.marker_file, module_path,
                                         outdir=args.outdir)
        if selected:
            print(f"   selected {len(selected)} cell types: {', '.join(selected)}")
            cm.sc_marker_hm(sdata.rna, hm_v1, args.marker_file,
                            celltype=selected, celltype_col="cell type",
                            gene_col="official gene symbol", outdir=args.outdir)
        else:
            print("   no marker overlap — skipping the swarmplot")

    if not args.no_spatial:
        print("\n" + "=" * 70)
        print("5. Spatial plots via sc.pl.spatial (over the H&E image)")
        print("=" * 70)
        cm.add_module_expression(sdata.adata, sdata.rna, state)
        cm.plot_spatial_modules(sdata.adata, img_key=args.img_key,
                                outdir=args.outdir)
        cm.plot_spatial_expression(sdata.adata, img_key=args.img_key,
                                   outdir=args.outdir)

        print("\n" + "=" * 70)
        print("6. Module sizes")
        print("=" * 70)
        for m in sorted({int(x) for x in state.hm_labels if x != 0}):
            genes = state.module_genes(m)
            print(f"   module {m:>3}: {len(genes):>5} genes, "
                  f"top genes: {', '.join(genes[:5])}")

        print("\n" + "=" * 70)
        print("7. Interactive exploration")
        print("=" * 70)
        print("   - Load the saved AnnData in a notebook and colour it yourself:")
        print("       import scanpy as sc")
        print("       adata = sc.read_h5ad('<your>.h5ad')   # or keep sdata.adata")
        print("       sc.pl.spatial(adata, color='leiden', img_key='hires')")
        print("       sc.pl.spatial(adata, color=['module_1_expr', ...],")
        print("                      img_key='hires', ncols=4)")
        print("   - Map a gene module to tissue anatomy with cellxgene:")
        print("       awk -F, '$1==1 {print $2}' spatial_out/HM_ModGene.csv")

    print(f"\nDone. Outputs in {args.outdir}/")


if __name__ == "__main__":
    main()
