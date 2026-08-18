#!/usr/bin/env python
"""End-to-end spatial-transcriptomics demo for clusmap.

Treats a 10x Visium dataset as bulk data (gene x spot), finds gene modules with
the same clustering engine, then renders each module's expression profile across
the tissue (optionally overlaid on the H&E image).

Usage
-----
    /opt/anaconda3/envs/clusmap/bin/python demo_spatial.py [DATA_DIR]

The default DATA_DIR is ``spatial/data/V1_Mouse_Brain_Sagittal_Posterior``
(the 10x Visium mouse-brain dataset downloaded with
``sc.datasets.visium_sge("V1_Mouse_Brain_Sagittal_Posterior")``).

Outputs (in ``spatial_out/`` by default):
    spatial_module_expression.png   one subplot per module (mean expression, H&E)
    spatial_modules.png             discrete spot -> module map (heatmap colours)
    spatial_module_scores.tsv       spots x modules score table (for your own tools)
    heatmap.png / HM_ModGene.csv / module_state.pkl   the usual clusmap artifacts
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
                    help="Visium output directory (default: %(default)s)")
    ap.add_argument("--outdir", default="spatial_out", help="where to write outputs")
    ap.add_argument("--n-hvgs", type=int, default=2000,
                    help="top highly-variable genes to cluster (default: %(default)s)")
    ap.add_argument("--image", choices=["lowres", "hires", "none"], default="lowres",
                    help="H&E backdrop resolution (default: %(default)s)")
    ap.add_argument("--min-score", type=float, default=None,
                    help="min module score for a spot to be assigned (default: none)")
    ap.add_argument("--no-heatmap", action="store_true",
                    help="skip the (large) gene x spot clustermap")
    args = ap.parse_args(argv)

    if not os.path.isdir(args.path):
        sys.exit(f"Data directory not found: {args.path}\n"
                 f"Download it first, e.g. sc.datasets.visium_sge('V1_Mouse_Brain_Sagittal_Posterior')")

    print("=" * 70)
    print("1. Import spatial data (spot = pseudo-bulk)")
    print("=" * 70)
    image = None if args.image == "none" else args.image
    sdata = cm.import_spatial(args.path, image=image)
    print(f"   {sdata.rna.shape[0]} genes x {sdata.rna.shape[1]} spots")
    print(f"   coords: {sdata.coords.shape[0]} spots, image={sdata.image}")

    print("\n" + "=" * 70)
    print(f"2. Cluster genes into modules (top {args.n_hvgs} HVGs)")
    print("=" * 70)
    rna = cm.select_hvgs(sdata.rna, n_top=args.n_hvgs)
    rna = cm.preprocess(rna)
    state = cm.gen_mod(rna, deepSplit=1, minClusterSize=30, outdir=args.outdir)

    if not args.no_heatmap:
        print("\n" + "=" * 70)
        print("3. Gene x spot clustermap (same engine as bulk)")
        print("=" * 70)
        cm.bulk_hm(rna, state, outdir=args.outdir)

    print("\n" + "=" * 70)
    print("4. Per-spot module expression + spatial rendering")
    print("=" * 70)
    scores = cm.spatial_module_scores(rna, state, method="mean", norm="zscore")
    scores.to_csv(os.path.join(args.outdir, "spatial_module_scores.tsv"), sep="\t")
    print(f"   spots x modules score table -> {args.outdir}/spatial_module_scores.tsv")

    cm.plot_spatial_expression(
        scores, sdata.coords, image=sdata.image, scale_factors=sdata.scale_factors,
        outdir=args.outdir,
    )

    assign = cm.assign_spots_to_modules(scores, min_score=args.min_score)
    cm.plot_spatial_modules(
        sdata.coords, assign, state=state, image=sdata.image,
        scale_factors=sdata.scale_factors, outdir=args.outdir,
    )

    print("\n" + "=" * 70)
    print("5. Module sizes")
    print("=" * 70)
    for m in sorted({int(x) for x in state.hm_labels if x != 0}):
        genes = state.module_genes(m)
        n_spots = int((assign == m).sum())
        print(f"   module {m:>3}: {len(genes):>5} genes, "
              f"{n_spots:>5} spots, top genes: {', '.join(genes[:5])}")

    print("\n" + "=" * 70)
    print("6. Explore a single module interactively (cellxgene)")
    print("=" * 70)
    print("   The grid figure above shows every module at once. To zoom into one")
    print("   module, load the gene list in cellxgene (https://cellxgene.cziscience.com):")
    print("     - open the Visium .h5ad for this dataset in cellxgene")
    print("     - paste a module's genes from spatial_out/HM_ModGene.csv (hm_mod == N)")
    print("     - e.g. `awk -F, '$1==1 {print $2}' spatial_out/HM_ModGene.csv`")
    print("   or colour the spots yourself from spatial_out/spatial_module_scores.tsv.")

    print(f"\nDone. Outputs in {args.outdir}/")


if __name__ == "__main__":
    main()
