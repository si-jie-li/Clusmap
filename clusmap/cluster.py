"""Hierarchical clustering + dynamic tree cut -> ModuleState."""
from __future__ import annotations

import os
import pickle
import time
from typing import Optional

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from dynamicTreeCut import cutreeHybrid

from .state import ModuleState


def gen_mod(
    rna_df: pd.DataFrame,
    *,
    method: str = "average",
    metric: str = "correlation",
    outdir: Optional[str] = ".",
    link: Optional[np.ndarray] = None,
    deepSplit: int = 1,
    minClusterSize: int = 30,
    pamStage: bool = False,
    save_raw: bool = False,
    **cutree_kwargs,
) -> ModuleState:
    """Cluster genes and cut the tree into modules.

    Only ``rna_df`` is required. ``deepSplit`` (0 coarse .. 4 fine) and
    ``minClusterSize`` are the two knobs users usually touch.

    Returns a :class:`ModuleState` (carries linkage, labels, gene order and
    supports split/merge edits). The legacy ``(link, mod, rna_df)`` values are
    available as ``state.linkage`` / ``{'labels': state.raw_labels}`` /
    ``rna_df``.

    Saved to ``outdir``: ``module_state.pkl`` (the full state — linkage +
    labels + gene order; reload with ``ModuleState.load``), ``module_state.json``
    (human-readable summary), ``ModGene.csv`` (original module numbering) and
    ``HM_ModGene.csv`` (heatmap numbering). Set ``save_raw=True`` to also dump
    the raw ``cutreeHybrid`` dict (``modules.pkl``) and the bare linkage matrix
    (``linkage.pkl``) — both are otherwise redundant with ``module_state.pkl``.
    """
    # de-duplicate labels so everything downstream stays aligned
    rna_df = rna_df.copy()
    rna_df.index = pd.Index(rna_df.index.astype(str)).where(
        ~pd.Index(rna_df.index.astype(str)).duplicated(),
        pd.Index(rna_df.index.astype(str)) + "_dup")
    rna_df = rna_df[~rna_df.index.duplicated(keep="first")]

    distance = pdist(rna_df.values, metric=metric)
    if link is None:
        print(">>> Hierarchical clustering ...", end="\r")
        t0 = time.time()
        link = linkage(distance, method=method)
        print(f"Hierarchical clustering took {time.time() - t0:.2f}s")

    print(">>> Dynamic tree cut ...", end="\r")
    t0 = time.time()
    args = {"deepSplit": deepSplit, "minClusterSize": minClusterSize, "pamStage": pamStage}
    args.update({k: v for k, v in cutree_kwargs.items()
                 if k in cutreeHybrid.__code__.co_varnames})
    mod = cutreeHybrid(link, distance, **args)
    print(f"Dynamic tree cut took {time.time() - t0:.2f}s")

    state = ModuleState.from_cutree(link, mod, rna_df.index, metric=metric,
                                    method=method, data=rna_df)
    print(f"Found {state.n_modules} modules "
          f"({int((state.raw_labels == 0).sum())} genes unassigned).")

    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        # original (un-renumbered) mapping, kept for reference
        pd.DataFrame({"module": state.raw_labels, "gene": rna_df.index}) \
            .sort_values(["module", "gene"]) \
            .to_csv(os.path.join(outdir, "ModGene.csv"), index=False)
        state.save(outdir)   # module_state.pkl (linkage+labels) + .json + HM_ModGene.csv
        if save_raw:         # redundant with module_state.pkl; off by default
            with open(os.path.join(outdir, "modules.pkl"), "wb") as fh:
                pickle.dump(mod, fh)
            with open(os.path.join(outdir, "linkage.pkl"), "wb") as fh:
                pickle.dump(link, fh)
        print(f"Clustering outputs saved to {outdir}/")

    return state
