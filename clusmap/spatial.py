"""Spatial-transcriptomics support: treat a spot as a pseudo-bulk sample.

This module lets you feed a 10x Visium dataset through the exact same clusmap
engine (gene x spot instead of gene x sample): hierarchical clustering + dynamic
tree cut to find gene modules, then annotate them as usual. The one genuinely
new capability is *spatial* — rendering each module's expression profile across
the tissue so you can map it back onto the H&E image.

Two spatial views are provided:

* ``plot_spatial_expression`` — one subplot per module, each showing that
  module's per-spot expression (mean of its genes, or the first-PC eigengene)
  as a continuous colour scale, optionally overlaid on the H&E image. This is
  the primary "where is this module active?" view.
* ``plot_spatial_modules`` — a single figure colouring each spot by the module
  it most strongly expresses (discrete, using the heatmap's module colours).

Spot -> module is scored via ``spatial_module_scores`` (spots x modules) and
assigned via ``assign_spots_to_modules``.

Nothing here re-implements the engine: ``spatial_module_scores`` delegates to
``analysis.module_eigengenes``, and the colors for the discrete view come from
``plot.module_color_map`` so they match ``bulk_hm`` exactly.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Union

import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix, issparse

from .state import ModuleState


# --------------------------------------------------------------------------- #
# container
# --------------------------------------------------------------------------- #
@dataclass
class SpatialDataset:
    """A spatial dataset as gene x spot bulk data plus physical coordinates.

    Attributes
    ----------
    rna:
        genes x spots expression matrix (columns = spot barcodes). Feeds straight
        into ``preprocess`` / ``gen_mod`` / ``bulk_hm`` like any bulk matrix.
    coords:
        index = spot barcode, columns ``x`` / ``y`` in full-resolution pixel
        coordinates (as written by Space Ranger ``tissue_positions``).
    scale_factors:
        ``scalefactors_json.json`` contents (``tissue_lowres_scalef``,
        ``tissue_hires_scalef``, ``spot_diameter_fullres``, ...).
    image:
        path to ``tissue_lowres_image.png`` (or an image array from AnnData),
        used as the H&E backdrop. ``None`` disables the overlay.
    image_hires:
        path to ``tissue_hires_image.png`` (full resolution) if present.
    """
    rna: pd.DataFrame
    coords: pd.DataFrame
    scale_factors: Dict[str, float]
    image: Optional[Union[str, np.ndarray]] = None
    image_hires: Optional[Union[str, np.ndarray]] = None


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _make_unique(names) -> list:
    """De-duplicate gene names in place (``X``, ``X_1``, ``X_2``, ...)."""
    seen, out = set(), []
    for n in map(str, names):
        base, k = n, 0
        while base in seen:
            k += 1
            base = f"{n}_{k}"
        seen.add(base)
        out.append(base)
    return out


def _first(cols, candidates) -> Optional[str]:
    for c in candidates:
        if c in cols:
            return c
    return None


def _decode(v) -> str:
    return v.decode() if isinstance(v, bytes) else str(v)


# --------------------------------------------------------------------------- #
# 10x matrix readers
# --------------------------------------------------------------------------- #
def _read_10x_h5(path: str):
    """Read ``filtered_feature_bc_matrix.h5`` -> (csc X, genes, barcodes)."""
    import h5py
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        shape = tuple(int(x) for x in m["shape"][:])          # (features, barcodes)
        X = csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]),
                       shape=shape)
        genes = [_decode(n) for n in m["features"]["name"][:]]
        barcodes = [_decode(b) for b in m["barcodes"][:]]
    return X, genes, barcodes


def _read_10x_mtx(path: str):
    """Read a ``filtered_feature_bc_matrix/`` mtx dir -> (csc X, genes, barcodes)."""
    import gzip
    from scipy.io import mmread

    def _col(p, i):
        with gzip.open(p, "rt") as fh:
            return [ln.rstrip("\n").split("\t")[i] for ln in fh]

    X = mmread(os.path.join(path, "matrix.mtx.gz"))
    if not issparse(X):
        X = csc_matrix(X)
    genes = _col(os.path.join(path, "features.tsv.gz"), 1)
    barcodes = _col(os.path.join(path, "barcodes.tsv.gz"), 0)
    return X, genes, barcodes


# --------------------------------------------------------------------------- #
# coordinate / scale-factor readers
# --------------------------------------------------------------------------- #
def _read_coords(path: str) -> pd.DataFrame:
    """Parse ``tissue_positions[_list].csv`` -> DataFrame (index=barcode, x, y).

    Handles both the 6-column no-header ``_list`` variant and the header'd
    ``tissue_positions.csv``; only ``in_tissue == 1`` spots are kept.
    """
    raw = pd.read_csv(path, header=None)
    first = str(raw.iloc[0, 0]).strip().lower()
    if first in ("barcode", "barcodes", "spot", "spot_id", "cell_id"):
        df = pd.read_csv(path)
        df.columns = [c.strip() for c in df.columns]
        bc = df.iloc[:, 0].astype(str)
        xcol = _first(df.columns, ["pxl_col_in_fullres", "pixel_x", "x", "col"])
        ycol = _first(df.columns, ["pxl_row_in_fullres", "pixel_y", "y", "row"])
        icol = _first(df.columns, ["in_tissue", "in tissue", "tissue"])
        in_tissue = df[icol].astype(int) if icol else pd.Series(1, index=df.index)
        x = df[xcol].astype(float)
        y = df[ycol].astype(float)
    else:
        # no header: barcode, in_tissue, array_row, array_col, pxl_col, pxl_row
        bc = raw.iloc[:, 0].astype(str)
        in_tissue = raw.iloc[:, 1].astype(int)
        x = raw.iloc[:, 4].astype(float)
        y = raw.iloc[:, 5].astype(float)
    coords = pd.DataFrame({"x": x.values, "y": y.values}, index=bc.values)
    coords = coords[in_tissue.values == 1]
    coords = coords[~coords.index.duplicated(keep="first")]
    return coords


def _read_scalefactors(path: Union[str, Path]) -> Dict[str, float]:
    p = Path(path)
    if not p.exists():
        return {}
    with open(p) as fh:
        return json.load(fh)


def _first_library(uns_spatial: dict) -> Optional[dict]:
    if "scalefactors" in uns_spatial or "images" in uns_spatial:
        return uns_spatial
    for v in uns_spatial.values():
        if isinstance(v, dict) and ("scalefactors" in v or "images" in v):
            return v
    return None


# --------------------------------------------------------------------------- #
# importers
# --------------------------------------------------------------------------- #
def import_spatial(path, *, counts_file=None, coords_file=None,
                   image: str = "lowres") -> SpatialDataset:
    """Load a 10x Visium output directory into a :class:`SpatialDataset`.

    Expects the Space Ranger layout under ``path``::

        path/
          filtered_feature_bc_matrix.h5      (or filtered_feature_bc_matrix/ mtx)
          spatial/tissue_positions_list.csv  (or tissue_positions.csv)
          spatial/scalefactors_json.json
          spatial/tissue_lowres_image.png    (or tissue_hires_image.png)

    ``image`` selects the H&E backdrop ("lowres" default, "hires", or ``None``
    to disable). ``counts_file`` / ``coords_file`` override auto-detection.
    """
    path = Path(path)

    # expression matrix -> genes x spots (dense)
    if counts_file is None:
        h5 = path / "filtered_feature_bc_matrix.h5"
        mtx = path / "filtered_feature_bc_matrix"
        if h5.exists():
            X, genes, spots = _read_10x_h5(str(h5))
        elif mtx.is_dir():
            X, genes, spots = _read_10x_mtx(str(mtx))
        else:
            raise FileNotFoundError(
                f"No 10x matrix under {path} (expected filtered_feature_bc_matrix.h5 "
                "or filtered_feature_bc_matrix/).")
    else:
        cf = Path(counts_file)
        X, genes, spots = _read_10x_mtx(str(cf)) if cf.is_dir() else _read_10x_h5(str(cf))

    genes = _make_unique(genes)
    rna = pd.DataFrame(X.toarray(), index=genes, columns=spots).astype(float)

    # coordinates (keep only spots present in the matrix)
    if coords_file is None:
        coords_file = path / "spatial" / "tissue_positions_list.csv"
        if not coords_file.exists():
            coords_file = path / "spatial" / "tissue_positions.csv"
    coords = _read_coords(str(coords_file))
    keep = [s for s in spots if s in coords.index]
    coords = coords.loc[keep]
    rna = rna.loc[:, keep]

    # scale factors + images
    scale_factors = _read_scalefactors(path / "spatial" / "scalefactors_json.json")
    lowres = path / "spatial" / "tissue_lowres_image.png"
    hires = path / "spatial" / "tissue_hires_image.png"
    image_path = None
    if image in ("lowres", "hires"):
        chosen = lowres if image == "lowres" else hires
        if chosen.exists():
            image_path = str(chosen)
    return SpatialDataset(rna=rna, coords=coords, scale_factors=scale_factors,
                          image=image_path,
                          image_hires=str(hires) if hires.exists() else None)


def from_adata(adata, *, image: str = "lowres") -> SpatialDataset:
    """Convert a scanpy ``read_visium`` AnnData into a :class:`SpatialDataset`.

    Uses ``adata.obsm["spatial"]`` for pixel coordinates and ``adata.uns`` for
    scale factors / images, so users who already load data with scanpy
    (e.g. ``sc.datasets.visium_sge(...)``) can hand off without re-parsing files.
    """
    genes = _make_unique([str(g) for g in adata.var_names])
    spots = [str(b) for b in adata.obs_names]
    X = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
    rna = pd.DataFrame(X.T, index=genes, columns=spots).astype(float)

    if "spatial" in adata.obsm:
        arr = np.asarray(adata.obsm["spatial"])
        coords = pd.DataFrame({"x": arr[:, 0], "y": arr[:, 1]}, index=spots)
    else:
        coords = pd.DataFrame({"x": np.arange(len(spots)),
                               "y": np.arange(len(spots))}, index=spots)

    scale_factors: Dict[str, float] = {}
    img: Optional[np.ndarray] = None
    img_hires: Optional[np.ndarray] = None
    lib = _first_library(adata.uns.get("spatial", {}))
    if lib:
        scale_factors = dict(lib.get("scalefactors", {}) or {})
        images = lib.get("images", {}) or {}
        if image in ("lowres", "hires") and image in images:
            img = images[image]
        if "hires" in images:
            img_hires = images["hires"]
    return SpatialDataset(rna=rna, coords=coords, scale_factors=scale_factors,
                          image=img, image_hires=img_hires)


# --------------------------------------------------------------------------- #
# feature selection
# --------------------------------------------------------------------------- #
def select_hvgs(rna: pd.DataFrame, *, n_top: int = 2000,
                flavor: str = "dispersion") -> pd.DataFrame:
    """Subset to the top ``n_top`` highly-variable genes (genes x spots).

    A simple, pure-pandas heuristic so clustering the full 30k-gene spatial
    matrix stays tractable (``pdist`` is O(n^2) in genes). ``flavor``:

    * ``"dispersion"`` (default) — variance-to-mean ratio of the raw counts
      (index of dispersion), suited to unfiltered count matrices.
    * ``"variance"`` — plain per-gene variance.

    Returns ``n_top`` or fewer rows; pass ``n_top=None`` to return unchanged.
    """
    if n_top is None or n_top >= rna.shape[0]:
        return rna.copy()
    means = rna.mean(axis=1)
    var = rna.var(axis=1)
    if flavor == "dispersion":
        score = var / means.replace(0, np.nan)
    elif flavor == "variance":
        score = var
    else:
        raise ValueError(f"Unknown flavor {flavor!r} (dispersion|variance).")
    top = score.fillna(0.0).sort_values(ascending=False).head(n_top).index
    return rna.loc[top].copy()


# --------------------------------------------------------------------------- #
# module scores / assignment
# --------------------------------------------------------------------------- #
def _norm_scores(scores: pd.DataFrame, norm: Optional[str]) -> pd.DataFrame:
    if norm in (None, "none"):
        return scores
    if norm == "zscore":
        sd = scores.std(axis=0)
        return ((scores - scores.mean(axis=0)) / sd.replace(0, np.nan)).fillna(0.0)
    if norm == "minmax":
        rng = scores.max(axis=0) - scores.min(axis=0)
        return ((scores - scores.min(axis=0)) / rng.replace(0, np.nan)).fillna(0.0)
    raise ValueError(f"Unknown norm {norm!r} (zscore|minmax|None).")


def spatial_module_scores(rna: pd.DataFrame, state: ModuleState, *,
                          method: str = "mean",
                          norm: Optional[str] = "zscore") -> pd.DataFrame:
    """Per-spot expression of each gene module (spots x modules).

    ``method``:

    * ``"mean"`` (default) — mean expression of the module's genes at each spot
      (a direct "spatial expression profile" of the module).
    * ``"eigengene"`` — the module's first-PC eigengene (delegates to
      :func:`analysis.module_eigengenes`), a WGCNA-style representative profile.

    ``norm`` z-scores (default) or min-max scales each module across spots so
    modules are comparable on one colour scale; ``None`` returns raw values.
    Columns are 1-based heatmap module ids, matching ``state`` / the heatmap.
    """
    from .analysis import module_eigengenes

    rna = rna.loc[[g for g in state.genes if g in rna.index]]
    modules = sorted({int(x) for x in state.hm_labels if x != 0})
    if method == "eigengene":
        eig = module_eigengenes(rna, state, outdir=None)
        scores = eig.T.reindex(rna.columns)
    elif method == "mean":
        data: Dict[int, Any] = {}
        for m in modules:
            genes = [g for g in state.module_genes(m) if g in rna.index]
            data[m] = (rna.loc[genes].mean(axis=0) if genes
                       else pd.Series(np.nan, index=rna.columns))
        scores = pd.DataFrame(data, index=rna.columns)
    else:
        raise ValueError(f"Unknown method {method!r} (mean|eigengene).")
    scores.columns.name = "module"
    scores.index.name = "spot"
    return _norm_scores(scores, norm)


def assign_spots_to_modules(scores: pd.DataFrame, *,
                            min_score: Optional[float] = None) -> pd.Series:
    """Assign each spot to its highest-scoring module (Series: spot -> module id).

    ``min_score`` (optional) thresholds on the best score; spots below it become
    0 (unassigned / mixed). ``None`` always assigns the argmax module.
    """
    assign = scores.idxmax(axis=1).astype(int)
    assign.name = "module"
    if min_score is not None:
        assign[scores.max(axis=1) < min_score] = 0
    return assign


# --------------------------------------------------------------------------- #
# plotting helpers
# --------------------------------------------------------------------------- #
def _load_image(image):
    from matplotlib import image as mpimg
    if isinstance(image, (str, os.PathLike)):
        return mpimg.imread(str(image))
    return np.asarray(image)


def _image_scalef(image, scale_factors) -> float:
    """Pick the coordinate->image scale factor from the image identity.

    Lowres images are downsampled by ``tissue_lowres_scalef``; hires images are
    full-resolution (scale 1). Array images without a filename are assumed
    lowres (the ``from_adata(image="lowres")`` default).
    """
    name = str(image) if isinstance(image, (str, os.PathLike)) else ""
    if "hires" in name:
        return 1.0
    return (scale_factors or {}).get("tissue_lowres_scalef", 1.0)


def _prepare_coords(coords: pd.DataFrame, image, scale_factors):
    """Return (px, py, overlay) display coordinates for a spatial scatter."""
    x = coords["x"].astype(float).values
    y = coords["y"].astype(float).values
    if image is None:
        return x, y, False
    sf = _image_scalef(image, scale_factors)
    return x / sf, y / sf, True


def _draw_spatial_base(ax, coords, image, scale_factors, spot_size):
    """Draw the H&E backdrop (if any) + axis limits; return (px, py, spot_size)."""
    px, py, overlay = _prepare_coords(coords, image, scale_factors)
    if overlay:
        arr = _load_image(image)
        h, w = arr.shape[:2]
        ax.imshow(arr, aspect="equal")
        ax.set_xlim(0, w)
        ax.set_ylim(h, 0)                       # top-left origin, like the image
        if spot_size is None:
            sf = _image_scalef(image, scale_factors)
            diam = (scale_factors or {}).get("spot_diameter_fullres", 89.5) * sf
            spot_size = max(1.0, diam ** 2)
    else:
        ax.invert_yaxis()                        # image convention (y down)
        spot_size = 6.0 if spot_size is None else spot_size
    return px, py, spot_size


# --------------------------------------------------------------------------- #
# plotting: discrete module membership
# --------------------------------------------------------------------------- #
def plot_spatial_modules(coords, assignment, *, state=None, palette="hsv",
                         color_map=None, image=None, scale_factors=None,
                         spot_size=None, title=None, legend=True, ax=None,
                         outdir=None, save_format="png", dpi=300):
    """Colour each spot by the module it most strongly expresses.

    ``assignment`` is a spot -> module Series (see ``assign_spots_to_modules``).
    Colours come from ``color_map`` or, if omitted, ``module_color_map(state)``
    so they match the heatmap exactly (0 / unassigned -> white). Pass ``image`` +
    ``scale_factors`` to overlay the H&E image.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    if color_map is None:
        if state is None:
            raise ValueError("Pass either color_map or state (for module colors).")
        from .plot import module_color_map
        color_map = module_color_map(state, palette)

    common = coords.index.intersection(assignment.index)
    coords = coords.loc[common]
    assignment = assignment.loc[common]

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    else:
        fig = ax.get_figure()

    px, py, s = _draw_spatial_base(ax, coords, image, scale_factors, spot_size)
    mods = assignment.values.astype(int)
    colors = [color_map.get(m, "#ffffff") for m in mods]
    ax.scatter(px, py, c=colors, s=s, marker="o", linewidths=0, alpha=0.9)

    if legend:
        present = [m for m in sorted(set(mods)) if m != 0]
        if 0 in set(mods):
            present.append(0)
        handles = [Patch(facecolor=color_map.get(m, "#ffffff"),
                         label=("unassigned" if m == 0 else f"module {m}"))
                   for m in present]
        ax.legend(handles=handles, title="Modules", bbox_to_anchor=(1.02, 1),
                  loc="upper left", frameon=False, fontsize=8)
    ax.set_title(title if title is not None
                 else f"Spatial modules ({len(set(mods) - {0})} modules)")
    ax.set_xticks([])
    ax.set_yticks([])

    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        fname = f"spatial_modules.{save_format.lower()}"
        fig.savefig(os.path.join(outdir, fname), dpi=dpi, bbox_inches="tight")
        print(f"Spatial module map saved to {os.path.join(outdir, fname)}")
    return fig, ax


# --------------------------------------------------------------------------- #
# plotting: per-module expression grid
# --------------------------------------------------------------------------- #
def plot_spatial_expression(scores, coords, *, image=None, scale_factors=None,
                            ncols=3, cmap="viridis", spot_size=None,
                            vmin=None, vmax=None, colorbar=True, title=None,
                            figsize=None, outdir=None, save_format="png",
                            dpi=150):
    """One subplot per module showing its spatial expression profile.

    ``scores`` is the spots x modules table from ``spatial_module_scores``; each
    subplot scatters the spots, coloured by that module's score, over the
    optional H&E image. ``vmin``/``vmax`` set a shared colour scale (defaults to
    the global min/max so modules are directly comparable). Returns ``(fig, axes)``.
    """
    import matplotlib.pyplot as plt

    modules = list(scores.columns)
    n = len(modules)
    ncols = max(1, min(ncols, n))
    nrows = int(np.ceil(n / ncols))
    if figsize is None:
        figsize = (ncols * 3.4, nrows * 3.4)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)

    if vmin is None:
        vmin = float(scores.min().min())
    if vmax is None:
        vmax = float(scores.max().max())

    sc = None
    for i, m in enumerate(modules):
        ax = axes.flat[i]
        px, py, s = _draw_spatial_base(ax, coords, image, scale_factors, spot_size)
        sc = ax.scatter(px, py, c=scores[m].values.astype(float), cmap=cmap,
                        s=s, marker="o", linewidths=0, vmin=vmin, vmax=vmax)
        ax.set_title(f"Module {m}")
        ax.set_xticks([])
        ax.set_yticks([])

    for j in range(n, nrows * ncols):
        axes.flat[j].set_visible(False)

    if colorbar and sc is not None:
        fig.colorbar(sc, ax=axes.ravel().tolist(), shrink=0.8,
                     label="module expression (z-score)")
    fig.suptitle(title or "Spatial expression of gene modules", fontsize=13)

    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        fname = f"spatial_module_expression.{save_format.lower()}"
        fig.savefig(os.path.join(outdir, fname), dpi=dpi, bbox_inches="tight")
        print(f"Spatial module expression saved to {os.path.join(outdir, fname)}")
    return fig, axes
