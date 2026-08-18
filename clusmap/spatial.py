"""Spatial-transcriptomics support (scanpy-centric, 10x Visium).

A spatial dataset is treated as pseudo-bulk: a spot is a sample, so the clusmap
engine — gene modules via hierarchical clustering + dynamic tree cut, then the
clusterheatmap / annotation / analysis steps — applies unchanged to the
**gene x spot** matrix. The spatial part itself is scanpy-native:

* ``import_spatial`` reads either a 10x Visium output *folder*
  (``sc.read_visium``) or a pre-made ``.h5ad`` *file* (``sc.read_h5ad``) into an
  AnnData and, when the data is raw (no ``leiden`` / ``log1p`` / HVG already
  present), runs the standard scanpy pipeline
  (normalize -> log1p -> HVG -> scale -> PCA -> neighbors -> Leiden) so that
  ``adata.obs['leiden']`` holds **spot clusters** (a module formed by spots).
* ``run_leiden`` is that pipeline, callable on any AnnData.
* ``plot_spatial_modules`` renders the **Leiden** spot clusters over the H&E
  image with ``sc.pl.spatial(adata, color='leiden', img_key='hires')``.
* ``plot_spatial_expression`` renders the **mean expression of each gene
  module** (a module formed by genes, from ``gen_mod``) across the tissue, one
  subplot per module, again via ``sc.pl.spatial``.
* ``spatial_hm`` draws the usual clusmap clusterheatmap in **two versions**
  reusing :func:`clusmap.plot.bulk_hm`: v1 clusters rows *and* columns;
  v2 groups columns by Leiden cluster (clustered within each group) with
  ``col_cluster=False``. The Leiden annotation is the default column colour
  band; pass ``col_cat`` to add more bands.

The spatial metadata (pixel coordinates, H&E image, scale factors) all lives in
the AnnData (``obsm['spatial']`` + ``uns['spatial']``) so nothing is
re-implemented — plots are plain ``sc.pl.spatial`` calls.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.sparse import issparse

from .state import ModuleState

try:  # scanpy is an optional (spatial) dependency
    import scanpy as sc
    _HAVE_SCANPY = True
except Exception:  # pragma: no cover - exercised only without the extra
    sc = None
    _HAVE_SCANPY = False


# --------------------------------------------------------------------------- #
# container
# --------------------------------------------------------------------------- #
@dataclass
class SpatialDataset:
    """A spatial dataset as an AnnData plus the clusmap gene x spot matrix.

    Attributes
    ----------
    adata:
        The (preprocessed) AnnData. ``obs['leiden']`` holds the spot clusters,
        ``obsm['spatial']`` the pixel coordinates, ``uns['spatial']`` the H&E
        images / scale factors. This is what the ``sc.pl.spatial`` plots and
        the two-version heatmap read from.
    rna:
        genes x spots expression matrix (log1p of the highly-variable genes;
        columns = spot barcodes). Feeds straight into ``gen_mod`` / ``bulk_hm``
        like any bulk matrix.
    coords:
        index = spot barcode, columns ``x`` / ``y`` (pixel coordinates, taken
        from ``adata.obsm['spatial']``). Provided for convenience.
    scale_factors:
        ``scalefactors_json.json`` contents (``tissue_lowres_scalef``,
        ``tissue_hires_scalef``, ``spot_diameter_fullres``, ...).
    image:
        lowres H&E image array from ``uns['spatial']`` (or ``None``).
    image_hires:
        hires H&E image array from ``uns['spatial']`` (or ``None``).
    """

    adata: Any
    rna: pd.DataFrame
    coords: pd.DataFrame
    scale_factors: Dict[str, float]
    image: Optional[np.ndarray] = None
    image_hires: Optional[np.ndarray] = None

    @property
    def leiden(self) -> pd.Series:
        """Spot -> Leiden cluster (a module formed by spots)."""
        return self.adata.obs["leiden"].astype(str)


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _require_scanpy():
    if not _HAVE_SCANPY:
        raise ImportError(
            "scanpy is required for spatial support. Install it with: "
            "pip install 'clusmap[spatial]'")


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


def _as_int(x) -> int:
    try:
        return int(x)
    except (TypeError, ValueError):
        return 0


def _first(cols, candidates) -> Optional[str]:
    for c in candidates:
        if c in cols:
            return c
    return None


def _read_coords(path: str) -> pd.DataFrame:
    """Parse ``tissue_positions[_list].csv`` -> DataFrame (index=barcode, x, y).

    Handles both the 6-column no-header ``_list`` variant and the header'd
    ``tissue_positions.csv``; only ``in_tissue == 1`` spots are kept. Mainly
    useful for inspecting a Space Ranger output directory directly.
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


def _library_id(adata) -> Optional[str]:
    """Return the first Visium library id in ``uns['spatial']`` (or ``None``)."""
    spatial = adata.uns.get("spatial", {}) or {}
    if not spatial:
        return None
    first = next(iter(spatial))
    if isinstance(spatial[first], dict) and (
            "images" in spatial[first] or "scalefactors" in spatial[first]):
        return first          # nested per-library layout
    return None               # already the single-library layout


def _has_images(adata) -> bool:
    lib = _library_id(adata)
    if lib is None:
        spatial = adata.uns.get("spatial", {}) or {}
        return bool(spatial.get("images"))
    return bool((adata.uns["spatial"][lib] or {}).get("images"))


def _rna_from_adata(adata: Any) -> pd.DataFrame:
    """log1p gene x spot DataFrame from an AnnData (HVG subset when marked)."""
    _require_scanpy()
    adata = adata.copy()
    adata.var_names_make_unique()
    hvg = adata.var.get("highly_variable")
    if hvg is not None and bool(hvg.any()):
        adata = adata[:, hvg]
    X = adata.X
    X = X.toarray() if issparse(X) else np.asarray(X)
    return pd.DataFrame(X.T, index=list(adata.var_names),
                        columns=list(adata.obs_names)).astype(float)


# --------------------------------------------------------------------------- #
# preprocessing + Leiden (spot clusters)
# --------------------------------------------------------------------------- #
def run_leiden(adata, *, n_top_genes: int = 2000, resolution: float = 0.5,
               random_state: int = 0, n_neighbors: int = 10,
               n_pcs: int = 40) -> Any:
    """Run the standard scanpy pipeline, adding ``obs['leiden']`` to *adata*.

    If the data is already preprocessed (``log1p`` in ``uns``), the
    normalization steps are skipped; if ``highly_variable`` is already in
    ``var``, the existing HVG mask is respected; if ``leiden`` is already in
    ``obs`` the data is returned untouched. Returns the same AnnData (a copy),
    with ``obs['leiden']`` holding the spot clusters.
    """
    _require_scanpy()
    if "leiden" in adata.obs:
        return adata
    adata = adata.copy()
    adata.var_names_make_unique()

    if "log1p" not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    if "highly_variable" not in adata.var:
        if n_top_genes is None or n_top_genes >= adata.n_vars:
            adata.var["highly_variable"] = True
        else:
            sc.pp.highly_variable_genes(adata, flavor="seurat",
                                        n_top_genes=n_top_genes)

    hvg = adata[:, adata.var["highly_variable"].values].copy()
    sc.pp.scale(hvg, max_value=10)
    n_comps = min(hvg.n_obs - 1, hvg.n_vars)
    sc.tl.pca(hvg, svd_solver="arpack", n_comps=min(n_pcs, n_comps))
    sc.pp.neighbors(hvg, n_neighbors=min(n_neighbors, hvg.n_obs - 1),
                    n_pcs=min(n_pcs, n_comps))
    sc.tl.leiden(hvg, resolution=resolution, random_state=random_state)
    adata.obs["leiden"] = hvg.obs["leiden"]
    return adata


# --------------------------------------------------------------------------- #
# importers
# --------------------------------------------------------------------------- #
def import_spatial(path, *, n_top_genes: int = 2000, resolution: float = 0.5,
                   random_state: int = 0) -> SpatialDataset:
    """Load spatial data from a Visium folder or an ``.h5ad`` file.

    Robust to both input formats:

    * a **10x Visium output folder** (``sc.read_visium``) — e.g.
      ``"V1_Mouse_Brain_Sagittal_Posterior"`` (must contain
      ``filtered_feature_bc_matrix.h5`` + a ``spatial/`` directory);
    * a **pre-made ``.h5ad`` file** (``sc.read_h5ad``) — a Visium AnnData with
      ``obsm['spatial']`` / ``uns['spatial']`` (e.g. produced by
      ``sc.datasets.visium_sge`` and saved, or a custom dataset).

    If the input carries no preprocessed information (no ``leiden``), the
    scanpy pipeline (``run_leiden``) is run exactly as for raw data, so the
    returned dataset always has ``adata.obs['leiden']``.

    Returns a :class:`SpatialDataset` (``adata`` + the ``rna`` gene x spot
    matrix for the clusmap engine).
    """
    _require_scanpy()
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Spatial data not found: {path}\nPass either a 10x Visium output "
            "folder or a .h5ad file.")
    if path.name.lower().endswith((".h5ad", ".h5ad.gz")):
        adata = sc.read_h5ad(str(path))
    elif path.is_dir():
        adata = sc.read_visium(str(path))
    else:
        raise ValueError(
            f"Unrecognized spatial input {path}: pass a Visium output folder "
            "(a directory) or a .h5ad file (a 10x filtered_feature_bc_matrix.h5 "
            "belongs to a folder, not a standalone .h5ad).")
    return from_adata(adata, n_top_genes=n_top_genes, resolution=resolution,
                      random_state=random_state)


def from_adata(adata, *, n_top_genes: int = 2000, resolution: float = 0.5,
               random_state: int = 0) -> SpatialDataset:
    """Convert any (Visium) AnnData into a :class:`SpatialDataset`.

    Ensures preprocessing + Leiden (via :func:`run_leiden`) unless the data is
    already annotated, then builds the ``rna`` gene x spot matrix and the
    convenience ``coords`` / scale-factor fields.
    """
    _require_scanpy()
    adata = run_leiden(adata, n_top_genes=n_top_genes, resolution=resolution,
                       random_state=random_state)
    rna = _rna_from_adata(adata)

    spots = [str(b) for b in adata.obs_names]
    if "spatial" in adata.obsm:
        arr = np.asarray(adata.obsm["spatial"])
        coords = pd.DataFrame({"x": arr[:, 0], "y": arr[:, 1]}, index=spots)
    else:
        coords = pd.DataFrame({"x": np.arange(len(spots)),
                               "y": np.arange(len(spots))}, index=spots)

    scale_factors: Dict[str, float] = {}
    img: Optional[np.ndarray] = None
    img_hires: Optional[np.ndarray] = None
    lib = _library_id(adata)
    if lib is not None:
        entry = adata.uns["spatial"][lib]
    else:
        entry = adata.uns.get("spatial", {})
    if entry:
        scale_factors = dict(entry.get("scalefactors", {}) or {})
        images = entry.get("images", {}) or {}
        if "lowres" in images:
            img = np.asarray(images["lowres"])
        if "hires" in images:
            img_hires = np.asarray(images["hires"])
    return SpatialDataset(adata=adata, rna=rna, coords=coords,
                          scale_factors=scale_factors, image=img,
                          image_hires=img_hires)


# --------------------------------------------------------------------------- #
# feature selection (pure-pandas fallback / convenience)
# --------------------------------------------------------------------------- #
def select_hvgs(rna: pd.DataFrame, *, n_top: int = 2000,
                flavor: str = "dispersion") -> pd.DataFrame:
    """Subset to the top ``n_top`` highly-variable genes (genes x spots).

    A simple, pure-pandas heuristic so clustering a large matrix stays
    tractable (``pdist`` is O(n^2) in genes). ``flavor``:

    * ``"dispersion"`` (default) — variance-to-mean ratio of the raw counts.
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
# module scores / assignment (gene modules on the rna matrix)
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
    """Per-spot expression of each *gene* module (spots x modules).

    ``method``:

    * ``"mean"`` (default) — mean expression of the module's genes at each spot
      (the "spatial expression profile" of the gene module).
    * ``"eigengene"`` — the module's first-PC eigengene (delegates to
      :func:`clusmap.analysis.module_eigengenes`), a WGCNA-style profile.

    ``norm`` z-scores (default) or min-max scales each module across spots so
    modules are comparable; ``None`` returns raw values. Columns are 1-based
    heatmap module ids, matching ``state`` / the heatmap.
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
    """Assign each spot to its highest-scoring *gene* module (spot -> module id).

    ``min_score`` (optional) thresholds on the best score; spots below it become
    0 (unassigned / mixed). ``None`` always assigns the argmax module.
    """
    assign = scores.idxmax(axis=1).astype(int)
    assign.name = "module"
    if min_score is not None:
        assign[scores.max(axis=1) < min_score] = 0
    return assign


def add_module_expression(adata, rna_df, state, *, method: str = "mean",
                          norm: Optional[str] = "zscore") -> pd.DataFrame:
    """Add per-gene-module expression to ``adata.obs`` and return the scores.

    Adds one ``obs['module_{id}_expr']`` column per gene module (z-scored by
    default so modules are comparable) — the columns that
    :func:`plot_spatial_expression` renders with ``sc.pl.spatial``. The scores
    DataFrame (spots x modules) is also returned.
    """
    scores = spatial_module_scores(rna_df, state, method=method, norm=norm)
    for m in scores.columns:
        adata.obs[f"module_{int(m)}_expr"] = scores[int(m)].reindex(
            adata.obs_names).values.astype(float)
    return scores


# --------------------------------------------------------------------------- #
# plotting: sc.pl.spatial over the H&E image
# --------------------------------------------------------------------------- #
def _spatial_fig(adata, color, *, img_key="hires", library_id=None,
                 ncols=1, title=None, **kwargs):
    _require_scanpy()
    if "spatial" not in adata.obsm:
        raise ValueError(
            "adata.obsm['spatial'] is missing — this AnnData has no Visium "
            "spatial coordinates; sc.pl.spatial cannot draw it.")
    if not _has_images(adata):
        print("[spatial] no H&E images found in adata.uns['spatial'] — "
              "sc.pl.spatial will still draw the tissue coordinates.")
    lib = library_id or _library_id(adata)
    return sc.pl.spatial(adata, color=color, img_key=img_key, library_id=lib,
                         ncols=ncols, title=title, show=False,
                         return_fig=True, **kwargs)


def plot_spatial_modules(adata, *, img_key: str = "hires", library_id=None,
                         title: Optional[str] = None, outdir: Optional[str] = ".",
                         save_format: str = "png", dpi: int = 150,
                         **kwargs):
    """Plot the **Leiden** spot clusters over the H&E image.

    Uses ``sc.pl.spatial(adata, color='leiden', img_key='hires')`` — a module
    formed *by spots*. Saves ``spatial_modules.<fmt>`` when ``outdir`` is given
    and returns the matplotlib ``Figure``.
    """
    _require_scanpy()
    if "leiden" not in adata.obs:
        raise ValueError(
            "adata.obs['leiden'] is missing — run import_spatial / run_leiden "
            "first (or pass an already-annotated AnnData).")
    if title is None:
        title = f"Leiden clusters (n = {adata.obs['leiden'].nunique()})"
    fig = _spatial_fig(adata, "leiden", img_key=img_key,
                       library_id=library_id, title=title, **kwargs)
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        fname = f"spatial_modules.{save_format.lower()}"
        fig.savefig(os.path.join(outdir, fname), dpi=dpi, bbox_inches="tight")
        print(f"Leiden spatial map saved to {os.path.join(outdir, fname)}")
    return fig


def plot_spatial_expression(adata, *, modules: Optional[List[int]] = None,
                            ncols: int = 3, img_key: str = "hires",
                            library_id=None, outdir: Optional[str] = ".",
                            save_format: str = "png", dpi: int = 150,
                            **kwargs):
    """Plot the **mean expression profile of each gene module** over the tissue.

    One subplot per gene module (a module formed *by genes*, from ``gen_mod``),
    coloured by the module's mean per-spot expression
    (``obs['module_{id}_expr']`` from :func:`add_module_expression`). Saves
    ``spatial_module_expression.<fmt>`` when ``outdir`` is given and returns
    the matplotlib ``Figure``.
    """
    _require_scanpy()
    cols = [c for c in adata.obs.columns
            if c.startswith("module_") and c.endswith("_expr")]
    if not cols:
        raise ValueError(
            "No 'module_*_expr' columns in adata.obs — run "
            "add_module_expression(adata, rna, state) first.")
    if modules is not None:
        cols = [f"module_{int(m)}_expr" for m in modules]
        missing = [c for c in cols if c not in adata.obs]
        if missing:
            raise ValueError(f"Missing module-expression columns: {missing}")
    ids = [int(c.split("_")[1]) for c in cols]
    titles = [f"Module {m}" for m in ids]
    fig = _spatial_fig(adata, cols, img_key=img_key, library_id=library_id,
                       ncols=ncols, title=titles, **kwargs)
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        fname = f"spatial_module_expression.{save_format.lower()}"
        fig.savefig(os.path.join(outdir, fname), dpi=dpi, bbox_inches="tight")
        print(f"Spatial module expression saved to {os.path.join(outdir, fname)}")
    return fig


# --------------------------------------------------------------------------- #
# two-version clusterheatmap (reuses bulk_hm)
# --------------------------------------------------------------------------- #
def sort_columns_by_leiden(rna_df: pd.DataFrame, leiden,
                           *, within_cluster: bool = True,
                           metric: str = "correlation",
                           method: str = "average"):
    """Order ``rna_df`` columns by Leiden group, clustered within each group.

    Returns ``(rna_df_sorted, leiden_sorted)``: the same genes, with the spot
    columns rearranged so that spots of one Leiden cluster sit together. When
    ``within_cluster`` (default), columns inside a group are re-ordered by their
    own correlation dendrogram, so each group still shows internal structure.
    """
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import pdist

    leiden = pd.Series(leiden).reindex(rna_df.columns)
    groups = sorted({c for c in leiden.dropna().unique()}, key=_as_int)
    order: List[Any] = []
    for g in groups:
        cols = [c for c in rna_df.columns if leiden.loc[c] == g]
        if within_cluster and len(cols) > 1:
            sub = rna_df[cols].astype(float)
            d = np.nan_to_num(pdist(sub.T.values, metric=metric))
            Z = linkage(d, method=method)
            order.extend(cols[i] for i in leaves_list(Z))
        else:
            order.extend(cols)
    return rna_df[order], leiden.reindex(order)


def spatial_hm(adata, rna_df, state: ModuleState, *, leiden_col: str = "leiden",
               versions=(1, 2),
               col_cat: Optional[Dict[str, pd.Series]] = None,
               col_color_manual: Optional[Dict[str, Dict[str, str]]] = None,
               col_legend: Optional[list] = None,
               hm_args: Optional[Dict[str, Any]] = None,
               title: Optional[str] = None, mod_palette: str = "hsv",
               col_palette: str = "tab20", save_format: str = "png",
               outdir: Optional[str] = "."):
    """Gene x spot clusterheatmap in two versions (reuses ``bulk_hm``).

    The **Leiden** annotation is the default column colour band
    (``col_cat={'leiden': ...}``); pass ``col_cat`` to add more bands
    (e.g. a sample-level metadata Series, exactly like ``bulk_hm``).

    * **v1** — rows and columns both clustered (``col_cluster=True``).
    * **v2** — rows clustered, columns grouped by Leiden cluster and clustered
      *within* each group (``col_cluster=False`` with pre-sorted columns).

    Returns ``(hm_v1, hm_v2)`` (the seaborn ``ClusterGrid`` objects). When
    ``outdir`` is given the figures are saved as ``spatial_heatmap_v1.<fmt>``
    and ``spatial_heatmap_v2.<fmt>``.
    """
    from .plot import bulk_hm

    if leiden_col not in adata.obs:
        raise ValueError(
            f"adata.obs[{leiden_col!r}] is missing — run import_spatial / "
            "run_leiden first.")
    leiden = adata.obs[leiden_col].astype(str).reindex(rna_df.columns)

    col_cat = dict(col_cat or {})
    col_cat.setdefault(leiden_col, leiden)
    col_legend = [leiden_col] + list(col_legend or [])

    hms: Dict[int, Any] = {}
    if 1 in versions:
        hms[1] = bulk_hm(rna_df, state, title=title, mod_palette=mod_palette,
                         col_palette=col_palette, col_cat=col_cat,
                         col_color_manual=col_color_manual,
                         col_legend=col_legend, hm_args=dict(hm_args or {}),
                         outdir=None)
    if 2 in versions:
        rna2, leiden2 = sort_columns_by_leiden(rna_df, leiden)
        cc2 = dict(col_cat)
        cc2[leiden_col] = leiden2
        hm_args2 = dict(hm_args or {})
        hm_args2["col_cluster"] = False
        hms[2] = bulk_hm(rna2, state, title=title, mod_palette=mod_palette,
                         col_palette=col_palette, col_cat=cc2,
                         col_color_manual=col_color_manual,
                         col_legend=col_legend, hm_args=hm_args2, outdir=None)

    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        for v, hm in hms.items():
            dpi = getattr(hm, "_clusmap_dpi", 300)
            fname = f"spatial_heatmap_v{v}.{save_format.lower()}"
            kw = {"dpi": dpi} if save_format.lower() == "png" else {}
            hm.savefig(os.path.join(outdir, fname), bbox_inches="tight", **kw)
            print(f"Spatial heatmap v{v} saved to {os.path.join(outdir, fname)}")
    return hms.get(1), hms.get(2)
