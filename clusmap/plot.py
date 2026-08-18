"""Bulk clustermap rendering and per-cluster sample statistics."""
from __future__ import annotations

import os
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

from .state import ModuleState
from .util import canon_norm, normalize_rows


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _module_colors(state: ModuleState, palette: str):
    """Per-gene RGB colors (original order) + hex array, 0 -> white."""
    labels = state.raw_labels
    uniq = sorted(set(labels))
    pal = sns.color_palette(palette, n_colors=len(uniq))
    lut = dict(zip(uniq, pal))
    lut[0] = (1, 1, 1)
    rgb = [lut[x] for x in labels]
    hex_arr = np.array([mcolors.to_hex(c) for c in rgb])
    return rgb, hex_arr


def module_color_map(state: ModuleState, palette: str = "hsv") -> Dict[int, str]:
    """``{hm_module_id: hex}`` colour map, consistent with ``bulk_hm``'s colours.

    1-based heatmap module id -> hex colour (0 / unassigned -> ``"#ffffff"``).
    Reuses the same palette/ordering as ``bulk_hm`` so any spatial rendering
    coloured with this map matches the heatmap exactly.
    """
    _, hex_arr = _module_colors(state, palette)
    cmap: Dict[int, str] = {}
    for m, c in zip(state.hm_labels, hex_arr):
        cmap.setdefault(int(m), str(c))
    cmap[0] = "#ffffff"
    return cmap


def _shrink_xticklabels(ax, angle: int) -> float:
    """Shrink x tick font until adjacent labels stop overlapping. Returns fs."""
    fig = ax.get_figure()
    fig.canvas.draw()
    ticks = ax.get_xticks()
    if len(ticks) < 2 or not ax.get_xticklabels():
        return 10.0
    renderer = fig.canvas.get_renderer()

    def overlaps(fs):
        ax.tick_params(axis="x", labelsize=fs)
        fig.canvas.draw()
        t = ax.get_xticklabels()[0]
        _, h, _ = renderer.get_text_width_height_descent(
            t.get_text() or "X", t.get_fontproperties(), ismath=False)
        step = abs(ax.transData.transform((ticks[1], 0))[0]
                   - ax.transData.transform((ticks[0], 0))[0])
        return step * np.sin(np.deg2rad(angle)) < h

    fs = 10.0
    for _ in range(30):
        if not overlaps(fs) or fs < 0.5:
            break
        fs -= 0.5
    return fs


def _build_col_colors(rna_df, col_cat, col_color_manual, palette):
    """Return (col_colors_df, {band: {cat: hex}}) for clustermap col_colors."""
    if col_cat is None:
        return None, {}
    color_map: Dict[str, Dict[str, str]] = {}
    data = {}
    for band, v in col_cat.items():
        s = (v.copy() if isinstance(v, pd.Series)
             else pd.Series(list(v), index=rna_df.columns)).astype("category")
        cats = list(s.cat.categories)
        if col_color_manual and band in col_color_manual:
            cmap = {k: (mcolors.to_hex(c) if not str(c).startswith("#") else c)
                    for k, c in col_color_manual[band].items()}
            missing = [c for c in cats if c not in cmap]
            for c, rgb in zip(missing, sns.color_palette(palette, len(missing))):
                cmap[c] = mcolors.to_hex(rgb)
        else:
            cmap = {c: mcolors.to_hex(rgb)
                    for c, rgb in zip(cats, sns.color_palette(palette, len(cats)))}
        color_map[band] = cmap
        data[band] = s.map(cmap).astype("object").fillna("#ffffff").values
    return pd.DataFrame(data, index=rna_df.columns), color_map


def _dpi_from_fontsize(fs: float) -> int:
    for thr, dpi in ((8, 300), (6, 400), (4, 600), (2, 800)):
        if fs >= thr:
            return dpi
    return 1000


# --------------------------------------------------------------------------- #
# main entry
# --------------------------------------------------------------------------- #
def bulk_hm(
    rna_df: pd.DataFrame,
    state: ModuleState,
    *,
    norm_method: str = "z_score",
    title: Optional[str] = None,
    hm_args: Optional[Dict[str, Any]] = None,
    xtickangle: int = 90,
    yticksize: int = 5,
    goi_list: Optional[List[str]] = None,
    goi_size: int = 10,
    mod_palette: str = "hsv",
    mod_num_font: Optional[Dict[str, Any]] = None,
    mod_num_left: float = 0.01,
    mod_num_right: float = 0.03,
    col_palette: str = "tab20",
    col_color_manual: Optional[Dict[str, Dict[str, str]]] = None,
    col_cat: Optional[Dict[str, Union[pd.Series, list, np.ndarray]]] = None,
    col_legend: Optional[list] = None,
    row_band: bool = True,
    row_band_font: Optional[Dict[str, Any]] = None,
    outdir: Optional[str] = ".",
    save_format: str = "png",
):
    """Draw the module clustermap. Pass a :class:`ModuleState` (from ``gen_mod``).

    Only ``rna_df`` and ``state`` are required; every visual knob has a default.
    Returns the seaborn ``ClusterGrid`` (``hm``); the right-hand band axis is
    attached as ``hm.ax_band`` when ``row_band=True``.
    """
    nm = canon_norm(norm_method)
    link = state.linkage
    labels = state.raw_labels

    rgb_colors, hex_colors = _module_colors(state, mod_palette)

    hm_args = {"robust": True, "yticklabels": False, "xticklabels": "auto",
               "cmap": "turbo", "col_cluster": True, **(hm_args or {})}
    for k in ("z_score", "standard_scale"):   # we normalize ourselves
        hm_args.pop(k, None)

    df = normalize_rows(rna_df, nm)

    # auto figure / dendrogram sizing when showing every sample name
    if hm_args.get("xticklabels") is True:
        n = df.shape[1]
        if "figsize" not in hm_args:
            base_h = 13.0
            base_w = min(base_h * 2, 0.1 * n) if 0.1 * n > base_h else base_h
            hm_args["figsize"] = (base_w, base_h)
        if "dendrogram_ratio" not in hm_args:
            row_ratio = hm_args["figsize"][1] / hm_args["figsize"][0] * 0.2
            hm_args["dendrogram_ratio"] = (row_ratio, 0.1)

    col_colors_df, col_color_map = _build_col_colors(
        df, col_cat, col_color_manual, col_palette)

    hm = sns.clustermap(df, row_linkage=link, row_colors=rgb_colors,
                        col_colors=col_colors_df, **hm_args)
    ax, fig = hm.ax_heatmap, hm.fig
    ax.set_zorder(2)
    ax.set_ylabel("")

    n_mods = len({i for i in set(labels) if i != 0})
    ax.set_title(title if title is not None
                 else f"n = {df.shape[0]}\nmodules = {n_mods}\nNorm: {nm}")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=xtickangle, ha="right")

    fs = _shrink_xticklabels(ax, xtickangle) if hm_args.get("xticklabels") is True else 10.0

    # genes of interest on the right
    order = hm.dendrogram_row.reordered_ind
    ordered_idx = df.index[order]
    if goi_list:
        pos = [ordered_idx.get_loc(g) for g in goi_list if g in ordered_idx]
        lab = [g for g in goi_list if g in ordered_idx]
        ax.set_yticks(pos)
        ax.set_yticklabels(lab, rotation=0, fontsize=goi_size)

    # module numbers (1-based, heatmap order)
    _draw_module_numbers(ax, np.asarray(labels)[order], mod_num_font,
                         mod_num_left, mod_num_right)

    # right-side color band
    hm.ax_band = None
    if row_band:
        hm.ax_band = _draw_row_band(fig, ax, np.asarray(rgb_colors)[order],
                                    row_band_font)

    if hm_args.get("yticklabels") is True:
        ax.set_yticks(np.arange(len(df.index)) + 0.5)
        ax.set_yticklabels(ordered_idx, rotation=0, fontsize=yticksize)

    if col_colors_df is not None and col_legend:
        _draw_col_legend(fig, col_color_map, col_legend)

    _place_cbar(hm, hm_args)

    dpi = _dpi_from_fontsize(fs) if hm_args.get("xticklabels") is True else 300
    hm._clusmap_dpi = dpi
    fig._clusmap_dpi = dpi

    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        fname = f"heatmap.{save_format.lower()}"
        kw = {"dpi": dpi} if save_format.lower() == "png" else {}
        hm.savefig(os.path.join(outdir, fname), bbox_inches="tight", **kw)
        print(f"Heatmap saved to {os.path.join(outdir, fname)}")
        state.to_modgene_df(colors=hex_colors).to_csv(
            os.path.join(outdir, "HM_ModGene.csv"), index=False)
        print(f"Heatmap module-gene mapping saved to {os.path.join(outdir, 'HM_ModGene.csv')}")

    return hm


def _draw_module_numbers(ax, labels_ordered, font, left, right):
    font = {"fontsize": 5, "fontweight": "bold", **(font or {})}
    change = np.r_[True, labels_ordered[1:] != labels_ordered[:-1]]
    starts = np.where(change)[0]
    ends = np.r_[starts[1:], len(labels_ordered)]
    new_id = 0
    for lab, s, e in zip(labels_ordered[starts], starts, ends):
        if lab == 0:
            continue
        new_id += 1
        y = s + (e - s) / 2.0
        x = -(left if new_id % 2 == 1 else right) * ax.get_xlim()[1]
        ax.text(x, y, str(new_id), va="center", ha="center", clip_on=False, **font)


def _draw_row_band(fig, ax, colors_ordered, font):
    font = {"band_width": 0.03, "alpha": 0.35, "band_gap": 0, **(font or {})}
    pos = ax.get_position()
    n = len(colors_ordered)
    ax_band = fig.add_axes([pos.x1 + font["band_gap"], pos.y0, font["band_width"], pos.height])
    ax_band.set_axis_off()
    rgba = np.array([mpl.colors.to_rgba(c, alpha=font["alpha"])
                     for c in colors_ordered]).reshape(-1, 1, 4)
    ax_band.imshow(rgba, aspect="auto", origin="upper")
    ax_band.set_ylim(n - 0.5, -0.5)
    ax_band.set_zorder(0)
    return ax_band


def _draw_col_legend(fig, col_color_map, col_legend):
    fig.subplots_adjust(right=0.82)
    y_anchor = 0.98
    for band, cmap in col_color_map.items():
        if band not in col_legend:
            continue
        handles = [Patch(color=c, label=cat) for cat, c in cmap.items()]
        leg = fig.legend(handles=handles, title=band, loc="upper left",
                         bbox_to_anchor=(0.83, y_anchor), borderaxespad=0.0,
                         frameon=False, fontsize=8, title_fontsize=9)
        fig.canvas.draw()
        bbox = leg.get_window_extent(renderer=fig.canvas.get_renderer())
        y_anchor = bbox.transformed(fig.transFigure.inverted()).y0 - 0.01


def _place_cbar(hm, hm_args):
    cbar_pos = hm_args.get("cbar_pos", "DEFAULT")
    if cbar_pos is None:
        if getattr(hm, "cax", None) is not None:
            hm.cax.set_visible(False)
    elif cbar_pos == "DEFAULT":
        hm.cax.set_position([0.02, 0.8, 0.05, 0.18])
    else:
        hm.cax.set_position(cbar_pos)


# --------------------------------------------------------------------------- #
# per-cluster sample statistics
# --------------------------------------------------------------------------- #
def cluster_sample_stats(
    rna_df: pd.DataFrame,
    hm_mod_gene_path: str,
    *,
    hm=None,
    module_col: str = "hm_mod",
    gene_col: str = "gene",
    stats: Optional[List[str]] = None,
    outdir: Optional[str] = ".",
    save_heatmap: bool = True,
) -> dict:
    """Compute per-cluster x sample statistics; returns {stat: DataFrame}."""
    stats = stats or ["mean", "std"]
    funcs = {
        "mean": np.mean, "std": np.std, "median": np.median,
        "iqr": lambda x: np.percentile(x, 75) - np.percentile(x, 25),
        "mad": lambda x: np.median(np.abs(x - np.median(x))),
        "cv": lambda x: np.std(x) / np.mean(x) if np.mean(x) else np.nan,
        "nonzero_frac": lambda x: np.count_nonzero(x) / len(x),
        "min": np.min, "max": np.max, "n": len,
    }
    selected = {k: v for k, v in funcs.items() if k in stats}

    if hm is not None:
        rna_df = rna_df[hm.data2d.columns.tolist()]

    map_df = pd.read_csv(hm_mod_gene_path, header=0)
    if module_col not in map_df.columns or gene_col not in map_df.columns:
        raise ValueError(f"Mapping file must contain columns: {module_col}, {gene_col}")
    gene_to_mod = pd.Series(map_df[module_col].values,
                            index=map_df[gene_col].astype(str)).to_dict()
    labels = pd.Series([gene_to_mod.get(g, np.nan) for g in rna_df.index],
                       index=rna_df.index, name=module_col)
    frames = {stat: rna_df.groupby(labels).agg(func) for stat, func in selected.items()}

    result = dict(frames)   # always defined (fixes prior NameError when outdir=None)
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        for stat, df in frames.items():
            df.to_csv(os.path.join(outdir, f"cluster_sample_{stat}.tsv"), sep="\t")
            print(f"Saved {stat} to {os.path.join(outdir, f'cluster_sample_{stat}.tsv')}")
            if save_heatmap:
                _stat_heatmap(df, stat, outdir)
    return result


def _stat_heatmap(df, stat, outdir):
    n = df.shape[1]
    base_h = 11.0
    base_w = min(base_h * 2, 0.1 * n) if 0.1 * n > base_h else base_h
    fig, ax = plt.subplots(figsize=(base_w, base_h))
    z = df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, np.nan), axis=0)
    sns.heatmap(z, cmap="turbo", robust=True, ax=ax, xticklabels=True)
    ax.set_title(f"{stat} expression per cluster x sample")
    ax.set_ylabel("Cluster")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    _shrink_xticklabels(ax, 45)
    plt.tight_layout()
    out = os.path.join(outdir, f"cluster_sample_{stat}_heatmap.png")
    plt.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Saved heatmap for {stat} to {out}")
