"""Module annotation: pseudo-bulk heatmap, marker swarm, celltype enrichment, GO."""
from __future__ import annotations

import os
import traceback
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats_mod
from statsmodels.stats.multitest import multipletests

from .util import canon_norm, normalize_rows


# --------------------------------------------------------------------------- #
# pseudo-bulk from single cell
# --------------------------------------------------------------------------- #
def compute_pseudo_bulk(h5ad_file, celltype_key, *, layer=None, outdir=".",
                        log_base=None, plus_constant=0.1, chunk_size=None):
    """Mean expression per cell type (genes x cell types).

    If ``chunk_size`` is given the file is read in backed mode in chunks
    (memory-efficient for large atlases); otherwise it is loaded fully.
    """
    import anndata as ad
    from scipy.sparse import issparse, csr_matrix

    if chunk_size:
        adata = ad.read_h5ad(h5ad_file, backed="r")
        celltypes = sorted(adata.obs[celltype_key].unique().tolist())
        ct_to_idx = {ct: i for i, ct in enumerate(celltypes)}
        pb = np.zeros((adata.n_vars, len(celltypes)))
        for start in range(0, adata.n_obs, chunk_size):
            end = min(start + chunk_size, adata.n_obs)
            chunk = adata[start:end]
            data = csr_matrix(chunk.layers[layer] if layer is not None else chunk.X)
            col_idx = [ct_to_idx[ct] for ct in chunk.obs[celltype_key]]
            ind = csr_matrix((np.ones(end - start), (np.arange(end - start), col_idx)),
                             shape=(end - start, len(celltypes)))
            pb += (data.T @ ind).toarray()
            print(f"Processed {end}/{adata.n_obs} cells", end="\r")
        counts = adata.obs[celltype_key].value_counts()[celltypes].values
        means = pb / counts
        var_names = adata.var_names
    else:
        adata = ad.read_h5ad(h5ad_file)
        expr = adata.layers[layer] if layer is not None else adata.X
        if not issparse(expr):
            expr = csr_matrix(expr)
        codes, celltypes = pd.factorize(adata.obs[celltype_key], sort=True)
        ind = csr_matrix((np.ones_like(codes), (np.arange(len(codes)), codes)),
                         shape=(len(codes), len(celltypes)))
        means = expr.T @ ind
        means = (means.toarray() if issparse(means) else means) / \
            np.array(ind.sum(axis=0)).ravel()
        var_names = adata.var_names

    df = pd.DataFrame(means, index=var_names, columns=list(celltypes))
    if log_base is not None:
        df = np.log(df + plus_constant) / np.log(log_base)
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        df.to_csv(os.path.join(outdir, "pseudo_bulk.tsv"), sep="\t")
        print(f"\nPseudo-bulk saved to {os.path.join(outdir, 'pseudo_bulk.tsv')}")
    return df


def pseudo_bulk_hm(rna_df, hm, *, pseudo_bulk_df=None, pseudo_bulk_path=None,
                   pseudo_bulk_sep="\t", log_base=None, plus_constant=0.1,
                   norm_method="z_score", xtickangle=0, cbar_pos=None,
                   left_gap=0.08, col_width=0.1, pseudo_hm_args=None,
                   title="Pseudo Bulk Heatmap", title_size=10,
                   outdir=".", save_format="png"):
    """Append a pseudo-bulk heatmap to an existing bulk clustermap.

    ``log_base`` is applied consistently whether the data comes from a file or a
    DataFrame (set to ``None`` if already log-scaled). Missing genes are 0-filled.
    """
    fig = hm.fig
    heat_pos = hm.ax_heatmap.get_position()
    ordered = list(rna_df.index[hm.dendrogram_row.reordered_ind])

    if pseudo_bulk_df is None and pseudo_bulk_path is not None:
        pseudo_bulk_df = pd.read_csv(pseudo_bulk_path, sep=pseudo_bulk_sep,
                                     index_col=0, header=0)
    pseudo_bulk_df = pseudo_bulk_df.copy()
    if log_base is not None:
        pseudo_bulk_df = np.log(pseudo_bulk_df + plus_constant) / np.log(log_base)

    missing = [g for g in ordered if g not in pseudo_bulk_df.index]
    for g in missing:
        pseudo_bulk_df.loc[g] = 0.0
    df = normalize_rows(pseudo_bulk_df.loc[ordered], norm_method)
    nm = canon_norm(norm_method)

    width = col_width * df.shape[1]
    ax = fig.add_axes([heat_pos.x1 + left_gap, heat_pos.y0, width, heat_pos.height])
    args = {"robust": True, "cmap": "turbo", "cbar": True, **(pseudo_hm_args or {})}
    if cbar_pos is not None:
        args["cbar_ax"] = fig.add_axes(cbar_pos)
    sns.heatmap(df, ax=ax, **args)
    ax.set_xlabel(""); ax.set_ylabel(""); ax.set_yticks([])
    ax.set_xticklabels(df.columns, rotation=xtickangle, ha="right")
    ax.set_title(f"{title}\nNorm: {nm}", fontsize=title_size)

    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        fname = f"pb_heatmap.{save_format.lower()}"
        kw = {"dpi": getattr(fig, "_clusmap_dpi", 300)} if save_format.lower() == "png" else {}
        hm.savefig(os.path.join(outdir, fname), bbox_inches="tight", **kw)
        print(f"Heatmap with pseudo-bulk saved to {os.path.join(outdir, fname)}")
    return hm, ax


# --------------------------------------------------------------------------- #
# marker swarm plot
# --------------------------------------------------------------------------- #
def sc_marker_hm(rna_df, hm, marker_file, *, sep="\t", celltype_col="celltype",
                 gene_col="gene", celltype=None, left_gap=0.08, col_width=0.1,
                 title="Marker Gene Swarm Plot", spine=True, axhspan_args=None,
                 swarm_args=None, marker_show=False, outdir=".", save_format="png"):
    """Append a marker-gene swarm plot (PanglaoDB-style) beside the heatmap."""
    fig = hm.fig
    heat_pos = hm.ax_heatmap.get_position()
    marker_df = pd.read_csv(marker_file, sep=sep, header=0)

    order = hm.dendrogram_row.reordered_ind
    ordered = list(rna_df.index[order])
    ordered_rev = ordered[::-1]
    gene_list = rna_df.index.astype(str).tolist()
    order_dict = {g: i for i, g in enumerate(ordered_rev)}

    if celltype is None:
        celltype = marker_df[celltype_col].unique().tolist()
    else:
        marker_df = marker_df[marker_df[celltype_col].isin(celltype)]

    matrix = pd.DataFrame(np.nan, columns=celltype, index=gene_list)
    low_map = {str(g).lower(): g for g in gene_list}
    for cell in celltype:
        markers = {g.lower() for g in
                   marker_df[marker_df[celltype_col] == cell][gene_col].astype(str)}
        hits = [low_map[g] for g in low_map if g in markers]
        matrix.loc[hits, cell] = [order_dict[h] for h in hits]

    row_colors_ordered = np.asarray(hm.row_colors)[order]
    width = col_width * matrix.shape[1]
    ax = fig.add_axes([heat_pos.x1 + left_gap, heat_pos.y0, width, heat_pos.height])
    axhspan_args = {"alpha": 0.1, **(axhspan_args or {})}
    for i in range(len(gene_list)):
        ax.axhspan(len(gene_list) - 0.5 - i, len(gene_list) - 0.5 - (i + 1),
                   facecolor=row_colors_ordered[i], zorder=0, **axhspan_args)

    long_df = matrix.reset_index().melt(id_vars="index", var_name="celltype",
                                        value_name="marker_pos")
    swarm = sns.swarmplot(data=long_df, x="celltype", y="marker_pos", zorder=1,
                          ax=ax, **(swarm_args or {}))
    swarm.set_title(title)
    ax.set_ylim(-0.5, len(gene_list) - 0.5)
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    if marker_show:
        mk = long_df.dropna(subset=["marker_pos"])
        ypos = sorted(mk["marker_pos"].unique())
        ax.set_yticks(ypos)
        ax.set_yticklabels([mk[mk["marker_pos"] == p]["index"].iloc[0] for p in ypos],
                           fontsize="small")
        ax.tick_params(axis="y", left=True)
    else:
        ax.set_yticks([])
    if not spine:
        sns.despine(ax=ax)
        for s in ax.spines.values():
            s.set_visible(False)

    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        fname = f"heatmap_swarm.{save_format.lower()}"
        kw = {"dpi": getattr(fig, "_clusmap_dpi", 300)} if save_format.lower() == "png" else {}
        hm.savefig(os.path.join(outdir, fname), bbox_inches="tight", **kw)
        print(f"Heatmap with swarm plot saved to {os.path.join(outdir, fname)}")
    return hm, swarm


# --------------------------------------------------------------------------- #
# celltype enrichment selection
# --------------------------------------------------------------------------- #
def celltype_selection(marker_path, module_path, *, marker_sep="\t", module_sep=",",
                       marker_col_celltype="cell type",
                       marker_col_gene="official gene symbol",
                       module_col_gene="gene", module_col_module="hm_mod",
                       p_val_threshold=0.05, max_celltypes_per_module=5,
                       outdir=".", hm_cmap="Reds", heatmap_top_n_celltypes=30):
    """Fisher-exact enrichment of cell-type markers in each module (FDR-corrected).

    Returns the list of selected (most correlated) cell types.
    """
    marker_df = pd.read_csv(marker_path, sep=marker_sep, header=0)
    marker_df.columns = marker_df.columns.str.strip()
    module_df = pd.read_csv(module_path, sep=module_sep, header=0)

    universe = set(module_df[module_col_gene].dropna().astype(str).str.upper())
    modules = module_df[module_col_module].dropna().unique()
    mod_genes = {m: set(module_df[module_df[module_col_module] == m][module_col_gene]
                        .dropna().astype(str).str.upper()) for m in modules}

    celltypes = marker_df[marker_col_celltype].dropna().unique()
    ct_markers = {ct: set(marker_df[marker_df[marker_col_celltype] == ct][marker_col_gene]
                          .dropna().astype(str).str.upper()) & universe for ct in celltypes}
    valid_cts = [ct for ct in celltypes if ct_markers[ct]]

    p_tab = pd.DataFrame(index=modules, columns=valid_cts, dtype=float)
    N = len(universe)
    for m, mg in mod_genes.items():
        for ct in valid_cts:
            cg = ct_markers[ct]
            a = len(mg & cg)
            b, c = len(mg) - a, len(cg) - a
            d = N - a - b - c
            _, p = stats_mod.fisher_exact([[a, b], [c, d]], alternative="greater")
            p_tab.loc[m, ct] = p

    _, fdr, _, _ = multipletests(p_tab.values.flatten(), method="fdr_bh")
    p_tab.iloc[:, :] = fdr.reshape(p_tab.shape)

    selected = set()
    for m in p_tab.index:
        sig = p_tab.loc[m].sort_values()
        selected.update(sig[sig < p_val_threshold].head(max_celltypes_per_module).index)
    selected = list(selected)

    if outdir is not None:
        outdir = os.path.join(outdir, "celltype_selection")
        os.makedirs(outdir, exist_ok=True)
        p_tab.to_csv(os.path.join(outdir, "p_val_table.tsv"), sep="\t")
        with open(os.path.join(outdir, "selected_celltypes.txt"), "w") as fh:
            fh.write("\n".join(selected) + "\n")
        rows = []
        for m in p_tab.index:
            row = {"module": m}
            for i, (ct, p) in enumerate(p_tab.loc[m].sort_values()
                                        .head(max_celltypes_per_module).items()):
                row[f"Rank_{i+1}"] = f"{ct} ({p:.2e})"
            rows.append(row)
        pd.DataFrame(rows).to_csv(os.path.join(outdir, "module_top_celltypes.tsv"),
                                  sep="\t", index=False)
        if selected:
            top = p_tab[selected].min(axis=0).nsmallest(heatmap_top_n_celltypes).index.tolist()
            data = -np.log10(p_tab[top].astype(float) + 1e-300)
            ylabels = [f"{m} ({p_tab.loc[m].idxmin()})" for m in data.index]
            plt.figure(figsize=(max(10, len(top) * 0.3), max(8, len(ylabels) * 0.2)))
            ax = sns.heatmap(data, cmap=hm_cmap, yticklabels=ylabels)
            ax.set_yticklabels(ax.get_yticklabels(),
                               size=max(4, min(10, 800 / max(1, len(ylabels)))))
            plt.title(f"-log10(P) of Module vs Top {len(top)} Cell Types")
            plt.xlabel("Selected Cell Types"); plt.ylabel("Modules (top cell type)")
            plt.xticks(rotation=45, ha="right",
                       size=max(4, min(10, 800 / max(1, len(top)))))
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, "module_celltype_heatmap.pdf"))
            plt.close()
    return selected


# --------------------------------------------------------------------------- #
# GO enrichment
# --------------------------------------------------------------------------- #
def mod_GO(modgene_path, *, bg=None, module_col="hm_mod", gene_col="gene",
           organism="Mouse", library=None, GO_category=("BP", "MF"),
           outdir=".", mod=None, barplot=None, dotplot=None,
           barplot_args=None, dotplot_args=None, enrichr_args=None):
    """Over-representation (Enrichr) GO analysis per module.

    ``module_col`` defaults to ``hm_mod`` to match ``HM_ModGene.csv``.
    """
    import gseapy as gp

    barplot = (bg is not None) if barplot is None else barplot
    dotplot = (bg is None) if dotplot is None else dotplot
    barplot_args, dotplot_args = barplot_args or {}, dotplot_args or {}
    enrichr_args = enrichr_args or {}
    if bg is not None and enrichr_args.get("background") is None:
        enrichr_args["background"] = bg.index.astype(str).tolist()

    mod_gene_df = pd.read_csv(modgene_path, header=0)
    os.makedirs(outdir, exist_ok=True)
    if mod is not None:
        mod_gene_df = mod_gene_df[mod_gene_df[module_col].astype(int).isin(mod)]

    def run(go_cat, lib):
        print(f"Using {go_cat} library: {lib}")
        for m, sub in mod_gene_df.groupby(module_col, sort=True):
            genes = sub[gene_col].astype(str).tolist()
            try:
                enr = gp.enrichr(gene_list=genes, gene_sets=lib, organism=organism,
                                 outdir=None, no_plot=True, **enrichr_args)
            except ValueError as e:
                if "No enrich terms" in str(e):
                    print(f"[WARN] no terms: module={m}, {go_cat}")
                    continue
                traceback.print_exc(); continue
            sub_out = os.path.join(outdir, go_cat)
            os.makedirs(sub_out, exist_ok=True)
            enr.results.to_csv(os.path.join(sub_out, f"module_{m}_GO.csv"), index=False)
            print(f"Module {m} {go_cat}: {enr.res2d.shape[0]} terms")
            if barplot:
                try:
                    gp.plot.barplot(enr.res2d, title=f"Module {m}\n{lib}",
                                    ofname=os.path.join(sub_out, f"module_{m}_barplot.png"),
                                    **barplot_args)
                except Exception:
                    traceback.print_exc()
            if dotplot:
                try:
                    gp.plot.dotplot(enr.res2d, title=f"Module {m}\n{lib}",
                                    ofname=os.path.join(sub_out, f"module_{m}_dotplot.png"),
                                    **{k: v for k, v in dotplot_args.items() if k != "DynamicSize"})
                except Exception:
                    traceback.print_exc()

    if library is None:
        libs = gp.get_library_name(organism)
        cat_prefix = {"BP": "GO_Biological_Process_", "MF": "GO_Molecular_Function_",
                      "CC": "GO_Cellular_Component_"}
        for cat in GO_category:
            pref = cat_prefix[cat]
            picks = [l for l in libs if l.startswith(pref) and l[len(pref):].isdigit()]
            if picks:
                run(cat, max(picks))
    else:
        for lib in library:
            tag = ("BP" if "Biological_Process" in lib else
                   "MF" if "Molecular_Function" in lib else
                   "CC" if "Cellular_Component" in lib else
                   "KEGG" if "KEGG" in lib else "Custom")
            run(f"{tag}_{lib}", lib)
