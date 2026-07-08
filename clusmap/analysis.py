"""Higher-level module analyses (WGCNA-style), all built on ModuleState.

* ``module_eigengenes`` – the first principal component of each module
  (a "representative expression profile", modules x samples).
* ``hub_genes`` – intramodular connectivity (kME = correlation of each gene with
  its own module eigengene); the top-kME genes are the module's hubs.
* ``module_trait_correlation`` – correlate eigengenes with sample traits.
* ``project_modules`` – assign a new dataset's genes to existing modules.
* ``module_report`` – a single self-contained HTML bringing module size, hub
  genes, GO, cell types and motifs together.
"""
from __future__ import annotations

import glob
import html
import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

from .state import ModuleState


# --------------------------------------------------------------------------- #
# eigengenes
# --------------------------------------------------------------------------- #
def _zscore_rows(mat: np.ndarray) -> np.ndarray:
    mu = mat.mean(axis=1, keepdims=True)
    sd = mat.std(axis=1, ddof=0, keepdims=True)
    sd[sd == 0] = 1.0
    return (mat - mu) / sd


def module_eigengenes(rna_df: pd.DataFrame, state: ModuleState, *,
                      outdir: Optional[str] = ".") -> pd.DataFrame:
    """First-PC eigengene per module (modules x samples).

    Each module's genes are z-scored across samples, then the first right
    singular vector is taken as the eigengene and sign-aligned so it correlates
    positively with the module's mean profile (WGCNA convention).
    """
    rna_df = rna_df.loc[state.genes]
    samples = list(rna_df.columns)
    rows: Dict[int, np.ndarray] = {}
    for m in sorted({int(x) for x in state.hm_labels if x != 0}):
        genes = state.module_genes(m)
        if not genes:
            continue
        Xz = _zscore_rows(rna_df.loc[genes].values.astype(float))
        # first right singular vector (samples direction)
        _, _, Vt = np.linalg.svd(Xz, full_matrices=False)
        eig = Vt[0]
        if np.corrcoef(eig, Xz.mean(axis=0))[0, 1] < 0:
            eig = -eig
        rows[m] = eig
    eig_df = pd.DataFrame(rows, index=samples).T
    eig_df.index.name = "module"
    if outdir is not None and len(eig_df):
        os.makedirs(outdir, exist_ok=True)
        eig_df.to_csv(os.path.join(outdir, "module_eigengenes.tsv"), sep="\t")
        print(f"Eigengenes saved to {os.path.join(outdir, 'module_eigengenes.tsv')}")
    return eig_df


# --------------------------------------------------------------------------- #
# hub genes (kME)
# --------------------------------------------------------------------------- #
def hub_genes(rna_df: pd.DataFrame, state: ModuleState, *,
              eigengenes: Optional[pd.DataFrame] = None, top_n: int = 20,
              outdir: Optional[str] = ".") -> pd.DataFrame:
    """Rank genes by intramodular connectivity (kME).

    kME(gene) = Pearson correlation of the gene's expression with its module
    eigengene. Returns a tidy table (gene, module, kME) sorted by module then
    descending kME; the top ``top_n`` per module are the hub genes.
    """
    if eigengenes is None:
        eigengenes = module_eigengenes(rna_df, state, outdir=None)
    rna_df = rna_df.loc[state.genes]
    recs = []
    for m in eigengenes.index:
        eig = eigengenes.loc[m].values.astype(float)
        for g in state.module_genes(int(m)):
            x = rna_df.loc[g].values.astype(float)
            if x.std() == 0:
                kme = np.nan
            else:
                kme = float(np.corrcoef(x, eig)[0, 1])
            recs.append((g, int(m), kme))
    df = pd.DataFrame(recs, columns=["gene", "module", "kME"])
    df = df.sort_values(["module", "kME"], ascending=[True, False]).reset_index(drop=True)
    df["hub_rank"] = df.groupby("module").cumcount() + 1
    df["is_hub"] = df["hub_rank"] <= top_n
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        df.to_csv(os.path.join(outdir, "hub_genes.tsv"), sep="\t", index=False)
        print(f"Hub genes saved to {os.path.join(outdir, 'hub_genes.tsv')}")
    return df


def top_hubs(hub_df: pd.DataFrame, top_n: int = 10) -> Dict[int, List[str]]:
    """Convenience: {module: [top hub genes]} from a hub_genes table."""
    return {int(m): sub.head(top_n)["gene"].tolist()
            for m, sub in hub_df[hub_df["kME"].notna()].groupby("module")}


# --------------------------------------------------------------------------- #
# module-trait correlation
# --------------------------------------------------------------------------- #
def module_trait_correlation(eigengenes: pd.DataFrame, traits: pd.DataFrame, *,
                             outdir: Optional[str] = ".", save_heatmap: bool = True):
    """Correlate module eigengenes with sample traits (FDR-corrected).

    ``traits`` is a samples x traits DataFrame; non-numeric columns are
    one-hot encoded. Returns (correlation, fdr_pvalue) DataFrames (modules x traits).
    """
    traits = traits.reindex(eigengenes.columns)        # align to eigengene samples
    num = traits.apply(pd.to_numeric, errors="coerce")
    non_numeric = [c for c in traits.columns if num[c].isna().all()]
    if non_numeric:
        dummies = pd.get_dummies(traits[non_numeric].astype("object"), prefix_sep="=")
        num = pd.concat([num.drop(columns=non_numeric), dummies], axis=1)
    num = num.dropna(axis=1, how="all")

    cor = pd.DataFrame(index=eigengenes.index, columns=num.columns, dtype=float)
    pval = cor.copy()
    for m in eigengenes.index:
        e = eigengenes.loc[m].values.astype(float)
        for t in num.columns:
            v = num[t].values.astype(float)
            mask = ~np.isnan(v)
            if mask.sum() < 3 or np.nanstd(v[mask]) == 0:
                cor.loc[m, t], pval.loc[m, t] = np.nan, np.nan
            else:
                r, p = stats.pearsonr(e[mask], v[mask])
                cor.loc[m, t], pval.loc[m, t] = r, p

    flat = pval.values.flatten().astype(float)
    ok = ~np.isnan(flat)
    if ok.any():
        flat[ok] = multipletests(flat[ok], method="fdr_bh")[1]
    fdr = pd.DataFrame(flat.reshape(pval.shape), index=pval.index, columns=pval.columns)

    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        cor.to_csv(os.path.join(outdir, "module_trait_cor.tsv"), sep="\t")
        fdr.to_csv(os.path.join(outdir, "module_trait_fdr.tsv"), sep="\t")
        if save_heatmap and cor.notna().any().any():
            _trait_heatmap(cor, fdr, os.path.join(outdir, "module_trait_heatmap.png"))
        print(f"Module-trait correlation saved to {outdir}/")
    return cor, fdr


def _trait_heatmap(cor, fdr, path):
    import matplotlib.pyplot as plt
    import seaborn as sns
    annot = cor.round(2).astype(str) + fdr.map(
        lambda p: "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "")
    fig, ax = plt.subplots(figsize=(max(6, cor.shape[1] * 0.9), max(4, cor.shape[0] * 0.4)))
    sns.heatmap(cor.astype(float), cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                annot=annot.values, fmt="", annot_kws={"size": 7}, ax=ax,
                cbar_kws={"label": "Pearson r"})
    ax.set_title("Module–trait correlation (*FDR<.05 **<.01 ***<.001)")
    ax.set_ylabel("Module"); ax.set_xlabel("Trait")
    plt.tight_layout(); plt.savefig(path, dpi=300); plt.close(fig)


# --------------------------------------------------------------------------- #
# module preservation
# --------------------------------------------------------------------------- #
def _module_stats(ref: np.ndarray, test: np.ndarray):
    """Density + connectivity stats for one gene set, comparing ref vs test.

    ref/test are genes x samples arrays for the *same* genes in two datasets.
    Returns (meanCor_test, meanAdj_test, cor_kIM, cor_kME).
    """
    cr = np.corrcoef(ref)
    ct = np.corrcoef(test)
    n = ct.shape[0]
    off = ~np.eye(n, dtype=bool)
    mean_cor = float(np.nanmean(ct[off]))
    mean_adj = float(np.nanmean(np.abs(ct[off])))
    # intramodular connectivity kIM = row-sum of |adjacency| (excl. self)
    kim_ref = np.nansum(np.abs(cr) - np.eye(n), axis=1)
    kim_test = np.nansum(np.abs(ct) - np.eye(n), axis=1)
    cor_kim = float(np.corrcoef(kim_ref, kim_test)[0, 1]) if np.std(kim_ref) and np.std(kim_test) else np.nan
    # kME = correlation of each gene with the set's first-PC eigengene
    def _kme(mat):
        _, _, Vt = np.linalg.svd(_zscore_rows(mat), full_matrices=False)
        eig = Vt[0]
        return np.array([np.corrcoef(mat[i], eig)[0, 1] for i in range(mat.shape[0])])
    kme_ref, kme_test = _kme(ref), _kme(test)
    cor_kme = float(np.corrcoef(kme_ref, kme_test)[0, 1]) if np.std(kme_ref) and np.std(kme_test) else np.nan
    return mean_cor, mean_adj, cor_kim, cor_kme


def module_preservation(ref_rna: pd.DataFrame, state: ModuleState, test_rna: pd.DataFrame, *,
                        n_perm: int = 200, min_genes: int = 5, seed: int = 0,
                        outdir: Optional[str] = ".") -> pd.DataFrame:
    """Permutation-based module preservation between a reference and a test dataset.

    For each module it measures, in the *test* data, whether the module's genes
    stay densely co-expressed (density) and keep the same hub structure
    (connectivity), then z-scores each statistic against random gene sets of the
    same size (``n_perm`` permutations). ``Zsummary`` summarises both:

        Z >= 10  strongly preserved · 2 <= Z < 10  moderate · Z < 2  not preserved

    The two datasets must share genes and have comparable sample axes
    (replicate, independent cohort, or a split of one dataset). This follows the
    Langfelder et al. (2011) preservation framework in spirit; it is a compact
    re-implementation, not a byte-for-byte WGCNA port.
    """
    rng = np.random.default_rng(seed)
    shared = [g for g in state.genes if g in test_rna.index and g in ref_rna.index]
    ref_rna = ref_rna.loc[shared]
    test_rna = test_rna.loc[shared]
    universe = list(shared)

    recs = []
    for m in sorted({int(x) for x in state.hm_labels if x != 0}):
        genes = [g for g in state.module_genes(m) if g in test_rna.index]
        if len(genes) < min_genes:
            continue
        obs = _module_stats(ref_rna.loc[genes].values, test_rna.loc[genes].values)
        # null: random gene sets of the same size from the shared universe
        null = np.full((n_perm, 4), np.nan)
        n = len(genes)
        for i in range(n_perm):
            pick = rng.choice(len(universe), size=n, replace=False)
            sub = [universe[j] for j in pick]
            null[i] = _module_stats(ref_rna.loc[sub].values, test_rna.loc[sub].values)

        def _z(idx):
            mu, sd = np.nanmean(null[:, idx]), np.nanstd(null[:, idx])
            return (obs[idx] - mu) / sd if sd > 0 else np.nan

        z_density = np.nanmean([_z(0), _z(1)])          # meanCor, meanAdj
        z_conn = np.nanmean([_z(2), _z(3)])             # cor.kIM, cor.kME
        zsummary = np.nanmean([z_density, z_conn])
        quality = ("strong" if zsummary >= 10 else
                   "moderate" if zsummary >= 2 else "none")
        recs.append({"module": m, "n_genes": n, "meanCor_test": round(obs[0], 3),
                     "cor_kIM": round(obs[2], 3) if obs[2] == obs[2] else np.nan,
                     "cor_kME": round(obs[3], 3) if obs[3] == obs[3] else np.nan,
                     "Zdensity": round(z_density, 2), "Zconnectivity": round(z_conn, 2),
                     "Zsummary": round(zsummary, 2), "preservation": quality})
    df = pd.DataFrame(recs)
    if outdir is not None and len(df):
        os.makedirs(outdir, exist_ok=True)
        df.to_csv(os.path.join(outdir, "module_preservation.tsv"), sep="\t", index=False)
        _preservation_plot(df, os.path.join(outdir, "module_preservation.png"))
        print(f"Module preservation saved to {os.path.join(outdir, 'module_preservation.tsv')}")
    return df


def _preservation_plot(df, path):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(max(5, len(df) * 0.4), 4))
    colors = ["#16a34a" if z >= 10 else "#f59e0b" if z >= 2 else "#dc2626"
              for z in df["Zsummary"]]
    ax.bar(df["module"].astype(str), df["Zsummary"], color=colors)
    ax.axhline(2, ls="--", c="#999", lw=1); ax.axhline(10, ls="--", c="#999", lw=1)
    ax.set_xlabel("Module"); ax.set_ylabel("Zsummary")
    ax.set_title("Module preservation (green≥10 strong, amber≥2 moderate, red<2 none)")
    plt.tight_layout(); plt.savefig(path, dpi=300); plt.close(fig)


# --------------------------------------------------------------------------- #
# cross-dataset projection
# --------------------------------------------------------------------------- #
def project_modules(new_rna_df: pd.DataFrame, eigengenes: pd.DataFrame, *,
                    min_cor: float = 0.5, outdir: Optional[str] = ".") -> pd.DataFrame:
    """Assign each gene in a *new* dataset to the best-matching reference module.

    Requires the new dataset to share sample columns with the reference
    eigengenes. Each gene is assigned to the module whose eigengene it correlates
    with most strongly (>= ``min_cor``), else module 0 (unassigned).
    """
    shared = [c for c in eigengenes.columns if c in new_rna_df.columns]
    if len(shared) < 3:
        raise ValueError("New dataset must share >=3 sample columns with the reference "
                         f"eigengenes (shared: {len(shared)}).")
    E = eigengenes[shared].values.astype(float)
    recs = []
    for g in new_rna_df.index:
        x = new_rna_df.loc[g, shared].values.astype(float)
        if np.nanstd(x) == 0:
            recs.append((g, 0, np.nan)); continue
        cors = [np.corrcoef(x, E[i])[0, 1] for i in range(E.shape[0])]
        best = int(np.nanargmax(cors)); r = cors[best]
        recs.append((g, int(eigengenes.index[best]) if r >= min_cor else 0, float(r)))
    df = pd.DataFrame(recs, columns=["gene", "assigned_module", "best_cor"])
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        df.to_csv(os.path.join(outdir, "projected_modules.tsv"), sep="\t", index=False)
        print(f"Projection saved to {os.path.join(outdir, 'projected_modules.tsv')}")
    return df


# --------------------------------------------------------------------------- #
# HTML report
# --------------------------------------------------------------------------- #
def _sparkline(values, w=160, h=28) -> str:
    v = np.asarray(values, float)
    if v.size == 0 or np.nanstd(v) == 0:
        return ""
    lo, hi = np.nanmin(v), np.nanmax(v)
    xs = np.linspace(2, w - 2, v.size)
    ys = h - 2 - (v - lo) / (hi - lo) * (h - 4)
    pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in zip(xs, ys))
    zero = h - 2 - (0 - lo) / (hi - lo) * (h - 4) if lo < 0 < hi else None
    base = f'<line x1="2" y1="{zero:.1f}" x2="{w-2}" y2="{zero:.1f}" stroke="#ccc"/>' if zero else ""
    return (f'<svg width="{w}" height="{h}" style="vertical-align:middle">{base}'
            f'<polyline points="{pts}" fill="none" stroke="#2563eb" stroke-width="1.5"/></svg>')


def _read_go_top(outdir, module, n=5):
    for cat in ("BP", "MF", "CC"):
        f = os.path.join(outdir, cat, f"module_{module}_GO.csv")
        if os.path.exists(f):
            df = pd.read_csv(f)
            term = next((c for c in ("Term", "term") if c in df.columns), None)
            if term:
                return cat, df[term].head(n).tolist()
    return None, []


def module_report(outdir: str, rna_df: pd.DataFrame, state: ModuleState, *,
                  eigengenes: Optional[pd.DataFrame] = None,
                  hub_df: Optional[pd.DataFrame] = None,
                  report_name: str = "module_report.html") -> str:
    """Bundle per-module size, hub genes, eigengene sparkline, GO, cell types and
    motifs (whatever exists in ``outdir``) into one self-contained HTML file."""
    if eigengenes is None:
        eigengenes = module_eigengenes(rna_df, state, outdir=None)
    if hub_df is None:
        hub_df = hub_genes(rna_df, state, eigengenes=eigengenes, outdir=None)
    hubs = top_hubs(hub_df, top_n=10)

    # optional annotations on disk
    ct_map = {}
    ct_file = os.path.join(outdir, "celltype_selection", "module_top_celltypes.tsv")
    if os.path.exists(ct_file):
        ct = pd.read_csv(ct_file, sep="\t")
        rank_cols = [c for c in ct.columns if c.startswith("Rank_")]
        for _, r in ct.iterrows():
            ct_map[int(r["module"])] = ", ".join(str(r[c]) for c in rank_cols
                                                 if pd.notna(r[c]))[:160]
    motif_map = {}
    for mf in glob.glob(os.path.join(outdir, "**", "module_motif_summary.tsv"), recursive=True):
        try:
            for _, r in pd.read_csv(mf, sep="\t", header=None).iterrows():
                motif_map.setdefault(int(r[0]), ", ".join(
                    str(x) for x in r[1:6] if pd.notna(x)))
        except Exception:
            pass
    pres_map = {}
    pres_file = os.path.join(outdir, "module_preservation.tsv")
    if os.path.exists(pres_file):
        for _, r in pd.read_csv(pres_file, sep="\t").iterrows():
            pres_map[int(r["module"])] = (float(r["Zsummary"]), str(r["preservation"]))

    sizes = {int(m): len(state.module_genes(int(m))) for m in eigengenes.index}
    rows_html = []
    for m in eigengenes.index:
        m = int(m)
        cat, go = _read_go_top(outdir, m)
        go_txt = (f"<b>GO ({cat}):</b> " + "; ".join(map(html.escape, go))) if go else ""
        ct_txt = f"<b>Cell types:</b> {html.escape(ct_map[m])}" if m in ct_map else ""
        mo_txt = f"<b>Motifs:</b> {html.escape(motif_map[m])}" if m in motif_map else ""
        extras = "<br>".join(t for t in (go_txt, ct_txt, mo_txt) if t)
        if m in pres_map:
            z, q = pres_map[m]
            color = {"strong": "#16a34a", "moderate": "#f59e0b"}.get(q, "#dc2626")
            pres_cell = f'<span style="color:{color};font-weight:600">Z={z:.1f} ({q})</span>'
        else:
            pres_cell = "—"
        rows_html.append(f"""
        <tr>
          <td class="m">{m}</td>
          <td>{sizes.get(m, 0)}</td>
          <td>{pres_cell}</td>
          <td>{_sparkline(eigengenes.loc[m].values)}</td>
          <td class="hub">{html.escape(', '.join(hubs.get(m, [])))}</td>
          <td class="ann">{extras}</td>
        </tr>""")

    doc = f"""<!doctype html><html><head><meta charset="utf-8">
<title>clusmap module report</title>
<style>
 body{{font:14px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;margin:2rem;color:#222}}
 h1{{margin:0 0 .2rem}} .sub{{color:#666;margin-bottom:1rem}}
 table{{border-collapse:collapse;width:100%}}
 th,td{{border-bottom:1px solid #eee;padding:8px 10px;text-align:left;vertical-align:top}}
 th{{background:#f8fafc;position:sticky;top:0}}
 td.m{{font-weight:700;color:#2563eb}} td.hub{{font-family:monospace;font-size:12px;max-width:280px}}
 td.ann{{font-size:12px;color:#444;max-width:360px}}
</style></head><body>
<h1>clusmap module report</h1>
<div class="sub">{rna_df.shape[0]} genes · {len(eigengenes)} modules · edits: {html.escape(str(state.history)) or 'none'}</div>
<table><thead><tr><th>Module</th><th>Size</th><th>Preservation</th><th>Eigengene</th>
<th>Hub genes (top kME)</th><th>Annotations</th></tr></thead>
<tbody>{''.join(rows_html)}</tbody></table>
</body></html>"""

    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, report_name)
    with open(path, "w") as fh:
        fh.write(doc)
    print(f"Module report written to {path}")
    return path
