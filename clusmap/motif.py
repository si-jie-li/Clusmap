"""Wrap-up of AME motif-enrichment results (from motif_analysis/ SLURM jobs).

The heavy lifting (promoter FASTA extraction + AME against JASPAR/HOCOMOCO/
CIS-BP) runs on an HPC cluster via the scripts in ``motif_analysis/``; this
module collects the per-module ``ame.tsv`` outputs into tidy tables.
"""
from __future__ import annotations

import glob
import os
import re
from typing import Literal, Optional, Tuple

import pandas as pd

_MODULE_RE = re.compile(r"Module(\d+)", re.IGNORECASE)
_NUMERIC = ["rank", "adj_p-value", "E-value", "%TP", "pos"]


def module_motif(
    root_dir: str,
    outdir: str = ".",
    *,
    mode: Literal["top_n", "p", "e", "pe"] = "top_n",
    top_n: int = 10,
    p_thresh: Optional[float] = 0.05,
    e_thresh: Optional[float] = 1.0,
    name_col: str = "motif_alt_ID",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Collect & filter AME results from ``root_dir/Module*/ame.tsv``.

    Writes a long-format table, a per-module summary, and adj-p / E-value
    matrices. ``mode`` selects filtering: ``top_n`` keeps the best N per module;
    ``p`` / ``e`` / ``pe`` apply adj-p and/or E-value thresholds.
    """
    files = sorted(glob.glob(os.path.join(root_dir, "Module*", "ame.tsv")))
    if not files:
        raise FileNotFoundError(f"No ame.tsv under {root_dir}/Module*/ame.tsv")
    os.makedirs(outdir, exist_ok=True)

    frames = []
    for f in files:
        m = _MODULE_RE.search(f)
        df = pd.read_csv(f, sep="\t", dtype=str)
        df.columns = [c.strip() for c in df.columns]
        for col in _NUMERIC:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        df["module"] = int(m.group(1)) if m else None
        keep = ["module", "rank", "motif_ID", "motif_alt_ID",
                "adj_p-value", "E-value", "%TP", "pos"]
        frames.append(df[[c for c in keep if c in df.columns]])

    raw = pd.concat(frames, ignore_index=True).sort_values(["module", "rank"])

    all_df = raw.copy()
    if mode == "top_n":
        all_df = all_df.groupby("module", as_index=False).head(top_n)
    elif mode == "p":
        all_df = all_df[all_df["adj_p-value"] < p_thresh]
    elif mode == "e":
        all_df = all_df[all_df["E-value"] < e_thresh]
    elif mode == "pe":
        all_df = all_df[(all_df["adj_p-value"] < p_thresh) & (all_df["E-value"] < e_thresh)]
    else:
        raise ValueError(f"Unknown mode: {mode}")

    all_df.to_csv(os.path.join(outdir, "module_motif_longformat.tsv"), sep="\t", index=False)

    # per-module compact summary
    summary = (all_df.groupby("module")[name_col].apply(list).reset_index())
    motif_wide = pd.DataFrame(summary[name_col].tolist())
    summary_df = pd.concat([summary[["module"]], motif_wide], axis=1)
    summary_df.to_csv(os.path.join(outdir, "module_motif_summary.tsv"),
                      sep="\t", index=False, header=False)

    # adj-p / E-value matrices (using original, unfiltered values of passing motifs)
    valid = all_df[name_col].dropna().unique()
    base = raw[raw[name_col].isin(valid)]
    base.pivot_table(index="module", columns=name_col, values="adj_p-value",
                     aggfunc="min").to_csv(
        os.path.join(outdir, "module_motif_adjpval_matrix.tsv"), sep="\t")
    base.pivot_table(index="module", columns=name_col, values="E-value",
                     aggfunc="min").to_csv(
        os.path.join(outdir, "module_motif_eval_matrix.tsv"), sep="\t")

    return all_df, summary_df
