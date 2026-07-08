import os
import re
import glob
import pandas as pd
from typing import Optional, Literal, Tuple
import numpy as np
import csv

def module_motif(root_dir: str, outdir: str = ".",
    mode: Literal["top_n", "p", "e", "pe"] = "top_n",
    top_n: int = 10, p_thresh: Optional[float] = 0.05, e_thresh: Optional[float] = 1) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Collect and filter AME results from: root_dir/Module*/ame.tsv
    TF is from the meme file, may be TF isoform or family name, or do not find in the meme used when doing the enrichment analysis
    ONLY CISBP database have motif name in motif_alt_ID column
    Parameters:
    - root_dir (str): Directory containing per-module folders (e.g. .../AME_out_HOCOMOCO)
    - outdir (str, optional): Directory to save output files. Default to ".".
    - mode (select from {"top_n","p","e","pe"}, optional): Filtering mode. Default to "top_n".
    - top_n (int, optional): If mode="top_n", keep top_n motifs per module (based on rank in ame.tsv). Default to 10.
    - p_thresh (float or None, optional): If mode="p" or "pe", p-value threshold. Default to 0.05.
    - e_thresh (float or None, optional): If mode="e" or "pe", E-value threshold. Default to 1.

    Saves:
    - module_motif_summary.tsv: Per-module summary, one row per module.

    Returns:
    - pd.DataFrame: all_df long-format table.
    - pd.DataFrame: Per-module summary table.
    """

    # extrct ame.tsv files
    ame_files = sorted(glob.glob(os.path.join(root_dir, "Module*", "ame.tsv")))
    if not ame_files:
        raise FileNotFoundError(f"No ame.tsv files found under: {root_dir}/Module*/ame.tsv")
    os.makedirs(outdir, exist_ok=True)
    summary_out_path = os.path.join(outdir, "module_motif_summary_lite.tsv")
    pval_out_path = os.path.join(outdir, "module_motif_adjpval_matrix.tsv")
    eval_out_path = os.path.join(outdir, "module_motif_eval_matrix.tsv")


    # regex to extract module number
    module_re = re.compile(r"Module(\d+)", re.IGNORECASE)

    dfs = []
    for f in ame_files:
        m = module_re.search(f)
        module = int(m.group(1)) if m else None  # if can not capture the number, set to None

        df = pd.read_csv(f, sep="\t", dtype=str,index_col=None)  # read as string first
        df.columns = [c.strip() for c in df.columns]  # safety: strip column spaces, necessary!

        numeric_cols = ["rank", "adj_p-value", "E-value", "%TP", "pos"]
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")  # Convert numeric columns (coerce errors to NaN)

        df["module"] = module
        cols = ["module","rank","motif_ID","motif_alt_ID","adj_p-value","E-value","%TP","pos"]
        df = df[[c for c in cols if c in df.columns]]
        dfs.append(df) # list of dataframes

    raw_df = pd.concat(dfs, ignore_index=True)  # keep unfiltered dataframe for matrix generation
    raw_df = raw_df.sort_values(["module", "rank"], ascending=[True, True])
    
    all_df = raw_df.copy()
    if mode == "top_n":
        all_df = all_df.groupby("module", as_index=False).head(top_n)
        print(f"First {top_n} remaind, rank basing on confidence")

    elif mode == "p":
        all_df = all_df[all_df["adj_p-value"] < p_thresh]
        print(f"Motif with adj_p-value < {p_thresh} remaind")
    elif mode == "e":
        all_df = all_df[all_df["E-value"] < e_thresh]
        print(f"Motif with E value < {e_thresh} remaind")

    elif mode == "pe":
        all_df = all_df[(all_df["adj_p-value"] < p_thresh) & (all_df["E-value"] < e_thresh)]
        print(f"Motif with adj_p-value < {p_thresh} AND E value < {e_thresh} remaind")

    else:
        raise ValueError(f"Unknown mode: {mode}")

    # --- 1. Original output: per-module motif list summary ---
    grouped = all_df.groupby("module")["motif_alt_ID"].apply(list).reset_index()
    motif_df = pd.DataFrame(grouped["motif_alt_ID"].tolist())
    summary_df = pd.concat([grouped[["module"]], motif_df], axis=1)
    summary_df.to_csv(summary_out_path, sep="\t", index=False, header=False)

    # --- 2. New output: adj_p-value and E-value matrix ---
    valid_motifs = all_df["motif_alt_ID"].dropna().unique()  # all_df is filtered by p or e or top_n
    matrix_base = raw_df[raw_df["motif_alt_ID"].isin(valid_motifs)]  # raw_df is unfiltered, keep those motif-module pairs that pass the filter in all_df, but keep their original p and E values for matrix generation

    # Use min aggfunc to robustly handle any duplicate motif names in the same module
    pval_matrix = matrix_base.pivot_table(index="module", columns="motif_alt_ID", values="adj_p-value", aggfunc="min").fillna("NaN")
    eval_matrix = matrix_base.pivot_table(index="module", columns="motif_alt_ID", values="E-value", aggfunc="min").fillna("NaN")

    pval_matrix.to_csv(pval_out_path, sep="\t")
    eval_matrix.to_csv(eval_out_path, sep="\t")

    return all_df, summary_df

# -----------------------------
if __name__ == "__main__":
    root = "/export/home/lisijie/motif_analysis/AME_out_CIS-BP"

    # 1) top N per module
    module_motif(root, "./wrap_result_CIS-BP303")

    # 2) E-value threshold
    # module_motif(root, "./wrap_result_JASPAR0303", mode="e", e_thresh=1e-2)

    # 3) p-value AND E-value
    # collect_ame_hits(root, out_long, mode="pe", p_thresh=1e-6, e_thresh=1e-5)
