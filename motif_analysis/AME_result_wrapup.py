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

    Parameters:
    - root_dir (str): Directory containing per-module folders (e.g. .../AME_out_HOCOMOCO)
    - outdir (str, optional): Directory to save output files. Default to ".".
    - mode (select from {"top_n","p","e","pe"}, optional): Filtering mode. Default to "top_n".
    - top_n (int, optional): If mode="top_n", keep top_n motifs per module (based on rank in ame.tsv). Default to 10.
    - p_thresh (float or None, optional): If mode="p" or "pe", p-value threshold. Default to 0.05.
    - e_thresh (float or None, optional): If mode="e" or "pe", E-value threshold. Default to 1.

    Saves:
    - module_motif_longformat.tsv: all_df long-format table; module, motif_ID, motif_alt_ID, adj_p-value, E-value, %TP, pos, rank
    - module_motif_summary.tsv: Per-module summar, one row per module, compact string of motif hits.

    Returns:
    - pd.DataFrame: all_df long-format table.
    - pd.DataFrame: Per-module summary table.
    """

    # extrct ame.tsv files
    ame_files = sorted(glob.glob(os.path.join(root_dir, "Module*", "ame.tsv")))
    if not ame_files:
        raise FileNotFoundError(f"No ame.tsv files found under: {root_dir}/Module*/ame.tsv")
    os.makedirs(outdir, exist_ok=True)
    long_out_path = os.path.join(outdir, "module_motif_longformat.tsv")
    summary_out_path = os.path.join(outdir, "module_motif_summary.tsv")

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
        df = df[cols]
        dfs.append(df) # list of dataframes

    all_df = pd.concat(dfs, ignore_index=True)  # long format of all modules
    
    all_df = all_df.sort_values(["module", "rank"], ascending=[True, True])
    if mode == "top_n":
        all_df = all_df.groupby("module", as_index=False).head(top_n)

    elif mode == "p":
        all_df = all_df[all_df["p-value"] < p_thresh]

    elif mode == "e":
        all_df = all_df[all_df["E-value"] < e_thresh]

    elif mode == "pe":
        all_df = all_df[(all_df["p-value"] < p_thresh) & (all_df["E-value"] < e_thresh)]

    else:
        raise ValueError(f"Unknown mode: {mode}")
    all_df.to_csv(long_out_path, sep="\t", index=False)

    def _fmt_hit(row) -> str:
        motif_id = row["motif_ID"]
        name = row["motif_alt_ID"]
        e = row["E-value"]
        ap = row["adj_p-value"]
        tp = row["%TP"]
        pos = row["pos"]
        return f"{motif_id}(TF={name},E={e},adjP={ap},%TP={tp},pos={pos})"


    tmp = all_df.copy()
    tmp["motif"] = tmp.apply(_fmt_hit, axis=1)

    grouped = (
        tmp.groupby("module", as_index=False)   # group by module number, but "module" column do not set as index
           .agg(
               pos=("pos", "first"),        # pos should be constant within a module for a given run, thus take the first pos value of that module group 
               n_motif=("motif_ID", "count"),  # count how many different motif_ID, number of motif in this result
               motif_list=("motif", list)
           )
           .sort_values("module")
    )
    motif_df = pd.DataFrame(grouped["motif_list"].tolist())
    summary_df = pd.concat([grouped[["module", "pos", "n_motif"]],motif_df], axis=1)
    summary_df.to_csv(summary_out_path, sep="\t", index=False,header=False)

    return all_df, summary_df

# -----------------------------
if __name__ == "__main__":
    root = "/export/home/lisijie/motif_analysis/AME_out_CIS-BP"

    # 1) top N per module
    # module_motif(root, "./wrap_result_CIS-BP")


    # 2) E-value threshold
    module_motif(root, "./wrap_result_CIS-BP", mode="e", e_thresh=1e-2)

    # 3) p-value AND E-value
    # collect_ame_hits(root, out_long, mode="pe", p_thresh=1e-6, e_thresh=1e-5)
