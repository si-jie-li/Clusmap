"""Input/output and preprocessing for clusmap.

The goal of this module is to make data loading *robust* and *low-effort*:
the user normally only needs to pass a file path. Format (.h5ad / .tsv / .csv /
.txt / .xlsx / pickled DataFrame) is auto-detected, and a standard
"gene x sample" matrix with a header row and a gene-symbol index column is
parsed without any extra parameters.
"""
from __future__ import annotations

import os
from typing import Optional

import numpy as np
import pandas as pd


def set_working_directory(path: str = ".") -> None:
    """Set the process working directory (kept for backward compatibility)."""
    os.chdir(path)
    print(f"Working directory set to: {os.getcwd()}")


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
_DELIMITERS = {".tsv": "\t", ".txt": "\t", ".csv": ",", ".gene": None}


def _sniff_delimiter(file_path: str, fallback: str = "\t") -> str:
    """Guess a delimiter from the first non-empty line of a text file."""
    with open(file_path, "r") as fh:
        for line in fh:
            if line.strip():
                break
        else:
            return fallback
    # order matters: tab first (most common in genomics), then comma, then ws
    for cand in ("\t", ",", ";"):
        if cand in line:
            return cand
    if "  " in line or " " in line:
        return r"\s+"
    return fallback


def _from_h5ad(file_path: str, layer: Optional[str] = None) -> pd.DataFrame:
    """Load a bulk matrix from an .h5ad file as genes x samples.

    AnnData convention is cells/samples x genes, so we transpose. ``var_names``
    become the gene index, ``obs_names`` the sample columns.
    """
    import anndata as ad
    from scipy.sparse import issparse

    adata = ad.read_h5ad(file_path)
    X = adata.layers[layer] if layer is not None else adata.X
    if issparse(X):
        X = X.toarray()
    df = pd.DataFrame(np.asarray(X).T, index=adata.var_names.astype(str),
                      columns=adata.obs_names.astype(str))
    return df


def import_data(
    file_path: str,
    *,
    file_delimiter: Optional[str] = None,
    header_path: Optional[str] = None,
    header_delimiter: str = "\t",
    str_col_num: Optional[int] = None,
    float_col_num: Optional[int] = None,
    index_col: int = 1,
    layer: Optional[str] = None,
) -> pd.DataFrame:
    """Import a bulk expression matrix as a genes x samples ``DataFrame``.

    For the common case you only need ``file_path``. Everything else is
    inferred. The index is set to gene symbols (needed downstream for marker
    mapping and GO).

    Parameters
    ----------
    file_path :
        Path to the data. ``.h5ad``, ``.tsv``, ``.csv``, ``.txt``, ``.gene``,
        ``.xlsx``/``.xls`` and pickled DataFrames (``.pkl``) are auto-detected.
    file_delimiter :
        Override the delimiter for text files. ``None`` (default) auto-sniffs.
    header_path, header_delimiter :
        Optional external header file (one line of sample names).
    str_col_num :
        Number of leading non-numeric (annotation) columns. ``None`` (default)
        auto-detects how many leading columns are non-numeric.
    float_col_num :
        Number of value columns to keep. ``None`` keeps all remaining columns.
    index_col :
        1-based column used as the gene index. Defaults to column 1.
    layer :
        For ``.h5ad`` inputs, which layer to read (``None`` -> ``adata.X``).
    """
    ext = os.path.splitext(file_path)[1].lower()

    # --- direct, fully-structured formats ---------------------------------- #
    if ext == ".h5ad":
        rna_df = _from_h5ad(file_path, layer=layer)
        print(f"Data imported from h5ad: {rna_df.shape[0]} genes x {rna_df.shape[1]} samples")
        return _finalize(rna_df)
    if ext in (".pkl", ".pickle"):
        rna_df = pd.read_pickle(file_path)
        return _finalize(rna_df)
    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(file_path, header=0)
        return _parse_table(df, str_col_num, float_col_num, index_col,
                            header_path, header_delimiter)

    # --- delimited text ---------------------------------------------------- #
    if file_delimiter is None:
        file_delimiter = _DELIMITERS.get(ext) or _sniff_delimiter(file_path)

    header_tf: Optional[int] = 0 if header_path is None else None
    df = pd.read_csv(file_path, sep=file_delimiter, header=header_tf,
                     engine="python" if file_delimiter == r"\s+" else "c")

    header_list = None
    if header_path is not None:
        with open(header_path, "r") as fh:
            header_list = fh.readline().strip().split(header_delimiter)

    return _parse_table(df, str_col_num, float_col_num, index_col,
                        None, header_delimiter, header_list=header_list)


def _parse_table(df, str_col_num, float_col_num, index_col, header_path,
                 header_delimiter, header_list=None):
    """Slice a raw table into a numeric genes x samples matrix."""
    # auto-detect number of leading annotation columns if not given
    if str_col_num is None:
        str_col_num = 0
        for c in range(df.shape[1]):
            if pd.api.types.is_numeric_dtype(pd.to_numeric(df.iloc[:, c],
                                                           errors="coerce")) \
                    and pd.to_numeric(df.iloc[:, c], errors="coerce").notna().mean() > 0.9:
                break
            str_col_num += 1
        str_col_num = max(str_col_num, 1)

    if float_col_num is None:
        float_col_num = df.shape[1] - str_col_num

    rna_df = df.iloc[:, str_col_num:float_col_num + str_col_num].apply(
        pd.to_numeric, errors="coerce")
    rna_df.index = df.iloc[:, index_col - 1].astype(str).tolist()

    if header_list is not None:
        rna_df.columns = header_list

    print(f"Data imported: {rna_df.shape[0]} genes x {rna_df.shape[1]} samples")
    return _finalize(rna_df)


def _finalize(rna_df: pd.DataFrame) -> pd.DataFrame:
    """Coerce to float, make labels string + unique, drop all-NaN rows."""
    rna_df = rna_df.astype(float)
    rna_df.index = rna_df.index.astype(str)
    rna_df.columns = rna_df.columns.astype(str)
    rna_df = rna_df[~rna_df.index.duplicated(keep="first")]
    rna_df = rna_df.dropna(how="all")
    return rna_df


# --------------------------------------------------------------------------- #
# Preprocessing
# --------------------------------------------------------------------------- #
def preprocess(
    rna_df: pd.DataFrame,
    *,
    start: Optional[int] = None,
    end: Optional[int] = None,
    min_expr: float = 10,
    log_base: Optional[int] = 2,
    plus_constant: float = 0.1,
    std_filter: bool = True,
) -> pd.DataFrame:
    """Slice, low-expression filter, log-transform and drop zero-variance rows.

    All steps have sensible defaults; pass ``min_expr=0`` / ``log_base=None`` to
    skip the respective step.
    """
    if start is not None or end is not None:
        rna_df = rna_df.iloc[slice(start, end), :]

    if min_expr and min_expr > 0:
        before = rna_df.shape[0]
        rna_df = rna_df[rna_df.max(axis=1) > min_expr]
        print(f"Filtered {before - rna_df.shape[0]} genes with max expression <= {min_expr}.")

    if log_base is not None:
        rna_df = np.log(rna_df + plus_constant) / np.log(log_base)
        print(f"Log transform applied (base {log_base}, +{plus_constant}).")

    if std_filter:
        before = rna_df.shape[0]
        rna_df = rna_df.loc[rna_df.std(axis=1) > 0, :]
        print(f"Dropped {before - rna_df.shape[0]} zero-variance genes.")

    print(f"Data preprocessed: {rna_df.shape[0]} genes x {rna_df.shape[1]} samples")
    return rna_df


def extract_color_cat(rna_df: pd.DataFrame, categories_dict: Optional[dict] = None,
                      case_sensitive: bool = False) -> dict:
    """Build column-category vectors by matching keywords in sample names.

    Returns ``{category: [value-or-None per column]}`` aligned to
    ``rna_df.columns`` -- ready to pass to ``bulk_hm(col_cat=...)``.
    """
    if categories_dict is None:
        return {"default": rna_df.columns.tolist()}

    result = {}
    for cat_name, keywords in categories_dict.items():
        extracted = []
        for col in rna_df.columns:
            hay = str(col) if case_sensitive else str(col).lower()
            found = None
            for kw in keywords:
                needle = kw if case_sensitive else kw.lower()
                if needle in hay:
                    found = kw
                    break
            extracted.append(found)
        result[cat_name] = extracted
    return result
