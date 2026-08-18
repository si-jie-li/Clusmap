"""Input/output and preprocessing for clusmap.

The goal of this module is to make data loading *robust* and *low-effort*:
the user normally only needs to pass a file path. Format (.h5ad / .tsv / .csv /
.txt / .xlsx / pickled DataFrame) is auto-detected, and a standard
"gene x sample" matrix with a header row and a gene-symbol index column is
parsed without any extra parameters.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

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


# --------------------------------------------------------------------------- #
# Orientation / structure auto-detection for delimited text
# --------------------------------------------------------------------------- #
# Gene symbols that Excel/ pandas would wrongly coerce to a number (dates like
# 1-Mar -> timestamp, or SEPT2/MARCH1 -> 2/1). We treat these as *non-numeric*
# so the numeric-block detector never mistakes a gene-name column for data.
_DATELIKE_GENE = re.compile(r"^\d{1,2}[-/][A-Za-z]{3}$")
_NUMFAKE_GENE = re.compile(
    r"^(SEPT|MARCH|DEC|FEB|OCT|NOV|JUN|JUL|AUG|APR)\d+$", re.IGNORECASE)


@dataclass
class TableStructure:
    """What auto-detection inferred about a messy delimited table.

    All indices are into the raw (header=None) string grid. ``transposed``
    means the raw numeric block was samples×genes and was flipped to
    genes×samples. ``confidence`` is ``high``/``medium``/``low``.
    """
    data_top: int
    data_bottom: int
    data_left: int
    data_right: int
    header_row: int = -1          # row holding sample names (``-1`` = none found)
    gene_col: int = -1            # column holding gene names (``-1`` = none)
    transposed: bool = False
    confidence: str = "high"
    notes: List[str] = field(default_factory=list)

    @property
    def n_genes(self) -> int:
        return self.data_bottom - self.data_top + 1

    @property
    def n_samples(self) -> int:
        return self.data_right - self.data_left + 1


def _is_number_cell(value: str) -> bool:
    """True if ``value`` reads as a number and is NOT a date-faked gene name."""
    s = str(value).strip()
    if not s or s.lower() in ("nan", "na", "null", "none", "."):
        return False
    if _DATELIKE_GENE.match(s) or _NUMFAKE_GENE.match(s):
        return False
    try:
        float(s)
        return True
    except ValueError:
        return False


def _numeric_mask(raw: pd.DataFrame) -> np.ndarray:
    """Boolean mask (R×C) of which cells in the string grid are numeric."""
    return np.vectorize(_is_number_cell)(raw.values.astype(str))


def _find_dense_block(mask: np.ndarray) -> Tuple[int, int, int, int]:
    """Locate the largest contiguous numeric rectangle in the boolean mask.

    Returns ``(top, bottom, left, right)`` inclusive. Strategy: for every pair
    of rows, find columns that are numeric across the whole band, then keep the
    band×columnspan with the best score (density × area). O(R²·C) but the
    published-figure files this targets are small (thousands of cells).
    """
    nrows, ncols = mask.shape
    best = (0, 0, 0, 0, -1)   # top, bottom, left, right, score
    # prefix sums per column for fast "all numeric in this band" queries
    cum = np.zeros((nrows + 1, ncols), dtype=int)
    cum[1:] = mask.astype(int)
    np.cumsum(cum, axis=0, out=cum)
    for top in range(nrows):
        for bottom in range(top, nrows):
            band = bottom - top + 1
            # columns numeric across the entire band
            col_count = cum[bottom + 1] - cum[top]
            full = (col_count == band)
            # longest contiguous run of such columns
            run_best = run_left = run_right = 0
            run_start = -1
            for c in range(ncols):
                if full[c]:
                    if run_start == -1:
                        run_start = c
                    cur_len = c - run_start + 1
                    if cur_len > run_best:
                        run_best = cur_len
                        run_left, run_right = run_start, c
                else:
                    run_start = -1
            if run_best == 0:
                continue
            area = band * run_best
            # density inside the candidate rectangle
            density = mask[top:bottom + 1, run_left:run_right + 1].mean()
            score = density * area
            if score > best[4]:
                best = (top, bottom, run_left, run_right, score)
    return best[0], best[1], best[2], best[3]


def _infer_header_and_genecol(top: int, left: int) -> Tuple[int, int]:
    """Pick the sample-header row and gene-name column adjacent to the block."""
    header_row = top - 1 if top > 0 else -1
    gene_col = left - 1 if left > 0 else -1
    return header_row, gene_col


# Words that appear in the corner cell diagonally above-left of the data block
# (the "label of labels", e.g. the ``gene`` in ``gene,0,5,10``). If the block's
# top row is really a numeric header row (sample names that happen to be
# numbers), this word is the giveaway that the row belongs to the header.
_HEADER_CORNER = {"", "gene", "genes", "sample", "samples", "id", "probe_id",
                  "target", "name", "symbol", "description"}


def _refine_block_edges(
    raw: pd.DataFrame, top: int, bottom: int, left: int,
) -> Tuple[int, int]:
    """Pull an absorbed numeric header row / gene column out of the block.

    A header row whose sample names are numbers (time points ``0,5,10``) is
    otherwise indistinguishable from a data row, so the dense-block finder
    merges it in. Detect it: the cell diagonally above-left of the block is a
    non-numeric header-corner word while the column below it holds gene names.
    """
    gene_col = left - 1
    if top == 0 and gene_col >= 0:
        corner = str(raw.values[top, gene_col]).strip().lower()
        if corner in _HEADER_CORNER:
            below = raw.values[top + 1:bottom + 1, gene_col].astype(str)
            if below.size and np.vectorize(_is_number_cell)(below).mean() < 0.3:
                return top + 1, left
    return top, left


def _infer_orientation(
    raw: pd.DataFrame,
    top: int, bottom: int, left: int, right: int,
    header_row: int, gene_col: int,
) -> Tuple[bool, str, List[str]]:
    """Decide whether the raw block is samples×genes (needs transpose).

    Evidence, in priority order:
      1. annotation position: a gene-name column hanging off the left labels
         the *rows* (genes already = rows); a gene-name row above labels the
         *columns* (genes = columns -> transpose).
      2. shape: the gene axis is usually much longer than the sample axis.
    Returns ``(transposed, confidence, notes)``.
    """
    notes: List[str] = []
    nrows = bottom - top + 1
    ncols = right - left + 1

    # --- evidence 1: does the left-adjacent column look like gene names? ---
    left_is_labels = False
    if gene_col >= 0:
        col_vals = raw.values[top:bottom + 1, gene_col].astype(str)
        frac_non_numeric = 1.0 - np.vectorize(_is_number_cell)(col_vals).mean()
        left_is_labels = frac_non_numeric >= 0.6

    # --- evidence 2: does the row above the block look like gene names? ---
    top_is_labels = False
    if header_row >= 0:
        row_vals = raw.values[header_row, left:right + 1].astype(str)
        frac_non_numeric = 1.0 - np.vectorize(_is_number_cell)(row_vals).mean()
        top_is_labels = frac_non_numeric >= 0.6

    if left_is_labels and not top_is_labels:
        notes.append("gene-name column on the left -> rows are genes")
        return False, "high", notes
    if top_is_labels and not left_is_labels:
        notes.append("gene-name row above -> columns are genes (will transpose)")
        return True, "high", notes

    # --- evidence 3: shape fallback ---
    if nrows > ncols * 1.5:
        notes.append(f"block taller than wide ({nrows}×{ncols}) -> rows are genes")
        return False, "medium", notes
    if ncols > nrows * 1.5:
        notes.append(f"block wider than tall ({nrows}×{ncols}) -> columns are genes (will transpose)")
        return True, "medium", notes

    # --- ambiguous ---
    notes.append(
        f"orientation ambiguous ({nrows}×{ncols}); defaulting to rows=genes, "
        "override with str_col_num/index_col if wrong")
    return False, "low", notes


def _auto_detect_structure(
    raw: pd.DataFrame,
) -> Optional[Tuple[TableStructure, pd.DataFrame, List[str], List[str]]]:
    """Full detection pipeline on a header=None string grid.

    Returns ``(structure, numeric_block, gene_names, sample_names)`` or ``None``
    if no numeric block was found. ``numeric_block`` is the raw string slice of
    the expression matrix (still in its original orientation).
    """
    mask = _numeric_mask(raw)
    if mask.mean() < 0.05:
        return None
    top, bottom, left, right = _find_dense_block(mask)
    if bottom < top or right < left:
        return None

    top, left = _refine_block_edges(raw, top, bottom, left)
    header_row, gene_col = _infer_header_and_genecol(top, left)
    transposed, confidence, notes = _infer_orientation(
        raw, top, bottom, left, right, header_row, gene_col)

    # gene names: column left of block (rows=genes) or row above (cols=genes)
    if gene_col >= 0 and not transposed:
        gene_names = raw.values[top:bottom + 1, gene_col].astype(str).tolist()
    elif transposed and header_row >= 0:
        gene_names = raw.values[header_row, left:right + 1].astype(str).tolist()
    else:
        gene_names = [f"gene_{i}" for i in range(
            bottom - top + 1 if not transposed else right - left + 1)]

    # sample names: row above block (rows=genes) or column left (cols=genes)
    if (header_row >= 0) and not transposed:
        sample_names = raw.values[header_row, left:right + 1].astype(str).tolist()
    elif transposed and gene_col >= 0:
        sample_names = raw.values[top:bottom + 1, gene_col].astype(str).tolist()
    else:
        sample_names = [f"sample_{j}" for j in range(
            right - left + 1 if not transposed else bottom - top + 1)]

    block = raw.iloc[top:bottom + 1, left:right + 1]

    structure = TableStructure(
        data_top=top, data_bottom=bottom, data_left=left, data_right=right,
        header_row=header_row, gene_col=gene_col, transposed=transposed,
        confidence=confidence, notes=notes)
    return structure, block, gene_names, sample_names


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

    # When the user pins any structural parameter, honour it exactly (legacy
    # path). Otherwise auto-detect orientation, header row and gene column.
    use_legacy = (str_col_num is not None or float_col_num is not None
                  or header_path is not None or index_col != 1)

    if use_legacy:
        header_tf: Optional[int] = 0 if header_path is None else None
        df = pd.read_csv(file_path, sep=file_delimiter, header=header_tf,
                         engine="python" if file_delimiter == r"\s+" else "c")
        header_list = None
        if header_path is not None:
            with open(header_path, "r") as fh:
                header_list = fh.readline().strip().split(header_delimiter)
        return _parse_table(df, str_col_num, float_col_num, index_col,
                            None, header_delimiter, header_list=header_list)

    # --- auto-detect path: read the whole grid as raw strings ----------------
    raw = pd.read_csv(file_path, sep=file_delimiter, header=None, dtype=str,
                      na_filter=False,
                      engine="python" if file_delimiter == r"\s+" else "c")
    detected = _auto_detect_structure(raw)
    if detected is None:
        raise ValueError(
            "Could not find a contiguous numeric expression block in "
            f"{file_path!r}. If the file has a non-standard layout, pin it "
            "explicitly with str_col_num / index_col / float_col_num.")
    structure, block, gene_names, sample_names = detected

    print(f"[auto-detect] header_row={structure.header_row}, "
          f"gene_col={structure.gene_col}, transposed={structure.transposed}, "
          f"confidence={structure.confidence}")
    for note in structure.notes:
        print(f"             - {note}")

    rna_df = _build_from_detection(structure, block, gene_names, sample_names)
    rna_df.attrs["clusmap_detected"] = structure
    print(f"Data imported: {rna_df.shape[0]} genes x {rna_df.shape[1]} samples")
    return _finalize(rna_df)


def _build_from_detection(
    structure: TableStructure, block: pd.DataFrame,
    gene_names: List[str], sample_names: List[str],
) -> pd.DataFrame:
    """Assemble a genes×samples DataFrame from a detection result.

    ``block`` is the raw string slice of the numeric region in its *original*
    orientation; transpose it if the detector decided the genes ran along the
    columns.
    """
    values = block.apply(pd.to_numeric, errors="coerce")
    if structure.transposed:
        values = values.T
    values = values.astype(float)
    values.index = [str(g) for g in gene_names]
    values.columns = [str(s) for s in sample_names]
    return values


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
