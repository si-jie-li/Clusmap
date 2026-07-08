"""Small shared helpers."""
from __future__ import annotations

import numpy as np
import pandas as pd

_NORM_ALIAS = {
    "zscore": "z_score", "z-score": "z_score", "row_zscore": "z_score",
    "center": "centralize", "centering": "centralize", "mean_center": "centralize",
    "min-max scaling": "min_max", "minmax": "min_max", "max scaling": "max",
    "standard_scale": "min_max",
}
_VALID_NORM = {"z_score", "centralize", "none", "min_max", "max"}


def canon_norm(norm_method: str) -> str:
    nm = (norm_method or "z_score").strip().lower()
    nm = _NORM_ALIAS.get(nm, nm)
    if nm not in _VALID_NORM:
        raise ValueError(
            f"Invalid norm_method={norm_method!r}. Choose from {sorted(_VALID_NORM)} "
            f"(aliases: {sorted(_NORM_ALIAS)})."
        )
    return nm


def normalize_rows(df: pd.DataFrame, norm_method: str) -> pd.DataFrame:
    """Row-wise normalize a genes x samples frame. Returns a new frame.

    Unlike the old code this never relies on clustermap's internal ``z_score``
    flag, so the exact same transform is applied in every figure.
    """
    nm = canon_norm(norm_method)
    df = df.copy()
    if nm == "none":
        return df
    if nm == "centralize":
        return df.sub(df.mean(axis=1), axis=0)
    if nm == "z_score":
        mean = df.mean(axis=1)
        std = df.std(axis=1, ddof=0).replace(0, np.nan)
        return df.sub(mean, axis=0).div(std, axis=0).fillna(0.0)
    if nm == "min_max":
        lo, hi = df.min(axis=1), df.max(axis=1)
        return df.sub(lo, axis=0).div((hi - lo).replace(0, np.nan), axis=0).fillna(0.0)
    if nm == "max":
        return df.div(df.max(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    return df
