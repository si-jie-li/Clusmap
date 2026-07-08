"""Figure building + selection logic for the interactive heatmap.

This module is deliberately free of any Streamlit import so the plotting and
selection-mapping logic can be unit-tested headlessly. The Streamlit UI lives in
``clusmap/app.py`` and calls into here.
"""
from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.colors as mcolors

from .state import ModuleState
from .util import normalize_rows


def module_color_lut(state: ModuleState, palette: str = "hsv") -> Dict[int, str]:
    """Map heatmap module id (1-based, 0=unassigned) -> hex color."""
    ids = sorted({int(i) for i in state.hm_labels})
    non_zero = [i for i in ids if i != 0]
    pal = sns.color_palette(palette, n_colors=max(len(non_zero), 1))
    lut = {i: mcolors.to_hex(pal[k]) for k, i in enumerate(non_zero)}
    lut[0] = "#ffffff"
    return lut


def ordered_frame(state: ModuleState, rna_df: pd.DataFrame, norm_method: str = "z_score"):
    """Return (normalized values in heatmap row order, ordered gene names,
    ordered hm module ids)."""
    order = state.order
    genes = state.genes[order]
    mods = state.hm_labels[order]
    df = normalize_rows(rna_df.loc[state.genes], norm_method)
    values = df.values[order]
    return values, genes, mods


def build_figure(state: ModuleState, rna_df: pd.DataFrame, *,
                 norm_method: str = "z_score", palette: str = "hsv",
                 cmap: str = "turbo", highlight_genes: Sequence[str] = (),
                 row_height: float = 14.0):
    """Build the Plotly figure: thin module band + expression heatmap.

    Imported lazily so the rest of the package does not require plotly.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    values, genes, mods = ordered_frame(state, rna_df, norm_method)
    lut = module_color_lut(state, palette)
    n = len(genes)

    fig = make_subplots(rows=1, cols=2, column_widths=[0.05, 0.95],
                        shared_yaxes=True, horizontal_spacing=0.01)

    # --- module band (col 1): crisp discrete colored strip ----------------- #
    # stepped colorscale so each module id maps to exactly one flat colour
    # (no interpolation between adjacent rows), matching bulk_hm's row band.
    uniq = sorted(set(int(m) for m in mods))
    code = {m: k for k, m in enumerate(uniq)}
    nseg = len(uniq)
    band_z = np.array([[code[int(m)] + 0.5] for m in mods])
    colorscale = []
    for k, m in enumerate(uniq):
        colorscale.append([k / nseg, lut[m]])
        colorscale.append([(k + 1) / nseg, lut[m]])
    fig.add_trace(go.Heatmap(z=band_z, x=["module"], y=list(range(n)),
                             zmin=0, zmax=nseg, colorscale=colorscale, showscale=False,
                             hovertext=[[f"module {int(m)}"] for m in mods],
                             hoverinfo="text", name="module"), row=1, col=1)

    # --- expression heatmap (col 2) ---------------------------------------- #
    # robust colour limits (2nd/98th percentile) === bulk_hm's robust=True,
    # so outliers don't wash out the palette.
    finite = values[np.isfinite(values)]
    if finite.size:
        zmin, zmax = (float(x) for x in np.percentile(finite, [2, 98]))
        if zmin == zmax:
            zmin, zmax = float(finite.min()), float(finite.max() or 1.0)
    else:
        zmin, zmax = None, None
    fig.add_trace(go.Heatmap(z=values, x=list(rna_df.columns), y=list(range(n)),
                             zmin=zmin, zmax=zmax,
                             colorscale=cmap, colorbar=dict(title=norm_method, len=0.4,
                                                            y=0.8),
                             customdata=np.array(genes).reshape(-1, 1),
                             hovertemplate="gene=%{customdata[0]}<br>"
                                           "sample=%{x}<br>val=%{z:.2f}<extra></extra>",
                             name="expr"), row=1, col=2)

    # gene labels live on the y axis as indices; show names sparsely or on hover
    fig.update_yaxes(autorange="reversed", showticklabels=False)
    fig.update_xaxes(showticklabels=True, tickangle=90, row=1, col=2)
    fig.update_xaxes(showticklabels=False, row=1, col=1)

    # highlight searched genes with horizontal lines
    name_to_row = {g: i for i, g in enumerate(genes)}
    for g in highlight_genes:
        i = name_to_row.get(g) or name_to_row.get(str(g))
        if i is None:
            low = {str(x).lower(): k for k, x in enumerate(genes)}
            i = low.get(str(g).lower())
        if i is not None:
            fig.add_hline(y=i, line=dict(color="black", width=1, dash="dot"),
                          row=1, col=2)

    fig.update_layout(height=max(400, min(1600, row_height * np.sqrt(n) * 4)),
                      margin=dict(l=10, r=10, t=30, b=120),
                      title=f"{n} genes · {state.n_modules} modules · {norm_method}")
    return fig, genes, mods


# --------------------------------------------------------------------------- #
# selection -> meaning (pure helpers, testable without plotly/streamlit)
# --------------------------------------------------------------------------- #
def rows_from_selection(points: List[dict]) -> List[int]:
    """Extract ordered-row indices (y) from a plotly selection/click payload."""
    rows = set()
    for p in points or []:
        y = p.get("y")
        if y is not None:
            rows.add(int(round(y)))
    return sorted(rows)


def genes_at_rows(genes: Sequence[str], rows: Sequence[int]) -> List[str]:
    return [genes[r] for r in rows if 0 <= r < len(genes)]


def modules_at_rows(mods: Sequence[int], rows: Sequence[int]) -> List[int]:
    """Distinct module ids touched by the selected rows (order preserved)."""
    seen, out = set(), []
    for r in rows:
        if 0 <= r < len(mods):
            m = int(mods[r])
            if m not in seen:
                seen.add(m); out.append(m)
    return out


def are_neighbors(mods_ordered: Sequence[int], a: int, b: int) -> bool:
    """True if modules a and b are adjacent blocks in heatmap order."""
    blocks = []
    for m in mods_ordered:
        if not blocks or blocks[-1] != m:
            blocks.append(int(m))
    if a not in blocks or b not in blocks:
        return False
    return abs(blocks.index(a) - blocks.index(b)) == 1
