"""Launch the interactive Streamlit editor with in-memory objects.

Recommended workflow — you've already loaded the matrix with the right
``import_data`` parameters and clustered with ``gen_mod``; hand those live
objects straight to the app so it never re-parses the file (which is what breaks
on non-standard header/index layouts):

    import clusmap as cm
    rna   = cm.import_data("weird_layout.tsv", str_col_num=2, index_col=2, header_path="h.txt")
    rna   = cm.preprocess(rna)
    state = cm.gen_mod(rna)
    cm.launch_app(rna, state)          # opens the browser app, already populated

Works from a Jupyter notebook or a plain script. Returns the Streamlit
subprocess handle; call ``.terminate()`` to stop it.

Only stdlib is imported at module load, so ``import clusmap`` never requires
Streamlit — it is resolved lazily when you actually launch.
"""
from __future__ import annotations

import os
import pickle
import subprocess
import sys
import tempfile
from pathlib import Path


def launch_app(rna_df, state, hm=None, *, outdir: str = "clusmap_out",
               port: int = 8501, headless: bool = False):
    """Pickle (rna_df, state) to a temp bundle and start ``streamlit run app.py``.

    Parameters
    ----------
    rna_df, state :
        The objects returned by ``import_data``/``preprocess`` and ``gen_mod``.
    hm :
        Optional ``bulk_hm`` result. If given, its figure is rendered as a
        static reference image inside the app (so you can compare the
        publication-quality heatmap against the interactive one).
    outdir :
        Default directory the app's Save button writes to.
    port :
        Port for the Streamlit server.
    headless :
        If True, don't auto-open a browser (useful on remote/HPC; open the
        printed URL yourself, or tunnel the port).
    """
    if state is None or rna_df is None:
        raise ValueError("Pass both rna_df and state (run gen_mod first).")

    bundle_dir = Path(tempfile.mkdtemp(prefix="clusmap_"))
    bundle_path = bundle_dir / "session.pkl"
    reference_png = None
    if hm is not None:
        reference_png = str(bundle_dir / "reference_heatmap.png")
        fig = getattr(hm, "fig", hm)            # ClusterGrid.fig, or a Figure
        fig.savefig(reference_png, bbox_inches="tight",
                    dpi=getattr(hm, "_clusmap_dpi", 150))
    with open(bundle_path, "wb") as fh:
        pickle.dump({"rna": rna_df, "state": state, "outdir": outdir,
                     "reference_png": reference_png}, fh)

    app_path = str(Path(__file__).with_name("app.py"))
    env = dict(os.environ, CLUSMAP_SESSION=str(bundle_path))
    cmd = [sys.executable, "-m", "streamlit", "run", app_path,
           "--server.port", str(port)]
    if headless:
        cmd += ["--server.headless", "true"]

    url = f"http://localhost:{port}"
    print(f"Launching clusmap interactive editor at {url}\n"
          f"  ({rna_df.shape[0]} genes × {rna_df.shape[1]} samples, "
          f"{state.n_modules} modules)\n"
          f"  Stop it with <returned process>.terminate(), or Ctrl-C in the cell.")
    return subprocess.Popen(cmd, env=env)
