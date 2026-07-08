"""clusmap — the whole pipeline as a no-code Streamlit app.

Launch:  streamlit run clusmap/app.py
(use the project env: /opt/anaconda3/envs/clusmap/bin/streamlit run clusmap/app.py)

Every pipeline step is a block with its tunable parameters, a Run button and its
outputs in place; the top panel runs everything in one click. Module curation is
folded into the heatmap block. Data is loaded in the Import block (path or
upload) — or handed off in-memory from a notebook via ``clusmap.launch_app``.
"""
from __future__ import annotations

import os

import streamlit as st

st.set_page_config(page_title="clusmap", layout="wide")
SS = st.session_state


def _maybe_load_handoff():
    """Load an in-memory (rna, state) bundle handed off by ``launch_app`` once.

    Lets you build ``rna`` with the right ``import_data`` params and cluster with
    ``gen_mod`` in a notebook, then jump straight to curation / downstream steps
    in the browser without re-parsing the file.
    """
    if SS.get("_handoff_done"):
        return
    SS["_handoff_done"] = True
    path = os.environ.get("CLUSMAP_SESSION")
    if path and os.path.exists(path):
        import pickle
        with open(path, "rb") as fh:
            bundle = pickle.load(fh)
        SS["rna"] = SS["rna_raw"] = bundle["rna"]
        SS["state"] = bundle["state"]
        SS["outdir"] = bundle.get("outdir", "clusmap_out")
        # the heatmap (ClusterGrid) isn't pickled; click Run on the
        # Clusterheatmap block to render it, then curate.


def main():
    st.title("clusmap")
    _maybe_load_handoff()
    from clusmap import gui
    gui.render_pipeline()


if __name__ == "__main__":
    main()
