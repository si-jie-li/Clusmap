"""Spec-driven Streamlit renderer — the whole clusmap pipeline, no code.

One page, one block per pipeline step (specs from :mod:`clusmap.gui_steps`):
each block shows its tunable parameters with the default as a gray placeholder
(common ones inline, rarely-used ones under an *Advanced* expander), a ▶ Run
button, and its outputs (tables / figures) in place. Status messages pop up in a
closable dialog rather than piling up at the bottom. A top control panel runs the
whole pipeline in one click and toggles the optional steps. Module curation is
folded into the heatmap block.
"""
from __future__ import annotations

import json
import os
import tempfile
from typing import Any, Dict, List

import streamlit as st

from .gui_steps import (STEPS, STEP_BY_KEY, ParamSpec, StepSpec,
                        ordered_module_blocks, curate_merge, curate_split,
                        curate_reassign, save_curation)

SS = st.session_state

# steps whose first file argument can be uploaded instead of typed
_UPLOAD_FIELD = {"import": "file_path", "markers": "marker_path",
                 "pseudobulk": "h5ad_file", "preservation": "test_file",
                 "project": "new_file", "trait": "trait_file"}


# --------------------------------------------------------------------------- #
# notifications -> closable popup dialog
# --------------------------------------------------------------------------- #
def _note(label: str, level: str, msg: str):
    SS.setdefault("_pending_notes", []).append((label, level, msg))


@st.dialog("clusmap")
def _notify_dialog(notes):
    for label, level, msg in notes:
        fn = getattr(st, level if level in ("success", "info", "warning", "error") else "info")
        fn(f"**{label}** — {msg}")
    if st.button("Close", use_container_width=True):
        st.rerun()


def _flush_notes():
    notes = SS.get("_pending_notes") or []
    if notes:
        SS["_pending_notes"] = []          # snapshot already captured; don't re-open next run
        _notify_dialog(notes)


# --------------------------------------------------------------------------- #
# parameter widgets — default shown as gray placeholder, value parsed back
# --------------------------------------------------------------------------- #
def _placeholder(p: ParamSpec) -> str:
    if p.kind == "list_str":
        return ", ".join(map(str, p.default or [])) or "(comma-separated)"
    if p.kind in ("intnone", "floatnone") and p.default is None:
        return "none"
    if p.kind == "json":
        return json.dumps(p.default) if p.default else "{}"
    return "" if p.default in (None, "") else str(p.default)


def _coerce(txt: str, p: ParamSpec):
    txt = (txt or "").strip()
    if txt == "":
        return p.default
    if p.kind == "list_str":
        return [x.strip() for x in txt.split(",") if x.strip()]
    if p.kind == "json":
        try:
            return json.loads(txt)
        except json.JSONDecodeError:
            return p.default
    if p.kind in ("int",):
        try: return int(float(txt))
        except ValueError: return p.default
    if p.kind in ("float",):
        try: return float(txt)
        except ValueError: return p.default
    if p.kind == "intnone":
        if txt.lower() in ("none", "null"): return None
        try: return int(float(txt))
        except ValueError: return p.default
    if p.kind == "floatnone":
        if txt.lower() in ("none", "null"): return None
        try: return float(txt)
        except ValueError: return p.default
    return txt          # str / path


def _widget(step_key: str, p: ParamSpec) -> Any:
    key = f"{step_key}.{p.name}"
    if p.kind == "bool":
        return st.toggle(p.name, value=bool(p.default), key=key, help=p.help)
    if p.kind == "enum":
        idx = p.choices.index(p.default) if p.default in (p.choices or []) else 0
        return st.selectbox(p.name, p.choices, index=idx, key=key, help=p.help)
    txt = st.text_input(p.name, value="", placeholder=_placeholder(p), key=key, help=p.help)
    return _coerce(txt, p)


def _collect_params(step: StepSpec) -> Dict[str, Any]:
    values: Dict[str, Any] = {pp.name: pp.default for pp in step.params}
    common = [pp for pp in step.params if not pp.advanced]
    advanced = [pp for pp in step.params if pp.advanced]

    def _render(plist):
        visible = [pp for pp in plist if pp.show_if is None or pp.show_if(values)]
        cols = st.columns(2)
        for i, pp in enumerate(visible):
            with cols[i % 2]:
                values[pp.name] = _widget(step.key, pp)
        # second pass: params unlocked by a value just set (e.g. organism == other)
        for pp in [pp for pp in plist if pp.show_if and pp.show_if(values) and pp not in visible]:
            with cols[len(visible) % 2]:
                values[pp.name] = _widget(step.key, pp)

    if common:
        _render(common)
    if advanced:
        with st.expander("⚙ Advanced parameters"):
            _render(advanced)
    return values


# --------------------------------------------------------------------------- #
# output rendering (data stays in-block; messages go to the popup)
# --------------------------------------------------------------------------- #
def _render_data(out: Dict[str, Any]):
    if not out:
        return
    if out.get("df") is not None:
        st.dataframe(out["df"], use_container_width=True, height=240)
    if out.get("image") and os.path.exists(out["image"]):
        st.image(out["image"], use_container_width=True)
    if out.get("html") and os.path.exists(out["html"]):
        with open(out["html"]) as fh:
            st.download_button("⬇ download HTML report", fh.read(),
                               file_name=os.path.basename(out["html"]),
                               mime="text/html", key=f"dl_{out['html']}")
    if out.get("log"):
        with st.expander("console output"):
            st.code(out["log"].strip() or "(no output)")


def _execute(step: StepSpec, params: Dict[str, Any]):
    """Run a step, store its data output, and queue its message for the popup."""
    missing = [n for n in step.needs if SS.get(n) is None]
    if missing:
        SS[f"status_{step.key}"] = "⏭️"
        _note(step.label, "info", f"Needs: {', '.join(missing)} — run earlier steps first.")
        return
    try:
        out = step.run(SS, params)
        SS[f"out_{step.key}"] = out
        level = out.get("level", "success")
        SS[f"status_{step.key}"] = {"warning": "⚠️", "info": "⏭️", "error": "❌"}.get(level, "✅")
        if out.get("text"):
            _note(step.label, level, out["text"])
    except Exception as e:                                  # noqa: BLE001
        SS[f"status_{step.key}"] = "❌"
        SS[f"out_{step.key}"] = {}
        _note(step.label, "error", f"{type(e).__name__}: {e}")


# --------------------------------------------------------------------------- #
# one step block
# --------------------------------------------------------------------------- #
def _render_step(step: StepSpec, run_all: bool):
    with st.container(border=True):
        head = st.columns([0.72, 0.10, 0.18])
        head[0].markdown(f"**{step.label}**", help=step.help)
        head[1].markdown(SS.get(f"status_{step.key}", " "))
        run_clicked = head[2].button("▶ Run", key=f"run_{step.key}", use_container_width=True)

        params = _collect_params(step)
        _maybe_uploader(step, params)

        enabled = True if step.group == "core" and step.default_on else bool(SS.get(f"enable_{step.key}"))
        if run_clicked or (run_all and enabled):
            with st.spinner(f"Running {step.label}…"):
                _execute(step, params)

        _render_data(SS.get(f"out_{step.key}"))

        if step.key == "bulk_hm" and SS.get("state") is not None:
            _render_curation()


def _maybe_uploader(step: StepSpec, params: Dict[str, Any]):
    field = _UPLOAD_FIELD.get(step.key)
    if not field:
        return
    up = st.file_uploader(f"…or upload the file for `{field}`", key=f"up_{step.key}",
                          help="For very large files, prefer the path box above "
                               "(read directly from disk, no upload size limit).")
    if up is not None:
        path = os.path.join(tempfile.mkdtemp(), up.name)
        with open(path, "wb") as fh:
            fh.write(up.getbuffer())
        params[field] = path
        st.caption(f"using uploaded file: {up.name}")


# --------------------------------------------------------------------------- #
# module curation — folded into the heatmap block (text-based, reliable)
# --------------------------------------------------------------------------- #
def _do_curation(fn, *a):
    out = fn(SS, *a)
    SS["out_bulk_hm"] = {**SS.get("out_bulk_hm", {}), **{k: out[k] for k in ("image",) if k in out}}
    _note("Curate modules", out.get("level", "success"), out.get("text", ""))


def _render_curation():
    state = SS["state"]
    blocks = ordered_module_blocks(state)
    with st.expander(f"✏️ Curate modules ({len(blocks)} modules — edits update the heatmap above)"):
        st.caption("Merge only adjacent modules; split re-clusters a module into two; "
                   "reassign moves named genes. Press **💾 Save** to write the updated "
                   "module_state.pkl + HM_ModGene.csv + heatmap to the output folder.")

        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Merge two neighbours**")
            mm = st.multiselect("modules to merge", blocks, max_selections=2,
                                key="cur_merge", label_visibility="collapsed")
            if st.button("Merge", use_container_width=True, key="btn_merge"):
                if len(mm) == 2:
                    _do_curation(curate_merge, mm[0], mm[1]); st.rerun()
                else:
                    _note("Curate modules", "info", "Pick exactly two neighbouring modules.")
        with c2:
            st.markdown("**Split a module**")
            sm = st.selectbox("module to split", blocks, key="cur_split",
                              label_visibility="collapsed")
            if st.button("Split", use_container_width=True, key="btn_split"):
                _do_curation(curate_split, sm); st.rerun()

        st.markdown("**Reassign genes**")
        c3, c4, c5 = st.columns([0.55, 0.2, 0.25])
        genes_txt = c3.text_input("genes (comma-separated)", key="cur_genes",
                                  placeholder="Sox2, Gapdh", label_visibility="collapsed")
        target = c4.number_input("→ module (0=unassign)", min_value=0,
                                 max_value=max(blocks or [1]), value=0, key="cur_target",
                                 label_visibility="collapsed")
        if c5.button("Reassign", use_container_width=True, key="btn_reassign"):
            genes = [g.strip() for g in genes_txt.split(",") if g.strip()]
            _do_curation(curate_reassign, genes, int(target)); st.rerun()

        st.markdown("**Find a gene's module**")
        c6, c7 = st.columns([0.7, 0.3])
        q = c6.text_input("search gene", key="cur_search", placeholder="gene name",
                          label_visibility="collapsed")
        if c7.button("Find", use_container_width=True, key="btn_find") and q:
            m = state.gene_module(q)
            _note("Find gene", "info" if m is None else "success",
                  f"'{q}' not found." if m is None else f"'{q}' is in module {m}.")

        st.divider()
        if st.button("💾 Save modules + refresh outputs", type="primary",
                     use_container_width=True, key="btn_save"):
            _do_curation(save_curation); st.rerun()


# --------------------------------------------------------------------------- #
# top control panel
# --------------------------------------------------------------------------- #
def _top_controls() -> bool:
    with st.container(border=True):
        st.markdown("**Run the whole pipeline**")
        top = st.columns([0.34, 0.4, 0.26])
        SS["outdir"] = top[0].text_input("Output folder", value=SS.get("outdir", "clusmap_out"))
        mode = top[1].radio("Annotation mode",
                            ["Markers (PanglaoDB)", "Single-cell pseudo-bulk", "None"],
                            horizontal=True,
                            help="How to characterise modules in a one-click run. "
                                 "You can still run either block manually below.")
        if SS.get("_annot_mode") != mode:
            SS["enable_markers"] = (mode == "Markers (PanglaoDB)")
            SS["enable_pseudobulk"] = (mode == "Single-cell pseudo-bulk")
            SS["_annot_mode"] = mode

        if mode == "Single-cell pseudo-bulk":
            pc = st.columns([0.7, 0.3])
            SS["_pb_h5ad"] = pc[0].text_input("Single-cell .h5ad", value=SS.get("_pb_h5ad", ""),
                                              placeholder="path/to/atlas.h5ad",
                                              help="Imported here so a one-click run can build pseudo-bulk.")
            SS["_pb_celltype"] = pc[1].text_input("cell-type obs key",
                                                  value=SS.get("_pb_celltype", "celltype"))

        st.caption("Include in the one-click run (core steps 1–4 always run):")
        opt = [s for s in STEPS if s.group in ("annotate", "downstream") or
               (s.group == "core" and not s.default_on)]
        cols = st.columns(4)
        for i, s in enumerate(opt):
            with cols[i % 4]:
                if s.key in ("markers", "pseudobulk"):
                    st.toggle(s.label.split(" · ")[0], value=bool(SS.get(f"enable_{s.key}")),
                              key=f"top_{s.key}", disabled=True,
                              help="Controlled by the Annotation mode above.")
                else:
                    SS[f"enable_{s.key}"] = st.toggle(
                        s.label.split(" · ")[0], value=bool(SS.get(f"enable_{s.key}", s.default_on)),
                        key=f"top_{s.key}", help=s.help)
        return st.button("⚡ Run pipeline", type="primary", use_container_width=True)


# --------------------------------------------------------------------------- #
# page
# --------------------------------------------------------------------------- #
def render_pipeline():
    # notes are collected during a run and shown by _flush_notes(), which clears
    # them. We do NOT reset here, so a note queued just before a curation st.rerun()
    # survives into the next run and still pops up.
    SS.setdefault("_pending_notes", [])
    for s in STEPS:                                # make sure enable flags exist
        SS.setdefault(f"enable_{s.key}", s.default_on)

    st.caption("Each block is a pipeline step. Tune parameters (hover the name for help; "
               "the gray text is the default), run one with **▶ Run**, or set up the top "
               "panel and press **⚡ Run pipeline**. Curate modules under the heatmap block.")

    run_all = _top_controls()

    groups = [("Core pipeline", "core"), ("Module annotation", "annotate"),
              ("Downstream analysis", "downstream")]
    for title, g in groups:
        st.subheader(title)
        for step in [s for s in STEPS if s.group == g]:
            _render_step(step, run_all)

    _flush_notes()
