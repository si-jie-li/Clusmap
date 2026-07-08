"""ModuleState: the single source of truth for clustering results.

Every downstream step (heatmap, GO, motif, pseudo-bulk) and every interactive
edit (split / merge / reassign / search) reads and writes this object, so module
identity stays consistent across the whole pipeline and across sessions.

Numbering convention
--------------------
* ``raw_labels``  : labels straight from ``cutreeHybrid`` (0 = unassigned).
* ``hm_labels``   : 1-based labels renumbered in *heatmap (dendrogram) order*,
  which is what the figures and ``HM_ModGene.csv`` use. 0 stays 0 (unassigned).

The heatmap row order is ``scipy...leaves_list(linkage)`` so it is fully
reproducible and independent of any drawing step.
"""
from __future__ import annotations

import json
import os
import pickle
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage, fcluster
from scipy.spatial.distance import pdist


@dataclass
class ModuleState:
    genes: np.ndarray                       # gene names, in original rna_df order
    linkage: np.ndarray                     # scipy linkage matrix
    raw_labels: np.ndarray                  # cutreeHybrid labels, original order
    metric: str = "correlation"
    method: str = "average"
    history: List[str] = field(default_factory=list)
    _data: Optional[np.ndarray] = None      # optional genes x samples values (for split)

    # ----------------------------------------------------------------- build #
    @classmethod
    def from_cutree(cls, link, mod, genes, *, metric="correlation",
                    method="average", data: Optional[pd.DataFrame] = None):
        return cls(
            genes=np.asarray(list(genes), dtype=object),
            linkage=np.asarray(link),
            raw_labels=np.asarray(mod["labels"] if isinstance(mod, dict) else mod),
            metric=metric, method=method,
            _data=None if data is None else np.asarray(data.values),
        )

    # --------------------------------------------------------------- queries #
    @property
    def order(self) -> np.ndarray:
        """Row order (positions into ``genes``) as drawn in the heatmap."""
        return leaves_list(self.linkage)

    @property
    def ordered_genes(self) -> np.ndarray:
        return self.genes[self.order]

    def _renumber(self) -> Dict[int, int]:
        """Map raw label -> 1-based hm id following first appearance in order."""
        seen, mapping, nxt = set(), {0: 0}, 0
        for lab in self.raw_labels[self.order]:
            if lab == 0 or lab in seen:
                continue
            seen.add(lab)
            nxt += 1
            mapping[lab] = nxt
        return mapping

    @property
    def hm_labels(self) -> np.ndarray:
        """1-based heatmap labels, aligned to original gene order."""
        mapping = self._renumber()
        return np.array([mapping[l] for l in self.raw_labels])

    def to_modgene_df(self, colors: Optional[Sequence[str]] = None) -> pd.DataFrame:
        """Emit the HM_ModGene table (hm_mod, gene[, hm_color]) in heatmap order."""
        order = self.order
        df = pd.DataFrame({
            "hm_mod": self.hm_labels[order],
            "gene": self.genes[order],
        })
        if colors is not None:
            df["hm_color"] = np.asarray(colors)[order]
        return df.sort_values(["hm_mod", "gene"]).reset_index(drop=True)

    def gene_module(self, gene: str) -> Optional[int]:
        """Heatmap module id a gene belongs to (None if absent, 0 = unassigned)."""
        idx = np.where(self.genes == gene)[0]
        if len(idx) == 0:  # case-insensitive fallback
            low = {str(g).lower(): i for i, g in enumerate(self.genes)}
            if str(gene).lower() not in low:
                return None
            idx = [low[str(gene).lower()]]
        return int(self.hm_labels[idx[0]])

    def module_genes(self, hm_id: int) -> List[str]:
        """Genes in a given heatmap module id, in heatmap order."""
        order = self.order
        labs = self.hm_labels[order]
        return list(self.genes[order][labs == hm_id])

    @property
    def n_modules(self) -> int:
        return len({l for l in self.raw_labels if l != 0})

    # ------------------------------------------------------------- mutations #
    def merge(self, hm_id_a: int, hm_id_b: int) -> "ModuleState":
        """Merge two modules (by heatmap id) into one. In place + returns self."""
        mapping = self._renumber()
        inv = {v: k for k, v in mapping.items()}
        raw_a, raw_b = inv[hm_id_a], inv[hm_id_b]
        self.raw_labels[self.raw_labels == raw_b] = raw_a
        self.history.append(f"merge {hm_id_a}+{hm_id_b}")
        return self

    def split(self, hm_id: int) -> "ModuleState":
        """Split a module into two by re-clustering its genes (needs data)."""
        if self._data is None:
            raise ValueError("split() needs expression data; build the state with "
                             "ModuleState.from_cutree(..., data=rna_df).")
        mapping = self._renumber()
        inv = {v: k for k, v in mapping.items()}
        raw = inv[hm_id]
        mask = self.raw_labels == raw
        if mask.sum() < 4:
            raise ValueError(f"Module {hm_id} too small to split ({mask.sum()} genes).")
        sub = self._data[mask]
        sub_link = linkage(pdist(sub, metric=self.metric), method=self.method)
        parts = fcluster(sub_link, t=2, criterion="maxclust")
        new_raw = (self.raw_labels.max() + 1)
        idx = np.where(mask)[0]
        self.raw_labels[idx[parts == 2]] = new_raw  # second part becomes a new module
        self.history.append(f"split {hm_id}")
        return self

    def reassign(self, genes: Sequence[str], hm_id: int) -> "ModuleState":
        """Move a set of genes into an existing module (or 0 to unassign)."""
        mapping = self._renumber()
        inv = {v: k for k, v in mapping.items()}
        raw = inv.get(hm_id, 0) if hm_id != 0 else 0
        low = {str(g).lower(): i for i, g in enumerate(self.genes)}
        for g in genes:
            i = low.get(str(g).lower())
            if i is not None:
                self.raw_labels[i] = raw
        self.history.append(f"reassign {len(genes)} genes -> {hm_id}")
        return self

    # ----------------------------------------------------------------- io    #
    def save(self, outdir: str) -> None:
        os.makedirs(outdir, exist_ok=True)
        with open(os.path.join(outdir, "module_state.pkl"), "wb") as fh:
            pickle.dump(self, fh)
        meta = {"n_modules": self.n_modules, "metric": self.metric,
                "method": self.method, "history": self.history}
        with open(os.path.join(outdir, "module_state.json"), "w") as fh:
            json.dump(meta, fh, indent=2)
        self.to_modgene_df().to_csv(os.path.join(outdir, "HM_ModGene.csv"), index=False)

    @staticmethod
    def load(path: str) -> "ModuleState":
        if os.path.isdir(path):
            path = os.path.join(path, "module_state.pkl")
        with open(path, "rb") as fh:
            return pickle.load(fh)
