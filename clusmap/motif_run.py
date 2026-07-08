"""Run the motif-enrichment (AME) pipeline from Python, wired to HM_ModGene.csv.

The heavy reference data (genome FASTA, JASPAR/HOCOMOCO/CIS-BP motif DBs) lives
on your machine or HPC cluster and is referenced through :mod:`clusmap.config`.
This module:

1. splits ``HM_ModGene.csv`` into per-module gene lists,
2. subsets a gene-named promoter FASTA into one FASTA per module (pure Python —
   no seqkit dependency),
3. runs ``ame`` for each module locally, via ``sbatch`` (SLURM), or over SSH to
   a cluster, and
4. collects the per-module ``ame.tsv`` outputs with
   :func:`clusmap.motif.module_motif`.

The bundled scripts in ``motif_analysis/`` (genome/promoter prep + the original
SLURM batch) remain the source of truth for one-time genome preparation; this
module is the per-run driver.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from typing import Dict, List, Optional

import pandas as pd

from . import config
from .motif import module_motif


# --------------------------------------------------------------------------- #
# FASTA handling (pure python)
# --------------------------------------------------------------------------- #
def _gene_of(header: str) -> str:
    """Gene symbol from a promoter FASTA header like 'GENE::chr1:100-600(+)'."""
    h = header[1:].strip() if header.startswith(">") else header.strip()
    return h.split("::", 1)[0]


def _index_promoters(promoter_fasta: str) -> Dict[str, List[str]]:
    """Map upper-cased gene symbol -> list of FASTA records (header+seq blocks)."""
    index: Dict[str, List[str]] = {}
    header, seq = None, []
    with open(promoter_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    index.setdefault(_gene_of(header).upper(), []).append(
                        header + "".join(seq))
                header, seq = line, []
            else:
                seq.append(line)
    if header is not None:
        index.setdefault(_gene_of(header).upper(), []).append(header + "".join(seq))
    return index


def prepare_module_fastas(
    modgene_path: str,
    promoter_fasta: str,
    outdir: str,
    *,
    module_col: str = "hm_mod",
    gene_col: str = "gene",
    min_genes: int = 50,
    skip_module_0: bool = True,
) -> Dict[int, str]:
    """Write one promoter FASTA per module; return {module: fasta_path}.

    Modules with fewer than ``min_genes`` matched promoters are skipped (AME is
    unreliable on tiny gene sets), matching the original ``batch_AME.sh`` policy.
    """
    df = pd.read_csv(modgene_path)
    index = _index_promoters(promoter_fasta)
    fasta_dir = os.path.join(outdir, "mod_fasta")
    os.makedirs(fasta_dir, exist_ok=True)

    written: Dict[int, str] = {}
    for module, sub in df.groupby(module_col, sort=True):
        m = int(module)
        if skip_module_0 and m == 0:
            continue
        records = []
        for g in sub[gene_col].astype(str):
            records.extend(index.get(g.upper(), []))
        if len(records) < min_genes:
            print(f"[skip] module {m}: only {len(records)} promoters (< {min_genes}).")
            continue
        path = os.path.join(fasta_dir, f"Module{m}.fa")
        with open(path, "w") as fh:
            fh.write("".join(records))
        written[m] = path
        print(f"[ok] module {m}: {len(records)} promoters -> {path}")
    return written


# --------------------------------------------------------------------------- #
# AME invocation
# --------------------------------------------------------------------------- #
def _ame_cmd(input_fasta: str, out_dir: str, motif_dbs: List[str], evalue: float) -> List[str]:
    return [
        "ame", "--oc", out_dir, "--control", "--shuffle--", "--kmer", "2",
        "--method", "fisher", "--scoring", "totalhits",
        "--evalue-report-threshold", str(evalue), "--noseq",
        input_fasta, *motif_dbs,
    ]


def _slurm_script(module: int, input_fasta: str, out_dir: str, motif_dbs: List[str],
                  evalue: float, cfg: dict) -> str:
    cmd = " ".join(_ame_cmd(input_fasta, out_dir, motif_dbs, evalue))
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=AME_m{module}",
        "#SBATCH --nodes=1",
        f"#SBATCH --cpus-per-task={cfg.get('slurm_cpus', 1)}",
        f"#SBATCH --time={cfg.get('slurm_time', '04:00:00')}",
    ]
    if cfg.get("slurm_partition"):
        lines.append(f"#SBATCH --partition={cfg['slurm_partition']}")
    lines += ["set -euo pipefail", f'mkdir -p "{out_dir}"', cmd, ""]
    return "\n".join(lines)


def run_ame(
    module_fastas: Dict[int, str],
    outdir: str,
    *,
    motif_dbs: Optional[List[str]] = None,
    run_mode: Optional[str] = None,
    evalue: Optional[float] = None,
    dry_run: bool = False,
) -> str:
    """Run AME for each module. Returns the AME output root (``outdir/AME_out``).

    ``run_mode`` is ``local`` (needs ``ame`` on PATH), ``slurm`` (submits one
    ``sbatch`` job per module), or ``ssh`` (rsyncs FASTAs to the cluster and
    submits there). Defaults come from :mod:`clusmap.config`.
    """
    cfg = config.load()
    run_mode = run_mode or cfg["run_mode"]
    evalue = cfg["ame_evalue"] if evalue is None else evalue
    if motif_dbs is None:
        motif_dbs = list((cfg.get("motif_db") or {}).values())
    if not motif_dbs:
        raise RuntimeError("No motif databases configured. "
                           "`clusmap-config set motif_db.JASPAR /path/to/db.meme`")

    ame_root = os.path.join(outdir, "AME_out")
    os.makedirs(ame_root, exist_ok=True)

    if run_mode == "local":
        if not dry_run and shutil.which("ame") is None:
            raise RuntimeError("`ame` (MEME-suite) not found on PATH. Install MEME, "
                               "use the clusmap Docker image, or run_mode='slurm'/'ssh'.")
        for m, fa in module_fastas.items():
            out_dir = os.path.join(ame_root, f"Module{m}")
            cmd = _ame_cmd(fa, out_dir, motif_dbs, evalue)
            print("[ame]", " ".join(cmd))
            if not dry_run:
                subprocess.run(cmd, check=True)

    elif run_mode == "slurm":
        script_dir = os.path.join(outdir, "slurm")
        os.makedirs(script_dir, exist_ok=True)
        for m, fa in module_fastas.items():
            out_dir = os.path.join(ame_root, f"Module{m}")
            script = os.path.join(script_dir, f"ame_module{m}.slurm")
            with open(script, "w") as fh:
                fh.write(_slurm_script(m, fa, out_dir, motif_dbs, evalue, cfg))
            print("[sbatch]", script)
            if not dry_run:
                subprocess.run(["sbatch", script], check=True)

    elif run_mode == "ssh":
        _run_ssh(module_fastas, ame_root, motif_dbs, evalue, cfg, dry_run)

    else:
        raise ValueError(f"Unknown run_mode {run_mode!r} (local|slurm|ssh).")

    return ame_root


def _run_ssh(module_fastas, ame_root, motif_dbs, evalue, cfg, dry_run):
    host = cfg.get("ssh_host")
    remote = cfg.get("remote_workdir")
    if not host or not remote:
        raise RuntimeError("run_mode='ssh' needs ssh_host and remote_workdir in config.")
    remote_fasta = f"{remote}/mod_fasta"
    if not dry_run:
        subprocess.run(["ssh", host, f"mkdir -p {remote_fasta} {remote}/AME_out"], check=True)
    for m, fa in module_fastas.items():
        dst = f"{host}:{remote_fasta}/Module{m}.fa"
        print("[scp]", fa, "->", dst)
        if not dry_run:
            subprocess.run(["scp", fa, dst], check=True)
        out_dir = f"{remote}/AME_out/Module{m}"
        cmd = " ".join(_ame_cmd(f"{remote_fasta}/Module{m}.fa", out_dir, motif_dbs, evalue))
        sbatch = (f"sbatch --job-name=AME_m{m} --cpus-per-task={cfg.get('slurm_cpus',1)} "
                  f"--time={cfg.get('slurm_time','04:00:00')} "
                  f"--wrap='mkdir -p {out_dir}; {cmd}'")
        print("[ssh-sbatch]", sbatch)
        if not dry_run:
            subprocess.run(["ssh", host, sbatch], check=True)
    print(f"[ssh] jobs submitted on {host}. After they finish, fetch results:\n"
          f"      rsync -av {host}:{remote}/AME_out/ {ame_root}/")


# --------------------------------------------------------------------------- #
# end-to-end driver
# --------------------------------------------------------------------------- #
def motif_pipeline(
    modgene_path: str,
    outdir: str = "motif_out",
    *,
    promoter_fasta: Optional[str] = None,
    motif_dbs: Optional[List[str]] = None,
    run_mode: Optional[str] = None,
    min_genes: Optional[int] = None,
    evalue: Optional[float] = None,
    wrap_mode: str = "top_n",
    dry_run: bool = False,
):
    """Full motif pipeline: split -> per-module FASTA -> AME -> wrap-up.

    For ``local`` run mode this returns the parsed result tables. For
    ``slurm``/``ssh`` the jobs are submitted asynchronously; re-run
    :func:`clusmap.motif.module_motif` on the AME output once they finish.
    """
    cfg = config.load()
    promoter_fasta = promoter_fasta or cfg.get("promoter_fasta")
    if not promoter_fasta:
        raise RuntimeError("No promoter FASTA. Pass promoter_fasta=... or "
                           "`clusmap-config set promoter_fasta /path/to/promoters.fa` "
                           "(generate it with motif_analysis/promoter_extract.sh + "
                           "promoter_seq_extract.slurm).")
    min_genes = cfg["min_genes_for_ame"] if min_genes is None else min_genes

    os.makedirs(outdir, exist_ok=True)
    module_fastas = prepare_module_fastas(modgene_path, promoter_fasta, outdir,
                                          min_genes=min_genes)
    if not module_fastas:
        raise RuntimeError("No modules had enough promoters to run AME.")

    ame_root = run_ame(module_fastas, outdir, motif_dbs=motif_dbs,
                       run_mode=run_mode, evalue=evalue, dry_run=dry_run)

    mode = run_mode or cfg["run_mode"]
    if mode == "local" and not dry_run:
        return module_motif(ame_root, os.path.join(outdir, "motif_results"), mode=wrap_mode)
    print(f"\nJobs dispatched ({mode}). When AME finishes, wrap up with:\n"
          f"  python -c \"from clusmap.motif import module_motif; "
          f"module_motif('{ame_root}', '{outdir}/motif_results', mode='{wrap_mode}')\"")
    return ame_root
