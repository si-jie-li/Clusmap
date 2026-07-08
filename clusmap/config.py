"""User configuration for paths that must NOT ship inside the package.

The genome FASTA and motif databases (JASPAR / HOCOMOCO / CIS-BP) are multi-GB
files that live once on the user's machine or HPC cluster. clusmap only stores
*paths* to them, plus how to run the motif jobs (locally, via SLURM, or over SSH
to a cluster).

Resolution order (highest priority first):
1. explicit argument passed in code,
2. environment variable ``CLUSMAP_<KEY>`` (e.g. ``CLUSMAP_GENOME_FASTA``),
3. ``~/.clusmap/config.yaml``,
4. built-in default (usually ``None``).

Set values from the shell::

    clusmap-config set genome_fasta /data/mm10.fa
    clusmap-config set motif_db.JASPAR /db/JASPAR2026_vertebrates.meme
    clusmap-config show
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Optional

CONFIG_DIR = Path(os.environ.get("CLUSMAP_CONFIG_DIR", Path.home() / ".clusmap"))
CONFIG_PATH = CONFIG_DIR / "config.yaml"

DEFAULTS: Dict[str, Any] = {
    "genome_fasta": None,            # e.g. /data/genomes/mm10.fa
    "tss_bed": None,                 # precomputed TSS bed (optional)
    "promoter_fasta": None,          # precomputed promoter FASTA (gene-named headers)
    "promoter_bp": 500,              # +/- bp around TSS when extracting promoters
    "motif_db": {},                  # {name: /path/to/db.meme}
    "run_mode": "local",             # local | slurm | ssh
    "ame_evalue": 1.0,
    "min_genes_for_ame": 50,
    # HPC / SSH
    "ssh_host": None,                # user@cluster
    "remote_workdir": None,          # scratch dir on the cluster
    "slurm_partition": None,
    "slurm_cpus": 1,
    "slurm_time": "04:00:00",
}


def _load_file() -> Dict[str, Any]:
    if not CONFIG_PATH.exists():
        return {}
    import yaml
    with open(CONFIG_PATH) as fh:
        return yaml.safe_load(fh) or {}


def load() -> Dict[str, Any]:
    """Merged config dict (defaults <- file <- environment)."""
    cfg = {**DEFAULTS, **_load_file()}
    for key in cfg:
        env = os.environ.get(f"CLUSMAP_{key.upper()}")
        if env is not None:
            cfg[key] = env
    return cfg


def get(key: str, default: Any = None) -> Any:
    """Get one config value, honouring a nested ``motif_db.NAME`` key."""
    cfg = load()
    if "." in key:
        head, tail = key.split(".", 1)
        return (cfg.get(head) or {}).get(tail, default)
    return cfg.get(key, default)


def set_(key: str, value: Any) -> None:
    """Persist one value to ``~/.clusmap/config.yaml`` (supports ``motif_db.NAME``)."""
    import yaml
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    cfg = _load_file()
    if "." in key:
        head, tail = key.split(".", 1)
        cfg.setdefault(head, {})
        cfg[head][tail] = _coerce(value)
    else:
        cfg[key] = _coerce(value)
    with open(CONFIG_PATH, "w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=True)


def _coerce(v: Any):
    if not isinstance(v, str):
        return v
    for cast in (int, float):
        try:
            return cast(v)
        except ValueError:
            pass
    if v.lower() in {"true", "false"}:
        return v.lower() == "true"
    if v.lower() in {"null", "none", ""}:
        return None
    return v


def builtin_motif(organism: str):
    """Resolve (promoter_fasta, [motif_dbs]) for a built-in organism.

    Looks first in the packaged reference dir ``clusmap/data/motif/<organism>/``
    (``promoter.fa`` + ``*.meme``) — this is where human/mouse references will be
    bundled at publish time — then falls back to user config keys
    (``promoter_<organism>`` / ``promoter_fasta`` and ``motif_db``).
    Returns ``(None, [])`` if nothing is configured.
    """
    org = (organism or "").lower()
    base = Path(__file__).parent / "data" / "motif" / org
    prom = base / "promoter.fa"
    dbs = sorted(str(p) for p in base.glob("*.meme")) if base.exists() else []
    if prom.exists() and dbs:
        return str(prom), dbs
    # fall back to user-configured paths
    cprom = get(f"promoter_{org}") or get("promoter_fasta")
    cdbs = list((get("motif_db") or {}).values())
    return cprom, cdbs


def require(*keys: str) -> Dict[str, Any]:
    """Return the values for ``keys``, raising a clear error if any is unset."""
    cfg = load()
    missing = [k for k in keys if not get(k)]
    if missing:
        raise RuntimeError(
            "Missing clusmap config: " + ", ".join(missing) +
            f".\nSet them with `clusmap-config set <key> <value>` "
            f"(config file: {CONFIG_PATH})."
        )
    return {k: get(k) for k in keys}


# --------------------------------------------------------------------------- #
# tiny CLI: clusmap-config show | get KEY | set KEY VALUE
# --------------------------------------------------------------------------- #
def _cli(argv: Optional[list] = None) -> None:
    import argparse
    import json
    p = argparse.ArgumentParser(prog="clusmap-config", description="clusmap path/HPC config")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("show")
    g = sub.add_parser("get"); g.add_argument("key")
    s = sub.add_parser("set"); s.add_argument("key"); s.add_argument("value")
    args = p.parse_args(argv)

    if args.cmd == "show":
        print(json.dumps(load(), indent=2))
        print(f"\n(config file: {CONFIG_PATH})")
    elif args.cmd == "get":
        print(get(args.key))
    elif args.cmd == "set":
        set_(args.key, args.value)
        print(f"set {args.key} = {get(args.key)}  -> {CONFIG_PATH}")


if __name__ == "__main__":
    _cli()
