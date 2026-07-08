"""Interactive command-line chat for the clusmap agent.

Run:
    python -m clusmap.agent                       # uses Claude, ANTHROPIC_API_KEY from env
    python -m clusmap.agent --provider anthropic --api-key sk-...
    python -m clusmap.agent --outdir results

You then just talk to it: "Cluster mouse-thymus.tsv into about 8 modules and
draw the heatmap coloured by organ." It asks for what it needs and runs the
pipeline through the tool layer.
"""
from __future__ import annotations

import argparse
import os
import sys

from .tools import TOOLS, ToolSession
from .backends import make_backend


def chat(provider: str = "anthropic", api_key: str | None = None,
         outdir: str = "clusmap_out", model: str | None = None) -> None:
    session = ToolSession(outdir=outdir)
    os.makedirs(outdir, exist_ok=True)
    kwargs = {"api_key": api_key}
    if model:
        kwargs["model"] = model
    backend = make_backend(provider, TOOLS, **kwargs)

    def on_tool(name, args):
        print(f"  · running {name}({', '.join(f'{k}={v!r}' for k, v in args.items())})",
              file=sys.stderr)
        return session.run(name, args)

    print(f"clusmap agent ({provider}) — outputs in {outdir}/. Type 'exit' to quit.\n")
    while True:
        try:
            user = input("you > ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            break
        if user.lower() in {"exit", "quit", ":q"}:
            break
        if not user:
            continue
        try:
            reply = backend.send(user, on_tool)
        except Exception as e:
            print(f"[agent error] {e}", file=sys.stderr)
            continue
        print(f"\nclusmap > {reply}\n")


def main(argv=None):
    p = argparse.ArgumentParser(description="Conversational clusmap agent.")
    p.add_argument("--provider", default=os.environ.get("CLUSMAP_PROVIDER", "anthropic"),
                   help="LLM provider (anthropic|claude|openai). Default anthropic.")
    p.add_argument("--api-key", default=None,
                   help="API key for the provider. Defaults to the provider's env var "
                        "(e.g. ANTHROPIC_API_KEY).")
    p.add_argument("--model", default=None, help="Override the model id.")
    p.add_argument("--outdir", default="clusmap_out", help="Where outputs are written.")
    args = p.parse_args(argv)
    chat(provider=args.provider, api_key=args.api_key, outdir=args.outdir, model=args.model)


if __name__ == "__main__":
    main()
