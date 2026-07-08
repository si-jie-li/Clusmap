"""Pluggable LLM backends for the clusmap agent.

The agent core only depends on the small ``LLMBackend`` interface: it feeds a
user message in and gets the assistant's reply out, while the backend handles
the provider-specific tool-calling loop and executes tools via the ``on_tool``
callback. Pick your provider by constructing the matching backend with its API
key — the tool layer and runner stay identical.

Claude is implemented here against the official ``anthropic`` SDK. To add
another provider, implement ``send`` for its tool-calling API (see
``OpenAIBackend`` for the contract).
"""
from __future__ import annotations

from typing import Callable, List, Protocol

ToolFn = Callable[[str, dict], str]   # (tool_name, args) -> result string

SYSTEM_PROMPT = (
    "You are the clusmap assistant. You help a biologist run a bulk RNA-seq "
    "module-discovery pipeline by calling tools, gathering the few parameters "
    "you need conversationally, and reporting results.\n\n"
    "Guidelines:\n"
    "- The usual order is: load_data -> preprocess -> cluster -> make_heatmap -> "
    "(celltype_swarm / go_enrichment / edits) -> save_state.\n"
    "- Ask only for parameters you genuinely need; otherwise use sensible defaults. "
    "When the user says e.g. 'about 8 modules', pick a deepSplit, run cluster, report the "
    "count, and offer to adjust rather than interrogating them up front.\n"
    "- After load_data, confirm the matrix looks like genes x samples before proceeding.\n"
    "- Be concise. Report what a tool produced (counts, file paths) and propose the next step."
)


class LLMBackend(Protocol):
    def send(self, user_text: str, on_tool: ToolFn) -> str:
        """Send one user turn, run the tool loop to completion, return reply text."""
        ...


# --------------------------------------------------------------------------- #
# Claude (Anthropic SDK)
# --------------------------------------------------------------------------- #
class AnthropicBackend:
    def __init__(self, tools: List[dict], api_key: str | None = None,
                 model: str = "claude-opus-4-8", max_tokens: int = 8000):
        import anthropic
        # api_key=None lets the SDK resolve ANTHROPIC_API_KEY / profile from env
        self.client = anthropic.Anthropic(api_key=api_key) if api_key else anthropic.Anthropic()
        self.tools = tools
        self.model = model
        self.max_tokens = max_tokens
        self.messages: list = []

    def send(self, user_text: str, on_tool: ToolFn) -> str:
        self.messages.append({"role": "user", "content": user_text})
        while True:
            resp = self.client.messages.create(
                model=self.model,
                max_tokens=self.max_tokens,
                system=SYSTEM_PROMPT,
                tools=self.tools,
                thinking={"type": "adaptive"},
                messages=self.messages,
            )
            self.messages.append({"role": "assistant", "content": resp.content})

            if resp.stop_reason != "tool_use":
                return "".join(b.text for b in resp.content if b.type == "text").strip()

            results = []
            for block in resp.content:
                if block.type == "tool_use":
                    out = on_tool(block.name, block.input or {})
                    results.append({"type": "tool_result", "tool_use_id": block.id,
                                    "content": out})
            self.messages.append({"role": "user", "content": results})


# --------------------------------------------------------------------------- #
# Extension point for other providers
# --------------------------------------------------------------------------- #
class OpenAIBackend:
    """Stub. Implement `send` against your provider's tool-calling API.

    Contract (mirror AnthropicBackend):
      1. Append the user message to an internal history.
      2. Call the chat/completions endpoint with `self.tools` converted to the
         provider's function/tool schema (our schemas are plain JSON Schema, so
         most providers accept them with a thin wrapper, e.g.
         {"type": "function", "function": {"name", "description", "parameters": input_schema}}).
      3. While the model returns tool calls: run each via on_tool(name, json.loads(args))
         and append the results in the provider's tool-result message shape.
      4. When the model returns a normal message, return its text.
    """
    def __init__(self, tools, api_key=None, model="gpt-4o", **_):
        raise NotImplementedError(
            "OpenAIBackend is a documented stub. Default to AnthropicBackend, or "
            "implement send() following the contract in this class's docstring."
        )

    def send(self, user_text: str, on_tool: ToolFn) -> str:  # pragma: no cover
        raise NotImplementedError


BACKENDS = {"anthropic": AnthropicBackend, "claude": AnthropicBackend,
            "openai": OpenAIBackend}


def make_backend(provider: str, tools: List[dict], **kwargs) -> LLMBackend:
    provider = (provider or "anthropic").lower()
    if provider not in BACKENDS:
        raise ValueError(f"Unknown provider {provider!r}. Options: {sorted(BACKENDS)}")
    return BACKENDS[provider](tools, **kwargs)
