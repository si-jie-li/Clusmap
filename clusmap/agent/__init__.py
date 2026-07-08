"""Conversational, provider-agnostic agent front-end for clusmap.

    from clusmap.agent import chat
    chat(provider="anthropic")          # talk to it in the terminal

Or from the shell:  python -m clusmap.agent
"""
from .tools import TOOLS, ToolSession
from .backends import LLMBackend, AnthropicBackend, make_backend
from .runner import chat, main

__all__ = ["TOOLS", "ToolSession", "LLMBackend", "AnthropicBackend",
           "make_backend", "chat", "main"]
