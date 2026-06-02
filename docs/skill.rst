LLM Reference (SKILL.md)
========================

``SKILL.md`` is a structured reference document designed to be loaded into an
LLM-powered coding assistant (Claude Code, ChatGPT, Cursor, etc.) to give it
accurate, up-to-date knowledge of the One Codex Python client. Without it,
assistants may hallucinate incorrect method names, unsupported filter
operators, or non-existent API patterns.

| :download:`Download SKILL.md <SKILL.md>`
| **Raw URL (stable):** ``https://onecodex.github.io/onecodex/SKILL.md``

.. rubric:: Usage

**Claude Code** — place this file in your project root. Claude Code will
automatically read files named ``SKILL.md`` in the working directory.

**Claude.ai / ChatGPT / other chat interfaces** — paste the contents directly
into the conversation, or attach the file at the start of a session before
asking questions.

**OpenAI Codex CLI** — pass the file as part of your instructions file or
prepend it to your prompt.

.. rubric:: Contents

.. include:: SKILL.md
   :parser: myst_parser.sphinx_
