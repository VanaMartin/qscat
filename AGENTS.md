# OpenCode repository instructions

Read `CLAUDE.md` in full before working in this repository and follow its
operating manual, lifecycle, boundaries, conventions, and verification rules.
It is the canonical repository guidance for all coding assistants.

Apply its **Search and edit loop** to every coding/review task and specialist
handoff: resolve committed-main anchors against current source, inspect local
changes relative to the indexed commit, and check main-index status at handoff.

OpenCode discovers the existing skills in `.claude/skills/` automatically.
Native specialist agents and slash commands live in `.opencode/agents/` and
`.opencode/commands/`; project permissions live in `opencode.json`.
See [.opencode/README.md](.opencode/README.md) for usage.
