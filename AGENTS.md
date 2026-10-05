# OpenCode repository instructions

Read `CLAUDE.md` in full before working in this repository and follow its
operating manual, lifecycle, boundaries, conventions, and verification rules.
It is the canonical repository guidance for all coding assistants.

Apply its **Search and edit loop** to every coding/review task and specialist
handoff: resolve indexed anchors against current source and check index freshness
at edit checkpoints.

OpenCode discovers the existing skills in `.claude/skills/` automatically.
Native specialist agents and slash commands live in `.opencode/agents/` and
`.opencode/commands/`; project permissions live in `opencode.json`.
See [.opencode/README.md](.opencode/README.md) for usage.
