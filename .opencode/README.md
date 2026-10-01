# OpenCode project setup

The repository supports OpenCode V2 with:

- `opencode.json`: the routine shell-command allowlist adapted from
  `.claude/settings.json`. Other shell commands ask for approval.
- `.opencode/agents/`: the six specialist agents adapted from `.claude/agents/`.
  Ask the primary agent to use a specialist by name, for example:
  `Use physics-reviewer to review this method and its validation evidence.`
- `.opencode/commands/`: `/review-ready` and `/tidy-history`, adapted from the
  Claude commands and using the same tracked procedures.
- `AGENTS.md`: an automatically loaded adapter directing OpenCode to the
  canonical operating manual in `CLAUDE.md`.
- `.claude/skills/`: the existing skills, discovered automatically by OpenCode.

Agents inherit the session's model. Read-only specialists cannot use edit tools;
the three structured-report agents return their JSON for the caller to save.
Shell access for specialists that need it asks for approval, except the Rust
engineer's routine development commands. Shell approval does not enforce a
read-only filesystem; only approve inspection commands for read-only agents.

Maintain the Claude and OpenCode specialist prompts together when changing their
roles. Shared workflow instructions belong in `.claude/skills/` or `CLAUDE.md`.
