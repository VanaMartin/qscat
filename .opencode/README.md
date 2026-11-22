# OpenCode project setup

The repository supports OpenCode V2 with:

- `opencode.json`: the routine shell-command allowlist adapted from
  `.claude/settings.json` and specialist lookup permissions. Search MCPs are owned
  by the native plugin, not duplicate process definitions. Other shell commands
  ask for approval.
- `.opencode/agents/`: the six specialist agents adapted from `.claude/agents/`.
  Ask the primary agent to use a specialist by name, for example:
  `Use physics-reviewer to review this method and its validation evidence.`
- `.opencode/commands/`: `/review-ready` and `/tidy-history`, adapted from the
  Claude commands and using the same tracked procedures.
- `AGENTS.md`: an automatically loaded adapter directing OpenCode to the
  canonical operating manual in `CLAUDE.md`.
- `.claude/skills/`: the existing skills, discovered automatically by OpenCode.

Agents inherit the session's model. Read-only specialists cannot use edit tools;
the three structured-report agents return only JSON for the caller to check and
save. The caller derives counts and keeps unresolved symbol decisions separate.
Shell access for specialists that need it asks for approval, except the Rust
engineer's routine development commands. Shell approval does not enforce a
read-only filesystem; only approve inspection commands for read-only agents.

Maintain the Claude and OpenCode specialist prompts together when changing their
roles. Shared workflow instructions belong in `.claude/skills/` or `CLAUDE.md`.

## Native repository search

Install [opencode-repo-search](https://github.com/VanaMartin/opencode-repo-search)
once per host if it is not already configured:

```sh
opencode plugin add 'github:VanaMartin/opencode-repo-search#latest'
```

Replace an existing target rather than adding a duplicate. The floating Git tag
can include release candidates; immutable version/full-commit pins are optional.
Reload the location after installing. The package owns its locked isolated runtime
outside the numerical workspace; this repository contains no search engine or
search dependency environment to install.

`search-package/runtime.json` and `search-package/profiles.json` supply native
defaults with two independent corpora: `lancedb` for code and `lancedb-articles`
for committed literature notes. Source remains fetched `origin/main`; the default
local policy is separately pinned by hash. Repository-private settings in Git
common metadata can override these tracked templates and be inherited by linked
worktrees. Read [`docs/agent-search.md`](../docs/agent-search.md) before changing
policy, scope or index identity.

Use `repo-knowledge-search`, or the existing `knowledge-search` guidance adapter.
Every coding/review agent follows `CLAUDE.md`'s **Search and edit loop**. Native
specialists and `explore` can load either skill and use Code Mode plus both
corpora's query/inspection tools. Read-only specialists leave maintenance to the
caller; exact search remains available without an index.

Check `/mcps` and `/search-status`; connected transport alone is not readiness.
Require `ready_for_queries: true` and inspect integrity, upstream currency and
refresh errors separately. Request `/search-sync` for an immediate source refresh,
or add `--corpus articles` to select literature. Note excerpts retain provenance
and literal locators, not independently verified full-text extraction.

Use OpenCode's targeted native plugin check/update for package releases. A sync
does not install updates, and restarting OpenCode does not necessarily stop the
separate search daemon. The operations guide covers disconnecting every affected
location, stopping the old daemon, reloading and verifying real queries after a
runtime-changing update. Preserve databases and caches; never reinstall an
old repository-local runtime to repair a mismatch.
