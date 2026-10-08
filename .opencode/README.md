# OpenCode project setup

The repository supports OpenCode V2 with:

- `opencode.json`: the routine shell-command allowlist adapted from
  `.claude/settings.json`, plus local code and article-note MCP mounts. Other shell
  commands ask for approval.
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

For repository and processed-article retrieval, use `knowledge-search` and read
[`docs/agent-search.md`](../docs/agent-search.md). `search/profiles.json` describes
the corpus/chunking policies consumed by the MCP/indexer; it is not native
OpenCode configuration.
Every coding/review agent follows `CLAUDE.md`'s **Search and edit loop**. Specialist
prompts point to that shared contract; lookup permissions include Code Mode and
the two LanceDB query/inspection tools. Read-only specialists leave index
maintenance to the caller. The built-in `explore` adapter also permits these
lookups; exact searches and source reads remain available when indexes are absent.

OpenCode mounts `tools.agent_search.mcp` from the isolated `search/` Python
environment. `uv run --frozen` installs its locked dependencies on a fresh clone.
The lightweight stdio proxy registers with one host-shared local search service,
used by both corpus mounts and compatible worktrees/CLI commands. That service
caches one loaded CPU embedding model; proxies import no database/model runtime.
The startup hook fetches and indexes committed upstream `main` in the background,
sharing one logical corpus across clones and branches of the same upstream.
Source bytes and indexing policy come from pinned Git objects. It refreshes
populated corpora and checks upstream every five minutes while mounted.
`table_details` reports the indexed commit, integrity, upstream currency, and
refresh state. Queries await first population; later refreshes retain the last
complete main snapshot. The reader supports vector, FTS, and hybrid lookup.

Use `uv run --project .opencode/search --no-sync python -m tools.agent_search status`
to inspect main-index status, and substitute `sync` for an immediate upstream-main
refresh. `bootstrap` is an alias for this ensure-main operation. Resolve retrieved
anchors locally and inspect `git diff <indexed-commit>` plus relevant untracked
files for branch work. See the guide for identity, migration, and verification.

`lancedb-articles` mounts the same reader with `--corpus articles`. It maintains
a separate fetched-main snapshot of tracked literature notes and supports
vector/FTS/hybrid queries with optional `paper_id` prefilters (note filename stems).
Results identify the note, literal Source/Pagination declarations and page clauses;
`tracked_note`/`note_fragment` describe excerpts, not fresh full-text verification.
No PDF extraction or page-offset inference occurs. Use `--corpus articles` on the
CLI's `status`/`sync` commands and read the note before turning a hit into a citation.

`python -m tools.agent_search.service status` in the isolated environment inspects
the shared PID, model count and mounts; `stop` releases the service. Disconnect
older mounts before stopping it after implementation/dependency changes, then
reconnect updated mounts. Incompatible search implementations are rejected.
The guide covers socket/log locations, lease expiry, idle shutdown and recovery.
