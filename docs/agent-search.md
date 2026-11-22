# Repository and literature search for coding agents

Search is provided by the native
[opencode-repo-search plugin](https://github.com/VanaMartin/opencode-repo-search),
not a Python search engine in this repository. The package owns the MCP launchers,
isolated locked Python environment, parsers, indexes and shared search daemon.
qModeling owns its corpus policy and agent workflow.

Use `repo-knowledge-search` when the plugin supplies it. The tracked
`.claude/skills/knowledge-search/SKILL.md` remains a guidance adapter for existing
specialists and other clients; it does not install or start another runtime.

## Install and connect

Requirements are OpenCode V2 >=2.0.15, <3, Git and
[uv](https://docs.astral.sh/uv/) on macOS or Linux. Search dependencies stay outside
the numerical workspace; do not add them to the application's environment.

For a new host, install the plugin once:

```sh
opencode plugin add 'github:VanaMartin/opencode-repo-search#latest'
```

If an entry already exists, replace its target in its owning OpenCode configuration
instead of adding a duplicate. Preserve object-form options. `latest` is a movable
Git tag for accepted releases, which can include release candidates; it does not
mean stable or install updates automatically. Use an immutable version tag or full
commit instead when reproducibility is required.

Reload the project location after installing or changing the package. This
repository supplies native defaults under `.opencode/search-package/`, so a clone
with those files needs no repeated initialization. The plugin registers:

| Corpus | MCP | Material |
|---|---|---|
| `repository` | `lancedb` | Owned code, tests, projects, validation, benchmarks and durable documentation |
| `articles` | `lancedb-articles` | Tracked literature notes, excluding their inventory README |

`opencode.json` contains permissions, not duplicate search process definitions.
Without the native plugin, search is unavailable rather than falling back to an
outdated repository-local engine.

In OpenCode, inspect the connections and their actual readiness:

```text
/mcps
/search-status
```

Connected transport is not a usable index. Require `ready_for_queries: true`, and
inspect integrity, upstream currency and refresh errors independently. Initial
indexing can take several minutes. Use `/search-sync` only when an explicit
upstream refresh is wanted; it does not install a plugin update.

## Settings, policy and linked worktrees

`.opencode/search-package/runtime.json` selects `origin/main` and the native profile
at `.opencode/search-package/profiles.json`. The profile uses
`repo-search-profiles-v1`, with independent enabled code and literature corpora.
Their include/exclude rules and embedding settings retain qModeling's search
boundaries. Reference trees, working plans, generated results and defective test
fixtures stay outside owned-code retrieval. PDFs are not ingested.

The repository defaults explicitly select **local policy**: the package validates
the effective profile and pins an immutable snapshot by SHA-256. This permits
reviewing a policy on a branch without pretending it has been published on main.
It never permits indexing working-tree source. A profile edit changes its hash
and can give every enabled corpus a new index identity and require rebuilding.
Use the corpus commands to review changes before applying them.

Effective configuration selection is explicit CLI overrides, worktree-private
settings, repository-private settings, then tracked settings. Status names the
selected physical path, scope, inherited state and overrides. Repository-private
settings live in `repo-search/package.json` inside the resolved Git common
directory; never hardcode an absolute `.git` path for linked worktrees.

For an existing setup, preview and explicitly adopt repository scope when desired:

```text
/search-init git-ignore --scope repository --adopt --dry-run
/search-init git-ignore --scope repository --adopt
```

Adoption retains policy, upstream and compatible index identity. It does not stage,
commit or push. Existing repository-private settings must agree with the adopted
files. Unchanged adopted files defer to repository settings; changing a preserved
private source file creates a worktree override. Do not reinitialize every worktree.
An already adopted installation needs no further setup after moving a session.

Shared policy is an alternative explicit choice: `policy_source: upstream` reads
the profile from the same fetched source commit. Publish that profile on the
configured upstream branch before using it. Do not change policy source or
database selection merely to repair a runtime mismatch.

## Source and publication contract

Every searchable source file comes from Git blobs at one exact commit fetched from
the configured upstream branch. The checked-out branch, local main, staged edits,
dirty files and untracked files never enter that committed-source snapshot.
Shared policy comes from the same commit; local policy provenance is separately
identified by hash. Fetching does not reset, rebase or change the task checkout.

Each corpus publishes a complete immutable generation through an atomic manifest
pointer. Failed refreshes retain the last complete publication. Different corpora
retain separate ownership, tables and publication state. Clones/worktrees with
matching upstream, branch, profile and local-policy identity share the corresponding
database and daemon by default; unrelated upstreams remain independent.

`table_details` and query freshness identify the source commit, policy provenance,
integrity, `ready_for_queries`, `upstream.state`, refresh phase and errors. A
complete snapshot can be behind upstream or have unknown remote currency. With
no complete snapshot, use local source and report unavailable retrieval. An empty
result establishes neither absence nor completeness.

## Queries and evidence

Discover the actual tool catalog before calling it. Each corpus exposes
`query_table` and `table_details`. Queries support `vector`, `fts` and `hybrid`;
start with `top_k=5` for code or `8` for notes. Literature queries can prefilter by
`paper_id`, the note filename stem. Inspection selections must match the mounted
corpus/database rather than silently inspecting something else.

Read current source after at most two exploratory queries. Resolve paths and
qualified symbols/headings against the task checkout; indexed line ranges are
snapshot hints, not edit offsets. Inspect the relevant branch delta:

```sh
git diff --name-status <indexed-commit> -- <scope>
git diff --cached --name-status -- <scope>
git diff --name-status -- <scope>
git ls-files --others --exclude-standard -- <scope>
```

Read changed bodies, omit deleted hits and resolve ambiguous anchors explicitly.
At handoff, name the indexed commit, current integrity/currency and relevant local
divergence. Read-only specialists leave explicit maintenance to the caller. Never
use an append-only ingestion connector to refresh these managed corpora.

Literature results preserve note hashes, Source/Pagination declarations and literal
page/locator clauses. `tracked_note` means a tracked excerpt, not independently
verified full text. Oversized splits are `note_fragment`; read adjacent note lines
to recover a complete equation or qualification. `extraction_page` is null: the
index processes no PDF pages and infers no page mapping. A note can discuss another
paper's equations, so read the owning note before attributing a claim. Follow the
tracked [literature conventions](../reference/literature/README.md) for citations.

## Package updates and daemon handover

Use OpenCode's native updater with the exact configured target:

```sh
opencode plugin list
opencode plugin check 'github:VanaMartin/opencode-repo-search#latest'
opencode plugin update 'github:VanaMartin/opencode-repo-search#latest'
```

Do not omit the target and update unrelated plugins. Advisory checks/notices do not
authorize installation, revision changes, service restarts or index deletion.
Versioned tags and full-commit pins do not advance to a different release.

There are separate OpenCode, MCP-proxy and search-daemon lifetimes. Restarting a
terminal session or OpenCode alone does not guarantee that the search daemon stops.
For an update that changes runtime compatibility:

1. Record the effective package launcher, database, fingerprint, daemon PID and
   proxy PIDs from `/search-status` and `table_details.service`.
2. Disconnect affected MCPs in every loaded location sharing the database. Identify
   surviving orphaned proxies by command, root and database; stop only those verified
   processes. Otherwise their heartbeats can restart the old service.
3. Append `service stop --root /path/to/worktree` to the old package launcher.
   Verify actual process exit and that the old daemon is not recreated. Service
   health/stop remain available across implementation mismatches.
4. Update the native package, reload affected locations and reconnect their MCPs.
5. Verify matching runtime fingerprints, correct database and policy, both corpora's
   readiness/integrity/currency, and a real query in each corpus. Compatible mounts
   should share one daemon and embedding model.

The compatibility fingerprint includes source, lockfile, dependencies and Python
patch version. A Python update alone can therefore require a handover. A generic
code/dependency mismatch is not evidence of Git or index corruption.

Use the exact launcher supplied by the package context or `/search-*` commands for
direct CLI operations, followed by the command, `--root` and optional `--corpus`.
Do not guess cache generation paths, copy the runtime into the project, or retain
an independent dependency lock. Runtime sockets/logs live under
`~/.cache/repo-search/<database-path-hash>/`; database and prepared-runtime paths are
reported separately. Preserve published indexes, old databases and caches during
updates. An old index-format migration may require a new index, not destructive
in-place reuse or deleting its predecessor.

## Verification and ownership

The package's own suite tests parsers, embeddings, Git snapshots, atomic
publication, query modes, daemon lifecycle and actual MCP transport. qModeling's
checks enforce its native configuration, corpus boundaries, agent permissions and
the absence of a duplicate executable search engine:

```sh
uv run --no-sync pytest tests/test_search_integration.py tests/test_docs_portability.py -q
```

After configuration changes, verify live `/mcps` and `/search-status` in the intended
worktree and test code plus literature retrieval. Configuration tests cannot prove
that a daemon is compatible or an index is query-ready.

Further documentation:

- [Package operations and worktrees](https://vanamartin.github.io/opencode-repo-search/)
- [OpenCode V2 plugins and native updates](https://opencode.ai/v2/docs/plugins)
- [OpenCode V2 MCP connections](https://opencode.ai/v2/docs/mcp-servers)
