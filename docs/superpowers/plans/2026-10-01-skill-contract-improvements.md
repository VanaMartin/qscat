# Skill-contract improvement plan

## Goal

Make repository skills consistent at stage boundaries and require evidence that
matches the claim being made. Keep shared procedures usable from both Claude and
OpenCode.

## 1. Repair contradictory guidance

- [x] Qualify conservation checks for closed Hermitian versus ECS/open systems;
  replace universal monotonic convergence with method-appropriate evidence and
  explicit tolerances.
- [x] Verify history reconstruction on the original base before integration;
  separately verify rebased changes, with explicit unstaging and isolated backups.
- [x] Make structured-report persistence caller-owned in both agent runtimes and
  update the consuming skill workflows.
- [x] Require dynamic-reference checks and provenance preservation before deletion.
- [x] Align container guidance with current stages, extras, version reporting,
  resource limits, and test tiers.
- [x] Require actual-observable evidence for tuner completion and distinguish
  production convergence, proxy evidence, and deferral.
- [x] Make lifecycle porting and independent physics-review gates explicit.
- [x] Check local references, runtime loading, and affected documentation contracts;
  review the final diff for contradictions and update the changelog.

### Step 1 verification record

- Documentation portability: 196 tests passed.
- Disposable Git exercise: preserved additions, deletions, executable file modes,
  and reconstruction tree identity; integrating upstream changed the whole tree
  while preserving the reconstructed branch delta.
- OpenCode loaded all six specialists and the three caller-owned report contracts.
- Skill frontmatter, local links, command references, and search-profile JSON checked.
- Independent workflow and physics reviews completed; findings corrected in the
  adapters, Docker guidance, and tuner boundary/refinement requirements.

## 2. Mechanically validate evidence and reports

Design a shared verification runner and versioned report schemas. Record scope,
source fingerprint, command, environment capabilities, exit status, result counts,
and exclusions. Guard test parallelism and validate report coverage before synthesis.

## 3. Add artifact/configuration contracts

Define separate integrity, authority/containment, offline-clone, semantic-round-trip,
and provenance checks. Add focused verification and a contract reviewer. Coordinate
test lifecycle work with the existing scaffolding/WIP proposal.

## 4. Review claims and publication handoffs

Add claim propagation and PR-body/final-diff reconciliation. Introduce focused
test-evidence and claim reviewers with caller-owned verification and persistence.

## 5. Evaluate workflow improvements

Exercise historical failure scenarios: inherited branch commits, failed or empty
test runs, mock/live divergence, lost configuration settings, probe-only convergence,
stale withdrawn claims, and duplicate regression coverage. Compare existing and
revised workflows using observable actions and outputs.

## Scope and checkpoints

Each numbered step is a separate implementation slice. Step 1 changes guidance
and agent contracts; later slices add executable enforcement. Permanent guidance
must explain its contracts without depending on this plan. Review changed guidance
against implementation and durable architecture decisions before marking a step done.

## Additional retrieval work

- [x] Measure owned code and durable documentation; inspect the installed LanceDB
  connector and embedding model rather than assuming advanced search is supported.
- [x] Add an exact/semantic routing skill, corpus policies, and a durable design
  note for repository and processed-article retrieval.
- [x] Configure a separate machine-local OpenCode article server/table selection.
- [x] Research existing code-indexing implementations and define stable anchors,
  content-keyed embedding reuse, and moving location metadata.
- [x] Verify the shared search/edit contract and lookup access across primary,
  built-in exploration, and both sets of specialist agent adapters.
- [x] Evaluate CocoIndex's current dependency compatibility; exercise the
  server-compatible checkpoint writer for line-only moves, changes/deletions,
  tokenizer limits, and worktree isolation. Managed CocoIndex watching is deferred
  to a coordinated reader/writer upgrade.
- [x] Populate the code corpus and verify live MCP retrieval and current freshness.
- [x] Implement managed code/note metadata, vector/BM25/hybrid query modes,
  committed-source reconciliation and health reporting; article lookup adds paper
  prefilters. Processed full-text ingestion remains a separate provenance contract.
- [ ] Evaluate retrieval quality and latency before selecting an optimized model
  or approximate vector index.

### Committed-main and article-note follow-on

- Code membership and indexing policy now come from fetched upstream-main Git
  objects. Clones share an upstream-owned corpus; mount and bounded periodic
  refreshes publish immutable generations through an atomic manifest pointer.
  The earlier worktree/empty-only population records below describe superseded
  stages, not current ownership or refresh behavior.
- The same writer/reader supports a separate article-note corpus. It indexes
  tracked literature notes, preserves literal edition and page declarations,
  marks oversized fragments, and returns note provenance without claiming fresh
  full-text verification. Paper-filtered hybrid lookup explicitly prefilters both
  vector/BM25 candidates to avoid the installed SDK's reversed filter forwarding.
- Article-note integration and publication regression checks: 24 passed.
  Documentation portability: 196 passed. Main advanced during the work; the
  maintainer incorporated three new notes, reaching 21 files / 1,633 chunks at
  `0f9768e`. Full-corpus note/blob hashes, source spans, payload/model keys and
  tokenizer budgets passed. An unchanged refresh changed zero rows/embeddings.
- Six fixed known-paper/page smoke cases hit the expected source and literal page
  within eight results in vector, FTS and hybrid modes. Warm full-query-function
  median times were 11.1 / 5.9 / 13.6 ms respectively. This is not exhaustive
  recall or evidence for an optimized model/ANN index. Local evidence is under
  `.superpowers/audit/article-notes/`.

### Shared-runtime follow-on

- Code/article MCPs are lightweight proxies into one local service per database;
  compatible worktrees and CLI commands share a strong loaded-transformer cache.
  Corpus publication and refresh ownership remain independent. Runtime diagnostics
  expose process/model/mount counts, and leased mounts support crash recovery,
  source-checkout failover and bounded idle shutdown.
- Unchanged refreshes check metadata and selected blobs before loading rows or
  vectors. Tests prohibit both row materialization and model loading on this path,
  including excluded-only main advances.
- Combined search/runtime and documentation checks: 228 passed (32 integration,
  196 portability) in 81.58 seconds. The cross-worktree test additionally ran MCP
  processes from separate checkout directories containing identical tooling bytes.
  Ruff lint/format, whitespace checks and the strict Sphinx HTML build passed.
- Live code/article connections were reconnected to the shared implementation and
  both report complete/current `0f9768e` snapshots. Refresh-only operation loaded
  zero models; vector/hybrid lookups brought the shared count to one.
- A controlled aggregate-RSS snapshot compared reconstructed per-corpus model
  processes (including the old no-op row-materialization step) with the shared
  service plus proxies, on the same real code/note snapshots. Two mounts measured
  1,154 → 411 MiB; four mounts across two worktrees measured 2,318 → 560 MiB.
  These macOS RSS sums include shared pages and are not unique physical memory or
  a guaranteed reduction for other workloads. All six code/article query-mode
  comparisons returned identical hits/scores; simultaneous forced refreshes and
  queries passed, with zero new embeddings/parsed files/row changes. Evidence is
  under `.superpowers/audit/search-runtime/`.

### Initial code population

- CocoIndex dependency inspection found that its current LanceDB extra requires
  LanceDB >=0.34 and Python >=3.11; the connected server uses LanceDB 0.21.2 and
  Python 3.10. Initial ingestion uses a compatible checkpoint writer in the
  isolated Python 3.12 `.opencode/search` environment.
- Implemented Python AST, Rust Tree-sitter, and Markdown-heading chunks; separate
  embedding payloads and moving provenance; exact payload/model vector reuse;
  per-file catch-up, ID reconciliation, worktree ownership, and freshness states.
  Native BM25 and hybrid retrieval are available through the companion CLI.
- Real parser/embedding/LanceDB integration checks: 8 passed, covering line moves,
  content edits, renames/deletions, parse failures, explicit untracked additions,
  ownership, tokenizer limits, and external mutation recovery. Ruff passed.
- Documentation portability: 196 passed. Strict Sphinx build passed.
- Code corpus populated: 518 selected files, 507 searchable files, 11 empty files,
  and 11,549 chunks. CPU MiniLM baseline, 384 dimensions; maximum full payload
  236 tokens against the 256-token limit. Exact vector scans and native BM25 FTS.
- Full-corpus verification found and corrected a trailing-blank-line locator bug;
  all stored excerpts now fit their source spans, and token/hash/ID checks pass.
  The correction reused 11,544 vectors; five new embeddings came from edited
  tooling/tests. An unchanged catch-up then reparsed zero files, generated zero
  embeddings, and changed zero rows (5.44 s, including model/environment startup).
- Connected MCP table inspection reports 11,549 rows; live vector queries return
  provenance-bearing results. Separate Python 3.10/PyArrow 19 reader compatibility
  was also exercised against the Python 3.12 writer's disposable table.
- Five known-source smoke cases: vector 4/5, FTS 5/5, hybrid 4/5 hit-at-5;
  respective warm-process median query times 10.5, 5.7, and 11.8 ms. Vector/hybrid
  missed the predeclared grid-evidence sources in the top five, though they found
  neighboring convergence code. This is a small smoke set, not exhaustive recall
  or evidence to select a different model/ANN index. Raw vectors occupy 17.7 MB.
- Local evidence is saved under `.superpowers/audit/skill-contracts/` in
  `index-population.json`, `index-location-refresh.json`, `index-idempotence.json`,
  `index-evaluation.json`, and `index-cli-query.json`. Current code manifest is
  content-addressed and worktree-owned; the scientific-article corpus is pending.

### MCP mount bootstrap

- Added `tools.agent_search.mcp` with a lifespan hook that runs initial population
  asynchronously under the indexer's writer lock. Absent/empty indexes are
  populated; later mounts of a populated table reuse it. Queries await the mount
  job and inspection exposes progress/failures. CLI `bootstrap` fires the same
  empty-only operation explicitly.
- Project `opencode.json` mounts this reader through the isolated environment's
  frozen lockfile, with a startup allowance for first-install dependency setup.
  The existing owned base table is reused; new clones select a deterministic
  root-scoped table. CLI and MCP share that selection policy.
- Integration checks: 13 passed, including real stdio MCP handshakes and queries,
  two concurrent mounts, newly cloned checkout isolation, empty-schema adoption,
  initial parse failure reporting, and live policy changes followed by checkpoint
  catch-up. Ruff, formatting, whitespace, 196 documentation portability tests,
  and the strict Sphinx build passed.
- Live OpenCode reports both MCP connections connected. Its code reader reports
  bootstrap complete with zero new embeddings on a populated remount, current
  freshness for 519 files / 11,595 chunks, and successful hybrid lookup. The
  updated corpus was caught up; `index-mount-hook-refresh.json` records it.

### Retrieval contract verification

- Both runtimes' six specialist prompts point to the shared source-resolution loop.
- OpenCode's effective loaded permissions allow query/inspection and Code Mode for
  all six specialists and built-in exploration; read-only edit/ingestion denials
  remain effective. Build, plan, and general also have lookup access.
- A disposable LanceDB 0.21.2 exercise confirmed that an unrelated edit can move a
  method's span and file fingerprint while retaining its logical/payload keys;
  updating only location metadata preserved stored text/vector values. A scoped
  deletion retained the unrelated file. This verifies SDK operations, not a
  working CocoIndex integration or retrieval quality.
- Strict Sphinx build passed with this worktree's `qscat` package selected;
  skill/agent frontmatter, local links, profile JSON, and whitespace checks passed.
