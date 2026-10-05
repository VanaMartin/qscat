# Repository and article search for coding agents

Use local exact search for known identifiers and LanceDB for conceptual discovery.
Read current source before acting on an indexed result. The `knowledge-search`
skill defines this routing; `.opencode/search/profiles.json` records the two corpus
policies. The MCP reader and checkpoint indexer in `tools/agent_search/` consume
the repository profile. OpenCode's project configuration mounts the code reader
automatically; its startup lifecycle bootstraps an absent/empty code table.

## Corpus boundaries

The repository profile includes owned source, tests, example configs, durable
physics/API/architecture documentation, and operating guidance. Exclude working
plans/specs, generated outputs, build directories, fixtures containing deliberate
defects, and legacy reference trees from ordinary retrieval. These can be searched
explicitly with local tools when the task requires them.

Keep processed scientific articles in `qscat_articles`, separate from repository
search in `qscat_knowledge` and its worktree-scoped variants. Start with tracked
`reference/literature/*.md` notes.
Add processed full text only from an explicit source root with provenance. Notes
and extractions are different source kinds; an extraction is not automatically a
verified reference claim. Local article data does not replace the tracked notes.

## Basic external connector

The basic external LanceDB MCP server, retained for the article connection, uses
environment variables `LANCEDB_URI`, `TABLE_NAME`,
`EMBEDDING_FUNCTION`, and `MODEL_NAME`. It creates a two-column `doc`/`vector` table
on first ingestion and appends strings. It does not read filenames or URLs passed
as strings. Supply actual chunk content, not a path to be embedded as if it were
the document.

Its current `query_table` always performs vector search even when `query_type`
requests another mode. `table_details` opens the configured table regardless of
its arguments. A connected server therefore proves neither a populated table nor
working hybrid retrieval. Table inspection and queries fail until a table exists.

Use the CPU `all-MiniLM-L6-v2` baseline already configured. Its installed model has
384 dimensions and a 256-token input limit. Budget with the model's tokenizer:
target 160-token content units, cap at 220, and ensure the complete embedded text
including the provenance header and special tokens fits 256. Overlap up to 24
tokens only when splitting a large unit; do not duplicate whole adjacent symbols.
These are starting settings to evaluate, not demonstrated optimal chunk sizes.

With the basic connector, put a compact chunk ID, path, and qualified anchor in
each `doc` string and keep snapshot/hash/location provenance in a local manifest
keyed by that ID. This provides leads but cannot substitute for server-side
filtering or incremental reconciliation. Avoid appending successive repository
snapshots to the same unfiltered table: stale matches would compete with current
code. Rebuild an explicitly owned derived table or use immutable snapshot tables
with explicit server selection; never erase another corpus to refresh one checkout.

## Stable documents and moving locations

Code has three distinct keys:

| Key | Depends on | Used for |
|---|---|---|
| Logical document anchor | repository, relative path, language/kind, qualified symbol or heading, duplicate disambiguator | Resolving the unit in current source |
| Chunk/content key | exact chunk text plus any deliberately embedded parent context | Deduplicating content and reconciling changed fragments |
| Embedding cache key | complete embedding payload, model revision, document/query parameters | Reusing expensive vectors |

Line ranges, whole-file hashes, snapshot/worktree IDs, and timestamps belong to
location/freshness metadata. Keep them outside the embedded text and its cache
key. An edit earlier in a file then changes positions and the file fingerprint,
while an unchanged function's payload and embedding remain reusable. Context that
is actually embedded, such as a signature or module name, remains part of the
payload: changing it can legitimately require a new vector.

Store module headers, functions/methods, compact classes/result holders, and
fragments of large units. Keep each fragment attached to its parent anchor. Reuse
content keys rather than making a global chunk ordinal the sole row identity;
inserting an earlier function must not renumber unrelated documents. Identical
fragments still need distinct source locations. Repeated declarations and headings
need an explicit disambiguator; an ambiguous match must be resolved against source.
Renames/moves create new path/name anchors and retire the old ones, with embedding
reuse only when the complete payload still matches.

### Worked editing example

Suppose `SparseLU.refactor` moves from line 300 to line 320 because a helper is
added above it. The logical anchor is still the file plus `SparseLU.refactor`.
The indexer reparses that file, updates its source fingerprint and the method's
span, and reuses the vector if the payload is identical. These line numbers are
illustrative, not source locators for this checkout.

The agent resolves `SparseLU.refactor` locally and reads its current body before
editing. If that body subsequently changes, its changed fragments get new content
keys and embeddings. A symbol-level content match does not prove its dependencies
or behavior are unchanged: imports, callees, configuration, and tests still need
the usual change-impact review.

### Index maintenance during editing

Use one incremental writer with per-file change detection, embedding caching,
and idempotent reconciliation. A filesystem watcher can update saved files; a
checkpoint catch-up covers changes made by shell commands, Git operations, and
other agents. Verify selected worktree and catch-up status before relying on
freshness. Parser errors in partially edited files must leave their results
explicitly stale/unavailable until a successful reconciliation.

Agents follow `CLAUDE.md`'s **Search and edit loop**: discover, resolve/read live
source, edit/validate, then check catch-up at a coherent checkpoint. Read-only
specialists keep their assigned evidence/report scope and leave maintenance to
the caller. If the writer is unavailable, local source remains usable and the
caller records the affected unsynchronised paths. Agents do not hand-maintain
vectors or derived document copies.

Keep worktrees isolated through separate selected manifests/tables or overlays.
An overlay must mask base locations for changed/deleted paths before adding the
working versions. Explicit task-owned new files can enter that overlay; the
tracked-files-only base policy still excludes arbitrary untracked outputs.

## MCP mount-time population

`opencode.json` configures the project-local `lancedb` server to run
`python -m tools.agent_search.mcp` through `uv run --project .opencode/search
--python 3.12 --frozen`. Opening a fresh clone in OpenCode installs the locked
search environment and starts this reader. The startup timeout allows dependency
setup on a new machine; embedding-model loading and indexing run after mount in a
background task so the MCP handshake and tool catalog are available immediately.

The server's lifespan invokes the indexer's `bootstrap` operation under the
exclusive writer lock. It populates an absent or empty selected table, and can
adopt an empty table with the basic connector's schema. A populated table is
skipped, including when its source is stale: the editing caller still owns
checkpoint `sync`. Concurrent mounts recheck emptiness after acquiring the same
lock, so only one performs initial ingestion. Ownership checks remain effective.

`table_details` exposes the selected table, root, bootstrap `pending`/`running`/
`complete`/`failed` state, row count, schema, and source freshness. `query_table`
awaits the mount job and supports `vector`, `fts`, and `hybrid`; results include
freshness, source metadata, excerpts, and scores. Bootstrap errors are reported by
inspection and queries. A stale result remains a lead requiring live-source
resolution; absent or incomplete indexes are not queried. The reader exposes
lookup/inspection tools; the single indexer owns derived writes.

The same operation can be fired manually:

```bash
uv run --project .opencode/search --no-sync python -m tools.agent_search bootstrap
```

Other MCP clients can mount the same module from the repository root with this
command and startup allowance. The implementation uses the MCP server's lifespan
rather than a client-specific plugin hook.

## Populate, refresh, and inspect the code index

The repository includes a checkpoint writer compatible with the connected
LanceDB 0.21.2 server. Its isolated Python 3.12 environment keeps indexing
dependencies outside the numerical workspace. Run these commands from the
repository root:

```bash
uv sync --project .opencode/search --python 3.12
uv run --project .opencode/search --no-sync python -m tools.agent_search plan
uv run --project .opencode/search --no-sync python -m tools.agent_search sync
uv run --project .opencode/search --no-sync python -m tools.agent_search status
uv run --project .opencode/search --no-sync python -m tools.agent_search query \
  "sparse factorization symbolic analysis reuse" --mode hybrid --limit 5
```

`sync` is the supported caller-owned catch-up operation at edit/test checkpoints.
It fingerprints selected files, reparses only changed files under the same profile
and chunker version, reuses exact payload/model embeddings, and reconciles removed
rows by ID. Python AST, Rust Tree-sitter items, and Markdown headings supply
anchors; other text/config formats use bounded file fragments. It checks the
complete payload with the installed MiniLM tokenizer before embedding. Empty
files are represented in the freshness manifest without searchable chunks.

The default database is `~/.local/share/opencode/lancedb`, or `LANCEDB_URI` when
set. MCP and CLI use the same automatic selection: reuse `qscat_knowledge` when
its manifest already belongs to this resolved worktree root; otherwise select
`qscat_knowledge_<16-character root SHA256>`. A newly cloned copy on the same
machine therefore selects an independent table. `--db`, `--table`, and `--root`
select explicit locations; the MCP module also accepts `TABLE_NAME`. Every table
has a sibling `<table>.manifest.json` and an exclusive writer lock. Ownership
checks reject reuse by another worktree. Changing the embedding model/configuration
requires a new table.

Tracked files form the base corpus. Include an explicit task-owned new source with
`sync --include-new relative/path.py`; this selection persists in the manifest.
Excluded sources remain excluded. Ordinary untracked outputs never enter through
a recursive directory scan. Parsed chunks, vectors, and manifests are derived
local data; source files are the editing surface.

`status` reports `current`, `stale`, `incomplete`, or `absent`, and lists changed
paths. Parse failures preserve the previous committed rows; the changed source
fingerprint reports their staleness. Interrupted publication and external table
mutations are detected through manifest state, table version, and row count.
Retry `sync` after a coherent source checkpoint. CLI queries verify ownership and
freshness first. The project MCP reader verifies ownership and includes freshness
in each result; its mount job handles initial population, while caller checkpoints
handle subsequent source changes.

The writer stores a tokenizer-bounded `text` embedding source separately from the
`doc` excerpt returned by the existing MCP reader. `doc` carries path, anchor,
declaration disambiguator, lines, worktree root, file hash, and chunk ID. Those
moving locations are not embedded. Structured columns carry the same provenance
for metadata-aware consumers. The companion CLI supports vector, native BM25 FTS,
and reciprocal-rank-fused hybrid search, also supported by the code MCP reader.
The external article connector remains vector-only. Exact scans remain the vector
baseline.

Run the isolated integration checks with:

```bash
uv run --project .opencode/search --no-sync python -m pytest tests/test_agent_search.py -q
```

They exercise real embeddings and LanceDB writes, unchanged-symbol line moves,
edits, deletions/renames, parse failure, explicit new-file selection, worktree
ownership, external mutations, and all three CLI query modes. MCP protocol checks
exercise concurrent mounts, automatic bootstrap, fresh-clone isolation, empty
schema adoption, and visible bootstrap failures. The numerical
workspace can run the parser-only checks without installing search dependencies.

## Target index design

The metadata-aware design needs the following capabilities in the MCP server or
a companion indexer:

1. **Symbol/section chunking.** Use Python AST boundaries, Rust item boundaries,
   Markdown headings, and configuration units. Attach parent class/module or
   heading context. Split oversized units with source line ranges.
2. **Structured identity.** Repository rows carry repository/worktree identity,
   snapshot ID, path, language, qualified symbol/heading anchor, derived line
   range, file hash, chunk/payload hashes, embedding model/revision/parameters,
   chunker version, and source kind. Article rows additionally
   carry paper identity/edition, printed and extraction pages, locators, note
   path, extraction version, and verification status. Missing article locators
   remain explicitly unverified; do not manufacture them.
3. **Incremental updates.** Reuse embeddings for unchanged content, replace all
   locations/chunks of changed files, remove deleted files, and publish the new
   snapshot manifest after the update succeeds. Update spans without reembedding
   identical payloads. Dirty worktrees need a distinct overlay or direct local
   search. Identify snapshots by their content manifest, not branch name alone.
4. **Hybrid retrieval.** Add BM25 to semantic search, fuse with reciprocal-rank
   fusion, and prefilter repository/snapshot or paper/version before ranking.
   Preserve identifiers, one-letter physics variables, Greek symbols, equation
   numbers, and DOI tokens in lexical text. Start with stemming and stop-word
   removal disabled for these corpora; retain exact local search for punctuation-
   sensitive identifiers and formulas.
5. **Bounded results.** Return source metadata, a short excerpt, and scores;
   suppress embeddings and full files. Retrieve a small candidate set, deduplicate,
   then read selected source regions. Group article hits by paper and repository
   hits by file so one long document cannot fill the answer.

### Vector index choice

Start with exact scans for both corpora. A few thousand 384-dimensional vectors
have a modest raw footprint; an ANN index is an optimization with a recall cost,
not a prerequisite for vector search. Add a vector index only after measuring
warm/cold latency and filtered recall against an exact-search baseline.

The inspected connector environment uses LanceDB 0.21.2. Its `create_index` supports
`IVF_FLAT`, `IVF_PQ`, `IVF_HNSW_SQ`, and `IVF_HNSW_PQ`; newer documentation also
describes indexes absent from that signature. Select by the installed SDK rather
than copying a new index name blindly. For a recall-first trial, evaluate
`IVF_FLAT` with cosine distance; tune partition/probe counts using the actual corpus.
Consider compression only when measured memory or latency warrants it. Keep index
and query metrics identical. Model changes require reembedding into a new table,
even when dimensions happen to match.

The existing connector cannot enable this target design through environment
variables alone. It needs query-mode handling, provenance fields, metadata
filters, replacement semantics, FTS/index maintenance, and health/score reporting.

## Existing implementations to build on

The following primary sources support the design; they are implementation
references, not measured retrieval results for qModeling.

- **CocoIndex's LanceDB code example** uses Tree-sitter-aware splitting, a row
  containing code/path/span/vector, content-dependent ID generation, and managed
  incremental writes. It embeds `chunk.text` and stores start/end lines separately.
  Its [source](https://github.com/cocoindex-io/cocoindex/blob/main/examples/code_embedding_lancedb/main.py)
  and [walkthrough](https://github.com/cocoindex-io/cocoindex/tree/main/examples/code_embedding_lancedb)
  are the preferred starting point for evaluating an incremental writer.
- **CocoIndex's memoization and LanceDB target docs** explain separating stable
  identity from freshness, reusing computations, and reconciling upserts/deletions.
  The target connector also owns compaction/index maintenance; keep external
  connections read-only while it owns writes.
  See [memoization keys](https://cocoindex.io/docs/advanced_topics/memoization_keys/),
  [ID generation](https://cocoindex.io/docs/common_resources/id_generation/), and
  [LanceDB connector](https://cocoindex.io/docs/connectors/lancedb/).
- **Continue's LanceDB implementation** separates cached content/vectors from
  workspace/branch index selection and handles compute/add/remove/delete changes.
  It corroborates content caching and location isolation, although its update
  granularity and IDs differ from the proposed symbol anchors.
  See [LanceDbIndex](https://github.com/continuedev/continue/blob/main/core/indexing/LanceDbIndex.ts)
  and [CodebaseIndexer](https://github.com/continuedev/continue/blob/main/core/indexing/CodebaseIndexer.ts).
- **Aider's repository map** combines signatures with dependency-based ranking
  to select useful context. Use our measured code map for structural/impact facts;
  vector similarity cannot establish caller coverage.
  See [repository map](https://aider.chat/docs/repomap.html).
- **LanceDB's update and reindexing docs** distinguish updating row data from
  maintaining FTS/ANN indexes. Reconciliation must retire removed source rows;
  compaction/index optimization alone does not find deleted source files.
  See [table updates](https://docs.lancedb.com/tables/update) and
  [reindexing](https://docs.lancedb.com/indexing/reindexing).

CocoIndex's current LanceDB extra requires Python >=3.11 and LanceDB >=0.34,
whereas the existing database uses LanceDB 0.21.2. The repository reader/writer
uses Python 3.12 with that database version; the external article connector uses
Python 3.10. Revisit CocoIndex for managed watching/reconciliation when upgrading the
reader/writer stack together. Adapt its source patterns, provenance, tokenizer
budget, and worktree selection to this repository. The example's generic chunk
sizes need tokenizer verification for MiniLM; its demo score conversion is not a
calibrated similarity measure.

The packaged `cocoindex-code` tool currently uses SQLite for its search target
(see its [indexer](https://github.com/cocoindex-io/cocoindex-code/blob/main/src/cocoindex_code/indexer.py)).
Its MCP refresh-on-search and agent workflow are useful references, but installing
it is not a LanceDB configuration. No CocoIndex writer is installed by these
profiles. Its watcher and refresh-on-search behavior are not supplied by our
checkpoint CLI; MCP metadata filters and article ingestion remain future work.

## Secondary OpenCode connection

Configure a second local MCP server named `lancedb-articles`, using the basic
external server command and database URI with `TABLE_NAME=qscat_articles`. This isolates
scientific retrieval without changing the primary server's table. Native OpenCode
V2 config uses `mcp.servers`:

```jsonc
{
  "$schema": "https://opencode.ai/config.json",
  "mcp": {
    "servers": {
      "lancedb-articles": {
        "type": "local",
        "command": ["uv", "--directory", "{env:LANCEDB_MCP_ROOT}", "run", "lancedb_mcp.py"],
        "environment": {
          "LANCEDB_URI": "{env:LANCEDB_URI}",
          "TABLE_NAME": "qscat_articles",
          "EMBEDDING_FUNCTION": "sentence-transformers",
          "MODEL_NAME": "all-MiniLM-L6-v2"
        }
      }
    }
  }
}
```

Set the machine-local paths in the environment available to the OpenCode service,
or configure absolute paths in global configuration. Preserve the original server
settings; do not put workstation paths in tracked project config. A second server
loads another model process, so measure total RSS alongside numerical workloads.
A future corpus-selecting server could multiplex both tables in one process.

Use `opencode mcp list` to confirm connection. Discover tools through OpenCode's
catalog; the article server has its own namespace. Read-only lookup needs query
and inspection access, not ingestion access. Specialist agents with deny-by-default
permissions also need explicit access to the search skill, MCP lookup tools, and
Code Mode when applicable; do not grant write/ingestion tools merely for retrieval.

## Evaluate before calling it optimal

Use repository questions with known relevant sources: outgoing-flux DA extraction,
ECS c-product, sparse symbolic reuse, resolved-config TD packet round trips, and
grid convergence versus proxy evidence. Use article questions with verified
printed pages and locators in tracked notes. Measure relevant-source recall, stale
or superseded hit rate, locator accuracy, latency, result tokens, and total memory.
Compare local search, vector-only retrieval, and hybrid retrieval on the same set.
Tune chunk size/model/ANN settings from those results, keeping paper and repository
evaluations separate.

## Software documentation

- [OpenCode V2 MCP configuration](https://opencode.ai/v2/docs/mcp-servers)
- [LanceDB hybrid search](https://docs.lancedb.com/search/hybrid-search)
- [LanceDB full-text search](https://docs.lancedb.com/search/full-text-search)
- [LanceDB vector indexes](https://docs.lancedb.com/indexing/vector-index)
- [Tracked literature and locator conventions](../reference/literature/README.md)
