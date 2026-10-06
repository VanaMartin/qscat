# Repository and article search for coding agents

Use local exact search for known identifiers and LanceDB for conceptual discovery.
Code retrieval is a shared baseline of **committed upstream `main`**. Read current
source and inspect its local delta before acting on a result. The `knowledge-search`
skill defines agent routing; `.opencode/search/profiles.json` describes the corpus
policies consumed by `tools.agent_search`.

## Corpus boundaries

The code corpus includes owned source, tests, example configs, durable
physics/API/architecture documentation, and operating guidance. It excludes
working plans/specs, generated outputs, build directories, fixtures containing
deliberate defects, legacy reference trees, and literature notes. Search those
locally when the task requires them.

Processed scientific articles belong to `qscat_articles`, with a separate
connection and provenance contract. Start with tracked `reference/literature/*.md`
notes; add processed full text only from an explicit source root with provenance.
A full-text extraction is not automatically a verified reference claim. The
main-only code contract does not ingest article files or populate the article table.

## The committed-main source contract

Every published code snapshot derives from one exact commit fetched from the
configured upstream's `main`. The default remote is `origin`; `--remote upstream`
selects another configured remote for a fork. Local `main`, the current branch,
staged files, dirty files, and untracked files are never indexing inputs.

The writer performs a bounded, noninteractive fetch of
`refs/heads/main` into `refs/remotes/<remote>/main`, resolves the commit and tree,
enumerates regular blobs with `git ls-tree`, and reads source through batched
`git cat-file` calls. Symlinks and submodules are excluded. File bytes do not pass
through checkout filters. The indexing profile is read from that same commit;
branch policy edits cannot change the shared corpus. A successful fetch is
required before publishing, so a stale offline clone cannot overwrite a newer
shared snapshot with its cached remote-tracking ref.

Fetching updates the remote-tracking ref, not the task's branch, staging area,
or working files. Published provenance records the upstream identity, main commit,
tree, profile hash, selected blob IDs and file SHA256s, chunker, and model spec.
Returned hits carry path, qualified anchor, declaration disambiguator, file/blob
hashes, and snapshot-specific line spans; their response names the main commit.

## Shared identity and publication

The default logical corpus is
`qscat_knowledge_main_<16-character identity SHA256>`. Its identity contains the
normalized upstream URL, `main`, and repository-relative committed profile path.
SSH and HTTPS URLs for the same upstream normalize to the same identity;
credentials are excluded. GitHub repository-name case is normalized. Forks and
other upstreams remain separate. Absolute checkout paths are not owners.

Clones and branches with the same upstream share a logical corpus and writer
lock. `--table` explicitly selects another logical name, retaining ownership
checks; use a new name when changing the embedding model/configuration. Ownership
rejects different upstream/profile identities and legacy manifest formats.

Each changed corpus is built as an immutable physical generation named
`<logical-name>_g_<unique-id>`. The writer constructs rows and the native BM25
index, validates table version and row count, then atomically replaces
`<logical-name>.manifest.json`. That manifest is the publication pointer.
Queries capture it once and pin its physical table version. A failed build or
interrupted publication leaves the previous completed snapshot available. Readers
in other processes can finish using their captured older generation.

Generations are retained because another mounted reader may still hold one.
Failed builds can leave unpublished generations, also excluded from search.
Reclaim obsolete generations only with readers disconnected and the published
manifest preserved; automatic cross-process generation reclamation is not supplied.

## Mount and ongoing refresh

`opencode.json` mounts the project-local reader with:

```bash
uv run --project .opencode/search --python 3.12 --frozen \
  python -m tools.agent_search.mcp
```

The frozen lockfile installs the isolated search environment in a fresh clone.
The startup allowance covers dependency setup. The MCP lifespan starts the
ensure-main job in the background so its handshake and tool catalog can become
available before embedding and indexing finish.

Every mount ensures main is indexed, including when the corpus is populated.
While mounted, the server checks upstream every 300 seconds; `--check-interval`
sets a positive interval in seconds. The shared writer lock and attempt timestamp
rate-limit concurrent mounts. Explicit CLI `sync` forces an immediate fetch and
reconciliation. This is a Git-upstream refresh loop, not a filesystem watcher.

With no prior snapshot, the first query waits for initial population or a reported
failure. With a completed snapshot, queries remain available while a refresh runs
or fails. No complete snapshot means unavailable semantic retrieval, so use local
source. A fetch failure exposes unknown upstream currency; a failed build after
observing a newer main exposes the old snapshot as behind. Neither is called current.

`table_details` and query responses distinguish:

| Field | Meaning |
|---|---|
| `source_commit` | Exact main commit represented by the returned snapshot |
| `integrity` | `complete`, `absent`, or `incomplete` published generation |
| `upstream.state` | `current`, `behind`, or `unknown` relative to a recent successful upstream check |
| `upstream.checked_at` | When the upstream commit was successfully observed |
| `refresh.state` | `pending`, `checking`, `indexing`, `complete`, or `failed` |
| `refresh.error` | Visible failure reason, when present |
| `state` | Summary: `current`, `stale`, `absent`, or `incomplete` |

The inspection tool's `maintainer_error` also exposes local failures when the
writer cannot create a persisted refresh record, such as an invalid database path.

`stale` can therefore mean a complete, usable older main snapshot or unknown remote
currency. Local branch divergence does not change main-index freshness. Once a
successful check ages beyond the interval, currency becomes unknown until checked
again. `status` and CLI queries inspect this information without fetching; the
mounted maintainer and explicit `sync` own refreshes.

## Commands and migration

Run from a checkout with the intended upstream configured:

```bash
uv sync --project .opencode/search --python 3.12
uv run --project .opencode/search --no-sync python -m tools.agent_search plan
uv run --project .opencode/search --no-sync python -m tools.agent_search sync
uv run --project .opencode/search --no-sync python -m tools.agent_search status
uv run --project .opencode/search --no-sync python -m tools.agent_search query \
  "sparse factorization symbolic analysis reuse" --mode hybrid --limit 5
```

`plan` fetches main and reports the prospective corpus without publishing it.
`bootstrap` is a compatibility alias for the same ensure-main operation as `sync`;
it refreshes populated corpora. `--root`, `--remote`, `--db`, `--table`, and
`--profile` provide explicit selections. The profile must be a committed,
repository-relative path. The database defaults to
`~/.local/share/opencode/lancedb`, or `LANCEDB_URI`. The MCP also accepts
`TABLE_NAME` for an explicit logical name. There is no untracked-file ingestion flag.

The main-only manifest format is `qmodeling-main-search-v2`. First migration builds
source membership and provenance afresh from fetched main. A healthy legacy
`qmodeling-search-v1` table from the same upstream can supply self-consistent
payload/model-keyed vectors. Legacy checkout locations, file membership, and
ownership are not copied into the new snapshot. Previously indexed branch content
can only contribute a cached vector when the exact payload/model pair is also
required by committed main.

## Agent workflow for branch work

Use main retrieval to discover candidate files, then resolve each anchor in current
source. Indexed line numbers are hints for the named commit, not current edit
locations. Include local branch additions and changes in discovery:

```bash
git diff --name-status <indexed-commit> -- <scope>
git diff --cached --name-status -- <scope>
git diff --name-status -- <scope>
git ls-files --others --exclude-standard -- <scope>
```

The first diff covers committed branch changes and tracked working-tree changes
relative to the indexed main commit. The next two expose staging and dirty state;
the last lists relevant untracked source. Read relevant diffs and source bodies.
Resolve renames/moves locally and omit deleted hits. If the named commit is missing
in this clone, fetch it through the configured upstream or report the unresolved
delta and continue with local discovery.

After editing, re-read affected source and review dependencies/tests. An unchanged
symbol's embedding says nothing about changed imports or callees. At handoff,
record the indexed commit, its main-index status, and relevant local divergence.
Branches do not refresh their own edits into the shared index. Read-only specialists
retain their assigned scope/report contract and leave explicit maintenance to the
caller. The basic append-only `ingest_docs` tool cannot refresh this corpus.

## Chunking, identity, and embedding reuse

Python AST functions/methods and class/module context, Rust Tree-sitter items, and
Markdown headings provide qualified anchors. Other text/config formats use bounded
file fragments. Repeated declarations and headings have a disambiguator; oversized
units are split while keeping their parent anchor. Empty files appear in the
manifest without searchable chunks.

Three keys have different jobs:

| Key | Depends on | Used for |
|---|---|---|
| Logical document anchor | Repository, path, language/kind, qualified anchor, declaration | Resolving a unit in source |
| Chunk/content key | Exact content and deliberately embedded context | Reconciling fragments |
| Embedding cache key | Complete payload and model revision/parameters | Reusing vectors |

Commits, blob/file hashes, timestamps, and line ranges are metadata outside the
embedded payload. If a helper moves `SparseLU.refactor` down several lines, the
writer reparses that changed main file, updates locations/hashes, and reuses the
method's vector if its full payload is unchanged. Renames/moves retire the old
locations; path context inside the payload can legitimately require new vectors.
Unchanged blob IDs skip parsing. A main commit changing only excluded content
advances commit provenance without rewriting rows or rebuilding the generation.

The CPU baseline is normalized `all-MiniLM-L6-v2`, with 384 dimensions and a
256-token limit. Budget with the installed model's tokenizer: target 160 content
tokens, cap at 220, and require the complete payload including headers and special
tokens to fit 256. Use up to 24 overlap tokens only when splitting oversized units.
These are starting settings, not demonstrated optimal chunk sizes.

The code MCP and CLI support vector, native BM25, and reciprocal-rank-fused hybrid
search. Lexical indexing starts with stemming and stop-word removal disabled,
preserving identifiers and small physics variables. Vector queries use cosine
distance and exact scans. Results are bounded to 1–50 and omit embeddings; local
exact search remains necessary for punctuation-sensitive identifiers and formulas.

## Verification and evaluation

The isolated environment exercises real parsers, embeddings, LanceDB writes, local
Git remotes, and stdio MCP connections:

```bash
uv run --project .opencode/search --no-sync python -m pytest tests/test_agent_search.py -q
```

Checks cover branch/staged/dirty/untracked exclusion, committed policy, main
advancement, unchanged refreshes, line-shift vector reuse, edits/deletions/renames,
shared-clone identity, atomic publication, pinned older readers, interrupted builds,
upstream outages, external mutation recovery, legacy migration, and periodic and
concurrent MCP mounts. The numerical workspace can run parser-only checks without
installing search dependencies.

Use known repository retrieval questions: outgoing-flux DA extraction, ECS
c-product, sparse symbolic reuse, resolved-config TD packet round trips, and grid
convergence versus proxy evidence. Measure relevant-source recall, stale/superseded
hit rate, locator accuracy, latency, result tokens, and total memory. Compare local,
vector, FTS, and hybrid retrieval on the same set. Article evaluation remains
separate and uses verified printed-page locators.

An ANN index is an optimization with a recall cost. Benchmark exact scans first;
add ANN only for measured latency/memory need and compare recall against exact
results. LanceDB 0.21.2 supports `IVF_FLAT`, `IVF_PQ`, `IVF_HNSW_SQ`, and
`IVF_HNSW_PQ`; newer names may not exist in the installed SDK. Keep metrics aligned
and reembed into a new logical corpus when model/configuration changes.

## Implementation research

These primary sources support the design; they are not qModeling retrieval results:

- [CocoIndex's LanceDB code example](https://github.com/cocoindex-io/cocoindex/tree/main/examples/code_embedding_lancedb)
  uses Tree-sitter-aware splitting, separate code/path/span/vector fields, and
  content-dependent IDs with managed incremental writes.
- CocoIndex documents [memoization keys](https://cocoindex.io/docs/advanced_topics/memoization_keys/),
  [ID generation](https://cocoindex.io/docs/common_resources/id_generation/), and
  [LanceDB reconciliation](https://cocoindex.io/docs/connectors/lancedb/).
- Continue's [LanceDbIndex](https://github.com/continuedev/continue/blob/main/core/indexing/LanceDbIndex.ts)
  and [CodebaseIndexer](https://github.com/continuedev/continue/blob/main/core/indexing/CodebaseIndexer.ts)
  separate cached vectors from workspace/branch selection and compute/add/remove/delete
  changes. Content reuse remains useful even though our source policy is main-only.
- [Aider's repository map](https://aider.chat/docs/repomap.html) combines signatures
  and dependency ranking. Use the measured code map for structural/impact facts;
  vector similarity cannot establish caller coverage.
- LanceDB's [table updates](https://docs.lancedb.com/tables/update) and
  [reindexing](https://docs.lancedb.com/indexing/reindexing) distinguish data changes
  from FTS/ANN maintenance. Index optimization alone cannot find deleted sources.

Current CocoIndex's LanceDB extra requires Python >=3.11 and LanceDB >=0.34.
The repository reader/writer uses Python 3.12 with LanceDB 0.21.2, while the
external article connector uses Python 3.10. Upgrade coordinated readers/writers
before adopting that integration. `cocoindex-code` currently uses SQLite, so
installing it is not a LanceDB configuration. No CocoIndex writer is installed.

## Article connection

The basic external connector used for articles takes `LANCEDB_URI`, `TABLE_NAME`,
`EMBEDDING_FUNCTION`, and `MODEL_NAME`. It appends actual text strings into a
`doc`/`vector` table; a supplied path or URL is embedded literally, not read.
Its query mode and inspection selection arguments are ignored: it is vector-only
and selects its configured table. Connection alone proves neither populated data
nor working hybrid retrieval.

A separate OpenCode V2 connection can use:

```jsonc
{
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

Keep workstation paths in machine-local environment/configuration. Use
`opencode mcp list` and the runtime tool catalog to verify connection and capabilities.
Read-only agents need query/inspection access and Code Mode access, not ingestion.
Each connector loads a model process; measure total memory with numerical workloads.
Metadata filters, corpus multiplexing, and article ingestion remain separate work.

## Software documentation

- [OpenCode V2 MCP configuration](https://opencode.ai/v2/docs/mcp-servers)
- [LanceDB hybrid search](https://docs.lancedb.com/search/hybrid-search)
- [LanceDB full-text search](https://docs.lancedb.com/search/full-text-search)
- [LanceDB vector indexes](https://docs.lancedb.com/indexing/vector-index)
- [Tracked literature and locator conventions](../reference/literature/README.md)
