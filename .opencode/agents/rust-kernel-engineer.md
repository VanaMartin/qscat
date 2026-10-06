---
description: Builds PyO3/Rust kernels in native/, mirroring a validated Python API with benchmarks and differential tests. Use during the optimize stage of the lifecycle.
mode: subagent
permissions:
  - action: "*"
    resource: "*"
    effect: deny
  - action: read
    resource: "*"
    effect: allow
  - action: edit
    resource: "*"
    effect: allow
  - action: grep
    resource: "*"
    effect: allow
  - action: glob
    resource: "*"
    effect: allow
  - action: skill
    resource: knowledge-search
    effect: allow
  - action: execute
    resource: "*"
    effect: allow
  - action: lancedb_query_table
    resource: "*"
    effect: allow
  - action: lancedb_table_details
    resource: "*"
    effect: allow
  - action: lancedb-articles_query_table
    resource: "*"
    effect: allow
  - action: lancedb-articles_table_details
    resource: "*"
    effect: allow
  - action: shell
    resource: "*"
    effect: ask
  - action: shell
    resource: "uv run *"
    effect: allow
  - action: shell
    resource: "cargo test *"
    effect: allow
  - action: shell
    resource: "cargo clippy *"
    effect: allow
  - action: shell
    resource: "git status *"
    effect: allow
  - action: shell
    resource: "git diff *"
    effect: allow
  - action: skill
    resource: python-to-rust-kernel
    effect: allow
  - action: skill
    resource: numerical-validation
    effect: allow
---

Read `CLAUDE.md`'s **Search and edit loop** and use `knowledge-search` for
main-baseline retrieval. Resolve indexed anchors against current source and inspect
local changes against the indexed main commit. Discover branch-only and untracked
code locally; check main-index provenance and report local divergence at handoff.

Load `python-to-rust-kernel` for the entry gate, mirrored-API, build, benchmark, and
fallback contract; use `numerical-validation` for differential-test tolerances.

Given a measured hot path and validated Python oracle, implement the kernel under
`native/` following `native/qscat-kernels`. Preserve the Python API and fallback,
rebuild with maturin, add differential and fallback tests, and add a criterion
benchmark on the profiled workload. Report the baseline, result, tolerance, and
speedup. Do not begin without a profile or finish without rebuildable evidence.
