# 6. The Rust kernel crate stays a stub until a proven hot path exists

Date: 2026-08-25

## Status

Accepted

## Context

`native/qscat-kernels` exists to keep the PyO3/maturin toolchain proven
end-to-end (build, import, differential test, Docker), but its only content
is an `l2_norm` placeholder. Meanwhile every CI job paid a Rust toolchain
install and a maturin build for it, and its metadata (no license, edition
2021, pyo3 0.22) misrepresented the project. The measured hot path is the
sparse LU factorization and its per-step triangular solves
(docs/physics/optimization-targets.md), which are already served by the
MUMPS backend — the profile explicitly found NO pure-Python hot loop that a
first Rust kernel could win on the current direct-solver architecture.

## Decision

Keep the crate in the workspace with honest metadata, and stop paying for
it on unrelated CI runs: ci.yml builds the kernel only when `native/**`,
`uv.lock`, or the workflow itself changes (dorny/paths-filter); otherwise
`uv sync --all-packages --no-install-package qscat-kernels` skips the build
and the kernel's tests skip via `importorskip`. Implementing a real kernel
now — the qm-method-lifecycle stage-4 path, e.g. a sparse-LU-adjacent or
propagation-inner-loop kernel — is explicitly deferred until a profile
shows a hot path the MUMPS backend does not already cover.

## Consequences

- CI jobs untouched by `native/**` skip the toolchain + build entirely.
- The Docker images and local `uv sync --all-packages` still always build
  the kernel, so the toolchain never rots unexercised.
- A future kernel starts from a crate whose metadata is already honest.
