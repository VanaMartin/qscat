---
name: containerize-and-run
description: Use when packaging a qModeling capability to run locally in Docker or prepare it for AWS — reproducible multi-stage builds via uv + maturin, CPU-only.
---

# containerize-and-run

## Overview

qModeling ships one multi-stage `docker/Dockerfile` with `build` / `test` /
`runtime` stages, verified working, CPU-only end to end. Containerizing new
work means reusing these stages, not inventing new Docker plumbing.

## When to Use

- Packaging a capability to run locally in Docker.
- Gating a change with a container-based test build before merge.
- Preparing a capability for eventual AWS deployment.

## Reuse the Existing Stages

`docker/Dockerfile` (see `docker/README.md`) defines:

- **`build`** — installs Rust and runs `uv sync --all-packages`, which
  installs `qscat` and builds/installs the Rust `qscat_kernels` extension in
  one step, on top of `ghcr.io/astral-sh/uv:python3.12-bookworm-slim`.
- **`test`** — `FROM build`, runs `uv run --no-sync pytest -q`. `--no-sync`
  is required here so `uv` does not re-resolve and prune the workspace
  members that `build` just compiled.
- **`runtime`** — a fresh slim base with only `/app` (including the built
  `.venv`) copied in from `build`; no compiler toolchain, no Rust source.
  Its default `CMD` imports `qscat` and `qscat_kernels` to prove the image is
  runnable.

For a new capability, add its dependencies to the relevant `pyproject.toml`
(root workspace or the specific package) rather than editing the Dockerfile;
the existing `build` stage picks up anything `uv sync --all-packages`
resolves.

## Reproducibility

- Local development: `uv sync --all-packages` (canonical setup command;
  resolves and installs the whole workspace, including Rust kernels).
- Container builds should prefer `uv sync --frozen` where a lockfile is
  committed, so the image installs exactly the locked versions rather than
  re-resolving — reproducible builds are the point of containerizing at all.
- Kernels built for a release/runtime image should use
  `maturin develop --release` (or `uv run maturin develop --release`) so the
  container ships optimized native code, not a debug build.

## The Gate: `docker build --target test`

Before trusting a container-packaged change, run the `test` stage as a gate:

```bash
docker build --target test -f docker/Dockerfile .
```

This fails the build if `pytest` fails inside the image — treat it the same
as a CI gate. Follow with the runtime build/run to confirm the shippable
image actually starts:

```bash
docker build --target runtime -t qmodeling:latest -f docker/Dockerfile .
docker run --rm qmodeling:latest
```

## AWS Deployment (Deferred)

AWS deployment is not implemented yet. When it is, it extends the `runtime`
stage (same CPU-only image) rather than introducing a parallel build path —
don't design a second Dockerfile or a GPU-oriented base image for it.

## Common Mistakes

- Running `uv sync` (without `--no-sync`) in the `test` stage — it can
  re-resolve and prune packages the `build` stage already compiled.
- Adding a new native crate under `native/` and forgetting it's picked up
  automatically by `uv sync --all-packages` via the workspace — no Dockerfile
  edit should be needed for a new kernel crate that follows the
  `python-to-rust-kernel` pattern.
- Building only the `runtime` target and skipping the `test` target gate
  before considering a container change done.
- Introducing GPU base images or CUDA — this repo is CPU-only by design (see
  `qscat-conventions`); `reference/libXcuda` is a read-only oracle, not a
  runtime dependency.
