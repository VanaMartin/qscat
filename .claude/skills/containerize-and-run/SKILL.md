---
name: containerize-and-run
description: Use when packaging a qModeling capability to run locally in Docker or prepare it for AWS — reproducible multi-stage builds via uv + maturin, CPU-only.
---

# containerize-and-run

## Overview

qModeling ships **two** Dockerfiles, layered, verified working, CPU-only end
to end. Containerizing new work means reusing these, not inventing new
Docker plumbing.

- **`docker/base.Dockerfile`** builds `qmodeling-base` — the architecture /
  BLAS-FFT vendor layer: OpenBLAS + LAPACK(E) + FFTW3 dev libs, pkg-config,
  **the Rust toolchain**, and uv/Python 3.12, on top of
  `ghcr.io/astral-sh/uv:python3.12-bookworm-slim`. Rebuilt rarely; a
  pkg-config sanity gate fails the build early if the CBLAS/LAPACKE/FFTW3
  ABIs aren't discoverable.
- **`docker/Dockerfile`** is `FROM ${BASE_IMAGE}` (default
  `qmodeling-base:latest`) and layers `build` → `test` → `runtime` on top.

## When to Use

- Packaging a capability to run locally in Docker.
- Gating a change with a container-based test build before merge.
- Preparing a capability for eventual AWS deployment.

## Reuse the Existing Stages

`docker/Dockerfile` (see `docker/README.md`) defines:

- **`build`** — `FROM ${BASE_IMAGE}` (Rust and system numerical libs already
  present there), runs `uv sync --all-packages`, which installs `qscat` and
  builds/installs the Rust `qscat_kernels` extension in one step.
- **`test`** — `FROM build`, runs `uv run --no-sync pytest -q`. `--no-sync`
  is required here so `uv` does not re-resolve and prune the workspace
  members that `build` just compiled. `docker/build.sh test` prints
  `5 passed` when this stage succeeds.
- **`runtime`** — a fresh slim base (not `qmodeling-base`) with only the
  non-dev shared libs installed, plus `/app` (including the built `.venv`
  and source tree) copied in from `build`; no compiler toolchain, no Rust
  source. Its default `CMD` imports `qscat` and `qscat_kernels` and prints
  `qscat 0.0.0 ready` to prove the image is runnable.

Build both layers with `docker/build.sh [test|runtime]` (default `test`) —
it always builds `qmodeling-base:latest` first, then the requested app
target, tagged `qmodeling:test` or `qmodeling:runtime`.

For a new capability, add its dependencies to the relevant `pyproject.toml`
(root workspace or the specific package) rather than editing either
Dockerfile; the existing `build` stage picks up anything
`uv sync --all-packages` resolves. Only touch `base.Dockerfile` when the
architecture/vendor choice itself changes (new system lib, new toolchain).

## Reproducibility

- Local development: `uv sync --all-packages` (canonical setup command;
  resolves and installs the whole workspace, including Rust kernels).
- Container builds should prefer `uv sync --frozen` where a lockfile is
  committed, so the image installs exactly the locked versions rather than
  re-resolving — reproducible builds are the point of containerizing at all.
- Kernels built for a release/runtime image should use
  `maturin develop --release` (or `uv run maturin develop --release`) so the
  container ships optimized native code, not a debug build.

## The Gate: `docker/build.sh test`

Before trusting a container-packaged change, run the `test` stage as a gate:

```bash
docker/build.sh test
```

This builds `qmodeling-base` then the app `test` target, and fails if
`pytest` fails inside the image — treat it the same as a CI gate (expect
`5 passed`). Follow with the runtime build/run to confirm the shippable
image actually starts:

```bash
docker/build.sh runtime
docker run --rm qmodeling:runtime
```

This prints `qscat 0.0.0 ready`.

## AWS Deployment (Deferred)

AWS deployment is not implemented yet. When it is, it extends the `runtime`
stage (same CPU-only image), pushing `qmodeling-base` to ECR once and
reusing it across app builds. The base layer (`base.Dockerfile`) is the
sanctioned way to retarget architecture/vendor choices (e.g. a future `-mkl`
or ARM/Graviton variant) by swapping `BASE_IMAGE` — don't introduce a
GPU-oriented base image, and don't add a third build path outside the
base/build/test/runtime layering for scope that isn't an architecture
change.

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
- Editing `base.Dockerfile` for ordinary code/dependency changes — that
  layer is architecture/vendor only and is rebuilt rarely; regular changes
  belong in `pyproject.toml` and flow through the `build` stage.
