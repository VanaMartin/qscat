# 3. Layered Docker with an OpenBLAS/FFTW3 vendor base

Date: 2026-07-14
Status: accepted

## Context
qModeling needs reproducible CPU containers for testing and (eventually) AWS
deployment. Two concerns bundled into a single Dockerfile would churn together:
the architecture/BLAS-FFT vendor choice (rarely changes) and the application
build/test/runtime steps (change on every commit). We also want to keep future
options open — an MKL-tuned x86-64 image, or an ARM/Graviton image — without
rewriting the application image.

## Decision
Split the image into two layers:

- `docker/base.Dockerfile` builds `qmodeling-base`: the architecture / BLAS-FFT
  vendor layer. It installs OpenBLAS + LAPACK(E) + FFTW3 dev libraries,
  pkg-config, a Rust toolchain, and uv/Python 3.12. All qModeling code targets
  the standard CBLAS / LAPACKE / FFTW3 ABIs (not a vendor-specific API), so
  this layer is swappable.
- `docker/Dockerfile` is `FROM ${BASE_IMAGE}` (default `qmodeling-base:latest`)
  and layers `build` → `test` → `runtime` on top, using the canonical
  `uv sync --all-packages` for setup.
- `docker/build.sh [test|runtime]` always builds the base first, then the
  requested app target, so the base is never stale relative to the app image.

The default vendor pick is OpenBLAS + LAPACK(E) + FFTW3 rather than Intel MKL,
because it is portable across architectures and keeps an ARM/Graviton base
open as a future variant. An MKL-backed x86-64 base is a planned future
variant, built by pointing `BASE_IMAGE` at a new tag — not implemented in
this phase.

## Consequences
The base image is rebuilt rarely (only on architecture/vendor changes) and can
be cached/pushed to ECR independently of ordinary application code changes,
which rebuild only the `build`/`test`/`runtime` layers. Verified working:
`docker/build.sh test` reports `5 passed`; `docker/build.sh runtime` starts a
container that prints `qscat 0.0.0 ready`.

The tradeoff is an extra layer of indirection (two Dockerfiles, a `BASE_IMAGE`
build-arg) versus a single flat Dockerfile. When an MKL or ARM/Graviton
variant is added, the `runtime` stage's non-dev shared-library set
(`libopenblas0`, `libfftw3-double3`, `liblapacke`) must be tied to whichever
base produced the `build` stage — the runtime stage cannot assume OpenBLAS if
the build stage came from an MKL base — so any additional base variant needs
either a matching runtime stage or an explicit runtime library selection
keyed off the base.
