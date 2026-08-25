# syntax=docker/dockerfile:1
# Build/test/runtime for qModeling. Architecture/vendor lives in BASE_IMAGE
# (see docker/base.Dockerfile). Swap BASE_IMAGE to retarget (e.g. an -mkl or ARM base)
# without editing this file.
ARG BASE_IMAGE=qmodeling-base:latest

FROM ${BASE_IMAGE} AS build
# The build context excludes .git (.dockerignore), so `git rev-parse`
# has no repo inside the image and qscat-run's manifest recorded
# git_sha="unknown" for every containerised run. The host passes its SHA
# in instead; qscat_run.artifacts reads QSCAT_GIT_SHA first.
ARG GIT_SHA=unknown
ENV QSCAT_GIT_SHA=${GIT_SHA}
WORKDIR /app
COPY . .
# --all-packages installs qscat AND builds/installs the Rust qscat_kernels.
# The `plot` extra (matplotlib) IS installed here so it flows into `runtime`,
# which copies this stage's venv verbatim (runtime has no toolchain, so it
# cannot re-sync/rebuild the Rust kernel) -- that is what makes qscat.viz's
# animation output usable in the runtime image (paired with ffmpeg below).
# NOTE: the mumps extra is deliberately NOT installed here -- shipping
# python-mumps into the production image (whose apt layer has no system MUMPS
# .so) would bloat it with a nonfunctional package. That extra is test-only.
RUN uv sync --all-packages --extra plot

FROM build AS test-deps
# Add the mumps extra here so the complex-symmetric MUMPS backend is available
# -- python-mumps links the system MUMPS provisioned in the base image. `plot`
# is re-stated because `uv sync` reconciles to exactly the given extras
# (omitting it would prune matplotlib). Neither this stage nor `test` below is
# copied by `runtime`, so the mumps extra stays test-only.
#
# This is the target `docker/run.sh` builds from: it needs the MUMPS backend
# for large production decks (e.g. H2P's ~1.15M-unknown grid) but must NOT pay
# for, or be blocked by, a full test-suite run on every compute invocation --
# that coupling was the reason this stage was split out of `test`.
RUN uv sync --all-packages --extra mumps --extra plot

FROM test-deps AS test
# Executes the suite -- reached only when explicitly requested
# (`docker/build.sh test`), never as a side effect of building `test-deps` or
# `runtime`. --no-sync skips a re-resolve/prune of the just-built workspace
# members; the base image's ffmpeg makes the ffmpeg-gated .mp4 viz test render
# rather than skip.
#
# This is where the `slow` tier runs unconditionally. CI runs it only on
# demand -- a `validate:*` label or a manual dispatch of
# .github/workflows/validation.yml (docs/adr/0005) -- so this image is the one
# place a full slow run is guaranteed. It runs BOTH tiers, but not the same way:
#
#   fast  -- parallel, `--dist loadfile` so a module's grid is built once per
#            file rather than once per worker.
#   slow  -- SERIAL, deliberately. These are the production-scale decks; the
#            H2+ DR example alone was measured at ~19 GB peak RSS, so running
#            them concurrently would OOM a container long before it saved
#            wall-clock.
RUN uv run --no-sync pytest -q -m "not slow" -n auto --dist loadfile
RUN uv run --no-sync pytest -q -m slow

FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS runtime
# Runtime needs the shared numerical libs (no -dev headers, no toolchain) plus
# ffmpeg -- the matplotlib FFMpegWriter backend for qscat.viz's .mp4 animation
# output (matplotlib itself arrives in the copied venv via the build stage's
# `--extra plot`; the .gif path needs only pillow, already a plot-extra dep).
RUN apt-get update && apt-get install -y --no-install-recommends \
      libopenblas0 libfftw3-double3 liblapacke \
      ffmpeg \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /app
# LOAD-BEARING: `uv sync --all-packages` installs `qscat` as an editable
# install whose .pth file points at /app/libs/qscat, so the source tree
# (not just .venv) must be copied. Do not reduce this to copying only
# .venv without first switching qscat to a non-editable/wheel install.
COPY --from=build /app /app
ENV PATH="/app/.venv/bin:${PATH}"
CMD ["python", "-c", "import qscat, qscat_kernels; print('qscat', qscat.__version__, 'ready')"]
