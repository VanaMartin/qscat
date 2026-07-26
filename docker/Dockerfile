# syntax=docker/dockerfile:1
# Build/test/runtime for qModeling. Architecture/vendor lives in BASE_IMAGE
# (see docker/base.Dockerfile). Swap BASE_IMAGE to retarget (e.g. an -mkl or ARM base)
# without editing this file.
ARG BASE_IMAGE=qmodeling-base:latest

FROM ${BASE_IMAGE} AS build
WORKDIR /app
COPY . .
# --all-packages installs qscat AND builds/installs the Rust qscat_kernels.
# NOTE: the mumps extra is deliberately NOT installed here -- `runtime` copies
# from this stage, and shipping python-mumps into the production image (whose
# apt layer has no system MUMPS .so) would bloat it with a nonfunctional
# package. The extra is added only in the `test` stage below.
RUN uv sync --all-packages

FROM build AS test
# Add the mumps extra ONLY here so the complex-symmetric MUMPS backend is
# exercised (not skipped) -- python-mumps links the system MUMPS provisioned in
# the base image. This stage is not copied by `runtime`, so the extra stays
# test-only. Then run the suite with --no-sync (no re-resolve/prune of the
# just-built workspace members).
RUN uv sync --all-packages --extra mumps \
    && uv run --no-sync pytest -q

FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS runtime
# Runtime needs only the shared numerical libs (no -dev headers, no toolchain)
RUN apt-get update && apt-get install -y --no-install-recommends \
      libopenblas0 libfftw3-double3 liblapacke \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /app
# LOAD-BEARING: `uv sync --all-packages` installs `qscat` as an editable
# install whose .pth file points at /app/libs/qscat, so the source tree
# (not just .venv) must be copied. Do not reduce this to copying only
# .venv without first switching qscat to a non-editable/wheel install.
COPY --from=build /app /app
ENV PATH="/app/.venv/bin:${PATH}"
CMD ["python", "-c", "import qscat, qscat_kernels; print('qscat', qscat.__version__, 'ready')"]
