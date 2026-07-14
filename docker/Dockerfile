# syntax=docker/dockerfile:1
# Build/test/runtime for qModeling. Architecture/vendor lives in BASE_IMAGE
# (see docker/base.Dockerfile). Swap BASE_IMAGE to retarget (e.g. an -mkl or ARM base)
# without editing this file.
ARG BASE_IMAGE=qmodeling-base:latest

FROM ${BASE_IMAGE} AS build
WORKDIR /app
COPY . .
# --all-packages installs qscat AND builds/installs the Rust qscat_kernels
RUN uv sync --all-packages

FROM build AS test
# --no-sync so uv does not re-resolve and prune the just-built workspace members
RUN uv run --no-sync pytest -q

FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS runtime
# Runtime needs only the shared numerical libs (no -dev headers, no toolchain)
RUN apt-get update && apt-get install -y --no-install-recommends \
      libopenblas0 libfftw3-double3 liblapacke \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY --from=build /app /app
ENV PATH="/app/.venv/bin:${PATH}"
CMD ["python", "-c", "import qscat, qscat_kernels; print('qscat', qscat.__version__, 'ready')"]
