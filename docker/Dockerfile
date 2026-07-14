# syntax=docker/dockerfile:1
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential curl ca-certificates && rm -rf /var/lib/apt/lists/*
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
WORKDIR /app
COPY . .
# --all-packages installs qscat AND builds/installs the Rust qscat_kernels
RUN uv sync --all-packages

FROM build AS test
# --no-sync so uv does not re-resolve and prune the just-built workspace members
RUN uv run --no-sync pytest -q

FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS runtime
WORKDIR /app
COPY --from=build /app /app
ENV PATH="/app/.venv/bin:${PATH}"
CMD ["python", "-c", "import qscat, qscat_kernels; print('qscat', qscat.__version__, 'ready')"]
