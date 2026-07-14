# syntax=docker/dockerfile:1
# qmodeling-base — architecture / BLAS-FFT vendor layer + toolchains.
# Rebuilt rarely. Default vendor: OpenBLAS + LAPACK(E) + FFTW3 (portable, ARM-friendly).
# Code targets the standard CBLAS / LAPACKE / FFTW3 ABIs, so this layer is swappable.
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim

# System numerical libraries (standard ABIs) + build toolchain
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential curl ca-certificates pkg-config \
      libopenblas-dev liblapacke-dev libfftw3-dev \
    && rm -rf /var/lib/apt/lists/*

# Rust toolchain (compiles qscat kernels)
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Sane default; kernels/tests can override
ENV OPENBLAS_NUM_THREADS=1

# Fail the build if the standard ABIs are not discoverable (sanity gate)
RUN pkg-config --exists openblas && pkg-config --exists fftw3 && pkg-config --exists lapacke \
    && echo "openblas $(pkg-config --modversion openblas), fftw3 $(pkg-config --modversion fftw3), lapacke $(pkg-config --modversion lapacke)"
