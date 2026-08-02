# syntax=docker/dockerfile:1
# qmodeling-base:gpu — GPU-OPTIMIZED base (CUDA).  STATUS: SCAFFOLD (deferred).
#
# GPU is an explicit "later" item on the roadmap (docs/superpowers/plans/
# 2026-08-02-qscat-hardening-roadmap.md, Part 4/5). This file exists so the base
# swap seam (docker/Dockerfile's ARG BASE_IMAGE + docker/build.sh gpu) is in
# place and a future GPU sparse backend for qscat.linalg.SparseLU (cuDSS /
# cuSOLVER — the fourth dispatch option alongside SuperLU/MUMPS/PARDISO) slots in
# exactly like MUMPS did, WITHOUT touching the app layers.
#
# No qscat GPU kernels exist yet (reference/libXcuda is the future-porting
# reference). Building/running this today gives a CUDA-capable CPU environment;
# the GPU code paths are TODO. Requires an NVIDIA GPU + nvidia-container-toolkit
# at run time (`docker run --gpus all ...`).
#
# CUDA ships cuBLAS/cuSPARSE/cuSOLVER in the -devel image; cuDSS (the sparse
# direct solver) would be added here when the GPU SparseLU backend lands.
FROM nvidia/cuda:12.4.1-devel-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential curl ca-certificates pkg-config \
      python3.12 python3.12-venv python3.12-dev \
      libfftw3-dev libmumps-seq-dev libscotch-dev \
    && rm -rf /var/lib/apt/lists/*

# uv (matches the CPU bases' Python/uv workflow)
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:${PATH}"

# Rust toolchain (compiles qscat kernels; future GPU kernels may add a CUDA
# build step here).
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# TODO(gpu): synthesize MUMPS .pc files (see base.Dockerfile), add cuDSS, and a
# pkg-config sanity gate — deferred until the GPU SparseLU backend is designed.
