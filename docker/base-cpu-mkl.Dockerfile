# syntax=docker/dockerfile:1
# qmodeling-base:cpu-mkl — x86-64 CPU-OPTIMIZED base (Intel MKL variant).
#
# STATUS: SCAFFOLD. This provisions Intel MKL + the build toolchain + system
# MUMPS so the layer is ready, but two pieces of wiring are deliberately left as
# tracked optimization work (roadmap Part 5), NOT claimed working here:
#   1. numpy/scipy must actually LINK MKL (the PyPI wheels bundle OpenBLAS). The
#      supported route is installing numpy/scipy from a MKL-linked channel
#      (conda-forge with `blas=*=mkl`, or Intel's channel) instead of the pip
#      wheels `uv sync` pulls. Until then this image builds/runs but uses the
#      bundled OpenBLAS, i.e. it is not yet faster than the portable base.
#   2. An MKL PARDISO backend for qscat.linalg.SparseLU (the eMoScat reference
#      solver) is the Part-5 hot-path target; the MUMPS backend below is the
#      current fast path.
# x86-64 only (MKL has no ARM build); the portable base.Dockerfile stays the
# ARM/Graviton path.
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim

# Base build toolchain + FFTW3 + system MUMPS (same as the portable base).
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential curl ca-certificates pkg-config gnupg \
      libfftw3-dev libmumps-seq-dev libscotch-dev \
    && rm -rf /var/lib/apt/lists/*

# Intel oneAPI MKL from Intel's apt repository (devel = headers + libs + .pc).
RUN curl -fsSL https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
      | gpg --dearmor -o /usr/share/keyrings/oneapi-archive-keyring.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
      > /etc/apt/sources.list.d/oneAPI.list \
    && apt-get update && apt-get install -y --no-install-recommends \
      intel-oneapi-mkl-devel \
    && rm -rf /var/lib/apt/lists/*
ENV MKLROOT=/opt/intel/oneapi/mkl/latest
ENV LD_LIBRARY_PATH="${MKLROOT}/lib/intel64:${LD_LIBRARY_PATH:-}"
ENV PKG_CONFIG_PATH="${MKLROOT}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

# Synthesize the MUMPS pkg-config .pc files Debian omits (see base.Dockerfile).
RUN triple="$(gcc -dumpmachine)"; \
    ver="$(dpkg-query -W -f='${Version}' libmumps-seq-dev | sed -E 's/^[0-9]+://; s/[-+].*$//')"; \
    mkdir -p "/usr/lib/${triple}/pkgconfig"; \
    for name in dmumps zmumps cmumps smumps; do \
      printf 'prefix=/usr\nlibdir=${prefix}/lib/%s\nincludedir=${prefix}/include\n\nName: %s_seq\nDescription: MUMPS sequential (Debian libmumps-seq-dev), %s\nVersion: %s\nLibs: -L${libdir} -l%s_seq\nCflags: -I${includedir}\n' \
        "$triple" "$name" "$name" "$ver" "$name" \
        > "/usr/lib/${triple}/pkgconfig/${name}_seq.pc"; \
    done

# Rust toolchain (compiles qscat kernels)
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# MKL threading: 1 by default; tune per workload.
ENV MKL_NUM_THREADS=1 OMP_NUM_THREADS=1

# Sanity gate: MKL + FFTW3 + MUMPS pkg-config must resolve.
RUN pkg-config --exists mkl-dynamic-lp64-seq && pkg-config --exists fftw3 \
    && pkg-config --exists dmumps_seq && pkg-config --exists zmumps_seq \
    && echo "mkl $(pkg-config --modversion mkl-dynamic-lp64-seq), mumps_seq $(pkg-config --modversion dmumps_seq)"
