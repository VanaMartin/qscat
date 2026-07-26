# syntax=docker/dockerfile:1
# qmodeling-base — architecture / BLAS-FFT vendor layer + toolchains.
# Rebuilt rarely. Default vendor: OpenBLAS + LAPACK(E) + FFTW3 (portable, ARM-friendly).
# Code targets the standard CBLAS / LAPACKE / FFTW3 ABIs, so this layer is swappable.
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim

# System numerical libraries (standard ABIs) + build toolchain.
# libmumps-seq-dev is the sequential (non-MPI) MUMPS direct sparse solver used
# by qscat.linalg's optional MUMPS backend (complex-symmetric SYM=2). Debian's
# sequential MUMPS is compiled against SCOTCH (fill-reducing ordering), so
# libscotch-dev is pulled in as its ordering library; METIS is NOT compiled
# into the Debian seq build, so libmetis-dev is deliberately not installed.
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential curl ca-certificates pkg-config \
      libopenblas-dev liblapacke-dev libfftw3-dev \
      libmumps-seq-dev libscotch-dev \
    && rm -rf /var/lib/apt/lists/*

# Debian's libmumps-seq-dev ships headers + libs but NO pkg-config .pc files,
# whereas python-mumps discovers MUMPS via pkg-config using the conda-forge
# names {d,z,c,s}mumps_seq. Synthesize those .pc files pointing at the Debian
# libs. The multiarch pkgconfig dir is on pkg-config's default search path, so
# no PKG_CONFIG_PATH is needed. Each shared lib records its dependencies
# (mumps_common_seq, SCOTCH, LAPACK, gfortran) via ELF NEEDED, so -l<name>_seq
# alone links transitively.
RUN triple="$(gcc -dumpmachine)"; \
    mkdir -p "/usr/lib/${triple}/pkgconfig"; \
    for name in dmumps zmumps cmumps smumps; do \
      printf 'prefix=/usr\nlibdir=${prefix}/lib/%s\nincludedir=${prefix}/include\n\nName: %s_seq\nDescription: MUMPS sequential (Debian libmumps-seq-dev), %s\nVersion: 5.5.0\nLibs: -L${libdir} -l%s_seq\nCflags: -I${includedir}\n' \
        "$triple" "$name" "$name" "$name" \
        > "/usr/lib/${triple}/pkgconfig/${name}_seq.pc"; \
    done

# Rust toolchain (compiles qscat kernels)
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Sane default; kernels/tests can override
ENV OPENBLAS_NUM_THREADS=1

# Fail the build if the standard ABIs are not discoverable (sanity gate).
# Includes MUMPS: the {d,z}mumps_seq .pc files must resolve so python-mumps can
# build against the system library.
RUN pkg-config --exists openblas && pkg-config --exists fftw3 && pkg-config --exists lapacke \
    && pkg-config --exists dmumps_seq && pkg-config --exists zmumps_seq \
    && echo "openblas $(pkg-config --modversion openblas), fftw3 $(pkg-config --modversion fftw3), lapacke $(pkg-config --modversion lapacke), mumps_seq $(pkg-config --modversion dmumps_seq)"
