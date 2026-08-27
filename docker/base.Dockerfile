# syntax=docker/dockerfile:1
# qmodeling-base — architecture / BLAS-FFT vendor layer + toolchains.
# Rebuilt rarely. Default vendor: OpenBLAS + LAPACK(E) + FFTW3 (portable, ARM-friendly).
# Code targets the standard CBLAS / LAPACKE / FFTW3 ABIs, so this layer is swappable.
# Tag kept for humans, digest is what builds: `docker buildx imagetools
# inspect ghcr.io/astral-sh/uv:python3.12-bookworm-slim` to refresh, and
# dependabot's docker ecosystem watches it.
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim@sha256:e5b65587bce7de595f299855d7385fe7fca39b8a74baa261ba1b7147afa78e58

# System numerical libraries (standard ABIs) + build toolchain.
# MUMPS (qscat.linalg's optional backend, complex-symmetric SYM=2) is NOT taken
# from Debian: `libmumps-seq` has no OpenMP (issue #4), so it is built from
# source below with OpenMP. We keep libscotch-dev (its fill-reducing ordering
# library) and add gfortran/cmake/git to build MUMPS. METIS is not used.
# ffmpeg is the matplotlib FFMpegWriter backend for qscat.viz's .mp4 animation
# output (the .gif path uses pillow, no system dep); it flows to the `build`/
# `test` stages so the ffmpeg-gated .mp4 viz test actually renders here.
RUN apt-get -o Acquire::Retries=5 update \
    && apt-get -o Acquire::Retries=5 install -y --no-install-recommends \
      build-essential gfortran cmake git curl ca-certificates pkg-config \
      libopenblas-dev libopenblas-openmp-dev libblas-dev liblapack-dev \
      liblapacke-dev libfftw3-dev \
      libscotch-dev \
      ffmpeg \
    && rm -rf /var/lib/apt/lists/*
# libblas-dev/liblapack-dev provide the generic reference BLAS/LAPACK dev
# symlinks that the MUMPS-from-source CMake find_package(BLAS/LAPACK) resolves at
# BUILD time (Debian's libmumps-seq-dev used to pull these in transitively; we
# dropped it). MUMPS links the libblas.so.3 / liblapack.so.3 SONAMEs, which at
# RUNTIME resolve to the OpenMP OpenBLAS selected below via update-alternatives
# (verified: the built libzmumps.so's NEEDED libopenblas.so.0 is the OpenMP one).

# Issue #3: leverage OpenMP-threaded BLAS/LAPACK. Debian ships pthread/openmp/
# serial OpenBLAS variants; the metapackage defaults to pthread. Select the
# OPENMP variant for the RUNTIME SONAMEs (libblas.so.3 / liblapack.so.3 /
# libopenblas.so.0) so at load time BLAS/LAPACK (numpy/scipy dense ops, SuperLU,
# and MUMPS's frontal factorizations) share the SAME libgomp thread pool as the
# OpenMP MUMPS (issue #4). We deliberately do NOT repoint the unversioned dev
# `.so` symlinks: they are what CMake find_library/find_package(BLAS) resolves
# at BUILD time, and leaving them at the default keeps the MUMPS-from-source
# BLAS discovery working; the OpenMP variant still loads at runtime via the
# shared SONAME.
RUN triple="$(gcc -dumpmachine)"; \
    for a in libblas.so.3 liblapack.so.3 libopenblas.so.0; do \
      update-alternatives --set "${a}-${triple}" "/usr/lib/${triple}/openblas-openmp/${a}"; \
    done

# Issue #4: build an OpenMP-enabled SEQUENTIAL MUMPS from source (Debian's
# libmumps-seq has no OpenMP). The scivision/mumps CMake wrapper FetchContent-
# downloads the MUMPS source and builds it non-MPI + OpenMP + SCOTCH ordering,
# linked against the OpenMP OpenBLAS selected above, as shared libs into
# /usr/local. All four arithmetics (s,d,c,z) are built; qscat's complex-
# symmetric ECS backend uses zmumps. libzmumps.so then NEEDs libgomp +
# libopenblas, so MUMPS shares the ONE OpenMP thread pool. Threads follow
# OMP_NUM_THREADS; ICNTL(48) tree parallelism defaults on in an OpenMP MUMPS.
RUN git clone --depth 1 https://github.com/scivision/mumps.git /tmp/mumps \
    && cmake -S /tmp/mumps -B /tmp/mumps/build \
         -DMUMPS_parallel=no -DENABLE_OPENMP=on -DENABLE_SCOTCH=on \
         -DBLAS_VENDOR=OpenBLAS -DBUILD_SHARED_LIBS=on \
         -DBUILD_SINGLE=on -DBUILD_DOUBLE=on -DBUILD_COMPLEX=on -DBUILD_COMPLEX16=on \
         -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release \
    && cmake --build /tmp/mumps/build -j"$(nproc)" \
    && cmake --install /tmp/mumps/build \
    && ldconfig \
    && rm -rf /tmp/mumps

# python-mumps discovers MUMPS via pkg-config using the conda-forge names
# {d,z,c,s}mumps_seq. Synthesize those .pc files pointing at the from-source
# OpenMP libs in /usr/local (/usr/local/lib/pkgconfig is on pkg-config's
# default path; /usr/local/lib is on ldconfig's). Each libNmumps.so records its
# deps (mumps_common, SCOTCH, OpenBLAS, libgomp, gfortran) via ELF NEEDED.
RUN mkdir -p /usr/local/lib/pkgconfig; \
    for name in dmumps zmumps cmumps smumps; do \
      printf 'prefix=/usr/local\nlibdir=${prefix}/lib\nincludedir=${prefix}/include\n\nName: %s_seq\nDescription: OpenMP MUMPS (scivision/mumps source build), %s\nVersion: 5.9.1\nLibs: -L${libdir} -l%s -lmumps_common\nCflags: -I${includedir}\n' \
        "$name" "$name" "$name" \
        > "/usr/local/lib/pkgconfig/${name}_seq.pc"; \
    done

# Rust toolchain (compiles qscat kernels)
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Issue #3: threading is ENABLED by default (no hard OPENBLAS_NUM_THREADS=1
# pin). Cap per-run with OMP_NUM_THREADS / OPENBLAS_NUM_THREADS for deterministic
# tests or to avoid oversubscription once process-level parallelism is added.

# Fail the build if the standard ABIs are not discoverable (sanity gate).
# Includes MUMPS: the {d,z}mumps_seq .pc files must resolve so python-mumps can
# build against the system library.
RUN pkg-config --exists openblas && pkg-config --exists fftw3 && pkg-config --exists lapacke \
    && pkg-config --exists dmumps_seq && pkg-config --exists zmumps_seq \
    && echo "openblas $(pkg-config --modversion openblas), fftw3 $(pkg-config --modversion fftw3), lapacke $(pkg-config --modversion lapacke), mumps_seq $(pkg-config --modversion dmumps_seq)"
