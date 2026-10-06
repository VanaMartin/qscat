# UKRmol+ CPU images

Build and run the external UKRmol+ scattering engine on a Docker-capable Linux
x86-64 host. Build context is the **repository root**. The source engine uses
the upstream MPI/compiler/BLAS/ScaLAPACK/Psi4 toolchain, separate from QSCAT's
uv/maturin application image because these external Fortran interfaces need
their compatible toolchain.

## Targets

| Target | Contents / default command |
|---|---|
| `pilot` | Digest-pinned reference engine plus Python/Perl dependencies; `bash` |
| `build` | UKRmol-in/out and bundled GBTOlib compiled and installed |
| `test` | Build plus 12 selected serial/two-rank upstream reference checks |
| `runtime` | Generic source-built **double** engine, copied only after `test`; `bash` |
| `co-pilot` | CO deck and Python tools on `pilot` |
| `co-source` (default) | CO deck and Python tools on `runtime` |

The generic targets contain no molecule-specific input. Experiment sources are
copied after the engine build/test stages, so changing the CO deck does not
recompile Fortran. Add another experiment layer or bind-mount its source onto
the generic image. The CO experiment is documented in
[`projects/ukrmol_co/README.md`](../../projects/ukrmol_co/README.md).

```bash
# Reference engine; no source build.
docker build --target pilot -f docker/ukrmol-plus/Dockerfile \
  -t qmodeling/ukrmol-plus:reference .

# Source engine; test is a dependency, including for runtime-only builds.
docker build --target runtime --build-arg BUILD_JOBS=4 \
  -f docker/ukrmol-plus/Dockerfile -t qmodeling/ukrmol-plus:source .

# Included CO experiment, source engine (default target).
docker build --build-arg BUILD_JOBS=4 -f docker/ukrmol-plus/Dockerfile \
  -t qmodeling/ukrmol-co:source .

# CO experiment, upstream reference binaries.
docker build --target co-pilot -f docker/ukrmol-plus/Dockerfile \
  -t qmodeling/ukrmol-co:reference .
```

`BUILD_JOBS` bounds compilation concurrency. The upstream gate uses one test
worker; each parallel deck launches two MPI ranks. BuildKit mounts a 1-GiB
`/dev/shm` tmpfs for that gate. These target-HF comparisons validate the selected
engine paths, not every engine feature or the accuracy of the CO model.

## Pinned inputs and build additions

| Input | Pin |
|---|---|
| Toolchain/Psi4 | `zdenekmasin/ukrmol_plus@sha256:1cef2ee6aeaca2a4d5b813ff31e7d28bef4436a4282398c8933a1c047a8ec76c` |
| UKRmol-in 3.3.0 | [Zenodo archive](https://zenodo.org/records/18340340/files/ukrmol-in-3.3.0.tar.gz?download=1), SHA256 `1b81e423b7e2a5a530c336cfc7f7ab50d6f2705abf771c71d6119fc065d9cac2` |
| UKRmol-out 3.3.0.1 | [Zenodo archive](https://zenodo.org/records/18538198/files/ukrmol-out-3.3.0.1.tar.gz?download=1), SHA256 `6fecd651d1b953004926bcd2485d18a1645f2a13fca93e176458b0d48f556a9b` |
| GBTOlib | Compatible sources bundled in the pinned UKRmol-in archive |
| Python additions | NumPy 2.3.4, SciPy 1.16.2, installed without dependency re-resolution |

`source-build.sh` downloads/checks both archives and builds via CMake with MPI,
64-bit default Fortran integers, the upstream numerical libraries and a runtime
library search path matching `/opt/ukrmolp/lib.double`. The outer archive's
internal root differs from its versioned filename; extraction strips that root.

The dependency layer adds full Perl (`Storable`, `Time::HiRes`), `numdiff` and
the system CA certificate bundle. `SSL_CERT_FILE` points the relocated Python
OpenSSL runtime at that bundle; its compiled-in certificate path refers to
the upstream build machine. This is needed for a cold-cache Zenodo download.
The build layer adds `curl`, certificates, `make`, `binutils` and `libc6-dev`,
and sets `OMPI_FC=/opt/compiler/bin/gfortran` and `OMPI_CC=/opt/compiler/bin/gcc`
for the relocated MPI wrappers. Reference **double** binaries/libraries are
removed before copying the installed source engine. **Quad remains the upstream
reference engine**; the double build is not a quad build.

Engine source pins and numerical comparisons make the workflow repeatable.
apt repositories and pip wheels are not fully mirrored/locked; rebuild image
IDs may differ. Record the resulting image ID and package/toolchain provenance.
The runtime retains the upstream toolchain/Psi4 base and is not a minimal image
built wholly from source. Use upstream archive/image licences for their contents;
the repository licence covers the added workflow files.

## Remote deployment

Choose `REMOTE_HOST`, `REMOTE_REPO` and a persistent `RUN_ROOT` under `/home`.
Clone or sync this checkout to the host, then run the build commands there.
Inspect physical-core topology and available resources before choosing CPUs:

```bash
ssh "$REMOTE_HOST" 'lscpu -e=CPU,CORE,SOCKET,NODE; free -h; docker version'
ssh "$REMOTE_HOST" "df -h '$RUN_ROOT'"
ssh "$REMOTE_HOST" "cd '$REMOTE_REPO' && docker build --target runtime \
  --build-arg BUILD_JOBS=4 -f docker/ukrmol-plus/Dockerfile \
  -t qmodeling/ukrmol-plus:source ."
```

On the host, prepare writable bind-mount directories for the image's `ukrmol`
user. Use the actual UID/GID from the image rather than assuming it matches the
host login. For a newly allocated artifact directory:

```bash
IMAGE=qmodeling/ukrmol-plus:source
ENGINE_UID=$(docker run --rm --entrypoint id "$IMAGE" -u)
ENGINE_GID=$(docker run --rm --entrypoint id "$IMAGE" -g)
sudo install -d -o "$ENGINE_UID" -g "$ENGINE_GID" \
  "$RUN_ROOT" "$RUN_ROOT/runs" "$RUN_ROOT/upstream"

docker run --rm --cpuset-cpus=0-3 --memory=16g --memory-swap=16g \
  --shm-size=1g -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
  -e MKL_NUM_THREADS=1 -v "$RUN_ROOT:/work" "$IMAGE" \
  bash -c 'mkdir /work/runs/example-001; cd /work/runs/example-001; \
    mkdir scratch; export TMPDIR="$PWD/scratch" PSI_SCRATCH="$PWD/scratch"; \
    command -v psi4; command -v scatci'
```

This example checks the initialized environment and writable persistent scratch.
Run the experiment's driver for an actual calculation. Use fresh names for each
attempt. `--cpuset-cpus` denotes host CPU IDs, not MPI ranks; one thread per
library avoids oversubscription when the driver launches MPI.

## Provenance and differential checks

```bash
docker image inspect --format '{{.Id}}' qmodeling/ukrmol-plus:source
docker run --rm qmodeling/ukrmol-plus:source \
  cat /opt/ukrmolp/source-provenance.json
docker run --rm qmodeling/ukrmol-plus:source \
  cat /opt/ukrmolp/source-tests/Temporary/LastTest.log
```

The original engine passed all 12 selected checks. The saved CO differential
check has identical printed resonance position/width and a maximum phase
difference of `1e-7 rad` over 491 energies. Its original image ID and source/test
evidence are in
[`source-build-results.json`](../../projects/ukrmol_co/source-build-results.json).
Revalidate a changed build with its own image ID; a recorded historical image
ID is not a downloadable registry reference.

The packaged targets were also rebuilt independently: 12/12 upstream checks
passed, generic `runtime` and reference `co-pilot` smoke checks passed, and
`co-source` completed a cold-cache CO calculation with the same printed
position/width and `1e-7 rad` maximum phase difference. The acquisition check
exposed and resolved the relocated Python CA-path issue. See
[`packaging-results.json`](../../projects/ukrmol_co/packaging-results.json).
