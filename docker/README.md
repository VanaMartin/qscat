# docker/

CPU-only images, split into two layers:

1. **`qmodeling-base`** (`docker/base.Dockerfile`) — owns architecture and BLAS/FFT
   vendor choices: OpenBLAS + LAPACK(E) + FFTW3 dev libs, pkg-config, a Rust toolchain,
   and uv/Python 3.12. Rebuilt rarely. Code targets the standard CBLAS / LAPACKE / FFTW3
   ABIs, so this layer is swappable — an `-mkl` (x86-64) or ARM/Graviton-tuned variant can
   replace it later without touching `docker/Dockerfile`, just by pointing `BASE_IMAGE` at
   the new tag.
2. **`docker/Dockerfile`** — `FROM ${BASE_IMAGE}` (default `qmodeling-base:latest`) and
   layers `build` → `test` → `runtime` on top. Setup uses the canonical
   `uv sync --all-packages`, which installs `qscat` and builds the Rust `qscat_kernels`
   in one step. The `runtime` stage starts fresh from the upstream uv/Python image and
   installs only the non-dev shared libs (`libopenblas0`, `libfftw3-double3`,
   `liblapacke`) needed at runtime — no compilers, no `-dev` headers.

## Build

```bash
# Build qmodeling-base, then the app test stage (fails the build if tests fail):
docker/build.sh test

# Build qmodeling-base, then the app runtime image, and run it:
docker/build.sh runtime
docker run --rm qmodeling:runtime
```

`docker/build.sh` always builds `qmodeling-base:latest` first, then the requested app
target (`test` by default), so the base is never stale relative to the app image.

## Running a `qscat-run` experiment

`docker/run.sh` builds the `test` image (has MUMPS — needed for the larger production
decks, e.g. H2P's ~1.15M-unknown full grid) and runs a `qscat-run` YAML config inside it
in one step:

```bash
docker/run.sh CONFIG [OUTPUT_DIR]     # OUTPUT_DIR defaults to runs/<config-stem>

# e.g.
docker/run.sh apps/qscat-run/examples/h2p-dr.yaml runs/h2p-dr
```

It mounts `CONFIG` into the container read-only at `/config.yaml` and `OUTPUT_DIR` out
at `/out`, then runs `uv run --no-sync qscat-run run /config.yaml --output /out`. This
replaces the molecule-specific `docker/run-n2.sh` (kept as-is for the older N2 LCP
benchmark harness) with a general entry point for any `qscat-run` config. See the
top-level `README.md` for the full `qscat-run` walkthrough.

## AWS deploy

AWS deployment pushes `qmodeling-base` to ECR once (it changes only when the
architecture/vendor choice changes) and reuses it across app builds — only the `build`/
`test`/`runtime` layers rebuild on ordinary code changes.
