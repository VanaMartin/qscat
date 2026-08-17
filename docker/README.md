# docker/

Images split into two layers: a rarely-rebuilt **base** (architecture / BLAS-FFT
vendor / GPU) and the **app** (`build` → `test` → `runtime`) that is `FROM` it. The
app layer never changes when you swap the base — it reads `ARG BASE_IMAGE`.

## Base variants

| variant | Dockerfile | what it is | status |
|---|---|---|---|
| `cpu` (default) | `base.Dockerfile` | OpenBLAS + LAPACK(E) + FFTW3 + system MUMPS; ARM/Graviton-friendly | working |
| `cpu-mkl` | `base-cpu-mkl.Dockerfile` | x86-64 Intel MKL variant (CPU-optimized) | **scaffold** — MKL installed, but numpy/scipy MKL-linking + an MKL PARDISO `SparseLU` backend are tracked optimization work (roadmap Part 5) |
| `gpu` | `base-gpu.Dockerfile` | CUDA base for future GPU kernels | **scaffold** — deferred; the swap seam + a future cuDSS `SparseLU` backend |

All bases target the standard CBLAS / LAPACKE / FFTW3 ABIs (and, for `cpu`/`cpu-mkl`,
provision system MUMPS with the pkg-config `.pc` files Debian omits), so the base is
swappable purely via `BASE_IMAGE`.

**`docker/Dockerfile`** — `FROM ${BASE_IMAGE}` (default `qmodeling-base:latest`) and
   layers `build` → `test-deps` → `test` → `runtime` on top. Setup uses the canonical
   `uv sync --all-packages`, which installs `qscat` and builds the Rust `qscat_kernels`
   in one step. `test-deps` adds the `mumps`/`plot` extras only — no test execution;
   `test` is `FROM test-deps` and runs the full `pytest -q` suite, so building/running
   test dependencies never silently pays for (or is blocked by) a full suite run. The
   `runtime` stage starts fresh from the upstream uv/Python image and installs only the
   non-dev shared libs (`libopenblas0`, `libfftw3-double3`, `liblapacke`) needed at
   runtime — no compilers, no `-dev` headers, no `mumps` extra.

## Build

```bash
# docker/build.sh [test|test-deps|runtime] [cpu|cpu-mkl|gpu]   (defaults: test cpu)

# Build the CPU base, then the app test stage (fails the build if tests fail):
docker/build.sh test

# Build the CPU base, then just the test/compute dependencies (mumps + plot),
# with no test execution -- this is what docker/run.sh uses:
docker/build.sh test-deps

# Build the CPU base, then the app runtime image, and run it:
docker/build.sh runtime
docker run --rm qmodeling:runtime-cpu

# Build against a different base (scaffolds):
docker/build.sh test cpu-mkl
docker/build.sh runtime gpu      # requires --gpus all at run time
```

`docker/build.sh` builds the selected base first, then the requested app target, so the
base is never stale relative to the app image. The `cpu` base is also tagged
`qmodeling-base:latest` (the app Dockerfile's default `BASE_IMAGE`).

## Running a `qscat-run` experiment

`docker/run.sh` builds the `test-deps` image (has MUMPS — needed for the larger
production decks, e.g. H2P's ~1.15M-unknown full grid — but does NOT run the test suite)
and runs a `qscat-run` YAML config inside it in one step:

```bash
docker/run.sh CONFIG [OUTPUT_DIR]     # OUTPUT_DIR defaults to runs/<config-stem>

# e.g.
docker/run.sh apps/qscat-run/examples/h2p-dr.yaml runs/h2p-dr
```

It mounts `CONFIG` into the container read-only at `/config.yaml` and `OUTPUT_DIR` out
at `/out`, then runs `uv run --no-sync qscat-run run /config.yaml --output /out` — a
general entry point for any `qscat-run` config (it supersedes the old molecule-specific
`run-n2.sh`). See the top-level `README.md` for the full `qscat-run` walkthrough.

## AWS deploy

AWS deployment pushes `qmodeling-base` to ECR once (it changes only when the
architecture/vendor choice changes) and reuses it across app builds — only the `build`/
`test`/`runtime` layers rebuild on ordinary code changes.
