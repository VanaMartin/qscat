---
name: containerize-and-run
description: Use when packaging or validating a qModeling capability in the layered CPU Docker images, including choosing compute versus test targets and checking runtime capabilities.
---

# containerize-and-run

## Reuse the CPU layers

Read `docker/README.md`, `docker/Dockerfile`, and `docker/build.sh` before changing
container behavior. `docker/base.Dockerfile` supplies the CPU architecture/vendor
layer: OpenBLAS, LAPACK(E), FFTW3, source-built OpenMP sequential MUMPS, ffmpeg,
Rust, and uv/Python.
The app Dockerfile consumes `BASE_IMAGE`; ordinary Python dependencies belong in
the owning `pyproject.toml`, while system libraries belong in the base layer.

The app stages are:

| Stage | Parent | Capabilities and purpose |
|---|---|---|
| `build` | `BASE_IMAGE` | `uv sync --all-packages --extra plot`; builds the workspace and includes plotting. |
| `test-deps` | `build` | Adds `mumps` and retains `plot`; installs compute/test dependencies without running tests. |
| `test` | `test-deps` | Runs the fast suite in parallel with `-m "not slow" -n auto --dist loadfile`, then the slow suite serially with `-m slow`. |
| `runtime` | Fresh slim uv/Python image | Copies `/app`, including the built environment and editable source, from `build`; adds runtime numerical libraries and ffmpeg. Has plotting but omits MUMPS and the build toolchain. |

Keep `--no-sync` on execution commands after the environment has been built;
reconciliation can prune workspace members or extras. Do not reduce the runtime
copy to `.venv` alone while installed packages still point to editable source.

## Choose the target from the claim

- **Verify a container change:** `docker/build.sh test cpu`. This runs both test
  tiers and fails on a failing tier. Record actual passed/skipped counts and
  capabilities; no fixed test count is an acceptance criterion.
  The current Dockerfile's parallel test `RUN` does not pin BLAS threads; host
  environment variables do not automatically enter that build step. For pinned
  verification with the existing images, build `test-deps` and run both tiers:

  ```bash
  docker/build.sh test-deps cpu
  docker run --rm -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e MKL_NUM_THREADS=1 qmodeling:test-deps-cpu uv run --no-sync pytest -q -m "not slow" -n auto --dist loadfile
  docker run --rm -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e MKL_NUM_THREADS=1 qmodeling:test-deps-cpu uv run --no-sync pytest -q -m slow
  ```

  Record which route ran; adding pins to the Dockerfile is a separate build change.
- **Run a calculation needing MUMPS:** `docker/run.sh CONFIG [OUTPUT_DIR]`. This
  builds the CPU `test-deps` image without running the suite, preserves the
  config's repository-relative path where possible, and mounts the output directory.
  A successful calculation is not evidence that either test tier ran.
- **Check the shipped runtime:**

  ```bash
  docker/build.sh runtime cpu
  docker run --rm qmodeling:runtime-cpu
  ```

  The startup command imports `qscat` and `qscat_kernels` and prints the installed
  `qscat.__version__`. Follow this import smoke check with the changed capability's
  actual execution path; an import alone does not validate a solver or renderer.

`docker/build.sh [test|test-deps|runtime] cpu` builds the base and tags the app as
`qmodeling:<target>-cpu`. `docker/run.sh` uses its own `qmodeling:test-deps` tag.
Use the tag produced by the command you actually ran.

## Reproducibility and resources

- Local workspace setup is `uv sync --all-packages`; plain `uv sync` can prune
  workspace members. Read the committed lockfile and actual build command before
  claiming a frozen dependency install. Current Docker stages use `uv sync`
  without `--frozen`; adding lock enforcement is a separate build change.
- `docker/build.sh` and `docker/run.sh` pass the host commit through `GIT_SHA`.
  Both build and fresh runtime stages set `QSCAT_GIT_SHA`; inspect the produced
  manifest before claiming container provenance is available.
- For a justified Rust kernel, follow `python-to-rust-kernel`. An explicit local
  optimized rebuild uses
  `uv run maturin develop --release --manifest-path native/qscat-kernels/Cargo.toml`.
  Record the actual build profile rather than assuming every workspace build is
  optimized.
- Parallelize only the fast tier with `--dist loadfile`; pin Linux BLAS threads
  to one per worker. Run production tests serially and allocate memory for the
  actual deck. Check competing processes before attributing a kill to grid size.
  See `docs/adr/0005-test-tiers-fast-and-slow.md`.
- Record missing optional backends and skipped tests. Runtime lacks MUMPS;
  production calculations requiring it belong in `test-deps`, not runtime.

GPU/CUDA and AWS deployment are deferred. Their scaffold files do not establish
supported execution. Infrastructure provisioning and public-site deployment are
owned by the separate infrastructure repository, not by these research images.
