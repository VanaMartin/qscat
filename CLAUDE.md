# CLAUDE.md — qModeling operating manual

## Purpose

qModeling is a Python-first quantum-mechanics (QM) research monorepo, home to
**QSCAT**. It resuscitates and extends prior C++/CUDA work (`eMoScat`,
`libXcuda`) as a maintainable, CPU-first Python codebase with Rust for the hot
paths that actually need it. GPU/CUDA and AWS deployment are recovered as
reference material but deliberately deferred — not part of this phase.

## The lifecycle

Every method or numerical capability moves through five stages:

1. **Design** — work out the math/algorithm (see `docs/physics/`), decide the
   API shape.
2. **Toy model** — a small, readable Python implementation under `projects/`.
3. **Validate** — check it against analytic benchmarks, conservation laws,
   convergence studies, or a reference implementation (`validation/`).
4. **Optimize in Rust** — once validated and proven a hot path, move it to a
   PyO3/maturin kernel under `native/`, keeping the Python version as the
   differential oracle.
5. **Promote to `qscat`** — validated, reusable code graduates into
   `libs/qscat`, the standard library other projects depend on.

**Invariants that hold at every stage:** the code is CPU-runnable on a laptop
(no GPU required), and it is containerizable (builds and runs under
`docker/`). See the `qm-method-lifecycle` skill for the enforced workflow.

## Repo map

```
libs/       qscat — the standard library: validated, reusable QM code
            (units, linalg, dvr, ecs, special, evolution submodules)
            - qscat.dvr: FEM-DVR-ECS radial grid (`FemDvrEcsGrid`), kinetic-
              energy assembly (`kinetic`), and diagonal-potential Hamiltonian
              + eigensolver helpers (`hamiltonian`, `eigen`) — see
              docs/physics/femdvr-ecs.md.
            - qscat.ecs: exterior-complex-scaling coordinate map (`ecs_map`),
              the single source of the `z(x) = x` / `R0 + (x-R0)e^{i theta}`
              transform used by qscat.dvr's complex tail.
native/     Rust kernels (qscat-kernels crate) built with PyO3/maturin,
            mirroring validated Python APIs for hot paths
projects/   per-problem research and toy models — lifecycle stages 1-2
validation/ analytic benchmarks, golden datasets, convergence studies
reference/  read-only oracles: eMoScat (C++/CUDA snapshot), libXcuda
            (CUDA submodule) — for porting reference only, never imported
docs/       specs/plans (docs/superpowers), physics notes (docs/physics),
            and ADRs (docs/adr)
docker/     layered CPU images: base (architecture/vendor) + app (build/
            test/runtime)
.claude/    skills and agents (below) plus repo permissions settings
```

## Tech decisions

- **Python >=3.12**, managed with `uv`; the exact interpreter is pinned in
  `.python-version` (`3.12`) for reproducibility and Docker parity.
- **uv workspace**: `pyproject.toml` declares `members = ["libs/*", "native/*"]`
  with `package = false` at the root. The canonical setup command is
  `uv sync --all-packages` — a plain `uv sync` prunes the `qscat` and
  `qscat_kernels` workspace members.
- **Rust + PyO3/maturin** for compiled kernels (`native/qscat-kernels`),
  built and installed into the active environment via
  `uv run maturin develop --manifest-path native/qscat-kernels/Cargo.toml`.
- **Testing**: pytest, hypothesis for property-based tests, pytest-benchmark
  for kernel benchmarks; every Rust kernel gets a Python differential test
  against its oracle.
- **Quality tooling**: ruff (lint), mypy --strict (types), cargo clippy (Rust
  lint). Pre-commit (`.pre-commit-config.yaml`) currently enforces only
  **ruff**, **ruff-format**, and **cargo fmt --check**; mypy and clippy are
  run manually / in CI, not (yet) pre-commit hooks.
- **Units**: atomic units throughout (`libs/qscat/qscat/units.py`); no ad hoc
  unit conversions scattered through method code.
- **Docker is layered, not a single Dockerfile**:
  - `docker/base.Dockerfile` builds `qmodeling-base` — the architecture /
    BLAS-FFT vendor layer: OpenBLAS + LAPACK(E) + FFTW3 dev libs, pkg-config,
    a Rust toolchain, and uv/Python 3.12. Code targets the standard CBLAS /
    LAPACKE / FFTW3 ABIs, so this layer is swappable. Default vendor choice
    (OpenBLAS/FFTW3) keeps ARM/Graviton open; an MKL x86-64 variant is a
    planned future alternative, not implemented yet.
  - `docker/Dockerfile` is `FROM ${BASE_IMAGE}` and layers `build` → `test` →
    `runtime` on top, using `uv sync --all-packages` for setup.
  - `docker/build.sh [test|runtime]` builds the base image then the
    requested app target. Verified working: `test` prints `5 passed`;
    `runtime` prints `qscat 0.0.0 ready`.
- **GPU/CUDA and AWS deployment are deferred.** `reference/libXcuda` is kept
  as future-GPU-kernel reference only.

## Skills & agents

| Name | Kind | Use when |
|---|---|---|
| `qm-method-lifecycle` | skill | Adding or porting any QM method/capability — enforces design → toy → validate → optimize-in-Rust → promote-to-qscat. |
| `numerical-validation` | skill | Validating quantum/numerical code where exact equality doesn't apply: analytic benchmarks, convergence studies, conservation checks, differential testing. |
| `python-to-rust-kernel` | skill | A validated Python method has a proven hot path — scaffolding a PyO3/maturin crate under `native/`, mirroring the API, benchmarking, keeping Python as the oracle. |
| `containerize-and-run` | skill | Packaging a capability to run in Docker locally or prep it for AWS — reproducible multi-stage `uv` + `maturin` builds, CPU-only. |
| `qscat-conventions` | skill | Unsure how the project names or measures things — atomic units, FEM-DVR-ECS notation, tolerance defaults, standard-library layout. |
| `port-scout` | agent | Before porting anything from `reference/eMoScat` or `reference/libXcuda` — read-only archaeologist that extracts the math/algorithm, not the C++. |
| `physics-reviewer` | agent | Before promoting a method into `qscat` — reviews for physical/numerical correctness (units, conservation, boundary conditions, ECS handling, convergence), not style. |
| `rust-kernel-engineer` | agent | During the optimize-in-Rust stage — builds PyO3/Rust kernels in `native/` mirroring a validated Python API, with benchmarks and differential tests. |

For general engineering workflow — brainstorming, TDD, systematic debugging,
writing/executing plans, requesting code review — use the **Superpowers**
skills (`superpowers:brainstorming`, `superpowers:test-driven-development`,
`superpowers:systematic-debugging`, `superpowers:writing-plans`,
`superpowers:executing-plans`, `superpowers:requesting-code-review`, etc.).
These are general-purpose and complement the qModeling-specific skills above.

## Reference oracles

`reference/` holds read-only prior codebases: `reference/eMoScat` (C++/CUDA
FEM-DVR-ECS electron-molecule scattering code, snapshot copy) and
`reference/libXcuda` (recovered CUDA layer, git submodule). Never edit,
build, or import these as dependencies — they exist only to be read. Always
use the `port-scout` agent before porting a method out of them; it extracts
the underlying math/algorithm so it can be reimplemented cleanly in Python.

## Common commands

```bash
# Setup — installs qscat and builds the Rust qscat_kernels in one step.
# Plain `uv sync` is NOT sufficient: it prunes the qscat/qscat_kernels members.
uv sync --all-packages

# Run the full test suite
uv run pytest

# Right after a maturin build in the same step, skip the (already-satisfied)
# sync to avoid a redundant resolve:
uv run --no-sync pytest

# Rebuild the Rust kernel extension in place during kernel iteration
uv run maturin develop --manifest-path native/qscat-kernels/Cargo.toml

# Lint / type-check
uv run ruff check .
# mypy is type-clean over the qscat library (type stubs for the Rust
# qscat_kernels extension are pending, so repo-wide strict mypy isn't
# claimed to pass yet):
uv run mypy libs/qscat

# Build and test the CPU Docker images (base, then test target)
docker/build.sh test
```
