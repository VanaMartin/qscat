# MUMPS Complex-Symmetric Sparse Solver Backend — Design Spec (sub-project #8)

**Date:** 2026-07-26
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending
**Lifecycle:** `qm-method-lifecycle` **stage 4 (optimize the hot path)** for the sparse
factorization. The first optimization sub-project. Builds on #5 (`qscat.linalg.SparseLU`),
consumed transparently by #6 (the driven solve) and #7 (the CN propagation).

## Context

The **sparse LU factorization is the proven hot path** for the exact 2-D solvers. On the full
N₂ 2-D grid (N = 143,380) the numeric factorization is ~128 s and ~13.6 GB with ×93 fill-in;
the back-substitution is already cheap (~0.44 s). `qscat.linalg.SparseLU` currently wraps
`scipy.sparse.linalg.splu` (SuperLU), a **general unsymmetric** solver.

But every matrix we factor is **complex-symmetric** (`A = Aᵀ ≠ A†`, verified to ~1e-13
throughout): the driven matrix `E_tot·I − H` (#6) and the Crank-Nicolson matrix
`A = I + iH·dt/2` (#7) are both symmetric because `H = Hᵀ`. SuperLU cannot exploit this — it
factors both triangles when one suffices. A solver with a **complex-symmetric mode** does
roughly half the work and uses roughly half the memory. This is exactly why eMoScat's original
implementation used a symmetric sparse direct solver (Intel MKL PARDISO) and reportedly
factored **all of N₂/NO/F₂ in under an hour** — the performance bar for this line of work
(see the `sparse-lu-optimization-target` memory).

**Chosen backend: MUMPS** (open-source, CeCILL-C), for its complex-symmetric mode
(`SYM=2`, `zmumps`), nested-dissection ordering, and — decisively — because it runs on **both
x86 and ARM** (Graviton, Apple Silicon), so one backend covers the dev box and every AWS
instance type. (MKL PARDISO is faster on x86 but x86-only; a PARDISO variant behind the same
dispatch is a possible later addition.)

**Scope (deliberately tight):** the benchmark harness + the MUMPS complex-symmetric backend
behind `SparseLU`'s dispatch, differential-tested against SuperLU. Nothing else.

## Verified environment findings (not assumptions)

- `python-mumps` (the maintained Kwant binding) is **not** a bundled wheel — it builds against
  a **system** MUMPS library (`dmumps_seq`, found via pkg-config). So MUMPS is a **vendor
  library** that must be provisioned first, exactly like OpenBLAS/LAPACKE/FFTW3 already are in
  `docker/base.Dockerfile`.
- **Docker/Debian is the clean path:** `libmumps-seq-dev` is in apt main. This is the primary
  validated environment (and where CI/AWS runs live).
- **macOS (Apple Silicon)** MUMPS is not in homebrew-core (needs a science tap or conda-forge)
  — more involved. Since the heavy runs belong in Docker/AWS and the Mac can run the
  container, Mac-native MUMPS is a documented convenience, not a blocker.
- MUMPS's fill-reducing quality depends on its ordering library (METIS/SCOTCH). The benchmark
  must **report which ordering MUMPS actually used**, since the fill-in win rides on it.

## Design — four layers, strict boundaries

1. **Vendor layer.** Add `libmumps-seq-dev` (+ its ordering lib if separate) to
   `docker/base.Dockerfile` beside the BLAS/FFTW libs, and extend the existing pkg-config
   sanity gate to check for it. This is the "base image decides architecture" layer.
2. **Binding.** `python-mumps` as an **optional** dependency (`qscat[mumps]` extra), never
   core. **Core `qscat` stays numpy/scipy-only** and universally installable via the SuperLU
   path.
3. **Dispatch inside `SparseLU`.** A `backend` selector — `"auto"` (MUMPS if the binding +
   library import, else SuperLU), `"scipy"` (force SuperLU), `"mumps"` (force MUMPS, error if
   absent). The MUMPS path auto-detects `A = Aᵀ` (cheap, O(nnz)) and uses **`SYM=2`
   complex-symmetric** factorization; a non-symmetric `A` falls to `SYM=0`. **SuperLU remains
   the always-available fallback AND the differential oracle.**
   **Correctness trap (a V1 differential-test target): MUMPS `SYM=2` takes only the UPPER
   TRIANGLE of `A`** (its convention) — supplying the full matrix double-counts off-diagonals
   and silently corrupts the factorization. The dispatch must honor the binding's
   symmetric-input contract (upper-triangular entries only for the symmetric path).
4. **Benchmark + differential test.** A `pytest-benchmark` harness measuring factor + solve
   (time and peak memory) on the real N₂ 2-D matrices, MUMPS-symmetric vs SuperLU, reported
   against the "<1 hr" bar; and a differential test asserting the two backends' solves agree
   to round-off on the same complex-symmetric matrix.

## Interface (backward-compatible)

```python
class SparseLU:
    def __init__(
        self, A, *,
        ordering: _Ordering = "COLAMD",   # SuperLU-path only (MUMPS uses its own analysis)
        backend: Literal["auto", "scipy", "mumps"] = "auto",
        symmetric: bool | None = None,    # None => auto-detect A == A.T; overridable
    ) -> None: ...
    def solve(self, b): ...              # unchanged signature
    @property
    def backend_used(self) -> str: ...    # "scipy" | "mumps" -- which actually ran
    @property
    def fill_factor(self) -> float: ...   # surfaces MUMPS's stat when MUMPS ran
    def memory_bytes(self) -> int: ...    # surfaces MUMPS's stat when MUMPS ran
    @property
    def ordering_used(self) -> str: ...   # SuperLU permc_spec, or MUMPS's chosen ordering
```

Every current call — `SparseLU(A)`, `SparseLU(A, ordering=...)` — keeps working: `backend`
defaults to `"auto"`, so with MUMPS absent the behaviour is bit-identical to today, and with
MUMPS present the symmetric factorization engages transparently. `make_sparse_cn_stepper`
(#7) and the driven solve (#6) both go through `SparseLU`, so both accelerate with **no change
to their code**.

## Validation

**V1 — differential correctness (the gate).** On several complex-symmetric matrices
(including one at real N₂ 2-D scale), the MUMPS `SYM=2` solve agrees with the SuperLU solve to
round-off (`‖x_mumps − x_scipy‖/‖x‖ < 1e-10`) and both give `‖Ax − b‖/‖b‖ < 1e-10`. "Faster"
must never mean "wrong."

**V2 — the N₂ physics results are unchanged.** With `backend="mumps"`, #6's exact anchors
still match Houfek (the group-E numbers) and #7's σ_TD still matches σ_TI — i.e. swapping the
factorization backend changes no physics. Run the existing anchor checks through both backends.

**V3 — the measured speedup + memory (the deliverable).** The benchmark reports MUMPS-symmetric
vs SuperLU factor time, solve time, peak memory, fill-in, and the ordering MUMPS used, on the
real N₂ 2-D matrix. **The benchmark measures the win; it does not assume it.** The expectation
from exploiting symmetry + nested dissection is ~2× faster / ~half memory, but if MUMPS-seq
(without a strong ordering lib) does not win, that is a real, reportable finding — and the
SuperLU fallback means nothing regresses. Report the numbers against the "<1 hr all-models"
bar.

**V4 — fallback + absence.** With MUMPS not installed, `backend="auto"` silently uses SuperLU
and every result is bit-identical to today; `backend="mumps"` raises a clear error naming the
missing binding/library. The core (numpy/scipy-only) test suite passes with no MUMPS present.

## Out of scope (each a later sub-project)

- **`complex64` factorization + iterative refinement** (Lever D).
- **Symbolic/numeric reuse across the TI energy sweep** (Lever C — one analysis, N numeric
  factorizations). The dispatch is designed not to preclude it, but it is not built here.
- **A Rust/PyO3 kernel** for the non-LU hot paths (assembly, the propagation loop). The LU
  itself delegates to MUMPS — reimplementing a competitive sparse direct solver is not the win.
- **MKL PARDISO** as a second dispatch backend (x86 speed).
- **The publishing pipeline** — multi-arch Docker MUMPS-variant images, PyPI core + `[mumps]`
  extra, conda-forge, maturin/cibuildwheel platform wheels. A dedicated "packaging" sub-project.
- Iterative solvers / preconditioning (Lever E).

## Verification

- `uv run pytest libs/qscat -q` (core, no MUMPS) → all pass, results bit-identical to `main`.
- In the Docker `test` image (MUMPS present): V1 differential + V2 physics-unchanged + V3
  benchmark all pass/run; V4 fallback verified.
- `uv run mypy libs/qscat` → 0; `uv run ruff check .` → clean.
- The N₂ harness `python -m validation.n2.experiment` → **23 PASS / 0 PENDING / 6 NOTE /
  0 FAIL**, exit 0, unchanged (with and without MUMPS).
- `docker/build.sh test` → passes with MUMPS provisioned in the base.
- The benchmark numbers (factor/solve/memory/fill/ordering, MUMPS vs SuperLU, vs the bar)
  recorded in `docs/physics/` (or a benchmark note) and `CLAUDE.md` updated (the MUMPS backend,
  the `qscat[mumps]` extra, the base-image vendor addition).
