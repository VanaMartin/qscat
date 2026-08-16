# Sparse shift-invert eigensolver — interior complex eigenvalues of an ECS Hamiltonian

Date: 2026-08-16

## Purpose

Add a **sparse shift-invert eigensolver** to `qscat.linalg`: given a large sparse
complex-symmetric matrix `A` and a complex shift `sigma`, return the `k`
eigenvalues nearest `sigma` and their eigenvectors, reusing the factorization
machinery `qscat.linalg.SparseLU` already provides.

qscat today can only diagonalize densely (`qscat.dvr.eigen` → `np.linalg.eig`),
which caps every eigenvalue computation at 1-D grid sizes. The electronic pole
walk behind `qscat.core.lcp` works precisely because a fixed-`R` electronic
Hamiltonian is small. Nothing in the library can ask for eigenvalues of the 2-D
tensor Hamiltonian `hamiltonian_nd` assembles (143k unknowns for the N₂ working
deck, ~1.15M for H₂⁺ at full size).

This spec covers **only the primitive and its 1-D validation**. It is stage 1 of
a longer program whose eventual target is the *exact* (non-Born-Oppenheimer)
resonance states of the 2-D model — the true S-matrix poles, against which the
existing BO/LCP levels (`qscat.core.lcp.resonance_levels`) can be measured. That
target motivates the interface but is deliberately out of scope here: the
primitive is validated where an exact oracle already exists — in 1-D, against
`np.linalg.eig` — before it is trusted at a scale where no oracle does.

## Why shift-invert, and why it is cheap here

Resonances are **interior** eigenvalues: they sit in the lower half of the
complex energy plane, surrounded by the discretized rotated continuum, not at
either end of the spectrum. Krylov methods converge to extremal eigenvalues, so
the standard remedy is the shift-invert spectral transform — run Arnoldi on
`(A - sigma*I)^-1`, whose extremal eigenvalues are `A`'s eigenvalues nearest
`sigma`.

That transform needs one sparse solve per matrix-vector product, which is
exactly what `SparseLU` does, and the shifted matrix `A - sigma*I` has the **same
sparsity pattern for every `sigma`** (a diagonal shift). `SparseLU.refactor`
therefore reuses the symbolic analysis across shifts — the same reuse that makes
the time-independent energy sweep ~5× cheaper on the N₂ deck
(`docs/physics/ti-energy-sweep-reuse.md`), and on the MUMPS backend the one that
skips the SCOTCH ordering (`docs/physics/mumps-sparse-backend.md`). A resonance
hunt is a sweep of shifts through the complex plane; it gets the same discount.

## What already exists

- `qscat.linalg.SparseLU` — cached sparse LU, SuperLU/MUMPS backends,
  `refactor(A_new)` for a same-pattern re-factorization, and fill/memory
  diagnostics.
- `qscat.linalg.c_product` — the bilinear (non-conjugated) ECS inner product, the
  correct normalization for complex-scaled eigenvectors.
- `qscat.dvr.eigen(H)` — dense `np.linalg.eig`, sorted by ascending `Re E`,
  returning Euclidean-normalized (`v†v = 1`) eigenvectors with a docstring note
  that ECS observables require re-normalizing under the c-product. **This is the
  oracle for stage 1.**
- `qscat.ecs.find_resonance_pole` / `match_angle_stable` — two-spectrum
  angle-stability matchers that consume already-computed spectra. They are
  unchanged by this work and are what a later stage will feed.
- `qscat.exceptions` — `QscatError` base with `ConvergenceError`, `BackendError`,
  `GridError`, `ModelError`.

## Design

### Interface

New module `libs/qscat/qscat/linalg/eigs.py`, exported from `qscat.linalg`.

```python
class ShiftInvertEigs:
    def __init__(
        self,
        A: sp.spmatrix,
        *,
        k: int = 6,
        ordering: _Ordering = "COLAMD",
        backend: _Backend = "auto",
        symmetric: bool | None = None,
        ncv: int | None = None,
        tol: float = 0.0,
        maxiter: int | None = None,
    ) -> None: ...

    def near(
        self, sigma: complex, *, k: int | None = None
    ) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]: ...

    # diagnostics, delegated to the held SparseLU
    @property
    def shape(self) -> tuple[int, int]: ...
    @property
    def backend_used(self) -> str: ...
    @property
    def ordering_used(self) -> str: ...
    @property
    def fill_factor(self) -> float: ...
    def memory_bytes(self) -> int: ...
    @property
    def n_factorizations(self) -> int: ...
```

A **class**, not a function, because the object owns the factorization: the first
`near` constructs the `SparseLU`, every later `near` calls `refactor`. The
constructor's `ordering` / `backend` / `symmetric` mirror `SparseLU`'s verbatim
and are forwarded unchanged; `k` / `ncv` / `tol` / `maxiter` are ARPACK controls,
with `k` overridable per call.

### Mechanics of `near(sigma)`

1. Form `M = A - sigma * I` — **scipy's sign convention**, not the driven
   solver's `E*I - H`. Both are legal shift-inverts, but ARPACK is told the shift
   separately, so passing `sigma*I - A` as `OPinv` returns eigenvalues reflected
   about `sigma`. This is the single easiest way to get a plausible-looking wrong
   answer here, so it gets its own regression test.
   Adding the explicit identity also forces the diagonal to be structurally
   present, which keeps the pattern constant across shifts even if `A` has
   structural zeros on its diagonal.
2. First call: `self._lu = SparseLU(M, ordering=..., backend=..., symmetric=...)`.
   Later calls: `self._lu.refactor(M)`. `n_factorizations` counts both and is how
   the reuse is asserted in tests.
3. `OPinv = LinearOperator(shape, matvec=self._lu.solve, dtype=complex128)`.
4. `vals, vecs = scipy.sparse.linalg.eigs(A, k=k, sigma=sigma, OPinv=OPinv,
   ncv=ncv, tol=tol, maxiter=maxiter)`. SciPy applies the spectral transform and
   returns eigenvalues of `A`, not of the transformed operator.
5. Sort by `|E - sigma|` ascending and return.

### Conventions, and why they differ from `dvr.eigen`

- **Ordering: nearest-shift-first**, not ascending `Re E`. A shift-invert result
  is a local window around `sigma`; its meaningful order is distance from the
  shift, and callers want the front of the array. `dvr.eigen` returns a *whole*
  spectrum, where ascending `Re E` is the meaningful order. Both are documented
  at their call sites; neither should be changed to match the other.
- **Eigenvectors: Euclidean-normalized**, exactly as `dvr.eigen` returns them,
  with the same docstring instruction to re-normalize under `c_product` for ECS
  observables. Consistency with the existing function outweighs the convenience
  of normalizing here, and the correct normalization depends on the region the
  caller integrates over (the LCP code normalizes over the *real* region only).

### Errors

- `k >= n` → `ValueError` naming `qscat.dvr.eigen` as the dense route for small
  problems (ARPACK requires `k < n`).
- ARPACK non-convergence → `ConvergenceError` (from `qscat.exceptions`) wrapping
  `scipy.sparse.linalg.ArpackNoConvergence`, carrying whatever partial
  eigenvalues ARPACK produced and naming the three remedies: raise `k`, raise
  `ncv`, or move `sigma`.
- A non-square `A`, or a `refactor` pattern mismatch, is left to `SparseLU`'s
  existing raises — this class adds no second copy of those checks.

## Validation

`libs/qscat/tests/test_shift_invert_eigs.py`. Every case is fast (seconds);
none is `@slow`.

1. **Synthetic differential.** A random sparse complex-symmetric matrix
   (`n ≈ 200`, fixed seed) against `np.linalg.eig`: the `k` eigenvalues nearest
   `sigma` must agree to `rtol = 1e-9`. That is the repo's cross-architecture
   floor for sparse-solve comparisons — tighter tolerances have failed CI on a
   different BLAS, and this comparison runs through a sparse factorization.
   Eigenvectors are compared after normalizing each to unit **c-norm**
   (`sqrt(vᵀv)`, not the Euclidean norm), requiring `|vᵀw| ≈ 1`; that is the
   right notion of "same eigenvector up to scale" for a complex-symmetric
   operator, and it is well-defined as long as `vᵀv ≠ 0`, which holds away from
   an exceptional point.
2. **Physics differential.** The real 1-D electronic FEM-DVR-ECS Hamiltonian
   built from `qscat.dvr` + `qscat.model.N2` at a fixed `R`: the resonance pole
   found through the sparse path must equal the one `dvr.eigen` finds, to the
   same `rtol = 1e-9`. This is the differential oracle that makes the primitive
   trustworthy; the test lives with the linalg tests but imports `qscat.model`,
   which is allowed for tests (the *library* layering rule constrains
   `qscat.core`, not the test suite).
3. **Sign-convention pin.** A spectrum deliberately asymmetric about `sigma`, so
   that `A - sigma*I` and `sigma*I - A` cannot agree; asserts the returned
   eigenvalues are the true nearest ones.
4. **Refactor reuse.** Two `near()` calls on one object return exactly what two
   freshly-constructed objects return, and `n_factorizations == 2` with a single
   symbolic analysis — reuse must be invisible in results and visible in
   diagnostics.
5. **Robustness probe.** The behaviour that will bite at 2-D scale, asserted
   qualitatively here and tabulated in the physics note:
   - **shift distance** — sweep `|sigma - E_pole|` and record the largest
     distance at which the pole still appears among the `k` returned;
   - **`k` / `ncv`** — the smallest values that still capture the pole at a
     realistic shift;
   - **continuum-adjacent shift** — place `sigma` on the rotated continuum and
     confirm what returns is continuum (it moves when the ECS angle changes,
     whereas the pole does not). This is the 1-D rehearsal of the
     angle-stability selection a later stage will need.

Gates: `uv run ruff check .`, `uv run mypy libs/qscat/qscat` clean, numpydoc-style
docstrings so the strict (`-W`) Sphinx build stays warning-free.

## Documentation

- New physics note `docs/physics/shift-invert-eigensolver.md`, sitting beside
  `mumps-sparse-backend.md` and `ti-energy-sweep-reuse.md` (the precedent for
  linear-algebra notes in that directory): the spectral transform, the
  same-pattern reuse argument, the convention differences against `dvr.eigen`,
  and the robustness table from validation case 5.
- `CHANGELOG.md` under Unreleased → Added.

## Out of scope

Each is a later stage, not a deferred detail of this one:

- **Anything 2-D** — no use of this primitive on `hamiltonian_nd`, no 2-D
  resonance search, no `qscat.core` API.
- **Two-angle selection at 2-D** — `find_resonance_pole` / `match_angle_stable`
  are unchanged; generalizing them to two ECS angles (`theta_r`, `theta_R`) is
  the next stage's problem.
- **Changes to `qscat.core.lcp`** — the electronic pole walk keeps its dense
  path, which remains the oracle. Swapping it onto the sparse primitive is a
  separate, optional change with its own justification (it would need a measured
  speedup to be worth touching validated physics code).
- **Complex-symmetric Lanczos.** ARPACK's general non-Hermitian iteration ignores
  the c-product symmetry of an ECS Hamiltonian; a complex-symmetric Lanczos
  (three-term recurrence) could be cheaper. That is an optimization to consider
  only if ARPACK proves to be the 2-D bottleneck, with ARPACK as its oracle.
- **MUMPS parameter tuning** for near-singular shifted matrices.

## Known risk

Shift-invert works *because* `A - sigma*I` is near-singular — that ill
conditioning is the amplification mechanism, and ARPACK is designed for it.
But MUMPS `SYM=2` may warn or lose accuracy when `sigma` sits pathologically
close to an eigenvalue. The shift-distance probe (validation case 5) is where
this surfaces. If it does, the response is a documented floor on
`|sigma - E|` relative to the spectrum scale — recorded as a limitation, not
hidden behind an automatic nudge.

## Success criteria

- `ShiftInvertEigs` reproduces `np.linalg.eig` on both the synthetic and the real
  1-D electronic Hamiltonian, eigenvalues and eigenvectors, within stated
  tolerances.
- Reuse across shifts is asserted, not assumed.
- The robustness behaviour (shift distance, `k`/`ncv`, continuum-adjacent shift)
  is measured and written down, so stage 2 starts from numbers rather than
  guesses.
