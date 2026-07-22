# N-Dimensional Sparse Hamiltonian Library — Design Spec (sub-project #5)

**Date:** 2026-07-22
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending
**Lifecycle:** `qm-method-lifecycle` stages 1–3 for a *library* capability. Builds on #1
(`qscat.dvr`, `qscat.ecs`). Consumed by #6 (the N₂ 2-D exact cross section), which gets
its own spec.

## Context

Sub-projects #1–#4 built a validated 1-D FEM-DVR-ECS stack and used it to solve the
electron–N₂ vibrational-excitation problem twice (time-independent and time-dependent)
under the **Local Complex Potential** approximation. Both agree with each other and land
within a documented factor-3 cross-model bound of Karel Houfek's golden data, with two
honestly-documented structural limitations of the LCP model itself.

The next step is to drop the approximation and solve the **2-D model exactly** — the same
model Houfek solved. Port-scout archaeology
(`.superpowers/sdd/n2-2d-exact-extraction.md`) established that eMoScat's 2-D Hamiltonian
is assembled as a **Kronecker sum of 1-D FEM-DVR-ECS kinetic operators plus a diagonal
potential** (`source/FemDvrEcs2d/OperatorRowCompressed2d.cpp:89`), and solved by sparse
LU.

That construction is **not specific to 2-D, and not specific to molecular scattering**. It
is the general recipe for a separable-kinetic, diagonal-potential Hamiltonian on a tensor
product of DVR grids in any number of dimensions. This sub-project builds it as reusable
library code, validated on its own terms, *before* any N₂ physics depends on it.

**Why the split matters:** in sub-projects #1–#3, having the grid and kinetic assembly
independently validated against analytic benchmarks meant that when cross sections
disagreed with reference data, we could attribute the disagreement to physics rather than
re-litigating the numerics. Preserving that property is the main design driver here.

## Scope

Two layers, with a strict dependency direction (`dvr` → `linalg`, never the reverse).

### Layer 0 — `qscat.linalg` (currently an empty stub)

Pure sparse linear algebra. No grids, no physics, no ECS.

```python
def kron_sum(ops: Sequence[sparse.spmatrix]) -> sparse.csr_matrix:
    """Σ_d  I ⊗ … ⊗ ops[d] ⊗ … ⊗ I   for arbitrary len(ops).

    Each ops[d] is square; the result is square of size prod(n_d). Acts on
    C-order-raveled arrays of shape (n_0, …, n_{D-1}) — i.e. LAST axis fastest.
    D == 1 returns ops[0] unchanged.
    """

class SparseLU:
    """Cached sparse LU factorization: factor once, solve many right-hand sides.

    Wraps scipy.sparse.linalg.splu. Exposes fill-in and memory diagnostics
    because at production problem sizes those, not flops, decide feasibility.
    """
    def __init__(self, A: sparse.spmatrix, *, ordering: str = "COLAMD") -> None: ...
    def solve(self, b: NDArray) -> NDArray:      # (N,) or (N, k)
        ...
    @property
    def fill_factor(self) -> float: ...          # (L.nnz + U.nnz) / A.nnz
    @property
    def memory_bytes(self) -> int: ...

def c_product(a: NDArray, b: NDArray) -> complex:
    """The ECS inner product: Σ a_i b_i, WITHOUT conjugation."""
```

`kron_sum` taking arbitrary sparse matrices — rather than grids — is deliberate: a future
angular DVR, finite-difference, or B-spline dimension composes with FEM-DVR-ECS
dimensions at no extra cost.

`c_product` exists because the conjugate-vs-no-conjugate distinction under ECS has already
been a recurring correctness trap (sub-project #3's S-matrix; the reference code's use of
`cblas_zdotc`, which is formally wrong and survives only because it zeroes every channel
function on the scaled tail). Naming it makes the choice explicit at every call site.

### Layer 1 — `qscat.dvr`, dimension-general operator assembly

```python
def kinetic_sparse(grid: FemDvrEcsGrid, mass: float) -> sparse.csr_matrix:
    """Sparse sibling of the existing dense `kinetic()`, assembled by COO scatter."""

class TensorGrid:
    """Tensor product of D FEM-DVR-ECS grids. C-order: last axis fastest."""
    def __init__(self, grids: Sequence[FemDvrEcsGrid]) -> None: ...
    @property
    def grids(self) -> tuple[FemDvrEcsGrid, ...]: ...
    @property
    def ndim(self) -> int: ...
    @property
    def shape(self) -> tuple[int, ...]: ...      # (n_0, …, n_{D-1})
    @property
    def size(self) -> int: ...                   # prod(shape)
    def points(self) -> tuple[NDArray, ...]:     # D broadcastable complex arrays
        ...
    def real_mask(self) -> NDArray[np.bool_]:    # flat, True where ALL coords unscaled
        ...
    def outer(self, vectors: Sequence[NDArray]) -> NDArray:
        """Separable state ⊗_d vectors[d] → flat vector of length `size`."""

def kinetic_nd(tgrid: TensorGrid, masses: Sequence[float]) -> sparse.csr_matrix:
    """= kron_sum([kinetic_sparse(g, m) for g, m in zip(tgrid.grids, masses)])"""

def potential_nd(tgrid: TensorGrid, V: Callable[..., ArrayLike]) -> NDArray:
    """V evaluated at the D-dimensional COMPLEX (ECS-scaled) points, flattened."""

def hamiltonian_nd(tgrid, masses, V) -> sparse.csr_matrix:
    """kinetic_nd(tgrid, masses) + diags(potential_nd(tgrid, V))"""
```

## Design decisions

**1. Index ordering: numpy-native C-order (last axis fastest).** `psi[i_0, …, i_{D-1}].ravel()`
pairs with `kron_sum`. eMoScat uses the opposite convention (first coordinate fastest,
`idx = i_r + i_R·N_r`, `FemDvrEcsGrid2d.cpp:169`). The two are physically identical and
differ only in basis ordering; we take the numpy-native one so that `reshape(tgrid.shape)`
does the obvious thing. **Documented in the `TensorGrid` docstring** so nobody diffs raw
index dumps against the reference and concludes something is broken.

**2. The existing dense `kinetic()` stays, as the differential oracle for
`kinetic_sparse()`.** Same pattern the repo uses for Rust kernels: the validated slow
implementation is retained specifically to test the fast one against. Their difference must
be at round-off (`< 1e-12` relative).

**3. `real_mask` is a safety feature, not a convenience.** Under ECS, any driving term or
channel projection must be confined to the unscaled region or the resulting matrix element
is meaningless. Making the mask a first-class property of `TensorGrid` means the physics
layer cannot silently omit it.

**4. The potential is diagonal.** This is the DVR approximation and is what makes the whole
Kronecker-sum construction valid; it is already the convention in `qscat.dvr.operators`.
Non-diagonal potentials are explicitly out of scope.

**5. Complex symmetric, not Hermitian.** ECS makes `H = Hᵀ ≠ H†`. Every routine must use
general (not Hermitian) algorithms. Verified empirically at production size: `max|H − Hᵀ| =
1.7e-13`.

## Validation

The point of this sub-project is that the library is proven correct *independently of any
physics application*.

**V1 — `kron_sum` against dense `np.kron`.** Random small sparse matrices, D = 1, 2, 3, 4,
with unequal dimensions per axis (catches transposed-index bugs that square cases hide).
Exact to round-off.

**V2 — `kinetic_sparse` against dense `kinetic`.** Differential test over several grid
specs including ECS tails. Relative difference < 1e-12, and the sparsity pattern must match
the analytic prediction `nnz = q²·tnel − 4q + 3 − tnel` (eMoScat `KineticEnergy.cpp:95`).

**V3 — analytic separable benchmarks at D = 1, 2, 3.** The generality is *exercised*, not
asserted:
- **D-dimensional particle in a box**: eigenvalues are sums of 1-D eigenvalues
  `Σ_d n_d²π²/(2 m_d L_d²)`. Known in closed form.
- **D-dimensional harmonic oscillator**: `Σ_d ω_d(n_d + ½)`.
Both with unequal per-axis extents/masses. Lowest eigenvalues to ~1e-10.

**V4 — D = 1 reproduces the existing 1-D stack bit-for-bit.** `TensorGrid([g])` +
`hamiltonian_nd` must equal `qscat.dvr.hamiltonian(g, V, m)` to round-off. This makes every
result already validated in sub-projects #1–#4 a regression test on the new code.

**V5 — `SparseLU` correctness and reuse.** Residual `‖Ax − b‖/‖b‖ < 1e-12` on a
complex-symmetric ECS matrix; multi-RHS solve agrees with looped single solves; one
factorization reused across many right-hand sides gives identical results.

**V6 — production-scale smoke test.** Assemble the real N₂ 2-D Hamiltonian
(N = 143,380) and check dimension, nnz, and complex symmetry — **without** factorizing
(too heavy for CI). Marked `slow`.

## Performance evidence (measured, not estimated)

Spike on the production `N2-model.json` deck, this laptop:

| quantity | value |
|---|---|
| N | 335 × 428 = **143,380** |
| nnz | **3,276,450** (22.9/row) — matches eMoScat's formula exactly |
| assembly | 0.07 s, 0.05 GB |
| `max\|H − Hᵀ\|` | 1.7e-13 |
| `splu` (COLAMD) | **128 s**, fill-in ×93, L+U = 3.05e8 nnz ≈ 4.87 GB |
| peak RSS | **13.6 GB** |
| back-substitution | 440 ms, residual 1.25e-15 |

Consequences that shape sub-project #6:
- All final channels at one energy **share one factorization**, so the 6 benchmark anchors
  are only 3 distinct energies ≈ 6.4 min.
- 13.6 GB peak RSS is too much for the Docker harness; a reduced grid with a documented
  tolerance will be needed there.
- A full σ(E) curve at ~128 s/energy is hours — a separate undertaking, not part of #6.

**Open tuning question, deferred to implementation:** whether `permc_spec="MMD_AT_PLUS_A"`
(appropriate for our structurally symmetric pattern) beats the default COLAMD's ×93
fill-in, and how much a shorter electronic box helps. `SparseLU` exposes `ordering` and the
fill/memory diagnostics precisely so this can be measured rather than guessed. This is a
performance choice; it cannot affect correctness.

## Out of scope

- **All N₂ physics** — driving term, T-matrix, cross sections. That is sub-project #6.
- Non-diagonal potentials; coupled partial waves.
- Iterative solvers, preconditioning, matrix-free operators. Direct sparse LU is
  demonstrably sufficient at the sizes we need.
- Sparse *time propagation*. `qscat.evolution.make_cn_stepper` is dense-only; generalizing
  it to sparse is needed for a future TD 2-D route, not for #6.
- Rust. Assembly is 0.07 s; the cost is entirely inside SuperLU. There is no hot path here
  to optimize.
- GPU/CUDA (deferred repo-wide). The port-scout confirmed eMoScat's 2-D path was never
  GPU either.

## Verification

- `uv run pytest libs/qscat -q` → all pass, including V1–V5 (V6 marked slow).
- `uv run mypy libs/qscat` → 0 errors.
- `uv run ruff check .` → clean.
- The N₂ harness (`uv run python -m validation.n2.experiment`) still reports
  **19 PASS / 0 PENDING / 2 NOTE / 0 FAIL**, exit 0 — this sub-project must not disturb
  any existing physics result.
- `docker/build.sh test` → passes.
- `CLAUDE.md` updated with `qscat.linalg` and the `qscat.dvr` tensor additions.

Note for anything runnable: `uv run python` does not inherit the repo-root `pythonpath`
that pytest gets from `pyproject.toml`, so scripts importing `projects.*`/`validation.*`
need `PYTHONPATH=.`. (Conda `base` may be active in the shell; `uv run` correctly ignores
it and uses `.venv/bin/python3` at 3.12.7. Bare `python` would not.)
