# Sparse shift-invert eigensolver Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `qscat.linalg.ShiftInvertEigs` — the `k` eigenvalues/eigenvectors of a large sparse complex-symmetric matrix nearest a complex shift — and validate it in 1-D against the dense `np.linalg.eig` oracle.

**Architecture:** A class holding a `SparseLU` factorization of `A - sigma*I`. `near(sigma)` re-forms the shifted matrix, factors it (first call) or `refactor`s it (later calls, reusing the symbolic analysis since a diagonal shift preserves the sparsity pattern), hands the solve to ARPACK via a `LinearOperator` as `OPinv`, and returns eigenvalues of `A` sorted by distance from `sigma`.

**Tech Stack:** numpy, scipy (`scipy.sparse`, `scipy.sparse.linalg.eigs` = ARPACK znaupd), `qscat.linalg.SparseLU` (SuperLU/MUMPS), `qscat.exceptions`, pytest.

**Spec:** `docs/superpowers/specs/2026-08-16-shift-invert-eigensolver-design.md`

## Global Constraints

- Python >= 3.12; atomic units throughout; no new runtime dependencies.
- `qscat.linalg` knows nothing about grids, potentials, or physics — the class takes a matrix, never a model or a grid.
- Differential tolerance against dense: `rtol = 1e-9` (the repo's cross-architecture floor for comparisons that run through a sparse factorization; tighter has failed CI on a different BLAS).
- Eigenvectors are returned Euclidean-normalized, exactly like `qscat.dvr.eigen`; callers re-normalize under `c_product` for ECS observables.
- Eigenvalues are returned sorted by `|E - sigma|` ascending — NOT by `Re E` (which is `dvr.eigen`'s convention for a full spectrum).
- Gates before each commit: `uv run ruff check .`, `uv run ruff format` on touched files only (the repo is not format-clean at HEAD), `uv run mypy libs/qscat/qscat`, and the touched tests.
- Run tests in the foreground (`uv run --no-sync pytest ...`); backgrounded pytest reports exit 0 with empty output in this environment.

---

### Task 1: `ShiftInvertEigs` core — construction, `near()`, sorting, guards

**Files:**
- Create: `libs/qscat/qscat/linalg/eigs.py`
- Modify: `libs/qscat/qscat/linalg/__init__.py` (import + `__all__` + the "Public API" docstring list)
- Test: `libs/qscat/tests/test_shift_invert_eigs.py`

**Interfaces:**
- Consumes: `SparseLU(A, *, ordering, backend, symmetric)`, `.refactor(A_new)`, `.solve(b)` from `qscat.linalg.sparse_lu`; `_Backend` / `_Ordering` type aliases from the same module; `ConvergenceError` from `qscat.exceptions`.
- Produces: `ShiftInvertEigs(A, *, k=6, ordering="COLAMD", backend="auto", symmetric=None, ncv=None, tol=0.0, maxiter=None)` with `near(sigma: complex, *, k: int | None = None) -> tuple[NDArray[complex128], NDArray[complex128]]`. Tasks 2-4 rely on these exact names.

- [ ] **Step 1: Write the failing tests**

Create `libs/qscat/tests/test_shift_invert_eigs.py`:

```python
"""Tests for `qscat.linalg.ShiftInvertEigs`.

The oracle throughout is dense `np.linalg.eig` on the same matrix: this class
must return exactly the eigenpairs the dense solver finds nearest the shift.
Matrices are COMPLEX SYMMETRIC (A == A.T, not Hermitian) -- what exterior
complex scaling produces and what every real use of this class will be.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from qscat.linalg import ShiftInvertEigs, c_product


def _complex_symmetric(n: int, seed: int) -> sp.csc_matrix:
    """A well-conditioned, sparse, complex-symmetric test matrix."""
    rng = np.random.default_rng(seed)
    nnz = 5 * n
    rows = rng.integers(0, n, size=nnz)
    cols = rng.integers(0, n, size=nnz)
    vals = rng.standard_normal(nnz) + 1j * rng.standard_normal(nnz)
    m = sp.coo_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex).tocsr()
    m = m + m.T                                              # symmetric, NOT conjugated
    m = m + sp.identity(n, format="csr", dtype=complex) * (10.0 + 3.0j)
    return sp.csc_matrix(m)


def _nearest_dense(A: sp.csc_matrix, sigma: complex, k: int) -> np.ndarray:
    """The k dense eigenvalues nearest `sigma`, ascending in |E - sigma|."""
    w = np.linalg.eigvals(A.toarray())
    return w[np.argsort(np.abs(w - sigma))][:k]


def test_eigenvalues_match_dense_nearest_the_shift() -> None:
    A = _complex_symmetric(200, seed=0)
    sigma = 9.0 + 2.0j
    k = 6
    vals, _ = ShiftInvertEigs(A, k=k).near(sigma)
    assert np.allclose(vals, _nearest_dense(A, sigma, k), rtol=1e-9, atol=1e-12)


def test_eigenvalues_are_sorted_by_distance_from_the_shift() -> None:
    A = _complex_symmetric(200, seed=1)
    sigma = 11.0 - 1.0j
    vals, _ = ShiftInvertEigs(A, k=5).near(sigma)
    d = np.abs(vals - sigma)
    assert np.all(np.diff(d) >= 0.0)


def test_eigenvectors_match_dense_up_to_scale() -> None:
    """Compared under the c-norm sqrt(v@v): the right notion of 'same vector up
    to scale' for a complex-symmetric operator (v@v != 0 away from an
    exceptional point)."""
    A = _complex_symmetric(150, seed=2)
    sigma = 8.5 + 1.5j
    vals, vecs = ShiftInvertEigs(A, k=3).near(sigma)
    w, V = np.linalg.eig(A.toarray())
    for i, val in enumerate(vals):
        j = int(np.argmin(np.abs(w - val)))
        u, v = V[:, j], vecs[:, i]
        u = u / np.sqrt(c_product(u, u))
        v = v / np.sqrt(c_product(v, v))
        assert abs(abs(c_product(u, v)) - 1.0) < 1e-6


def test_shift_sign_convention_is_A_minus_sigma_I() -> None:
    """A spectrum deliberately asymmetric about sigma: passing `sigma*I - A` as
    OPinv instead of `A - sigma*I` returns eigenvalues reflected about sigma,
    which this pins down."""
    diag = np.array([0.0, 1.0, 2.0, 3.0, 10.0, 11.0, 12.0], dtype=complex)
    A = sp.csc_matrix(sp.diags(diag, dtype=complex))
    sigma = 9.0 + 0.0j
    vals, _ = ShiftInvertEigs(A, k=2).near(sigma)
    assert np.allclose(np.sort_complex(vals), np.array([10.0, 11.0], dtype=complex))


def test_k_too_large_raises_and_points_at_the_dense_route() -> None:
    A = _complex_symmetric(20, seed=3)
    with pytest.raises(ValueError, match="qscat.dvr.eigen"):
        ShiftInvertEigs(A, k=19).near(1.0 + 0.0j)


def test_non_square_raises() -> None:
    A = sp.csc_matrix((5, 7), dtype=complex)
    with pytest.raises(ValueError, match="square"):
        ShiftInvertEigs(A)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run --no-sync pytest libs/qscat/tests/test_shift_invert_eigs.py -q`
Expected: collection error / `ImportError: cannot import name 'ShiftInvertEigs' from 'qscat.linalg'`.

- [ ] **Step 3: Write the implementation**

Create `libs/qscat/qscat/linalg/eigs.py`:

```python
"""Sparse shift-invert eigensolver: interior eigenvalues near a complex shift.

Krylov methods converge to EXTREMAL eigenvalues, but a resonance is an INTERIOR
one -- it sits in the lower half of the complex plane surrounded by the
discretized (rotated) continuum. The standard remedy is the shift-invert
spectral transform: run Arnoldi on `(A - sigma*I)^-1`, whose extremal
eigenvalues are `A`'s eigenvalues nearest `sigma`.

That transform costs one sparse solve per matrix-vector product, which is what
`SparseLU` provides -- and because `A - sigma*I` has the SAME sparsity pattern
for every `sigma` (a diagonal shift), `SparseLU.refactor` reuses the symbolic
analysis across shifts. A shift sweep therefore gets the same discount as the
time-independent energy sweep (see `docs/physics/ti-energy-sweep-reuse.md`).

See `docs/physics/shift-invert-eigensolver.md`.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from qscat.exceptions import ConvergenceError

from .sparse_lu import SparseLU, _Backend, _Ordering

__all__ = ["ShiftInvertEigs"]


class ShiftInvertEigs:
    """The `k` eigenpairs of a sparse matrix nearest a complex shift.

    Holds the factorization: the first `near` call factors `A - sigma*I`, and
    every later call `refactor`s it, reusing the symbolic analysis (the shifted
    matrix's sparsity pattern does not depend on `sigma`). On the MUMPS backend
    that skips the SCOTCH ordering; on scipy it re-runs `splu` (correct, no
    reuse) -- the same trade `SparseLU.refactor` documents.

    `ordering`, `backend` and `symmetric` are forwarded verbatim to `SparseLU`.
    `k`, `ncv`, `tol` and `maxiter` are ARPACK controls; `k` may be overridden
    per call.

    Parameters
    ----------
    A : scipy.sparse.spmatrix
        Square matrix. Promoted to complex CSC internally.
    k : int, optional
        Number of eigenpairs to return per call (default 6).
    ordering : {"COLAMD", "NATURAL", "MMD_ATA", "MMD_AT_PLUS_A"}, optional
        SuperLU column ordering, forwarded to `SparseLU`.
    backend : {"auto", "scipy", "mumps"}, optional
        Factorization backend, forwarded to `SparseLU`.
    symmetric : bool or None, optional
        Complex-symmetry flag, forwarded to `SparseLU` (auto-detected if None).
    ncv : int or None, optional
        Krylov subspace size. ARPACK's default is used when None.
    tol : float, optional
        ARPACK relative tolerance; 0.0 means machine precision.
    maxiter : int or None, optional
        Maximum ARPACK restarts.

    Notes
    -----
    Eigenvectors are returned with numpy's Euclidean (`v^dagger v = 1`)
    normalization, exactly as `qscat.dvr.eigen` returns them. For
    exterior-complex-scaling observables re-normalize under the bilinear
    `c_product`.
    """

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
    ) -> None:
        if A.shape[0] != A.shape[1]:
            raise ValueError(f"matrix must be square, got shape {A.shape}")
        self._A: sp.csc_matrix[np.complex128] = sp.csc_matrix(A, dtype=np.complex128)
        self._n = int(self._A.shape[0])
        self._eye = sp.identity(self._n, format="csc", dtype=np.complex128)
        self._k = int(k)
        self._ordering: _Ordering = ordering
        self._backend: _Backend = backend
        self._symmetric = symmetric
        self._ncv = ncv
        self._tol = float(tol)
        self._maxiter = maxiter
        self._lu: SparseLU | None = None
        self._n_factorizations = 0

    def near(
        self, sigma: complex, *, k: int | None = None
    ) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
        """The `k` eigenpairs of `A` nearest `sigma`, nearest first.

        Returns `(energies, vectors)` with `energies` sorted by
        `|E - sigma|` ascending and `vectors[:, i]` the eigenvector of
        `energies[i]`.

        Raises
        ------
        ValueError
            If `k` is not in `1 <= k < n - 1` (an ARPACK constraint).
        qscat.exceptions.ConvergenceError
            If ARPACK fails to converge.
        """
        k_eff = self._k if k is None else int(k)
        if k_eff < 1:
            raise ValueError(f"k must be >= 1, got {k_eff}")
        if k_eff >= self._n - 1:
            raise ValueError(
                f"k={k_eff} must be < n-1 = {self._n - 1} (an ARPACK constraint); "
                "for a matrix this small use the dense qscat.dvr.eigen instead"
            )

        shift = complex(sigma)
        # scipy's convention: OPinv solves (A - sigma*I) x = b. NOT the driven
        # solver's (E*I - H); passing that instead silently returns eigenvalues
        # reflected about sigma. Adding the identity explicitly also forces the
        # diagonal to be structurally present, so the pattern -- and hence the
        # symbolic analysis reused by `refactor` -- is identical for every shift.
        shifted = (self._A - shift * self._eye).tocsc()
        if self._lu is None:
            self._lu = SparseLU(
                shifted,
                ordering=self._ordering,
                backend=self._backend,
                symmetric=self._symmetric,
            )
        else:
            self._lu.refactor(shifted)
        self._n_factorizations += 1

        op_inv = spla.LinearOperator(
            (self._n, self._n), matvec=self._lu.solve, dtype=np.complex128
        )
        try:
            vals, vecs = spla.eigs(
                self._A,
                k=k_eff,
                sigma=shift,
                OPinv=op_inv,
                ncv=self._ncv,
                tol=self._tol,
                maxiter=self._maxiter,
            )
        except spla.ArpackNoConvergence as exc:
            raise ConvergenceError(
                f"ARPACK did not converge at sigma={shift!r} with k={k_eff}: "
                f"{np.size(exc.eigenvalues)} of {k_eff} eigenvalues converged. "
                "Raise k, raise ncv, or move sigma closer to the target."
            ) from exc

        order = np.argsort(np.abs(vals - shift))
        return (
            np.asarray(vals[order], dtype=np.complex128),
            np.asarray(vecs[:, order], dtype=np.complex128),
        )

    @property
    def shape(self) -> tuple[int, int]:
        return (self._n, self._n)

    @property
    def n_factorizations(self) -> int:
        """How many times a shifted matrix has been factored (analysis + refactors)."""
        return self._n_factorizations
```

- [ ] **Step 4: Wire the export**

In `libs/qscat/qscat/linalg/__init__.py`: add `from .eigs import ShiftInvertEigs`, add `"ShiftInvertEigs"` to `__all__`, and add one bullet to the module docstring's "Public API" list:

```
  - `ShiftInvertEigs` -- the k eigenpairs nearest a complex shift (sparse
    shift-invert Arnoldi on top of `SparseLU`, reusing its analysis across
    shifts).
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run --no-sync pytest libs/qscat/tests/test_shift_invert_eigs.py -q`
Expected: 6 passed.

- [ ] **Step 6: Gates and commit**

```bash
uv run ruff check libs/qscat && uv run ruff format libs/qscat/qscat/linalg/eigs.py libs/qscat/tests/test_shift_invert_eigs.py
uv run mypy libs/qscat/qscat
git add libs/qscat/qscat/linalg/eigs.py libs/qscat/qscat/linalg/__init__.py libs/qscat/tests/test_shift_invert_eigs.py
git commit -m "feat(linalg): ShiftInvertEigs -- k eigenpairs nearest a complex shift"
```

---

### Task 2: Factorization reuse across shifts, diagnostics, and non-convergence

**Files:**
- Modify: `libs/qscat/qscat/linalg/eigs.py`
- Test: `libs/qscat/tests/test_shift_invert_eigs.py` (append)

**Interfaces:**
- Consumes: `ShiftInvertEigs` from Task 1; `SparseLU.backend_used`, `.ordering_used`, `.fill_factor`, `.memory_bytes()`.
- Produces: `backend_used`, `ordering_used`, `fill_factor` properties and `memory_bytes()` on `ShiftInvertEigs`; the guarantee that repeated `near()` calls agree with fresh objects.

- [ ] **Step 1: Write the failing tests** (append to the test file)

```python
def test_repeated_shifts_reuse_the_factorization_object() -> None:
    """Reuse must be invisible in the results and visible in the diagnostics."""
    A = _complex_symmetric(200, seed=4)
    s1, s2 = 9.0 + 2.0j, 12.0 - 1.0j
    solver = ShiftInvertEigs(A, k=4)
    v1, _ = solver.near(s1)
    v2, _ = solver.near(s2)
    assert solver.n_factorizations == 2
    fresh1, _ = ShiftInvertEigs(A, k=4).near(s1)
    fresh2, _ = ShiftInvertEigs(A, k=4).near(s2)
    assert np.allclose(v1, fresh1, rtol=1e-9, atol=1e-12)
    assert np.allclose(v2, fresh2, rtol=1e-9, atol=1e-12)


def test_diagnostics_delegate_to_the_factorization() -> None:
    A = _complex_symmetric(100, seed=5)
    solver = ShiftInvertEigs(A, k=3)
    with pytest.raises(RuntimeError, match="near"):
        _ = solver.backend_used                     # nothing factored yet
    solver.near(9.0 + 1.0j)
    assert solver.backend_used in {"scipy", "mumps"}
    assert solver.ordering_used
    assert solver.fill_factor > 0.0
    assert solver.memory_bytes() > 0
    assert solver.shape == (100, 100)


def test_non_convergence_raises_convergence_error() -> None:
    """maxiter=1 is far too few restarts: ARPACK bails, and we translate."""
    from qscat.exceptions import ConvergenceError

    A = _complex_symmetric(300, seed=6)
    solver = ShiftInvertEigs(A, k=8, ncv=10, maxiter=1)
    with pytest.raises(ConvergenceError, match="ARPACK"):
        solver.near(0.0 + 0.0j)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run --no-sync pytest libs/qscat/tests/test_shift_invert_eigs.py -q -k "reuse or diagnostics or convergence"`
Expected: FAIL — `AttributeError: 'ShiftInvertEigs' object has no attribute 'backend_used'`.

- [ ] **Step 3: Add the diagnostics**

Append to `ShiftInvertEigs` in `libs/qscat/qscat/linalg/eigs.py`:

```python
    def _require_lu(self) -> SparseLU:
        if self._lu is None:
            raise RuntimeError(
                "no factorization yet -- call near(sigma) before reading diagnostics"
            )
        return self._lu

    @property
    def backend_used(self) -> str:
        """Which factorization engine ran (`"scipy"` or `"mumps"`)."""
        return self._require_lu().backend_used

    @property
    def ordering_used(self) -> str:
        return self._require_lu().ordering_used

    @property
    def fill_factor(self) -> float:
        return self._require_lu().fill_factor

    def memory_bytes(self) -> int:
        """Factor memory. NOT cheap -- a method, not a property, so the cost is opt-in."""
        return self._require_lu().memory_bytes()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run --no-sync pytest libs/qscat/tests/test_shift_invert_eigs.py -q`
Expected: 9 passed.

If `test_non_convergence_raises_convergence_error` does NOT raise (ARPACK converged anyway), tighten it: raise `k` to 20 with `ncv=22, maxiter=1`. Do not delete the test — the translation path must be covered.

- [ ] **Step 5: Gates and commit**

```bash
uv run ruff check libs/qscat && uv run mypy libs/qscat/qscat
uv run --no-sync pytest libs/qscat/tests/test_shift_invert_eigs.py -q
git add libs/qscat/qscat/linalg/eigs.py libs/qscat/tests/test_shift_invert_eigs.py
git commit -m "feat(linalg): reuse the analysis across shifts + delegated diagnostics"
```

---

### Task 3: Physics differential test — the real 1-D electronic ECS Hamiltonian

**Files:**
- Test: `libs/qscat/tests/test_shift_invert_eigs.py` (append)

**Interfaces:**
- Consumes: `ShiftInvertEigs.near`; `qscat.core.grids.electronic_grid`; `qscat.dvr.kinetic`, `kinetic_sparse`, `eigen`; `qscat.ecs.find_resonance_pole`; `qscat.model.N2` (`N2.surface(points, R)` is the electronic potential surface at fixed nuclear `R`).
- Produces: nothing importable — this is the differential gate that makes the primitive trustworthy on a real ECS spectrum.

Context for the implementer: the electronic Hamiltonian at fixed `R` is `T(mass=1) + diag(V_surface)`, exactly as `qscat.core.lcp._h_el` builds it. The dense route (`eigen` + `find_resonance_pole` across two ECS angles) is the repo's validated pole finder; `docs/physics/n2-resonance.md` records `E_res ≈ 2.445 eV`, `Gamma ≈ 0.455 eV` at equilibrium, i.e. a pole near `0.090 - 0.0084j` Ha — which is where the seed shift below comes from.

- [ ] **Step 1: Write the failing test** (append)

```python
def _n2_electronic(angle_deg: float, R: float = 2.02):
    """Dense and sparse builds of the SAME fixed-R electronic ECS Hamiltonian."""
    from qscat.core.grids import electronic_grid
    from qscat.dvr import kinetic, kinetic_sparse
    from qscat.model import N2

    g = electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=angle_deg)
    V = np.asarray(N2.surface(g.points, R), dtype=np.complex128)
    H_dense = kinetic(g, 1.0) + np.diag(V)
    H_sparse = sp.csc_matrix(kinetic_sparse(g, 1.0) + sp.diags(V))
    return g, H_dense, H_sparse


def test_sparse_reproduces_the_dense_electronic_resonance_pole() -> None:
    """The differential oracle: on the real N2 electronic ECS Hamiltonian, the
    sparse shift-invert must return the same pole the dense path finds."""
    from qscat.dvr import eigen
    from qscat.ecs import find_resonance_pole

    _, Hd_a, Hs_a = _n2_electronic(35.0)
    _, Hd_b, _ = _n2_electronic(44.0)
    # Dense, two-angle: the repo's validated pole finder.
    E_pole, residual = find_resonance_pole(
        eigen(Hd_a)[0], eigen(Hd_b)[0], (0.0, 0.3, -0.1, 0.0)
    )
    assert residual < 1e-3                       # a genuine angle-stable pole

    # Sparse, one angle, seeded from a DELIBERATELY OFFSET physical guess.
    sigma = 0.10 - 0.010j
    vals, _ = ShiftInvertEigs(Hs_a, k=8).near(sigma)
    hit = vals[int(np.argmin(np.abs(vals - E_pole)))]
    assert abs(hit - E_pole) <= 1e-9 * abs(E_pole) + 1e-12


def test_sparse_eigenvectors_match_dense_on_the_electronic_hamiltonian() -> None:
    from qscat.dvr import eigen

    _, Hd, Hs = _n2_electronic(35.0)
    sigma = 0.10 - 0.010j
    vals, vecs = ShiftInvertEigs(Hs, k=4).near(sigma)
    w, V = eigen(Hd)
    for i, val in enumerate(vals):
        j = int(np.argmin(np.abs(w - val)))
        assert abs(w[j] - val) <= 1e-9 * abs(val) + 1e-12
        u = V[:, j] / np.sqrt(c_product(V[:, j], V[:, j]))
        v = vecs[:, i] / np.sqrt(c_product(vecs[:, i], vecs[:, i]))
        assert abs(abs(c_product(u, v)) - 1.0) < 1e-6
```

- [ ] **Step 2: Run to verify it fails, then passes**

Run: `uv run --no-sync pytest libs/qscat/tests/test_shift_invert_eigs.py -q -k electronic`
Expected first: FAIL only if something is wrong — these tests exercise Task 1 code, so they should PASS immediately. If they fail, the failure is the finding: report it before changing tolerances. Legitimate adjustments: widening the `find_resonance_pole` window if it raises "window contains no eigenvalues"; raising `k` if the pole is not among the 8 returned (record the value that works — it is a data point for Task 4's probe). Loosening `1e-9` is NOT a legitimate adjustment without a stated reason.

- [ ] **Step 3: Commit**

```bash
uv run ruff check libs/qscat && uv run --no-sync pytest libs/qscat/tests/test_shift_invert_eigs.py -q
git add libs/qscat/tests/test_shift_invert_eigs.py
git commit -m "test(linalg): differential gate on the N2 electronic ECS Hamiltonian"
```

---

### Task 4: Robustness probe + physics note + CHANGELOG

**Files:**
- Test: `libs/qscat/tests/test_shift_invert_eigs.py` (append)
- Create: `docs/physics/shift-invert-eigensolver.md`
- Modify: `CHANGELOG.md` (Unreleased → Added)

**Interfaces:**
- Consumes: `_n2_electronic` from Task 3; `ShiftInvertEigs`; `qscat.ecs.match_angle_stable`.
- Produces: the measured numbers stage 2 starts from — how far the shift may sit from the pole, the smallest workable `k`, and what a continuum-adjacent shift returns.

- [ ] **Step 1: Write the probe tests** (append)

```python
@pytest.mark.parametrize("offset", [0.001, 0.01, 0.05, 0.2])
def test_pole_is_found_from_an_offset_shift(offset: float) -> None:
    """How far may the seed shift sit from the pole and still find it at k=8?
    The largest passing offset is the number quoted in the physics note."""
    from qscat.dvr import eigen
    from qscat.ecs import find_resonance_pole

    _, Hd_a, Hs_a = _n2_electronic(35.0)
    _, Hd_b, _ = _n2_electronic(44.0)
    E_pole, _ = find_resonance_pole(
        eigen(Hd_a)[0], eigen(Hd_b)[0], (0.0, 0.3, -0.1, 0.0)
    )
    vals, _ = ShiftInvertEigs(Hs_a, k=8).near(E_pole + offset * (1.0 + 1.0j))
    assert np.min(np.abs(vals - E_pole)) <= 1e-9 * abs(E_pole) + 1e-12


def test_continuum_adjacent_shift_returns_angle_unstable_eigenvalues() -> None:
    """The 1-D rehearsal of the selection stage 2 needs: a shift parked on the
    rotated continuum returns eigenvalues that MOVE when the ECS angle changes,
    while the pole does not. Nothing here selects; it demonstrates the signal."""
    from qscat.ecs import match_angle_stable

    _, _, Hs_a = _n2_electronic(35.0)
    _, _, Hs_b = _n2_electronic(44.0)
    sigma = 0.6 - 0.25j                            # deep in the rotated continuum
    va, _ = ShiftInvertEigs(Hs_a, k=10).near(sigma)
    vb, _ = ShiftInvertEigs(Hs_b, k=10).near(sigma)
    window = (0.0, 1.5, -1.0, 0.0)
    stable, _, _ = match_angle_stable(va, vb, window, rel_tol=1e-4, atol=1e-8)
    assert stable.size == 0                        # continuum: nothing is stable
```

- [ ] **Step 2: Run the probes and RECORD the outcomes**

Run: `uv run --no-sync pytest libs/qscat/tests/test_shift_invert_eigs.py -q -k "offset or continuum" -v`

Write down: which offsets pass, the smallest `k` that still captures the pole at `offset=0.05` (try `k` = 2, 4, 6, 8 by hand), and whether the continuum probe returned an empty stable set. If a parametrized offset fails, do NOT delete it — change the parametrize list to the measured boundary and say so in the note. A failing large offset is the expected, reportable result: it is the primitive's working range.

- [ ] **Step 3: Write the physics note**

Create `docs/physics/shift-invert-eigensolver.md` covering, in this order:
1. **Why shift-invert** — resonances are interior eigenvalues; Krylov converges to extremal ones; the transform `(A - sigma*I)^-1`.
2. **Why it is cheap here** — `A - sigma*I` keeps the pattern for every `sigma`, so `SparseLU.refactor` reuses the analysis (MUMPS skips SCOTCH; scipy re-runs `splu`); cross-reference `docs/physics/ti-energy-sweep-reuse.md` and `docs/physics/mumps-sparse-backend.md`.
3. **Conventions** — the `A - sigma*I` sign convention and what goes wrong with the other one; nearest-shift-first ordering vs `dvr.eigen`'s `Re E` ordering; Euclidean eigenvector normalization and the `c_product` re-normalization callers must do.
4. **Validation** — the dense oracle, the `rtol = 1e-9` floor and why, and the N2 electronic-Hamiltonian result.
5. **Measured working range** — a small table of the Step 2 numbers: shift offset vs found/not-found, smallest workable `k`, and the continuum-adjacent behaviour.
6. **Limits** — near-singularity of `A - sigma*I` is the mechanism, not a bug, but MUMPS `SYM=2` may warn or lose accuracy for a shift pathologically close to an eigenvalue; state whatever was observed, or state that it was not observed at the offsets tested. No 2-D claims: this is validated in 1-D only.

- [ ] **Step 4: CHANGELOG**

Add under `## [Unreleased]` → `### Added`:

```markdown
- `qscat.linalg.ShiftInvertEigs`: the `k` eigenpairs of a sparse complex-symmetric
  matrix nearest a complex shift (shift-invert Arnoldi with `SparseLU` as the
  inner solve). Because `A - sigma*I` keeps its sparsity pattern for every shift,
  repeated `near(sigma)` calls reuse the symbolic analysis via
  `SparseLU.refactor` — the eigenvalue analogue of the time-independent energy
  sweep. Validated in 1-D against dense `np.linalg.eig` on both synthetic
  complex-symmetric matrices and the real N₂ electronic FEM-DVR-ECS Hamiltonian,
  eigenvalues and eigenvectors; see docs/physics/shift-invert-eigensolver.md for
  the measured working range (shift offset, `k`, continuum-adjacent behaviour).
```

- [ ] **Step 5: Full gates and commit**

```bash
uv run ruff check . && uv run mypy libs/qscat/qscat
uv run --no-sync pytest libs/qscat/tests -m "not slow" -n 8 -q
uv run --no-sync sphinx-build -b html -W --keep-going docs docs/_build/html
git add libs/qscat/tests/test_shift_invert_eigs.py docs/physics/shift-invert-eigensolver.md CHANGELOG.md
git commit -m "docs(linalg): shift-invert eigensolver note + measured working range"
```

Expected: the full fast suite green (374 passed + the new tests), Sphinx build succeeded.

---

## Self-review notes

- **Spec coverage:** interface → Task 1; reuse + diagnostics + `ConvergenceError` → Task 2; synthetic differential + sign pin + `k` guard → Task 1; physics differential → Task 3; robustness probe, note, CHANGELOG → Task 4. Every spec section maps to a task.
- **Out-of-scope items stay out:** no 2-D use, no `qscat.core` edits, no two-angle selection API (Task 4 only *demonstrates* the signal with the existing `match_angle_stable`), no Lanczos, no MUMPS tuning.
- **Naming is consistent across tasks:** `ShiftInvertEigs`, `near(sigma, k=)`, `n_factorizations`, `backend_used`, `ordering_used`, `fill_factor`, `memory_bytes()`, `_n2_electronic`.
