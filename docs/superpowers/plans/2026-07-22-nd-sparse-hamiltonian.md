# N-Dimensional Sparse Hamiltonian Library Implementation Plan (sub-project #5)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the dimension-general tensor-grid + Kronecker-sum sparse Hamiltonian layer in `qscat`, validated against analytic benchmarks at D = 1, 2, 3, so the exact 2-D N₂ solver (sub-project #6) rests on independently-proven numerics.

**Architecture:** Two layers with a strict dependency direction (`qscat.dvr` → `qscat.linalg`, never the reverse). **Layer 0** (`qscat.linalg`) is pure sparse linear algebra with no knowledge of grids or physics: `kron_sum`, `SparseLU`, `c_product`. **Layer 1** (`qscat.dvr`) adds grid-aware assembly: `kinetic_sparse`, `TensorGrid`, `kinetic_nd`/`potential_nd`/`hamiltonian_nd`. No N₂ physics is touched in this sub-project.

**Tech Stack:** Python 3.12, numpy >= 2, scipy >= 1.14 (`scipy.sparse`, `scipy.sparse.linalg.splu`), pytest.

**Design spec:** `docs/superpowers/specs/2026-07-22-nd-sparse-hamiltonian-design.md`
**Source archaeology:** `.superpowers/sdd/n2-2d-exact-extraction.md`

## Global Constraints

- Python `>=3.12`. Run everything through `uv` — **never** bare `python`/`pip`/conda. Conda `base` may be active in the shell; `uv run` correctly ignores it and uses `.venv/bin/python3` (3.12.7).
- Atomic units throughout (Hartree, bohr). No ad hoc unit conversions.
- Package-absolute imports only (`from qscat.linalg import kron_sum`). No `sys.path` hacks, no `importlib.util.spec_from_file_location`.
- `uv run mypy libs/qscat` must stay at **0 errors** (strict mode). Annotate every public signature.
- `uv run ruff check .` must stay clean. Line length 100; lint rules `E, F, I, UP, B, NPY`.
- **The existing N₂ harness must not regress:** `uv run python -m validation.n2.experiment` stays at **19 PASS / 0 PENDING / 2 NOTE / 0 FAIL**, exit 0.
- **Index ordering is numpy-native C-order: LAST axis fastest.** `psi[i_0, …, i_{D-1}].ravel()` pairs with `kron_sum`. eMoScat uses the opposite convention (first coordinate fastest); this is a deliberate, documented divergence.
- ECS makes `H = Hᵀ ≠ H†` — **complex symmetric, not Hermitian**. Never use a Hermitian-only algorithm (`eigh`, `vdot`-based norms) on these matrices.
- Commit trailer on every commit: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`

---

## File Structure

| File | Responsibility |
|---|---|
| `libs/qscat/pyproject.toml` | Add the missing `scipy` runtime dependency (Task 1) |
| `pyproject.toml` | Register the `slow` pytest marker (Task 6) |
| `libs/qscat/qscat/linalg/kron.py` | **Create.** `kron_sum` — Kronecker sum over arbitrary D |
| `libs/qscat/qscat/linalg/inner.py` | **Create.** `c_product` — the no-conjugate ECS inner product |
| `libs/qscat/qscat/linalg/sparse_lu.py` | **Create.** `SparseLU` — cached factorization + diagnostics |
| `libs/qscat/qscat/linalg/__init__.py` | **Modify.** Export the above |
| `libs/qscat/qscat/dvr/kinetic.py` | **Modify.** Add `kinetic_sparse` alongside the dense `kinetic` |
| `libs/qscat/qscat/dvr/tensor.py` | **Create.** `TensorGrid`, `kinetic_nd`, `potential_nd`, `hamiltonian_nd` |
| `libs/qscat/qscat/dvr/__init__.py` | **Modify.** Export the above |
| `libs/qscat/tests/test_kron_sum.py` | **Create.** V1 |
| `libs/qscat/tests/test_sparse_lu.py` | **Create.** V5 |
| `libs/qscat/tests/test_kinetic_sparse.py` | **Create.** V2 |
| `libs/qscat/tests/test_tensor_grid.py` | **Create.** `TensorGrid` geometry/mask/outer |
| `libs/qscat/tests/test_hamiltonian_nd.py` | **Create.** V3 (analytic benchmarks), V4 (D=1 regression) |
| `libs/qscat/tests/test_nd_scale.py` | **Create.** V6, `@pytest.mark.slow` |
| `docs/physics/nd-tensor-hamiltonian.md` | **Create.** The method note |
| `CLAUDE.md` | **Modify.** Document `qscat.linalg` + the `qscat.dvr` tensor additions |

---

### Task 1: `qscat.linalg.kron_sum` and `c_product`

**Files:**
- Modify: `libs/qscat/pyproject.toml`
- Create: `libs/qscat/qscat/linalg/kron.py`, `libs/qscat/qscat/linalg/inner.py`
- Modify: `libs/qscat/qscat/linalg/__init__.py`
- Test: `libs/qscat/tests/test_kron_sum.py`

**Interfaces:**
- Consumes: nothing (first task).
- Produces:
  - `kron_sum(ops: Sequence[sparse.spmatrix]) -> sparse.csr_matrix`
  - `c_product(a: npt.ArrayLike, b: npt.ArrayLike) -> complex`

**Background you need:** `qscat` currently declares only `numpy>=2` as a runtime dependency, but `qscat/evolution/crank_nicolson.py` already imports `scipy.linalg`. That is a latent packaging bug — it works only because scipy is in the repo's dev group. This task fixes it, since `qscat.linalg` will depend on scipy properly.

- [ ] **Step 1: Fix the scipy dependency**

In `libs/qscat/pyproject.toml`, change:

```toml
dependencies = ["numpy>=2"]
```

to:

```toml
dependencies = ["numpy>=2", "scipy>=1.14"]
```

- [ ] **Step 2: Write the failing tests — `libs/qscat/tests/test_kron_sum.py`**

```python
"""Tests for `qscat.linalg.kron_sum` (V1) and `qscat.linalg.c_product`.

`kron_sum` is checked against dense `np.kron` at D = 1, 2, 3, 4 with UNEQUAL
per-axis dimensions -- square cases hide transposed-index bugs.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from qscat.linalg import c_product, kron_sum


def _dense_kron_sum(mats: list[np.ndarray]) -> np.ndarray:
    """Reference: sum_d I x ... x mats[d] x ... x I, built with dense np.kron."""
    sizes = [m.shape[0] for m in mats]
    total = int(np.prod(sizes))
    out = np.zeros((total, total), dtype=complex)
    for d, m in enumerate(mats):
        term = np.eye(1, dtype=complex)
        for e, n in enumerate(sizes):
            term = np.kron(term, m if e == d else np.eye(n, dtype=complex))
        out += term
    return out


def _random_mats(rng: np.random.Generator, sizes: list[int]) -> list[np.ndarray]:
    return [
        rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)) for n in sizes
    ]


@pytest.mark.parametrize("sizes", [[4], [3, 4], [2, 3, 4], [2, 3, 2, 3]])
def test_kron_sum_matches_dense_np_kron(sizes: list[int]) -> None:
    rng = np.random.default_rng(0)
    mats = _random_mats(rng, sizes)
    got = kron_sum([sp.csr_matrix(m) for m in mats]).toarray()
    want = _dense_kron_sum(mats)
    assert got.shape == (int(np.prod(sizes)),) * 2
    assert np.allclose(got, want, rtol=0, atol=1e-12)


def test_kron_sum_single_operator_is_identity_operation() -> None:
    rng = np.random.default_rng(1)
    (m,) = _random_mats(rng, [5])
    assert np.allclose(kron_sum([sp.csr_matrix(m)]).toarray(), m, rtol=0, atol=1e-14)


def test_kron_sum_acts_on_c_order_ravel() -> None:
    """LAST axis fastest: (A (x) I + I (x) B) vec(psi) == vec(A@psi + psi@B.T)."""
    rng = np.random.default_rng(2)
    A, B = _random_mats(rng, [3, 4])
    psi = rng.standard_normal((3, 4)) + 1j * rng.standard_normal((3, 4))
    got = kron_sum([sp.csr_matrix(A), sp.csr_matrix(B)]) @ psi.ravel()
    want = (A @ psi + psi @ B.T).ravel()
    assert np.allclose(got, want, rtol=0, atol=1e-12)


def test_kron_sum_rejects_empty() -> None:
    with pytest.raises(ValueError, match="at least one"):
        kron_sum([])


def test_kron_sum_rejects_non_square() -> None:
    with pytest.raises(ValueError, match="square"):
        kron_sum([sp.csr_matrix(np.zeros((2, 3)))])


def test_c_product_does_not_conjugate() -> None:
    """The whole point: c_product != vdot for complex vectors."""
    a = np.array([1j, 2.0])
    assert c_product(a, a) == pytest.approx(3.0 + 0j)   # (1j)^2 + 4 = 3
    assert np.vdot(a, a) == pytest.approx(5.0 + 0j)     # |1j|^2 + 4 = 5, NOT what we want


def test_c_product_is_symmetric() -> None:
    rng = np.random.default_rng(3)
    a = rng.standard_normal(6) + 1j * rng.standard_normal(6)
    b = rng.standard_normal(6) + 1j * rng.standard_normal(6)
    assert c_product(a, b) == pytest.approx(c_product(b, a))
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `uv run pytest libs/qscat/tests/test_kron_sum.py -q`
Expected: FAIL — `ImportError: cannot import name 'c_product' from 'qscat.linalg'`

- [ ] **Step 4: Implement `libs/qscat/qscat/linalg/kron.py`**

```python
"""Kronecker sum of operators over an arbitrary number of dimensions.

The construction behind every separable-kinetic Hamiltonian on a tensor
product of grids:

    kron_sum([A_0, ..., A_{D-1}]) = sum_d  I x ... x A_d x ... x I

Pure linear algebra -- this module knows nothing about grids, potentials or
physics, and accepts ANY square sparse matrices. That is deliberate: a future
angular-DVR, finite-difference or B-spline dimension composes with FEM-DVR-ECS
dimensions at no extra cost.

Index convention: numpy-native C order, i.e. the LAST axis is fastest, so the
result acts on `psi.ravel()` for `psi` of shape `(n_0, ..., n_{D-1})`. Note
eMoScat uses the OPPOSITE convention (first coordinate fastest,
`idx = i_r + i_R*N_r`, `FemDvrEcsGrid2d.cpp:169`); the two are physically
identical and differ only in basis ordering.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import scipy.sparse as sp

__all__ = ["kron_sum"]


def kron_sum(ops: Sequence[sp.spmatrix]) -> sp.csr_matrix:
    """Assemble `sum_d I x ... x ops[d] x ... x I` as a CSR matrix.

    Each `ops[d]` must be square. The result is square with dimension
    `prod(n_d)`. `D == 1` returns `ops[0]` unchanged (as CSR).
    """
    mats = list(ops)
    if not mats:
        raise ValueError("kron_sum requires at least one operator")
    for d, m in enumerate(mats):
        if m.shape[0] != m.shape[1]:
            raise ValueError(f"operator {d} is not square: shape {m.shape}")

    sizes = [int(m.shape[0]) for m in mats]
    if len(mats) == 1:
        return sp.csr_matrix(mats[0])

    total: sp.csr_matrix | None = None
    for d, m in enumerate(mats):
        left = int(np.prod(sizes[:d])) if d else 1
        right = int(np.prod(sizes[d + 1 :])) if d < len(mats) - 1 else 1
        term = sp.kron(
            sp.identity(left, format="csr", dtype=complex),
            sp.kron(m, sp.identity(right, format="csr", dtype=complex), format="csr"),
            format="csr",
        )
        total = term if total is None else total + term

    assert total is not None  # len(mats) >= 2 guarantees at least one iteration
    return sp.csr_matrix(total)
```

- [ ] **Step 5: Implement `libs/qscat/qscat/linalg/inner.py`**

```python
"""The c-product: the bilinear inner product exterior complex scaling requires.

Under ECS the Hamiltonian is complex SYMMETRIC (`H = H^T`), not Hermitian, so
the natural pairing is the bilinear `sum_i a_i b_i` with NO complex conjugate
-- not `numpy.vdot`'s sesquilinear `sum_i conj(a_i) b_i`.

Getting this wrong is a recurring, quiet failure mode: it produces
plausible-looking complex "cross sections" with the wrong phase rather than an
obvious error. It has already bitten this repo once (sub-project #3's S-matrix,
where the Hermitian convention gave negative sigma), and the reference C++ code
uses `cblas_zdotc` here -- formally wrong, and correct in practice only because
it zeroes every channel function on the complex-scaled tail. Naming the
operation makes the choice explicit at every call site.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["c_product"]


def c_product(a: npt.ArrayLike, b: npt.ArrayLike) -> complex:
    """`sum_i a_i b_i` -- the bilinear (NOT conjugated) inner product."""
    av = np.asarray(a, dtype=np.complex128).ravel()
    bv = np.asarray(b, dtype=np.complex128).ravel()
    if av.shape != bv.shape:
        raise ValueError(f"shape mismatch: {av.shape} vs {bv.shape}")
    return complex(np.dot(av, bv))
```

- [ ] **Step 6: Update `libs/qscat/qscat/linalg/__init__.py`**

```python
"""Linear-algebra helpers: dimension-general Kronecker sums, cached sparse
factorizations, and the exterior-complex-scaling c-product.

Pure linear algebra -- nothing here knows about grids, potentials or physics,
so it composes with any discretization.

Public API:
  - `kron_sum` -- `sum_d I x ... x A_d x ... x I` for arbitrary D.
  - `c_product` -- the bilinear (non-conjugated) ECS inner product.

See `docs/physics/nd-tensor-hamiltonian.md`.
"""

from __future__ import annotations

from .inner import c_product
from .kron import kron_sum

__all__ = ["kron_sum", "c_product"]
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `uv run pytest libs/qscat/tests/test_kron_sum.py -q`
Expected: PASS (10 tests — 4 parametrized + 6)

- [ ] **Step 8: Type-check and lint**

Run: `uv run mypy libs/qscat && uv run ruff check .`
Expected: `Success: no issues found` and no ruff diagnostics.

- [ ] **Step 9: Commit**

```bash
git add libs/qscat/pyproject.toml libs/qscat/qscat/linalg libs/qscat/tests/test_kron_sum.py
git commit -m "$(cat <<'EOF'
feat(linalg): dimension-general kron_sum and the ECS c-product

kron_sum assembles sum_d I x ... x A_d x ... x I for arbitrary D over any
square sparse matrices, in numpy-native C order (last axis fastest).
Validated against dense np.kron at D=1,2,3,4 with unequal per-axis
dimensions, which square-only cases would not catch.

c_product names the bilinear, non-conjugated inner product that exterior
complex scaling requires, so the choice is explicit at every call site
rather than an implicit np.dot-vs-vdot trap.

Also declares qscat's scipy dependency, which crank_nicolson.py has been
importing without it being listed.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: `qscat.linalg.SparseLU`

**Files:**
- Create: `libs/qscat/qscat/linalg/sparse_lu.py`
- Modify: `libs/qscat/qscat/linalg/__init__.py`
- Test: `libs/qscat/tests/test_sparse_lu.py`

**Interfaces:**
- Consumes: nothing from Task 1 (independent), but lands in the same subpackage.
- Produces:
  - `SparseLU(A: sparse.spmatrix, *, ordering: str = "COLAMD")`
  - `.solve(b: npt.NDArray) -> npt.NDArray` — accepts `(N,)` or `(N, k)`
  - `.fill_factor: float`, `.memory_bytes: int`, `.shape: tuple[int, int]`

**Background you need:** `scipy.sparse.linalg.splu` requires CSC input and emits a `SparseEfficiencyWarning` otherwise, so convert explicitly. `ordering` maps to scipy's `permc_spec`; valid values are `"NATURAL"`, `"MMD_ATA"`, `"MMD_AT_PLUS_A"`, `"COLAMD"`. The diagnostics exist because at production sizes fill-in and memory — not flops — decide feasibility: a measured spike on the real N₂ 2-D deck gave ×93 fill-in and 13.6 GB peak RSS.

- [ ] **Step 1: Write the failing tests — `libs/qscat/tests/test_sparse_lu.py`**

```python
"""Tests for `qscat.linalg.SparseLU` (V5).

Exercised on a COMPLEX SYMMETRIC matrix (H = H^T, not Hermitian), which is what
exterior complex scaling produces and what every real use of this class will be.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
import scipy.sparse as sp
from qscat.linalg import SparseLU


def _complex_symmetric(n: int, seed: int) -> sp.csc_matrix:
    """A well-conditioned, sparse, complex-symmetric test matrix."""
    rng = np.random.default_rng(seed)
    nnz = 5 * n
    rows = rng.integers(0, n, size=nnz)
    cols = rng.integers(0, n, size=nnz)
    vals = rng.standard_normal(nnz) + 1j * rng.standard_normal(nnz)
    m = sp.coo_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex).tocsr()
    m = m + m.T                                             # complex SYMMETRIC, no conjugate
    m = m + sp.identity(n, format="csr", dtype=complex) * (10.0 + 3.0j)  # diagonally dominant
    return sp.csc_matrix(m)


def test_solve_residual_is_at_round_off() -> None:
    n = 200
    A = _complex_symmetric(n, seed=0)
    rng = np.random.default_rng(10)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    x = SparseLU(A).solve(b)
    assert np.linalg.norm(A @ x - b) / np.linalg.norm(b) < 1e-12


def test_matrix_is_complex_symmetric_not_hermitian() -> None:
    """Guard the fixture itself -- a Hermitian matrix would not exercise the point."""
    A = _complex_symmetric(50, seed=1)
    assert abs(A - A.T).max() < 1e-14
    assert abs(A - A.conj().T).max() > 1e-3


def test_multi_rhs_matches_looped_single_solves() -> None:
    n = 150
    A = _complex_symmetric(n, seed=2)
    rng = np.random.default_rng(11)
    B = rng.standard_normal((n, 4)) + 1j * rng.standard_normal((n, 4))
    lu = SparseLU(A)
    together = lu.solve(B)
    assert together.shape == (n, 4)
    for j in range(4):
        assert np.allclose(together[:, j], lu.solve(B[:, j]), rtol=0, atol=1e-12)


def test_factorization_is_reused_across_solves() -> None:
    """One factorization, many right-hand sides -- the whole reason this class exists."""
    n = 120
    A = _complex_symmetric(n, seed=3)
    lu = SparseLU(A)
    rng = np.random.default_rng(12)
    for _ in range(5):
        b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        assert np.linalg.norm(A @ lu.solve(b) - b) / np.linalg.norm(b) < 1e-12


def test_diagnostics_are_reported() -> None:
    A = _complex_symmetric(200, seed=4)
    lu = SparseLU(A)
    assert lu.shape == (200, 200)
    assert lu.fill_factor >= 1.0
    assert lu.memory_bytes > 0


def test_ordering_is_configurable_and_changes_fill() -> None:
    """Both orderings must solve correctly; fill-in generally differs."""
    n = 300
    A = _complex_symmetric(n, seed=5)
    rng = np.random.default_rng(13)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    for ordering in ("COLAMD", "MMD_AT_PLUS_A", "NATURAL"):
        lu = SparseLU(A, ordering=ordering)
        assert np.linalg.norm(A @ lu.solve(b) - b) / np.linalg.norm(b) < 1e-12


def test_rejects_non_square() -> None:
    with pytest.raises(ValueError, match="square"):
        SparseLU(sp.csc_matrix(np.zeros((3, 4))))


def test_accepts_csr_input_without_warning() -> None:
    """CSR in must be converted internally, not warned about."""
    A = sp.csr_matrix(_complex_symmetric(80, seed=6))
    rng = np.random.default_rng(14)
    b = rng.standard_normal(80) + 1j * rng.standard_normal(80)
    with warnings.catch_warnings():
        warnings.simplefilter("error", sp.SparseEfficiencyWarning)
        x = SparseLU(A).solve(b)
    assert np.linalg.norm(A @ x - b) / np.linalg.norm(b) < 1e-12
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest libs/qscat/tests/test_sparse_lu.py -q`
Expected: FAIL — `ImportError: cannot import name 'SparseLU' from 'qscat.linalg'`

- [ ] **Step 3: Implement `libs/qscat/qscat/linalg/sparse_lu.py`**

```python
"""Cached sparse LU factorization: factor once, solve many right-hand sides.

A thin, typed wrapper over `scipy.sparse.linalg.splu` that adds the two things
the bare function lacks for our use: an explicit CSC conversion (splu warns
otherwise), and fill-in / memory diagnostics.

Those diagnostics are not decoration. At the sizes this library targets, the
factorization -- not the solve -- is the whole cost, and fill-in decides whether
a problem fits in RAM at all. A measured spike on the production N2 2-D deck
(N = 143,380, nnz = 3,276,450) gave x93 fill-in, 3.05e8 nonzeros in L+U, and
13.6 GB peak RSS with the default COLAMD ordering, against a 440 ms
back-substitution. Choosing an ordering is therefore a real decision, and
`ordering` + `fill_factor` + `memory_bytes` exist so it can be MEASURED rather
than guessed. It cannot affect correctness -- only speed and memory.

Reusing one factorization across right-hand sides is the point: in a scattering
calculation every final channel at a given energy shares the same matrix.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

__all__ = ["SparseLU"]


class SparseLU:
    """LU factorization of a square sparse matrix, reusable across solves.

    `ordering` is scipy's `permc_spec`: one of `"NATURAL"`, `"MMD_ATA"`,
    `"MMD_AT_PLUS_A"`, `"COLAMD"` (the default). For a structurally symmetric
    pattern -- which a Kronecker-sum Hamiltonian has -- `"MMD_AT_PLUS_A"` is
    often the better choice; measure with `fill_factor` before assuming.
    """

    def __init__(self, A: sp.spmatrix, *, ordering: str = "COLAMD") -> None:
        if A.shape[0] != A.shape[1]:
            raise ValueError(f"matrix must be square, got shape {A.shape}")
        csc = sp.csc_matrix(A)
        self._shape: tuple[int, int] = (int(csc.shape[0]), int(csc.shape[1]))
        self._nnz: int = int(csc.nnz)
        self._ordering = ordering
        self._lu = spla.splu(csc, permc_spec=ordering)

    @property
    def shape(self) -> tuple[int, int]:
        return self._shape

    @property
    def ordering(self) -> str:
        return self._ordering

    @property
    def fill_factor(self) -> float:
        """`(L.nnz + U.nnz) / A.nnz` -- how much denser the factors are."""
        if self._nnz == 0:
            return 1.0
        return float(self._lu.L.nnz + self._lu.U.nnz) / float(self._nnz)

    @property
    def memory_bytes(self) -> int:
        """Bytes actually held by the L and U factors (data + index arrays)."""
        total = 0
        for factor in (self._lu.L, self._lu.U):
            csc = factor.tocsc()
            total += csc.data.nbytes + csc.indices.nbytes + csc.indptr.nbytes
        return int(total)

    def solve(self, b: npt.NDArray[np.complex128]) -> npt.NDArray[np.complex128]:
        """Solve `A x = b` for one `(N,)` or several `(N, k)` right-hand sides."""
        rhs = np.asarray(b)
        if rhs.shape[0] != self._shape[0]:
            raise ValueError(
                f"right-hand side has leading dimension {rhs.shape[0]}, "
                f"expected {self._shape[0]}"
            )
        out: npt.NDArray[np.complex128] = self._lu.solve(
            rhs.astype(np.complex128, copy=False)
        )
        return out
```

- [ ] **Step 4: Export it from `libs/qscat/qscat/linalg/__init__.py`**

Add `SparseLU` to the imports, `__all__`, and the module docstring's Public API list:

```python
from .inner import c_product
from .kron import kron_sum
from .sparse_lu import SparseLU

__all__ = ["kron_sum", "c_product", "SparseLU"]
```

and add this bullet to the docstring's Public API section:

```
  - `SparseLU` -- cached sparse LU factorization (factor once, solve many),
    with fill-in and memory diagnostics.
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest libs/qscat/tests/test_sparse_lu.py -q`
Expected: PASS (8 tests)

- [ ] **Step 6: Type-check and lint**

Run: `uv run mypy libs/qscat && uv run ruff check .`
Expected: clean.

- [ ] **Step 7: Commit**

```bash
git add libs/qscat/qscat/linalg libs/qscat/tests/test_sparse_lu.py
git commit -m "$(cat <<'EOF'
feat(linalg): SparseLU -- cached sparse factorization with fill diagnostics

Factor once, solve many right-hand sides: in a scattering calculation every
final channel at a given energy shares the same matrix. Tested on complex
SYMMETRIC (not Hermitian) matrices, which is what exterior complex scaling
produces.

fill_factor and memory_bytes are exposed because at the sizes this targets
the factorization is the entire cost and fill-in decides whether a problem
fits in RAM: the production N2 2-D deck factorizes to x93 fill-in and
13.6 GB with default COLAMD, against a 440 ms back-substitution. The
ordering choice is a measurement, not a guess.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: `qscat.dvr.kinetic_sparse`

**Files:**
- Modify: `libs/qscat/qscat/dvr/kinetic.py`
- Modify: `libs/qscat/qscat/dvr/__init__.py`
- Test: `libs/qscat/tests/test_kinetic_sparse.py`

**Interfaces:**
- Consumes: nothing from Tasks 1–2.
- Produces: `kinetic_sparse(grid: FemDvrEcsGrid, mass: float) -> sparse.csr_matrix`

**Background you need:** The existing dense `kinetic()` scatter-**adds** each element's retained sub-block into a dense array, and adjacent elements deliberately share one bridge global index so the `+=` accumulates the bridge coupling. The sparse version reproduces this exactly by emitting COO triplets: **`coo_matrix` sums duplicate `(row, col)` entries when converted to CSR**, which is precisely the same accumulation. Do not try to special-case bridges.

The nonzero count is predicted analytically by eMoScat's formula (`KineticEnergy.cpp:95`):

```
nnz = q**2 * tnel - 4*q + 3 - tnel
```

where `q = grid.nq` and `tnel = len(grid.spec.elements)`. This is the union of the per-element `q x q` blocks, overlapping by one index at each of the `tnel - 1` bridges, minus the two dropped Dirichlet endpoints (each removing `2q - 1` entries). It has been verified against the real N₂ grids: electronic `q=8, tnel=48 -> 2995`; nuclear `q=14, tnel=33 -> 6382`.

- [ ] **Step 1: Write the failing tests — `libs/qscat/tests/test_kinetic_sparse.py`**

```python
"""Tests for `qscat.dvr.kinetic_sparse` (V2).

The existing DENSE `kinetic()` -- already validated in sub-project #1 against
analytic particle-in-a-box and bound-state theta-independence -- is retained
specifically as the differential oracle for this sparse implementation.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec, kinetic, kinetic_sparse

CASES = {
    # name: (quadrature, elements)
    "all-real, uniform": (8, [ElementSpec(1.0) for _ in range(4)]),
    "all-real, graded": (6, [ElementSpec(0.5), ElementSpec(1.0), ElementSpec(2.0)]),
    "with ECS tail": (
        10,
        [ElementSpec(1.0), ElementSpec(1.0), ElementSpec(2.0, 35.0), ElementSpec(3.0, 35.0)],
    ),
    "single element": (7, [ElementSpec(1.5)]),
}


def _grid(name: str) -> FemDvrEcsGrid:
    q, els = CASES[name]
    return FemDvrEcsGrid(GridSpec(quadrature=q, elements=list(els), x_min=0.0))


@pytest.mark.parametrize("name", list(CASES))
@pytest.mark.parametrize("mass", [1.0, 12766.36])
def test_sparse_matches_dense_oracle(name: str, mass: float) -> None:
    grid = _grid(name)
    dense = kinetic(grid, mass)
    got = kinetic_sparse(grid, mass)
    assert got.shape == dense.shape
    scale = np.abs(dense).max()
    assert np.abs(got.toarray() - dense).max() <= 1e-12 * scale


# The formula below assumes the two dropped Dirichlet endpoints live in
# DIFFERENT elements, which requires tnel >= 2. For a single element both drops
# land in the same block and the count is instead (q-2)^2 -- i.e. the matrix is
# simply dense at n = q-2. The single-element grid is still covered by the
# differential test above, which is the check that actually matters.
@pytest.mark.parametrize("name", [n for n in CASES if n != "single element"])
def test_sparsity_matches_analytic_nnz_formula(name: str) -> None:
    """eMoScat KineticEnergy.cpp:95 -- nnz = q^2*tnel - 4q + 3 - tnel (tnel >= 2)."""
    grid = _grid(name)
    m = kinetic_sparse(grid, 1.0)
    m.eliminate_zeros()
    q = grid.nq
    tnel = len(grid.spec.elements)
    assert m.nnz == q**2 * tnel - 4 * q + 3 - tnel


def test_single_element_grid_is_dense() -> None:
    """tnel == 1: both Dirichlet drops fall in one block, giving a dense (q-2)^2."""
    grid = _grid("single element")
    m = kinetic_sparse(grid, 1.0)
    m.eliminate_zeros()
    assert grid.n == grid.nq - 2
    assert m.nnz == (grid.nq - 2) ** 2


def test_returns_csr_and_is_actually_sparse() -> None:
    grid = _grid("all-real, uniform")
    m = kinetic_sparse(grid, 1.0)
    assert isinstance(m, sp.csr_matrix)
    assert m.nnz < grid.n**2


def test_complex_symmetric_under_ecs() -> None:
    """ECS gives H = H^T but H != H^dagger; the kinetic term must already be so."""
    grid = _grid("with ECS tail")
    m = kinetic_sparse(grid, 1.0)
    assert abs(m - m.T).max() < 1e-12
    assert abs(m - m.conj().T).max() > 1e-6
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest libs/qscat/tests/test_kinetic_sparse.py -q`
Expected: FAIL — `ImportError: cannot import name 'kinetic_sparse' from 'qscat.dvr'`

- [ ] **Step 3: Implement `kinetic_sparse` in `libs/qscat/qscat/dvr/kinetic.py`**

Add these imports at the top of the file (keeping the existing ones):

```python
import scipy.sparse as sp
```

Change `__all__` to:

```python
__all__ = ["kinetic", "kinetic_sparse"]
```

Append this function to the end of the file:

```python
def kinetic_sparse(grid: FemDvrEcsGrid, mass: float) -> sp.csr_matrix:
    """Sparse (CSR) FEM-DVR kinetic-energy matrix -- the sparse sibling of `kinetic`.

    Identical mathematics to the dense `kinetic()`, which is retained as this
    function's differential oracle. The only structural difference is that
    per-element blocks are emitted as COO triplets instead of scatter-added into
    a dense array: `coo_matrix` SUMS duplicate `(row, col)` entries on
    conversion, which reproduces the dense version's `+=` bridge accumulation
    exactly. No bridge special-casing is needed or wanted.

    Nonzero count is `nq**2 * tnel - 4*nq + 3 - tnel` (eMoScat
    `KineticEnergy.cpp:95`) -- the union of the per-element `nq x nq` blocks,
    overlapping by one index at each bridge, less the two dropped Dirichlet
    endpoints.
    """
    n = grid.n
    nq = grid.nq

    _, wl = gll_nodes_weights(nq)  # reference GLL weights on (-1, 1)

    rows: list[npt.NDArray[np.intp]] = []
    cols: list[npt.NDArray[np.intp]] = []
    vals: list[npt.NDArray[np.complex128]] = []

    for k, (local, global_idx) in enumerate(grid.element_maps):
        hz = grid.hz[k]

        wze = hz * wl
        dBF = grid.dLp.T / hz

        norm = np.ones(nq, dtype=complex)
        norm[local] = 1.0 / np.sqrt(grid.weights[global_idx])
        dBF_n = dBF * norm[:, np.newaxis]

        T_local = (1.0 / (2.0 * mass)) * np.einsum("l,al,bl->ab", wze, dBF_n, dBF_n)
        block = T_local[np.ix_(local, local)]

        gi, gj = np.meshgrid(global_idx, global_idx, indexing="ij")
        rows.append(gi.ravel())
        cols.append(gj.ravel())
        vals.append(block.ravel())

    coo = sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n),
        dtype=complex,
    )
    return sp.csr_matrix(coo)  # duplicate (row, col) entries are summed here
```

- [ ] **Step 4: Export it from `libs/qscat/qscat/dvr/__init__.py`**

Change the import line and `__all__`:

```python
from .kinetic import kinetic, kinetic_sparse
```

```python
__all__ = [
    "ElementSpec",
    "GridSpec",
    "FemDvrEcsGrid",
    "kinetic",
    "kinetic_sparse",
    "hamiltonian",
    "eigen",
    "gll_nodes_weights",
    "diff_matrix",
]
```

and update the docstring's `kinetic` bullet to:

```
  - `kinetic`, `kinetic_sparse` -- assemble the FEM-DVR kinetic-energy matrix
    on a grid, dense or sparse (CSR). The dense one is the sparse one's
    differential oracle.
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest libs/qscat/tests/test_kinetic_sparse.py -q`
Expected: PASS (14 tests — 8 parametrized differential + 3 nnz-formula + 1 single-element + 2 structural)

- [ ] **Step 6: Confirm nothing regressed and type-check**

Run: `uv run pytest libs/qscat -q && uv run mypy libs/qscat && uv run ruff check .`
Expected: all pass, mypy clean, ruff clean.

- [ ] **Step 7: Commit**

```bash
git add libs/qscat/qscat/dvr libs/qscat/tests/test_kinetic_sparse.py
git commit -m "$(cat <<'EOF'
feat(dvr): kinetic_sparse -- sparse FEM-DVR kinetic assembly

Same mathematics as the dense kinetic(), emitted as COO triplets instead of
scatter-added into a dense array. coo_matrix sums duplicate (row, col)
entries on conversion to CSR, which reproduces the dense version's bridge
accumulation exactly -- no bridge special-casing.

The dense kinetic() is deliberately retained as the differential oracle.
Sparsity is checked against eMoScat's analytic nnz formula
q^2*tnel - 4q + 3 - tnel, verified independently against both production N2
grids (2995 and 6382 nonzeros).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: `qscat.dvr.TensorGrid`

**Files:**
- Create: `libs/qscat/qscat/dvr/tensor.py`
- Modify: `libs/qscat/qscat/dvr/__init__.py`
- Test: `libs/qscat/tests/test_tensor_grid.py`

**Interfaces:**
- Consumes: `FemDvrEcsGrid` (existing).
- Produces:
  - `TensorGrid(grids: Sequence[FemDvrEcsGrid])`
  - `.grids -> tuple[FemDvrEcsGrid, ...]`, `.ndim -> int`, `.shape -> tuple[int, ...]`, `.size -> int`
  - `.points() -> tuple[npt.NDArray[np.complex128], ...]` — D **broadcastable** arrays
  - `.real_mask() -> npt.NDArray[np.bool_]` — flat, length `size`
  - `.outer(vectors: Sequence[npt.ArrayLike]) -> npt.NDArray[np.complex128]` — flat, length `size`

**Background you need:** `points()` returns arrays shaped for broadcasting, not a full meshgrid — for D=2 that is `(n_0, 1)` and `(1, n_1)`. This keeps memory down and lets a potential be written naturally as `V(r, R) = v0(R) + v_int(r, R)`.

`real_mask()` is a **safety feature**. Under ECS, any driving term or channel projection must be confined to the unscaled region or the resulting matrix element is meaningless. `grid.real_points` holds the unscaled coordinate and `grid.R0` is the ECS pivot, so a point is unscaled iff `real_points <= R0`. The mask is the elementwise AND across all dimensions.

- [ ] **Step 1: Write the failing tests — `libs/qscat/tests/test_tensor_grid.py`**

```python
"""Tests for `qscat.dvr.TensorGrid`: geometry, broadcasting, the ECS real-region
mask, and separable-state construction.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec, TensorGrid


def _real_grid(q: int, n_el: int, length: float = 1.0) -> FemDvrEcsGrid:
    return FemDvrEcsGrid(
        GridSpec(quadrature=q, elements=[ElementSpec(length) for _ in range(n_el)])
    )


def _ecs_grid(q: int, n_real: int, n_cplx: int) -> FemDvrEcsGrid:
    els = [ElementSpec(1.0) for _ in range(n_real)]
    els += [ElementSpec(1.0, 35.0) for _ in range(n_cplx)]
    return FemDvrEcsGrid(GridSpec(quadrature=q, elements=els))


def test_shape_size_and_ndim() -> None:
    ga, gb, gc = _real_grid(6, 2), _real_grid(5, 3), _real_grid(4, 2)
    tg = TensorGrid([ga, gb, gc])
    assert tg.ndim == 3
    assert tg.shape == (ga.n, gb.n, gc.n)
    assert tg.size == ga.n * gb.n * gc.n
    assert tg.grids == (ga, gb, gc)


def test_points_are_broadcastable_not_meshgrid() -> None:
    ga, gb = _real_grid(6, 2), _real_grid(5, 3)
    tg = TensorGrid([ga, gb])
    pa, pb = tg.points()
    assert pa.shape == (ga.n, 1)
    assert pb.shape == (1, gb.n)
    # broadcasting them together reproduces the full grid
    assert np.broadcast_shapes(pa.shape, pb.shape) == (ga.n, gb.n)
    assert np.allclose(pa.ravel(), ga.points)
    assert np.allclose(pb.ravel(), gb.points)


def test_points_for_d1_is_plain_1d() -> None:
    g = _real_grid(6, 3)
    (p,) = TensorGrid([g]).points()
    assert p.shape == (g.n,)
    assert np.allclose(p, g.points)


def test_real_mask_is_and_across_dimensions() -> None:
    ga = _ecs_grid(6, 3, 2)
    gb = _ecs_grid(5, 2, 2)
    tg = TensorGrid([ga, gb])
    mask = tg.real_mask()
    assert mask.shape == (tg.size,)
    assert mask.dtype == np.bool_
    ma = ga.real_points <= ga.R0
    mb = gb.real_points <= gb.R0
    assert np.array_equal(mask, np.outer(ma, mb).ravel())
    # a genuine mixture -- otherwise the test proves nothing
    assert 0 < int(mask.sum()) < mask.size


def test_real_mask_all_true_when_no_ecs() -> None:
    tg = TensorGrid([_real_grid(6, 2), _real_grid(5, 3)])
    assert tg.real_mask().all()


def test_outer_builds_separable_state_in_c_order() -> None:
    ga, gb = _real_grid(6, 2), _real_grid(5, 3)
    tg = TensorGrid([ga, gb])
    rng = np.random.default_rng(0)
    a = rng.standard_normal(ga.n) + 1j * rng.standard_normal(ga.n)
    b = rng.standard_normal(gb.n) + 1j * rng.standard_normal(gb.n)
    psi = tg.outer([a, b])
    assert psi.shape == (tg.size,)
    assert np.allclose(psi.reshape(tg.shape), np.outer(a, b))


def test_outer_rejects_wrong_length() -> None:
    tg = TensorGrid([_real_grid(6, 2), _real_grid(5, 3)])
    with pytest.raises(ValueError, match="expected 2"):
        tg.outer([np.ones(3)])


def test_rejects_empty_grid_list() -> None:
    with pytest.raises(ValueError, match="at least one"):
        TensorGrid([])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest libs/qscat/tests/test_tensor_grid.py -q`
Expected: FAIL — `ImportError: cannot import name 'TensorGrid' from 'qscat.dvr'`

- [ ] **Step 3: Implement `libs/qscat/qscat/dvr/tensor.py`**

```python
"""Tensor products of FEM-DVR-ECS grids, and the N-dimensional Hamiltonian
assembled on them.

A separable-kinetic, diagonal-potential Hamiltonian on a tensor product of DVR
grids is a Kronecker sum plus a diagonal:

    H = sum_d  I x ... x T_d x ... x I  +  diag(V(x_0, ..., x_{D-1}))

Nothing about that is specific to two dimensions, or to molecular scattering.
This module is the dimension-general form; `qscat.linalg.kron_sum` does the
Kronecker algebra and knows nothing about grids.

Index convention: numpy-native C order, LAST axis fastest, so a state of shape
`tgrid.shape` ravels to the vector the Hamiltonian acts on. eMoScat uses the
opposite convention (first coordinate fastest); the two are physically
identical and differ only in basis ordering.

See `docs/physics/nd-tensor-hamiltonian.md`.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.linalg import kron_sum

from .grid import FemDvrEcsGrid
from .kinetic import kinetic_sparse

__all__ = ["TensorGrid", "kinetic_nd", "potential_nd", "hamiltonian_nd"]


class TensorGrid:
    """Tensor product of D FEM-DVR-ECS grids (C order: last axis fastest)."""

    def __init__(self, grids: Sequence[FemDvrEcsGrid]) -> None:
        tup = tuple(grids)
        if not tup:
            raise ValueError("TensorGrid requires at least one grid")
        self._grids = tup

    @property
    def grids(self) -> tuple[FemDvrEcsGrid, ...]:
        return self._grids

    @property
    def ndim(self) -> int:
        return len(self._grids)

    @property
    def shape(self) -> tuple[int, ...]:
        return tuple(g.n for g in self._grids)

    @property
    def size(self) -> int:
        return int(np.prod(self.shape))

    def _broadcast_shape(self, d: int) -> tuple[int, ...]:
        """Shape that puts axis `d`'s data on axis `d` and 1 everywhere else."""
        return (1,) * d + (-1,) + (1,) * (self.ndim - d - 1)

    def points(self) -> tuple[npt.NDArray[np.complex128], ...]:
        """The D complex (ECS-scaled) coordinate arrays, shaped to broadcast.

        For D = 2 these are `(n_0, 1)` and `(1, n_1)`, so a potential can be
        written naturally as `V(r, R)` without materializing a full meshgrid.
        For D = 1 the single array is plain 1-D of shape `(n_0,)`.
        """
        return tuple(
            np.asarray(g.points, dtype=np.complex128).reshape(self._broadcast_shape(d))
            for d, g in enumerate(self._grids)
        )

    def real_mask(self) -> npt.NDArray[np.bool_]:
        """Flat boolean mask, True where EVERY coordinate is in the unscaled region.

        Under exterior complex scaling a driving term or channel projection is
        only meaningful on the unscaled region, so anything of that kind must be
        masked with this before use. Making it a property of the grid is
        deliberate: the physics layer should not have to remember.
        """
        mask: npt.NDArray[np.bool_] | None = None
        for d, g in enumerate(self._grids):
            md = np.asarray(g.real_points <= g.R0, dtype=bool).reshape(
                self._broadcast_shape(d)
            )
            mask = md if mask is None else (mask & md)
        assert mask is not None  # ndim >= 1 guaranteed by __init__
        return np.asarray(np.broadcast_to(mask, self.shape).ravel(), dtype=bool)

    def outer(self, vectors: Sequence[npt.ArrayLike]) -> npt.NDArray[np.complex128]:
        """Separable state `⊗_d vectors[d]`, flattened to length `size`."""
        vecs = list(vectors)
        if len(vecs) != self.ndim:
            raise ValueError(f"expected {self.ndim} vectors, got {len(vecs)}")
        for d, v in enumerate(vecs):
            got = np.asarray(v).shape
            if got != (self.shape[d],):
                raise ValueError(
                    f"vector {d} has shape {got}, expected {(self.shape[d],)}"
                )
        out = np.asarray(vecs[0], dtype=np.complex128)
        for v in vecs[1:]:
            out = np.multiply.outer(out, np.asarray(v, dtype=np.complex128))
        return np.asarray(out.ravel(), dtype=np.complex128)


def kinetic_nd(tgrid: TensorGrid, masses: Sequence[float]) -> sp.csr_matrix:
    """`sum_d I x ... x T_d x ... x I`, with `T_d` built at mass `masses[d]`."""
    ms = list(masses)
    if len(ms) != tgrid.ndim:
        raise ValueError(f"expected {tgrid.ndim} masses, got {len(ms)}")
    return kron_sum(
        [kinetic_sparse(g, m) for g, m in zip(tgrid.grids, ms, strict=True)]
    )


def potential_nd(
    tgrid: TensorGrid, V: Callable[..., npt.ArrayLike]
) -> npt.NDArray[np.complex128]:
    """Evaluate `V` at the D-dimensional COMPLEX points, flattened.

    `V` is called as `V(x_0, ..., x_{D-1})` with the broadcastable arrays from
    `TensorGrid.points()`. It MUST NOT coerce its arguments to a real dtype:
    the points are complex on the ECS tail, and discarding the imaginary part
    silently destroys the analytic continuation the method depends on.
    """
    vals = np.asarray(V(*tgrid.points()), dtype=np.complex128)
    return np.asarray(
        np.broadcast_to(vals, tgrid.shape).ravel(), dtype=np.complex128
    )


def hamiltonian_nd(
    tgrid: TensorGrid, masses: Sequence[float], V: Callable[..., npt.ArrayLike]
) -> sp.csr_matrix:
    """`H = kinetic_nd(tgrid, masses) + diag(potential_nd(tgrid, V))` as CSR.

    Complex symmetric (`H = H^T`), NOT Hermitian, whenever any grid has an ECS
    tail. Use general algorithms only.
    """
    T = kinetic_nd(tgrid, masses)
    V_diag = potential_nd(tgrid, V)
    return sp.csr_matrix(T + sp.diags(V_diag, format="csr"))
```

- [ ] **Step 4: Export from `libs/qscat/qscat/dvr/__init__.py`**

Add the import (after the `.spec` import to keep alphabetical-ish order consistent with the file):

```python
from .tensor import TensorGrid, hamiltonian_nd, kinetic_nd, potential_nd
```

and add to `__all__`: `"TensorGrid"`, `"kinetic_nd"`, `"potential_nd"`, `"hamiltonian_nd"`.

Add to the module docstring's Public API list:

```
  - `TensorGrid` -- tensor product of D FEM-DVR-ECS grids (C order, last axis
    fastest), with the ECS real-region mask and separable-state construction.
  - `kinetic_nd`, `potential_nd`, `hamiltonian_nd` -- the N-dimensional
    Kronecker-sum Hamiltonian assembled on a `TensorGrid`, sparse (CSR).
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest libs/qscat/tests/test_tensor_grid.py -q`
Expected: PASS (8 tests)

- [ ] **Step 6: Type-check and lint**

Run: `uv run mypy libs/qscat && uv run ruff check .`
Expected: clean.

- [ ] **Step 7: Commit**

```bash
git add libs/qscat/qscat/dvr libs/qscat/tests/test_tensor_grid.py
git commit -m "$(cat <<'EOF'
feat(dvr): TensorGrid + N-dimensional Kronecker-sum Hamiltonian

TensorGrid composes D FEM-DVR-ECS grids in numpy-native C order (last axis
fastest), exposing broadcastable complex points, separable-state
construction, and the ECS real-region mask.

real_mask() is a safety feature rather than a convenience: under exterior
complex scaling, a driving term or channel projection confined to the
unscaled region is the difference between a meaningful matrix element and a
meaningless one, and the reference C++ gets away with a Hermitian dot
product only because it zeroes those functions on the tail.

kinetic_nd/potential_nd/hamiltonian_nd assemble H = sum_d I x .. T_d .. x I
+ diag(V) sparsely for any D. Analytic validation lands in the next task.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 5 (CRUX): Analytic validation at D = 1, 2, 3

**Files:**
- Test: `libs/qscat/tests/test_hamiltonian_nd.py`

**Interfaces:**
- Consumes: `TensorGrid`, `kinetic_nd`, `potential_nd`, `hamiltonian_nd` (Task 4); `hamiltonian` (existing dense, for V4).
- Produces: nothing new — this task is pure validation, and it is the reason the whole sub-project exists.

**Background you need:** This is where "works in arbitrary dimension" stops being a claim and becomes a measurement. Two analytic benchmarks, each run at **D = 1, 2 and 3**, with **unequal per-axis extents and masses** so that a transposed-axis bug cannot pass.

1. **Particle in a D-dimensional box.** With Dirichlet walls at both ends — exactly what `FemDvrEcsGrid` does by dropping both endpoints — the eigenvalues are `E = sum_d n_d^2 * pi^2 / (2 * m_d * L_d^2)`, `n_d = 1, 2, ...`. Exact, no potential involved: this tests `kinetic_nd` alone.
2. **D-dimensional harmonic oscillator.** `V = sum_d 0.5 * m_d * omega_d^2 * x_d^2` on a box wide enough that the walls do not matter; eigenvalues `E = sum_d omega_d * (n_d + 0.5)`. This tests `potential_nd` and `hamiltonian_nd` together.

Plus **V4**: at D = 1 the new code must reproduce the existing, already-validated 1-D stack to round-off, which turns every result from sub-projects #1–#4 into a regression test on this code.

A note on tolerances: the box benchmark is spectrally exact and should reach ~1e-10 or better. The oscillator is limited by the finite box and finite basis; the plan's tolerances are **starting points**. After the tests pass, tighten each one to just above the accuracy actually achieved and leave a comment recording the measured value — the repo already does this (see `libs/qscat/tests/test_crank_nicolson.py`'s tolerance comment). A loose tolerance that never fails is not a test.

- [ ] **Step 1: Write the failing tests — `libs/qscat/tests/test_hamiltonian_nd.py`**

```python
"""Analytic validation of the N-dimensional Hamiltonian (V3) and the D=1
regression against the existing 1-D stack (V4).

The generality is EXERCISED, not asserted: every benchmark runs at D = 1, 2
and 3, with unequal per-axis extents and masses so a transposed-axis bug
cannot pass.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from qscat.dvr import (
    ElementSpec,
    FemDvrEcsGrid,
    GridSpec,
    TensorGrid,
    hamiltonian,
    hamiltonian_nd,
    kinetic,
    kinetic_nd,
    potential_nd,
)

Floats = npt.NDArray[np.float64]


def _box_grid(length: float, n_el: int, q: int, x_min: float = 0.0) -> FemDvrEcsGrid:
    """All-real grid tiling [x_min, x_min + length] with `n_el` equal elements."""
    return FemDvrEcsGrid(
        GridSpec(
            quadrature=q,
            elements=[ElementSpec(length / n_el) for _ in range(n_el)],
            x_min=x_min,
        )
    )


def _lowest_dense(H: sp.csr_matrix, k: int) -> Floats:
    """Lowest `k` eigenvalues by real part, via a dense eigensolve."""
    vals = np.linalg.eigvals(H.toarray())
    return np.asarray(np.sort(vals.real)[:k], dtype=np.float64)


def _lowest_shift_invert(H: sp.csr_matrix, k: int, sigma: float) -> Floats:
    """Lowest `k` eigenvalues near `sigma`, via sparse shift-invert."""
    vals = spla.eigs(H.tocsc(), k=k, sigma=sigma, return_eigenvectors=False)
    return np.asarray(np.sort(vals.real)[:k], dtype=np.float64)


def _analytic_box_levels(
    lengths: tuple[float, ...], masses: tuple[float, ...], k: int
) -> Floats:
    """Lowest `k` levels of sum_d n_d^2 pi^2 / (2 m_d L_d^2)."""
    per_axis = [
        [n**2 * np.pi**2 / (2.0 * m * L**2) for n in range(1, k + 2)]
        for L, m in zip(lengths, masses, strict=True)
    ]
    sums = [sum(c) for c in itertools.product(*per_axis)]
    return np.asarray(np.sort(np.array(sums))[:k], dtype=np.float64)


def _analytic_oscillator_levels(omegas: tuple[float, ...], k: int) -> Floats:
    """Lowest `k` levels of sum_d omega_d (n_d + 1/2)."""
    per_axis = [[w * (n + 0.5) for n in range(k + 1)] for w in omegas]
    sums = [sum(c) for c in itertools.product(*per_axis)]
    return np.asarray(np.sort(np.array(sums))[:k], dtype=np.float64)


# --------------------------------------------------------------------------
# V3a: particle in a D-dimensional box  (kinetic_nd alone -- no potential)
# --------------------------------------------------------------------------

def test_box_d1() -> None:
    L, m = (1.0,), (1.0,)
    tg = TensorGrid([_box_grid(L[0], 4, 8)])
    got = _lowest_dense(kinetic_nd(tg, m), 4)
    want = _analytic_box_levels(L, m, 4)
    assert np.allclose(got, want, rtol=1e-10, atol=0)


def test_box_d2_unequal_extents_and_masses() -> None:
    L, m = (1.0, 1.3), (1.0, 2.0)
    tg = TensorGrid([_box_grid(L[0], 4, 8), _box_grid(L[1], 4, 8)])
    got = _lowest_dense(kinetic_nd(tg, m), 5)
    want = _analytic_box_levels(L, m, 5)
    assert np.allclose(got, want, rtol=1e-10, atol=0)


def test_box_d3_unequal_extents_and_masses() -> None:
    L, m = (1.0, 1.3, 0.9), (1.0, 2.0, 1.5)
    tg = TensorGrid([_box_grid(L[d], 2, 6) for d in range(3)])
    got = _lowest_dense(kinetic_nd(tg, m), 4)
    want = _analytic_box_levels(L, m, 4)
    assert np.allclose(got, want, rtol=1e-9, atol=0)


# --------------------------------------------------------------------------
# V3b: harmonic oscillator  (potential_nd + hamiltonian_nd)
# --------------------------------------------------------------------------

def _oscillator_V(
    masses: tuple[float, ...], omegas: tuple[float, ...]
) -> Callable[..., npt.ArrayLike]:
    """V = sum_d 0.5 m_d omega_d^2 x_d^2, for the broadcastable coords of points()."""

    def V(*coords: npt.NDArray[np.complex128]) -> npt.ArrayLike:
        total: npt.NDArray[np.complex128] = np.zeros((), dtype=np.complex128)
        for x, m, w in zip(coords, masses, omegas, strict=True):
            total = total + 0.5 * m * w**2 * x**2
        return total

    return V


def test_oscillator_d1() -> None:
    m, w = (1.0,), (1.0,)
    tg = TensorGrid([_box_grid(16.0, 8, 10, x_min=-8.0)])
    H = hamiltonian_nd(tg, m, _oscillator_V(m, w))
    got = _lowest_dense(H, 4)
    want = _analytic_oscillator_levels(w, 4)
    assert np.allclose(got, want, rtol=1e-6, atol=0)


def test_oscillator_d2_unequal_frequencies() -> None:
    m, w = (1.0, 1.0), (1.0, 1.7)
    tg = TensorGrid([_box_grid(14.0, 4, 8, x_min=-7.0) for _ in range(2)])
    H = hamiltonian_nd(tg, m, _oscillator_V(m, w))
    got = _lowest_dense(H, 4)
    want = _analytic_oscillator_levels(w, 4)
    assert np.allclose(got, want, rtol=1e-4, atol=0)


def test_oscillator_d3_unequal_frequencies() -> None:
    m, w = (1.0, 1.0, 1.0), (1.0, 1.3, 1.7)
    tg = TensorGrid([_box_grid(12.0, 2, 10, x_min=-6.0) for _ in range(3)])
    H = hamiltonian_nd(tg, m, _oscillator_V(m, w))
    want = _analytic_oscillator_levels(w, 3)
    got = _lowest_shift_invert(H, 3, sigma=float(want[0]) - 0.3)
    assert np.allclose(got, want, rtol=1e-3, atol=0)


# --------------------------------------------------------------------------
# V4: D = 1 reproduces the existing 1-D stack bit-for-bit
# --------------------------------------------------------------------------

def test_d1_reproduces_existing_1d_hamiltonian_with_ecs() -> None:
    """Makes every sub-project #1-#4 result a regression test on this code."""
    els = [ElementSpec(1.0), ElementSpec(1.0), ElementSpec(2.0, 35.0), ElementSpec(2.0, 35.0)]
    g = FemDvrEcsGrid(GridSpec(quadrature=9, elements=els, x_min=0.0))
    mass = 12766.36

    def V(z: npt.NDArray[np.complex128]) -> npt.ArrayLike:
        return 1.0 / (1.0 + z**2)

    H_dense = hamiltonian(g, V, mass)
    H_nd = hamiltonian_nd(TensorGrid([g]), [mass], V).toarray()

    scale = np.abs(H_dense).max()
    assert np.abs(H_nd - H_dense).max() <= 1e-13 * scale


def test_d1_kinetic_nd_matches_dense_kinetic() -> None:
    g = FemDvrEcsGrid(
        GridSpec(quadrature=7, elements=[ElementSpec(1.0) for _ in range(3)])
    )
    dense = kinetic(g, 2.5)
    got = kinetic_nd(TensorGrid([g]), [2.5]).toarray()
    assert np.abs(got - dense).max() <= 1e-13 * np.abs(dense).max()


# --------------------------------------------------------------------------
# Structural: ECS makes H complex symmetric, never Hermitian
# --------------------------------------------------------------------------

def test_hamiltonian_nd_is_complex_symmetric_under_ecs() -> None:
    els = [ElementSpec(1.0), ElementSpec(1.0), ElementSpec(2.0, 35.0)]
    g = FemDvrEcsGrid(GridSpec(quadrature=8, elements=els))
    tg = TensorGrid([g, g])
    H = hamiltonian_nd(tg, [1.0, 1.0], lambda a, b: 1.0 / (1.0 + a**2 + b**2))
    assert abs(H - H.T).max() < 1e-10
    assert abs(H - H.conj().T).max() > 1e-6


def test_potential_nd_preserves_complex_points_on_the_ecs_tail() -> None:
    """A potential that coerced to float would silently kill the ECS tail."""
    els = [ElementSpec(1.0), ElementSpec(1.0), ElementSpec(2.0, 35.0)]
    g = FemDvrEcsGrid(GridSpec(quadrature=8, elements=els))
    tg = TensorGrid([g])
    vals = potential_nd(tg, lambda z: z**2)
    assert np.abs(vals.imag).max() > 1e-3
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest libs/qscat/tests/test_hamiltonian_nd.py -q`
Expected: FAIL (collection error or assertion failures, depending on what Task 4 left incomplete). If Task 4 is complete these may already pass — that is fine and expected for a validation-only task; proceed to Step 3 and treat any failure as a real defect to fix.

- [ ] **Step 3: Make them pass**

If a benchmark fails, the defect is in Task 4's code, not the test. Debug in this order:

1. **Wrong index order** — symptom: D=2 with *equal* dimensions passes but unequal fails, or `test_kron_sum_acts_on_c_order_ravel` passes while a box test does not. Check `_broadcast_shape` and that `kron_sum` puts `left = prod(sizes[:d])` on the left.
2. **Masses applied to the wrong axis** — symptom: D=1 fine, D=2 levels wrong but plausible. Check the `zip(tgrid.grids, ms, strict=True)` pairing in `kinetic_nd`.
3. **Oscillator levels systematically high** — the box is too narrow or the basis too small. Widen the box or raise the element count; do NOT loosen the tolerance to hide it.
4. **`eigs` fails to converge in the D=3 oscillator** — move `sigma` closer to the target level, or raise `k`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest libs/qscat/tests/test_hamiltonian_nd.py -q -v`
Expected: PASS (10 tests)

- [ ] **Step 5: Tighten every tolerance to the measured accuracy**

Run with reported values and tighten each `rtol` to just above what is actually achieved, adding a comment recording the measured number. For example, if `test_oscillator_d2_unequal_frequencies` achieves 3e-6, set `rtol=1e-5` with a comment `# measured 3e-6 at this basis size`.

Re-run to confirm still green: `uv run pytest libs/qscat/tests/test_hamiltonian_nd.py -q`

- [ ] **Step 6: Full suite, types, lint, and the N₂ harness**

```bash
uv run pytest libs/qscat -q
uv run mypy libs/qscat
uv run ruff check .
uv run python -m validation.n2.experiment
```
Expected: all tests pass; mypy 0 errors; ruff clean; the harness still reports **19 PASS / 0 PENDING / 2 NOTE / 0 FAIL** and exits 0.

- [ ] **Step 7: Commit**

```bash
git add libs/qscat/tests/test_hamiltonian_nd.py
git commit -m "$(cat <<'EOF'
test(dvr): analytic validation of the N-dimensional Hamiltonian at D=1,2,3

Two analytic benchmarks, each exercised at D = 1, 2 and 3 with unequal
per-axis extents, masses and frequencies so a transposed-axis bug cannot
pass: the D-dimensional particle in a box (exact, tests kinetic_nd alone)
and the D-dimensional harmonic oscillator (tests potential_nd and
hamiltonian_nd together).

Also pins D=1 to the existing dense 1-D stack to round-off, which makes
every result already validated in sub-projects #1-#4 a regression test on
this code, and checks that ECS leaves H complex symmetric rather than
Hermitian and that potential_nd keeps the imaginary part of the tail.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Production-scale smoke test, docs, and `CLAUDE.md`

**Files:**
- Modify: `pyproject.toml` (register the `slow` marker)
- Test: `libs/qscat/tests/test_nd_scale.py`
- Create: `docs/physics/nd-tensor-hamiltonian.md`
- Modify: `CLAUDE.md`

**Interfaces:**
- Consumes: everything from Tasks 1–5.
- Produces: nothing consumed by later tasks in this plan (sub-project #6 consumes the library itself).

**Background you need:** V6 checks that the library survives contact with a production-scale problem — the real N₂ 2-D grid dimensions — **without factorizing**, which needs 13.6 GB and 128 s and has no place in a routine test run. Keep it grid-only with a generic analytic potential: `libs/qscat` must not import from `validation/` or `projects/`, and the N₂-specific assembly belongs to sub-project #6.

Measured reference values from the design-spec spike, for the real `N2-model.json` deck:

| quantity | value |
|---|---|
| electronic grid | `q=8`, 33 real + 15 ECS elements → `n = 335` |
| nuclear grid | `q=14`, 23 real + 10 ECS elements → `n = 428` |
| N | 143,380 |
| nnz | 3,276,450 (22.9/row) |
| `max\|H − Hᵀ\|` | 1.7e-13 |

- [ ] **Step 1: Register the `slow` marker in `pyproject.toml`**

In the `[tool.pytest.ini_options]` block, add:

```toml
markers = [
    "slow: heavier checks (production-scale assembly); deselect with -m 'not slow'",
]
```

- [ ] **Step 2: Write the test — `libs/qscat/tests/test_nd_scale.py`**

```python
"""V6: the library survives a production-scale assembly.

Grid dimensions are the real eMoScat N2 2-D deck
(`reference/eMoScat/input/experimental/N2-model.json`), but the potential here
is a generic analytic function: `libs/qscat` must not depend on `validation/`
or `projects/`, and the N2-specific assembly belongs to sub-project #6.

Deliberately does NOT factorize. A measured spike put that at 128 s and
13.6 GB peak RSS -- real, and no business in a routine test run.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pytest
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec, TensorGrid, hamiltonian_nd

MU = 12766.36


def _ecs_tail(base: float, n: int, *, skip: int = 2, alpha: float = 0.2) -> list[float]:
    """eMoScat `uniform_increment`/`exp`: `skip` elements at `base`, then growing."""
    return [
        base if i < skip else base * float(np.exp(alpha * (i - skip + 1)))
        for i in range(n)
    ]


def _build(
    segments: list[tuple[float, float]],
    order: int,
    tail_base: float,
    n_tail: int,
    angle: float = 35.0,
) -> FemDvrEcsGrid:
    els: list[ElementSpec] = []
    start = 0.0
    for end, length in segments:
        k = round((end - start) / length)
        els += [ElementSpec((end - start) / k) for _ in range(k)]
        start = end
    els += [ElementSpec(h, angle) for h in _ecs_tail(tail_base, n_tail)]
    return FemDvrEcsGrid(GridSpec(quadrature=order, elements=els, x_min=0.0))


@pytest.mark.slow
def test_production_scale_2d_assembly() -> None:
    g_el = _build(
        [(1.0, 0.2), (5.0, 1.0), (7.0, 2.0), (10.0, 3.0), (98.0, 4.0)], 8, 4.0, 15
    )
    g_nu = _build([(1.5, 0.5), (3.0, 0.15), (4.0, 0.5), (12.0, 1.0)], 14, 1.0, 10)

    assert (g_el.n, g_nu.n) == (335, 428)

    tg = TensorGrid([g_el, g_nu])
    assert tg.size == 143_380

    def V(
        r: npt.NDArray[np.complex128], R: npt.NDArray[np.complex128]
    ) -> npt.ArrayLike:
        return 1.0 / (1.0 + r**2) + 1.0 / (1.0 + R**2)

    H = hamiltonian_nd(tg, [1.0, MU], V)
    assert H.shape == (143_380, 143_380)
    # matches eMoScat's own nnz formula, independently re-derived
    assert H.nnz == 3_276_450
    assert abs(H - H.T).max() < 1e-10          # complex symmetric, never Hermitian
```

- [ ] **Step 3: Run it**

Run: `uv run pytest libs/qscat/tests/test_nd_scale.py -q -m slow`
Expected: PASS (1 test), in a few seconds — assembly alone is ~0.1 s.

If `nnz` differs, the ECS tail element count is off; check `_ecs_tail` against the counts in the table above (`n = 335` and `n = 428` are the gate — they must match before `nnz` can).

- [ ] **Step 4: Confirm the default run deselects it**

Run: `uv run pytest libs/qscat -q -m "not slow"`
Expected: all fast tests pass, the scale test deselected, no `PytestUnknownMarkWarning`.

- [ ] **Step 5: Write `docs/physics/nd-tensor-hamiltonian.md`**

Cover, in prose with formulas:
- The construction: `H = sum_d I x ... x T_d x ... x I + diag(V)`, and the two conditions that make it valid — **separable kinetic energy** and a **diagonal (DVR) potential**.
- Why it is dimension-general, and that `kron_sum` accepts any square sparse operators so non-FEM-DVR dimensions compose.
- The index convention (C order, last axis fastest) and the explicit divergence from eMoScat's first-fastest convention.
- ECS consequences: `H = Hᵀ ≠ H†`; the c-product; why `real_mask` exists and what goes wrong without it.
- The validation actually performed: box and oscillator at D = 1, 2, 3; the D=1 regression against the 1-D stack; the sparsity formula `q²·tnel − 4q + 3 − tnel`.
- Measured cost at production scale (the table above, plus the 128 s / ×93 / 13.6 GB factorization figures and the 440 ms back-substitution), and the consequence: all channels at one energy share a factorization, so cost scales with the number of *energies*, not channels.
- Link to `docs/superpowers/specs/2026-07-22-nd-sparse-hamiltonian-design.md` and `.superpowers/sdd/n2-2d-exact-extraction.md`.

- [ ] **Step 6: Update `CLAUDE.md`**

In the repo-map block under `libs/`, add a `qscat.linalg` entry and extend the `qscat.dvr` entry:

```
            - qscat.linalg: dimension-general sparse linear algebra --
              `kron_sum` (Kronecker sum over arbitrary D), `SparseLU` (cached
              factorization with fill-in/memory diagnostics), and `c_product`
              (the bilinear, non-conjugated ECS inner product) -- see
              docs/physics/nd-tensor-hamiltonian.md.
```

and append to the existing `qscat.dvr` bullet:

```
              Also `kinetic_sparse` and the N-dimensional tensor layer
              (`TensorGrid`, `kinetic_nd`, `potential_nd`, `hamiltonian_nd`):
              H = sum_d I x .. T_d .. x I + diag(V) for any D, sparse (CSR),
              validated on analytic box/oscillator benchmarks at D = 1, 2, 3.
```

- [ ] **Step 7: Full verification**

```bash
uv run pytest libs/qscat -q
uv run pytest -q -m "not slow"
uv run mypy libs/qscat
uv run ruff check .
uv run python -m validation.n2.experiment
docker/build.sh test
```
Expected: all green; mypy 0; ruff clean; harness **19 PASS / 0 PENDING / 2 NOTE / 0 FAIL**, exit 0; docker test target passes.

- [ ] **Step 8: Commit**

```bash
git add pyproject.toml libs/qscat/tests/test_nd_scale.py docs/physics/nd-tensor-hamiltonian.md CLAUDE.md
git commit -m "$(cat <<'EOF'
test+docs: production-scale assembly check and the N-dim Hamiltonian note

A slow-marked smoke test assembles at the real eMoScat N2 2-D grid
dimensions (N = 143,380) and pins nnz = 3,276,450 -- eMoScat's own formula,
independently re-derived -- plus complex symmetry. It deliberately does not
factorize: measured at 128 s and 13.6 GB, that has no place in a routine
run. The potential is generic, keeping libs/qscat free of any dependency on
validation/ or projects/.

Documents the construction, the C-order index convention and its deliberate
divergence from eMoScat, the ECS consequences, and the measured costs.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Final verification

- [ ] `uv run pytest libs/qscat -q` — all pass, including the slow scale test.
- [ ] `uv run pytest -q` — the whole repo suite passes; nothing in `projects/` or `validation/` regressed.
- [ ] `uv run mypy libs/qscat` — 0 errors.
- [ ] `uv run ruff check .` — clean.
- [ ] `uv run python -m validation.n2.experiment` — **19 PASS / 0 PENDING / 2 NOTE / 0 FAIL**, exit 0.
- [ ] `docker/build.sh test` — passes.
- [ ] V1–V6 from the design spec are each covered by a named test.
- [ ] No N₂ physics, no T-matrix, no cross sections anywhere in `libs/qscat` — that is sub-project #6.
- [ ] `CLAUDE.md` and `docs/physics/nd-tensor-hamiltonian.md` describe what was actually built and measured.
