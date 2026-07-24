# N₂ Exact 2-D Time-Dependent VE Cross-Section Implementation Plan (sub-project #7)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Solve the exact 2-D electron–N₂ problem in the time domain — launch a Gaussian electron wavepacket, propagate it under `H_2D` by sparse Crank-Nicolson, and read the cross section off the correlation functions — with #6's exact σ_TI as the differential oracle. Numeric time-series outputs are the primary deliverable; figures are the visual layer.

**Architecture:** `projects/n2_2d_td_cross_section/`, built on #6's exact `H_2D`/grids/channel functions and #4's TD pattern. The generic sparse Crank-Nicolson propagator promotes to `qscat.evolution`. The energy transform is Tannor–Weeks (genuine wavepacket + η-deconvolution), not the #4 doorway form — the exact model has no Γ, so there is no doorway.

**Tech Stack:** Python 3.12, numpy, `scipy.sparse`, `qscat.linalg.SparseLU`, `qscat.evolution`, `scipy.special.spherical_jn`.

**Design spec:** `docs/superpowers/specs/2026-07-24-n2-2d-td-cross-section-design.md`
**Source archaeology:** `.superpowers/sdd/n2-2d-exact-extraction.md` (sections 4.1, 5.1–5.3)

## Global Constraints

- Python `>=3.12`. Everything through `uv` — never bare `python`/`pip`/conda. `PYTHONPATH=.` for `uv run python` scripts.
- **Atomic units throughout.** `μ = 12766.36`, electron mass `1`, `l = 2`. Cross sections in bohr².
- Package-absolute imports only. `projects/` may import `projects.*`; **`projects/` must NOT import `validation/`** (reverse only). Harness/oracle glue lives in `validation/`.
- `uv run mypy libs/qscat` stays at **0 errors**. `uv run ruff check .` stays clean (line length 100; rules `E, F, I, UP, B, NPY`; per-file `E741` ignore for a bare `l` is the established precedent).
- The existing N₂ harness must not regress: **23 PASS / 0 PENDING / 4 NOTE / 0 FAIL**, exit 0. Group F *adds* rows.
- **`H_2D` is complex symmetric, never Hermitian.** Use `qscat.linalg.c_product`, never `np.vdot`. Under the ECS contour `‖Ψ(t)‖` decays (outgoing flux absorbed) — this is correct, not a bug.
- **The σ prefactor is the S-matrix form `σ = π|S − δ|²/(2E)`**, which equals #6's `4π³|T|²/(2E)` via `S = 1 − 2πiT`. Do NOT write `4π³|S − δ|²`.
- Index order is numpy-native C-order (last axis fastest). Axis 0 = electronic `r`, axis 1 = nuclear `R`.
- **#6's exact σ_TI is the oracle:** `S_TD(E) → S_TI(E)`. Convergence to it (not a loose bound) is the crux gate.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`

---

## File Structure

| File | Responsibility |
|---|---|
| `libs/qscat/qscat/evolution/crank_nicolson.py` | **Modify.** Add `make_sparse_cn_stepper` |
| `libs/qscat/qscat/evolution/__init__.py` | **Modify.** Export it |
| `libs/qscat/tests/test_crank_nicolson.py` | **Modify.** Add sparse-vs-dense tests |
| `projects/n2_2d_td_cross_section/__init__.py` | **Create.** Empty package marker |
| `projects/n2_2d_td_cross_section/wavepacket.py` | **Create.** Gaussian wavepacket + 2-D initial state |
| `projects/n2_2d_td_cross_section/td_propagation.py` | **Create.** Propagate + two-cadence sampling (the data engine) |
| `projects/n2_2d_td_cross_section/td_cross_section.py` | **Create.** Tannor–Weeks transform → σ(E) (crux) |
| `projects/n2_2d_td_cross_section/observation.py` | **Create.** Numeric outputs + PNG snapshots |
| `projects/n2_2d_td_cross_section/convergence.py` | **Create.** Box/dt/T sizing; the TD working grid |
| `projects/n2_2d_td_cross_section/test_*.py` | **Create.** Per-task tests |
| `validation/n2/td_exact2d.py` | **Create.** Group F: σ_TD at anchors vs #6 σ_TI |
| `validation/n2/experiment.py` | **Modify.** Add group F |
| `docs/physics/n2-2d-td-cross-section.md` | **Create.** Method, numeric outputs, convergence, σ(E)-vs-σ_TI |
| `CLAUDE.md` | **Modify.** Add the sub-project + `qscat.evolution.make_sparse_cn_stepper` |

---

### Task 1: Sparse Crank-Nicolson propagator → `qscat.evolution`

**Files:**
- Modify: `libs/qscat/qscat/evolution/crank_nicolson.py`, `libs/qscat/qscat/evolution/__init__.py`
- Test: `libs/qscat/tests/test_crank_nicolson.py`

**Interfaces:**
- Consumes: `qscat.linalg.SparseLU`.
- Produces: `make_sparse_cn_stepper(H: sparse.spmatrix, dt: float) -> Callable[[NDArray], NDArray]`.

**Background you need.** The dense `make_cn_stepper` (same file) builds `A = I + iH dt/2`, `B = I − iH dt/2`, LU-factors `A` once, and each `stepper(psi)` returns `lu_solve(lu, B @ psi)`. The sparse sibling is identical mathematics with `A` factored by `qscat.linalg.SparseLU` (one back-substitution per step). `H` here is the ~10⁴–10⁵-dimension sparse `H_2D`; dense factorization is infeasible, which is the whole reason this exists. The dense stepper is the differential oracle.

- [ ] **Step 1: Write the failing tests — append to `libs/qscat/tests/test_crank_nicolson.py`**

```python
import scipy.sparse as sp
from qscat.evolution import make_sparse_cn_stepper


def test_sparse_cn_matches_dense_cn() -> None:
    """The sparse stepper must equal the dense one to round-off on the same H."""
    rng = np.random.default_rng(7)
    n = 40
    a = sp.random(n, n, density=0.1, format="csr", random_state=rng, dtype=float)
    b = sp.random(n, n, density=0.1, format="csr", random_state=rng, dtype=float)
    H_sp = (a + 1j * b).tocsr()
    H_dense = H_sp.toarray()
    dt = 0.05
    psi0 = (rng.standard_normal(n) + 1j * rng.standard_normal(n)).astype(np.complex128)
    dense_step = make_cn_stepper(H_dense, dt)
    sparse_step = make_sparse_cn_stepper(H_sp, dt)
    assert np.linalg.norm(sparse_step(psi0) - dense_step(psi0)) < 1e-11


def test_sparse_cn_matches_exact_exp_for_hermitian() -> None:
    rng = np.random.default_rng(8)
    n = 6
    a = sp.random(n, n, density=0.3, format="csr", random_state=rng, dtype=float)
    b = sp.random(n, n, density=0.3, format="csr", random_state=rng, dtype=float)
    m = (a + 1j * b).toarray()
    H = m + m.conj().T  # Hermitian
    dt = 1e-3
    step = make_sparse_cn_stepper(sp.csr_matrix(H), dt)
    psi0 = (rng.standard_normal(n) + 1j * rng.standard_normal(n)).astype(np.complex128)
    exact = expm(-1j * H * dt) @ psi0
    assert np.linalg.norm(step(psi0) - exact) < 5e-8


def test_sparse_cn_unitary_for_hermitian() -> None:
    rng = np.random.default_rng(9)
    n = 8
    a = sp.random(n, n, density=0.4, format="csr", random_state=rng, dtype=float)
    H_dense = a.toarray()
    H = sp.csr_matrix(H_dense + H_dense.T)  # real symmetric => Hermitian
    step = make_sparse_cn_stepper(H, 0.1)
    psi = (rng.standard_normal(n) + 1j * rng.standard_normal(n)).astype(np.complex128)
    n0 = np.vdot(psi, psi).real
    for _ in range(50):
        psi = step(psi)
    assert abs(np.vdot(psi, psi).real - n0) < 1e-10


def test_sparse_cn_decays_for_absorbing_H() -> None:
    """Negative-imaginary diagonal (ECS/optical) => norm strictly decreases."""
    H = sp.diags([1.0 - 0.1j, 2.0 - 0.2j, 3.0 - 0.05j], format="csr")
    step = make_sparse_cn_stepper(H, 0.05)
    psi = np.ones(3, dtype=np.complex128)
    n0 = np.vdot(psi, psi).real
    for _ in range(100):
        psi = step(psi)
    assert np.vdot(psi, psi).real < n0
```

- [ ] **Step 2: Run → fail.** `uv run pytest libs/qscat/tests/test_crank_nicolson.py -q` → `ImportError: cannot import name 'make_sparse_cn_stepper'`.

- [ ] **Step 3: Implement in `libs/qscat/qscat/evolution/crank_nicolson.py`.** Add the imports and function:

```python
import scipy.sparse as sp

from qscat.linalg import SparseLU
```

Extend `__all__` to `["make_cn_stepper", "make_sparse_cn_stepper"]`, then append:

```python
def make_sparse_cn_stepper(
    H: sp.spmatrix, dt: float
) -> Callable[[npt.NDArray[np.complexfloating[Any, Any]]], npt.NDArray[np.complex128]]:
    """Sparse Crank-Nicolson stepper -- the sparse sibling of `make_cn_stepper`.

    Same Cayley form `(I + i H dt/2) psi_{n+1} = (I - i H dt/2) psi_n`, but
    `A = I + i H dt/2` is factored once with `qscat.linalg.SparseLU` and each
    step is a single sparse back-substitution. For the ~1e4-1e5-dimension
    sparse Hamiltonians this targets, dense factorization is infeasible; the
    dense `make_cn_stepper` is retained as this function's differential oracle.

    `H` must be square and sparse. Complex symmetric (ECS) `H` is fine -- no
    Hermiticity is assumed.
    """
    n = H.shape[0]
    ident = sp.identity(n, format="csr", dtype=np.complex128)
    A = (ident + 0.5j * dt * H).tocsc()
    B = (ident - 0.5j * dt * H).tocsr()
    lu = SparseLU(A)

    def stepper(
        psi: npt.NDArray[np.complexfloating[Any, Any]],
    ) -> npt.NDArray[np.complex128]:
        result: npt.NDArray[np.complex128] = lu.solve(B @ psi)
        return result

    return stepper
```

Add `make_sparse_cn_stepper` to `libs/qscat/qscat/evolution/__init__.py`'s import, `__all__`, and the docstring's Public API list (one bullet: "the sparse sibling for large sparse H, factoring once with `SparseLU`").

- [ ] **Step 4: Run → pass** (4 new tests). Then `uv run mypy libs/qscat` → 0; `uv run ruff check .` → clean.

- [ ] **Step 5: Commit.**

```bash
git add libs/qscat/qscat/evolution libs/qscat/tests/test_crank_nicolson.py
git commit -m "$(cat <<'EOF'
feat(evolution): make_sparse_cn_stepper -- sparse Crank-Nicolson propagator

The sparse sibling of make_cn_stepper: same Cayley form, but A = I + iH dt/2
is factored once with qscat.linalg.SparseLU and each step is one sparse
back-substitution. Dense factorization is infeasible at the 1e4-1e5
dimensions the 2-D time-dependent solver needs; the dense stepper is kept as
the differential oracle and the sparse one is tested against it to round-off,
plus exp-match / unitarity / absorbing-decay.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Gaussian wavepacket + 2-D initial state

**Files:**
- Create: `projects/n2_2d_td_cross_section/__init__.py` (empty), `projects/n2_2d_td_cross_section/wavepacket.py`
- Test: `projects/n2_2d_td_cross_section/test_wavepacket.py`

**Interfaces:**
- Consumes: `qscat.dvr.TensorGrid`, `qscat.linalg.c_product`, `projects.n2_ti_cross_section.vibrational.vibrational_states`.
- Produces:
  - `gaussian_coeffs(grid, r0, p0, sigma) -> NDArray[complex128]` — masked DVR coefficients of `g(r)` on the electronic grid.
  - `initial_state(tgrid, chi_v, *, r0, p0, sigma) -> NDArray[complex128]` — `g(r) ⊗ χ_v(R)`, flat, masked.

**Background you need.** The incident wavepacket (`eMoScat input.cpp:240`) is
`g(r) = (2π σ²)^{-1/4} exp(−(r−r₀)²/(4σ²)) exp(i p₀ r)`, L²-normalized in the continuum. On the FEM-DVR grid a *function* becomes coefficients via `c_j = g(r_j)·sqrt(w_j)` (the same convention as #6's `channel_vector`), using `TensorGrid.sqrt_weights()[0]` for the electronic axis. `g` is evaluated on **real** electronic points only and the state is masked to the unscaled region (`real_mask`) — a wavepacket launched into the ECS tail is meaningless. `p₀ < 0` is inward (toward the molecule at small r).

**Hermitian norm vs c-product — get this distinction right, it matters throughout #7.** The initial state is renormalized to unit **Hermitian L2 norm** (`‖ψ‖₂ = sqrt(Σ|ψ_j|²)`, the true remaining-probability norm), NOT the c-product self-pairing. For a wavepacket with a momentum phase, `c_product(g,g) = Σ g_j² = Σ envelope_j² e^{2i p₀ r_j} w_j` is a small oscillatory *complex* number, nowhere near 1 — it is not a norm. The Hermitian L2 norm is what decays monotonically under the absorbing ECS contour (provable: for CN, `A†A − B†B = −2 Im(H) dt ≥ 0` when `Im(H) ≤ 0`, so `‖A⁻¹B‖₂ ≤ 1`), and is the physically meaningful depletion diagnostic. The **c-product** (no conjugate) is used ONLY for the correlation functions and the S-matrix, where the ECS analytic continuation requires it.

- [ ] **Step 1: Write the failing tests — `projects/n2_2d_td_cross_section/test_wavepacket.py`**

```python
"""Incident Gaussian electron wavepacket and the 2-D initial state g(r) x chi_v(R)."""

from __future__ import annotations

import numpy as np
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_2d_td_cross_section.wavepacket import gaussian_coeffs, initial_state
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states

# A box big enough to hold a wavepacket launched near r0=20 (test-scale, not production).
TG = TensorGrid(
    [
        n2_electronic_grid(r_max=40.0, order=8, n_complex=6),
        n2_nuclear_grid(quadrature=10, r_max=22.0, n_complex=5),
    ]
)
EPS, CHI = vibrational_states(TG.grids[1], MU, 4)


def test_gaussian_is_localized_near_r0() -> None:
    g = n2_electronic_grid(r_max=40.0, order=8, n_complex=6)
    r0 = 20.0
    coeffs = gaussian_coeffs(g, r0=r0, p0=-0.35, sigma=4.0)
    # density peak (|coeff|^2, which is the DVR density weight) sits near r0
    r_real = g.real_points
    peak_r = r_real[np.argmax(np.abs(coeffs[: r_real.size]) ** 2)]
    assert abs(peak_r - r0) < 3.0


def test_gaussian_carries_inward_momentum() -> None:
    """<p> < 0 for p0 < 0 -- the wavepacket moves toward the molecule."""
    g = n2_electronic_grid(r_max=40.0, order=8, n_complex=6)
    coeffs = gaussian_coeffs(g, r0=20.0, p0=-0.35, sigma=4.0)
    # phase increments negatively in r: consecutive real-region coeffs rotate clockwise
    real = g.real_points <= g.R0
    c = coeffs[: g.real_points.size][real[: g.real_points.size]]
    phases = np.unwrap(np.angle(c[np.abs(c) > 1e-6]))
    assert phases[-1] - phases[0] < 0.0


def test_gaussian_masked_to_unscaled_region() -> None:
    g = n2_electronic_grid(r_max=40.0, order=8, n_complex=6)
    coeffs = gaussian_coeffs(g, r0=20.0, p0=-0.35, sigma=4.0)
    tail = g.real_points > g.R0
    assert np.all(coeffs[tail] == 0.0)


def test_initial_state_is_unit_hermitian_norm_and_separable() -> None:
    psi = initial_state(TG, CHI[0], r0=20.0, p0=-0.35, sigma=4.0)
    assert psi.shape == (TG.size,)
    # Hermitian L2 norm == 1 (the physical probability norm), NOT the c-product
    assert abs(float(np.linalg.norm(psi)) - 1.0) < 1e-10
    # separable: reshape factorizes to outer(g_coeff, chi) up to scale
    block = psi.reshape(TG.shape)
    u, s, vh = np.linalg.svd(block)
    assert s[0] / s.sum() > 0.999  # essentially rank-1


def test_initial_state_masked() -> None:
    psi = initial_state(TG, CHI[0], r0=20.0, p0=-0.35, sigma=4.0)
    assert np.all(psi[~TG.real_mask()] == 0.0)
```

- [ ] **Step 2: Run → fail** (`ModuleNotFoundError: ... wavepacket`).

- [ ] **Step 3: Implement `projects/n2_2d_td_cross_section/wavepacket.py`.**

```python
"""Incident Gaussian electron wavepacket and the 2-D initial state.

`g(r) = (2 pi sigma^2)^{-1/4} exp(-(r-r0)^2/(4 sigma^2)) exp(i p0 r)`
(eMoScat `input.cpp:240`), converted to FEM-DVR coefficients on the unscaled
electronic region (`c_j = g(r_j) sqrt(w_j)`, same convention as #6's
`channel_vector`). `p0 < 0` launches the packet inward, toward the molecule;
the ECS tail (not evaluated here) absorbs whatever leaves during propagation.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.dvr import FemDvrEcsGrid, TensorGrid
from qscat.linalg import c_product

__all__ = ["gaussian_coeffs", "initial_state"]


def gaussian_coeffs(
    grid: FemDvrEcsGrid, *, r0: float, p0: float, sigma: float
) -> npt.NDArray[np.complex128]:
    """DVR coefficients of `g(r)` on `grid`, zero on the ECS tail."""
    r = grid.real_points
    envelope = (2.0 * np.pi * sigma**2) ** -0.25 * np.exp(-((r - r0) ** 2) / (4.0 * sigma**2))
    g_vals = envelope * np.exp(1j * p0 * r)
    sqrt_w = np.sqrt(np.asarray(grid.weights, dtype=np.complex128))
    coeffs = (g_vals * sqrt_w).astype(np.complex128)
    coeffs[r > grid.R0] = 0.0  # unscaled region only
    return coeffs


def initial_state(
    tgrid: TensorGrid,
    chi_v: npt.NDArray[np.complex128],
    *,
    r0: float,
    p0: float,
    sigma: float,
) -> npt.NDArray[np.complex128]:
    """`Psi(0) = g(r) chi_v(R)`, flat, masked, renormalized to unit c-product norm."""
    g_coeff = gaussian_coeffs(tgrid.grids[0], r0=r0, p0=p0, sigma=sigma)
    chi = np.asarray(chi_v, dtype=np.complex128)
    chi = chi / np.sqrt(c_product(chi, chi))
    psi = tgrid.outer([g_coeff, chi])
    psi[~tgrid.real_mask()] = 0.0
    # Hermitian L2 (probability) norm -- NOT the c-product, which for a phased
    # wavepacket is a small oscillatory complex number, not a norm. The c-product
    # is used only for correlations/S-matrix (module `td_propagation`/`correlation`).
    norm = float(np.linalg.norm(psi))
    return np.asarray(psi / norm, dtype=np.complex128)
```

- [ ] **Step 4: Run → pass** (5 tests). `uv run ruff check .` → clean.

- [ ] **Step 5: Commit.**

```bash
git add projects/n2_2d_td_cross_section
git commit -m "$(cat <<'EOF'
feat(n2-2d-td): incident Gaussian electron wavepacket + 2-D initial state

g(r) = (2 pi sigma^2)^-1/4 exp(-(r-r0)^2/4sigma^2) exp(i p0 r), converted to
FEM-DVR coefficients on the unscaled electronic region (c_j = g(r_j) sqrt(w_j),
#6's convention) and tensored with the neutral vibrational state chi_v(R).
p0 < 0 launches inward; masked to the unscaled region and renormalized to unit
c-product norm so norm decay under the ECS contour is measured cleanly.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Propagation engine + two-cadence sampling

**Files:**
- Create: `projects/n2_2d_td_cross_section/td_propagation.py`
- Test: `projects/n2_2d_td_cross_section/test_td_propagation.py`

**Interfaces:**
- Consumes: Task 1 (`make_sparse_cn_stepper`), Task 2 (`initial_state`), #6 (`build_h2d`, `channel_vector`... but see below), `qscat.linalg.c_product`.
- Produces:
  - `PropagationResult` dataclass: `t` `(n_t,)`, `c` `(n_t, n_channels)` complex, `norm` `(n_t,)`, `snapshots: list[Snapshot]`.
  - `Snapshot` dataclass: `time: float`, `rho_R: NDArray[float64]`, `rho_r: NDArray[float64]`, `psi: NDArray[complex128] | None`.
  - `propagate(tgrid, psi0, out_channels, *, dt, n_steps, sample_period, snapshot_times, keep_psi_at) -> PropagationResult`.

**Background you need — this is the numeric-output data engine.** Two cadences: `c_{v'}(t_n)` is recorded at EVERY step (it is the transform's raw material — fine), while densities/norm snapshots are recorded on a coarse schedule (`sample_period` steps, or explicit `snapshot_times`). `out_channels` is a list of *test-function* states `Φ_{v'}` (built in Task 4, but for this task's tests pass simple separable states). Nuclear density `ρ(R,t) = Σ_r |Ψ(r,R,t)|²` over electronic real points, reshaping the flat state to `tgrid.shape` (axis 0 = r, axis 1 = R) and summing axis 0 restricted to the electronic unscaled region; electronic density symmetrically. The full `Ψ` is kept only at times in `keep_psi_at` (memory). **Norm is the Hermitian L2 norm `np.linalg.norm(psi)`** (real, ≥0, monotone non-increasing under the absorbing contour — see Task 2's note), the physical depletion diagnostic; the c-product is used only for the correlations `c_{v'}(t) = c_product(Φ_{v'}, Ψ(t))`.

- [ ] **Step 1: Write the failing tests — `projects/n2_2d_td_cross_section/test_td_propagation.py`**

```python
"""Two-cadence propagation: fine c(t), coarse density/norm snapshots."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.dvr import TensorGrid
from qscat.linalg import c_product

from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_2d_td_cross_section.td_propagation import propagate
from projects.n2_2d_td_cross_section.wavepacket import initial_state
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states

TG = TensorGrid(
    [
        n2_electronic_grid(r_max=40.0, order=8, n_complex=6),
        n2_nuclear_grid(quadrature=10, r_max=22.0, n_complex=5),
    ]
)
EPS, CHI = vibrational_states(TG.grids[1], MU, 4)
PSI0 = initial_state(TG, CHI[0], r0=20.0, p0=-0.35, sigma=4.0)
# simple separable test-function channels chi_v'(R) x (uniform electronic weight)
OUT = [TG.outer([np.ones(TG.shape[0]), CHI[v]]) for v in range(3)]


def test_correlation_recorded_every_step_at_t0_zero() -> None:
    res = propagate(TG, PSI0, OUT, dt=1.0, n_steps=20, sample_period=5)
    assert res.t.shape == (21,)
    assert res.c.shape == (21, 3)
    # c_{v'}(0) = c_product(Phi_v', Psi0)
    for k in range(3):
        assert res.c[0, k] == c_product(OUT[k], PSI0)


def test_norm_decays_under_absorbing_contour() -> None:
    res = propagate(TG, PSI0, OUT, dt=1.0, n_steps=60, sample_period=10)
    assert res.norm[0] == pytest.approx(1.0, abs=1e-9)
    assert res.norm[-1] < res.norm[0]  # ECS absorbs outgoing flux
    assert np.all(np.diff(res.norm) <= 1e-12)  # monotone non-increasing


def test_snapshots_on_coarse_cadence_and_densities_nonneg() -> None:
    res = propagate(TG, PSI0, OUT, dt=1.0, n_steps=20, sample_period=5)
    assert [s.time for s in res.snapshots] == [0.0, 5.0, 10.0, 15.0, 20.0]
    for s in res.snapshots:
        assert s.rho_R.shape == (TG.shape[1],)
        assert s.rho_r.shape == (TG.shape[0],)
        assert np.all(s.rho_R >= 0.0) and np.all(s.rho_r >= 0.0)
        assert s.psi is None  # not requested


def test_full_psi_kept_only_at_requested_times() -> None:
    res = propagate(
        TG, PSI0, OUT, dt=1.0, n_steps=20, sample_period=5, keep_psi_at=[0.0, 10.0]
    )
    kept = {s.time: s.psi for s in res.snapshots}
    assert kept[0.0] is not None and kept[10.0] is not None
    assert kept[5.0] is None
    assert kept[0.0].shape == (TG.size,)


def test_explicit_snapshot_times() -> None:
    res = propagate(TG, PSI0, OUT, dt=0.5, n_steps=40, snapshot_times=[0.0, 5.0, 20.0])
    assert [s.time for s in res.snapshots] == [0.0, 5.0, 20.0]
```

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement `projects/n2_2d_td_cross_section/td_propagation.py`.**

```python
"""Time-domain propagation engine: propagate Psi(0) under H_2D and sample.

Two cadences (the numeric-output design):
  * `c_{v'}(t_n)` -- recorded at EVERY step; the raw material of the
    Tannor-Weeks energy transform and the literal "formation from the
    correlation functions".
  * density/norm snapshots -- recorded on a COARSE schedule (`sample_period`
    steps, or explicit `snapshot_times`), so the wavefunction is observed at
    static points without storing every step.

`H_2D` is time-independent, so the sparse Crank-Nicolson factorization is built
once and reused; under the ECS contour `||Psi||` decays as outgoing flux is
absorbed (the resonance depletes). Norm here is the magnitude of the ECS
c-product self-pairing, not the Hermitian norm.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.dvr import TensorGrid
from qscat.evolution import make_sparse_cn_stepper
from qscat.linalg import c_product

from projects.n2_2d_cross_section.hamiltonian2d import build_h2d

__all__ = ["Snapshot", "PropagationResult", "propagate"]


@dataclass(frozen=True)
class Snapshot:
    time: float
    rho_R: npt.NDArray[np.float64]  # nuclear density, sum_r |Psi|^2 (unscaled)
    rho_r: npt.NDArray[np.float64]  # electronic density, sum_R |Psi|^2 (unscaled)
    psi: npt.NDArray[np.complex128] | None  # full state, only if requested


@dataclass(frozen=True)
class PropagationResult:
    t: npt.NDArray[np.float64]  # (n_t,)  sample times n*dt
    c: npt.NDArray[np.complex128]  # (n_t, n_channels)  c_{v'}(t_n)
    norm: npt.NDArray[np.float64]  # (n_t,)  |c_product(Psi,Psi)|^{1/2}
    snapshots: list[Snapshot]


def _densities(
    tgrid: TensorGrid, psi: npt.NDArray[np.complex128]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    block = psi.reshape(tgrid.shape)
    dens = np.abs(block) ** 2
    r_real = tgrid.grids[0].real_points <= tgrid.grids[0].R0
    R_real = tgrid.grids[1].real_points <= tgrid.grids[1].R0
    rho_R = dens[r_real, :].sum(axis=0)
    rho_r = dens[:, R_real].sum(axis=1)
    return rho_R.astype(np.float64), rho_r.astype(np.float64)


def propagate(
    tgrid: TensorGrid,
    psi0: npt.NDArray[np.complex128],
    out_channels: list[npt.NDArray[np.complex128]],
    *,
    dt: float,
    n_steps: int,
    sample_period: int = 0,
    snapshot_times: list[float] | None = None,
    keep_psi_at: list[float] | None = None,
) -> PropagationResult:
    """Propagate and sample. See module docstring for the two cadences."""
    H = build_h2d(tgrid)
    step = make_sparse_cn_stepper(H, dt)

    n_t = n_steps + 1
    t = np.arange(n_t, dtype=np.float64) * dt
    n_ch = len(out_channels)
    c = np.empty((n_t, n_ch), dtype=np.complex128)
    norm = np.empty(n_t, dtype=np.float64)

    if snapshot_times is not None:
        snap_set = {round(x / dt) for x in snapshot_times}
    elif sample_period > 0:
        snap_set = set(range(0, n_t, sample_period)) | {n_t - 1}
    else:
        snap_set = {0, n_t - 1}
    keep_set = {round(x / dt) for x in (keep_psi_at or [])}

    psi = psi0.astype(np.complex128).copy()
    snapshots: list[Snapshot] = []
    for n in range(n_t):
        for k in range(n_ch):
            c[n, k] = c_product(out_channels[k], psi)  # correlation: c-product
        norm[n] = float(np.linalg.norm(psi))  # Hermitian L2: physical, monotone
        if n in snap_set:
            rho_R, rho_r = _densities(tgrid, psi)
            snapshots.append(
                Snapshot(float(t[n]), rho_R, rho_r, psi.copy() if n in keep_set else None)
            )
        if n < n_steps:
            psi = step(psi)

    return PropagationResult(t=t, c=c, norm=norm, snapshots=snapshots)
```

- [ ] **Step 4: Run → pass** (5 tests). `uv run ruff check .` → clean.

- [ ] **Step 5: Commit** ("feat(n2-2d-td): two-cadence propagation engine — fine c(t), coarse density/norm snapshots").

---

### Task 4 (CRUX): Tannor–Weeks energy transform → σ(E), gated by TD ≈ TI

**Files:**
- Create: `projects/n2_2d_td_cross_section/correlation.py`, `projects/n2_2d_td_cross_section/td_cross_section.py`
- Test: `projects/n2_2d_td_cross_section/test_td_cross_section.py`

**Interfaces:**
- Consumes: Task 2/3; #6 (`channel_vector` for `Φ_{v'}`, `riccati_bessel_en`, `MU`, `ELL`); #6's oracle `projects.n2_2d_cross_section.cross_section_2d.ve_cross_section_2d`; `vibrational_states`.
- Produces:
  - `outgoing_channel(tgrid, chi_v, *, r0_out, p0_out, sigma_out) -> NDArray` — the test function `Φ_{v'} = g_out(r) ⊗ χ_{v'}(R)`, masked. It is **energy-independent** (the k'-dependence lives only in `eta_outgoing`), so it takes no `k`.
  - `eta_incident(grid, k, l, *, r0, p0, sigma) -> complex`; `eta_outgoing(grid, kp, l, *, r0_out, p0_out, sigma_out) -> complex`.
  - `td_ve_cross_section_2d(tgrid, eps, chi, v_init, vprimes, E, *, dt, n_steps, wp_in, wp_out) -> NDArray` — σ via propagate + correlate + Tannor–Weeks transform.

**Background you need — the transform, and why it is the crux.** From the correlation functions `c_{v'}(t)` (Task 3),

```
S_{v→v'}(E) = [2π · conj(eta_out_{v'}(E)) · eta_in_v(E)]^{-1}  * Σ_n w_n e^{i E_tot t_n} c_{v'}(t_n) dt
sigma       = π |S − δ|^2 / (2E)        (elastic subtracts 1; = 4π³|T|²/(2E))
```

(`eMoScat TestFunction2d.cpp:298-307`), `E_tot = E + ε_{v_init}`, `w_n` composite (Simpson, trapezoid fallback) — reuse #4's `_quadrature_weights` pattern. The η factors deconvolve the wavepackets' own spectral content: `eta_in(E) = c_product(g_in, F_{E,l})` (incident Gaussian coeffs · energy-normalized regular free coeffs, on the electronic real region) and `eta_out(E) = c_product(g_out, F^out_{E',l})`. The outgoing test function `Φ_{v'} = g_out(r) ⊗ χ_{v'}(R)` (masked, coefficient form).

**This is the exact differential oracle, and the gate.** `S_TD(E) = S_TI(E)`, so **σ_TD must converge to #6's `ve_cross_section_2d`** at converged `dt`/`n_steps`. If it does not, debug in this order (do NOT loosen the tolerance):
1. **The η definitions** — `eta_in`/`eta_out` must be the SAME energy-normalized `F` (`riccati_bessel_en`) that #6 projects onto, coefficient-converted with the SAME `sqrt(w_r)`.
2. **The 2π and the conjugation** — `conj(eta_out)` not `eta_out`; the `1/(2π …)` prefactor.
3. **`E_tot` in the phase** — `E + ε_{v_init}`, not `E`.
4. **c-product everywhere** (no `vdot`); the mask on `g_out`.
5. **Finite-T truncation** — `‖Ψ(T)‖` must have decayed; lengthen `T` before suspecting a bug.
6. **The prefactor** — `π|S−δ|²/(2E)`, and elastic subtracts δ.

If it will not converge to #6 on any affordable `dt`/`T`, report **BLOCKED** with σ_TD vs σ_TI and the norm-decay profile — exactly as #4's crux task did.

- [ ] **Step 1: Write the crux tests — `projects/n2_2d_td_cross_section/test_td_cross_section.py`.** Build the shared 2-D setup ONCE at module scope (a modest but wavepacket-capable grid). Tests:
  - **V2a (TD≈TI):** at an anchor (e.g. E=0.1, v'=1) and one more (E=0.15, v'=1), σ_TD (converged dt/n_steps) agrees with #6's `ve_cross_section_2d` to `rtol ≤ 0.10`. σ_TD real and ≥0.
  - **V2b (free-particle sanity):** with the interaction scaled to zero via #6's `lam_scale=0` path is not available here (TD has no lam_scale); instead assert closed channels give 0 and that a v' with `E_tot−ε_{v'}≤0` is 0.
  - **V4 (convergence/depletion):** σ_TD at (E=0.1, v'=1) with `n_steps` and `2·n_steps` agree to `rtol ≤ 0.05`; and `‖Ψ(T)‖ < 0.2·‖Ψ(0)‖`.

  (Write concrete asserts with the real values you measure; set the tolerances at the measured level per the "tighten to measured" rule.)

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement `correlation.py`** (`outgoing_channel`, `eta_incident`, `eta_outgoing` — coefficient overlaps via `c_product` on the electronic real region, reusing `riccati_bessel_en` and `gaussian_coeffs`) **and `td_cross_section.py`** (`td_ve_cross_section_2d`: build `Φ_{v'}` for each v', `propagate` from `initial_state`, then the Tannor–Weeks transform per the formula above; accept scalar or array `E`, matching #6's return shape convention). Reuse #4's `_quadrature_weights`.

- [ ] **Step 4: Run.** Tune `dt` (~0.5–2 a.u.), `n_steps` (T long enough for norm decay), and the wavepacket parameters (`wp_in`/`wp_out`) until V2a/V4 pass. Follow the debug order above on failure. **Commit** once green, recording the converged (dt, n_steps, T, norm-decay, σ_TD-vs-σ_TI ratios) in the message.

---

### Task 5: The full σ(E) curve + convergence / box sizing

**Files:**
- Create: `projects/n2_2d_td_cross_section/convergence.py`
- Test: `projects/n2_2d_td_cross_section/test_td_convergence.py`

**Interfaces:**
- Consumes: Task 4.
- Produces: `TD_WORKING_GRID: dict`, `td_working_tgrid() -> TensorGrid`, `sigma_curve(E_grid, vprimes, *, dt, n_steps, ...) -> NDArray` (one propagation → whole curve), `usable_window(E_grid, ...) -> tuple[float, float]`.

**Background you need.** The wavepacket needs room: a packet launched at `r₀` needs `r_max` comfortably beyond `r₀ + several σ`, so the TD box is larger than #6's TI working grid (`r_max=16`). Size `r_max`, `r₀`, `p₀`, `σ`, `dt`, `T` so that (a) σ_TD matches σ_TI at the anchors, and (b) the **usable energy window** — where `|η_in(E)|` is large enough that the deconvolution is not noise-dominated — covers the resonance region. **The whole σ(E) curve comes from ONE propagation** (the transform is cheap once `c(t)` is stored) — this is TD's signature advantage and the "cross section formation".

- [ ] **Step 1:** Write `convergence.py` sweeping box/`dt`/`T`/wavepacket params, timing each, and reporting σ_TD-vs-σ_TI and the usable window. Save the table to `.superpowers/sdd/task-5-td-convergence-table.md`.
- [ ] **Step 2:** Run the study (background + poll — it is minutes-scale per config). Record the table.
- [ ] **Step 3:** Choose `TD_WORKING_GRID` — smallest box/T where σ_TD matches σ_TI to ~1% across the usable window, each value commented with the measured number. Expose `td_working_tgrid()`.
- [ ] **Step 4:** Write `test_td_convergence.py` (`@pytest.mark.slow`): σ_TD(E) matches σ_TI(E) across the usable window within the measured tolerance; the usable window is reported and non-empty; refining `T`/`dt` past the working values changes the curve < ~1%.
- [ ] **Step 5:** Commit with the convergence table and the chosen grid in the message.

---

### Task 6: Numeric-output layer + PNG snapshots

**Files:**
- Create: `projects/n2_2d_td_cross_section/observation.py`
- Test: `projects/n2_2d_td_cross_section/test_observation.py`

**Interfaces:**
- Consumes: Tasks 3–5.
- Produces:
  - `save_numeric_outputs(result, sigma_E, E_grid, path) -> None` — writes the sampled arrays (`t`, `c(t)`, `norm(t)`, per-snapshot `ρ(R,t)`/`ρ(r,t)`, `E_grid`, `σ(E)`) to a documented `.npz`.
  - `plot_snapshots(result, path) -> None`, `plot_correlation(result, path) -> None`, `plot_sigma_vs_ti(E_grid, sigma_td, sigma_ti, usable, path) -> None`.

**Background you need.** The numeric arrays are the primary deliverable and are what the tests assert on; the PNGs are a thin `matplotlib.use("Agg")` layer over them. `plot_sigma_vs_ti` overlays σ_TD(E) on #6's exact σ_TI(E) and shades the usable window (outside it is η-deconvolution noise — plot it faded or clipped, honestly labeled). Snapshots show `ρ(R,t)` and `ρ(r,t)` at the sampled times (the incoming packet → transient anion at the molecule → decay). Commit the three PNGs into `docs/physics/figures/`.

- [ ] **Step 1:** Tests — the `.npz` round-trips (load → arrays match `result`); re-transforming the saved `c(t)` reproduces the saved `σ(E)` (self-consistency, V5); densities real/≥0; the three plot functions produce non-empty PNG files. No pass/fail on the physics shape.
- [ ] **Step 2:** Implement `observation.py`.
- [ ] **Step 3:** Generate the real figures at `TD_WORKING_GRID` and commit them to `docs/physics/figures/` (`n2-2d-td-snapshots.png`, `n2-2d-td-correlation.png`, `n2-2d-td-sigma.png`). Record the numeric `.npz` output location.
- [ ] **Step 4:** Commit.

---

### Task 7: Harness group F + docs

**Files:**
- Create: `validation/n2/td_exact2d.py`
- Modify: `validation/n2/experiment.py`, `CLAUDE.md`
- Create: `docs/physics/n2-2d-td-cross-section.md`

**Interfaces:**
- Consumes: Task 4/5; `validation.n2.exact2d` (#6's σ_TI at the anchors); `validation.n2.reference`.
- Produces: `TdExact2dResult` dataclass + `compute_td_exact2d_results()` (`lru_cache`d).

**Background you need.** Mirror `validation/n2/exact2d.py` (#6's group-E module) and the C5/D/E harness pattern. Group F computes σ_TD at the gated anchors and compares to #6's σ_TI (tight rtol — the exact differential oracle) and to Houfek (factor-3, for context). Because TD needs the big box, **decide from measurement** whether group F runs the full `TD_WORKING_GRID` or a documented reduced config to stay Docker-affordable (time it; the reduced config's looser tolerance must be stated in the row detail and the docs). Guard in try/except → labeled FAIL. Preserve exit 0 and the existing 23/0/4/0.

- [ ] **Step 1:** Implement `validation/n2/td_exact2d.py` (reuse #6's σ_TI via `exact2d`; group anchors by energy; one propagation per energy).
- [ ] **Step 2:** Wire group F into `experiment.py` (PASS/FAIL on σ_TD-vs-σ_TI rtol at gated anchors; detail prints TD/TI and TD/Houfek). Time it; pick full-vs-reduced grid; document.
- [ ] **Step 3:** Write `docs/physics/n2-2d-td-cross-section.md`: the method (wavepacket, sparse CN, Tannor–Weeks transform, S_TD=S_TI); the numeric-output design (two cadences, the `.npz`); the validation ladder (sparse-CN-vs-dense, TD≈TI, convergence, norm decay); the convergence/box study and `TD_WORKING_GRID`; the full σ(E) curve overlaid on σ_TI with the honest usable window; the three figures; and the framing (this is the time-domain route to the same exact cross section, and the sparse LU is the eventual optimize-in-Rust target, NOT done here). Cross-reference `docs/physics/n2-2d-cross-section.md`.
- [ ] **Step 4:** Update `CLAUDE.md` — add `projects/n2_2d_td_cross_section`, `qscat.evolution.make_sparse_cn_stepper`, group F, and the new doc.
- [ ] **Step 5:** Full verification:
```
uv run pytest projects/n2_2d_td_cross_section libs/qscat validation/n2 -q
uv run pytest -q -m "not slow"
uv run mypy libs/qscat
uv run ruff check .
uv run python -m validation.n2.experiment      # group F present; exit 0; no regression of 23/0/4/0
docker/build.sh test
```
- [ ] **Step 6:** Commit.

---

## Final verification

- [ ] `make_sparse_cn_stepper` promoted to `qscat.evolution`, matches the dense stepper to round-off, mypy-clean.
- [ ] σ_TD converges to #6's exact σ_TI at the gated anchors (the exact differential oracle); convergence in dt/T; norm decays.
- [ ] The full σ(E) curve comes from one propagation and overlays σ_TI across the honestly-reported usable window (η-deconvolution noise outside it not passed off as signal).
- [ ] The numeric outputs (`c(t)`, `ρ(R,t)`, `ρ(r,t)`, `‖Ψ(t)‖`, `σ(E)`) are saved as inspectable arrays and are self-consistent (re-transforming saved `c(t)` reproduces saved `σ(E)`); the three figures committed.
- [ ] Harness group F wired and guarded; existing 23 PASS / 0 PENDING / 4 NOTE / 0 FAIL not regressed; docker green.
- [ ] No model parameter tuned to match reference data; nothing beyond the sparse CN stepper promoted to `qscat`; no `projects/` → `validation/` import.
- [ ] The sparse-LU optimization is documented as the deferred next lifecycle stage, not attempted here.
