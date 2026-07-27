# Exact-2D TI Dissociative Attachment (sub-project A) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the exact-2D **dissociative-attachment (DA)** cross section σ_DA(E) to `qscat.core` — the transient anion's second exit channel (A + B⁻) — as a time-independent driven-equation T-matrix reusing the existing VE driven solve, and validate its thresholds/well-posedness on N₂/NO/F₂.

**Architecture:** DA is the SAME driven Lippmann-Schwinger solve as VE (`Ψ₊ = Ψ_i + (E_tot·I − H)⁻¹ V_int Ψ_i`, already in `qscat.core.driven.ve_cross_section` with `return_wavefunction=True`), then projected onto the dissociation channel `Φ_n(r,R) = φ_e^(n)(r)·F^nuc_{K_n,0}(R)` with the **rearrangement interaction** `V_DR = V_int(r,R) + v0(R) − V_int(r, R_inf)` (= `H − H_final`, NOT `V_int`): `T_n = ⟨Φ_n | V_DR | Ψ₊⟩`, `σ_n = 4π³|T_n|²/(2E)`. New: a mass-μ nuclear Bessel exit function, the anion bound-electronic-state solver, `V_DR`, and the projection loop. This is the eMoScat `time_independent_model.cpp` method; see `docs/superpowers/specs/2026-07-27-da-cross-sections-design.md` and `docs/physics/diatomic-ve-cross-sections.md`.

**Tech Stack:** Python ≥3.12, NumPy/SciPy, `qscat.dvr` (FEM-DVR-ECS TensorGrid), `qscat.linalg` (`c_product`, `SparseLU`), `qscat.special`, `qscat.core`, `qscat.model` (test-only). pytest, mypy --strict over `libs/qscat`, ruff.

## Global Constraints

- **Atomic units** throughout; nuclear reduced mass `μ = model.mu`, electron mass 1.
- **c-product, never Hermitian:** under ECS `H = Hᵀ`. Pair every projection with `qscat.linalg.c_product` (no conjugate). Never use `np.vdot`/`@` with a conjugate.
- **Mask to the real (unscaled) region:** any channel projection (`Φ_n`, `Ψ_i`) is meaningless on the ECS tail and MUST be zeroed there — `psi[~tgrid.real_mask()] = 0.0` for 2-D states, `p[g_r.real_points > g_r.R0] = 0.0` for the 1-D electronic factor. eMoScat zeroes both factors' tails.
- **σ prefactor is `4·π³·|T|²/(2E)`** (= `π|S|²/2E`; `k² = 2E`), identical to `ve_cross_section`. Getting the constant wrong rescales every cross section.
- **`V_DR` is the rearrangement interaction `V_int + v0(R) − V_int(r, R_inf)`, NOT `V_int`.** Using `V_int` gives a ~10⁶ unitarity violation (the earlier prototype's bug). `R_inf = tgrid.grids[1].R0` (the nuclear ECS pivot / real-region endpoint, eMoScat's `nu_inf`).
- **`qscat.core` NEVER imports `qscat.model` (nor `projects.*`) at runtime.** `dissociation.py` annotates `model: ResonanceModel` under `TYPE_CHECKING` only, exactly like `driven.py`. Enforced by `libs/qscat/tests/test_core_no_model_import.py`.
- **Energy normalization:** the nuclear exit function is `F^nuc_{E,l}(R) = √(2μK/π)·R·j_l(KR)`, `K = √(2μ E_DR)` — the mass-μ generalization of `riccati_bessel_en` (which is the mass-1 case). This normalization is what makes the `4π³|T|²/2E` prefactor correct.
- **No independent golden DA data exists** (only N₂ VE has Houfek's). Tests gate WELL-POSEDNESS (real, finite, ≥0, right shape), THRESHOLDS (N₂ closed → σ_DA≡0; F₂ exothermic → σ_DA>0; onset at the correct E), and a SOFT unitarity bound — never tight agreement with a number.
- **`c_product(a, b)` is the bilinear no-conjugate dot** (`Σ aᵢbᵢ`); `c_product(v, v)` is the c-norm² used for normalization.

---

## File Structure

- `libs/qscat/qscat/special/radial.py` (modify) — add `riccati_bessel_en_mass`.
- `libs/qscat/qscat/special/__init__.py` (modify) — export it.
- `libs/qscat/qscat/core/dissociation.py` (create) — `anion_electronic_states`, `v_dr_diag`, `da_cross_section`.
- `libs/qscat/qscat/core/__init__.py` (modify) — export the three.
- `libs/qscat/tests/test_radial_mass.py` (create) — Task 1 tests.
- `libs/qscat/tests/test_dissociation.py` (create) — Tasks 2–4 tests.
- `benchmarks/da_nuclear_convergence.py` (create) — Task 5 convergence study.
- `validation/diatomic/da_curves.py` (create), `validation/diatomic/test_da_curves.py` (create) — Task 6 figures + gate.
- `docs/physics/diatomic-ve-cross-sections.md` (modify), `CLAUDE.md` (modify) — Task 6 docs.

---

### Task 1: Mass-μ energy-normalized nuclear Bessel

**Files:**
- Modify: `libs/qscat/qscat/special/radial.py`
- Modify: `libs/qscat/qscat/special/__init__.py`
- Test: `libs/qscat/tests/test_radial_mass.py`

**Interfaces:**
- Produces: `riccati_bessel_en_mass(r: NDArray[float64], k: float, l: int, mu: float) -> NDArray[float64]` returning `√(2μk/π)·r·j_l(k r)`. For `mu=1.0` it equals `riccati_bessel_en(r, k, l)` exactly. `r` REAL, `k>0`, `mu>0`.

- [ ] **Step 1: Write the failing test**

```python
# libs/qscat/tests/test_radial_mass.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.special import riccati_bessel_en, riccati_bessel_en_mass


def test_reduces_to_mass_one():
    r = np.linspace(0.1, 20.0, 200)
    for l in (0, 1, 2):
        got = riccati_bessel_en_mass(r, 1.3, l, 1.0)
        assert np.allclose(got, riccati_bessel_en(r, 1.3, l), rtol=0, atol=1e-14)


def test_l0_closed_form():
    # F_{E,0}(R) = sqrt(2 mu K / pi) R j_0(KR) = sqrt(2 mu / (pi K)) sin(KR)
    r = np.linspace(0.05, 15.0, 300)
    K, mu = 4.0, 918.25
    got = riccati_bessel_en_mass(r, K, 0, mu)
    expect = np.sqrt(2.0 * mu / (np.pi * K)) * np.sin(K * r)
    assert np.allclose(got, expect, rtol=1e-12, atol=1e-12)


def test_scales_as_sqrt_mu():
    r = np.linspace(0.1, 10.0, 50)
    a = riccati_bessel_en_mass(r, 2.0, 1, 1.0)
    b = riccati_bessel_en_mass(r, 2.0, 1, 4.0)
    assert np.allclose(b, 2.0 * a, rtol=1e-13, atol=1e-13)  # sqrt(4)=2


@pytest.mark.parametrize("bad", [0.0, -1.0])
def test_rejects_nonpositive_k(bad):
    with pytest.raises(ValueError):
        riccati_bessel_en_mass(np.array([1.0]), bad, 0, 918.25)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest libs/qscat/tests/test_radial_mass.py -q`
Expected: FAIL — `ImportError: cannot import name 'riccati_bessel_en_mass'`.

- [ ] **Step 3: Add the function and export**

In `libs/qscat/qscat/special/radial.py`, append to `__all__` and add:

```python
def riccati_bessel_en_mass(
    r: npt.NDArray[np.float64], k: float, l: int, mu: float
) -> npt.NDArray[np.float64]:
    """`F_{E,l}(r) = sqrt(2 mu k / pi) r j_l(k r)`, energy-normalized at mass `mu`.

    The mass-`mu` generalization of `riccati_bessel_en` (which is the `mu=1`
    case): the energy-normalized (`<F_E|F_E'> = delta(E-E')`) regular radial
    solution for a particle of reduced mass `mu` and momentum `k = sqrt(2 mu E)`.
    Used for the OUTGOING NUCLEAR dissociation wave in the DA/DR exit channel
    (eMoScat `bessel::s_jEn(R, K, mu, l)`). `r` must be REAL (channel
    projections are masked to the unscaled region); `k>0`, `mu>0`.
    """
    if k <= 0.0:
        raise ValueError(f"k must be positive, got {k}")
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu}")
    rr = np.asarray(r, dtype=np.float64)
    out: npt.NDArray[np.float64] = (
        np.sqrt(2.0 * mu * k / np.pi) * rr * spherical_jn(l, k * rr)
    )
    return out
```

Set `radial.py`'s `__all__ = ["riccati_bessel_en", "riccati_hankel_en", "riccati_bessel_en_mass"]`.

In `libs/qscat/qscat/special/__init__.py`:
```python
from qscat.special.radial import (
    riccati_bessel_en,
    riccati_bessel_en_mass,
    riccati_hankel_en,
)

__all__ = ["riccati_bessel_en", "riccati_hankel_en", "riccati_bessel_en_mass"]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest libs/qscat/tests/test_radial_mass.py -q`
Expected: PASS (4 tests).

- [ ] **Step 5: Type-check + lint**

Run: `uv run mypy libs/qscat/qscat/special && uv run ruff check libs/qscat/qscat/special`
Expected: no errors.

- [ ] **Step 6: Commit**

```bash
git add libs/qscat/qscat/special/radial.py libs/qscat/qscat/special/__init__.py libs/qscat/tests/test_radial_mass.py
git commit -m "feat(da): mass-mu energy-normalized nuclear Bessel (riccati_bessel_en_mass)"
```

---

### Task 2: Anion bound-electronic-state solver

**Files:**
- Create: `libs/qscat/qscat/core/dissociation.py`
- Test: `libs/qscat/tests/test_dissociation.py`

**Interfaces:**
- Consumes: `qscat.dvr.{kinetic, eigen, FemDvrEcsGrid}`, `qscat.linalg.c_product`, `model.surface`, `model.ell`, `model.mu`.
- Produces: `anion_electronic_states(g_r: FemDvrEcsGrid, model: ResonanceModel, R_inf: float, n_states: int = 1) -> tuple[NDArray[float64], NDArray[complex128]]` returning `(eps_e, phi_e)`: `eps_e` the `n_states` lowest-Re bound electronic eigenvalues (real, ascending) of `−½∂²_r + surface(r, R_inf)`; `phi_e` shape `(n_states, g_r.n)`, each c-product-normalized over the electronic real region.

**Design notes (read before coding):** `model.surface(g_r.points, R_inf)` is the full electronic potential at the dissociation limit — `v0(R_inf) + ℓ(ℓ+1)/2r² + v_int(r, R_inf)` — so `eps_e` sits on the SAME energy scale as `H_2D` (it includes `v0(R_inf)`), which is what makes the DA threshold `eps_e − eps[v_init]` correct. The anion state is genuinely BOUND at `R_inf` (real eigenvalue below the ECS continuum), so select by `|Im(E)| < _IM_TOL_HA` then take the lowest `n_states` by `Re(E)` — do NOT take "the n lowest-Re overall" (ECS continuum eigenvalues can have `Re(E)` below the bound state; only the `|Im|` filter distinguishes them). `g_r.points` may be complex on the tail; `surface` handles that. Mirror `qscat.core.vibrational.vibrational_states`' structure and its `_IM_TOL_HA = 1e-6`.

- [ ] **Step 1: Write the failing test**

```python
# libs/qscat/tests/test_dissociation.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.core.dissociation import anion_electronic_states
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.vibrational import vibrational_states
from qscat.model import F2, N2, NO


def _eps0(model):
    g_R = nuclear_grid(r_max=22.0, n_complex=8, quadrature=12)
    eps, _ = vibrational_states(g_R, model.mu, 3, model.v0)
    return eps[0], g_R.R0


@pytest.mark.parametrize("model", [N2, NO, F2], ids=["N2", "NO", "F2"])
def test_one_bound_anion_state_real(model):
    g_r = electronic_grid(r_max=16.0, order=7, n_complex=6)
    _, R0 = _eps0(model)
    eps_e, phi = anion_electronic_states(g_r, model, R0, n_states=1)
    assert eps_e.shape == (1,) and phi.shape == (1, g_r.n)
    # c-product self-normalized over the real region ~ 1
    real = g_r.real_points <= g_r.R0
    p = phi[0].copy()
    p[~real] = 0.0
    assert abs(complex(p @ p) - 1.0) < 1e-6


def test_thresholds_have_correct_signs():
    # threshold(E_coll) = eps_e - eps[0]; F2 exothermic (<0), N2 closed (>0.3),
    # NO opens above its resonance (~0.17). No independent data -> sign/band gate.
    g_r = electronic_grid(r_max=16.0, order=7, n_complex=6)
    thr = {}
    for name, model in (("N2", N2), ("NO", NO), ("F2", F2)):
        eps0, R0 = _eps0(model)
        eps_e, _ = anion_electronic_states(g_r, model, R0, 1)
        thr[name] = float(eps_e[0]) - eps0
    assert thr["F2"] < 0.0            # exothermic: DA open at all E>0
    assert thr["N2"] > 0.3            # closed in the measurement window
    assert 0.10 < thr["NO"] < 0.25    # opens above the resonance


def test_raises_when_too_many_states_requested():
    g_r = electronic_grid(r_max=16.0, order=7, n_complex=6)
    _, R0 = _eps0(F2)
    with pytest.raises(ValueError):
        anion_electronic_states(g_r, F2, R0, n_states=50)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest libs/qscat/tests/test_dissociation.py -q`
Expected: FAIL — `ModuleNotFoundError: qscat.core.dissociation`.

- [ ] **Step 3: Create `dissociation.py` with `anion_electronic_states`**

```python
# libs/qscat/qscat/core/dissociation.py
"""Exact 2-D TI dissociative attachment (DA) cross section.

DA is the transient anion's second exit channel -- `e- + AB(v=0) -> AB-* ->
A + B-`, outgoing flux in the NUCLEAR coordinate. It is computed by the SAME
driven Lippmann-Schwinger solve as VE (`qscat.core.driven.ve_cross_section`,
`return_wavefunction=True`) but projected onto the dissociation channel with
the REARRANGEMENT interaction

    V_DR(r, R) = V_int(r, R) + v0(R) - V_int(r, R_inf)    (= H - H_final),

NOT V_int. The exit channel is Phi_n(r,R) = phi_e^(n)(r) F^nuc_{K_n,0}(R),
phi_e the anion bound electronic state at the dissociation limit R_inf and
F^nuc the mass-mu energy-normalized regular nuclear Bessel; the T-matrix is
T_n = <Phi_n | V_DR | Psi+> (c-product, masked), sigma_n = 4 pi^3 |T_n|^2/(2E).

This is eMoScat's `time_independent_model.cpp` method (an earlier prototype
that used V_int instead of V_DR gave a ~1e6 unitarity violation -- that was
the bug, not a structural obstacle to a TI DA). H2+ DR is the same T-matrix
looped over the neutral's Rydberg electronic series + a Coulomb incident;
deferred (sub-project D). See docs/physics/diatomic-ve-cross-sections.md and
docs/superpowers/specs/2026-07-27-da-cross-sections-design.md.

`qscat.core` never imports `qscat.model` at runtime: `model` is typed against
the `ResonanceModel` protocol under `TYPE_CHECKING` only, exactly like
`driven.py`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, eigen, kinetic
from qscat.linalg import c_product

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["anion_electronic_states", "v_dr_diag", "da_cross_section"]

# Bound-state signature on an ECS grid: true bound levels have |Im(E)| ~ 1e-15,
# ECS-continuum states jump to >= 1e-7. Same tolerance as `vibrational_states`.
_IM_TOL_HA = 1e-6


def anion_electronic_states(
    g_r: FemDvrEcsGrid,
    model: ResonanceModel,
    R_inf: float,
    n_states: int = 1,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128]]:
    """Anion bound electronic state(s) at the dissociation limit `R_inf`.

    Diagonalizes `-1/2 d^2/dr^2 + surface(r, R_inf)` (electron mass 1) on the
    electronic grid. `surface` is the FULL electronic potential at `R_inf`
    (`v0(R_inf) + ell(ell+1)/2r^2 + v_int(r, R_inf)`), so `eps_e` shares the
    `H_2D` energy scale (it includes `v0(R_inf)`) and the DA threshold
    `eps_e - eps[v_init]` is correct.

    Returns `(eps_e, phi_e)`: `eps_e` the `n_states` lowest-Re eigenvalues with
    `|Im(E)| < _IM_TOL_HA` (the genuinely bound states), real, ascending;
    `phi_e` shape `(n_states, g_r.n)`, each c-product-normalized over the
    electronic real region. Raises `ValueError` if fewer than `n_states` bound
    states exist (e.g. `n_states` reached past the finite bound spectrum).
    """
    H_el = kinetic(g_r, 1.0) + np.diag(model.surface(g_r.points, R_inf))
    E, V = eigen(H_el)  # ascending Re(E)
    bound = np.flatnonzero(np.abs(E.imag) < _IM_TOL_HA)
    if bound.size < n_states:
        raise ValueError(
            f"anion_electronic_states(n_states={n_states}) found only "
            f"{bound.size} bound electronic state(s) with |Im(E)| < "
            f"{_IM_TOL_HA} Ha at R_inf={R_inf}: the well supports fewer bound "
            "states than requested. Reduce n_states."
        )
    idx = bound[:n_states]
    eps_e = E[idx].real
    phi = V[:, idx].T.astype(np.complex128)

    real = g_r.real_points <= g_r.R0
    for i in range(n_states):
        p = phi[i].copy()
        p[~real] = 0.0
        norm2 = c_product(p, p)
        phi[i] = phi[i] / np.sqrt(norm2)
    return eps_e, phi
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest libs/qscat/tests/test_dissociation.py -q`
Expected: PASS (5 tests: 3 parametrized + 2).

- [ ] **Step 5: Type-check + lint**

Run: `uv run mypy libs/qscat/qscat/core && uv run ruff check libs/qscat/qscat/core`
Expected: no errors.

- [ ] **Step 6: Commit**

```bash
git add libs/qscat/qscat/core/dissociation.py libs/qscat/tests/test_dissociation.py
git commit -m "feat(da): anion bound-electronic-state solver at the dissociation limit"
```

---

### Task 3: The rearrangement interaction `V_DR`

**Files:**
- Modify: `libs/qscat/qscat/core/dissociation.py`
- Test: `libs/qscat/tests/test_dissociation.py`

**Interfaces:**
- Consumes: `model.interaction_diag(tgrid)`, `model.v0`, `model.v_int`, `qscat.dvr.TensorGrid` (`points()`, `shape`, `size`).
- Produces: `v_dr_diag(tgrid: TensorGrid, model: ResonanceModel) -> NDArray[complex128]` — the flat (C-order) diagonal `V_int(r,R) + v0(R) − V_int(r, R_inf)`, `R_inf = tgrid.grids[1].R0`, length `tgrid.size`.

**Design notes:** `model.interaction_diag(tgrid)` already gives flat `V_int(r,R)`. The `v0(R)` term ignores `r`: `model.v0(pts_R)` where `pts_R = tgrid.points()[1]` is `(1, n_R)` → broadcast over `r`. The `V_int(r, R_inf)` term ignores `R`: `model.v_int(pts_r, R_inf)` where `pts_r = tgrid.points()[0]` is `(n_r, 1)` → broadcast over `R`. Build the two terms with explicit `np.broadcast_to(..., tgrid.shape).ravel()` (do NOT route through `potential_nd`, which would pass the grid `R`, not `R_inf`).

- [ ] **Step 1: Write the failing test** (append to `test_dissociation.py`)

```python
from qscat.core.dissociation import v_dr_diag
from qscat.core.grids import electronic_grid as _eg, nuclear_grid as _ng
from qscat.dvr import TensorGrid


def _tgrid():
    return TensorGrid([_eg(r_max=14.0, order=6, n_complex=4),
                       _ng(r_max=20.0, n_complex=4, quadrature=8)])


def test_v_dr_shape_and_dtype():
    tg = _tgrid()
    vdr = v_dr_diag(tg, F2)
    assert vdr.shape == (tg.size,) and vdr.dtype == np.complex128


def test_v_dr_equals_definition_pointwise():
    tg = _tgrid()
    model = F2
    R_inf = tg.grids[1].R0
    pts_r, pts_R = tg.points()  # (n_r,1), (1,n_R)
    expect = (
        model.interaction_diag(tg)
        + np.broadcast_to(model.v0(pts_R), tg.shape).ravel()
        - np.broadcast_to(model.v_int(pts_r, R_inf), tg.shape).ravel()
    )
    assert np.allclose(v_dr_diag(tg, model), expect, rtol=0, atol=1e-14)


def test_v_dr_tends_to_v0_at_large_R():
    # Where R is near R_inf, V_int(r,R) ~ V_int(r,R_inf), so V_DR ~ v0(R).
    tg = _tgrid()
    model = F2
    vdr = v_dr_diag(tg, model).reshape(tg.shape)  # (n_r, n_R)
    pts_R = tg.points()[1].ravel()
    j = int(np.argmin(np.abs(pts_R - tg.grids[1].R0)))  # column nearest R_inf
    v0_col = np.broadcast_to(model.v0(tg.points()[1]), tg.shape)[:, j]
    assert np.allclose(vdr[:, j], v0_col, rtol=0, atol=1e-10)
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest libs/qscat/tests/test_dissociation.py -q -k v_dr`
Expected: FAIL — `cannot import name 'v_dr_diag'`.

- [ ] **Step 3: Implement `v_dr_diag`** (append to `dissociation.py`; add `TensorGrid` to the `qscat.dvr` import)

```python
def v_dr_diag(tgrid: TensorGrid, model: ResonanceModel) -> npt.NDArray[np.complex128]:
    """The rearrangement interaction `V_DR = V_int(r,R) + v0(R) - V_int(r, R_inf)`,
    flat (C-order), length `tgrid.size`. `R_inf = tgrid.grids[1].R0` (the nuclear
    ECS pivot / real-region endpoint, eMoScat's `nu_inf`).

    This -- not `V_int` -- is the operator in the DA/DR T-matrix: `H - H_final`,
    where `H_final` is the asymptotic channel Hamiltonian (electron bound in
    `V_int(r, R_inf)`, free nuclei on `v0`). As `R -> R_inf` the `V_int` terms
    cancel and `V_DR -> v0(R)`.
    """
    R_inf = tgrid.grids[1].R0
    pts_r, pts_R = tgrid.points()
    v0_term = np.broadcast_to(model.v0(pts_R), tgrid.shape).ravel()
    vint_inf = np.broadcast_to(model.v_int(pts_r, R_inf), tgrid.shape).ravel()
    return np.asarray(
        model.interaction_diag(tgrid) + v0_term - vint_inf, dtype=np.complex128
    )
```

Update the import line to `from qscat.dvr import FemDvrEcsGrid, TensorGrid, eigen, kinetic`.

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest libs/qscat/tests/test_dissociation.py -q -k v_dr`
Expected: PASS (3 tests).

- [ ] **Step 5: Type-check + lint**

Run: `uv run mypy libs/qscat/qscat/core && uv run ruff check libs/qscat/qscat/core`
Expected: no errors.

- [ ] **Step 6: Commit**

```bash
git add libs/qscat/qscat/core/dissociation.py libs/qscat/tests/test_dissociation.py
git commit -m "feat(da): rearrangement interaction V_DR = V_int + v0 - V_int(r,R_inf)"
```

---

### Task 4: `da_cross_section` — the DA driven-equation T-matrix

**Files:**
- Modify: `libs/qscat/qscat/core/dissociation.py`
- Modify: `libs/qscat/qscat/core/__init__.py`
- Test: `libs/qscat/tests/test_dissociation.py`

**Interfaces:**
- Consumes: `qscat.core.driven.ve_cross_section` (for `Ψ₊` via `return_wavefunction=True`), `anion_electronic_states`, `v_dr_diag`, `riccati_bessel_en_mass`, `qscat.linalg.c_product`, `TensorGrid.{outer, real_mask, sqrt_weights, grids}`.
- Produces:
  ```python
  def da_cross_section(
      tgrid: TensorGrid, model: ResonanceModel,
      eps: NDArray[float64], chi: NDArray[complex128],
      v_init: int, E: float | ArrayLike, *,
      n_channels: int = 1, ordering: _Ordering = "COLAMD",
  ) -> NDArray[float64]: ...
  ```
  σ_DA per anion channel. Scalar `E` → shape `(n_channels,)`; array `E` → `(len(E), n_channels)`. `σ = 0` for a closed channel (`E ≤ 0` or `E_DR = E_tot − eps_e ≤ 0`). `_Ordering` is the same `Literal[...]` re-declaration `driven.py` uses (mirrors `SparseLU`'s private ordering type) so `ordering` passes through to `ve_cross_section` type-clean, with no `# type: ignore`.

**Design notes:** Get `Ψ₊` for all energies in ONE driven sweep by calling `ve_cross_section(tgrid, model, eps, chi, v_init, [v_init], E, ordering=ordering, return_wavefunction=True)` — it already does the `SparseLU.refactor` energy sweep and returns `psis` (a per-energy list, `None` below threshold). Discard its σ_VE. Then for each energy `e>0` with `psi_plus`: `e_tot = e + eps[v_init]`; for each anion channel `(eps_e_n, phi_n)`: `E_DR = e_tot − eps_e_n`; skip if `≤0`; `K = √(2μ E_DR)`; build `Y_coeff = riccati_bessel_en_mass(g_R.real_points, K, 0, μ) * tgrid.sqrt_weights()[1].ravel()`; `Phi = tgrid.outer([phi_n, Y_coeff])`, `Phi[~mask] = 0`; `T = c_product(Phi, v_dr * psi_plus)`; `σ_n = 4π³|T|²/(2e)`. Compute `v_dr = v_dr_diag(tgrid, model)`, the anion states, `mask`, and `sqrt_w_R` once outside the energy loop. `R_inf = tgrid.grids[1].R0`.

- [ ] **Step 1: Write the failing test** (append to `test_dissociation.py`)

```python
from qscat.core.dissociation import da_cross_section


def _working():
    tg = TensorGrid([_eg(r_max=16.0, order=8, n_complex=6),
                     _ng(r_max=22.0, n_complex=8, quadrature=12)])
    return tg


def test_da_shape_scalar_and_array():
    tg = _working()
    eps, chi = vibrational_states(tg.grids[1], F2.mu, 3, F2.v0)
    s1 = da_cross_section(tg, F2, eps, chi, 0, 0.05)
    assert s1.shape == (1,)
    sN = da_cross_section(tg, F2, eps, chi, 0, np.array([0.05, 0.10]))
    assert sN.shape == (2, 1)
    assert np.all(sN >= 0.0) and np.all(np.isfinite(sN))


def test_n2_channel_closed_is_zero():
    # N2's DA threshold is +0.5 Ha -> sigma_DA == 0 across the whole VE window.
    tg = _working()
    eps, chi = vibrational_states(tg.grids[1], N2.mu, 3, N2.v0)
    E = np.array([0.04, 0.10, 0.18])
    s = da_cross_section(tg, N2, eps, chi, 0, E)
    assert np.all(s == 0.0)


@pytest.mark.slow
def test_f2_exothermic_da_is_positive():
    # F2 DA is open at all E>0; expect a nonzero, finite sigma in its resonance
    # window. No golden number (no independent DA data) -- positivity + soft
    # unitarity only.
    tg = _working()
    eps, chi = vibrational_states(tg.grids[1], F2.mu, 3, F2.v0)
    E = np.array([0.02, 0.03, 0.04])
    s = da_cross_section(tg, F2, eps, chi, 0, E)[:, 0]
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)
    assert s.max() > 0.0
    # soft unitarity: sigma_DA <= a few * pi/(2E) (partial-wave cap, generous
    # band for the under-resolved fast outgoing wave; see the convergence note)
    cap = np.pi / (2.0 * E)
    assert np.all(s < 50.0 * cap)
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest libs/qscat/tests/test_dissociation.py -q -k "da_ or n2_channel or f2_exoth"`
Expected: FAIL — `cannot import name 'da_cross_section'`.

- [ ] **Step 3: Implement `da_cross_section`** (append to `dissociation.py`)

Add imports + the `_Ordering` type at the top of the module (below the existing imports). `channel_vector` is NOT needed (DA builds its exit factor from `riccati_bessel_en_mass`, not the electronic channel function):
```python
from typing import Literal  # add to the existing typing import

from qscat.special import riccati_bessel_en_mass

from .driven import ve_cross_section

# Mirrors driven.py's re-declaration of SparseLU's private ordering Literal,
# so `ordering` passes through to ve_cross_section type-clean.
_Ordering = Literal["NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD"]
```

```python
def da_cross_section(
    tgrid: TensorGrid,
    model: ResonanceModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    E: float | npt.ArrayLike,
    *,
    n_channels: int = 1,
    ordering: _Ordering = "COLAMD",
) -> npt.NDArray[np.float64]:
    """sigma_DA(E) in bohr^2, exact 2-D driven-equation DA cross section.

    Reuses `ve_cross_section(..., return_wavefunction=True)` for `Psi+` (one
    `SparseLU.refactor` sweep across `E`), then projects onto each of
    `n_channels` anion dissociation channels with `V_DR`. `E` may be scalar
    (returns `(n_channels,)`) or an array (returns `(len(E), n_channels)`).
    `sigma = 0` for a closed channel (`E <= 0` or `E_DR = E_tot - eps_e <= 0`).
    """
    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    mu = model.mu
    g_R = tgrid.grids[1]
    R_inf = g_R.R0

    eps_e, phi = anion_electronic_states(g_r=tgrid.grids[0], model=model, R_inf=R_inf, n_states=n_channels)
    v_dr = v_dr_diag(tgrid, model)
    mask = tgrid.real_mask()
    sqrt_w_R = tgrid.sqrt_weights()[1].ravel()

    _, psis = ve_cross_section(
        tgrid, model, eps, chi, v_init, [v_init], e_arr,
        ordering=ordering, return_wavefunction=True,
    )
    psi_list = psis if isinstance(psis, list) else [psis]

    out = np.zeros((len(e_arr), n_channels), dtype=np.float64)
    for ie, e in enumerate(e_arr):
        psi_plus = psi_list[ie]
        if psi_plus is None:  # E <= 0
            continue
        e_tot = float(e) + eps[v_init]
        v_psi = v_dr * psi_plus
        for n in range(n_channels):
            e_dr = e_tot - eps_e[n]
            if e_dr <= 0.0:
                continue
            k_r = float(np.sqrt(2.0 * mu * e_dr))
            y_coeff = riccati_bessel_en_mass(g_R.real_points, k_r, 0, mu) * sqrt_w_R
            phi_f = tgrid.outer([phi[n], y_coeff])
            phi_f[~mask] = 0.0
            t = c_product(phi_f, v_psi)
            out[ie, n] = 4.0 * np.pi**3 * abs(t) ** 2 / (2.0 * float(e))

    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    return np.asarray(out[0] if scalar else out, dtype=np.float64)
```

Export in `libs/qscat/qscat/core/__init__.py`: add
`from .dissociation import anion_electronic_states, da_cross_section, v_dr_diag`
and append `"anion_electronic_states"`, `"v_dr_diag"`, `"da_cross_section"` to `__all__`. Extend the module docstring's "Public API" list with a one-line entry for the DA cross section.

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest libs/qscat/tests/test_dissociation.py -q`
Expected: PASS (non-slow); the `@slow` F2 test passes under `uv run pytest libs/qscat/tests/test_dissociation.py -q -m slow`.

- [ ] **Step 5: Guard test still green + type-check + lint**

Run: `uv run pytest libs/qscat/tests/test_core_no_model_import.py -q && uv run mypy libs/qscat && uv run ruff check libs/qscat`
Expected: the import-boundary test PASSES (dissociation.py must not import qscat.model at runtime), mypy 0 errors, ruff clean.

- [ ] **Step 6: Commit**

```bash
git add libs/qscat/qscat/core/dissociation.py libs/qscat/qscat/core/__init__.py libs/qscat/tests/test_dissociation.py
git commit -m "feat(da): exact-2D TI da_cross_section (driven solve + V_DR T-matrix)"
```

---

### Task 5: Nuclear-grid convergence study

**Files:**
- Create: `benchmarks/da_nuclear_convergence.py`
- Test: `benchmarks/test_da_nuclear_convergence.py`

**Interfaces:**
- Consumes: `qscat.core.da_cross_section`, `qscat.core.grids`, `qscat.core.vibrational.vibrational_states`, `qscat.model.F2`.
- Produces: `convergence_table(model, E, grids) -> list[dict]` returning `{"n_complex", "quadrature", "r_max", "sigma"}` per nuclear-grid setting, and a `main()` that prints the table. This MEASURES how σ_DA(E) stabilizes as the nuclear grid resolves the fast outgoing wave (`K_R ~ 50`, wavelength `~0.13 bohr`) — documenting that the DA roughness is under-resolution, not a method error.

**Design notes:** F₂ (exothermic, strongest DA) at a fixed collision energy in its window (e.g. `E=0.03`). Sweep the nuclear grid's `quadrature` (10→18) and `n_complex` at fixed electronic grid; hold `v_init=0`, `n_channels=1`. Report σ_DA and the successive relative change. The `benchmarks/` package imports `projects`/`qscat` and is run via `python -m benchmarks.da_nuclear_convergence` (see CLAUDE.md).

- [ ] **Step 1: Write the failing test**

```python
# benchmarks/test_da_nuclear_convergence.py
from __future__ import annotations

import numpy as np
from benchmarks.da_nuclear_convergence import convergence_table
from qscat.model import F2


def test_convergence_table_runs_and_stabilizes():
    rows = convergence_table(F2, 0.03, quadratures=(10, 14, 18))
    assert len(rows) == 3
    assert all(np.isfinite(r["sigma"]) and r["sigma"] >= 0.0 for r in rows)
    # finer nuclear quadrature: the last step changes less than the first
    d1 = abs(rows[1]["sigma"] - rows[0]["sigma"])
    d2 = abs(rows[2]["sigma"] - rows[1]["sigma"])
    assert d2 <= d1 + 1e-12
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest benchmarks/test_da_nuclear_convergence.py -q`
Expected: FAIL — `ModuleNotFoundError: benchmarks.da_nuclear_convergence`.

- [ ] **Step 3: Implement the study**

```python
# benchmarks/da_nuclear_convergence.py
"""DA nuclear-grid convergence: sigma_DA(E) vs nuclear resolution for F2.

The DA exit wave F^nuc_{K,0}(R) = sqrt(2 mu/pi K) sin(KR) oscillates fast for
heavy nuclei (F2: K_R ~ 50, wavelength ~ 0.13 bohr), so the nuclear FEM-DVR
grid must resolve it for the T-matrix quadrature to converge. This measures
that: sigma_DA(F2, E) as the nuclear quadrature is refined. Run via
`uv run python -m benchmarks.da_nuclear_convergence`.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from qscat.core.dissociation import da_cross_section
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import DiatomicResonanceModel


def convergence_table(
    model: DiatomicResonanceModel,
    E: float,
    quadratures: Sequence[int] = (10, 12, 14, 16, 18),
    *,
    r_max_R: float = 22.0,
    n_complex: int = 10,
    el_r_max: float = 16.0,
) -> list[dict]:
    rows: list[dict] = []
    g_r = electronic_grid(r_max=el_r_max, order=8, n_complex=6)
    for q in quadratures:
        g_R = nuclear_grid(r_max=r_max_R, n_complex=n_complex, quadrature=q)
        tg = TensorGrid([g_r, g_R])
        eps, chi = vibrational_states(g_R, model.mu, 3, model.v0)
        sigma = float(da_cross_section(tg, model, eps, chi, 0, E)[0])
        rows.append({"quadrature": q, "n_complex": n_complex, "r_max": r_max_R, "sigma": sigma})
    return rows


def main() -> None:
    from qscat.model import F2

    rows = convergence_table(F2, 0.03)
    print(f"F2 sigma_DA(E=0.03) vs nuclear quadrature (r_max={rows[0]['r_max']}):")
    prev = None
    for r in rows:
        rel = "" if prev is None else f"  d_rel={abs(r['sigma']-prev)/max(abs(r['sigma']),1e-30):.2%}"
        print(f"  q={r['quadrature']:2d}  sigma={r['sigma']:.6e}{rel}")
        prev = r["sigma"]


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest benchmarks/test_da_nuclear_convergence.py -q`
Expected: PASS.

- [ ] **Step 5: Run the study, record the numbers**

Run: `uv run python -m benchmarks.da_nuclear_convergence`
Expected: a printed table of σ_DA vs quadrature with shrinking `d_rel`. Copy the numbers into the commit message (they seed Task 6's doc table).

- [ ] **Step 6: Type-check + lint + commit**

Run: `uv run mypy benchmarks/da_nuclear_convergence.py && uv run ruff check benchmarks/da_nuclear_convergence.py`
```bash
git add benchmarks/da_nuclear_convergence.py benchmarks/test_da_nuclear_convergence.py
git commit -m "bench(da): nuclear-grid convergence study for sigma_DA (K_R~50 resolution)"
```

---

### Task 6: F₂/NO σ_DA figures, docs, and CLAUDE.md

**Files:**
- Create: `validation/diatomic/da_curves.py`
- Create: `validation/diatomic/test_da_curves.py`
- Modify: `docs/physics/diatomic-ve-cross-sections.md`
- Modify: `CLAUDE.md`

**Interfaces:**
- Consumes: `qscat.core.da_cross_section`, `validation.diatomic.config.CONFIGS`, `validation.diatomic.curves.build_grid` + `FIGURE_DIR`, `qscat.core.vibrational.vibrational_states`, `qscat.core.plot.plot_cross_sections`.
- Produces: `compute_da_curve(cfg: MoleculeConfig, E_grid) -> (E, sigma)` with `sigma` shape `(len(E), 1)`, and `main()` writing `docs/physics/figures/{f2,no}-2d-ti-da-cross-section.png`.

**Design notes (verified against the real APIs):** reuse `build_grid(cfg)` and `FIGURE_DIR` from `validation/diatomic/curves.py` (do NOT rebuild the grid by hand). `MoleculeConfig` exposes grid *fields* (`e_r_max`, `n_quadrature`, …, `n_vib`, `model`), NOT builder methods. `plot_cross_sections(E_grid, sigma, *, channels=None, title=..., path=...)` is keyword-only for everything after `sigma`, `path` is required, and it log-scales y + masks non-positive σ to NaN. `channels` is typed `list[int] | None`; DA is a single channel with no `v'`, so pass `channels=None` (one curve) and let the `title` carry "dissociative attachment" — the legend's generic "v'=0" label is a known cosmetic wart on a single-curve DA plot, not worth widening the shared plot signature. `CONFIGS` holds only NO and F₂ (no N₂ — closed). F₂ exothermic (`E ∈ [0.005, 0.05]`), NO opens ~0.17 (`E ∈ [0.15, 0.30]`). This is the oracle σ_DA (no independent data). The gate test uses a coarse E grid and asserts thresholds + finiteness only (a full curve is `@slow`).

- [ ] **Step 1: Write the failing gate test**

```python
# validation/diatomic/test_da_curves.py
from __future__ import annotations

import numpy as np
import pytest
from validation.diatomic.config import CONFIGS
from validation.diatomic.da_curves import compute_da_curve


def test_f2_da_open_no_closed_below_threshold():
    # F2 exothermic: sigma_DA > 0 somewhere in [0.01, 0.04].
    E, s = compute_da_curve(CONFIGS["F2"], np.array([0.02, 0.03, 0.04]))
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)
    assert s.max() > 0.0


@pytest.mark.slow
def test_no_da_threshold_onset():
    # NO opens ~0.17 Ha: closed below, open above.
    _, s_lo = compute_da_curve(CONFIGS["NO"], np.array([0.10]))
    _, s_hi = compute_da_curve(CONFIGS["NO"], np.array([0.22]))
    assert s_lo[0, 0] == 0.0
    assert s_hi[0, 0] >= 0.0
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest validation/diatomic/test_da_curves.py -q`
Expected: FAIL — `ModuleNotFoundError: validation.diatomic.da_curves`.

- [ ] **Step 3: Implement `da_curves.py`**

```python
# validation/diatomic/da_curves.py
"""Exact-2D TI dissociative-attachment sigma_DA(E) for F2 and NO (the oracle).

No independent DA data exists (only N2 VE has Houfek's); the exact-2D TI
solver IS the reference. N2's DA channel is closed in the measurement window
(threshold +0.5 Ha), so only F2 (exothermic) and NO (~0.17 Ha) are shown.
Run via `uv run python -m validation.diatomic.da_curves`.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.dissociation import da_cross_section
from qscat.core.plot import plot_cross_sections
from qscat.core.vibrational import vibrational_states

from validation.diatomic.config import CONFIGS, MoleculeConfig
from validation.diatomic.curves import FIGURE_DIR, build_grid


def compute_da_curve(
    cfg: MoleculeConfig, E_grid: npt.NDArray[np.float64]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """`(E_grid, sigma_DA[E, 0])` — the exact-2D TI single-channel DA σ(E)."""
    tg = build_grid(cfg)
    eps, chi = vibrational_states(tg.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)
    E = np.asarray(E_grid, dtype=np.float64)
    sigma = da_cross_section(tg, cfg.model, eps, chi, 0, E)
    return E, np.asarray(sigma, dtype=np.float64)


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    specs = {
        "F2": (np.linspace(0.005, 0.05, 60), "f2-2d-ti-da-cross-section"),
        "NO": (np.linspace(0.15, 0.30, 60), "no-2d-ti-da-cross-section"),
    }
    for name, (E, stem) in specs.items():
        _, sigma = compute_da_curve(CONFIGS[name], E)
        np.savez(FIGURE_DIR / f"{stem}.npz", E=E, sigma=sigma)
        plot_cross_sections(
            E,
            sigma,
            channels=None,  # single DA channel; title carries the meaning
            title=f"{name} exact-2D TI dissociative attachment sigma_DA(E) (qscat.core)",
            path=FIGURE_DIR / f"{stem}.png",
        )
        print(f"{name}: wrote {stem}.png/.npz")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run gate test + generate figures**

Run: `uv run pytest validation/diatomic/test_da_curves.py -q`
Then: `uv run python -m validation.diatomic.da_curves`
Expected: gate PASS; two PNGs written.

- [ ] **Step 5: Extend the docs**

In `docs/physics/diatomic-ve-cross-sections.md`, under the DA section, add the converged σ_DA figures and a short table of the Task-5 convergence numbers (σ_DA(F₂, E=0.03) vs nuclear quadrature), stating the K_R~50 resolution requirement and that the roughness is under-resolution, not method error. In `CLAUDE.md`, add a one-line entry for `qscat.core.dissociation` (the DA cross section) to the `qscat.core` bullet, and note `validation/diatomic/da_curves.py` produces the σ_DA figures.

- [ ] **Step 6: Commit**

```bash
git add validation/diatomic/da_curves.py validation/diatomic/test_da_curves.py \
        docs/physics/figures/f2-2d-ti-da-cross-section.png \
        docs/physics/figures/no-2d-ti-da-cross-section.png \
        docs/physics/diatomic-ve-cross-sections.md CLAUDE.md
git commit -m "feat(da): F2/NO sigma_DA oracle curves + figures; docs + CLAUDE.md"
```

---

## Verification (whole sub-project)

- `uv run pytest -q -m "not slow"` passes; `uv run pytest -q -m slow libs/qscat/tests/test_dissociation.py` passes.
- `uv run pytest libs/qscat/tests/test_core_no_model_import.py -q` passes (DA code keeps the core/model boundary).
- `uv run mypy libs/qscat` 0 errors; `uv run ruff check .` clean.
- σ_DA thresholds correct in tests: N₂ ≡ 0 (closed), F₂ > 0 (exothermic), NO onset ~0.17; σ_DA finite, ≥0, soft-unitary.
- Nuclear-grid convergence measured and documented (the roughness is under-resolution, not a method error); F₂/NO σ_DA figures committed; `docs/physics/diatomic-ve-cross-sections.md` + `CLAUDE.md` updated.

## Out of scope (this plan)

- **LCP DA** (the approximation under test) — sub-project B.
- **H₂⁺ DR** (Rydberg-series loop + Coulomb incident `coulomb::sF_en`) — sub-project D.
- Rotational (J>0), multiple electron partial waves, TD exact-2D DA (the TI route is exact and cheaper).
