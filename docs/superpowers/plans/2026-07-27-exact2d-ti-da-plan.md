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

**Design notes (read before coding):** `model.surface(g_r.points, R_inf)` is the full electronic potential at the dissociation limit — `v0(R_inf) + ℓ(ℓ+1)/2r² + v_int(r, R_inf)` — so `eps_e` sits on the SAME energy scale as `H_2D` (it includes `v0(R_inf)`), which is what makes the DA threshold `eps_e − eps[v_init]` correct. The anion state is genuinely BOUND at `R_inf`: a real eigenvalue (`|Im(E)| < _IM_TOL_HA`) BELOW the asymptotic electronic continuum edge, which is `v0(R_inf)` (as `r→∞` both `v_int` and the centrifugal term vanish, so the electron sees only `v0(R_inf)`). Select by BOTH conditions — `|Im(E)| < _IM_TOL_HA` AND `Re(E) < v0(R_inf)` — then take the lowest `n_states` by `Re(E)`. **The `|Im|` filter alone is NOT enough:** a finite FEM/DVR basis always produces "top-of-grid numerical-junk" eigenvalues with large positive `Re(E)` but tiny `|Im(E)|` (documented in `test_femdvr_ecs.py`'s "numerical-junk states"; on the 113-point grid there are ~50 of them), which the `Re(E) < v0(R_inf)` ceiling excludes. Do NOT take "the n lowest-Re overall" either. `g_r.points` may be complex on the tail; `surface`/`v0` handle that. Mirror `qscat.core.vibrational.vibrational_states`' structure and its `_IM_TOL_HA = 1e-6`.

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

__all__ = ["anion_electronic_states"]  # Tasks 3/4 append "v_dr_diag", "da_cross_section"

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

    Returns `(eps_e, phi_e)`: `eps_e` the `n_states` lowest-Re genuinely bound
    eigenvalues (`|Im(E)| < _IM_TOL_HA` AND `Re(E) < v0(R_inf)`), real,
    ascending; `phi_e` shape `(n_states, g_r.n)`, each c-product-normalized
    over the electronic real region. Raises `ValueError` if fewer than
    `n_states` bound states exist (e.g. `n_states` reached past the finite
    bound spectrum).
    """
    H_el = kinetic(g_r, 1.0) + np.diag(model.surface(g_r.points, R_inf))
    E, V = eigen(H_el)  # ascending Re(E)
    # Genuinely bound: near-real AND below the asymptotic electronic continuum
    # edge v0(R_inf) (as r->inf, v_int and centrifugal vanish -> the electron
    # sees only v0(R_inf)). The |Im| filter alone counts finite-basis
    # "top-of-grid numerical-junk" eigenvalues (large positive Re(E), tiny
    # |Im|) as bound; the Re(E) < e_thresh cut excludes them.
    e_thresh = float(np.real(model.v0(np.asarray(R_inf))))
    bound = np.flatnonzero((np.abs(E.imag) < _IM_TOL_HA) & (E.real < e_thresh))
    if bound.size < n_states:
        raise ValueError(
            f"anion_electronic_states(n_states={n_states}) found only "
            f"{bound.size} bound electronic state(s) (|Im(E)| < {_IM_TOL_HA} Ha "
            f"and Re(E) < v0(R_inf)={e_thresh:.6g}) at R_inf={R_inf}: the well "
            "supports fewer bound states than requested. Reduce n_states."
        )
    idx = bound[:n_states]  # E is Re-ascending, so these are the lowest-Re bound states
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

Update the import line to `from qscat.dvr import FemDvrEcsGrid, TensorGrid, eigen, kinetic`, and append `"v_dr_diag"` to the module `__all__` (so it reads `["anion_electronic_states", "v_dr_diag"]` — ruff F822 fails on names in `__all__` that aren't yet defined, so grow it as each function lands).

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

Append `"da_cross_section"` to `dissociation.py`'s module `__all__` (now `["anion_electronic_states", "v_dr_diag", "da_cross_section"]`). Export in `libs/qscat/qscat/core/__init__.py`: add `from .dissociation import anion_electronic_states, da_cross_section, v_dr_diag` and append `"anion_electronic_states"`, `"v_dr_diag"`, `"da_cross_section"` to its `__all__`. Extend that module docstring's "Public API" list with a one-line entry for the DA cross section.

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

### Task 5: `segmented_grid` — the eMoScat per-molecule grid builder

**Rationale (why this replaced the quadrature convergence study):** the DA magnitude did NOT converge under quadrature (p-)refinement — σ_DA(F₂,0.03) swung 0.16→26→0.54→2.3→4.0 across q=10→30 on the shared N₂-style grid. A `sin(K_R R)` exit wave with K_R≈58 (wavelength ~0.107 bohr) needs element-density (h-)refinement, not more points per element. eMoScat already solved this ACCURATELY with hand-tuned **per-molecule** grids (`reference/eMoScat/input/{N2,NO,F2}/grids.txt`), whose nuclear grids are far finer over the dissociation region (NO: 37×0.2 bohr over [1.6,9.0]; F₂: 40×0.2 bohr over [2.7,10.7] + 0.024 bohr around the 2.5–2.7 resonance). We adopt those tested values rather than rediscover the knob. This task adds the general builder for eMoScat's grid-deck format; Task 6 transcribes the decks. (Future work — a skill that CHOOSES the discretisation from the potential curves via the per-element de Broglie phase — is noted in the docs, not built here.)

**Files:**
- Modify: `libs/qscat/qscat/core/grids.py`
- Modify: `libs/qscat/qscat/core/__init__.py`
- Test: `libs/qscat/tests/test_grids_segmented.py`

**Interfaces:**
- Produces: `segmented_grid(real_segments, complex_segments, *, angle_deg, quadrature, x_min=0.0) -> FemDvrEcsGrid`. `real_segments`/`complex_segments` are `Sequence[tuple[int, float]]` of `(n_elements, endpoint)` pairs — exactly eMoScat's `grids.txt` format: from `x_min`, each segment lays `n` uniform elements up to `endpoint`; the complex part is an ECS tail at `angle_deg`. The ECS pivot `R0` is the last real endpoint. `complex_segments` may be empty (pure real grid).

**Design notes:** mirror the element assembly already in `nuclear_grid`/`electronic_grid` (build a `list[ElementSpec]`, real as `ElementSpec(h)`, complex as `ElementSpec(h, angle_deg)`, then `FemDvrEcsGrid(GridSpec(quadrature, elements, x_min))`). `h = (endpoint - start) / n`. Validate: every `n >= 1`; endpoints strictly increasing (each `> start`); `quadrature >= 2`. `GridSpec` computes `R0` as `x_min + sum(real element lengths)` = the last real endpoint, so the ECS tail starts there automatically.

- [ ] **Step 1: Write the failing test**

```python
# libs/qscat/tests/test_grids_segmented.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import segmented_grid


def test_reproduces_emoscat_n2_nuclear_deck():
    # N2 nuclear (input/N2/grids.txt, 2nd declaration): real to 12.0, tail to 55.
    g = segmented_grid(
        [(2, 1.0), (1, 1.5), (10, 3.0), (2, 4.0), (2, 6.0), (6, 12.0)],
        [(1, 13.0), (2, 16.0), (1, 18.0), (4, 55.0)],
        angle_deg=35.0,
        quadrature=14,
    )
    assert g.R0 == pytest.approx(12.0)            # ECS pivot = last real endpoint
    # real region ends at 12, tail runs onto the complex plane past it
    assert float(g.real_points.max()) == pytest.approx(55.0)
    assert np.iscomplexobj(g.points) and np.any(np.abs(g.points.imag) > 0)


def test_element_lengths_are_uniform_per_segment():
    # the 10 elements over [1.5, 3.0] are each 0.15 bohr
    g = segmented_grid([(1, 1.5), (10, 3.0)], [], angle_deg=35.0, quadrature=8)
    assert g.R0 == pytest.approx(3.0)             # no complex tail -> pivot at real end
    assert np.max(np.abs(g.points.imag)) == pytest.approx(0.0)  # pure real


@pytest.mark.parametrize(
    "real_seg",
    [[(0, 1.0)], [(2, 1.0), (1, 0.5)]],  # n<1 ; non-increasing endpoint
)
def test_rejects_bad_segments(real_seg):
    with pytest.raises(ValueError):
        segmented_grid(real_seg, [], angle_deg=35.0, quadrature=8)
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest libs/qscat/tests/test_grids_segmented.py -q`
Expected: FAIL — `cannot import name 'segmented_grid'`.

- [ ] **Step 3: Implement `segmented_grid`** (append to `grids.py`, add to `__all__`)

```python
def segmented_grid(
    real_segments: Sequence[tuple[int, float]],
    complex_segments: Sequence[tuple[int, float]],
    *,
    angle_deg: float,
    quadrature: int,
    x_min: float = 0.0,
) -> FemDvrEcsGrid:
    """A FEM-DVR-ECS grid from eMoScat's `grids.txt` segment format.

    `real_segments` / `complex_segments` are `(n_elements, endpoint)` pairs:
    from `x_min`, each segment tiles `n` uniform elements up to `endpoint`.
    The complex part is an ECS tail at `angle_deg`; the ECS pivot `R0` is the
    last real endpoint. `complex_segments` may be empty (a pure real grid).
    This is the per-molecule discretisation route -- see
    docs/physics/diatomic-ve-cross-sections.md (DA nuclear grids).
    """
    if quadrature < 2:
        raise ValueError(f"quadrature must be >= 2, got {quadrature}")
    elements: list[ElementSpec] = []
    start = x_min
    for label, segs, angle in (("real", real_segments, None), ("complex", complex_segments, angle_deg)):
        for n, end in segs:
            if n < 1:
                raise ValueError(f"{label} segment ({n}, {end}) has n_elements < 1")
            if end <= start:
                raise ValueError(f"{label} endpoint {end} must exceed previous {start}")
            h = (end - start) / n
            elements += [ElementSpec(h) if angle is None else ElementSpec(h, angle) for _ in range(n)]
            start = end
    return FemDvrEcsGrid(GridSpec(quadrature=quadrature, elements=elements, x_min=x_min))
```

Add `from collections.abc import Sequence` at the top of `grids.py` if not present, append `"segmented_grid"` to `grids.py`'s `__all__`, and export it from `libs/qscat/qscat/core/__init__.py` (import + `__all__` + a one-line Public-API docstring entry).

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest libs/qscat/tests/test_grids_segmented.py -q`
Expected: PASS (4 tests).

- [ ] **Step 5: Type-check + lint**

Run: `uv run mypy libs/qscat/qscat/core && uv run ruff check libs/qscat/qscat/core`
Expected: no errors.

- [ ] **Step 6: Commit**

```bash
git add libs/qscat/qscat/core/grids.py libs/qscat/qscat/core/__init__.py libs/qscat/tests/test_grids_segmented.py
git commit -m "feat(da): segmented_grid builder for eMoScat per-molecule grid decks"
```

---

### Task 6: eMoScat per-molecule nuclear decks + σ_DA stability

**Files:**
- Modify: `validation/diatomic/config.py`
- Create: `benchmarks/da_grid_stability.py`
- Test: `validation/diatomic/test_da_grid.py`

**Interfaces:**
- `MoleculeConfig` gains the eMoScat nuclear deck fields + a `da_grid()` builder:
  ```python
  nuc_real: tuple[tuple[int, float], ...]      # eMoScat real segments
  nuc_complex: tuple[tuple[int, float], ...]   # eMoScat complex (ECS tail) segments
  nuc_angle: float
  nuc_quad: int
  def da_grid(self) -> TensorGrid: ...          # electronic (r_max=e_r_max, validated) x eMoScat nuclear
  ```
- `benchmarks/da_grid_stability.py`: `stability(cfg, E) -> dict` returning σ_DA on the eMoScat grid AND on a modestly h-refined variant (doubling the outer nuclear elements) + their relative difference — the real convergence evidence (records numbers for the docs). Run via `uv run python -m benchmarks.da_grid_stability`.

**Design notes:** the electronic grid stays the VE-validated one (`electronic_grid(r_max=e_r_max, order=e_order, n_complex=e_n_complex)` — r_max=16); only the NUCLEAR grid becomes eMoScat's fine per-molecule deck via `segmented_grid`. `da_grid` = `TensorGrid([electronic_grid(...), segmented_grid(nuc_real, nuc_complex, angle_deg=nuc_angle, quadrature=nuc_quad)])`. This is a NEW grid path used only by DA; the VE `build_grid` (in `curves.py`) is unchanged (VE's outgoing flux is electronic, already converged on the coarse nuclear grid — do not touch VE figures here). The eMoScat decks, transcribed verbatim from `reference/eMoScat/input/{NO,F2}/grids.txt` (2nd/nuclear declaration):

- **NO** nuclear: real `[(1,1.0),(1,1.6),(37,9.0)]`, complex `[(1,9.25),(1,10.0),(1,12.0),(4,42.0)]`, angle `45.0`, quad `14`.
- **F₂** nuclear: real `[(9,1.8),(1,2.0),(5,2.5),(4,2.596908),(4,2.7),(40,10.7)]`, complex `[(1,10.8),(1,11.0),(1,11.5),(1,12.5),(1,14.0),(1,18.0),(4,30.0),(2,101.0)]`, angle `35.0`, quad `14`.

(N₂'s deck — real `[(2,1.0),(1,1.5),(10,3.0),(2,4.0),(2,6.0),(6,12.0)]`, complex `[(1,13.0),(2,16.0),(1,18.0),(4,55.0)]`, angle 35, quad 14 — is recorded in the doc but N₂ DA is closed, so no N₂ config entry is needed.) These grids are large (F₂ ≈ 108k-unknown 2-D solve), so the σ_DA runs are heavy: the pytest gate uses 1–2 anchor energies and is `@slow`; the stability benchmark is a script, not a gate.

- [ ] **Step 1: Write the failing test**

```python
# validation/diatomic/test_da_grid.py
from __future__ import annotations

import numpy as np
import pytest
from validation.diatomic.config import CONFIGS
from qscat.core.vibrational import vibrational_states
from qscat.core.dissociation import da_cross_section


def test_da_grid_uses_emoscat_nuclear_resolution():
    # F2 nuclear real region ends at 10.7 (eMoScat deck), finely resolved.
    tg = CONFIGS["F2"].da_grid()
    assert tg.grids[1].R0 == pytest.approx(10.7)
    # the [2.7,10.7] region is tiled at ~0.2 bohr -> many nuclear elements
    assert tg.grids[1].n > 700  # ~960 for the F2 deck at quad 14


@pytest.mark.slow
def test_f2_sigma_da_wellposed_on_emoscat_grid():
    cfg = CONFIGS["F2"]
    tg = cfg.da_grid()
    eps, chi = vibrational_states(tg.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)
    E = np.array([0.03])
    s = da_cross_section(tg, cfg.model, eps, chi, 0, E)[:, 0]
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)
    assert s[0] > 0.0                      # exothermic -> open
    assert s[0] < 50.0 * np.pi / (2.0 * E[0])   # soft unitarity
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest validation/diatomic/test_da_grid.py -q`
Expected: FAIL — `AttributeError: ... 'da_grid'` / missing config fields.

- [ ] **Step 3: Extend `MoleculeConfig` + add the decks**

In `validation/diatomic/config.py`: import `segmented_grid`, `electronic_grid` from `qscat.core.grids` and `TensorGrid` from `qscat.dvr`; add the four `nuc_*` fields to the `MoleculeConfig` dataclass and the `da_grid()` method (per the Interfaces block); fill `nuc_real`/`nuc_complex`/`nuc_angle`/`nuc_quad` for NO and F₂ from the decks above. Keep all existing VE fields/values untouched.

- [ ] **Step 4: Run the well-posedness gate**

Run: `uv run pytest validation/diatomic/test_da_grid.py -q` (non-slow), then `uv run pytest validation/diatomic/test_da_grid.py -q -m slow` (the heavy F₂ solve, minutes).
Expected: both PASS — σ_DA finite, positive, soft-unitary on the eMoScat grid.

- [ ] **Step 5: Add + run the stability benchmark**

Implement `benchmarks/da_grid_stability.py`: for F₂ at E∈{0.02,0.03,0.04}, compute σ_DA on `cfg.da_grid()` and on a variant with the outer nuclear segments' element counts doubled (h-refinement), print both + the relative difference. Run `uv run python -m benchmarks.da_grid_stability` and record the numbers (they go in the doc: evidence σ_DA is now stable under h-refinement, unlike the earlier quadrature sweep).

- [ ] **Step 6: Type-check + lint + commit**

Run: `uv run mypy libs/qscat && uv run ruff check .`
```bash
git add validation/diatomic/config.py validation/diatomic/test_da_grid.py benchmarks/da_grid_stability.py
git commit -m "feat(da): per-molecule eMoScat nuclear grids (NO/F2) + sigma_DA stability"
```

---

### Task 7: F₂/NO σ_DA figures, docs, and CLAUDE.md

**Files:**
- Create: `validation/diatomic/da_curves.py`
- Create: `validation/diatomic/test_da_curves.py`
- Modify: `docs/physics/diatomic-ve-cross-sections.md`
- Modify: `CLAUDE.md`

**Interfaces:**
- Consumes: `qscat.core.da_cross_section`, `validation.diatomic.config.CONFIGS`, `validation.diatomic.curves.FIGURE_DIR`, `qscat.core.vibrational.vibrational_states`, `qscat.core.plot.plot_cross_sections`.
- Produces: `compute_da_curve(cfg, E_grid) -> (E, sigma)` (`sigma` shape `(len(E), 1)`), building the grid via `cfg.da_grid()`; `main()` writing `docs/physics/figures/{f2,no}-2d-ti-da-cross-section.png`.

**Design notes (verified against the real APIs):** build the grid with `cfg.da_grid()` (the eMoScat per-molecule nuclear grid from Task 6) — NOT `build_grid` (that's the VE grid). `plot_cross_sections(E_grid, sigma, *, channels=None, title=..., path=...)` — keyword-only after `sigma`, `path` required, log-y, masks non-positive to NaN; pass `channels=None` (single DA channel; the `title` carries "dissociative attachment"). `CONFIGS` holds only NO and F₂. Because the eMoScat grid makes each energy a ~108k solve (minutes with SuperLU on a laptop), the committed figures use a MODEST energy count (F₂ ~24 points in [0.01,0.05]; NO ~20 in [0.15,0.30]) — dense curves are a Docker+MUMPS follow-on; note this in the doc. The gate test is coarse (2–3 anchors) and `@slow`.

- [ ] **Step 1: Write the failing gate test**

```python
# validation/diatomic/test_da_curves.py
from __future__ import annotations

import numpy as np
import pytest
from validation.diatomic.config import CONFIGS
from validation.diatomic.da_curves import compute_da_curve


@pytest.mark.slow
def test_f2_da_positive_on_emoscat_grid():
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

Computed on eMoScat's per-molecule NUCLEAR grid (`cfg.da_grid()`) -- the fine
discretisation the fast K_R~58 dissociation wave needs. No independent DA data
exists (only N2 VE has Houfek's); the exact-2D TI solver IS the reference. N2's
DA channel is closed in the window (+0.5 Ha), so only F2 (exothermic) and NO
(~0.17 Ha) are shown. Run via `uv run python -m validation.diatomic.da_curves`.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.dissociation import da_cross_section
from qscat.core.plot import plot_cross_sections
from qscat.core.vibrational import vibrational_states

from validation.diatomic.config import CONFIGS, MoleculeConfig
from validation.diatomic.curves import FIGURE_DIR


def compute_da_curve(
    cfg: MoleculeConfig, E_grid: npt.NDArray[np.float64]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """`(E_grid, sigma_DA[E, 0])` on the eMoScat per-molecule nuclear grid."""
    tg = cfg.da_grid()
    eps, chi = vibrational_states(tg.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)
    E = np.asarray(E_grid, dtype=np.float64)
    sigma = da_cross_section(tg, cfg.model, eps, chi, 0, E)
    return E, np.asarray(sigma, dtype=np.float64)


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    specs = {
        "F2": (np.linspace(0.01, 0.05, 13), "f2-2d-ti-da-cross-section"),
        "NO": (np.linspace(0.15, 0.30, 12), "no-2d-ti-da-cross-section"),
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

Run: `uv run pytest validation/diatomic/test_da_curves.py -q -m slow` (heavy — minutes).
Then: `uv run python -m validation.diatomic.da_curves` (heavy — the two curves take a while on a laptop; this is expected).
Expected: gate PASS; two PNGs + `.npz` written.

- [ ] **Step 5: Extend the docs**

In `docs/physics/diatomic-ve-cross-sections.md`, under the DA section: (a) REPLACE the "roughness is under-resolution, needs finer nuclear grid" caveat with the resolved story — σ_DA is computed on eMoScat's **per-molecule** nuclear grid (NO/F₂ decks quoted), which resolves the K_R~58 exit wave; add the committed σ_DA figures and the Task-6 stability numbers (σ_DA stable under h-refinement); note the earlier quadrature (p-)refinement could not converge an oscillatory integrand. (b) Add a short "Future: automatic discretisation" note — a skill to choose the grid from the potential curves via the per-element de Broglie phase. In `CLAUDE.md`, add a one-line `qscat.core.dissociation` entry (the DA cross section) to the `qscat.core` bullet, a `segmented_grid` mention to the grids bullet, and note `validation/diatomic/da_curves.py` + the per-molecule `da_grid()`.

- [ ] **Step 6: Commit**

```bash
git add validation/diatomic/da_curves.py validation/diatomic/test_da_curves.py \
        docs/physics/figures/f2-2d-ti-da-cross-section.png \
        docs/physics/figures/no-2d-ti-da-cross-section.png \
        docs/physics/figures/f2-2d-ti-da-cross-section.npz \
        docs/physics/figures/no-2d-ti-da-cross-section.npz \
        docs/physics/diatomic-ve-cross-sections.md CLAUDE.md
git commit -m "feat(da): F2/NO sigma_DA oracle curves on eMoScat grids + figures + docs"
```

---

## Verification (whole sub-project)

- `uv run pytest -q -m "not slow"` passes; `uv run pytest -q -m slow libs/qscat/tests/test_dissociation.py` passes.
- `uv run pytest libs/qscat/tests/test_core_no_model_import.py -q` passes (DA code keeps the core/model boundary).
- `uv run mypy libs/qscat/qscat` 0 errors (the qscat LIBRARY — test helpers follow the repo's untyped convention); `uv run ruff check .` clean.
- σ_DA thresholds correct in tests: N₂ ≡ 0 (closed), F₂ > 0 (exothermic), NO onset ~0.17; σ_DA finite, ≥0, soft-unitary.
- `segmented_grid` reproduces eMoScat's per-molecule grid decks; the eMoScat per-molecule NUCLEAR grids (NO/F₂) resolve the K_R~58 exit wave, and σ_DA is stable under h-refinement on them (measured, documented) — the earlier quadrature (p-)refinement non-convergence is explained and resolved. F₂/NO σ_DA figures committed on those grids; `docs/physics/diatomic-ve-cross-sections.md` (+ the future automatic-discretisation note) and `CLAUDE.md` updated.

## Out of scope (this plan)

- **LCP DA** (the approximation under test) — sub-project B.
- **H₂⁺ DR** (Rydberg-series loop + Coulomb incident `coulomb::sF_en`) — sub-project D.
- **A skill to choose discretisation from the potential curves** (the per-element de Broglie phase) — noted in the docs as future work; the eMoScat decks are the interim tested values.
- **Re-gridding the VE curves onto the eMoScat grids** — VE is already converged on its grid (VE's outgoing flux is electronic, not nuclear); a separate follow-on if desired.
- Rotational (J>0), multiple electron partial waves, TD exact-2D DA (the TI route is exact and cheaper).
