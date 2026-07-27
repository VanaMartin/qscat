# LCP Dissociative Attachment (sub-project B) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the **local-complex-potential (LCP) approximation** of the dissociative-attachment cross section σ_DA(E) to `qscat.core` — the 1-D nuclear reduction (the approximation *under test*) — and validate it against the exact-2D DA oracle (`da_cross_section`) on F₂/NO (+ N₂ sanity), showing where the LCP holds and where it departs.

**Architecture:** The LCP replaces the full electron–nuclear problem with a 1-D nuclear wavepacket on a *local complex potential* `V_d(R) − iΓ(R)/2`, where `(V_d(R), Γ(R))` is the fixed-R electronic resonance pole. A doorway `d_{v₀}(R) = √(Γ(R)/2π)·χ_{v₀}` is propagated (order-3 Padé) under `H_res = T_nuc + diag(V_d − iΓ/2)`; DA is the outgoing flux at the dissociation boundary. Two model-independent additions to `qscat.core`: `local_complex_potential(model, …)` (the R-dependent `V_d/Γ` via `qscat.ecs.find_resonance_pole`, seeded from the anion bound state `anion_electronic_states` already gives at R_inf) and `lcp_da_cross_section(…)` (doorway propagation + boundary-flux energy transform). This is eMoScat's `ModelLCP` method (`SMatrix.cpp`); it is entirely 1-D nuclear (cheap — no 2-D solves).

**Tech Stack:** Python ≥3.12, NumPy/SciPy, `qscat.dvr` (FemDvrEcsGrid, `kinetic_sparse`, `eigen`), `qscat.ecs.find_resonance_pole`, `qscat.evolution.make_pade_stepper`, `qscat.core.dissociation` (`anion_electronic_states`, `da_cross_section` the oracle), `qscat.core.vibrational`, `qscat.model`. pytest, mypy --strict over `libs/qscat/qscat`, ruff.

## Global Constraints

- **Atomic units**; nuclear reduced mass `μ = model.mu`, electron mass 1.
- **c-product, never Hermitian** for every projection/correlation (`np.dot`, no conjugate) — under ECS `H = Hᵀ`. The DVR basis is `1/√w`-normalized, so the c-product is a plain coefficient dot.
- **`V_d(R) = Re(E_pole(R))`, `Γ(R) = max(0, −2·Im(E_pole(R)))`**, where `E_pole(R)` is the two-angle-matched resonance pole of the electronic Hamiltonian `−½∂²_r + model.surface(r, R)` at fixed `R`. Because `model.surface` **includes** `v0(R)`, `V_d` is `Re(E_pole)` directly (NOT `v0(R) + E_res` — that would double-count `v0`; the N₂ project's `vres.py` adds `v0` separately only because its `v_eff_el` excludes it). At `R → R_inf` the pole closes to the bound anion, so `V_d(R_inf) = ε_e` (the exact-DA threshold from `anion_electronic_states`) and `Γ(R_inf) → 0`.
- **DA threshold `ε_e = V_d(R_inf)`**; the DA channel is open only where `E_tot − ε_e > 0` (`E_tot = E + eps[v_init]`), `K = √(2μ(E_tot − ε_e))`.
- **σ prefactor `4·π³·|S|²/(2E)`** (identical to VE and the exact DA).
- **Propagator: order-3 diagonal Padé** (`make_pade_stepper(H_res, dt, 3)`, eMoScat's LCP setting), `H_res` **sparse** (`kinetic_sparse(grid, μ) + sp.diags(V_d − 0.5j·Γ)`).
- **`qscat.core` never imports `qscat.model`/`projects` at runtime** — `model: ResonanceModel` under `TYPE_CHECKING` only (enforced by `test_core_no_model_import.py`).
- **LCP is the approximation under test, the exact-2D `da_cross_section` is the ORACLE** — the deliverable is the *comparison* (where LCP agrees/departs), NOT tuning LCP to match. No independent DA data exists.

## File Structure

- `libs/qscat/qscat/core/lcp.py` (create) — `local_complex_potential`, `lcp_da_cross_section`.
- `libs/qscat/qscat/core/__init__.py` (modify) — export both.
- `libs/qscat/tests/test_lcp.py` (create) — Tasks 1–2 tests.
- `validation/diatomic/config.py` (modify) — per-molecule LCP electronic-grid angles + nuclear grid for the comparison.
- `validation/diatomic/lcp_da_curves.py` (create), `validation/diatomic/test_lcp_da_curves.py` (create) — Task 3 LCP-vs-exact comparison + gate.
- `docs/physics/diatomic-ve-cross-sections.md` (modify), `CLAUDE.md` (modify) — Task 4 docs.

---

### Task 1: `local_complex_potential` — model-independent V_d(R)/Γ(R)

**Files:**
- Create: `libs/qscat/qscat/core/lcp.py`
- Test: `libs/qscat/tests/test_lcp.py`

**Interfaces:**
- Consumes: `qscat.dvr.{FemDvrEcsGrid, eigen, kinetic}`, `qscat.ecs.find_resonance_pole`, `qscat.core.dissociation.anion_electronic_states`, `model.surface`/`model.v0`/`model.mu`.
- Produces:
  ```python
  def local_complex_potential(
      model: ResonanceModel,
      nuclear_grid: FemDvrEcsGrid,
      elec_grid_a: FemDvrEcsGrid,
      elec_grid_b: FemDvrEcsGrid,
      *,
      re_half_width: float = 0.05,
      im_half_width: float = 0.05,
      resid_tol: float = 1e-3,
  ) -> tuple[NDArray[complex128], NDArray[float64]]: ...   # (V_d, Gamma), shape (nuclear_grid.n,)
  ```
  `V_d[i]`, `Γ[i]` at nuclear grid point `nuclear_grid.points[i]`. `V_d` complex128 (evaluated at the possibly-ECS-rotated point), `Γ` float64 ≥ 0.

**Design notes (read before coding):** This is the N₂ `resonance_curve`/`vres_on_grid` continuation walk (`projects/n2_resonance/pole.py`, `projects/n2_ti_cross_section/vres.py`) made **model-independent** by (a) building `H_el(R) = kinetic(elec_grid, 1.0) + np.diag(model.surface(elec_grid.points, R))` instead of the N₂ `v_eff_el`, and (b) **seeding the continuation from the anion bound state at R_inf** rather than a per-molecule equilibrium window. Steps:
1. `R_inf = nuclear_grid.R0`. Get the seed pole: `eps_e, _ = anion_electronic_states(elec_grid_a, model, R_inf, 1)` → the bound anion energy `ε_e` (real). Seed `window = (ε_e − re_half_width, ε_e + re_half_width, −im_half_width, im_half_width)`.
2. Sort the **real** nuclear grid points descending in R (start at the outer edge nearest R_inf, where the pole ≈ bound anion). Walk **inward** (decreasing R): at each `R`, `E_a, _ = eigen(H_el(R, elec_grid_a))`, `E_b, _ = eigen(H_el(R, elec_grid_b))`, `E_pole, resid = find_resonance_pole(E_a, E_b, window)`. Accept if `resid < resid_tol`: set `V_d = E_pole.real`, `Γ = max(0, −2·E_pole.imag)`, recenter `window` on `E_pole` (±half-widths). On the FIRST rejected step (raise or `resid ≥ tol`), stop advancing and **freeze** `E_pole`'s resonance-shift at the last accepted value for all remaining (smaller-R) real points — small-R breakdown, exactly as `vres.py` documents (physically irrelevant: `v0(R)` is Ha above threshold there, χ negligible). Catch `(ValueError, np.linalg.LinAlgError)` as a rejected step.
3. `V_d[i] = ` the accepted/frozen `V_d` for real points; `Γ[i]` likewise (frozen region keeps its last Γ). **Frozen-region caveat:** freeze the *shift* `s = V_d_lastgood − v0(R_lastgood)`, and set `V_d[i] = v0(R_i) + s` at each frozen point's own `R_i` (so the rapidly-rising `v0` is still tracked; only the electronic shift is frozen). Store this consistently for the accepted region too (`s_i = V_d_i − v0(R_i)`), so a single formula `V_d[i] = v0(R_i) + s_i` holds everywhere.
4. **Complex tail** (`R.imag ≠ 0`, i.e. `R > R0`): clamp `Γ = 0`, `V_d = v0(R_complex) + s_asymptote`, `s_asymptote = s` at the largest real R (`= ε_e − v0(R_inf)`). Mirrors `vres.py`.

Build the two electronic grids at two ECS angles (e.g. 35° and 44°, N₂'s) via `qscat.core.grids.electronic_grid(r_max=…, angle_deg=…, …)` — the caller passes them in (per-molecule extents live in config, Task 3).

- [ ] **Step 1: Write the failing test**

```python
# libs/qscat/tests/test_lcp.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.lcp import local_complex_potential
from qscat.core.dissociation import anion_electronic_states
from qscat.model import F2, N2


def _elec_grids():
    return (electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=35.0),
            electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=44.0))


def test_vd_gamma_shapes_and_gamma_nonneg():
    g_R = nuclear_grid(r_max=22.0, n_complex=6, quadrature=10)
    ga, gb = _elec_grids()
    Vd, Gamma = local_complex_potential(F2, g_R, ga, gb)
    assert Vd.shape == (g_R.n,) and Gamma.shape == (g_R.n,)
    assert Vd.dtype == np.complex128 and Gamma.dtype == np.float64
    assert np.all(Gamma >= 0.0)


def test_gamma_closes_and_vd_matches_anion_at_large_R():
    # As R -> R_inf the pole closes to the bound anion: Gamma -> ~0 and
    # V_d(R_inf) == eps_e (the exact-DA threshold from anion_electronic_states).
    g_R = nuclear_grid(r_max=22.0, n_complex=6, quadrature=10)
    ga, gb = _elec_grids()
    Vd, Gamma = local_complex_potential(F2, g_R, ga, gb)
    R = g_R.points
    real = R.imag == 0.0
    i_outer = np.flatnonzero(real)[np.argmax(R[real].real)]  # largest real R
    eps_e, _ = anion_electronic_states(ga, F2, g_R.R0, 1)
    assert Gamma[i_outer] < 1e-3                               # closed at the edge
    assert abs(Vd[i_outer].real - eps_e[0]) < 5e-3            # == anion asymptote


def test_gamma_positive_in_resonance_region():
    # At smaller R (inside the crossing) the anion is a real resonance: Gamma>0.
    g_R = nuclear_grid(r_max=22.0, n_complex=6, quadrature=10)
    ga, gb = _elec_grids()
    Vd, Gamma = local_complex_potential(F2, g_R, ga, gb)
    R = g_R.points.real
    band = (R > 1.5) & (R < 2.5)
    assert Gamma[band].max() > 1e-4                           # genuine width somewhere
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest libs/qscat/tests/test_lcp.py -q`
Expected: FAIL — `ModuleNotFoundError: qscat.core.lcp`.

- [ ] **Step 3: Implement `local_complex_potential`**

Create `libs/qscat/qscat/core/lcp.py` with the module docstring (cite eMoScat `ModelLCP`, the research-program "approximation under test" framing, and the `V_d = Re(E_pole)` no-double-count note), the `TYPE_CHECKING` `ResonanceModel` import, and:

```python
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, eigen, kinetic
from qscat.ecs import find_resonance_pole

from .dissociation import anion_electronic_states

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["local_complex_potential"]  # Task 2 appends "lcp_da_cross_section"


def _h_el(model: ResonanceModel, R: complex, g: FemDvrEcsGrid) -> npt.NDArray[np.complex128]:
    return kinetic(g, 1.0) + np.diag(model.surface(g.points, R))


def local_complex_potential(
    model: ResonanceModel,
    nuclear_grid: FemDvrEcsGrid,
    elec_grid_a: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    *,
    re_half_width: float = 0.05,
    im_half_width: float = 0.05,
    resid_tol: float = 1e-3,
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64]]:
    """Local complex potential `(V_d(R), Gamma(R))` for the LCP DA/VE model.

    `V_d(R) = Re(E_pole(R))`, `Gamma(R) = max(0, -2 Im(E_pole(R)))`, `E_pole`
    the two-angle-matched resonance pole of `-1/2 d^2/dr^2 + model.surface(r,R)`
    (surface includes v0(R), so V_d = Re(E_pole) directly). Seeded from the
    bound anion at `R_inf = nuclear_grid.R0` (`anion_electronic_states`) and
    continued INWARD; small-R breakdown freezes the electronic shift; the
    complex tail clamps Gamma=0. See module docstring.
    """
    R_inf = nuclear_grid.R0
    eps_e, _ = anion_electronic_states(elec_grid_a, model, R_inf, 1)
    window = (eps_e[0] - re_half_width, eps_e[0] + re_half_width, -im_half_width, im_half_width)

    pts = nuclear_grid.points
    real_idx = np.flatnonzero(pts.imag == 0.0)
    order = np.argsort(pts[real_idx].real)[::-1]          # descending R: outer -> inner
    walk = real_idx[order]
    R_real = pts[walk].real

    shift = np.empty(walk.size, dtype=np.float64)          # s = V_d - v0(R)
    gamma_w = np.empty(walk.size, dtype=np.float64)
    last_s: float | None = None
    last_g = 0.0
    broken = False
    for j, gidx in enumerate(walk):
        R = float(R_real[j])
        if not broken:
            try:
                E_pole, resid = find_resonance_pole(
                    eigen(_h_el(model, R, elec_grid_a))[0],
                    eigen(_h_el(model, R, elec_grid_b))[0],
                    window,
                )
            except (ValueError, np.linalg.LinAlgError):
                resid = np.inf
            else:
                if resid < resid_tol:
                    v0R = float(np.real(model.v0(np.asarray(R))))
                    last_s = E_pole.real - v0R
                    last_g = max(0.0, -2.0 * E_pole.imag)
                    window = (E_pole.real - re_half_width, E_pole.real + re_half_width,
                              E_pole.imag - im_half_width, E_pole.imag + im_half_width)
                    shift[j], gamma_w[j] = last_s, last_g
                    continue
            broken = True
        if last_s is None:
            raise RuntimeError("local_complex_potential: pole finder failed at the seed edge")
        shift[j], gamma_w[j] = last_s, last_g

    Vd = np.empty(nuclear_grid.n, dtype=np.complex128)
    Gamma = np.zeros(nuclear_grid.n, dtype=np.float64)
    Vd[walk] = model.v0(R_real) + shift
    Gamma[walk] = gamma_w

    tail = np.flatnonzero(pts.imag != 0.0)
    if tail.size:
        assert last_s is not None
        s_asym = shift[0]                                  # shift at the largest real R
        Vd[tail] = model.v0(pts[tail]) + s_asym
    return Vd, Gamma
```

(Note: `shift[0]` is the outermost real R — the walk starts there — so `s_asym` is the asymptotic electronic binding `ε_e − v0(R_inf)`.)

Add `from collections.abc import Sequence`-style imports only if needed; keep `__all__` to `["local_complex_potential"]` (Task 2 grows it).

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest libs/qscat/tests/test_lcp.py -q`
Expected: PASS (3 tests). These do ~1 electronic diagonalization pair per real nuclear point (~hundreds); on the small grid it runs in seconds-to-tens-of-seconds.

- [ ] **Step 5: Differential check against the N₂ project oracle**

Add one more test comparing to the existing validated N₂ `vres_on_grid` on the N₂ nuclear grid (the differential oracle — same physics, different code path). The two use different electronic-shift bookkeeping (`vres.py`: `v0 + E_res` with `v_eff_el` excluding v0; here: `Re(E_pole)` with `surface` including v0), so compare the **observable** `V_d(R)` and `Γ(R)` themselves, which must agree:

```python
def test_matches_n2_vres_oracle():
    from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
    from projects.n2_ti_cross_section.vres import vres_on_grid
    g_R = n2_nuclear_grid()
    ga, gb = _elec_grids()
    Vd, Gamma = local_complex_potential(N2, g_R, ga, gb)
    Vd_ref, Gamma_ref = vres_on_grid(g_R)
    real = g_R.points.imag == 0.0
    # compare on the resonance region where both are well-defined (R in [1.5,3.5])
    R = g_R.points.real
    band = real & (R > 1.5) & (R < 3.5)
    assert np.allclose(Vd[band].real, Vd_ref[band].real, atol=5e-3)
    assert np.allclose(Gamma[band], Gamma_ref[band], atol=5e-3)
```

Run: `uv run pytest libs/qscat/tests/test_lcp.py::test_matches_n2_vres_oracle -q`
Expected: PASS. (This test imports `projects.*` — that is fine for a LIBRARY TEST, which may cross-check against a project oracle; the `qscat.core` PACKAGE still imports no project. If it reveals a real bookkeeping discrepancy, fix the code, not the tolerance — the two V_d/Γ observables must agree.)

- [ ] **Step 6: Type-check + lint + commit**

Run: `uv run mypy libs/qscat/qscat/core && uv run ruff check libs/qscat/qscat/core libs/qscat/tests/test_lcp.py`
```bash
git add libs/qscat/qscat/core/lcp.py libs/qscat/tests/test_lcp.py
git commit -m "feat(lcp): model-independent local complex potential V_d(R)/Gamma(R)"
```

---

### Task 2: `lcp_da_cross_section` — doorway propagation + boundary flux

**Files:**
- Modify: `libs/qscat/qscat/core/lcp.py`
- Modify: `libs/qscat/qscat/core/__init__.py`
- Test: `libs/qscat/tests/test_lcp.py`

**Interfaces:**
- Consumes: `qscat.dvr.kinetic_sparse`, `qscat.evolution.make_pade_stepper`, `scipy.sparse`.
- Produces:
  ```python
  def lcp_da_cross_section(
      nuclear_grid: FemDvrEcsGrid, mu: float,
      Vd: NDArray[complex128], Gamma: NDArray[float64],
      eps: NDArray[float64], chi: NDArray[complex128],
      v_init: int, eps_e: float, E: float | ArrayLike, *,
      dt: float = 1.0, n_steps: int = 1500, order: int = 3,
  ) -> NDArray[float64]: ...
  ```
  σ_DA(E) (bohr²). Scalar `E` → shape `()` (0-d) or `(1,)`; array `E` → `(len(E),)`. `σ=0` where `E ≤ 0` or `E_tot − eps_e ≤ 0` (closed).

**Design notes:** mirror `projects/n2_td_cross_section/td_cross_section.py`, but (a) the propagated observable is `ψ(X, t_n)` — the wavefunction coefficient at the **dissociation-boundary grid point** `X` (the outermost REAL point, `boundary_idx = argmax(points.real where imag==0)`, `X = points[boundary_idx].real`), recorded each step; (b) the energy transform is the DA one: `S_DA(E) = √(K/2πμ)·e^{−iKX}·Σ_n w_n·e^{i(E+eps[v_init])t_n}·ψ(X,t_n)·dt`, `K = √(2μ(E_tot − eps_e))`, `σ_DA = 4π³|S_DA|²/2E`. Doorway `d(R) = √(Γ(R)/2π)·χ_{v_init}`. Propagator: `H_res = kinetic_sparse(grid, mu) + sp.diags(Vd − 0.5j·Gamma)` (SPARSE), `step = make_pade_stepper(H_res, dt, order)`. Reuse the Simpson/trapezoidal `_quadrature_weights` (copy from `td_cross_section.py` — a small helper). Record `psi_X[n]` for `n = 0..n_steps` (t_0=0 included), stepping between records.

- [ ] **Step 1: Write the failing test** (append to `test_lcp.py`)

```python
from qscat.core.lcp import lcp_da_cross_section
from qscat.core.vibrational import vibrational_states


def _f2_lcp_inputs():
    g_R = nuclear_grid(r_max=22.0, n_complex=6, quadrature=10)
    ga, gb = _elec_grids()
    Vd, Gamma = local_complex_potential(F2, g_R, ga, gb)
    eps, chi = vibrational_states(g_R, F2.mu, 3, F2.v0)
    eps_e = float(Vd[np.flatnonzero(g_R.points.imag == 0.0)][np.argmax(g_R.points.real[g_R.points.imag == 0.0])].real)
    return g_R, Vd, Gamma, eps, chi, eps_e


def test_lcp_da_shape_and_nonneg():
    g_R, Vd, Gamma, eps, chi, eps_e = _f2_lcp_inputs()
    E = np.array([0.02, 0.03, 0.04])
    s = lcp_da_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, eps_e, E, dt=1.0, n_steps=400)
    assert s.shape == (3,)
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)


def test_lcp_da_closed_channel_is_zero():
    # F2 is exothermic (no closed region for E>0), so test the closed guard
    # directly: a threshold eps_e above E_tot forces E_DR<0 -> sigma must be 0.
    g_R, Vd, Gamma, eps, chi, _ = _f2_lcp_inputs()
    eps_e_high = float(eps[0]) + 1.0
    s = lcp_da_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, eps_e_high,
                             np.array([0.03]), dt=1.0, n_steps=400)
    assert s[0] == 0.0


@pytest.mark.slow
def test_lcp_da_f2_positive():
    g_R, Vd, Gamma, eps, chi, eps_e = _f2_lcp_inputs()
    E = np.array([0.02, 0.03, 0.04])
    s = lcp_da_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, eps_e, E, dt=1.0, n_steps=2000)
    assert s.max() > 0.0                     # exothermic F2: DA is open and nonzero
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest libs/qscat/tests/test_lcp.py -q -k lcp_da`
Expected: FAIL — `cannot import name 'lcp_da_cross_section'`.

- [ ] **Step 3: Implement `lcp_da_cross_section`** (append to `lcp.py`)

```python
import scipy.sparse as sp
from qscat.dvr import kinetic_sparse
from qscat.evolution import make_pade_stepper


def _quadrature_weights(n_t: int) -> npt.NDArray[np.float64]:
    """Composite Simpson weights (unscaled by dt); trapezoidal fallback for even n_t.
    Copied from projects/n2_td_cross_section/td_cross_section.py."""
    if n_t < 2:
        raise ValueError("need at least 2 time samples")
    if n_t % 2 == 1:
        w = np.ones(n_t); w[1:-1:2] = 4.0; w[2:-1:2] = 2.0; w /= 3.0
    else:
        w = np.full(n_t, 2.0); w[0] = w[-1] = 1.0; w /= 2.0
    return np.asarray(w, dtype=np.float64)


def lcp_da_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    eps_e: float,
    E: float | npt.ArrayLike,
    *,
    dt: float = 1.0,
    n_steps: int = 1500,
    order: int = 3,
) -> npt.NDArray[np.float64]:
    """LCP dissociative-attachment sigma_DA(E) (bohr^2): the 1-D nuclear
    doorway propagated on V_d - i Gamma/2, DA = flux at the dissociation
    boundary. eMoScat ModelLCP/SMatrix.cpp. The approximation under test vs the
    exact-2D `da_cross_section` oracle."""
    pts = nuclear_grid.points
    real_idx = np.flatnonzero(pts.imag == 0.0)
    b = int(real_idx[np.argmax(pts[real_idx].real)])   # outermost real point index
    X = float(pts[b].real)

    doorway = np.sqrt(Gamma / (2.0 * np.pi)).astype(np.complex128) * chi[v_init]
    H_res = kinetic_sparse(nuclear_grid, mu) + sp.diags(Vd - 0.5j * Gamma)
    step = make_pade_stepper(H_res.tocsc(), dt, order)

    n_t = n_steps + 1
    t = np.arange(n_t, dtype=np.float64) * dt
    psi_X = np.empty(n_t, dtype=np.complex128)
    psi = doorway.copy()
    for n in range(n_t):
        psi_X[n] = psi[b]
        if n < n_steps:
            psi = step(psi)

    w = _quadrature_weights(n_t)
    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    out = np.zeros(e_arr.size, dtype=np.float64)
    for ie, e in enumerate(e_arr):
        if e <= 0.0:
            continue
        e_tot = float(e) + eps[v_init]
        e_dr = e_tot - eps_e
        if e_dr <= 0.0:
            continue
        K = float(np.sqrt(2.0 * mu * e_dr))
        phase = np.exp(1j * e_tot * t)
        S = np.sqrt(K / (2.0 * np.pi * mu)) * np.exp(-1j * K * X) * np.sum(w * phase * psi_X) * dt
        out[ie] = 4.0 * np.pi**3 * abs(S) ** 2 / (2.0 * float(e))

    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    return np.asarray(out[0] if scalar else out, dtype=np.float64)
```

Append `"lcp_da_cross_section"` to `lcp.py`'s `__all__`. Export both names from `libs/qscat/qscat/core/__init__.py` (import + `__all__` + Public-API docstring entries).

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest libs/qscat/tests/test_lcp.py -q` (non-slow), then `uv run pytest libs/qscat/tests/test_lcp.py -q -m slow`.
Expected: all PASS.

- [ ] **Step 5: Boundary guard + type-check + lint**

Run: `uv run pytest libs/qscat/tests/test_core_no_model_import.py -q && uv run mypy libs/qscat/qscat && uv run ruff check libs/qscat`
Expected: import-boundary test PASSES (lcp.py imports no `qscat.model` at runtime), mypy 0, ruff clean.

- [ ] **Step 6: Commit**

```bash
git add libs/qscat/qscat/core/lcp.py libs/qscat/qscat/core/__init__.py libs/qscat/tests/test_lcp.py
git commit -m "feat(lcp): lcp_da_cross_section (doorway propagation + boundary-flux transform)"
```

---

### Task 3: Per-molecule LCP config + LCP-vs-exact comparison

**Files:**
- Modify: `validation/diatomic/config.py`
- Create: `validation/diatomic/lcp_da_curves.py`
- Test: `validation/diatomic/test_lcp_da_curves.py`

**Interfaces:**
- `MoleculeConfig` gains LCP electronic-angle-grid fields + builders:
  ```python
  lcp_angle_a: float = 35.0
  lcp_angle_b: float = 44.0
  def lcp_elec_grids(self) -> tuple[FemDvrEcsGrid, FemDvrEcsGrid]: ...   # electronic_grid at the two angles (r_max=e_r_max, order=e_order, n_complex=e_n_complex)
  def lcp_nuclear_grid(self) -> FemDvrEcsGrid: ...                        # the eMoScat da nuclear grid (reuse the da_grid()'s nuclear factor)
  ```
- `validation/diatomic/lcp_da_curves.py`: `compute_lcp_da_curve(cfg, E_grid, *, dt=1.0, n_steps=1500) -> (E, sigma_lcp)` and `main()` that overlays σ_DA^LCP vs σ_DA^exact (from `validation.diatomic.da_curves.compute_da_curve`) into `docs/physics/figures/{f2,no}-2d-da-lcp-vs-exact.png`.

**Design notes:** the LCP nuclear grid should be the SAME per-molecule eMoScat nuclear grid used for the exact DA (`da_grid().grids[1]`) so the comparison is grid-consistent — expose it as `lcp_nuclear_grid()` (build the nuclear factor via `segmented_grid(nuc_real, nuc_complex, angle_deg=nuc_angle, quadrature=nuc_quad)`). `eps_e` comes from `local_complex_potential`'s `V_d` at the outer real point (== the anion asymptote). `compute_lcp_da_curve`: build `Vd, Gamma = local_complex_potential(cfg.model, g_R, *cfg.lcp_elec_grids())`, `eps, chi = vibrational_states(g_R, cfg.model.mu, cfg.n_vib, cfg.model.v0)`, then `lcp_da_cross_section(...)`. The LCP is 1-D (cheap) EXCEPT `local_complex_potential` does ~2 electronic diagonalizations per nuclear real point (~hundreds) — a few minutes; the propagation itself is fast.

- [ ] **Step 1: Write the failing gate test**

```python
# validation/diatomic/test_lcp_da_curves.py
from __future__ import annotations

import numpy as np
import pytest
from validation.diatomic.config import CONFIGS
from validation.diatomic.lcp_da_curves import compute_lcp_da_curve


@pytest.mark.slow
def test_f2_lcp_da_positive_and_finite():
    E, s = compute_lcp_da_curve(CONFIGS["F2"], np.array([0.02, 0.03, 0.04]), n_steps=1500)
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)
    assert s.max() > 0.0


@pytest.mark.slow
def test_f2_lcp_within_order_of_magnitude_of_exact():
    # The scientific check: the LCP approximation should be the SAME ORDER as
    # the exact oracle across F2's DA window (a cross-model band, NOT exact
    # agreement -- the LCP is the approximation under test).
    from validation.diatomic.da_curves import compute_da_curve
    E = np.array([0.02, 0.03, 0.04])
    _, s_lcp = compute_lcp_da_curve(CONFIGS["F2"], E, n_steps=1500)
    _, s_exact = compute_da_curve(CONFIGS["F2"], E)
    ratio = s_lcp / s_exact[:, 0]
    assert np.all((ratio > 0.1) & (ratio < 10.0))   # within 10x either way
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest validation/diatomic/test_lcp_da_curves.py -q`
Expected: FAIL — `ModuleNotFoundError`.

- [ ] **Step 3: Extend config + implement `lcp_da_curves.py`**

Add the `lcp_*` fields/methods to `MoleculeConfig` (verify the actual field names for the nuclear deck — `nuc_real`/`nuc_complex`/`nuc_angle`/`nuc_quad` from sub-project A Task 6 — and reuse `segmented_grid`). Implement `compute_lcp_da_curve` + `main()` per the Interfaces. `main()` overlays LCP vs exact via `plot_cross_sections` (pass the exact curve as `reference=(E, sigma_exact)`), writing the comparison figure. Read `validation/diatomic/config.py`, `validation/diatomic/da_curves.py`, and `qscat/core/plot.py` first and match the real signatures.

- [ ] **Step 4: Run the gate**

Run: `uv run pytest validation/diatomic/test_lcp_da_curves.py -q -m slow` (heavy — the V_d/Γ pole search + propagation; minutes).
Expected: PASS — LCP σ_DA finite/positive and within ~10× of the exact oracle across F₂'s window.

- [ ] **Step 5: Type-check + lint + commit**

Run: `uv run mypy libs/qscat/qscat && uv run ruff check .`
```bash
git add validation/diatomic/config.py validation/diatomic/lcp_da_curves.py validation/diatomic/test_lcp_da_curves.py
git commit -m "feat(lcp): per-molecule LCP config + F2/NO LCP-vs-exact DA comparison"
```

---

### Task 4: Comparison figures, docs, and CLAUDE.md

**Files:**
- Create: `docs/physics/figures/{f2,no}-2d-da-lcp-vs-exact.png` (generated)
- Modify: `docs/physics/diatomic-ve-cross-sections.md`
- Modify: `CLAUDE.md`

**Design notes:** generate the F₂ and NO LCP-vs-exact overlays via `validation.diatomic.lcp_da_curves.main()` (heavy — controller/background). N₂ is a closed-channel sanity (both σ_DA ≈ 0) — mention, no figure needed. The doc section states the scientific finding: where the LCP tracks the exact DA and where it departs (the LCP's fixed-R local-complex-potential ansatz vs the exact non-adiabatic 2-D dynamics), quantified by the committed overlays.

- [ ] **Step 1: Generate the figures**

Run: `uv run python -m validation.diatomic.lcp_da_curves` (heavy, minutes — the pole search dominates). Expected: two overlay PNGs written.

- [ ] **Step 2: Write the doc section**

In `docs/physics/diatomic-ve-cross-sections.md`, add an "LCP DA vs exact-2D DA" subsection: the LCP method (1-D doorway on `V_d − iΓ/2`, boundary flux), the model-independent `V_d/Γ` via `find_resonance_pole` seeded from the anion state, the committed overlays, and the quantitative agreement/departure (the scientific point — the exact solver is the oracle, the LCP is under test). Note N₂'s closed-channel sanity.

- [ ] **Step 3: Update CLAUDE.md**

Add a one-line `qscat.core.lcp` entry (`local_complex_potential`, `lcp_da_cross_section` — the LCP DA approximation) to the `qscat.core` bullet, and note `validation/diatomic/lcp_da_curves.py`.

- [ ] **Step 4: Commit**

```bash
git add docs/physics/figures/f2-2d-da-lcp-vs-exact.png docs/physics/figures/no-2d-da-lcp-vs-exact.png \
        docs/physics/diatomic-ve-cross-sections.md CLAUDE.md
git commit -m "docs(lcp): F2/NO LCP-vs-exact DA comparison figures + docs + CLAUDE.md"
```

---

## Verification (whole sub-project)

- `uv run pytest -q -m "not slow"` passes; the `@slow` LCP tests pass.
- `uv run pytest libs/qscat/tests/test_core_no_model_import.py -q` passes (LCP code keeps the core/model boundary).
- `uv run mypy libs/qscat/qscat` 0 errors; `uv run ruff check .` clean.
- `local_complex_potential` matches the N₂ `vres_on_grid` oracle (V_d/Γ) in the resonance region; `Γ→0`, `V_d(R_inf)=ε_e` at the edge.
- `lcp_da_cross_section` respects the DA threshold, is finite/≥0, and (F₂) positive; the LCP-vs-exact comparison is within a documented cross-model band. Committed overlays + docs + CLAUDE.md updated.

## Out of scope (this plan)

- **Promoting the LCP VE** (`lcp_ve_cross_section`, the correlation observable) — `local_complex_potential` enables it, but B is DA-focused; a natural follow-on.
- **Tuning the LCP to match the exact oracle** — forbidden by the research program; the deliverable is the comparison.
- **H₂⁺ DR** (Coulomb) — sub-project D.
- Rotational (J>0), multiple partial waves.
