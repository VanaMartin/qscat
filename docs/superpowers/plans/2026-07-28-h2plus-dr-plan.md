# H₂⁺ ionic model + dissociative recombination (sub-project D) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the H₂⁺ **ionic** model (Coulomb special functions, the Coulomb-tail potential, the large discretization) into `qscat`, and add the **dissociative-recombination (DR)** cross section — validated analytically + on a small laptop proxy, with the full 1300-bohr deck Docker/MUMPS-ready.

**Architecture:** H₂⁺ is the first ion: the electron sees a −1/r Coulomb tail, so channel functions become **Coulomb** (charge z=−1), the process is DR (`e⁻+H₂⁺(v)→H+H` via a Rydberg resonance dissociating on the neutral curve), and there is a **Rydberg series** of exit channels. Additions: `qscat.special.coulomb` (energy-normalized Coulomb F/G/H via mpmath, the charge-z generalization of `riccati_bessel_en`); a narrowed `ResonanceModel` protocol with a `charge` attribute + a new **ionic** model form + `H2P` registry entry; a charge-aware `channel_vector`; and `dr_cross_section` — `da_cross_section` generalized to a Coulomb incident + a Rydberg-channel loop. See `docs/superpowers/specs/2026-07-28-h2plus-dr-design.md`.

**Tech Stack:** Python ≥3.12, NumPy/SciPy, **mpmath** (new qscat dep — complex-arg Coulomb functions), `qscat.dvr`, `qscat.linalg`, `qscat.special`, `qscat.core` (`da_cross_section`, `anion_electronic_states`, `v_dr_diag`, `riccati_bessel_en_mass`, `channel_vector`, `segmented_grid`), `qscat.model`. pytest, mypy --strict over `libs/qscat/qscat`, ruff.

## Global Constraints

- **Atomic units**; electron mass 1, nuclear reduced mass μ = model.mu (918.25 for H₂⁺), electronic partial wave ℓ = model.ell (1 for H₂⁺), Coulomb charge z = model.charge (−1 for H₂⁺, 0 for neutrals).
- **c-product, never Hermitian** (`qscat.linalg.c_product`, no conjugate) for every projection. eMoScat uses a conjugated dot (`zdotc`) in the DR T-matrix, but the ECS-correct choice is the c-product (as everywhere in qscat) — a validation check confirms the rotated-tail difference is negligible.
- **Coulomb functions:** `η = m·z/k`, `ρ = k·x`, energy-normalization `√(2m/(πk))`. `H⁺ = G + iF` (outgoing), `H⁻ = G − iF`. Do NOT replicate eMoScat's `sH1` bug (it returned F, not G+iF). At z=0, `F_en` MUST equal `riccati_bessel_en` (mass m).
- **`V_DR = V_int(r,R) + v0(R) − V_int(r, R_inf)`** (the rearrangement interaction, → 0 as R→∞), `R_inf = tgrid.grids[1].R0`. Same form as the DA `v_dr_diag`.
- **σ prefactor `4·π³·|T|²/(2E)`** (= `4π³|T|²/k_el²`), identical to VE/DA.
- **Potentials evaluated on the COMPLEX ECS coordinate** (`tgrid.points()`), never the real part — fixes eMoScat's `// FIXME`.
- **`qscat.core` never imports `qscat.model`/`projects` at runtime** — `model: ResonanceModel` under `TYPE_CHECKING` (enforced by `test_core_no_model_import.py`).
- **Full size is non-laptop** — tests use a reduced proxy grid; the full deck is Docker/MUMPS-ready but not run in the suite.

## File Structure

- `libs/qscat/pyproject.toml` (modify) — add `mpmath` to dependencies.
- `libs/qscat/qscat/special/coulomb.py` (create), `special/__init__.py` (modify) — Coulomb functions.
- `libs/qscat/qscat/model/diatomic.py` (modify) — narrow protocol (drop `lam`, add `charge`); `DiatomicResonanceModel.charge=0`.
- `libs/qscat/qscat/model/ionic.py` (create), `model/library.py` (modify), `model/__init__.py` (modify) — the ionic model + `H2P`.
- `libs/qscat/qscat/core/channels.py` (modify) — charge-aware `channel_vector`.
- `libs/qscat/qscat/core/dissociation.py` (modify) — `dr_cross_section` (+ export in `core/__init__.py`).
- `libs/qscat/qscat/core/grids.py` (modify) — an exp-growth ECS-tail nuclear/electronic builder if the existing ones don't cover the H₂⁺ deck.
- `validation/h2plus/` (create) — `config.py` (H₂⁺ deck + reduced proxy), `dr.py` (the DR driver), `test_dr.py`.
- Tests: `libs/qscat/tests/test_coulomb.py`, `test_ionic_model.py`; extend `test_channels`/`test_dissociation`.
- `docs/physics/h2plus-dr.md` (create), `CLAUDE.md` (modify).

---

### Task 1: Coulomb special functions (`qscat.special.coulomb`)

**Files:**
- Modify: `libs/qscat/pyproject.toml` (add `mpmath` dep)
- Create: `libs/qscat/qscat/special/coulomb.py`
- Modify: `libs/qscat/qscat/special/__init__.py`
- Test: `libs/qscat/tests/test_coulomb.py`

**Interfaces:**
- Produces: `coulomb_f_en(x, k, z, m, l)`, `coulomb_g_en(...)`, `coulomb_h1_en(...)` — energy-normalized regular / irregular / outgoing Coulomb functions, accepting REAL or COMPLEX `x` arrays, returning `complex128`. `coulomb_f_en(x, k, 0, m, l)` equals `riccati_bessel_en(x, k, l)` at m=1 (and its mass-m generalization).

**Design notes:** `mpmath.coulombf(l, eta, rho)` / `coulombg(l, eta, rho)` accept complex `rho`; loop over the (small-per-call, but vectorize with a comprehension) points, `eta = m*z/k`, `rho = k*x`, prefactor `√(2m/(πk))`. `H⁺ = G + iF`. Cast each `mpmath.mpc` to `complex`. Guard `k>0`, `m>0`. This is slow (mpmath) — acceptable (H₂⁺ isn't laptop-scale; optimize later). Add `mpmath` to `[project.dependencies]` in `libs/qscat/pyproject.toml` and `uv sync --all-packages`.

- [ ] **Step 1: Add mpmath + write the failing test**

Add `"mpmath>=1.3"` to `dependencies` in `libs/qscat/pyproject.toml`; run `uv sync --all-packages`.

```python
# libs/qscat/tests/test_coulomb.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.special import coulomb_f_en, coulomb_h1_en, riccati_bessel_en, riccati_hankel_en


def test_f_en_reduces_to_riccati_bessel_at_zero_charge():
    r = np.linspace(0.5, 40.0, 60)
    for l in (0, 1, 2):
        got = coulomb_f_en(r, 1.3, 0.0, 1.0, l)
        assert np.allclose(got.real, riccati_bessel_en(r, 1.3, l), rtol=1e-8, atol=1e-9)
        assert np.allclose(got.imag, 0.0, atol=1e-9)


def test_h1_en_reduces_to_riccati_hankel_at_zero_charge():
    r = np.linspace(0.5, 30.0, 40)
    got = coulomb_h1_en(r, 1.0, 0.0, 1.0, 1)      # G + iF -> Riccati-Hankel h1
    assert np.allclose(got, riccati_hankel_en(r, 1.0, 1), rtol=1e-7, atol=1e-8)


def test_attractive_coulomb_known_value():
    # mpmath.coulombf(1, -0.5, 2.0) = 0.972687664241193 (eta = m z/k = -0.5)
    # coulomb_f_en(x, k, z, m, 1) with k*x=2, m z/k=-0.5: pick k=1, x=2, z=-0.5, m=1
    got = coulomb_f_en(np.array([2.0]), 1.0, -0.5, 1.0, 1)
    expect = np.sqrt(2.0 / np.pi) * 0.972687664241193
    assert abs(got[0] - expect) < 1e-9


def test_accepts_complex_ecs_argument():
    r = np.array([3.0 + 0.4j, 10.0 + 2.0j])       # ECS-rotated points
    got = coulomb_f_en(r, 1.0, -1.0, 1.0, 1)
    assert got.shape == (2,) and np.all(np.isfinite(got))


@pytest.mark.parametrize("bad", [(0.0, 1.0), (1.0, 0.0)])
def test_rejects_nonpositive_k_or_m(bad):
    with pytest.raises(ValueError):
        coulomb_f_en(np.array([1.0]), bad[0], -1.0, bad[1], 1)
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest libs/qscat/tests/test_coulomb.py -q`
Expected: FAIL — `cannot import name 'coulomb_f_en'`.

- [ ] **Step 3: Implement `coulomb.py`**

```python
# libs/qscat/qscat/special/coulomb.py
"""Energy-normalized Coulomb functions (the charge-z generalization of the
Riccati-Bessel/Hankel radial functions in `radial.py`).

For a scattering particle of mass `m` in a Coulomb field of charge `z` (the
Sommerfeld parameter `eta = m z / k`, momentum `k`), the energy-normalized
regular / irregular / outgoing radial solutions are

    F_en(x) = sqrt(2 m/(pi k)) F_l(eta, k x)
    G_en(x) = sqrt(2 m/(pi k)) G_l(eta, k x)
    H1_en(x) = sqrt(2 m/(pi k)) (G_l + i F_l)(eta, k x)   [outgoing, H+ = G + iF]

with F_l/G_l the standard regular/irregular Coulomb functions (mpmath, which
accepts COMPLEX arguments -- needed for ECS-rotated x). At z=0 (eta=0),
F_l(0, rho) = rho j_l(rho), so `coulomb_f_en(., ., 0, m, l)` reduces to
`riccati_bessel_en` at mass m -- the differential-test hook. eMoScat's
`sH1` wrapper had a copy-paste bug (returned F, not G+iF); we define H+ = G + iF
correctly. (eMoScat coulomb.cpp / coulcc.f; the DR incident wave uses F_en.)
"""

from __future__ import annotations

import mpmath
import numpy as np
import numpy.typing as npt

__all__ = ["coulomb_f_en", "coulomb_g_en", "coulomb_h1_en"]


def _fg(x: npt.NDArray[np.complex128], k: float, z: float, m: float, l: int
        ) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
    if k <= 0.0:
        raise ValueError(f"k must be positive, got {k}")
    if m <= 0.0:
        raise ValueError(f"m must be positive, got {m}")
    eta = m * z / k
    xs = np.asarray(x, dtype=np.complex128).ravel()
    pref = np.sqrt(2.0 * m / (np.pi * k))
    f = np.empty(xs.size, dtype=np.complex128)
    g = np.empty(xs.size, dtype=np.complex128)
    for i, xv in enumerate(xs):
        rho = mpmath.mpc(k * xv)
        f[i] = complex(mpmath.coulombf(l, eta, rho))
        g[i] = complex(mpmath.coulombg(l, eta, rho))
    return (pref * f).reshape(np.shape(x)), (pref * g).reshape(np.shape(x))


def coulomb_f_en(x: npt.ArrayLike, k: float, z: float, m: float, l: int) -> npt.NDArray[np.complex128]:
    """`sqrt(2m/pi k) F_l(m z/k, k x)`; reduces to `riccati_bessel_en(x,k,l)` at z=0, m=1."""
    f, _ = _fg(np.asarray(x, dtype=np.complex128), k, z, m, l)
    return np.asarray(f, dtype=np.complex128)


def coulomb_g_en(x: npt.ArrayLike, k: float, z: float, m: float, l: int) -> npt.NDArray[np.complex128]:
    """Energy-normalized irregular Coulomb function `sqrt(2m/pi k) G_l`."""
    _, g = _fg(np.asarray(x, dtype=np.complex128), k, z, m, l)
    return np.asarray(g, dtype=np.complex128)


def coulomb_h1_en(x: npt.ArrayLike, k: float, z: float, m: float, l: int) -> npt.NDArray[np.complex128]:
    """Energy-normalized OUTGOING Coulomb function `sqrt(2m/pi k)(G_l + i F_l)`."""
    f, g = _fg(np.asarray(x, dtype=np.complex128), k, z, m, l)
    return np.asarray(g + 1j * f, dtype=np.complex128)
```

Export the three from `libs/qscat/qscat/special/__init__.py` (import + `__all__`).

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest libs/qscat/tests/test_coulomb.py -q`
Expected: PASS (6 tests). (Note: `riccati_hankel_en` = `j_l + i y_l`; confirm `G_l(0,ρ)=−ρ y_l(ρ)` sign works out so `G+iF → h1`; if the `h1` test fails on a sign, the reduction is `G_en+iF_en` vs `riccati_hankel_en` — adjust the TEST's expected only if the math genuinely differs, and document.)

- [ ] **Step 5: Type-check + lint + commit**

Run: `uv run mypy libs/qscat/qscat/special && uv run ruff check libs/qscat/qscat/special libs/qscat/tests/test_coulomb.py`
```bash
git add libs/qscat/pyproject.toml libs/qscat/qscat/special/coulomb.py libs/qscat/qscat/special/__init__.py libs/qscat/tests/test_coulomb.py
git commit -m "feat(h2+): energy-normalized Coulomb functions (mpmath, charge-z generalization)"
```

---

### Task 2: `charge` protocol + the ionic model (`qscat.model`)

**Files:**
- Modify: `libs/qscat/qscat/model/diatomic.py` (narrow protocol; `DiatomicResonanceModel.charge`)
- Create: `libs/qscat/qscat/model/ionic.py`
- Modify: `libs/qscat/qscat/model/library.py`, `libs/qscat/qscat/model/__init__.py`
- Test: `libs/qscat/tests/test_ionic_model.py`; extend `libs/qscat/tests/test_model.py`

**Interfaces:**
- `ResonanceModel` protocol: **remove `lam`** (the engine never calls it — verified: core reads only `mu/ell/v0/v_int/surface/hamiltonian/interaction_diag`), **add** `@property def charge(self) -> int`.
- `DiatomicResonanceModel`: add `charge: int = 0` (a defaulted field, LAST — existing fields have no defaults; registry entries keep working). Its `lam` method stays (a concrete detail its own `v_int` uses).
- Produces: `IonicResonanceModel` (frozen dataclass implementing the narrowed protocol) with the H₂⁺ form + `charge`, and an `H2P` registry instance.

**Design notes (extracted formulas — atomic units):**
- `v0(R) = V0·(e^{−2α(R−R0)} − 2 e^{−α(R−R0)})`, H₂⁺: `V0=0.1027, R0=2.0, α=0.69` (the ion Morse; the initial vibrational state lives here).
- `v_int(r,R) = −a1·(1 − tanh Q(R))·S(R)·(e^{−r²/3}/r)`, `Q(R)=(a2 − R − a3 R⁴)/7`, `S(R)=tanh(R/a4)⁴`, H₂⁺: `a1=1.6435, a2=6.2, a3=0.0125, a4=1.15` (the σ-capture interaction).
- `surface(r,R) = v0(R) + v_int(r,R) + ℓ(ℓ+1)/(2r²) + charge/r` (H₂⁺ charge=−1 → the `−1/r` electron–core Coulomb attraction).
- `mu=918.25, ell=1, charge=−1`. `hamiltonian`/`interaction_diag` identical to `DiatomicResonanceModel`'s (`hamiltonian_nd(tgrid,[1.0,mu],surface)` / `potential_nd(tgrid,v_int)`). All complex-safe (`np.asarray(..., complex128)`; never coerce to real — ECS tails).

- [ ] **Step 1: Write the failing test**

```python
# libs/qscat/tests/test_ionic_model.py
from __future__ import annotations

import numpy as np
from qscat.model import H2P, ResonanceModel


def test_h2p_is_a_resonance_model_with_charge():
    assert isinstance(H2P, ResonanceModel)
    assert H2P.charge == -1 and H2P.ell == 1 and H2P.mu == 918.25


def test_h2p_v0_is_the_ion_morse():
    R = np.array([2.0, 3.0, 8.0])
    V0, R0, a = 0.1027, 2.0, 0.69
    expect = V0 * (np.exp(-2 * a * (R - R0)) - 2 * np.exp(-a * (R - R0)))
    assert np.allclose(H2P.v0(R).real, expect, atol=1e-12)
    assert abs(H2P.v0(np.array([2.0]))[0].real + 0.1027) < 1e-12   # min -V0 at R0


def test_h2p_v_int_matches_sigma_capture():
    r, R = np.array([1.5]), np.array([2.5])
    a1, a2, a3, a4 = 1.6435, 6.2, 0.0125, 1.15
    Q = (a2 - R - a3 * R**4) / 7.0
    S = np.tanh(R / a4) ** 4
    E = np.exp(-r**2 / 3.0) / r
    expect = -a1 * (1 - np.tanh(Q)) * S * E
    assert np.allclose(H2P.v_int(r, R).real, expect, atol=1e-12)


def test_h2p_surface_has_coulomb_tail():
    # surface - (v0 + v_int + centrifugal) == charge/r == -1/r
    r, R = np.array([2.0]), np.array([3.0])
    cent = H2P.ell * (H2P.ell + 1) / (2.0 * r**2)
    tail = H2P.surface(r, R) - (H2P.v0(R) + H2P.v_int(r, R) + cent)
    assert np.allclose(tail, -1.0 / r, atol=1e-12)


def test_diatomic_models_are_neutral():
    from qscat.model import N2, F2
    assert N2.charge == 0 and F2.charge == 0
    assert isinstance(N2, ResonanceModel)
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest libs/qscat/tests/test_ionic_model.py -q`
Expected: FAIL — `cannot import name 'H2P'`.

- [ ] **Step 3: Narrow the protocol; add `charge`; create the ionic model**

In `model/diatomic.py`: delete the `lam` method from the `ResonanceModel` Protocol; add `@property def charge(self) -> int: ...` (with a docstring: "the Coulomb charge z for the channel functions; 0 for neutral targets, −1 for a singly-charged cation like H₂⁺"). Keep `DiatomicResonanceModel.lam` (concrete). Add `charge: int = 0` as the LAST field of `DiatomicResonanceModel`.

Create `model/ionic.py` with `IonicResonanceModel` (frozen dataclass, fields `mu, ell, charge, V0, R0, alpha, a1, a2, a3, a4`) implementing `v0`/`v_int`/`surface`/`hamiltonian`/`interaction_diag` per the Design notes (complex-safe). Add `H2P = IonicResonanceModel(mu=918.25, ell=1, charge=-1, V0=0.1027, R0=2.0, alpha=0.69, a1=1.6435, a2=6.2, a3=0.0125, a4=1.15)` to `model/library.py`; export `IonicResonanceModel`, `H2P` from `model/__init__.py`.

In `test_model.py`: if the protocol change or a `lam`-in-protocol assertion breaks an existing test, update it (the engine doesn't use `lam`; `DiatomicResonanceModel.lam` still exists for the N2 checks).

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest libs/qscat/tests/test_ionic_model.py libs/qscat/tests/test_model.py -q`
Expected: PASS.

- [ ] **Step 5: Type-check + lint + commit**

Run: `uv run mypy libs/qscat/qscat/model && uv run ruff check libs/qscat/qscat/model libs/qscat/tests/test_ionic_model.py`
```bash
git add libs/qscat/qscat/model/ libs/qscat/tests/test_ionic_model.py libs/qscat/tests/test_model.py
git commit -m "feat(h2+): charge-aware ResonanceModel protocol + H2+ ionic model (H2P)"
```

---

### Task 3: Charge-aware `channel_vector`

**Files:**
- Modify: `libs/qscat/qscat/core/channels.py`
- Test: `libs/qscat/tests/test_core_channels.py` (or wherever `channel_vector` is tested — else add)

**Interfaces:**
- `channel_vector(tgrid, k, chi_v, l, *, charge=0)` — when `charge==0`, the existing fast `riccati_bessel_en` path (neutrals bit-unchanged); when `charge!=0`, the electronic factor is `coulomb_f_en(r, k, charge, 1.0, l)` (mass-1 electron). Everything else (the `sqrt_weights` factor, `chi` c-normalization, `real_mask`) identical.

**Design notes:** the ONLY change is the electronic radial function: `f_vals = riccati_bessel_en(g_r.real_points, k, l) if charge == 0 else coulomb_f_en(g_r.real_points, k, float(charge), 1.0, l)`. `coulomb_f_en` returns complex; `f_coeff` becomes complex (fine — `psi` is already complex). Keep `charge=0` the default so all existing callers (VE/DA) are unchanged.

- [ ] **Step 1: Write the failing test** (append to the channels test)

```python
def test_channel_vector_charge_zero_unchanged():
    # charge=0 must be identical to the default (no Coulomb path)
    from qscat.core.channels import channel_vector
    from qscat.core.grids import electronic_grid, nuclear_grid
    from qscat.dvr import TensorGrid
    import numpy as np
    tg = TensorGrid([electronic_grid(r_max=14.0, order=6, n_complex=4),
                     nuclear_grid(r_max=20.0, n_complex=4, quadrature=8)])
    chi = np.zeros(tg.grids[1].n, dtype=np.complex128); chi[0] = 1.0
    a = channel_vector(tg, 0.6, chi, 1)
    b = channel_vector(tg, 0.6, chi, 1, charge=0)
    assert np.array_equal(a, b)


def test_channel_vector_coulomb_is_finite_and_differs():
    from qscat.core.channels import channel_vector
    from qscat.core.grids import electronic_grid, nuclear_grid
    from qscat.dvr import TensorGrid
    import numpy as np
    tg = TensorGrid([electronic_grid(r_max=14.0, order=6, n_complex=4),
                     nuclear_grid(r_max=20.0, n_complex=4, quadrature=8)])
    chi = np.zeros(tg.grids[1].n, dtype=np.complex128); chi[0] = 1.0
    free = channel_vector(tg, 0.6, chi, 1)
    coul = channel_vector(tg, 0.6, chi, 1, charge=-1)
    assert np.all(np.isfinite(coul)) and not np.allclose(coul, free)
```

- [ ] **Step 2–5:** implement the `charge`-dispatch (per Design notes), run the tests (PASS), confirm no VE/DA regression (`uv run pytest libs/qscat/tests/test_dissociation.py libs/qscat/tests/test_core_driven.py -q -m "not slow"`), mypy + ruff, commit:
```bash
git add libs/qscat/qscat/core/channels.py libs/qscat/tests/test_core_channels.py
git commit -m "feat(h2+): charge-aware channel_vector (Coulomb incident when charge != 0)"
```

---

### Task 4: `dr_cross_section` — the Rydberg-loop DR generalization

**Files:**
- Modify: `libs/qscat/qscat/core/dissociation.py`
- Modify: `libs/qscat/qscat/core/__init__.py`
- Test: `libs/qscat/tests/test_dissociation.py`

**Interfaces:**
- Produces:
  ```python
  def dr_cross_section(tgrid, model, eps, chi, v_init, E, *, n_channels=3, ordering="COLAMD"): ...
  ```
  σ_DR per Rydberg channel, shape `(n_channels,)` scalar E / `(len(E), n_channels)` array E. `σ_n=0` where `E≤0` or `E_tot−E_ryd(n)≤0` (closed).

**Design notes:** `dr_cross_section` is `da_cross_section` **generalized** (read `da_cross_section` first and mirror it):
- The Rydberg exit states are `anion_electronic_states(tgrid.grids[0], model, R_inf, n_states=n_channels)` → `(eps_e[n], phi_e[n])` — the SAME bound-electronic-state solver (they're bound below the −1/r continuum). Their energies are `E_ryd(n)=eps_e[n]`.
- `V_DR = v_dr_diag(tgrid, model)` — SAME.
- The incident is **Coulomb**: `psi_i = channel_vector(tgrid, k, chi[v_init], model.ell, charge=model.charge)` (charge=−1 → Coulomb). Get `Ψ₊` via the driven sweep — either reuse `ve_cross_section(..., return_wavefunction=True)` IF it can pass `charge` to `channel_vector` (it currently can't — so either (a) add a `charge` pass-through to `ve_cross_section`/`_sigma_at_one_energy`, or (b) replicate the ~12-line driven loop inside `dr_cross_section` using the charged `channel_vector`). **Prefer (b)** for a focused change: build `H`, `v_diag=interaction_diag`, per-E `SparseLU`/`refactor`, `psi_plus = psi_i + lu.solve(v_diag*psi_i)`.
- Per open channel n: `E_DR = E_tot − eps_e[n]`; `K = √(2μ E_DR)`; `Y_coeff = riccati_bessel_en_mass(g_R.real_points, K, 0, μ) * sqrt_weights()[1].ravel()`; `Phi = tgrid.outer([phi_e[n], Y_coeff])`, `Phi[~mask]=0`; `T = c_product(Phi, V_DR·psi_plus)`; `σ_n = 4π³|T|²/2E`.

Export `dr_cross_section` from `core/__init__.py`.

- [ ] **Step 1: Write the failing test** (append; uses a REDUCED H₂⁺ proxy grid — the fine points come from Task 5's config, but a self-contained small grid here keeps the lib test independent)

```python
def _h2p_proxy():
    # small ionic proxy: electronic to ~60 bohr (holds a couple Rydberg states +
    # the incident), nuclear to ~14. Big enough for well-posedness, laptop-fast.
    from qscat.core.grids import electronic_grid, nuclear_grid
    from qscat.dvr import TensorGrid
    return TensorGrid([electronic_grid(r_max=60.0, order=8, n_complex=6),
                       nuclear_grid(r_max=22.0, n_complex=6, quadrature=10)])


@pytest.mark.slow
def test_dr_wellposed_and_threshold_ordered():
    from qscat.core.dissociation import dr_cross_section
    from qscat.core.vibrational import vibrational_states
    from qscat.model import H2P
    tg = _h2p_proxy()
    eps, chi = vibrational_states(tg.grids[1], H2P.mu, 3, H2P.v0)
    E = np.array([0.01, 0.03])
    s = dr_cross_section(tg, H2P, eps, chi, 0, E, n_channels=2)
    assert s.shape == (2, 2)
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)
```

- [ ] **Step 2–6:** run→fail; implement `dr_cross_section` (mirror `da_cross_section`, add the Rydberg loop + charged incident); run→pass (the `@slow` proxy test — the electronic Coulomb `channel_vector` via mpmath makes it heavier, minutes); import-guard + mypy + ruff; commit:
```bash
git add libs/qscat/qscat/core/dissociation.py libs/qscat/qscat/core/__init__.py libs/qscat/tests/test_dissociation.py
git commit -m "feat(h2+): dr_cross_section (Rydberg-channel loop + Coulomb incident)"
```

---

### Task 5: H₂⁺ discretization + config (`validation/h2plus`)

**Files:**
- Create: `validation/h2plus/__init__.py`, `validation/h2plus/config.py`
- Modify: `libs/qscat/qscat/core/grids.py` (only if the exp-growth ECS tail isn't already expressible)
- Test: `validation/h2plus/test_config.py`

**Interfaces:**
- `validation/h2plus/config.py`: `full_grid()` (the real H₂⁺ deck — electronic real→1300 + exp-ECS tail 5°; nuclear real→14 + exp-ECS tail 22°; order 8) and `proxy_grid()` (the reduced laptop grid: electronic r_max~60, nuclear~14), both `TensorGrid`; plus energy grid + `n_channels`.

**Design notes:** the electronic real region is 5 segments `[(10,1.0),(10,4.0),(16,20.0),(20,100.0),(120,1300.0)]` (element sizes 0.1/0.3/1.0/4.0/10.0 — `(n_elements, endpoint)` pairs); the nuclear real region `[(5,1.0),(20,4.0),(67,14.0)]` (sizes 0.2/0.15/0.15). **The ECS tail is EXP-GROWTH, which `segmented_grid` (uniform-per-segment) cannot express — add a new builder** to `grids.py`:
```python
def fem_grid_exp_tail(real_segments, *, angle_deg, quadrature, tail_n, tail_alpha=0.2, tail_skip=2, x_min=0.0) -> FemDvrEcsGrid: ...
```
It assembles the real segments exactly like `segmented_grid` (uniform `ElementSpec(h)` per `(n,endpoint)`), then appends `tail_n` complex `ElementSpec(h_i, angle_deg)` whose lengths are the existing `_ecs_tail(base, tail_n, skip=tail_skip, alpha=tail_alpha)` sequence (`base` = the last real element length: `skip` elements at `base`, then `base·e^{alpha·(i−skip+1)}`) — reusing the module's `_ecs_tail` helper verbatim. Export it. Then: electronic full grid = `fem_grid_exp_tail([(10,1.0),(10,4.0),(16,20.0),(20,100.0),(120,1300.0)], angle_deg=5.0, quadrature=8, tail_n=25)`; nuclear = `fem_grid_exp_tail([(5,1.0),(20,4.0),(67,14.0)], angle_deg=22.0, quadrature=8, tail_n=25)`. `proxy_grid` shrinks the electronic real region to `[(10,1.0),(10,4.0),(16,20.0),(10,60.0)]` (drop 100/1300; ~60 bohr) with a smaller `tail_n` so a laptop SuperLU solve is feasible; nuclear proxy `[(5,1.0),(20,4.0),(40,14.0)]`. `full_grid` is Docker/MUMPS-sized.

- [ ] **Steps:** TDD — first the `fem_grid_exp_tail` builder + a test (real endpoint = R0; the tail element lengths grow as `_ecs_tail`; complex on the tail); then `config.py` `full_grid()`/`proxy_grid()` with tests (`full_grid` electronic `R0≈1300`, `proxy_grid` electronic `R0≈60`, both nuclear `R0≈14`, complex dtype); mypy/ruff; commit:
```bash
git add validation/h2plus/ libs/qscat/qscat/core/grids.py
git commit -m "feat(h2+): H2+ discretization (full 1300-bohr deck + laptop proxy grid)"
```

---

### Task 6: Small-proxy DR validation + Docker-ready + docs

**Files:**
- Create: `validation/h2plus/dr.py`, `validation/h2plus/test_dr.py`
- Create: `docs/physics/h2plus-dr.md`
- Modify: `CLAUDE.md`, `docs/superpowers/specs/2026-07-27-da-cross-sections-design.md` (mark sub-project D delivered)

**Interfaces:**
- `validation/h2plus/dr.py`: `compute_dr(grid, *, energies, n_channels)` → σ_DR array, using `H2P` + `dr_cross_section`; a `main()` for a Docker smoke run on `full_grid()` (a couple of energies) writing a `.npz`.

**Design notes:** the acceptance is a **well-posedness/threshold** gate on `proxy_grid()` (finite, ≥0, per-channel thresholds respected, the c-product-vs-conjugated-dot difference negligible — compute both and assert `|Δσ|/σ < ~1e-3` on the proxy), NOT a converged σ_DR (the real grid is 1300 bohr). Document the Docker/MUMPS path for the full deck (a smoke `main()` run, not a suite gate). The doc states the ionic physics (Coulomb tail, Rydberg series, DR vs DA), the reduced-proxy caveat, and the eMoScat `sH1`-bug / `zdotc`-vs-c-product findings.

- [ ] **Step 1: Write the failing gate test** (`@slow`, proxy grid): σ_DR finite/≥0; channel `n` closed below `eps_e[n]−eps[0]`; the c-product vs conjugated-dot difference is < 1e-3 on the proxy.
- [ ] **Step 2–4:** implement `dr.py` (+ the conjugated-dot comparison helper), run the proxy gate (PASS, `@slow`), and a `main()` smoke path.
- [ ] **Step 5: Docs + Docker-ready:** write `docs/physics/h2plus-dr.md`; add a `qscat.special.coulomb` + `qscat.model.H2P` + `dr_cross_section` line to `CLAUDE.md`; note the `docker/` deck for the full run; mark sub-project D delivered in the DA/DR spec.
- [ ] **Step 6: Commit:**
```bash
git add validation/h2plus/ docs/physics/h2plus-dr.md CLAUDE.md docs/superpowers/specs/2026-07-27-da-cross-sections-design.md
git commit -m "feat(h2+): small-proxy DR validation + Docker-ready + docs"
```

---

## Verification (whole sub-project)

- `uv run pytest -q -m "not slow"` pass; the `@slow` H₂⁺ proxy tests pass.
- `uv run pytest libs/qscat/tests/test_core_no_model_import.py -q` pass (DR keeps the core/model boundary).
- `uv run mypy libs/qscat/qscat` 0 errors; `uv run ruff check .` clean.
- Coulomb functions match z→0→Bessel + the known value + accept complex args; the ionic `H2P` model matches the extracted formulas and is a `ResonanceModel` with `charge=−1`; neutrals keep `charge=0` and `channel_vector(charge=0)` is bit-unchanged.
- `dr_cross_section` is well-posed on the proxy, respects per-channel thresholds, and the c-product choice is justified (≪1e-3 vs the conjugated dot on the proxy); the full H₂⁺ deck builds and smoke-runs (Docker/MUMPS path documented).
- `docs/physics/h2plus-dr.md` + `CLAUDE.md` updated; DA/DR spec's sub-project D marked delivered.

## Out of scope (this plan)

- **A converged full-size σ_DR(E) curve/figure** (needs the 1300-bohr Docker/MUMPS run) — proxy gate + Docker-ready deck delivered; converged run is a follow-on.
- **The π channel** (`p_pi_potential`); **optimizing the Coulomb functions** (Rust/COULCC port); **rotational/coupled-channel DR.**
