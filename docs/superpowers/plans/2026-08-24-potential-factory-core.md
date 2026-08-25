# Potential Factory — Core (ansatz, tracker, fitter, round-trip oracle) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A toy-stage factory under `projects/potential_factory/` that takes a tiered `Target` (neutral curve, pole curves, energy-dependent width) and fits a `FlexibleDiatomicModel` to it in stages, proven by round-tripping the existing N₂/NO/F₂ models' own calculated curves back to their published constants and curves.

**Architecture:** Three layers, each a single file with one job: `ansatz.py` (the analytic, ECS-safe surface family — a strict superset of `DiatomicResonanceModel`), `tracker.py` (fixed-`R` resonance-pole solver with a spurious-pole gate, the c-product Hellmann–Feynman gradient, Newton in `(λ, α)`, and continuation over `R`), `fit.py` + `report.py` (staged T0→T1→T3 least squares that stop-and-report). `extract.py` produces a `Target` *from* an existing model with the repo's own forward models, which is what makes the round-trip an exact oracle. No new physics: every forward quantity is computed by existing `qscat.core` code.

**Tech Stack:** Python 3.12, numpy, scipy (`least_squares`, `brentq`), pytest (`@pytest.mark.slow` for the round trips); `qscat.dvr`, `qscat.ecs`, `qscat.core.{lcp,grids,dissociation,vibrational,nrm}`, `qscat.linalg.c_product`.

**Spec:** `docs/superpowers/specs/2026-08-24-potential-factory-design.md` (read it first; the survey `docs/physics/potential-factory-options.md` explains why the ansatz has two free functions of `R`).

## Global Constraints

- Atomic units throughout; eV appears only in messages/docs (`EV = 27.211386` Hartree⁻¹).
- `r` and `R` may be COMPLEX (ECS tails): every potential term converts to `complex128` and never coerces to real (same rule as `qscat.model.diatomic`).
- `projects/potential_factory` may import `qscat.model` and `qscat.core`; **`qscat.core` must not import it** (`libs/qscat/tests/test_core_no_model_import.py` stays green).
- All intra-repo imports package-absolute (`from projects.potential_factory.ansatz import ...`); tests live beside the code as `projects/potential_factory/test_*.py`; run with `uv run --no-sync pytest projects/potential_factory -q`.
- Lint/format before every commit: `uv run ruff check . && uv run ruff format .`.
- Every `git add` lists files explicitly (never `-a`).
- A fake angle-stable pole is accepted by `qscat.ecs.find_resonance_pole` near threshold (spike, 2026-08-24: `(0.11 eV, 0.26 eV)` for `l ≤ 1`); every pole this package returns must pass the gate in Task 3.
- Nothing is fitted to an experimental observable.

---

### Task 1: Package scaffold and the ansatz `FlexibleDiatomicModel`

**Files:**
- Create: `projects/potential_factory/__init__.py` (empty)
- Create: `projects/potential_factory/ansatz.py`
- Test: `projects/potential_factory/test_ansatz.py`

**Interfaces:**
- Consumes: `qscat.model.DiatomicResonanceModel`, `qscat.dvr.TensorGrid`, `qscat.dvr.hamiltonian_nd`, `qscat.dvr.potential_nd`.
- Produces:
  - `y_p(R, R_e, p) -> ndarray[complex128]`
  - `SmoothR(f_inf, f_0, f_1, R_f, coeffs=(), R_e=2.0, p=3)` callable `(R) -> complex128 array`
  - `FlexibleDiatomicModel(mu, ell, D_e, R_e, betas, p, lam, alpha, shell, alpha_b, r_b, charge=0)` implementing `ResonanceModel` (`v0`, `v_int`, `surface`, `hamiltonian`, `interaction_diag`, plus `lam_R(R)`, `alpha_R(R)`, `shell_R(R)` helpers)
  - `from_diatomic(model: DiatomicResonanceModel) -> FlexibleDiatomicModel`

- [ ] **Step 1: Write the failing tests**

```python
# projects/potential_factory/test_ansatz.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.model import F2, N2, NO, ResonanceModel

from projects.potential_factory.ansatz import FlexibleDiatomicModel, SmoothR, from_diatomic, y_p

R_REAL = np.linspace(1.2, 6.0, 25)
R_CPLX = 12.0 + np.linspace(0.1, 4.0, 9) * np.exp(1j * np.deg2rad(35.0))
r_REAL = np.linspace(0.05, 12.0, 40)
r_CPLX = 16.0 + np.linspace(0.1, 6.0, 7) * np.exp(1j * np.deg2rad(40.0))


@pytest.mark.parametrize("model", [N2, NO, F2], ids=["N2", "NO", "F2"])
def test_from_diatomic_reproduces_published_form_to_roundoff(model):
    flex = from_diatomic(model)
    for R in (R_REAL, R_CPLX):
        np.testing.assert_allclose(flex.v0(R), model.v0(R), rtol=0, atol=1e-14)
        np.testing.assert_allclose(flex.lam_R(R), model.lam(R), rtol=0, atol=1e-14)
        for r in (r_REAL, r_CPLX):
            rr, RR = np.meshgrid(r, R, indexing="ij")
            np.testing.assert_allclose(flex.v_int(rr, RR), model.v_int(rr, RR), rtol=0, atol=1e-14)
            np.testing.assert_allclose(flex.surface(rr, RR), model.surface(rr, RR), rtol=0, atol=1e-13)


def test_flexible_model_is_a_resonance_model():
    assert isinstance(from_diatomic(N2), ResonanceModel)


def test_emo_with_one_beta_is_morse():
    flex = from_diatomic(N2)
    assert flex.betas == (N2.alpha0,)
    assert flex.D_e == N2.D0 and flex.R_e == N2.R0


def test_y_p_is_zero_at_R_e_and_tends_to_one():
    assert y_p(2.0, 2.0, 3) == 0.0
    assert abs(y_p(1e6, 2.0, 3) - 1.0) < 1e-12


def test_smooth_r_reduces_to_houfek_sigmoid():
    lam0 = (N2.lambda_c - N2.lambda_inf) * (1 + np.exp(N2.lambda_1 * (N2.R_c - N2.R_lambda)))
    s = SmoothR(f_inf=N2.lambda_inf, f_0=lam0, f_1=N2.lambda_1, R_f=N2.R_lambda)
    np.testing.assert_allclose(s(R_REAL), N2.lam(R_REAL), atol=1e-14)


def test_shell_term_adds_a_barrier():
    base = from_diatomic(N2)
    shell = SmoothR(f_inf=0.5, f_0=0.0, f_1=1.0, R_f=0.0)
    with_shell = FlexibleDiatomicModel(
        **{**base.__dict__, "shell": shell, "alpha_b": 2.0, "r_b": 3.0}
    )
    r = np.array([3.0])
    R = np.array([2.0])
    assert with_shell.v_int(r, R).real > base.v_int(r, R).real
    assert abs(with_shell.v_int(r, R) - base.v_int(r, R) - 0.5) < 1e-14


def test_complex_inputs_are_not_coerced_to_real():
    flex = from_diatomic(N2)
    out = flex.v0(R_CPLX)
    assert out.dtype == np.complex128 and np.any(out.imag != 0.0)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run --no-sync pytest projects/potential_factory/test_ansatz.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'projects.potential_factory'`

- [ ] **Step 3: Write the ansatz**

```python
# projects/potential_factory/__init__.py
"""Potential factory (toy stage): fit a FlexibleDiatomicModel to a tiered Target."""
```

```python
# projects/potential_factory/ansatz.py
"""The analytic, ECS-safe surface family the factory emits.

`FlexibleDiatomicModel` is a strict superset of `qscat.model.DiatomicResonanceModel`:
an Expanded-Morse-Oscillator neutral curve (Le Roy; one beta == Morse), a Gaussian
electron-molecule well whose depth `lam(R)` AND range `alpha(R)` are smooth
functions of `R`, and an optional repulsive Gaussian shell `shell(R) exp(-alpha_b
(r - r_b)^2)`. `from_diatomic` embeds the published models exactly, so N2/NO/F2
are points of this parameter space (the round-trip oracle relies on that).

Every term is entire in `r`; in `R` the only singularities are the poles of
`y_p` at `|R| = R_e`, which an ECS tail pivoted at `R_0 > R_e` never reaches.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.dvr import TensorGrid, hamiltonian_nd, potential_nd
from qscat.model import DiatomicResonanceModel

__all__ = ["y_p", "SmoothR", "FlexibleDiatomicModel", "from_diatomic"]


def y_p(R: npt.ArrayLike, R_e: float, p: int) -> npt.NDArray[np.complex128]:
    """Le Roy's dimensionless radial variable `(R^p - R_e^p) / (R^p + R_e^p)`."""
    Rc = np.asarray(R, dtype=np.complex128)
    out = (Rc**p - R_e**p) / (Rc**p + R_e**p)
    return np.asarray(out, dtype=np.complex128)


@dataclass(frozen=True)
class SmoothR:
    """`f(R) = f_inf + f_0 / (1 + exp(f_1 (R - R_f))) * (1 + sum_i coeffs[i] y_p(R)^(i+1))`.

    With `coeffs == ()` this is exactly Houfek's sigmoid
    `lambda_inf + lambda_0 / (1 + exp(lambda_1 (R - R_lambda)))`.
    """

    f_inf: float
    f_0: float
    f_1: float
    R_f: float
    coeffs: tuple[float, ...] = ()
    R_e: float = 2.0
    p: int = 3

    def __call__(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        Rc = np.asarray(R, dtype=np.complex128)
        sig = self.f_0 / (1.0 + np.exp(self.f_1 * (Rc - self.R_f)))
        poly = np.ones_like(Rc)
        if self.coeffs:
            y = y_p(Rc, self.R_e, self.p)
            poly = poly + sum(c * y ** (i + 1) for i, c in enumerate(self.coeffs))
        return np.asarray(self.f_inf + sig * poly, dtype=np.complex128)


@dataclass(frozen=True)
class FlexibleDiatomicModel:
    """EMO neutral curve + Gaussian well with `lam(R)`, `alpha(R)` + optional shell."""

    mu: float
    ell: int
    D_e: float
    R_e: float
    betas: tuple[float, ...]
    p: int
    lam: SmoothR
    alpha: SmoothR
    shell: SmoothR | None
    alpha_b: float
    r_b: float
    charge: int = 0

    # -- neutral curve -------------------------------------------------------
    def beta_R(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        Rc = np.asarray(R, dtype=np.complex128)
        y = y_p(Rc, self.R_e, self.p)
        out = sum(b * y**i for i, b in enumerate(self.betas))
        return np.asarray(out + 0.0 * Rc, dtype=np.complex128)

    def v0(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """EMO: `D_e [ (1 - exp(-beta(R)(R - R_e)))^2 - 1 ]`, minimum -D_e at R_e, 0 at infinity."""
        Rc = np.asarray(R, dtype=np.complex128)
        e = np.exp(-self.beta_R(Rc) * (Rc - self.R_e))
        return np.asarray(self.D_e * ((1.0 - e) ** 2 - 1.0), dtype=np.complex128)

    # -- interaction ---------------------------------------------------------
    def lam_R(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        return self.lam(R)

    def alpha_R(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        return self.alpha(R)

    def shell_R(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        Rc = np.asarray(R, dtype=np.complex128)
        if self.shell is None:
            return np.zeros_like(Rc)
        return self.shell(Rc)

    def v_int(self, r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        rr = np.asarray(r, dtype=np.complex128)
        Rc = np.asarray(R, dtype=np.complex128)
        well = -self.lam_R(Rc) * np.exp(-self.alpha_R(Rc) * rr**2)
        if self.shell is None:
            return np.asarray(well, dtype=np.complex128)
        barrier = self.shell_R(Rc) * np.exp(-self.alpha_b * (rr - self.r_b) ** 2)
        return np.asarray(well + barrier, dtype=np.complex128)

    def surface(self, r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        rr = np.asarray(r, dtype=np.complex128)
        out = self.v0(R) + self.ell * (self.ell + 1) / (2.0 * rr**2) + self.v_int(rr, R)
        return np.asarray(out, dtype=np.complex128)

    def hamiltonian(self, tgrid: TensorGrid) -> sp.csr_matrix:
        return hamiltonian_nd(tgrid, [1.0, self.mu], self.surface)

    def interaction_diag(self, tgrid: TensorGrid) -> npt.NDArray[np.complex128]:
        return potential_nd(tgrid, self.v_int)

    def with_shell(self, shell: SmoothR | None, alpha_b: float, r_b: float) -> FlexibleDiatomicModel:
        return replace(self, shell=shell, alpha_b=alpha_b, r_b=r_b)


def from_diatomic(model: DiatomicResonanceModel) -> FlexibleDiatomicModel:
    """Embed a published `DiatomicResonanceModel` exactly (to round-off)."""
    lam0 = (model.lambda_c - model.lambda_inf) * (
        1.0 + np.exp(model.lambda_1 * (model.R_c - model.R_lambda))
    )
    return FlexibleDiatomicModel(
        mu=model.mu,
        ell=model.ell,
        D_e=model.D0,
        R_e=model.R0,
        betas=(model.alpha0,),
        p=3,
        lam=SmoothR(f_inf=model.lambda_inf, f_0=float(lam0), f_1=model.lambda_1, R_f=model.R_lambda, R_e=model.R0),
        alpha=SmoothR(f_inf=model.alpha_c, f_0=0.0, f_1=1.0, R_f=0.0, R_e=model.R0),
        shell=None,
        alpha_b=1.0,
        r_b=0.0,
        charge=model.charge,
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run --no-sync pytest projects/potential_factory/test_ansatz.py -q`
Expected: 9 passed (3 parametrized + 6)

- [ ] **Step 5: Lint, then commit**

```bash
uv run ruff check projects/potential_factory && uv run ruff format projects/potential_factory
git add projects/potential_factory/__init__.py projects/potential_factory/ansatz.py projects/potential_factory/test_ansatz.py
git commit -m "feat(factory): FlexibleDiatomicModel — EMO v0 + Gaussian well with lam(R), alpha(R) + optional shell; embeds N2/NO/F2 exactly"
```

---

### Task 2: Flat parameter names (`params` / `with_params`)

**Files:**
- Modify: `projects/potential_factory/ansatz.py` (append)
- Test: `projects/potential_factory/test_ansatz.py` (append)

**Interfaces:**
- Produces: `params(model) -> dict[str, float]` with keys `D_e, R_e, beta0..betaN, lam.f_inf, lam.f_0, lam.f_1, lam.R_f, lam.c0.., alpha.f_inf, ..., shell.f_inf, ..., alpha_b, r_b` (shell keys only when `shell is not None`); `with_params(model, updates: Mapping[str, float]) -> FlexibleDiatomicModel` (unknown key → `KeyError`); `pack(model, names) -> ndarray`, `unpack(model, names, x) -> FlexibleDiatomicModel`. `mu`, `ell`, `p`, `charge` are never parameters.

- [ ] **Step 1: Write the failing tests**

```python
# append to projects/potential_factory/test_ansatz.py
from projects.potential_factory.ansatz import pack, params, unpack, with_params  # noqa: E402


def test_params_round_trip_and_unknown_key():
    flex = from_diatomic(NO)
    p = params(flex)
    assert p["D_e"] == NO.D0 and p["beta0"] == NO.alpha0 and p["lam.f_inf"] == NO.lambda_inf
    assert "shell.f_inf" not in p
    back = with_params(flex, {"lam.f_inf": 7.0, "beta0": 1.5})
    assert back.lam.f_inf == 7.0 and back.betas == (1.5,) and back.D_e == NO.D0
    with pytest.raises(KeyError):
        with_params(flex, {"nope": 1.0})


def test_pack_unpack_are_inverse():
    flex = from_diatomic(F2)
    names = ["lam.f_inf", "lam.f_0", "alpha.f_inf"]
    x = pack(flex, names)
    assert x.shape == (3,)
    again = unpack(flex, names, x * 1.0)
    assert params(again) == params(flex)
    moved = unpack(flex, names, x + 0.1)
    assert abs(moved.lam.f_inf - flex.lam.f_inf - 0.1) < 1e-15
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run --no-sync pytest projects/potential_factory/test_ansatz.py -q -k "params or pack"`
Expected: FAIL with `ImportError: cannot import name 'pack'`

- [ ] **Step 3: Implement**

```python
# append to projects/potential_factory/ansatz.py
from collections.abc import Mapping, Sequence  # noqa: E402  (move to the top when editing)


def _smooth_params(prefix: str, s: SmoothR) -> dict[str, float]:
    out = {f"{prefix}.f_inf": s.f_inf, f"{prefix}.f_0": s.f_0, f"{prefix}.f_1": s.f_1, f"{prefix}.R_f": s.R_f}
    for i, c in enumerate(s.coeffs):
        out[f"{prefix}.c{i}"] = c
    return out


def _smooth_update(s: SmoothR, prefix: str, upd: Mapping[str, float]) -> SmoothR:
    kw = {"f_inf": s.f_inf, "f_0": s.f_0, "f_1": s.f_1, "R_f": s.R_f}
    coeffs = list(s.coeffs)
    for key, val in upd.items():
        if not key.startswith(prefix + "."):
            continue
        field = key[len(prefix) + 1 :]
        if field in kw:
            kw[field] = float(val)
        elif field.startswith("c") and field[1:].isdigit() and int(field[1:]) < len(coeffs):
            coeffs[int(field[1:])] = float(val)
        else:
            raise KeyError(key)
    return replace(s, coeffs=tuple(coeffs), **kw)


def params(model: FlexibleDiatomicModel) -> dict[str, float]:
    """Flat, ordered name -> value map of every FITTABLE parameter."""
    out: dict[str, float] = {"D_e": model.D_e, "R_e": model.R_e}
    for i, b in enumerate(model.betas):
        out[f"beta{i}"] = b
    out.update(_smooth_params("lam", model.lam))
    out.update(_smooth_params("alpha", model.alpha))
    if model.shell is not None:
        out.update(_smooth_params("shell", model.shell))
        out["alpha_b"] = model.alpha_b
        out["r_b"] = model.r_b
    return out


def with_params(model: FlexibleDiatomicModel, updates: Mapping[str, float]) -> FlexibleDiatomicModel:
    """A copy with the named parameters replaced. Unknown names raise KeyError."""
    known = params(model)
    for key in updates:
        if key not in known:
            raise KeyError(key)
    betas = list(model.betas)
    for i in range(len(betas)):
        if f"beta{i}" in updates:
            betas[i] = float(updates[f"beta{i}"])
    return replace(
        model,
        D_e=float(updates.get("D_e", model.D_e)),
        R_e=float(updates.get("R_e", model.R_e)),
        betas=tuple(betas),
        lam=_smooth_update(model.lam, "lam", updates),
        alpha=_smooth_update(model.alpha, "alpha", updates),
        shell=None if model.shell is None else _smooth_update(model.shell, "shell", updates),
        alpha_b=float(updates.get("alpha_b", model.alpha_b)),
        r_b=float(updates.get("r_b", model.r_b)),
    )


def pack(model: FlexibleDiatomicModel, names: Sequence[str]) -> npt.NDArray[np.float64]:
    p = params(model)
    return np.array([p[n] for n in names], dtype=np.float64)


def unpack(model: FlexibleDiatomicModel, names: Sequence[str], x: npt.ArrayLike) -> FlexibleDiatomicModel:
    xs = np.asarray(x, dtype=np.float64)
    return with_params(model, dict(zip(names, xs.tolist(), strict=True)))
```

(Move the `Mapping, Sequence` import to the module's import block and add the four names to `__all__`.)

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run --no-sync pytest projects/potential_factory/test_ansatz.py -q`
Expected: 11 passed

- [ ] **Step 5: Commit**

```bash
uv run ruff check projects/potential_factory && uv run ruff format projects/potential_factory
git add projects/potential_factory/ansatz.py projects/potential_factory/test_ansatz.py
git commit -m "feat(factory): flat parameter names — params/with_params/pack/unpack"
```

---

### Task 3: `ElectronicPair` and the gated fixed-`R` pole solver

**Files:**
- Create: `projects/potential_factory/tracker.py`
- Test: `projects/potential_factory/test_tracker.py`

**Interfaces:**
- Consumes: `qscat.core.grids.electronic_grid`, `qscat.dvr.kinetic`, `qscat.dvr.eigen`, `qscat.ecs.find_resonance_pole`.
- Produces:
  - `Pole(energy: complex, residual: float)` with `.gamma` (= `max(0, -2 Im)`) and `.shift` (= `Re`)
  - `ElectronicPair(angles=(35.0, 44.0), r_max=16.0, order=7, n_complex=6)` with `.grid_a`, `.grid_b`, `.hamiltonians(v_fn) -> (Ha, Hb)`, and `.pole(v_fn, window, *, resid_tol=1e-3, gate_frac=0.05, e_floor=0.006) -> Pole | None`
  - `Window = tuple[float, float, float, float]` (re_lo, re_hi, im_lo, im_hi)
  - `DEFAULT_RES_WINDOW = (0.006, 0.8, -0.6, 0.0)`; `DEFAULT_BOUND_WINDOW = (-0.5, -1e-6, -1e-6, 1e-6)`

The gate is the spec's "angle-stable AND residual ≪ Γ AND Re E above a threshold floor": accept iff `resid < resid_tol` and (`gamma == 0` or `resid < gate_frac * gamma`) and (`Re E >= e_floor` or the pole is bound: `Re E < 0` and `|Im E| < 1e-6`).

- [ ] **Step 1: Write the failing tests**

```python
# projects/potential_factory/test_tracker.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.model import N2

from projects.potential_factory.tracker import DEFAULT_RES_WINDOW, ElectronicPair, Pole

EV = 27.211386


@pytest.fixture(scope="module")
def pair():
    return ElectronicPair()


def _v_n2(R):
    def v(r):
        return N2.surface(r, R)

    return v


def test_n2_pole_at_R0_matches_documented_values(pair):
    # docs/physics/n2-resonance.md: E_res(R0) = 2.445 eV, Gamma(R0) = 0.455 eV.
    p = pair.pole(_v_n2(N2.R0), DEFAULT_RES_WINDOW)
    assert isinstance(p, Pole)
    shift_eV = (p.shift - N2.v0(N2.R0).real) * EV
    assert abs(shift_eV - 2.445) < 0.02
    assert abs(p.gamma * EV - 0.455) < 0.02
    assert p.residual < 0.05 * p.gamma


def test_gate_rejects_the_spike_fake_pole(pair):
    # Bare l=1 Gaussian well with lam=0.5, alpha=0.3 has NO resonance; the
    # 2026-08-24 spike saw find_resonance_pole return (0.11 eV, 0.26 eV) here.
    def v(r):
        rr = np.asarray(r, dtype=np.complex128)
        return -0.5 * np.exp(-0.3 * rr**2) + 1.0 / rr**2

    assert pair.pole(v, DEFAULT_RES_WINDOW) is None


def test_bound_state_is_accepted_through_bound_window(pair):
    from projects.potential_factory.tracker import DEFAULT_BOUND_WINDOW

    # N2 at R = 3.0 bohr: the anion is bound (Gamma == 0, shift < 0).
    p = pair.pole(_v_n2(3.0), DEFAULT_BOUND_WINDOW)
    assert p is not None and p.gamma == 0.0 and p.shift < N2.v0(3.0).real
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run --no-sync pytest projects/potential_factory/test_tracker.py -q`
Expected: FAIL with `ModuleNotFoundError` for `tracker`

- [ ] **Step 3: Implement**

```python
# projects/potential_factory/tracker.py
"""Fixed-R resonance-pole solving with a spurious-pole gate, the c-product
Hellmann-Feynman gradient, Newton in (lam, alpha), and continuation over R."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.core.grids import electronic_grid
from qscat.dvr import FemDvrEcsGrid, eigen, kinetic
from qscat.ecs import find_resonance_pole

__all__ = ["Pole", "Window", "ElectronicPair", "DEFAULT_RES_WINDOW", "DEFAULT_BOUND_WINDOW"]

Window = tuple[float, float, float, float]
VFn = Callable[[npt.NDArray[np.complex128]], npt.NDArray[np.complex128]]

DEFAULT_RES_WINDOW: Window = (0.006, 0.8, -0.6, 0.0)
DEFAULT_BOUND_WINDOW: Window = (-0.5, -1e-6, -1e-6, 1e-6)
_BOUND_IM_TOL = 1e-6


@dataclass(frozen=True)
class Pole:
    energy: complex
    residual: float

    @property
    def gamma(self) -> float:
        return max(0.0, -2.0 * self.energy.imag)

    @property
    def shift(self) -> float:
        return float(self.energy.real)


class ElectronicPair:
    """Two electronic FEM-DVR-ECS grids differing only in the ECS angle, with
    their kinetic matrices precomputed once."""

    def __init__(
        self,
        angles: tuple[float, float] = (35.0, 44.0),
        r_max: float = 16.0,
        order: int = 7,
        n_complex: int = 6,
    ) -> None:
        self.grid_a: FemDvrEcsGrid = electronic_grid(r_max=r_max, order=order, n_complex=n_complex, angle_deg=angles[0])
        self.grid_b: FemDvrEcsGrid = electronic_grid(r_max=r_max, order=order, n_complex=n_complex, angle_deg=angles[1])
        self._Ta = kinetic(self.grid_a, 1.0)
        self._Tb = kinetic(self.grid_b, 1.0)

    def hamiltonians(self, v_fn: VFn) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
        return (
            self._Ta + np.diag(v_fn(self.grid_a.points)),
            self._Tb + np.diag(v_fn(self.grid_b.points)),
        )

    def pole(
        self,
        v_fn: VFn,
        window: Window,
        *,
        resid_tol: float = 1e-3,
        gate_frac: float = 0.05,
        e_floor: float = 0.006,
    ) -> Pole | None:
        """The angle-stable pole in `window`, or None if no candidate passes the gate."""
        Ha, Hb = self.hamiltonians(v_fn)
        try:
            E, resid = find_resonance_pole(eigen(Ha)[0], eigen(Hb)[0], window)
        except ValueError:
            return None
        p = Pole(complex(E), float(resid))
        if p.residual >= resid_tol:
            return None
        bound = p.energy.real < 0.0 and abs(p.energy.imag) < _BOUND_IM_TOL
        if bound:
            return p
        if p.gamma <= 0.0 or p.residual >= gate_frac * p.gamma:
            return None
        if p.energy.real < e_floor:
            return None
        return p
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run --no-sync pytest projects/potential_factory/test_tracker.py -q`
Expected: 3 passed (≈20 s; the dense eigensolves dominate)

- [ ] **Step 5: Commit**

```bash
uv run ruff check projects/potential_factory && uv run ruff format projects/potential_factory
git add projects/potential_factory/tracker.py projects/potential_factory/test_tracker.py
git commit -m "feat(factory): ElectronicPair.pole — angle-stable pole with the residual<<Gamma / threshold-floor gate"
```

---

### Task 4: The c-product Hellmann–Feynman gradient

**Files:**
- Modify: `projects/potential_factory/tracker.py`
- Test: `projects/potential_factory/test_tracker.py` (append)

**Interfaces:**
- Consumes: `qscat.linalg.c_product`.
- Produces: `pole_vector(H, E) -> ndarray[complex128]` (the right eigenvector of the eigenvalue nearest `E`, c-normalised `ψᵀψ = 1`); `pole_sensitivity(H, E) -> ndarray[complex128]` = `ψ²` (so `dE/dV_i = psi_i**2`, exact for a complex-symmetric `H`).

- [ ] **Step 1: Write the failing test**

```python
# append to projects/potential_factory/test_tracker.py
from projects.potential_factory.tracker import pole_sensitivity  # noqa: E402


def test_c_product_gradient_matches_finite_difference(pair):
    Ha, _ = pair.hamiltonians(_v_n2(N2.R0))
    p = pair.pole(_v_n2(N2.R0), DEFAULT_RES_WINDOW)
    assert p is not None
    dEdV = pole_sensitivity(Ha, p.energy)
    # perturb one real-region diagonal entry and re-find the nearest eigenvalue
    i = int(np.argmin(np.abs(pair.grid_a.points - 1.5)))
    h = 1e-5
    Hp = Ha.copy()
    Hp[i, i] += h
    Ep = eigen_nearest(Hp, p.energy)
    fd = (Ep - p.energy) / h
    assert abs(fd - dEdV[i]) < 1e-3 * max(1.0, abs(dEdV[i]))


def eigen_nearest(H, E):
    from qscat.dvr import eigen

    vals = eigen(H)[0]
    return vals[int(np.argmin(np.abs(vals - E)))]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run --no-sync pytest projects/potential_factory/test_tracker.py -q -k gradient`
Expected: FAIL with `ImportError: cannot import name 'pole_sensitivity'`

- [ ] **Step 3: Implement**

```python
# append to projects/potential_factory/tracker.py
from qscat.dvr import eigen as _eigen  # noqa: E402  (fold into the import block)
from qscat.linalg import c_product  # noqa: E402


def pole_vector(H: npt.NDArray[np.complex128], E: complex) -> npt.NDArray[np.complex128]:
    """Right eigenvector of the eigenvalue of `H` nearest `E`, c-normalised (psi^T psi = 1)."""
    vals, vecs = _eigen(H)
    j = int(np.argmin(np.abs(vals - E)))
    psi = np.asarray(vecs[:, j], dtype=np.complex128)
    norm = c_product(psi, psi)
    if abs(norm) < 1e-12:
        raise ValueError("eigenvector is c-product self-orthogonal; cannot normalise")
    return psi / np.sqrt(norm)


def pole_sensitivity(H: npt.NDArray[np.complex128], E: complex) -> npt.NDArray[np.complex128]:
    """`dE/dV_i = psi_i^2` for complex-symmetric `H` (Hellmann-Feynman under the c-product)."""
    psi = pole_vector(H, E)
    return psi * psi
```

Add both names to `__all__`.

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run --no-sync pytest projects/potential_factory/test_tracker.py -q -k gradient`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
uv run ruff check projects/potential_factory && uv run ruff format projects/potential_factory
git add projects/potential_factory/tracker.py projects/potential_factory/test_tracker.py
git commit -m "feat(factory): pole_sensitivity — exact dE/dV from one c-normalised eigenvector"
```

---

### Task 5: Newton in `(λ, α)` for a target pole

**Files:**
- Modify: `projects/potential_factory/tracker.py`
- Test: `projects/potential_factory/test_tracker.py` (append)

**Interfaces:**
- Produces:
  - `WellParams(lam: float, alpha: float)` (frozen dataclass)
  - `well_potential(ell, lam, alpha, extra: VFn | None) -> VFn` — `-lam exp(-alpha r^2) + ell(ell+1)/2r^2 (+ extra(r))`
  - `solve_pole_params(pair, ell, target: complex, seed: WellParams, *, extra=None, window=None, max_iter=30, tol=1e-8) -> tuple[WellParams, Pole]`; raises `ConvergenceError` (from `qscat.exceptions`) if Newton does not converge. For a **bound** target (`target.imag == 0 and target.real < 0`) only `lam` is solved and `alpha` is held at the seed.

Newton step: with `psi` c-normalised on `grid_a`, `dE/dlam = sum_i psi_i^2 * (-exp(-alpha r_i^2))`, `dE/dalpha = sum_i psi_i^2 * (lam r_i^2 exp(-alpha r_i^2))`; solve the real 2×2 system `[[Re J_lam, Re J_alpha],[Im J_lam, Im J_alpha]] dx = -[Re f, Im f]` with `f = E - target`; damp by halving while `|f|` does not decrease; recentre the window on the current pole (`±0.05` Ha real, `[-0.05, 0]` imaginary; for bound: `DEFAULT_BOUND_WINDOW`).

- [ ] **Step 1: Write the failing tests**

```python
# append to projects/potential_factory/test_tracker.py
from projects.potential_factory.tracker import WellParams, solve_pole_params, well_potential  # noqa: E402


def test_newton_recovers_n2_well_parameters_from_a_perturbed_seed(pair):
    R = N2.R0
    lam_true = float(N2.lam(R).real)
    v_true = well_potential(N2.ell, lam_true, N2.alpha_c, None)
    target = pair.pole(v_true, DEFAULT_RES_WINDOW)
    assert target is not None
    sol, pole = solve_pole_params(pair, N2.ell, target.energy, WellParams(lam=4.5, alpha=0.5))
    assert abs(sol.lam - lam_true) < 1e-5
    assert abs(sol.alpha - N2.alpha_c) < 1e-5
    assert abs(pole.energy - target.energy) < 1e-8


def test_newton_bound_target_solves_lam_only(pair):
    R = 3.0
    lam_true = float(N2.lam(R).real)
    v_true = well_potential(N2.ell, lam_true, N2.alpha_c, None)
    from projects.potential_factory.tracker import DEFAULT_BOUND_WINDOW

    target = pair.pole(v_true, DEFAULT_BOUND_WINDOW)
    assert target is not None
    sol, _ = solve_pole_params(pair, N2.ell, complex(target.energy.real, 0.0), WellParams(lam=lam_true * 0.9, alpha=N2.alpha_c))
    assert abs(sol.lam - lam_true) < 1e-5 and sol.alpha == N2.alpha_c
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run --no-sync pytest projects/potential_factory/test_tracker.py -q -k newton`
Expected: FAIL with `ImportError`

- [ ] **Step 3: Implement**

```python
# append to projects/potential_factory/tracker.py
from qscat.exceptions import ConvergenceError  # noqa: E402  (fold into imports)


@dataclass(frozen=True)
class WellParams:
    lam: float
    alpha: float


def well_potential(ell: int, lam: float, alpha: float, extra: VFn | None) -> VFn:
    def v(r: npt.NDArray[np.complex128]) -> npt.NDArray[np.complex128]:
        rr = np.asarray(r, dtype=np.complex128)
        out = -lam * np.exp(-alpha * rr**2) + ell * (ell + 1) / (2.0 * rr**2)
        if extra is not None:
            out = out + extra(rr)
        return np.asarray(out, dtype=np.complex128)

    return v


def _is_bound_target(target: complex) -> bool:
    return target.imag == 0.0 and target.real < 0.0


def _window_about(E: complex, bound: bool) -> Window:
    if bound:
        return DEFAULT_BOUND_WINDOW
    return (E.real - 0.05, E.real + 0.05, min(E.imag - 0.05, -0.05), 0.0)


def solve_pole_params(
    pair: ElectronicPair,
    ell: int,
    target: complex,
    seed: WellParams,
    *,
    extra: VFn | None = None,
    window: Window | None = None,
    max_iter: int = 30,
    tol: float = 1e-8,
) -> tuple[WellParams, Pole]:
    """Newton on (lam, alpha) so that the gated pole equals `target`."""
    bound = _is_bound_target(target)
    lam, alpha = seed.lam, seed.alpha
    r = pair.grid_a.points
    win = window if window is not None else _window_about(target, bound)
    pole = pair.pole(well_potential(ell, lam, alpha, extra), win)
    if pole is None:
        raise ConvergenceError("no gated pole at the seed parameters")
    for _ in range(max_iter):
        f = pole.energy - target
        if abs(f) < tol:
            return WellParams(lam, alpha), pole
        Ha, _ = pair.hamiltonians(well_potential(ell, lam, alpha, extra))
        s = pole_sensitivity(Ha, pole.energy)
        g = np.exp(-alpha * r**2)
        dlam = complex(np.sum(s * (-g)))
        dalpha = complex(np.sum(s * (lam * r**2 * g)))
        if bound:
            step = np.array([-f.real / dlam.real, 0.0])
        else:
            J = np.array([[dlam.real, dalpha.real], [dlam.imag, dalpha.imag]])
            step = np.linalg.solve(J, -np.array([f.real, f.imag]))
        damp = 1.0
        while damp > 1e-4:
            lam_n, alpha_n = lam + damp * step[0], alpha + damp * step[1]
            if lam_n <= 0.0 or alpha_n <= 0.0:
                damp *= 0.5
                continue
            cand = pair.pole(well_potential(ell, lam_n, alpha_n, extra), _window_about(pole.energy, bound))
            if cand is not None and abs(cand.energy - target) < abs(f):
                lam, alpha, pole = lam_n, alpha_n, cand
                break
            damp *= 0.5
        else:
            raise ConvergenceError(f"Newton stalled at lam={lam:.6g}, alpha={alpha:.6g}, |f|={abs(f):.3g}")
    raise ConvergenceError(f"Newton did not converge in {max_iter} iterations")
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run --no-sync pytest projects/potential_factory/test_tracker.py -q -k newton`
Expected: 2 passed

- [ ] **Step 5: Commit**

```bash
uv run ruff check projects/potential_factory && uv run ruff format projects/potential_factory
git add projects/potential_factory/tracker.py projects/potential_factory/test_tracker.py
git commit -m "feat(factory): solve_pole_params — damped Newton in (lam, alpha) on the c-product gradient"
```

---

### Task 6: Continuation over `R` (`track_curve`)

**Files:**
- Modify: `projects/potential_factory/tracker.py`
- Test: `projects/potential_factory/test_tracker.py` (append)

**Interfaces:**
- Produces: `TrackResult(R: ndarray, lam: ndarray, alpha: ndarray, converged: ndarray[bool], poles: list[Pole | None])` and `track_curve(pair, ell, R_desc, target_of_R: Callable[[float], complex], seed: WellParams, *, extra_of_R: Callable[[float], VFn | None] | None = None) -> TrackResult`. `R_desc` strictly descending (raise `ValueError` otherwise). Each node is seeded from the previous accepted node; on `ConvergenceError` the node is flagged `converged=False`, its `(lam, alpha)` copied from the last accepted node (freeze-and-flag, mirroring `resonance_pole_walk`), and tracking continues.

- [ ] **Step 1: Write the failing test**

```python
# append to projects/potential_factory/test_tracker.py
from projects.potential_factory.tracker import track_curve  # noqa: E402


@pytest.mark.slow
def test_track_curve_recovers_n2_lam_and_alpha_over_R(pair):
    R_desc = np.linspace(3.0, 1.6, 15)

    def target(R):
        v = well_potential(N2.ell, float(N2.lam(R).real), N2.alpha_c, None)
        p = pair.pole(v, DEFAULT_BOUND_WINDOW) if R > 2.5 else pair.pole(v, DEFAULT_RES_WINDOW)
        assert p is not None
        return complex(p.energy.real, 0.0) if p.gamma == 0.0 else p.energy

    res = track_curve(pair, N2.ell, R_desc, target, WellParams(lam=float(N2.lam(3.0).real) * 0.95, alpha=N2.alpha_c))
    assert res.converged.all()
    np.testing.assert_allclose(res.lam, N2.lam(R_desc).real, rtol=1e-4)
    np.testing.assert_allclose(res.alpha, N2.alpha_c, rtol=1e-4)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run --no-sync pytest projects/potential_factory/test_tracker.py -q -k track_curve`
Expected: FAIL with `ImportError`

- [ ] **Step 3: Implement**

```python
# append to projects/potential_factory/tracker.py
@dataclass(frozen=True)
class TrackResult:
    R: npt.NDArray[np.float64]
    lam: npt.NDArray[np.float64]
    alpha: npt.NDArray[np.float64]
    converged: npt.NDArray[np.bool_]
    poles: list[Pole | None]


def track_curve(
    pair: ElectronicPair,
    ell: int,
    R_desc: npt.ArrayLike,
    target_of_R: Callable[[float], complex],
    seed: WellParams,
    *,
    extra_of_R: Callable[[float], VFn | None] | None = None,
) -> TrackResult:
    R = np.asarray(R_desc, dtype=np.float64)
    if R.size > 1 and np.any(np.diff(R) >= 0.0):
        raise ValueError("R_desc must be strictly descending")
    lam = np.empty(R.size)
    alpha = np.empty(R.size)
    ok = np.zeros(R.size, dtype=bool)
    poles: list[Pole | None] = []
    cur = seed
    for j, Rj in enumerate(R):
        extra = extra_of_R(float(Rj)) if extra_of_R is not None else None
        try:
            cur, pole = solve_pole_params(pair, ell, target_of_R(float(Rj)), cur, extra=extra)
            ok[j] = True
            poles.append(pole)
        except ConvergenceError:
            poles.append(None)
        lam[j], alpha[j] = cur.lam, cur.alpha
    return TrackResult(R=R, lam=lam, alpha=alpha, converged=ok, poles=poles)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run --no-sync pytest projects/potential_factory/test_tracker.py -q -k track_curve -m slow`
Expected: PASS (≈2–4 min)

- [ ] **Step 5: Commit**

```bash
uv run ruff check projects/potential_factory && uv run ruff format projects/potential_factory
git add projects/potential_factory/tracker.py projects/potential_factory/test_tracker.py
git commit -m "feat(factory): track_curve — (lam, alpha) continuation over R with freeze-and-flag"
```

---

### Task 7: `Target` and `extract_target` (the round-trip data path)

**Files:**
- Create: `projects/potential_factory/target.py`
- Create: `projects/potential_factory/extract.py`
- Test: `projects/potential_factory/test_target.py`

**Interfaces:**
- Produces (`target.py`):
  - `Curve` — `Curve.from_table(x, y)` (cubic `scipy.interpolate.CubicSpline`, `extrapolate=False` → NaN outside) or `Curve.from_callable(fn)`; `__call__(x) -> ndarray[float64]`; `.x` (table nodes or `None`)
  - `NeutralTarget(curve: Curve | None, constants: dict[str, float], R_range: tuple[float, float])` — constants may hold `R_e, D_e, omega_e, omega_e_x_e`
  - `ResonanceTarget(v_ion: Curve, gamma: Curve, ea: float, R_range)` — `v_ion(R)` absolute (includes `v0`), `gamma(R) ≥ 0`, `ea` = fragment electron affinity (Ha): asymptote `v_ion(∞) − v0(∞) = −ea`
  - `CouplingTarget(gamma_tilde: Callable[[eps, R], ndarray], eps_window: tuple[float, float], R_range, alpha_exponent: float)` with `CouplingTarget.from_alt_houfek(a0, a1, a2, b0, b1, alpha) -> CouplingTarget` implementing `2π ε^α (a0 + a1 R) e^{a2 R} e^{−(b0 + b1 R) ε}` and `CouplingTarget.from_table(eps, R, gamma_table)` (bilinear `RegularGridInterpolator`)
  - `Provenance(source: str, locator: str)`; `Target(name, mu, ell, charge, coordinates: tuple[str, ...], neutral, resonance, coupling, eigenphase=None, provenance: dict[str, Provenance])`
- Produces (`extract.py`): `extract_target(model, *, pair: ElectronicPair, R_desc, R_inf=10.0, eps_window=(0.002, 0.25), n_eps=12, name="model") -> Target` using `qscat.core.lcp.resonance_pole_walk` (T1), `qscat.core.dissociation.anion_electronic_states` (EA), `AsymptoticDiscreteState` + `v_dk_plus` + `gamma_from_coupling` (T3 table).

- [ ] **Step 1: Write the failing tests**

```python
# projects/potential_factory/test_target.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.core.dissociation import anion_electronic_states
from qscat.model import N2

from projects.potential_factory.extract import extract_target
from projects.potential_factory.target import CouplingTarget, Curve
from projects.potential_factory.tracker import ElectronicPair


def test_curve_from_table_interpolates_and_is_nan_outside():
    c = Curve.from_table(np.array([1.0, 2.0, 3.0]), np.array([1.0, 4.0, 9.0]))
    assert abs(c(2.5) - 6.25) < 0.05
    assert np.isnan(c(5.0))


def test_alt_houfek_coupling_has_the_threshold_exponent():
    ct = CouplingTarget.from_alt_houfek(a0=13.83669, a1=0.892095, a2=-0.935987, b0=3.015014, b1=0.718160, alpha=2.5, R_range=(1.8, 2.8))
    e1, e2 = 1e-4, 2e-4
    slope = np.log(ct.gamma_tilde(e2, 2.0) / ct.gamma_tilde(e1, 2.0)) / np.log(e2 / e1)
    assert abs(slope - 2.5) < 1e-3


@pytest.mark.slow
def test_extract_target_from_n2_is_self_consistent():
    pair = ElectronicPair()
    R_desc = np.linspace(3.0, 1.6, 8)
    t = extract_target(N2, pair=pair, R_desc=R_desc, n_eps=4)
    assert t.ell == 2 and t.coordinates == ("R",)
    assert t.neutral is not None and t.resonance is not None and t.coupling is not None
    np.testing.assert_allclose(t.neutral.curve(R_desc), N2.v0(R_desc).real, atol=1e-12)
    eps_e, _ = anion_electronic_states(pair.grid_a, N2, 10.0, 1)
    assert abs(-t.resonance.ea - (eps_e[0] - N2.v0(10.0).real)) < 1e-10
    g = t.resonance.gamma(R_desc)
    assert np.all(g >= 0.0) and g[-1] > 0.0 and g[0] == 0.0
    # threshold law of the extracted coupling: slope ~ l + 1/2 at small eps
    e1, e2 = 0.002, 0.004
    slope = np.log(t.coupling.gamma_tilde(e2, 1.9) / t.coupling.gamma_tilde(e1, 1.9)) / np.log(e2 / e1)
    assert abs(slope - 2.5) < 0.3
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run --no-sync pytest projects/potential_factory/test_target.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement `target.py`**

```python
# projects/potential_factory/target.py
"""What a molecule IS to the factory: tiered, provenance-carrying target data."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
from scipy.interpolate import CubicSpline, RegularGridInterpolator

__all__ = ["Curve", "NeutralTarget", "ResonanceTarget", "CouplingTarget", "Provenance", "Target"]

FArr = npt.NDArray[np.float64]


class Curve:
    """A real curve on one coordinate: a cubic spline through a table, or a callable."""

    def __init__(self, fn: Callable[[FArr], FArr], x: FArr | None) -> None:
        self._fn = fn
        self.x = x

    @classmethod
    def from_table(cls, x: npt.ArrayLike, y: npt.ArrayLike) -> Curve:
        xs = np.asarray(x, dtype=np.float64)
        ys = np.asarray(y, dtype=np.float64)
        order = np.argsort(xs)
        spl = CubicSpline(xs[order], ys[order], extrapolate=False)
        return cls(lambda t: np.asarray(spl(t), dtype=np.float64), xs[order])

    @classmethod
    def from_callable(cls, fn: Callable[[FArr], FArr]) -> Curve:
        return cls(lambda t: np.asarray(fn(np.asarray(t, dtype=np.float64)), dtype=np.float64), None)

    def __call__(self, x: npt.ArrayLike) -> FArr:
        return self._fn(np.asarray(x, dtype=np.float64))


@dataclass(frozen=True)
class NeutralTarget:
    curve: Curve | None
    constants: dict[str, float]
    R_range: tuple[float, float]


@dataclass(frozen=True)
class ResonanceTarget:
    v_ion: Curve
    gamma: Curve
    ea: float
    R_range: tuple[float, float]


@dataclass(frozen=True)
class CouplingTarget:
    gamma_tilde: Callable[[npt.ArrayLike, npt.ArrayLike], FArr]
    eps_window: tuple[float, float]
    R_range: tuple[float, float]
    alpha_exponent: float

    @classmethod
    def from_alt_houfek(
        cls, *, a0: float, a1: float, a2: float, b0: float, b1: float, alpha: float,
        R_range: tuple[float, float], eps_window: tuple[float, float] = (0.002, 0.22),
    ) -> CouplingTarget:
        """Alt & Houfek, PRA 103, 032829 (2021) Eq. (25)-(27): 2 pi eps^alpha A(R) exp(-B(R) eps)."""

        def g(eps: npt.ArrayLike, R: npt.ArrayLike) -> FArr:
            e = np.asarray(eps, dtype=np.float64)
            r = np.asarray(R, dtype=np.float64)
            A = (a0 + a1 * r) * np.exp(a2 * r)
            B = b0 + b1 * r
            return np.asarray(2.0 * np.pi * e**alpha * A * np.exp(-B * e), dtype=np.float64)

        return cls(g, eps_window, R_range, alpha)

    @classmethod
    def from_table(
        cls, eps: npt.ArrayLike, R: npt.ArrayLike, table: npt.ArrayLike, *, alpha: float,
    ) -> CouplingTarget:
        e = np.asarray(eps, dtype=np.float64)
        r = np.asarray(R, dtype=np.float64)
        tab = np.asarray(table, dtype=np.float64)  # shape (e.size, r.size)
        interp = RegularGridInterpolator((e, r), tab, bounds_error=False, fill_value=np.nan)

        def g(eps_q: npt.ArrayLike, R_q: npt.ArrayLike) -> FArr:
            eq = np.asarray(eps_q, dtype=np.float64)
            rq = np.asarray(R_q, dtype=np.float64)
            eb, rb = np.broadcast_arrays(eq, rq)
            pts = np.stack([eb.ravel(), rb.ravel()], axis=-1)
            return np.asarray(interp(pts).reshape(eb.shape), dtype=np.float64)

        return cls(g, (float(e.min()), float(e.max())), (float(r.min()), float(r.max())), alpha)


@dataclass(frozen=True)
class Provenance:
    source: str
    locator: str


@dataclass(frozen=True)
class Target:
    name: str
    mu: float
    ell: int
    charge: int
    coordinates: tuple[str, ...]
    neutral: NeutralTarget | None
    resonance: ResonanceTarget | None
    coupling: CouplingTarget | None
    eigenphase: object | None = None  # reserved: T2 tables (loader only, no loss in v1)
    provenance: dict[str, Provenance] = field(default_factory=dict)
```

- [ ] **Step 4: Implement `extract.py`**

```python
# projects/potential_factory/extract.py
"""Build a Target FROM an existing model with the repo's own forward models.

This is the round-trip data path: the extracted curves are exact properties of
the model, so fitting them back must recover the model."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.dissociation import anion_electronic_states
from qscat.core.lcp import resonance_pole_walk
from qscat.core.nrm.coupling import gamma_from_coupling, v_dk_plus
from qscat.core.nrm.discrete_state import AsymptoticDiscreteState
from qscat.model import ResonanceModel

from projects.potential_factory.target import CouplingTarget, Curve, NeutralTarget, Provenance, ResonanceTarget, Target
from projects.potential_factory.tracker import ElectronicPair

__all__ = ["extract_target"]


def extract_target(
    model: ResonanceModel,
    *,
    pair: ElectronicPair,
    R_desc: npt.ArrayLike,
    R_inf: float = 10.0,
    eps_window: tuple[float, float] = (0.002, 0.25),
    n_eps: int = 12,
    name: str = "model",
) -> Target:
    R = np.asarray(R_desc, dtype=np.float64)
    if R.size > 1 and np.any(np.diff(R) >= 0.0):
        raise ValueError("R_desc must be strictly descending")
    v0 = model.v0(R).real

    # T1: the two-angle pole walk, seeded on the bound anion at the outermost R.
    eps_e, _ = anion_electronic_states(pair.grid_a, model, float(R[0]), 1)
    seed = (eps_e[0] - 0.05, eps_e[0] + 0.05, -0.05, 1e-6)
    shift, gamma = resonance_pole_walk(model, R, pair.grid_a, pair.grid_b, seed)
    v_ion = v0 + shift
    eps_inf, _ = anion_electronic_states(pair.grid_a, model, R_inf, 1)
    ea = -(eps_inf[0] - model.v0(R_inf).real)

    # T3: the model's own Gamma~(eps, R) with the R-independent discrete state.
    phi_d = AsymptoticDiscreteState(pair.grid_a, model, R_inf)
    eps = np.geomspace(eps_window[0], eps_window[1], n_eps)
    table = np.empty((eps.size, R.size))
    R_asc = R[::-1]
    for i, e in enumerate(eps):
        table[i] = gamma_from_coupling(v_dk_plus(pair.grid_a, model, phi_d, R_asc, float(e)))

    rng = (float(R.min()), float(R.max()))
    return Target(
        name=name,
        mu=model.mu,
        ell=model.ell,
        charge=model.charge,
        coordinates=("R",),
        neutral=NeutralTarget(Curve.from_table(R, v0), {}, rng),
        resonance=ResonanceTarget(Curve.from_table(R, v_ion), Curve.from_table(R, gamma), float(ea), rng),
        coupling=CouplingTarget.from_table(eps, R_asc, table, alpha=model.ell + 0.5),
        provenance={"all": Provenance(f"extract_target({name})", "computed, not published")},
    )
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run --no-sync pytest projects/potential_factory/test_target.py -q -m "not slow"` then `-m slow`
Expected: 2 passed; then 1 passed (≈1–2 min)

- [ ] **Step 6: Commit**

```bash
uv run ruff check projects/potential_factory && uv run ruff format projects/potential_factory
git add projects/potential_factory/target.py projects/potential_factory/extract.py projects/potential_factory/test_target.py
git commit -m "feat(factory): tiered Target + extract_target (the round-trip data path)"
```

---

### Task 8: `fit_neutral` — EMO least squares (T0)

**Files:**
- Create: `projects/potential_factory/fit.py`
- Create: `projects/potential_factory/report.py`
- Test: `projects/potential_factory/test_fit.py`

**Interfaces:**
- Produces (`report.py`): `TierResult(name: str, status: Literal["met", "not met", "not attempted"], rms: float, max: float, detail: str)`; `Tolerances(v0_rms=2e-4, omega_e_rel=0.01, e_res_rms=1e-3, gamma_rel=0.10, gamma_floor=2e-3, coupling_log_rms=0.2)` (Hartree / relative; **placeholders until the sensitivity budget plan replaces them**).
- Produces (`fit.py`): `fit_neutral(target: NeutralTarget, seed: FlexibleDiatomicModel, *, n_beta: int = 1, tol: Tolerances) -> tuple[FlexibleDiatomicModel, TierResult]` — frees `D_e, R_e, beta0..beta{n_beta-1}` and least-squares the EMO to `target.curve` on 200 points over `R_range` (or to the constants when `curve is None`: `R_e`, `D_e` direct; `beta0 = omega_e * sqrt(mu / (2 D_e))` — the Morse relation `ω_e = β √(2D_e/μ)`).

- [ ] **Step 1: Write the failing test**

```python
# projects/potential_factory/test_fit.py
from __future__ import annotations

import numpy as np
from qscat.model import N2

from projects.potential_factory.ansatz import from_diatomic, with_params
from projects.potential_factory.fit import fit_neutral
from projects.potential_factory.report import Tolerances
from projects.potential_factory.target import Curve, NeutralTarget


def test_fit_neutral_recovers_morse_constants_from_the_curve():
    R = np.linspace(1.4, 5.0, 60)
    target = NeutralTarget(Curve.from_table(R, N2.v0(R).real), {}, (1.5, 4.8))
    seed = with_params(from_diatomic(N2), {"D_e": 0.6, "R_e": 2.2, "beta0": 1.0})
    model, res = fit_neutral(target, seed, tol=Tolerances())
    assert res.status == "met"
    assert abs(model.D_e - N2.D0) < 1e-8 and abs(model.R_e - N2.R0) < 1e-8 and abs(model.betas[0] - N2.alpha0) < 1e-8


def test_fit_neutral_from_constants_uses_the_morse_relation():
    omega_e = N2.alpha0 * np.sqrt(2.0 * N2.D0 / N2.mu)
    target = NeutralTarget(None, {"R_e": N2.R0, "D_e": N2.D0, "omega_e": omega_e}, (1.5, 4.8))
    model, res = fit_neutral(target, from_diatomic(N2), tol=Tolerances())
    assert res.status == "met" and abs(model.betas[0] - N2.alpha0) < 1e-12
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run --no-sync pytest projects/potential_factory/test_fit.py -q`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement `report.py` and `fit_neutral`**

```python
# projects/potential_factory/report.py
"""FitReport: the model, per-tier residuals, verdicts, and provenance."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Literal

__all__ = ["TierResult", "Tolerances", "FitReport"]

Status = Literal["met", "not met", "not attempted"]


@dataclass(frozen=True)
class TierResult:
    name: str
    status: Status
    rms: float
    max: float
    detail: str


@dataclass(frozen=True)
class Tolerances:
    """PLACEHOLDERS (Hartree / relative) until validation/factory's sensitivity
    budget replaces them -- see docs/superpowers/plans/2026-08-24-potential-factory-budget.md."""

    v0_rms: float = 2e-4
    omega_e_rel: float = 0.01
    e_res_rms: float = 1e-3
    gamma_rel: float = 0.10
    gamma_floor: float = 2e-3
    coupling_log_rms: float = 0.2


@dataclass
class FitReport:
    target_name: str
    parameters: dict[str, float]
    tiers: list[TierResult]
    ecs_bounds_deg: dict[str, float]
    crossing_R: float | None
    da_threshold_sign: int | None
    provenance: dict[str, dict[str, str]] = field(default_factory=dict)

    def to_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2))

    @classmethod
    def from_json(cls, path: str | Path) -> FitReport:
        d = json.loads(Path(path).read_text())
        d["tiers"] = [TierResult(**t) for t in d["tiers"]]
        return cls(**d)
```

```python
# projects/potential_factory/fit.py
"""Staged fitting: T0 (neutral) -> T1 (pole curves) -> T3 (energy-dependent width).
Each stage seeds the next; a stage that misses its tolerance stops and reports."""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from projects.potential_factory.ansatz import FlexibleDiatomicModel, pack, unpack, with_params
from projects.potential_factory.report import TierResult, Tolerances
from projects.potential_factory.target import NeutralTarget

__all__ = ["fit_neutral"]


def _with_n_beta(model: FlexibleDiatomicModel, n_beta: int) -> FlexibleDiatomicModel:
    """Grow (with zeros) or truncate the EMO beta expansion to `n_beta` terms."""
    betas = list(model.betas) + [0.0] * max(0, n_beta - len(model.betas))
    return replace(model, betas=tuple(betas[:n_beta]))


def fit_neutral(
    target: NeutralTarget,
    seed: FlexibleDiatomicModel,
    *,
    n_beta: int = 1,
    tol: Tolerances,
) -> tuple[FlexibleDiatomicModel, TierResult]:
    model = _with_n_beta(seed, n_beta)
    if target.curve is None:
        c = target.constants
        beta0 = c["omega_e"] * np.sqrt(model.mu / (2.0 * c["D_e"]))
        model = with_params(model, {"D_e": c["D_e"], "R_e": c["R_e"], "beta0": float(beta0)})
        return model, TierResult("T0", "met", 0.0, 0.0, "constants only; Morse relation for beta0")
    R = np.linspace(target.R_range[0], target.R_range[1], 200)
    y = target.curve(R)
    names = ["D_e", "R_e"] + [f"beta{i}" for i in range(n_beta)]

    def resid(x):
        m = unpack(model, names, x)
        return m.v0(R).real - y

    sol = least_squares(resid, pack(model, names), xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=2000)
    fitted = unpack(model, names, sol.x)
    r = resid(sol.x)
    rms, mx = float(np.sqrt(np.mean(r**2))), float(np.max(np.abs(r)))
    status = "met" if rms <= tol.v0_rms else "not met"
    return fitted, TierResult("T0", status, rms, mx, f"EMO with {n_beta} beta(s) on {R.size} points")
```

(`replace` is `dataclasses.replace`; import it at the top of `fit.py`.)

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run --no-sync pytest projects/potential_factory/test_fit.py -q`
Expected: 2 passed

- [ ] **Step 5: Commit**

```bash
uv run ruff check projects/potential_factory && uv run ruff format projects/potential_factory
git add projects/potential_factory/fit.py projects/potential_factory/report.py projects/potential_factory/test_fit.py
git commit -m "feat(factory): fit_neutral (EMO least squares) + TierResult/Tolerances/FitReport"
```

---

### Task 9: `fit_resonance` — track, smooth, constrain, re-verify (T1)

**Files:**
- Modify: `projects/potential_factory/fit.py`
- Test: `projects/potential_factory/test_fit.py` (append)

**Interfaces:**
- Consumes: `track_curve`, `WellParams`, `ElectronicPair`, `resonance_pole_walk` (re-verification), `anion_electronic_states` (EA constraint).
- Produces: `fit_resonance(target: ResonanceTarget, model: FlexibleDiatomicModel, *, pair: ElectronicPair, n_nodes: int = 24, tol: Tolerances, lam_coeffs: int = 0, alpha_coeffs: int = 0) -> tuple[FlexibleDiatomicModel, TierResult]`. Steps: (1) nodes `R_desc = linspace(R_range[1], R_range[0], n_nodes)`; target pole per node `E = (v_ion − v0_model)(R) − i γ(R)/2` (γ==0 → bound target); (2) `track_curve` seeded from `(model.lam_R(R_max), model.alpha_R(R_max))`; (3) least-squares the `lam`/`alpha` `SmoothR` parameters (`lam.f_inf, lam.f_0, lam.f_1, lam.R_f, lam.c*` and `alpha.f_inf, alpha.f_0, alpha.f_1, alpha.R_f, alpha.c*`, with `coeffs` grown to `lam_coeffs`/`alpha_coeffs` zeros) to the tracked `lam_j`, `alpha_j` on converged nodes; (4) EA constraint: `brentq` on `lam.f_inf` so that `anion_electronic_states(grid_a, model, R_inf=10)` gives `eps_e − v0(10) = −target.ea` (skip when `|v_ion(R_max) − v0(R_max) + ea| > 0.05` i.e. the table does not reach the asymptote — then report it); (5) re-verify with `resonance_pole_walk` on the smoothed model: residuals `rms/max` of `E_res` (Ha) and of `Γ` relative where `γ_target > tol.gamma_floor`; `status = "met"` iff both within tolerance.

- [ ] **Step 1: Write the failing test**

```python
# append to projects/potential_factory/test_fit.py
import pytest  # noqa: E402

from projects.potential_factory.extract import extract_target  # noqa: E402
from projects.potential_factory.fit import fit_resonance  # noqa: E402
from projects.potential_factory.tracker import ElectronicPair  # noqa: E402


@pytest.mark.slow
def test_fit_resonance_round_trips_n2_curves():
    pair = ElectronicPair()
    R_desc = np.linspace(3.0, 1.6, 12)
    target = extract_target(N2, pair=pair, R_desc=R_desc, n_eps=3)
    seed = with_params(from_diatomic(N2), {"lam.f_inf": 5.5, "alpha.f_inf": 0.5})
    model, res = fit_resonance(target.resonance, seed, pair=pair, n_nodes=12, tol=Tolerances())
    assert res.status == "met", res.detail
    np.testing.assert_allclose(model.lam_R(R_desc).real, N2.lam(R_desc).real, rtol=2e-3)
    np.testing.assert_allclose(model.alpha_R(R_desc).real, N2.alpha_c, rtol=2e-3)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run --no-sync pytest projects/potential_factory/test_fit.py -q -k resonance -m slow`
Expected: FAIL with `ImportError: cannot import name 'fit_resonance'`

- [ ] **Step 3: Implement**

```python
# append to projects/potential_factory/fit.py
from dataclasses import replace  # noqa: E402  (fold into imports)

from qscat.core.dissociation import anion_electronic_states  # noqa: E402
from qscat.core.lcp import resonance_pole_walk  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

from projects.potential_factory.ansatz import SmoothR  # noqa: E402
from projects.potential_factory.target import ResonanceTarget  # noqa: E402
from projects.potential_factory.tracker import ElectronicPair, WellParams, track_curve  # noqa: E402

_R_INF = 10.0


def _grow_coeffs(s: SmoothR, n: int) -> SmoothR:
    c = list(s.coeffs) + [0.0] * max(0, n - len(s.coeffs))
    return replace(s, coeffs=tuple(c[:n]))


def _smooth_names(prefix: str, s: SmoothR) -> list[str]:
    return [f"{prefix}.f_inf", f"{prefix}.f_0", f"{prefix}.f_1", f"{prefix}.R_f"] + [f"{prefix}.c{i}" for i in range(len(s.coeffs))]


def _fit_smooth(model: FlexibleDiatomicModel, prefix: str, R: np.ndarray, y: np.ndarray) -> FlexibleDiatomicModel:
    s = getattr(model, prefix)
    names = _smooth_names(prefix, s)

    def resid(x):
        m = unpack(model, names, x)
        return getattr(m, prefix)(R).real - y

    sol = least_squares(resid, pack(model, names), xtol=1e-12, ftol=1e-12, max_nfev=5000)
    return unpack(model, names, sol.x)


def _apply_ea_constraint(model: FlexibleDiatomicModel, pair: ElectronicPair, ea: float) -> FlexibleDiatomicModel:
    def g(f_inf: float) -> float:
        m = with_params(model, {"lam.f_inf": f_inf})
        eps_e, _ = anion_electronic_states(pair.grid_a, m, _R_INF, 1)
        return float(eps_e[0] - m.v0(_R_INF).real + ea)

    f0 = model.lam.f_inf
    lo, hi = 0.5 * f0, 1.5 * f0
    if g(lo) * g(hi) > 0:
        return model  # no bracket: leave lam.f_inf; the re-verification will report it
    return with_params(model, {"lam.f_inf": brentq(g, lo, hi, xtol=1e-10)})


def fit_resonance(
    target: ResonanceTarget,
    model: FlexibleDiatomicModel,
    *,
    pair: ElectronicPair,
    n_nodes: int = 24,
    tol: Tolerances,
    lam_coeffs: int = 0,
    alpha_coeffs: int = 0,
) -> tuple[FlexibleDiatomicModel, TierResult]:
    R_desc = np.linspace(target.R_range[1], target.R_range[0], n_nodes)
    model = replace(model, lam=_grow_coeffs(model.lam, lam_coeffs), alpha=_grow_coeffs(model.alpha, alpha_coeffs))

    def pole_target(R: float) -> complex:
        s = float(target.v_ion(R) - model.v0(R).real)
        g = float(target.gamma(R))
        return complex(s, 0.0) if g <= 0.0 else complex(s, -0.5 * g)

    seed = WellParams(float(model.lam_R(R_desc[0]).real), float(model.alpha_R(R_desc[0]).real))
    tr = track_curve(pair, model.ell, R_desc, pole_target, seed)
    ok = tr.converged
    if ok.sum() < 4:
        return model, TierResult("T1", "not met", np.inf, np.inf, f"only {int(ok.sum())} nodes tracked")
    model = _fit_smooth(model, "lam", tr.R[ok], tr.lam[ok])
    model = _fit_smooth(model, "alpha", tr.R[ok], tr.alpha[ok])
    reaches_asymptote = abs(float(target.v_ion(R_desc[0]) - model.v0(R_desc[0]).real) + target.ea) < 0.05
    if reaches_asymptote:
        model = _apply_ea_constraint(model, pair, target.ea)

    # re-verify on the SMOOTHED model -- this is where the residual lives
    eps_e, _ = anion_electronic_states(pair.grid_a, model, float(R_desc[0]), 1)
    seed_win = (eps_e[0] - 0.05, eps_e[0] + 0.05, -0.05, 1e-6)
    shift, gamma = resonance_pole_walk(model, R_desc, pair.grid_a, pair.grid_b, seed_win)
    e_err = (model.v0(R_desc).real + shift) - target.v_ion(R_desc)
    g_t = target.gamma(R_desc)
    mask = g_t > tol.gamma_floor
    g_rel = (gamma[mask] - g_t[mask]) / g_t[mask] if mask.any() else np.zeros(0)
    e_rms, e_max = float(np.sqrt(np.mean(e_err**2))), float(np.max(np.abs(e_err)))
    g_rms = float(np.sqrt(np.mean(g_rel**2))) if g_rel.size else 0.0
    g_max = float(np.max(np.abs(g_rel))) if g_rel.size else 0.0
    met = e_rms <= tol.e_res_rms and g_max <= tol.gamma_rel and bool(ok.all())
    detail = (
        f"tracked {int(ok.sum())}/{ok.size} nodes; E_res rms={e_rms:.2e} max={e_max:.2e} Ha; "
        f"Gamma rel rms={g_rms:.3f} max={g_max:.3f}; EA constraint {'applied' if reaches_asymptote else 'skipped (table short of asymptote)'}"
    )
    return model, TierResult("T1", "met" if met else "not met", e_rms, max(e_max, g_max), detail)
```

Add `fit_resonance` to `__all__`.

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run --no-sync pytest projects/potential_factory/test_fit.py -q -k resonance -m slow`
Expected: PASS (≈5 min)

- [ ] **Step 5: Commit**

```bash
uv run ruff check projects/potential_factory && uv run ruff format projects/potential_factory
git add projects/potential_factory/fit.py projects/potential_factory/test_fit.py
git commit -m "feat(factory): fit_resonance — track, smooth lam(R)/alpha(R), EA constraint, re-verify"
```

---

### Task 10: `fit_coupling` — the shell on `log Γ̃(ε, R)` (T3)

**Files:**
- Modify: `projects/potential_factory/fit.py`
- Test: `projects/potential_factory/test_fit.py` (append)

**Interfaces:**
- Consumes: `AsymptoticDiscreteState`, `v_dk_plus`, `gamma_from_coupling`.
- Produces: `model_gamma_tilde(model, pair, eps, R_asc) -> ndarray (eps.size, R.size)`; `fit_coupling(target: CouplingTarget, model, *, pair, n_eps=8, n_R=8, tol, r_b=3.0, alpha_b=2.0) -> tuple[FlexibleDiatomicModel, TierResult]` — installs a shell `SmoothR(f_inf=0, f_0=0, f_1=1, R_f=R_mid)` if absent, frees `shell.f_inf, shell.f_0, r_b` and least-squares `log(model Γ̃) − log(target Γ̃)` over the `(ε, R)` grid (`eps` geometric in `eps_window`, `R` uniform in `R_range`), then reports `rms/max` of the log residual; `status = "met"` iff `rms <= tol.coupling_log_rms`.

- [ ] **Step 1: Write the failing tests**

```python
# append to projects/potential_factory/test_fit.py
from projects.potential_factory.fit import fit_coupling, model_gamma_tilde  # noqa: E402


@pytest.mark.slow
def test_fit_coupling_round_trip_keeps_the_shell_negligible():
    pair = ElectronicPair()
    R_desc = np.linspace(3.0, 1.6, 8)
    target = extract_target(N2, pair=pair, R_desc=R_desc, n_eps=6)
    model, res = fit_coupling(target.coupling, from_diatomic(N2), pair=pair, n_eps=4, n_R=4, tol=Tolerances())
    assert res.status == "met", res.detail
    assert res.rms < 0.02
    assert abs(model.shell_R(2.0).real) < 1e-3


@pytest.mark.slow
def test_fit_coupling_moves_the_shell_for_a_steeper_falloff():
    pair = ElectronicPair()
    R_desc = np.linspace(3.0, 1.6, 8)
    base = extract_target(N2, pair=pair, R_desc=R_desc, n_eps=6).coupling
    eps = np.geomspace(*base.eps_window, 6)
    R_asc = R_desc[::-1]
    g = model_gamma_tilde(from_diatomic(N2), pair, eps, R_asc) * np.exp(-3.0 * eps)[:, None]
    from projects.potential_factory.target import CouplingTarget

    steeper = CouplingTarget.from_table(eps, R_asc, g, alpha=2.5)
    model, res = fit_coupling(steeper, from_diatomic(N2), pair=pair, n_eps=4, n_R=4, tol=Tolerances())
    assert model.shell is not None and abs(model.shell_R(2.0).real) > 1e-3
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run --no-sync pytest projects/potential_factory/test_fit.py -q -k coupling -m slow`
Expected: FAIL with `ImportError`

- [ ] **Step 3: Implement**

```python
# append to projects/potential_factory/fit.py
from qscat.core.nrm.coupling import gamma_from_coupling, v_dk_plus  # noqa: E402
from qscat.core.nrm.discrete_state import AsymptoticDiscreteState  # noqa: E402

from projects.potential_factory.target import CouplingTarget  # noqa: E402


def model_gamma_tilde(model: FlexibleDiatomicModel, pair: ElectronicPair, eps: np.ndarray, R_asc: np.ndarray) -> np.ndarray:
    """The model's own 2 pi |V_dk+(eps, R)|^2 with the R-independent discrete state."""
    phi_d = AsymptoticDiscreteState(pair.grid_a, model, _R_INF)
    out = np.empty((eps.size, R_asc.size))
    for i, e in enumerate(eps):
        out[i] = gamma_from_coupling(v_dk_plus(pair.grid_a, model, phi_d, R_asc, float(e)))
    return out


def fit_coupling(
    target: CouplingTarget,
    model: FlexibleDiatomicModel,
    *,
    pair: ElectronicPair,
    n_eps: int = 8,
    n_R: int = 8,
    tol: Tolerances,
    r_b: float = 3.0,
    alpha_b: float = 2.0,
) -> tuple[FlexibleDiatomicModel, TierResult]:
    eps = np.geomspace(target.eps_window[0], target.eps_window[1], n_eps)
    R_asc = np.linspace(target.R_range[0], target.R_range[1], n_R)
    y = np.log(target.gamma_tilde(eps[:, None], R_asc[None, :]))
    if model.shell is None:
        R_mid = 0.5 * (target.R_range[0] + target.R_range[1])
        model = model.with_shell(SmoothR(f_inf=0.0, f_0=0.0, f_1=1.0, R_f=R_mid, R_e=model.R_e), alpha_b, r_b)
    names = ["shell.f_inf", "shell.f_0", "r_b"]

    def resid(x):
        m = unpack(model, names, x)
        g = model_gamma_tilde(m, pair, eps, R_asc)
        return (np.log(np.maximum(g, 1e-300)) - y).ravel()

    x0 = pack(model, names)
    r0 = resid(x0)
    rms0 = float(np.sqrt(np.mean(r0**2)))
    if rms0 <= 0.5 * tol.coupling_log_rms:
        return model, TierResult("T3", "met", rms0, float(np.max(np.abs(r0))), "shell not needed; seed already within half the tolerance")
    sol = least_squares(resid, x0, diff_step=1e-3, xtol=1e-8, ftol=1e-8, max_nfev=60)
    fitted = unpack(model, names, sol.x)
    r = resid(sol.x)
    rms, mx = float(np.sqrt(np.mean(r**2))), float(np.max(np.abs(r)))
    status = "met" if rms <= tol.coupling_log_rms else "not met"
    return fitted, TierResult("T3", status, rms, mx, f"log Gamma~ over {n_eps}x{n_R} (eps, R); shell fitted")
```

Add both to `__all__`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run --no-sync pytest projects/potential_factory/test_fit.py -q -k coupling -m slow`
Expected: 2 passed (each `resid` costs `n_eps × n_R` dense solves; ≈5–10 min total)

- [ ] **Step 5: Commit**

```bash
uv run ruff check projects/potential_factory && uv run ruff format projects/potential_factory
git add projects/potential_factory/fit.py projects/potential_factory/test_fit.py
git commit -m "feat(factory): fit_coupling — shell fitted on log Gamma~(eps, R) against the target width"
```

---

### Task 11: `fit()` orchestration, ECS-boundedness check, and the end-to-end round trip

**Files:**
- Modify: `projects/potential_factory/fit.py`, `projects/potential_factory/report.py`
- Test: `projects/potential_factory/test_roundtrip.py`

**Interfaces:**
- Produces (`report.py`): `ecs_bounded(model, pair, R_tail: ndarray[complex]) -> dict[str, float]` — evaluates `|surface|` on the electronic tail points of `pair.grid_a` at `R = R_e` and on `R_tail` (nuclear tail points, given by the caller) at `r = 1.0`, and returns `{"electronic_deg": 45.0, "nuclear_deg": 90.0, "tail_growth": max|V_tail| / max|V_real|}`; raises `ValueError` if `tail_growth > 10`.
- Produces (`fit.py`): `fit(target: Target, seed: FlexibleDiatomicModel, *, pair: ElectronicPair, tol: Tolerances = Tolerances(), n_beta=1, n_nodes=24, lam_coeffs=0, alpha_coeffs=0) -> tuple[FlexibleDiatomicModel, FitReport]` running T0 → T1 → T3 in order over the tiers present, stopping after the first `"not met"` (later tiers `"not attempted"`), then filling `crossing_R` (the sign change of `v_ion − v0` from the target, `None` if none) and `da_threshold_sign` (`sign(−ea − (ε_0 − v0(R_e)))` using `qscat.core.vibrational.vibrational_states` on `qscat.core.grids.nuclear_grid()`; `+1` endothermic, `−1` exothermic).

- [ ] **Step 1: Write the failing test**

```python
# projects/potential_factory/test_roundtrip.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.model import F2, N2, NO

from projects.potential_factory.ansatz import from_diatomic, with_params
from projects.potential_factory.extract import extract_target
from projects.potential_factory.fit import fit
from projects.potential_factory.report import FitReport
from projects.potential_factory.tracker import ElectronicPair

CASES = {"N2": (N2, (3.0, 1.6)), "NO": (NO, (3.2, 1.7)), "F2": (F2, (4.0, 2.0))}


@pytest.mark.slow
@pytest.mark.parametrize("name", list(CASES))
def test_round_trip_recovers_the_published_model(name, tmp_path):
    model, (R_hi, R_lo) = CASES[name]
    pair = ElectronicPair()
    R_desc = np.linspace(R_hi, R_lo, 10)
    target = extract_target(model, pair=pair, R_desc=R_desc, n_eps=4, name=name)
    seed = with_params(
        from_diatomic(model),
        {"D_e": model.D0 * 1.2, "R_e": model.R0 * 1.05, "beta0": model.alpha0 * 0.9,
         "lam.f_inf": model.lambda_inf * 1.1, "alpha.f_inf": model.alpha_c * 1.3},
    )
    fitted, report = fit(target, seed, pair=pair, n_nodes=10)
    assert [t.status for t in report.tiers] == ["met", "met", "met"], [t.detail for t in report.tiers]
    assert abs(fitted.D_e - model.D0) < 1e-8 and abs(fitted.R_e - model.R0) < 1e-8
    np.testing.assert_allclose(fitted.lam_R(R_desc).real, model.lam(R_desc).real, rtol=3e-3)
    np.testing.assert_allclose(fitted.alpha_R(R_desc).real, model.alpha_c, rtol=3e-3)
    assert report.crossing_R is not None and abs(report.crossing_R - model.R_c) < 0.1
    report.to_json(tmp_path / "r.json")
    back = FitReport.from_json(tmp_path / "r.json")
    assert back.parameters == report.parameters and back.tiers == report.tiers
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run --no-sync pytest projects/potential_factory/test_roundtrip.py -q -m slow -k N2`
Expected: FAIL with `ImportError: cannot import name 'fit'`

- [ ] **Step 3: Implement**

```python
# append to projects/potential_factory/report.py
import numpy as np  # noqa: E402  (fold into imports)


def ecs_bounded(model, pair, R_tail) -> dict[str, float]:
    """Growth of |V| on the ECS tails relative to the real region; raises if > 10x."""
    r = pair.grid_a.points
    real = r.imag == 0.0
    v_r = np.abs(model.surface(r, model.R_e))
    v_R = np.abs(model.v0(np.asarray(R_tail, dtype=np.complex128)))
    ref = max(float(np.max(v_r[real][r[real].real > 0.3])), 1e-12)
    growth = max(float(np.max(v_r[~real])), float(np.max(v_R))) / ref
    if growth > 10.0:
        raise ValueError(f"potential grows {growth:.1f}x on the ECS tail; not absorbing")
    return {"electronic_deg": 45.0, "nuclear_deg": 90.0, "tail_growth": growth}
```

```python
# append to projects/potential_factory/fit.py
from qscat.core.grids import nuclear_grid  # noqa: E402
from qscat.core.vibrational import vibrational_states  # noqa: E402

from projects.potential_factory.ansatz import params  # noqa: E402
from projects.potential_factory.report import FitReport, ecs_bounded  # noqa: E402
from projects.potential_factory.target import Target  # noqa: E402


def _crossing(target: Target, model: FlexibleDiatomicModel) -> float | None:
    if target.resonance is None:
        return None
    R = np.linspace(target.resonance.R_range[0], target.resonance.R_range[1], 400)
    d = target.resonance.v_ion(R) - model.v0(R).real
    s = np.sign(d)
    idx = np.flatnonzero(s[:-1] * s[1:] < 0)
    if idx.size == 0:
        return None
    i = int(idx[0])
    return float(R[i] - d[i] * (R[i + 1] - R[i]) / (d[i + 1] - d[i]))


def _da_sign(target: Target, model: FlexibleDiatomicModel) -> int | None:
    if target.resonance is None:
        return None
    g = nuclear_grid()
    eps, _ = vibrational_states(g, model.mu, 1, model.v0)
    threshold = (-target.resonance.ea) - eps[0]  # DA threshold energy from v=0
    return int(np.sign(threshold)) if threshold != 0.0 else 0


def fit(
    target: Target,
    seed: FlexibleDiatomicModel,
    *,
    pair: ElectronicPair,
    tol: Tolerances = Tolerances(),
    n_beta: int = 1,
    n_nodes: int = 24,
    lam_coeffs: int = 0,
    alpha_coeffs: int = 0,
) -> tuple[FlexibleDiatomicModel, FitReport]:
    model = seed
    tiers: list[TierResult] = []
    halted = False
    stages = [
        ("T0", target.neutral, lambda m: fit_neutral(target.neutral, m, n_beta=n_beta, tol=tol)),
        ("T1", target.resonance, lambda m: fit_resonance(target.resonance, m, pair=pair, n_nodes=n_nodes, tol=tol, lam_coeffs=lam_coeffs, alpha_coeffs=alpha_coeffs)),
        ("T3", target.coupling, lambda m: fit_coupling(target.coupling, m, pair=pair, tol=tol)),
    ]
    for name, present, stage in stages:
        if present is None:
            tiers.append(TierResult(name, "not attempted", np.nan, np.nan, "no target data"))
            continue
        if halted:
            tiers.append(TierResult(name, "not attempted", np.nan, np.nan, "an earlier tier was not met"))
            continue
        model, res = stage(model)
        tiers.append(res)
        halted = res.status == "not met"
    R_tail = 12.0 + np.linspace(0.1, 6.0, 8) * np.exp(1j * np.deg2rad(35.0))
    report = FitReport(
        target_name=target.name,
        parameters=params(model),
        tiers=tiers,
        ecs_bounds_deg=ecs_bounded(model, pair, R_tail),
        crossing_R=_crossing(target, model),
        da_threshold_sign=_da_sign(target, model),
        provenance={k: {"source": v.source, "locator": v.locator} for k, v in target.provenance.items()},
    )
    return model, report
```

Add `fit` to `__all__`; `ecs_bounded` to `report.__all__`.

- [ ] **Step 4: Run the round trips**

Run: `uv run --no-sync pytest projects/potential_factory/test_roundtrip.py -q -m slow`
Expected: 3 passed (≈15–30 min total; run `-k N2` first while iterating)

- [ ] **Step 5: Run the whole package + the core import guard, then commit**

Run: `uv run --no-sync pytest projects/potential_factory libs/qscat/tests/test_core_no_model_import.py -q -m "not slow"`
Expected: all pass

```bash
uv run ruff check projects/potential_factory && uv run ruff format projects/potential_factory
git add projects/potential_factory/fit.py projects/potential_factory/report.py projects/potential_factory/test_roundtrip.py
git commit -m "feat(factory): fit() orchestration + FitReport with crossing/DA-sign/ECS checks; N2/NO/F2 round-trip oracle"
```

---

### Task 12: Theory note and repo map entry

**Files:**
- Create: `docs/physics/potential-factory.md`
- Modify: `docs/physics/open-directions.md` (the survey stays listed; add a line pointing at the new note) and `docs/physics/README.md` is untouched (the note joins the `resonances` section toctree — open `docs/physics/resonances.md` and add `potential-factory` to its toctree)
- Modify: `CLAUDE.md` — add a `projects/potential_factory` bullet under `projects/`

- [ ] **Step 1: Write the note**

The note follows the house form (Status, Relates to, Units; then the method, the equations actually used, the validation evidence with the measured round-trip numbers from Task 11's run, limitations). Required content:

```markdown
# Potential factory — fitting model surfaces to target curves

**Status:** toy stage (`projects/potential_factory/`); promoted only after the O₂ plan.
**Relates to:** `docs/physics/potential-factory-options.md` (why), the spec
`docs/superpowers/specs/2026-08-24-potential-factory-design.md` (what),
`docs/physics/nonlocal-resonance-model.md` (the T3 forward model), `docs/physics/n2-resonance.md`.
**Units:** atomic units.

## What it does
(Target tiers T0/T1/T3; ansatz `FlexibleDiatomicModel`, the EMO and the two free functions; the shell.)

## The forward models and the gradient
(`ElectronicPair.pole` and the gate; `pole_sensitivity = psi^2` under the c-product — state that this
is Hellmann–Feynman for a complex-symmetric matrix and give the 1e-3 finite-difference check from
Task 4; Newton in (λ, α); continuation.)

## Round-trip oracle — measured
(The table from `test_roundtrip.py`: per molecule, T0 constants recovered to 1e-8; λ(R), α(R) curves to
the measured rtol; E_res rms and Γ max-rel from the T1 `TierResult.detail`; T3 log-rms; wall time.)

## Limitations
(Sigmoid-constant degeneracy — why the gate is on curves; `Tolerances` are placeholders until the
budget plan; the fake near-threshold pole and the gate; `l = 0` needs the shell; polar targets out of scope.)
```

Fill every parenthesis with the actual measured numbers from the Task 11 run (copy `TierResult.detail` strings).

- [ ] **Step 2: Add the CLAUDE.md bullet**

Under `projects/`, after the `n2_2d_td_cross_section` entry:

```
            - `potential_factory`: the toy-stage POTENTIAL FACTORY — fits a
              `FlexibleDiatomicModel` (EMO `v0` + Gaussian well with `lam(R)`
              AND `alpha(R)` + optional shell; embeds N2/NO/F2 exactly) to a
              tiered `Target` (T0 neutral curve, T1 pole curves, T3 the
              published energy-dependent width) in stages that stop-and-report.
              Proven by round-tripping the existing models' OWN calculated
              curves (`extract_target`) back to their constants/curves. Nothing
              is fitted to experiment — see docs/physics/potential-factory.md.
```

- [ ] **Step 3: Build the docs locally if the plot/docs extras are installed (optional), then commit**

```bash
git add docs/physics/potential-factory.md docs/physics/resonances.md docs/physics/open-directions.md CLAUDE.md
git commit -m "docs(factory): theory note with the measured N2/NO/F2 round trip; repo-map entry"
```

---

## Self-review against the spec

- **Target format** (tiers, coordinates tuple, provenance, reserved `eigenphase`) → Task 7. The `rydberg` slot for ions is *not* added (spec: out of scope for v1; the `Target` dataclass gains it in the O₂ plan's promotion task if needed — recorded there).
- **Ansatz** (EMO, `λ(R)`, `α(R)`, shell, ECS analyticity, `ell` not fitted, Houfek models as exact points) → Tasks 1–2; the boundedness check → Task 11.
- **Forward models** (pole walk, EA, `Γ̃` via `v_dk_plus`, c-product gradient) → Tasks 3–7, 10.
- **Staged fitting, stop-and-report, EA hard constraint, re-verify after smoothing, `log Γ̃` loss, `FitReport` fields (parameters, per-tier residuals, verdicts, ECS bounds, `R_c`, DA sign, provenance)** → Tasks 8–11.
- **Validation (a) round trip on N₂/NO/F₂, gate on `v0` constants + curves** → Task 11. **(e) spurious-pole guard** → Task 3. **(b), (c), (d)** belong to the O₂ plan; the **tolerance budget** to the budget plan — `Tolerances` is explicitly a placeholder (Task 8).
- **`MoleculePreset` / serialised parameter file** → JSON in Task 11; the preset is the O₂ plan's job (it needs tuner grids for a *new* molecule; the round-trip molecules already have decks).
- Type consistency checked: `SmoothR(f_inf, f_0, f_1, R_f, coeffs, R_e, p)`, `FlexibleDiatomicModel(..., lam, alpha, shell, alpha_b, r_b)`, `Pole.energy/residual/gamma/shift`, `WellParams(lam, alpha)`, `TrackResult(R, lam, alpha, converged, poles)`, `TierResult(name, status, rms, max, detail)`, `Tolerances(v0_rms, omega_e_rel, e_res_rms, gamma_rel, gamma_floor, coupling_log_rms)`, `FitReport(target_name, parameters, tiers, ecs_bounds_deg, crossing_R, da_threshold_sign, provenance)` are used with these exact names throughout.
