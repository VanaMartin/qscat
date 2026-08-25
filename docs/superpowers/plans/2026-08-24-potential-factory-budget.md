# Potential Factory — Tolerance Budget (N₂/NO sensitivity study) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the placeholder `Tolerances` of the factory with a per-feature budget measured on N₂-like and NO-like: the largest deviation of each curve feature (`E_res`, `Γ`, `R_c`, `ω_e`, `B`) that keeps every exact observable within a stated observable-level tolerance.

**Architecture:** `validation/factory/` holds a perturbation catalogue (parameter-space moves realised on `FlexibleDiatomicModel`, so every perturbed model is a legitimate `ResonanceModel`), a curve-metric function (measures what a perturbation *did* to the curves — measured, not assumed), an observable-metric function (exact VE/DA cross sections, BO levels, scored with `qscat.core.assignment.peak_positions`/`peak_alignment`), and a bisection that finds the budget per feature. A CLI writes a CSV + Markdown table; a theory note records it; the factory's `Tolerances` defaults are set from the table and a test pins them to it.

**Tech Stack:** numpy, scipy, pytest (`@pytest.mark.slow`), `qscat.core.{driven,dissociation,lcp,vibrational,assignment,grids}`, `qscat_run.presets` (the per-molecule decks; `validation` may import `apps`, as `validation/diatomic/test_da_grid.py` already does), `projects.potential_factory`.

**Spec:** `docs/superpowers/specs/2026-08-24-potential-factory-design.md` — section "Tolerance budget — calibrated on N₂ and NO". Depends on the core plan `2026-08-24-potential-factory-core.md` being merged.

## Global Constraints

- The study is **oracle-based**: every metric compares the perturbed model's exact 2-D solution with the unperturbed model's exact 2-D solution. **No experimental data enters.**
- Observable-level tolerance (the one hand-chosen number set): peak positions within **0.25** of a resonance width; peak heights and integrated σ within **10 %**; BO level positions within **0.1** of their width; level widths within **10 %**. Constants in `validation/factory/budget.py` as `OBSERVABLE_TOL`.
- NO's DA enters through peak positions and the threshold only (never relative magnitude — exponentially sensitive near threshold, `docs/physics/nonlocal-resonance-model.md` §8).
- Converged decks only: `qscat_run.presets` `ti_grid()` for N₂ and NO (the eMoScat-derived, gated grids), default energies from the preset.
- Atomic units; eV only in tables' headers where stated.
- Runs are long (each cross-section sweep is minutes); every test that solves the 2-D problem is `@pytest.mark.slow`, and the CLI caches results as `.npz` under `validation/factory/results/` (gitignored except the final CSV/MD).

---

### Task 1: Perturbation catalogue and curve metrics

**Files:**
- Create: `validation/factory/__init__.py` (empty), `validation/factory/perturb.py`
- Test: `validation/factory/test_perturb.py`

**Interfaces:**
- Consumes: `projects.potential_factory.ansatz.{from_diatomic, with_params, params}`, `projects.potential_factory.tracker.ElectronicPair`, `projects.potential_factory.extract.extract_target`, `projects.potential_factory.fit.model_gamma_tilde`, `qscat.core.vibrational.vibrational_states`, `qscat.core.grids.nuclear_grid`.
- Produces:
  - `Perturbation(name: str, param: str, feature: str)` and `CATALOGUE: tuple[Perturbation, ...]` = `(("e_res_shift", "lam.f_inf", "E_res"), ("gamma_scale", "alpha.f_inf", "Gamma"), ("crossing_move", "lam.R_f", "R_c"), ("ladder", "beta0", "omega_e"), ("falloff", "shell.f_inf", "B"))`
  - `perturbed(model, pert, amount) -> FlexibleDiatomicModel` — multiplicative on the parameter (`value * (1 + amount)`), except `crossing_move` (additive, bohr) and `falloff` (additive, Hartree; installs a shell at `r_b = 3.0, alpha_b = 2.0` if absent)
  - `CurveDeltas(e_res_rms: float, gamma_rel_max: float, r_c_shift: float, omega_e_rel: float, b_rel: float)` and `curve_deltas(base, pert_model, *, pair, R_desc) -> CurveDeltas` — extracts both models' T1/T3 targets with `extract_target` and reports: rms of `ΔE_res` over nodes (Ha); max relative `ΔΓ` where `Γ_base > 2e-3`; shift of the crossing (bohr); relative change of `ε_1 − ε_0` from `vibrational_states(nuclear_grid(), mu, 2, v0)`; relative change of the log-slope of `Γ̃` in `ε` at `R_e` (`B ≈ −d ln Γ̃/dε` after removing `ε^{l+1/2}`).

- [ ] **Step 1: Write the failing tests**

```python
# validation/factory/test_perturb.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.model import N2

from projects.potential_factory.ansatz import from_diatomic
from projects.potential_factory.tracker import ElectronicPair
from validation.factory.perturb import CATALOGUE, curve_deltas, perturbed


def test_catalogue_covers_the_five_spec_features():
    assert {p.feature for p in CATALOGUE} == {"E_res", "Gamma", "R_c", "omega_e", "B"}


def test_perturbed_moves_only_the_named_parameter():
    base = from_diatomic(N2)
    p = next(p for p in CATALOGUE if p.name == "ladder")
    m = perturbed(base, p, 0.05)
    assert abs(m.betas[0] - 1.05 * base.betas[0]) < 1e-15 and m.D_e == base.D_e


def test_falloff_installs_a_shell():
    base = from_diatomic(N2)
    p = next(p for p in CATALOGUE if p.name == "falloff")
    m = perturbed(base, p, 0.01)
    assert m.shell is not None and abs(m.shell_R(2.0).real - 0.01) < 1e-15


@pytest.mark.slow
def test_curve_deltas_zero_perturbation_is_zero_and_ladder_moves_omega_e():
    base = from_diatomic(N2)
    pair = ElectronicPair()
    R_desc = np.linspace(3.0, 1.6, 8)
    z = curve_deltas(base, base, pair=pair, R_desc=R_desc)
    assert z.e_res_rms == 0.0 and z.gamma_rel_max == 0.0 and z.omega_e_rel == 0.0
    p = next(p for p in CATALOGUE if p.name == "ladder")
    d = curve_deltas(base, perturbed(base, p, 0.05), pair=pair, R_desc=R_desc)
    assert 0.03 < d.omega_e_rel < 0.07
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run --no-sync pytest validation/factory/test_perturb.py -q -m "not slow"`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement**

```python
# validation/factory/perturb.py
"""Parameter-space perturbations of a FlexibleDiatomicModel and the MEASURED
curve deviations they cause (nothing is assumed about what a parameter does)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.core.grids import nuclear_grid
from qscat.core.vibrational import vibrational_states

from projects.potential_factory.ansatz import FlexibleDiatomicModel, SmoothR, with_params
from projects.potential_factory.extract import extract_target
from projects.potential_factory.fit import model_gamma_tilde
from projects.potential_factory.tracker import ElectronicPair

__all__ = ["Perturbation", "CATALOGUE", "perturbed", "CurveDeltas", "curve_deltas"]


@dataclass(frozen=True)
class Perturbation:
    name: str
    param: str
    feature: str


CATALOGUE: tuple[Perturbation, ...] = (
    Perturbation("e_res_shift", "lam.f_inf", "E_res"),
    Perturbation("gamma_scale", "alpha.f_inf", "Gamma"),
    Perturbation("crossing_move", "lam.R_f", "R_c"),
    Perturbation("ladder", "beta0", "omega_e"),
    Perturbation("falloff", "shell.f_inf", "B"),
)

_ADDITIVE = {"crossing_move", "falloff"}


def perturbed(model: FlexibleDiatomicModel, pert: Perturbation, amount: float) -> FlexibleDiatomicModel:
    if pert.name == "falloff" and model.shell is None:
        model = model.with_shell(SmoothR(f_inf=0.0, f_0=0.0, f_1=1.0, R_f=model.R_e, R_e=model.R_e), 2.0, 3.0)
    from projects.potential_factory.ansatz import params

    cur = params(model)[pert.param]
    new = cur + amount if pert.name in _ADDITIVE else cur * (1.0 + amount)
    return with_params(model, {pert.param: new})


@dataclass(frozen=True)
class CurveDeltas:
    e_res_rms: float
    gamma_rel_max: float
    r_c_shift: float
    omega_e_rel: float
    b_rel: float


def _crossing(v_ion: npt.NDArray[np.float64], v0: npt.NDArray[np.float64], R: npt.NDArray[np.float64]) -> float:
    d = v_ion - v0
    s = np.sign(d)
    idx = np.flatnonzero(s[:-1] * s[1:] < 0)
    if idx.size == 0:
        return float("nan")
    i = int(idx[0])
    return float(R[i] - d[i] * (R[i + 1] - R[i]) / (d[i + 1] - d[i]))


def _omega_e(model: FlexibleDiatomicModel) -> float:
    eps, _ = vibrational_states(nuclear_grid(), model.mu, 2, model.v0)
    return float(eps[1] - eps[0])


def _b_slope(model: FlexibleDiatomicModel, pair: ElectronicPair) -> float:
    eps = np.geomspace(0.01, 0.2, 6)
    g = model_gamma_tilde(model, pair, eps, np.array([model.R_e]))[:, 0]
    y = np.log(g) - (model.ell + 0.5) * np.log(eps)
    slope, _ = np.polyfit(eps, y, 1)
    return float(-slope)


def curve_deltas(
    base: FlexibleDiatomicModel,
    pert_model: FlexibleDiatomicModel,
    *,
    pair: ElectronicPair,
    R_desc: npt.ArrayLike,
) -> CurveDeltas:
    R = np.asarray(R_desc, dtype=np.float64)
    tb = extract_target(base, pair=pair, R_desc=R, n_eps=3)
    tp = extract_target(pert_model, pair=pair, R_desc=R, n_eps=3)
    assert tb.resonance is not None and tp.resonance is not None
    e_b, e_p = tb.resonance.v_ion(R), tp.resonance.v_ion(R)
    g_b, g_p = tb.resonance.gamma(R), tp.resonance.gamma(R)
    mask = g_b > 2e-3
    g_rel = np.max(np.abs(g_p[mask] - g_b[mask]) / g_b[mask]) if mask.any() else 0.0
    rc_b = _crossing(e_b, base.v0(R).real, R)
    rc_p = _crossing(e_p, pert_model.v0(R).real, R)
    w_b, w_p = _omega_e(base), _omega_e(pert_model)
    b_b, b_p = _b_slope(base, pair), _b_slope(pert_model, pair)
    return CurveDeltas(
        e_res_rms=float(np.sqrt(np.mean((e_p - e_b) ** 2))),
        gamma_rel_max=float(g_rel),
        r_c_shift=float(rc_p - rc_b),
        omega_e_rel=float(abs(w_p - w_b) / w_b),
        b_rel=float(abs(b_p - b_b) / max(abs(b_b), 1e-12)),
    )
```

(Move the inner `params` import to the top of the module.)

- [ ] **Step 4: Run tests**

Run: `uv run --no-sync pytest validation/factory/test_perturb.py -q -m "not slow"` then `-m slow`
Expected: 3 passed; then 1 passed

- [ ] **Step 5: Commit**

```bash
uv run ruff check validation/factory && uv run ruff format validation/factory
git add validation/factory/__init__.py validation/factory/perturb.py validation/factory/test_perturb.py
git commit -m "feat(validation/factory): perturbation catalogue + measured curve deltas"
```

---

### Task 2: Observable metrics on the converged decks

**Files:**
- Create: `validation/factory/observables.py`
- Test: `validation/factory/test_observables.py`

**Interfaces:**
- Consumes: `qscat_run.presets` (`PRESETS["N2"]`/`["NO"]` — use whatever lookup `presets.py` exposes; read it first: the preset has `ti_grid()`, `default_energies`, `n_vib`), `qscat.core.vibrational.vibrational_states`, `qscat.core.driven.ve_cross_section`, `qscat.core.dissociation.{da_cross_section, anion_electronic_states}`, `qscat.core.lcp.resonance_levels`, `qscat.core.assignment.{peak_positions, peak_alignment}`, `qscat.core.grids.{electronic_grid, nuclear_grid}`.
- Produces:
  - `Observables(E: ndarray, sigma_ve: dict[int, ndarray], sigma_da: ndarray | None, level_E: ndarray, level_gamma: ndarray)` and `compute_observables(model, molecule: str, *, vprimes=(0, 1, 2), n_energies: int | None = None) -> Observables` — the TI VE sweep on the preset's `ti_grid()` over the preset's default energies (subsampled to `n_energies` when given), DA only when the preset lists `"da"` in `valid_observables`, BO levels via `resonance_levels` on two nuclear angles (`nuclear_grid(angle_deg=35)`, `nuclear_grid(angle_deg=40)`) and the electronic pair `electronic_grid(angle_deg=35/44)`, `n_levels=6`.
  - `ObservableDeltas(peak_pos_widths: float, peak_height_rel: float, integrated_rel: float, level_pos_widths: float, level_width_rel: float, da_threshold_shift: float | None)` and `observable_deltas(base: Observables, pert: Observables) -> ObservableDeltas` — peak positions via `peak_positions(E, σ_v')` for each `v'`, aligned with `peak_alignment(marks=pert_peaks, peaks=base_peaks, width=mean(base.level_gamma))` → max distance in widths; peak heights: max relative change of `σ` at the base peak indices; integrated: `|∫σ_pert − ∫σ_base| / ∫σ_base` (trapezoid) max over `v'`; levels: max `|ΔE_v| / Γ_v` and max `|ΔΓ_v| / Γ_v`; DA threshold shift: energy of the first `σ_DA > 1e-3 · max` (None when no DA).

- [ ] **Step 1: Write the failing tests**

```python
# validation/factory/test_observables.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.model import N2

from projects.potential_factory.ansatz import from_diatomic
from validation.factory.observables import Observables, compute_observables, observable_deltas


def test_observable_deltas_of_identical_inputs_is_zero():
    E = np.linspace(0.05, 0.15, 60)
    s = 1.0 + 50.0 * np.exp(-((E - 0.1) / 0.003) ** 2) + 30.0 * np.exp(-((E - 0.12) / 0.003) ** 2)
    obs = Observables(E=E, sigma_ve={0: s, 1: 0.5 * s}, sigma_da=None, level_E=np.array([0.1, 0.12]), level_gamma=np.array([0.004, 0.004]))
    d = observable_deltas(obs, obs)
    assert d.peak_pos_widths == 0.0 and d.peak_height_rel == 0.0 and d.integrated_rel == 0.0 and d.level_pos_widths == 0.0


def test_observable_deltas_sees_a_shifted_peak_in_width_units():
    E = np.linspace(0.05, 0.15, 400)
    peak = lambda c: 1.0 + 50.0 * np.exp(-((E - c) / 0.003) ** 2)  # noqa: E731
    a = Observables(E=E, sigma_ve={0: peak(0.100)}, sigma_da=None, level_E=np.array([0.1]), level_gamma=np.array([0.004]))
    b = Observables(E=E, sigma_ve={0: peak(0.102)}, sigma_da=None, level_E=np.array([0.102]), level_gamma=np.array([0.004]))
    d = observable_deltas(a, b)
    assert 0.4 < d.peak_pos_widths < 0.6 and 0.4 < d.level_pos_widths < 0.6


@pytest.mark.slow
def test_compute_observables_runs_on_the_n2_deck():
    obs = compute_observables(from_diatomic(N2), "N2", n_energies=12)
    assert set(obs.sigma_ve) == {0, 1, 2} and obs.sigma_da is None
    assert obs.level_E.size == 6 and np.all(obs.level_gamma > 0)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run --no-sync pytest validation/factory/test_observables.py -q -m "not slow"`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement**

```python
# validation/factory/observables.py
"""Exact observables of a model on its converged deck, and how much they moved."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.core.assignment import peak_alignment, peak_positions
from qscat.core.dissociation import anion_electronic_states, da_cross_section
from qscat.core.driven import ve_cross_section
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.lcp import resonance_levels
from qscat.core.vibrational import vibrational_states
from qscat.model import ResonanceModel
from qscat_run import presets  # read presets.py for the exact lookup name

__all__ = ["Observables", "compute_observables", "ObservableDeltas", "observable_deltas"]

FArr = npt.NDArray[np.float64]


@dataclass(frozen=True)
class Observables:
    E: FArr
    sigma_ve: dict[int, FArr]
    sigma_da: FArr | None
    level_E: FArr
    level_gamma: FArr


def compute_observables(
    model: ResonanceModel,
    molecule: str,
    *,
    vprimes: tuple[int, ...] = (0, 1, 2),
    n_energies: int | None = None,
) -> Observables:
    preset = presets.get(molecule)  # replace with the real accessor after reading presets.py
    tgrid = preset.ti_grid()
    E = np.asarray(preset.default_energies.values(), dtype=np.float64)  # ditto: EnergySpec -> array
    if n_energies is not None and n_energies < E.size:
        E = E[np.linspace(0, E.size - 1, n_energies).astype(int)]
    g_R = tgrid.grids[1]
    eps, chi = vibrational_states(g_R, model.mu, preset.n_vib, model.v0)
    sig = ve_cross_section(tgrid, model, eps, chi, 0, list(vprimes), E)
    sigma_ve = {v: np.asarray(sig[:, i], dtype=np.float64) for i, v in enumerate(vprimes)}
    sigma_da = None
    if "da" in preset.valid_observables:
        sigma_da = np.asarray(da_cross_section(tgrid, model, eps, chi, 0, E), dtype=np.float64)
    ga, gb = electronic_grid(angle_deg=35.0), electronic_grid(angle_deg=44.0)
    lv = resonance_levels(model, nuclear_grid(angle_deg=35.0), nuclear_grid(angle_deg=40.0), ga, gb, n_levels=6)
    return Observables(E=E, sigma_ve=sigma_ve, sigma_da=sigma_da, level_E=np.asarray(lv.energies.real), level_gamma=np.asarray(-2.0 * lv.energies.imag))
```

(`ResonanceLevels`' attribute for the complex levels must be read from `qscat/core/lcp.py:502` and used by its real name; the `sig` indexing must follow `ve_cross_section`'s documented return shape — check `driven.py`'s docstring: `(E.size, len(vprimes))` or the transpose — and fix the line accordingly.)

```python
@dataclass(frozen=True)
class ObservableDeltas:
    peak_pos_widths: float
    peak_height_rel: float
    integrated_rel: float
    level_pos_widths: float
    level_width_rel: float
    da_threshold_shift: float | None


def _da_threshold(E: FArr, s: FArr) -> float:
    idx = np.flatnonzero(s > 1e-3 * np.max(s))
    return float(E[idx[0]]) if idx.size else float("nan")


def observable_deltas(base: Observables, pert: Observables) -> ObservableDeltas:
    width = float(np.mean(base.level_gamma))
    pos, height, integ = 0.0, 0.0, 0.0
    for v, s_b in base.sigma_ve.items():
        s_p = pert.sigma_ve[v]
        pk_b = peak_positions(base.E, s_b)
        pk_p = peak_positions(pert.E, s_p)
        if pk_b.size and pk_p.size:
            al = peak_alignment(pk_p, pk_b, width=width)
            pos = max(pos, float(np.max(al.distances)))  # use PeakAlignment's real field name
            ib = np.searchsorted(base.E, pk_b)
            height = max(height, float(np.max(np.abs(s_p[ib] - s_b[ib]) / s_b[ib])))
        integ = max(integ, float(abs(np.trapz(s_p, pert.E) - np.trapz(s_b, base.E)) / np.trapz(s_b, base.E)))
    n = min(base.level_E.size, pert.level_E.size)
    lp = float(np.max(np.abs(pert.level_E[:n] - base.level_E[:n]) / base.level_gamma[:n]))
    lw = float(np.max(np.abs(pert.level_gamma[:n] - base.level_gamma[:n]) / base.level_gamma[:n]))
    da = None
    if base.sigma_da is not None and pert.sigma_da is not None:
        da = _da_threshold(pert.E, pert.sigma_da) - _da_threshold(base.E, base.sigma_da)
    return ObservableDeltas(pos, height, integ, lp, lw, da)
```

- [ ] **Step 4: Run tests**

Run: `uv run --no-sync pytest validation/factory/test_observables.py -q -m "not slow"` then `-m slow`
Expected: 2 passed; then 1 passed (minutes)

- [ ] **Step 5: Commit**

```bash
uv run ruff check validation/factory && uv run ruff format validation/factory
git add validation/factory/observables.py validation/factory/test_observables.py
git commit -m "feat(validation/factory): exact observables on the preset decks + width-unit deltas"
```

---

### Task 3: The budget search and CLI

**Files:**
- Create: `validation/factory/budget.py`, `validation/factory/sensitivity.py` (the `__main__`)
- Modify: `.gitignore` (add `validation/factory/results/*.npz`)
- Test: `validation/factory/test_budget.py`

**Interfaces:**
- Produces (`budget.py`): `OBSERVABLE_TOL = ObservableDeltas(peak_pos_widths=0.25, peak_height_rel=0.10, integrated_rel=0.10, level_pos_widths=0.10, level_width_rel=0.10, da_threshold_shift=None)`; `within(d: ObservableDeltas, tol=OBSERVABLE_TOL) -> bool` (DA threshold compared as `|shift| <= 0.25 * width` only when both are present — pass `width` through `within(d, tol, width)`); `BudgetRow(molecule, perturbation, feature, amount, curve: CurveDeltas, obs: ObservableDeltas)`; `find_budget(model, molecule, pert, *, pair, R_desc, amounts: Sequence[float], n_energies) -> list[BudgetRow]` — evaluates the catalogue at the given amounts (a geometric ladder, e.g. `0.005, 0.01, 0.02, 0.04, 0.08`), stopping at the first amount that is *not* `within`; the budget for that feature is the last `within` row's `curve` delta.
- Produces (`sensitivity.py`): `python -m validation.factory.sensitivity --molecule N2 [--n-energies 40] [--out validation/factory/results]` → `results/<mol>_rows.csv` (one `BudgetRow` per line) and `results/<mol>_budget.md` (the per-feature budget table), caching each `Observables` as `results/<mol>_<pert>_<amount>.npz`.

- [ ] **Step 1: Write the failing test**

```python
# validation/factory/test_budget.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.model import N2

from projects.potential_factory.ansatz import from_diatomic
from projects.potential_factory.tracker import ElectronicPair
from validation.factory.budget import OBSERVABLE_TOL, find_budget, within
from validation.factory.observables import ObservableDeltas
from validation.factory.perturb import CATALOGUE


def test_within_respects_every_field():
    ok = ObservableDeltas(0.1, 0.05, 0.05, 0.05, 0.05, None)
    bad = ObservableDeltas(0.3, 0.05, 0.05, 0.05, 0.05, None)
    assert within(ok, OBSERVABLE_TOL, width=0.004) and not within(bad, OBSERVABLE_TOL, width=0.004)


@pytest.mark.slow
def test_find_budget_ladder_stops_at_the_first_violation():
    pert = next(p for p in CATALOGUE if p.name == "ladder")
    rows = find_budget(from_diatomic(N2), "N2", pert, pair=ElectronicPair(), R_desc=np.linspace(3.0, 1.6, 8), amounts=(0.002, 0.2), n_energies=10)
    assert 1 <= len(rows) <= 2
    assert rows[0].amount == 0.002
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run --no-sync pytest validation/factory/test_budget.py -q -m "not slow"`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement**

```python
# validation/factory/budget.py
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from projects.potential_factory.ansatz import FlexibleDiatomicModel
from projects.potential_factory.tracker import ElectronicPair
from validation.factory.observables import ObservableDeltas, compute_observables, observable_deltas
from validation.factory.perturb import CurveDeltas, Perturbation, curve_deltas, perturbed

__all__ = ["OBSERVABLE_TOL", "within", "BudgetRow", "find_budget"]

OBSERVABLE_TOL = ObservableDeltas(
    peak_pos_widths=0.25, peak_height_rel=0.10, integrated_rel=0.10,
    level_pos_widths=0.10, level_width_rel=0.10, da_threshold_shift=None,
)


def within(d: ObservableDeltas, tol: ObservableDeltas = OBSERVABLE_TOL, *, width: float) -> bool:
    ok = (
        d.peak_pos_widths <= tol.peak_pos_widths
        and d.peak_height_rel <= tol.peak_height_rel
        and d.integrated_rel <= tol.integrated_rel
        and d.level_pos_widths <= tol.level_pos_widths
        and d.level_width_rel <= tol.level_width_rel
    )
    if d.da_threshold_shift is not None:
        ok = ok and abs(d.da_threshold_shift) <= 0.25 * width
    return ok


@dataclass(frozen=True)
class BudgetRow:
    molecule: str
    perturbation: str
    feature: str
    amount: float
    curve: CurveDeltas
    obs: ObservableDeltas


def find_budget(
    model: FlexibleDiatomicModel,
    molecule: str,
    pert: Perturbation,
    *,
    pair: ElectronicPair,
    R_desc: npt.ArrayLike,
    amounts: Sequence[float],
    n_energies: int | None = None,
) -> list[BudgetRow]:
    base_obs = compute_observables(model, molecule, n_energies=n_energies)
    width = float(np.mean(base_obs.level_gamma))
    rows: list[BudgetRow] = []
    for a in amounts:
        m = perturbed(model, pert, a)
        obs = observable_deltas(base_obs, compute_observables(m, molecule, n_energies=n_energies))
        rows.append(BudgetRow(molecule, pert.name, pert.feature, float(a), curve_deltas(model, m, pair=pair, R_desc=R_desc), obs))
        if not within(obs, width=width):
            break
    return rows
```

```python
# validation/factory/sensitivity.py
"""python -m validation.factory.sensitivity --molecule N2 -- the tolerance-budget study."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
from pathlib import Path

import numpy as np
from qscat.model import N2, NO

from projects.potential_factory.ansatz import from_diatomic
from projects.potential_factory.tracker import ElectronicPair
from validation.factory.budget import find_budget, within
from validation.factory.perturb import CATALOGUE

AMOUNTS = (0.005, 0.01, 0.02, 0.04, 0.08, 0.16)
MODELS = {"N2": (N2, (3.0, 1.6)), "NO": (NO, (3.2, 1.7))}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--molecule", choices=list(MODELS), required=True)
    ap.add_argument("--n-energies", type=int, default=None)
    ap.add_argument("--out", type=Path, default=Path("validation/factory/results"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    model, (R_hi, R_lo) = MODELS[args.molecule]
    flex = from_diatomic(model)
    pair = ElectronicPair()
    R_desc = np.linspace(R_hi, R_lo, 12)
    rows = []
    for pert in CATALOGUE:
        rows += find_budget(flex, args.molecule, pert, pair=pair, R_desc=R_desc, amounts=AMOUNTS, n_energies=args.n_energies)
    with (args.out / f"{args.molecule}_rows.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["perturbation", "feature", "amount", *asdict(rows[0].curve), *asdict(rows[0].obs)])
        for r in rows:
            w.writerow([r.perturbation, r.feature, r.amount, *asdict(r.curve).values(), *asdict(r.obs).values()])
    lines = ["| feature | perturbation | largest amount within tolerance | curve deviation at that amount |", "|---|---|---|---|"]
    for pert in CATALOGUE:
        ok = [r for r in rows if r.perturbation == pert.name and within(r.obs, width=1.0)]
        if ok:
            r = ok[-1]
            lines.append(f"| {r.feature} | {r.perturbation} | {r.amount:g} | {asdict(r.curve)} |")
        else:
            lines.append(f"| {pert.feature} | {pert.name} | none (smallest amount already violates) | - |")
    (args.out / f"{args.molecule}_budget.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
```

(`within(..., width=1.0)` in the summary is only a placeholder for the printout when `da_threshold_shift` is None; when NO's DA rows are present, thread the base width through by storing it on `BudgetRow` — add a `width: float` field and use it.)

- [ ] **Step 4: Run tests, then the study**

Run: `uv run --no-sync pytest validation/factory/test_budget.py -q -m "not slow"` → 1 passed; `-m slow` → 1 passed.
Run: `uv run --no-sync python -m validation.factory.sensitivity --molecule N2` then `--molecule NO` (each: hours on the laptop — run in Docker/`sadaharu` if available; the `.npz` cache makes reruns cheap).
Expected: `validation/factory/results/{N2,NO}_budget.md` written.

- [ ] **Step 5: Commit (code + the two CSV/MD result files; `.npz` ignored)**

```bash
printf 'validation/factory/results/*.npz\n' >> .gitignore
uv run ruff check validation/factory && uv run ruff format validation/factory
git add .gitignore validation/factory/budget.py validation/factory/sensitivity.py validation/factory/test_budget.py validation/factory/results/N2_rows.csv validation/factory/results/N2_budget.md validation/factory/results/NO_rows.csv validation/factory/results/NO_budget.md
git commit -m "feat(validation/factory): tolerance-budget search + CLI; N2/NO budget tables"
```

---

### Task 4: Set the factory's `Tolerances` from the budget and record it

**Files:**
- Modify: `projects/potential_factory/report.py` (`Tolerances` defaults)
- Create: `docs/physics/potential-factory-budget.md`
- Modify: `docs/physics/potential-factory.md` (replace the "placeholder" sentence with a pointer), `docs/physics/resonances.md` toctree
- Test: `projects/potential_factory/test_report.py`

**Interfaces:**
- Produces: `Tolerances` defaults = the **minimum over N₂ and NO** of the per-feature budget (`e_res_rms` from the `E_res` row's `curve.e_res_rms`, `gamma_rel` from the `Gamma` row's `curve.gamma_rel_max`, `omega_e_rel` from the `omega_e` row, `coupling_log_rms` from the `B` row's `curve.b_rel`; `v0_rms` derived from `omega_e_rel` as `omega_e_rel * D_e / 20` for N₂, stated as such), and `BUDGET_SOURCE = "validation/factory/results/{N2,NO}_budget.md (2026-..)"`.

- [ ] **Step 1: Write the failing test**

```python
# projects/potential_factory/test_report.py
from __future__ import annotations

import csv
from pathlib import Path

from projects.potential_factory.report import BUDGET_SOURCE, Tolerances


def _budget(molecule: str, feature: str, key: str) -> float:
    rows = list(csv.DictReader(Path(f"validation/factory/results/{molecule}_rows.csv").open()))
    ok = [r for r in rows if r["feature"] == feature]
    return float(ok[-2][key]) if len(ok) > 1 else float(ok[-1][key])  # last row within tolerance


def test_tolerances_defaults_equal_the_committed_budget():
    t = Tolerances()
    assert "validation/factory/results" in BUDGET_SOURCE
    assert t.e_res_rms == min(_budget("N2", "E_res", "e_res_rms"), _budget("NO", "E_res", "e_res_rms"))
    assert t.gamma_rel == min(_budget("N2", "Gamma", "gamma_rel_max"), _budget("NO", "Gamma", "gamma_rel_max"))
    assert t.omega_e_rel == min(_budget("N2", "omega_e", "omega_e_rel"), _budget("NO", "omega_e", "omega_e_rel"))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run --no-sync pytest projects/potential_factory/test_report.py -q`
Expected: FAIL (`ImportError: BUDGET_SOURCE` / values differ)

- [ ] **Step 3: Set the defaults**

Edit `Tolerances` in `projects/potential_factory/report.py`: replace the PLACEHOLDER docstring with `"""Per-feature budget measured by validation/factory/sensitivity.py; see docs/physics/potential-factory-budget.md."""`, add `BUDGET_SOURCE = "validation/factory/results/{N2,NO}_budget.md (<date of the run>)"` and set each default to the number from the tables (typed as literals, with a comment naming the row it came from).

- [ ] **Step 4: Write the note**

`docs/physics/potential-factory-budget.md` — house form; contents: the observable-level tolerance (the one hand-chosen set), the catalogue, both budget tables verbatim, the resulting `Tolerances`, the transfer caveat (`l = 2/1` boomerang systems → O₂ fine; re-measure for H₂), the NO-DA caveat, and the sentence that nothing in the study touched an experiment. Add it to `docs/physics/resonances.md`'s toctree; in `docs/physics/potential-factory.md` replace "`Tolerances` are placeholders…" with a pointer to this note.

- [ ] **Step 5: Run tests and commit**

Run: `uv run --no-sync pytest projects/potential_factory -q -m "not slow"`
Expected: all pass

```bash
git add projects/potential_factory/report.py projects/potential_factory/test_report.py docs/physics/potential-factory-budget.md docs/physics/potential-factory.md docs/physics/resonances.md
git commit -m "feat(factory): Tolerances set from the measured N2/NO budget; budget note"
```

---

## Self-review against the spec

- Procedure steps 1–4 of "Tolerance budget" → Tasks 1–3 (perturbations realised as parameter moves and *measured*; exact observables; width-unit metrics; the budget = largest deviation within the observable tolerance). The hand-chosen observable tolerance → `OBSERVABLE_TOL` (Task 3), the stated numbers copied verbatim from the spec.
- Both spec caveats → Task 4's note (transfer to O₂ / re-measure for H₂; NO DA via threshold and peaks only — enforced in `observable_deltas`, which never compares DA magnitudes).
- Oracle-only → Global Constraints; no loader for experimental data exists anywhere in this plan.
- Type consistency: `CurveDeltas(e_res_rms, gamma_rel_max, r_c_shift, omega_e_rel, b_rel)`, `ObservableDeltas(peak_pos_widths, peak_height_rel, integrated_rel, level_pos_widths, level_width_rel, da_threshold_shift)`, `BudgetRow(molecule, perturbation, feature, amount, curve, obs)`, `Perturbation(name, param, feature)` used consistently; two named unknowns (the `presets` accessor, `PeakAlignment`'s distance field, `ResonanceLevels`' level attribute, `ve_cross_section`'s return orientation) are flagged in Task 2 as "read the source and use the real name" rather than guessed.
