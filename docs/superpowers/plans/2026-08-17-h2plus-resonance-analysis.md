# H₂⁺ DR exact-resonance analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Mark the genuine (non-Born-Oppenheimer) resonance positions on the H₂⁺ DR cross sections, measure their shift from the BO quasi-bound levels, and explain the cross-section features that are not resonances as S-matrix zeros.

**Architecture:** Four computations feed one analysis. The library gains the complex DR transition amplitude (currently formed and discarded). A driver builds the BO Rydberg levels from `anion_electronic_states` + `vibrational_states`. `exact_resonance_states` — validated on N₂, never run on a Coulomb model — supplies the exact poles, box-probed before it is sized. A σ(E) sweep on the full deck runs on sadaharu under MUMPS. The exact poles and an independent fit of the swept amplitude are then compared against each other and against the peaks.

**Tech Stack:** numpy/scipy, `qscat.core` (`dissociation`, `resonance`, `plot`), `qscat.viz`, pytest, Docker + MUMPS on sadaharu, `apps/qscat-run` configs.

**Spec:** `docs/superpowers/specs/2026-08-17-h2plus-resonance-analysis-design.md`

## Global Constraints

- Atomic units throughout. Energies in Hartree; the electron energy `E` is measured from the `v_init` ion vibrational threshold, as `dr_cross_section` already defines it.
- Energy windows: `[0, 0.008]`, `[0.010, 0.018]`, `[0.020, 0.027]` Ha — the published panels trimmed short of each ion vibrational threshold. Do not widen them into the accumulation region without saying so.
- `qscat.core` must never import `qscat.model` at runtime (`tests/…test_core_no_model_import.py` enforces it); molecule specifics live in `validation/h2plus/`.
- The full deck (`H2P:emoscat`, ~1.15 M unknowns, electronic → 1300 bohr) is **MUMPS-only and sadaharu-only**. `H2P:proxy` is the laptop reduction and is never a source of quoted physics.
- ECS angle ceilings for this model: nuclear θ < 22.5° (enforced by `IonicResonanceModel.max_nuclear_ecs_angle_deg`), electronic θ < 45°.
- Differential tolerance against a dense/reference computation: `rtol = 1e-9` (the repo's cross-architecture floor for anything through a sparse factorization).
- Widths below ~1e-6 Ha are noise (recorded for the BO levels); positions stay quotable where widths do not.
- Gates before each commit: `uv run ruff check .`, `uv run mypy libs/qscat/qscat`, the touched tests. Foreground only — backgrounded pytest reports exit 0 with empty output here.
- `qscat.viz` tests **skip in GitHub CI** (matplotlib is only in the `plot` extra); run them locally or in Docker before believing a viz change.

---

### Task 1: Cost probe on sadaharu — measure before committing

**Files:**
- Create: `docs/superpowers/plans/2026-08-17-h2plus-cost-probe.md` (the recorded outcome)

**Interfaces:**
- Consumes: `apps/qscat-run/examples/h2p-dr-ti.yaml`, `docker/run.sh`.
- Produces: three numbers every later task depends on — seconds per full-deck energy point, peak RSS, and whether `SparseLU.refactor` reuse holds across an energy sweep at this size.

This task exists because a full-deck energy point at ~1.15 M unknowns has never been timed in this repository, and Task 6's sampling density is chosen from the answer. A coarse sweep through a Rydberg series is worse than no sweep — it aliases the peaks it is meant to locate — so the choice must be made from a measurement, not a guess.

- [ ] **Step 1: Write a two-energy probe config**

Create `apps/qscat-run/examples/h2p-dr-probe.yaml`:

```yaml
# Cost probe: two energies on the FULL deck, to measure per-energy cost and
# confirm the analysis-reuse path holds at ~1.15M unknowns. Not a physics run.
molecule: H2P
methods: [ti]
observables:
  - {kind: dr, channels: 3}
v_init: 0
energies: {min: 0.0120, max: 0.0121, step: 0.0001}
grid: {preset: emoscat}
artifacts:
  cross_section: true
backend: mumps
output_dir: runs/h2p-dr-probe
```

- [ ] **Step 2: Run it on sadaharu**

This is an operator step — the full deck does not run on the Mac. On sadaharu, from the repo root:

```bash
docker/build.sh test
/usr/bin/time -v docker/run.sh apps/qscat-run/examples/h2p-dr-probe.yaml runs/h2p-dr-probe
```

Record from the output: wall time, the `Maximum resident set size` line, and (from the run log) whether the second energy was faster than the first.

- [ ] **Step 3: Record the numbers and the decision**

Write `docs/superpowers/plans/2026-08-17-h2plus-cost-probe.md` with: seconds for energy 1, seconds for energy 2, peak RSS, and the resulting sampling budget for Task 6 computed as `(0.023 Ha of windows) / (points affordable in the time budget)`. State the budget explicitly, e.g. "at 40 s/point and a 12 h budget, 1080 points → 2.1e-5 Ha spacing".

If a point costs more than ~2 minutes, record the fallback decision from the spec: narrow to the two windows with the richest structure rather than thinning resolution everywhere.

- [ ] **Step 4: Commit**

```bash
git add apps/qscat-run/examples/h2p-dr-probe.yaml docs/superpowers/plans/2026-08-17-h2plus-cost-probe.md
git commit -m "perf(h2plus): measure the full-deck per-energy cost on sadaharu"
```

---

### Task 2: The DR transition amplitude (library)

**Files:**
- Modify: `libs/qscat/qscat/core/dissociation.py` (the `dr_cross_section` overloads and body)
- Test: `libs/qscat/tests/test_dissociation.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `dr_cross_section(..., return_amplitude=True)` returning `(sigma, amp)` with `amp` complex, shaped exactly like `sigma` (`(n_channels,)` for scalar `E`, `(len(E), n_channels)` for an array). Composable with `return_wavefunction`, in which case the return is `(sigma, psi, amp)`. Tasks 6, 7 and 9 consume it.

Note on naming: the quantity returned is the **transition amplitude** `t` the solver already forms, with `sigma = 4*pi**3 * abs(t)**2 / (2*E)`. It is *not* literally a unitary S-matrix element — the thesis's `S_DR` differs from it by the standard `S = -2*pi*i*T` factor and its own normalization. That factor is a fixed rotation and rescale, so it changes neither the zeros nor the shape of the real/imaginary crossings that Fig. 4.10 is about. Returning the amplitude as computed is honest; synthesizing an "S" would mean guessing a normalization.

- [ ] **Step 1: Write the failing test**

Append to `libs/qscat/tests/test_dissociation.py`:

```python
def test_dr_amplitude_reproduces_the_returned_sigma() -> None:
    """The amplitude and sigma must not drift apart: sigma is 4pi^3|t|^2/2E."""
    from qscat.core.dissociation import dr_cross_section

    tgrid, eps, chi = _h2p_small_problem()  # existing helper in this module
    energies = np.array([0.012, 0.014])
    sigma, amp = dr_cross_section(
        H2P, tgrid, eps=eps, chi=chi, E=energies, v_init=0,
        n_channels=3, return_amplitude=True,
    )
    assert amp.shape == sigma.shape
    assert amp.dtype == np.complex128
    recomputed = 4.0 * np.pi**3 * np.abs(amp) ** 2 / (2.0 * energies[:, None])
    open_channels = sigma > 0.0
    assert np.allclose(recomputed[open_channels], sigma[open_channels], rtol=1e-12)
    # A closed channel contributes exactly zero to both.
    assert np.all(amp[~open_channels] == 0.0)


def test_dr_amplitude_composes_with_the_wavefunction_return() -> None:
    from qscat.core.dissociation import dr_cross_section

    tgrid, eps, chi = _h2p_small_problem()
    sigma, psi, amp = dr_cross_section(
        H2P, tgrid, eps=eps, chi=chi, E=0.012, v_init=0, n_channels=3,
        return_wavefunction=True, return_amplitude=True,
    )
    assert psi is not None and psi.shape == (tgrid.size,)
    assert amp.shape == sigma.shape
```

If `_h2p_small_problem()` does not exist in that module, build the fixture from `qscat_run.presets` `H2P:proxy` and `vibrational_states`, mirroring whatever the neighbouring DR tests already do — read them first and follow that pattern rather than inventing a second one.

- [ ] **Step 2: Run to verify it fails**

Run: `uv run --no-sync pytest libs/qscat/tests/test_dissociation.py -q -k amplitude`
Expected: FAIL — `TypeError: dr_cross_section() got an unexpected keyword argument 'return_amplitude'`.

- [ ] **Step 3: Implement**

In `libs/qscat/qscat/core/dissociation.py`, add the parameter to the signature and the four overloads (`return_wavefunction` × `return_amplitude`), then in the body:

```python
    out = np.zeros((len(e_arr), n_channels), dtype=np.float64)
    amp = np.zeros((len(e_arr), n_channels), dtype=np.complex128)
```

inside the channel loop, right where `t` is formed:

```python
            t = c_product(phi_f, v_psi)
            amp[ie, n] = t
            out[ie, n] = 4.0 * np.pi**3 * abs(t) ** 2 / (2.0 * float(e))
```

and at the return:

```python
    sigma = np.asarray(out[0] if scalar else out, dtype=np.float64)
    amplitude = np.asarray(amp[0] if scalar else amp, dtype=np.complex128)
    psi_out = psi_list[0] if scalar else psi_list
    if return_wavefunction and return_amplitude:
        return sigma, psi_out, amplitude
    if return_amplitude:
        return sigma, amplitude
    if return_wavefunction:
        return sigma, psi_out
    return sigma
```

Extend the docstring with a `return_amplitude` paragraph stating the `sigma = 4*pi**3*abs(t)**2/(2*E)` relation and the `S = -2*pi*i*T` caveat above.

- [ ] **Step 4: Run to verify it passes**

Run: `uv run --no-sync pytest libs/qscat/tests/test_dissociation.py -q`
Expected: PASS, including the pre-existing DR tests.

- [ ] **Step 5: Gates and commit**

```bash
uv run ruff check . && uv run mypy libs/qscat/qscat
git add libs/qscat/qscat/core/dissociation.py libs/qscat/tests/test_dissociation.py
git commit -m "feat(core): return the DR transition amplitude from dr_cross_section"
```

---

### Task 3: Born-Oppenheimer Rydberg levels `ω_i^j`

**Files:**
- Create: `validation/h2plus/rydberg_levels.py`
- Test: `validation/h2plus/test_rydberg_levels.py`

**Interfaces:**
- Consumes: `qscat.core.anion_electronic_states(g_r, model, R_inf, n_states)` → `(eps, phi)`; `qscat.core.vibrational_states(grid, mu, n, v0)` → `VibrationalBasis` with `.energies`.
- Produces: `rydberg_levels(model, g_r, g_R, *, n_curves, n_vib) -> RydbergLevels` — a frozen dataclass with `curves` (shape `(n_curves, g_R.n)`, the `E_Ryn(R)` on the nuclear grid) and `energies` (shape `(n_curves, n_vib)`, real). Tasks 5, 9 and 10 consume it.

These levels are the thesis's dashed verticals. They are **real** in the BO picture: a Rydberg state is bound in its own electronic curve, and the width it really has comes from the coupling to the dissociative channel that BO discards. That absence is exactly what Task 5 measures against.

- [ ] **Step 1: Write the failing test**

```python
"""BO quasi-bound levels in the H2+ Rydberg curves."""
from __future__ import annotations

import numpy as np
from qscat.model import H2P
from validation.h2plus.config import proxy_grid
from validation.h2plus.rydberg_levels import rydberg_levels


def test_curves_are_ordered_and_below_the_ion() -> None:
    tg = proxy_grid()
    res = rydberg_levels(H2P, tg.grids[0], tg.grids[1], n_curves=3, n_vib=4)
    assert res.curves.shape == (3, tg.grids[1].n)
    # Rydberg curves are ordered in energy at every R.
    real = np.abs(tg.grids[1].points.imag) < 1e-12
    ordered = np.diff(res.curves[:, real].real, axis=0)
    assert np.all(ordered > 0.0)


def test_levels_are_real_and_ascending() -> None:
    tg = proxy_grid()
    res = rydberg_levels(H2P, tg.grids[0], tg.grids[1], n_curves=2, n_vib=5)
    assert res.energies.shape == (2, 5)
    assert np.all(np.isreal(res.energies))
    assert np.all(np.diff(res.energies, axis=1) > 0.0)  # vibrational ladder
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run --no-sync pytest validation/h2plus/test_rydberg_levels.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'validation.h2plus.rydberg_levels'`.

- [ ] **Step 3: Implement**

```python
"""Born-Oppenheimer quasi-bound levels in the H2+ Rydberg curves.

The neutral's bound electronic states at fixed R form the Rydberg series
`E_Ryn(R)`; the vibrational levels those curves support are the quasi-bound
states a DR cross-section peak is conventionally assigned to. They carry NO
width here: a Rydberg state is bound in its own curve, and the width it really
has comes from the coupling to the dissociative channel the BO picture drops.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.core import anion_electronic_states, vibrational_states
from qscat.dvr import FemDvrEcsGrid


@dataclass(frozen=True)
class RydbergLevels:
    curves: npt.NDArray[np.complex128]   # (n_curves, n_R)
    energies: npt.NDArray[np.float64]    # (n_curves, n_vib)


def rydberg_levels(
    model, g_r: FemDvrEcsGrid, g_R: FemDvrEcsGrid, *, n_curves: int, n_vib: int
) -> RydbergLevels:
    """`E_Ryn(R)` and the vibrational levels each curve supports."""
    pts = g_R.points
    curves = np.empty((n_curves, pts.size), dtype=np.complex128)
    for j, R in enumerate(pts):
        eps, _ = anion_electronic_states(
            g_r=g_r, model=model, R_inf=complex(R), n_states=n_curves
        )
        curves[:, j] = eps

    energies = np.empty((n_curves, n_vib), dtype=np.float64)
    for n in range(n_curves):
        curve = curves[n]

        def v_n(R, _curve=curve):
            # The curve is tabulated ON g_R.points, and vibrational_states
            # evaluates its potential at exactly those points, so this is a
            # lookup rather than an interpolation.
            return _curve

        basis = vibrational_states(g_R, model.mu, n_vib, v_n)
        energies[n] = np.asarray(basis.energies, dtype=np.float64)

    return RydbergLevels(curves=curves, energies=energies)
```

If `vibrational_states` calls `v0` with an argument that is not exactly `g_R.points`, replace the lookup with an interpolation of `curve` over `R` and say so in a comment — do not silently return a mismatched array.

- [ ] **Step 4: Run to verify it passes**

Run: `uv run --no-sync pytest validation/h2plus/test_rydberg_levels.py -q`
Expected: PASS.

- [ ] **Step 5: Check against the published positions**

The thesis marks these levels on Fig. 4.7 (p. 70). Compute them on the **full** deck's nuclear grid and print the ones falling in each of the three windows:

```bash
uv run python -c "
from qscat.model import H2P
from validation.h2plus.config import full_grid
from validation.h2plus.rydberg_levels import rydberg_levels
tg = full_grid()
r = rydberg_levels(H2P, tg.grids[0], tg.grids[1], n_curves=3, n_vib=14)
print(r.energies)
"
```

Compare the pattern against the published dashed verticals: consecutive panels are offset by roughly one ion vibrational quantum (~0.0095 Ha). **This is what settles the `ω_i^j` index convention** — record which index is the Rydberg curve and which the vibrational quantum number in the module docstring, as a measured fact.

If the positions do not reproduce, stop and report: either the curve construction or the level assignment is wrong, and every downstream seed depends on it.

- [ ] **Step 6: Commit**

```bash
uv run ruff check . && uv run --no-sync pytest validation/h2plus/test_rydberg_levels.py -q
git add validation/h2plus/rydberg_levels.py validation/h2plus/test_rydberg_levels.py
git commit -m "feat(h2plus): BO quasi-bound levels in the Rydberg curves"
```

---

### Task 4: Pole box-convergence probe

**Files:**
- Create: `validation/h2plus/pole_box_probe.py`

**Interfaces:**
- Consumes: `rydberg_levels` (Task 3); `qscat.core.exact_resonance_states(model, grid_base, grid_electronic, grid_nuclear, *, shifts, window, k, rel_tol, atol)`.
- Produces: the smallest converged electronic box for the pole search, recorded in the module docstring and consumed by Task 5.

The cross-section deck needs 1300 bohr because the incident Coulomb wave must be represented. A Rydberg resonance is a **closed-channel, localized** state whose extent goes as `n²` (~120 bohr at `n = 11`), so it may converge on a far smaller box. This task measures that instead of assuming it — and if it fails, that is the finding, and the campaign falls back to the `S(E)`-fitting route alone.

- [ ] **Step 1: Write the probe**

```python
"""Does the H2+ pole search converge on a small electronic box?

The cross-section deck needs 1300 bohr for the incident Coulomb wave. A Rydberg
resonance is closed-channel and localized, so it may not. Measured here rather
than assumed, because the answer sets the cost of the whole pole campaign.
"""

from __future__ import annotations

import time

import numpy as np
from qscat.core import exact_resonance_states
from qscat.core.grids import fem_grid_exp_tail, nuclear_grid
from qscat.dvr import TensorGrid
from qscat.model import H2P
from validation.h2plus.config import full_grid
from validation.h2plus.rydberg_levels import rydberg_levels

WINDOW = (0.0, 0.05, -0.01, 0.0)   # widened in Step 2 if it catches nothing


def main() -> None:
    tg_full = full_grid()
    lv = rydberg_levels(H2P, tg_full.grids[0], tg_full.grids[1], n_curves=3, n_vib=14)
    # Seeds: BO levels inside the middle window, where the structure is cleanest.
    flat = np.sort(lv.energies.ravel())
    seeds = [complex(e, -1e-4) for e in flat if 0.010 <= e <= 0.018][:4]
    print(f"seeds: {[f'{s:.6f}' for s in seeds]}", flush=True)

    for r_max in (150.0, 300.0, 600.0):
        el_a = fem_grid_exp_tail(r_max=r_max, angle_deg=35.0)
        el_b = fem_grid_exp_tail(r_max=r_max, angle_deg=42.0)
        nu_a = nuclear_grid(angle_deg=20.0, r_max=14.0, n_complex=6, quadrature=12)
        nu_b = nuclear_grid(angle_deg=15.0, r_max=14.0, n_complex=6, quadrature=12)
        t0 = time.perf_counter()
        res = exact_resonance_states(
            H2P,
            TensorGrid([el_a, nu_a]), TensorGrid([el_b, nu_a]), TensorGrid([el_a, nu_b]),
            shifts=seeds, k=8, window=WINDOW,
        )
        print(
            f"r_max={r_max:6.0f}  n2d={el_a.n * nu_a.n:8d}  "
            f"found={res.energies.size}  {time.perf_counter() - t0:.0f}s",
            flush=True,
        )
        for e, g, re_, rn in zip(
            res.energies, res.widths, res.residual_electronic,
            res.residual_nuclear, strict=True,
        ):
            print(f"    E={e.real:+.9f}  G={g:.3e}  res_el={re_:.1e}  res_nuc={rn:.1e}",
                  flush=True)


if __name__ == "__main__":
    main()
```

Check `fem_grid_exp_tail`'s actual signature before running (it is the H₂⁺ builder used by `validation/h2plus/config.py`) and pass the arguments it wants — the call above assumes `r_max` and `angle_deg` only.

- [ ] **Step 2: Run it and read the result**

Run: `uv run python -m validation.h2plus.pole_box_probe`

Three outcomes, each with a defined response:
- **Poles converge by 300 bohr** (position stable to ≲1e-6 Ha between 300 and 600, residuals ≤1e-5): use that box in Task 5 and record it.
- **Still moving at 600 bohr**: escalate to sadaharu and probe 1200 bohr; if it still moves, the pole campaign needs the full deck and Task 5 becomes a sadaharu job.
- **Nothing found in the window**: widen `WINDOW`'s imaginary extent (these states may be far narrower than N₂'s — try `im_lo = -1e-5`), then raise `k`. If still empty, report it: `exact_resonance_states` has never run on a Coulomb model and a negative result here redirects the campaign.

- [ ] **Step 3: Record and commit**

Write the outcome into the module docstring as measured fact (box chosen, positions, residuals, timings).

```bash
git add validation/h2plus/pole_box_probe.py
git commit -m "test(h2plus): box-convergence probe for the Rydberg pole search"
```

---

### Task 5: Pole campaign across the three windows

**Files:**
- Create: `validation/h2plus/exact_poles.py`

**Interfaces:**
- Consumes: `rydberg_levels` (Task 3); the box chosen by Task 4; `exact_resonance_states`.
- Produces: `poles.npz` with `energies` (complex), `widths`, `residual_electronic`, `residual_nuclear`, `seed_level` (the `(curve, vib)` index each pole was seeded from), plus `exact_poles(window)` returning the same as a dataclass. Tasks 7, 8, 9 and 10 consume it.

- [ ] **Step 1: Write the driver**

Seed every BO level in each window, run `exact_resonance_states` per window, cache to `.npz` (gitignored) so downstream figure work does not re-solve. Structure it exactly like `validation/n2/exact_resonance_figures.py`'s cache block — read that file and copy the pattern, including `main(cache=...)` and the "delete the file to force a recompute" comment.

Record per pole which BO level seeded it; Task 9's shift table is meaningless without that pairing, and pairing by nearest-energy after the fact fails exactly where the shift is largest.

- [ ] **Step 2: Run for all three windows**

Run: `uv run python -m validation.h2plus.exact_poles`

Expected: one pole per BO seed in the clean parts of each window, with `residual_electronic` and `residual_nuclear` both small. Fewer poles than seeds near the top of a window is expected — that is the Rydberg accumulation the spec warns about, and it must be reported per window, not silently dropped.

- [ ] **Step 3: Commit**

```bash
git add validation/h2plus/exact_poles.py
git commit -m "feat(h2plus): exact 2-D resonance poles across the three DR windows"
```

---

### Task 6: The σ(E) sweep on sadaharu

**Files:**
- Create: `apps/qscat-run/examples/h2p-dr-windows.yaml`

**Interfaces:**
- Consumes: Task 1's sampling budget; Task 2's `return_amplitude` (through qscat-run's DR observable — check whether the runner threads it, and if not, add the artifact there as part of this task).
- Produces: `sigma(E)` and the complex amplitude for DR₀…DR₂ across the three windows, as npz artifacts. Tasks 7, 9 and 10 consume them.

- [ ] **Step 1: Write three window configs**

One config per window (separate runs keep a failure from costing all three). For the first:

```yaml
# H2+ DR, window 1 of 3. Full deck -> sadaharu/MUMPS only.
molecule: H2P
methods: [ti]
observables:
  - {kind: dr, channels: 3}
v_init: 0
energies: {min: 0.0, max: 0.008, step: 2.0e-5}   # step from Task 1's budget
grid: {preset: emoscat}
artifacts:
  cross_section: true
  amplitude: true        # add to the runner if absent (see Step 2)
backend: mumps
output_dir: runs/h2p-dr-w1
```

Replace `step` with the value Task 1 recorded. Do not invent it.

- [ ] **Step 2: Thread the amplitude through qscat-run if needed**

Check `apps/qscat-run/qscat_run/` for where the `dr` observable calls `dr_cross_section`. If the amplitude is not persisted, add it as an artifact alongside the cross section, with a test in `apps/qscat-run/tests/` asserting the npz carries a complex array shaped like the σ array.

- [ ] **Step 3: Run on sadaharu**

```bash
docker/run.sh apps/qscat-run/examples/h2p-dr-w1.yaml runs/h2p-dr-w1
docker/run.sh apps/qscat-run/examples/h2p-dr-w2.yaml runs/h2p-dr-w2
docker/run.sh apps/qscat-run/examples/h2p-dr-w3.yaml runs/h2p-dr-w3
```

Copy the `runs/` npz artifacts back to the working machine for the analysis tasks.

- [ ] **Step 4: Commit the configs**

```bash
git add apps/qscat-run/examples/h2p-dr-w*.yaml
git commit -m "feat(h2plus): full-deck DR sweep configs for the three windows"
```

---

### Task 7: Independent pole extraction from the amplitude

**Files:**
- Create: `validation/h2plus/pole_fits.py`
- Test: `validation/h2plus/test_pole_fits.py`

**Interfaces:**
- Consumes: the swept amplitude from Task 6; the exact poles from Task 5.
- Produces: `fit_poles(E, amplitude) -> ndarray of complex` and a comparison table against Task 5.

This is what makes "the poles sit on the peaks" a check rather than a tautology: the exact poles come from an eigenproblem that never sees the cross section; these come from the cross section and never see the eigenproblem.

- [ ] **Step 1: Write the failing test on synthetic data**

```python
def test_fit_recovers_a_planted_pole() -> None:
    """A synthetic Breit-Wigner amplitude with a known pole must be recovered."""
    import numpy as np
    from validation.h2plus.pole_fits import fit_poles

    e_pole, gamma = 0.0123456, 4.0e-6
    E = np.linspace(e_pole - 20 * gamma, e_pole + 20 * gamma, 400)
    amp = 1.0 / (E - (e_pole - 0.5j * gamma)) + 0.3  # resonance + background
    got = fit_poles(E, amp)
    assert got.size >= 1
    best = got[int(np.argmin(np.abs(got.real - e_pole)))]
    assert abs(best.real - e_pole) < 0.05 * gamma
    assert abs(-2.0 * best.imag - gamma) < 0.1 * gamma
```

- [ ] **Step 2: Run to verify it fails, then implement**

Run: `uv run --no-sync pytest validation/h2plus/test_pole_fits.py -q` → FAIL (module missing).

Implement `fit_poles` by fitting `a/(E - E_pole) + b + c*E` to the complex amplitude in a window around each local maximum of `|amp|`, via `scipy.optimize.least_squares` on the real and imaginary parts jointly. Complex least squares on the amplitude — not on `|amp|²` — because the phase carries the pole position and squaring discards it.

- [ ] **Step 3: Apply to the real sweep and compare**

Report, per window: fitted pole vs Task 5's exact pole, in meV, and the widths. Agreement to within the sampling step is success; a systematic offset is a finding to investigate before publishing either number.

- [ ] **Step 4: Commit**

```bash
git add validation/h2plus/pole_fits.py validation/h2plus/test_pole_fits.py
git commit -m "feat(h2plus): extract DR poles from the swept amplitude"
```

---

### Task 8: Assignment by overlap

**Files:**
- Create: `validation/h2plus/assignment.py`

**Interfaces:**
- Consumes: `dr_cross_section(..., return_wavefunction=True)`; Task 5's poles and states; `qscat.linalg.c_product`.
- Produces: a table of `abs(c_product(psi_res, psi_plus(E)))` for each pole against `Ψ⁺` at energies on and off the corresponding feature — the Table 4.2 analogue with exact states.

- [ ] **Step 1: Write the driver**

For each pole, take three energies: at the pole's `Re E`, and at `Re E ± 5Γ`. Solve `dr_cross_section(..., E=those, return_wavefunction=True)` on the full deck (sadaharu) or the Task-4 box if it proved sufficient, then overlap each `Ψ⁺` with each pole's state under the c-product, normalizing both.

- [ ] **Step 2: Check the expected behaviour**

The overlap must **peak at the pole and fall off it** — that is the prediction from `Ψ_sc ≈ ψ_res·<ψ̃_res|V|Ψ_i>/(E − E_res)`. A flat overlap means the state is not what is driving that feature, and the assignment fails; report it rather than presenting the table as if it worked.

- [ ] **Step 3: Commit**

```bash
git add validation/h2plus/assignment.py
git commit -m "feat(h2plus): assign DR features to exact resonance states by overlap"
```

---

### Task 9: Marked cross-section plotting (library)

**Files:**
- Modify: `libs/qscat/qscat/core/plot.py` (`plot_cross_sections`)
- Test: `libs/qscat/tests/test_plot_markers.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `plot_cross_sections(..., markers=...)` where `markers` is `dict[str, tuple[Sequence[float], dict]]` — a label mapped to positions plus a matplotlib style dict, drawn as labelled verticals. Task 10 consumes it.

- [ ] **Step 1: Write the failing test**

```python
def test_markers_draw_one_vertical_per_position(tmp_path) -> None:
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from qscat.core import plot_cross_sections

    E = np.linspace(0.0, 0.01, 50)
    sigma = np.ones((50, 2))
    out = tmp_path / "m.png"
    plot_cross_sections(
        E, sigma, path=out,
        markers={"exact": ([0.002, 0.004], {"color": "k"}),
                 "BO": ([0.003], {"color": "r", "linestyle": "--"})},
    )
    assert out.exists() and out.stat().st_size > 0
    plt.close("all")
```

- [ ] **Step 2: Run to verify it fails, implement, re-run**

Implement by looping the mapping and calling `ax.axvline(x, **style)` for each position, labelling only the first of each group so the legend carries one entry per marker set. Keep it generic — no physics, matching the module's existing contract.

- [ ] **Step 3: Commit**

```bash
uv run ruff check . && uv run mypy libs/qscat/qscat
git add libs/qscat/qscat/core/plot.py libs/qscat/tests/test_plot_markers.py
git commit -m "feat(core): labelled vertical markers on plot_cross_sections"
```

---

### Task 10: Figures, note, and reference rows

**Files:**
- Create: `validation/h2plus/resonance_analysis.py`, `docs/physics/h2plus-resonance-analysis.md`
- Create: `docs/physics/figures/h2plus-dr-resonances-w{1,2,3}.png`, `h2plus-dr-smatrix-zeros.png`, `h2plus-exact-vs-bo-shift.png`
- Modify: `reference/literature/vana-2017-thesis.md`, `CHANGELOG.md`, `CLAUDE.md`, `docs/index.md`, `docs/physics/README.md`, `docs/api/core.md`

**Interfaces:**
- Consumes: everything above.
- Produces: the deliverables.

- [ ] **Step 1: The three annotated panels**

`σ_DR0` and `σ_DR1` on a log axis per window, exact poles as solid verticals, BO `ω_i^j` as dashed, thresholds shaded — the Fig. 4.7 analogue via Task 9's `markers`.

- [ ] **Step 2: The S-matrix zero panels**

For each cross-section minimum with no pole nearby: zoomed `σ_DR1` above `Re amp` / `Im amp` on a linear scale, showing both crossing zero — the Fig. 4.10 analogue. Identify candidates as local minima of `σ_DR1` whose distance to the nearest pole exceeds several widths.

- [ ] **Step 3: The shift figure and table**

Exact `Re E` minus its seeding BO level, in meV, per level, with `Γ` alongside. State plainly that BO carries no width here, so widths have no BO counterpart — they are new information, not a correction.

- [ ] **Step 4: The physics note**

`docs/physics/h2plus-resonance-analysis.md`: method, what each figure shows, the shift table, the assignment table, and a Limits section carrying the accumulation caveat, the ~1e-6 Ha width floor, the ECS angle ceiling, and whichever box Task 4 settled on. Wire it into `docs/index.md`'s Theory toctree and `docs/physics/README.md`, or the `-W` build fails on an orphaned document.

- [ ] **Step 5: Reference-note rows**

Add Fig. 4.3 (p. 64), Fig. 4.7 (p. 70), Fig. 4.10 (p. 73) and Table 4.2 (p. 73) to `reference/literature/vana-2017-thesis.md`'s "What this repository uses" table, each with what we take from it — the `mastering-references` rule is that every published fact the repo relies on carries a locator.

- [ ] **Step 6: Full gates and commit**

```bash
uv run ruff check . && uv run mypy libs/qscat/qscat
uv run --no-sync pytest libs/qscat/tests -m "not slow" -n 8 -q
uv run --no-sync pytest tests/ validation/h2plus -q
uv run --no-sync sphinx-build -b html -W --keep-going docs docs/_build/html
git add -A
git commit -m "docs(h2plus): exact resonances against the DR cross sections"
```

---

## Self-review notes

- **Spec coverage:** S-matrix element → Task 2; BO levels → Task 3; exact poles → Tasks 4–5; independent fits → Task 7; assignment → Task 8; figures A/B/C → Task 10 (with the plotting primitive in Task 9); cost strategy → Task 1; reference rows → Task 10 Step 5. Every spec section maps to a task.
- **Naming deviation, deliberate:** the spec and the approved question said `return_smatrix`; the plan implements `return_amplitude`, because what the solver forms is the transition amplitude `t` with `sigma = 4*pi**3*abs(t)**2/(2*E)`, not a unitary S element. The `S = -2*pi*i*T` factor is a fixed rotation that preserves the zeros Fig. 4.10 is about. Synthesizing an "S" would mean guessing a normalization.
- **Consistent names across tasks:** `rydberg_levels` / `RydbergLevels.energies`, `exact_poles`, `fit_poles`, `return_amplitude`, `markers`.
- **Two tasks can return negative results** and say so explicitly rather than being forced green: Task 4 (poles may need the full deck, or may not be findable at all on a Coulomb model) and Task 8 (a flat overlap means the assignment fails).
