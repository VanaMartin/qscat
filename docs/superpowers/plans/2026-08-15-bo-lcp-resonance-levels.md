# BO/LCP Resonance Levels Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Diagonalize the nuclear problem in the LCP complex potential
`V_res(R) = E_res(R) − (i/2)Γ(R)` to obtain the quasi-bound (Siegert) vibrational
levels `E_v − iΓ_v/2` of the molecular anion — the Born-Oppenheimer approximation
to the 2-D model's resonance energies — and expose them through `qscat-run`.

**Architecture:** Step 1 (the complex curve) already exists as
`qscat.core.lcp.local_complex_potential`. This plan adds step 2: a
complex-symmetric FEM-DVR-ECS diagonalization of `T(mu) + diag(V_res)` on the
nuclear grid, with physical levels selected by a two-nuclear-ECS-angle stability
criterion (the vectorized generalization of `qscat.ecs.find_resonance_pole`).
Because `E_res(R)` at real `R` is angle-independent, the expensive electronic pole
walk runs once and only the tail assembly plus one extra nuclear diagonalization
are repeated. A golden-rule comparator (`Γ = 0` levels plus
`⟨χ_v|Γ|χ_v⟩`) rides along for free and reproduces what eMoScat/the thesis
actually computed.

**Tech Stack:** Python ≥3.12, numpy/scipy, `qscat.dvr` (FEM-DVR-ECS grids,
`kinetic`, `eigen`), `qscat.ecs` (pole matching), `qscat.linalg.c_product`,
pytest, `apps/qscat-run` (YAML → artifacts CLI), `uv` workspace.

## Global Constraints

- **Atomic units throughout.** Energies in Hartree, lengths in bohr. No unit
  conversions in method code (`libs/qscat/qscat/units.py` is the only source).
- **`qscat.core` never imports `qscat.model` or `projects.*` at runtime.** `model`
  is typed against the `ResonanceModel` protocol under `TYPE_CHECKING` only,
  exactly as `driven.py`/`dissociation.py`/`lcp.py` already do. Enforced by
  `libs/qscat/tests/test_core_no_model_import.py`.
- **Complex-safe everywhere.** Never coerce an array to real — ECS tails are
  complex. Use `np.complex128`.
- **Bilinear c-product, never conjugated.** `qscat.linalg.c_product` for every
  inner product in the rotated region. Never `np.vdot`/`.conj()`.
- **Package-absolute imports only.** e.g. `from qscat.core.lcp import ...`. No
  `sys.path` manipulation, no bare intra-directory imports.
- **Naming follows the thesis.** `V_res`/`E_res`/`Γ`/`ω_j` in docstrings;
  qscat's `Vd` **is** the thesis's `E_res`. Call the states *complex-scaled (ECS)
  resonance eigenstates* or *quasi-bound levels* — **never** "Siegert
  pseudostates", which is a different construction (Hvizdoš et al. App. A).
- **Reduced masses are model-supplied.** `model.mu`. Never hard-code.
- Run tests with `uv run --no-sync pytest`. Lint with `uv run ruff check .` and
  `uv run ruff format`.
- Spec: `docs/superpowers/specs/2026-08-15-bo-lcp-resonance-levels-design.md`.

## File Structure

| File | Responsibility |
|---|---|
| `libs/qscat/qscat/ecs/pole.py` (modify) | Add `match_angle_stable`, the multi-state two-angle matcher; `find_resonance_pole` keeps its behaviour, both share one private pairing helper. |
| `libs/qscat/qscat/ecs/__init__.py` (modify) | Export `match_angle_stable`. |
| `libs/qscat/qscat/core/lcp.py` (modify) | Extract `_assemble_lcp`; add `ResonanceLevels`, `lcp_resonance_levels` (numeric core), `resonance_levels` (model convenience). |
| `libs/qscat/qscat/model/ionic.py` (modify) | Add `max_nuclear_ecs_angle_deg = 22.5` (the H₂⁺ divergence bound). |
| `libs/qscat/tests/test_ecs_pole_match.py` (create) | `match_angle_stable` unit tests. |
| `libs/qscat/tests/test_lcp_resonance_levels.py` (create) | Analytic Morse + constant-Γ oracles, convergence, golden-rule consistency, Γ support condition. |
| `apps/qscat-run/qscat_run/config.py` (modify) | `resonance_levels` observable kind; `artifacts.resonance_levels` flag; `grid.nuclear_angle_b`; `energies` optional for levels-only runs. |
| `apps/qscat-run/qscat_run/presets.py` (modify) | Angle-parameterized nuclear decks + `resolve_lcp_nuclear_grid_b`; `resonance_levels` in `VALIDITY`. |
| `apps/qscat-run/qscat_run/runner.py` (modify) | `ResonanceLevelsRun` dataclass; compute levels on the LCP path (both entry points). |
| `apps/qscat-run/qscat_run/artifacts.py` (modify) | `resonance_levels.csv` / `.npz` / `.png` writers. |
| `apps/qscat-run/examples/f2-resonance-levels.yaml` (create) | Levels-only example config. |
| `docs/physics/lcp-resonance-levels.md` (create) | The physics note. |
| `CLAUDE.md` (modify) | Repo-map entry for the new capability. |

---

### Task 1: `match_angle_stable` — the multi-state two-angle matcher

**Files:**
- Modify: `libs/qscat/qscat/ecs/pole.py`
- Modify: `libs/qscat/qscat/ecs/__init__.py`
- Test: `libs/qscat/tests/test_ecs_pole_match.py` (create)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `qscat.ecs.match_angle_stable(eigs_a, eigs_b, window, *, rel_tol=1e-4,
  atol=1e-8) -> tuple[NDArray[complex128], NDArray[float64], NDArray[intp]]`
  returning `(energies, residuals, indices_into_filtered_a_mapped_back_to_eigs_a)`.
  Task 3 uses all three.

**Background for the implementer:** an ECS-rotated Hamiltonian's *discretized
continuum* eigenvalues move when you change the rotation angle `theta`; a true
resonance pole does not. Diagonalizing the same physical Hamiltonian at two
angles and keeping eigenvalues that agree between the two therefore isolates the
physical states. `find_resonance_pole` already does this for a *single* pole
(global argmin). This task generalizes it to *every* stable state in a window,
which is what a vibrational level spectrum needs.

- [ ] **Step 1: Write the failing tests**

Create `libs/qscat/tests/test_ecs_pole_match.py`:

```python
"""Unit tests for `qscat.ecs.match_angle_stable` (multi-state two-angle matcher)."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.ecs import find_resonance_pole, match_angle_stable

WIDE = (-10.0, 10.0, -10.0, 10.0)


def test_keeps_states_present_in_both_spectra():
    # Three "stable" states appear in both spectra to ~1e-12; the rest are
    # "continuum" and differ between the two by O(0.1).
    stable = np.array([-1.0 - 0.01j, -0.5 - 0.02j, 0.25 - 0.05j])
    a = np.concatenate([stable, np.array([2.0 - 1.0j, 3.0 - 1.5j])])
    b = np.concatenate([stable + 1e-12, np.array([2.4 - 1.3j, 3.6 - 1.9j])])

    energies, residuals, idx = match_angle_stable(a, b, WIDE)

    assert energies.shape == (3,)
    assert np.allclose(energies, stable, atol=1e-9)
    assert np.all(residuals < 1e-9)
    # `idx` indexes into the ORIGINAL `eigs_a`, not a filtered copy.
    assert np.allclose(np.asarray(a)[idx], stable, atol=1e-9)


def test_window_excludes_out_of_range_states():
    a = np.array([-1.0 - 0.01j, 5.0 - 0.01j])
    b = a + 1e-12
    energies, _residuals, _idx = match_angle_stable(a, b, (-2.0, 0.0, -1.0, 1.0))
    assert energies.shape == (1,)
    assert np.isclose(energies[0].real, -1.0)


def test_results_are_sorted_by_ascending_real_part():
    stable = np.array([0.25 - 0.05j, -1.0 - 0.01j, -0.5 - 0.02j])
    energies, _residuals, _idx = match_angle_stable(stable, stable + 1e-12, WIDE)
    assert np.all(np.diff(energies.real) > 0)


def test_relative_tolerance_scales_with_magnitude():
    # |E| = 100, partner off by 1e-3 -> |dE|/|E| = 1e-5 < rel_tol=1e-4: accepted.
    a = np.array([100.0 + 0.0j])
    b = np.array([100.001 + 0.0j])
    assert match_angle_stable(a, b, WIDE, rel_tol=1e-4)[0].size == 1
    # Same absolute gap at |E| = 0.1 is a 1e-2 relative gap: rejected. The
    # absolute floor `atol` must not rescue it either.
    a2 = np.array([0.1 + 0.0j])
    b2 = np.array([0.101 + 0.0j])
    assert match_angle_stable(a2, b2, WIDE, rel_tol=1e-4, atol=1e-8)[0].size == 0


def test_empty_window_raises():
    a = np.array([1.0 + 0.0j])
    with pytest.raises(ValueError, match="contains no eigenvalues"):
        match_angle_stable(a, a, (5.0, 6.0, -1.0, 1.0))


def test_no_stable_pair_returns_empty_not_an_error():
    a = np.array([1.0 + 0.0j, 2.0 + 0.0j])
    b = np.array([1.5 + 0.0j, 2.5 + 0.0j])
    energies, residuals, idx = match_angle_stable(a, b, WIDE)
    assert energies.size == 0 and residuals.size == 0 and idx.size == 0


def test_find_resonance_pole_behaviour_is_unchanged():
    # The single-pole finder returns the globally closest pair WHATEVER its
    # residual -- it must NOT inherit match_angle_stable's tolerance cut.
    a = np.array([1.0 + 0.0j, 2.0 + 0.0j])
    b = np.array([1.5 + 0.0j, 2.5 + 0.0j])
    pole, residual = find_resonance_pole(a, b, WIDE)
    assert np.isclose(pole.real, 1.25) and np.isclose(residual, 0.5)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run --no-sync pytest libs/qscat/tests/test_ecs_pole_match.py -v`
Expected: FAIL — `ImportError: cannot import name 'match_angle_stable' from 'qscat.ecs'`.

- [ ] **Step 3: Implement `match_angle_stable`**

In `libs/qscat/qscat/ecs/pole.py`, replace the module's `__all__` and the body of
`find_resonance_pole` from the `fa = _filter_window(...)` line onward with the
following. `_filter_window` (already present, unchanged) keeps returning values;
the new `_window_indices` returns *indices* so callers can map back to `eigs_a`.

```python
__all__ = ["find_resonance_pole", "match_angle_stable"]


def _window_indices(
    E: npt.NDArray[np.complex128], window: tuple[float, float, float, float]
) -> npt.NDArray[np.intp]:
    re_lo, re_hi, im_lo, im_hi = window
    mask = (E.real >= re_lo) & (E.real <= re_hi) & (E.imag >= im_lo) & (E.imag <= im_hi)
    return np.flatnonzero(mask)


def _paired(
    eigs_a: npt.ArrayLike,
    eigs_b: npt.ArrayLike,
    window: tuple[float, float, float, float],
    caller: str,
) -> tuple[
    npt.NDArray[np.complex128],
    npt.NDArray[np.complex128],
    npt.NDArray[np.intp],
    npt.NDArray[np.intp],
    npt.NDArray[np.float64],
]:
    """Window-filter both spectra and pair each surviving `a` with its nearest `b`.

    Returns `(fa, fb, ia, nearest, dist)`: the windowed values, `ia` their
    indices into the ORIGINAL `eigs_a`, `nearest[k]` the index into `fb` closest
    to `fa[k]`, and `dist[k]` that distance. Shared by `find_resonance_pole`
    (which takes the global argmin) and `match_angle_stable` (which applies a
    tolerance cut) so there is exactly one implementation of the criterion.
    """
    a = np.asarray(eigs_a, dtype=np.complex128)
    b = np.asarray(eigs_b, dtype=np.complex128)
    ia = _window_indices(a, window)
    ib = _window_indices(b, window)
    if ia.size == 0 or ib.size == 0:
        raise ValueError(
            f"{caller}: window {window} contains no eigenvalues in "
            f"{'eigs_a' if ia.size == 0 else 'eigs_b'} "
            f"(found {ia.size} in A, {ib.size} in B) -- window too tight "
            "or spectrum too coarse."
        )
    fa, fb = a[ia], b[ib]
    diffs = np.abs(fa[:, None] - fb[None, :])
    nearest = np.asarray(np.argmin(diffs, axis=1), dtype=np.intp)
    dist = np.asarray(diffs[np.arange(fa.size), nearest], dtype=np.float64)
    return fa, fb, ia, nearest, dist


def match_angle_stable(
    eigs_a: npt.ArrayLike,
    eigs_b: npt.ArrayLike,
    window: tuple[float, float, float, float],
    *,
    rel_tol: float = 1e-4,
    atol: float = 1e-8,
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64], npt.NDArray[np.intp]]:
    """Every angle-stable eigenvalue shared by two ECS spectra, not just one.

    The multi-state generalization of `find_resonance_pole`: `eigs_a`/`eigs_b`
    are complex eigenvalue arrays of the same Hamiltonian at two different ECS
    rotation angles. An eigenvalue of `eigs_a` inside `window` is ACCEPTED when
    its nearest `eigs_b` partner satisfies

        |E_a - E_b| < max(rel_tol * |E_a|, atol)

    -- eMoScat's `DiscreteStates` criterion, vectorized. Discretized continuum
    eigenvalues rotate with the angle and fail it; bound and resonance states do
    not. Returns `(energies, residuals, indices)`, ascending in `Re E`:
    `energies` the midpoints `(E_a + E_b)/2` (matching `find_resonance_pole`'s
    convention), `residuals` the `|E_a - E_b|` per accepted state, and `indices`
    the positions in the ORIGINAL `eigs_a` -- so a caller holding the grid-`a`
    eigenvectors can pull the matching columns straight out.

    An empty result is a normal outcome (no stable state in `window`), NOT an
    error. `ValueError` is raised only when `window` catches nothing at all in
    one of the two spectra, mirroring `find_resonance_pole`.
    """
    fa, fb, ia, nearest, dist = _paired(eigs_a, eigs_b, window, "match_angle_stable")
    keep = dist < np.maximum(rel_tol * np.abs(fa), atol)
    energies = 0.5 * (fa[keep] + fb[nearest[keep]])
    residuals = dist[keep]
    indices = ia[keep]
    order = np.argsort(energies.real)
    return (
        np.asarray(energies[order], dtype=np.complex128),
        np.asarray(residuals[order], dtype=np.float64),
        np.asarray(indices[order], dtype=np.intp),
    )
```

Then rewrite `find_resonance_pole`'s body (keep its existing docstring verbatim)
so it reuses `_paired`:

```python
    fa, fb, _ia, nearest, dist = _paired(eigs_a, eigs_b, window, "find_resonance_pole")
    i = int(np.argmin(dist))
    ea, eb = fa[i], fb[nearest[i]]
    return complex(0.5 * (ea + eb)), float(np.abs(ea - eb))
```

Delete the now-unused `_filter_window` if nothing else references it (check with
`grep -rn "_filter_window" libs/`).

Then in `libs/qscat/qscat/ecs/__init__.py`:

```python
from .pole import find_resonance_pole, match_angle_stable

__all__ = ["ecs_map", "find_resonance_pole", "match_angle_stable"]
```

and add to that module's docstring, after the existing `find_resonance_pole`
paragraph:

```
`match_angle_stable` is its multi-state sibling: same acceptance criterion,
but it returns EVERY angle-stable eigenvalue in a window (with the indices
needed to recover the eigenvectors), which is what a level spectrum needs.
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run --no-sync pytest libs/qscat/tests/test_ecs_pole_match.py libs/qscat/tests/test_pole.py -v`
Expected: PASS — all new tests, and every pre-existing `find_resonance_pole` test
still green (its behaviour must be byte-identical).

- [ ] **Step 5: Run the wider suite for regressions**

Run: `uv run --no-sync pytest libs/qscat/tests -q -k "pole or lcp or resonance"`
Expected: PASS. `local_complex_potential` calls `find_resonance_pole` on every
walk step, so this is the real regression signal.

- [ ] **Step 6: Commit**

```bash
uv run ruff format libs/qscat/qscat/ecs libs/qscat/tests/test_ecs_pole_match.py
uv run ruff check libs/qscat/qscat/ecs libs/qscat/tests/test_ecs_pole_match.py
git add libs/qscat/qscat/ecs/pole.py libs/qscat/qscat/ecs/__init__.py libs/qscat/tests/test_ecs_pole_match.py
git commit -m "feat(ecs): match_angle_stable -- multi-state two-angle pole matcher"
```

---

### Task 2: Extract `_assemble_lcp` from `local_complex_potential`

**Files:**
- Modify: `libs/qscat/qscat/core/lcp.py:267-276` (the assembly tail of `local_complex_potential`)

**Interfaces:**
- Consumes: nothing.
- Produces: `_assemble_lcp(model, grid, shift, gamma_w) -> tuple[NDArray[complex128], NDArray[float64]]`
  — private to `lcp.py`, used by Task 4's `resonance_levels` to place ONE walk's
  result onto TWO nuclear grids.

This is a **pure refactor with zero behaviour change.** The existing
`libs/qscat/tests/test_lcp.py` is the gate: it must pass untouched.

- [ ] **Step 1: Confirm the existing tests pass before touching anything**

Run: `uv run --no-sync pytest libs/qscat/tests/test_lcp.py -q`
Expected: PASS. Record the count — it must be identical after the refactor.

- [ ] **Step 2: Add the helper**

In `libs/qscat/qscat/core/lcp.py`, immediately above `local_complex_potential`:

```python
def _assemble_lcp(
    model: ResonanceModel,
    grid: FemDvrEcsGrid,
    shift: npt.NDArray[np.float64],
    gamma_w: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64]]:
    """Place one `resonance_pole_walk` result onto `grid` as `(V_d, Gamma)`.

    `shift`/`gamma_w` are indexed by DESCENDING real `R` (the walk order).
    Real nodes get `V_d = v0(R) + shift`; the complex ECS tail gets the
    analytic continuation `v0(z) + shift[0]` (the shift at the largest real `R`,
    i.e. the asymptotic electronic shift) with `Gamma = 0`.

    Factored out of `local_complex_potential` so `resonance_levels` can run the
    expensive electronic walk ONCE and lay the same curve onto two nuclear grids
    that differ only in their ECS tail angle.
    """
    pts = grid.points
    real_idx = np.flatnonzero(pts.imag == 0.0)
    order = np.argsort(pts[real_idx].real)[::-1]  # descending R: outer -> inner
    walk = real_idx[order]

    Vd = np.empty(grid.n, dtype=np.complex128)
    Gamma = np.zeros(grid.n, dtype=np.float64)
    Vd[walk] = model.v0(pts[walk].real) + shift
    Gamma[walk] = gamma_w

    tail = np.flatnonzero(pts.imag != 0.0)
    if tail.size:
        Vd[tail] = model.v0(pts[tail]) + shift[0]
    return Vd, Gamma
```

- [ ] **Step 3: Rewrite `local_complex_potential`'s tail to call it**

Replace everything in `local_complex_potential` from `Vd = np.empty(...)` to the
final `return Vd, Gamma` with:

```python
    return _assemble_lcp(model, nuclear_grid, shift, gamma_w)
```

The lines above it (seed, `real_idx`/`order`/`walk`/`R_real`, the
`resonance_pole_walk` call) stay exactly as they are — `R_real` is still needed
as the walk's input.

- [ ] **Step 4: Run the tests to verify nothing changed**

Run: `uv run --no-sync pytest libs/qscat/tests/test_lcp.py -q`
Expected: PASS with the SAME count as Step 1. Any difference means the refactor
changed behaviour — revert and redo.

- [ ] **Step 5: Commit**

```bash
uv run ruff format libs/qscat/qscat/core/lcp.py
uv run ruff check libs/qscat/qscat/core/lcp.py
git add libs/qscat/qscat/core/lcp.py
git commit -m "refactor(lcp): extract _assemble_lcp (no behaviour change)"
```

---

### Task 3: `ResonanceLevels` + `lcp_resonance_levels` (the numeric core)

**Files:**
- Modify: `libs/qscat/qscat/core/lcp.py`
- Test: `libs/qscat/tests/test_lcp_resonance_levels.py` (create)

**Interfaces:**
- Consumes: `qscat.ecs.match_angle_stable` (Task 1).
- Produces:
  - `ResonanceLevels` frozen dataclass with fields `energies`, `widths`,
    `states`, `residuals`, `real_weight`, `golden_rule`.
  - `lcp_resonance_levels(nuclear_grid_a, nuclear_grid_b, mu, Vd_a, Vd_b, Gamma, *,
    window=None, n_levels=None, rel_tol=1e-4, atol=1e-8, golden_rule=True) -> ResonanceLevels`.
  Task 4 wraps it; Task 7 consumes the dataclass.

**Physics the implementer needs.** The nuclear Hamiltonian is
`H_N = T(mu) + diag(W)`, `W(R) = V_d(R) − iΓ(R)/2` — complex **symmetric** (not
Hermitian), which is why eigenvectors must be normalized with the bilinear
c-product `Σ c_i² = 1` and never with `‖c‖₂`. `qscat.dvr.kinetic(grid, mass)`
builds `T`; `qscat.dvr.eigen(H)` returns `(E, V)` already sorted by ascending
`Re E`, with numpy's `v†v = 1` normalization that we then replace.

The **golden-rule comparator** solves the same problem with `Γ` set to zero (NOT
with `Re(V_d)`: taking a real part would destroy the analytic continuation in the
ECS tail, where `V_d = v0(z)` is legitimately complex). Its first-order width is
`Γ_v = Σ_i c_i² Γ_i` — in a DVR the matrix element of a diagonal operator is just
that sum. Comparator levels are paired to complex levels by nearest `Re E`.

- [ ] **Step 1: Write the failing tests**

Create `libs/qscat/tests/test_lcp_resonance_levels.py`:

```python
"""Analytic oracles for `qscat.core.lcp.lcp_resonance_levels`.

Two exact benchmarks, no convergence hand-waving:

1. Gamma = 0 with a bare Morse curve -> the analytic Morse spectrum,
   E_n = -D (1 - alpha (n + 1/2) / sqrt(2 mu D))^2, with Im E ~ 0.
2. Gamma = Gamma_0 CONSTANT -> a constant imaginary term commutes with
   everything, so the spectrum must shift rigidly by exactly -i Gamma_0/2.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import segmented_grid
from qscat.core.lcp import lcp_resonance_levels

# A deep, heavy Morse well (N2-like) so the analytic infinite-domain formula
# holds to high accuracy: V(0) ~ 64 Ha is an effectively infinite inner wall.
MU, D0, ALPHA0, RE = 12766.36, 0.75102, 1.1535, 2.01943

# Lowest five analytic levels; the sixth (~ -0.6838) lies outside WINDOW.
WINDOW = (-0.76, -0.69, -1e-6, 1e-6)


def morse_levels(n: np.ndarray) -> np.ndarray:
    x = ALPHA0 * (n + 0.5) / np.sqrt(2.0 * MU * D0)
    return -D0 * (1.0 - x) ** 2


def morse(R: np.ndarray) -> np.ndarray:
    z = np.asarray(R, dtype=np.complex128)
    return D0 * (np.exp(-2 * ALPHA0 * (z - RE)) - 2 * np.exp(-ALPHA0 * (z - RE)))


def grid_pair(angle_a: float = 35.0, angle_b: float = 25.0, n_real: int = 30):
    """Two nuclear grids sharing every real node, differing only in tail angle."""
    real, complex_ = [(n_real, 6.0)], [(6, 16.0)]
    return (
        segmented_grid(real, complex_, angle_deg=angle_a, quadrature=10),
        segmented_grid(real, complex_, angle_deg=angle_b, quadrature=10),
    )


def test_zero_width_reproduces_the_analytic_morse_spectrum():
    ga, gb = grid_pair()
    Vd_a, Vd_b = morse(ga.points), morse(gb.points)
    Gamma = np.zeros(ga.n, dtype=np.float64)

    out = lcp_resonance_levels(ga, gb, MU, Vd_a, Vd_b, Gamma, window=WINDOW)

    assert out.energies.shape == (5,)
    np.testing.assert_allclose(
        out.energies.real, morse_levels(np.arange(5)), rtol=1e-5
    )
    assert np.all(np.abs(out.energies.imag) < 1e-8)
    assert np.all(out.widths < 1e-8)
    # Bound states live entirely in the real region.
    assert np.all(out.real_weight > 0.999)


def test_constant_width_shifts_the_spectrum_rigidly():
    ga, gb = grid_pair()
    Vd_a, Vd_b = morse(ga.points), morse(gb.points)
    g0 = 0.01
    zero = lcp_resonance_levels(
        ga, gb, MU, Vd_a, Vd_b, np.zeros(ga.n), window=WINDOW
    )
    shifted = lcp_resonance_levels(
        ga,
        gb,
        MU,
        Vd_a,
        Vd_b,
        np.full(ga.n, g0),
        window=(-0.76, -0.69, -0.5 * g0 - 1e-6, -0.5 * g0 + 1e-6),
    )
    # H(Gamma_0) = H(0) - i (Gamma_0/2) I exactly -- round-off, not tolerance.
    np.testing.assert_allclose(
        shifted.energies, zero.energies - 0.5j * g0, atol=1e-10
    )
    np.testing.assert_allclose(shifted.widths, g0, atol=1e-10)


def test_states_are_c_product_normalized():
    from qscat.linalg import c_product

    ga, gb = grid_pair()
    out = lcp_resonance_levels(
        ga, gb, MU, morse(ga.points), morse(gb.points), np.zeros(ga.n), window=WINDOW
    )
    for state in out.states:
        assert abs(c_product(state, state) - 1.0) < 1e-10


def test_n_levels_truncates_to_the_lowest():
    ga, gb = grid_pair()
    out = lcp_resonance_levels(
        ga,
        gb,
        MU,
        morse(ga.points),
        morse(gb.points),
        np.zeros(ga.n),
        window=WINDOW,
        n_levels=2,
    )
    assert out.energies.shape == (2,)
    np.testing.assert_allclose(out.energies.real, morse_levels(np.arange(2)), rtol=1e-5)


def test_mismatched_real_regions_raise():
    ga, _ = grid_pair()
    _, gb_wrong = grid_pair(n_real=20)  # different real discretization
    with pytest.raises(ValueError, match="real nodes"):
        lcp_resonance_levels(
            ga, gb_wrong, MU, morse(ga.points), morse(gb_wrong.points),
            np.zeros(ga.n), window=WINDOW,
        )


def test_non_positive_mass_raises():
    ga, gb = grid_pair()
    with pytest.raises(ValueError, match="mu must be positive"):
        lcp_resonance_levels(
            ga, gb, 0.0, morse(ga.points), morse(gb.points), np.zeros(ga.n),
            window=WINDOW,
        )
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run --no-sync pytest libs/qscat/tests/test_lcp_resonance_levels.py -v`
Expected: FAIL — `ImportError: cannot import name 'lcp_resonance_levels'`.

- [ ] **Step 3: Implement the dataclass and the numeric core**

In `libs/qscat/qscat/core/lcp.py`, add to the imports:

```python
import warnings
from dataclasses import dataclass

from qscat.ecs import find_resonance_pole, match_angle_stable
```

(the `find_resonance_pole` import already exists — extend that line).

Extend `__all__` with `"ResonanceLevels"`, `"lcp_resonance_levels"`.

Append to the module:

```python
@dataclass(frozen=True)
class ResonanceLevels:
    """The quasi-bound vibrational levels of the anion in the LCP curve.

    The thesis's `omega_j` (Vana 2017, Sec. 1.5/3.4), promoted from real levels
    of `Re V_res` to genuine complex eigenvalues. These are complex-scaled (ECS)
    resonance eigenstates -- NOT Siegert pseudostates, which carry an
    outgoing-wave condition at a finite radius and a surface-corrected
    orthogonality relation (Hvizdos et al., Phys. Rev. A 97, 022704 (2018),
    App. A). ECS rotates rather than truncates, so the plain bilinear c-product
    is the complete inner product here.

    - `energies`: `E_v - i Gamma_v/2` (Hartree), ascending in `Re E`.
    - `widths`: `Gamma_v = max(0, -2 Im E_v)` (Hartree). A level below the anion
      dissociation limit carries only the ELECTRONIC autodetachment width; one
      above it also carries a NUCLEAR (dissociative) width. Both come out of the
      one diagonalization.
    - `states`: shape `(n_levels, grid.n)` DVR COEFFICIENTS `c_i`
      (`psi(R_i) = c_i / sqrt(w_i)`), c-product-normalized: `sum_i c_i^2 = 1`.
    - `residuals`: the two-angle stability residual per level. A level whose
      residual does not fall under grid refinement is not converged.
    - `real_weight`: fraction of `|c|^2` inside the real region -- a diagnostic,
      not a normalization. Near 1 for a well-localized level.
    - `golden_rule`: `E_v^(0) - i <chi_v|Gamma|chi_v>/2`, the perturbative
      comparator (the `Gamma = 0` levels plus the first-order width). This is
      what eMoScat and the thesis actually computed. Agreement with `energies`
      means the level is perturbative; divergence means it is genuinely broad
      and the non-perturbative treatment is load-bearing. `nan` where no
      comparator level could be paired, and all-`nan` when `golden_rule=False`.
    """

    energies: npt.NDArray[np.complex128]
    widths: npt.NDArray[np.float64]
    states: npt.NDArray[np.complex128]
    residuals: npt.NDArray[np.float64]
    real_weight: npt.NDArray[np.float64]
    golden_rule: npt.NDArray[np.complex128]


def _default_window(
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    grid: FemDvrEcsGrid,
    atol: float,
) -> tuple[float, float, float, float]:
    """`Re` spanning the anion curve over the REAL nodes, `Im` down to `-max Gamma`.

    The real span is taken from the real nodes only (the ECS tail's continued
    `v0(z)` is complex and says nothing about where levels lie), and covers the
    whole curve so neither the well bottom nor levels above the neutral
    dissociation limit `v0(inf) = 0` are cut.
    """
    real = grid.points.imag == 0.0
    v = Vd[real].real
    return (float(v.min()), float(v.max()), -float(max(Gamma.max(), atol)), atol)


def _levels_from(
    grid_a: FemDvrEcsGrid,
    grid_b: FemDvrEcsGrid,
    mu: float,
    W_a: npt.NDArray[np.complex128],
    W_b: npt.NDArray[np.complex128],
    window: tuple[float, float, float, float],
    rel_tol: float,
    atol: float,
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64], npt.NDArray[np.complex128]]:
    """Diagonalize `T(mu) + diag(W)` on both grids, keep the angle-stable states.

    Returns `(energies, residuals, states)` with `states` shape
    `(n_levels, grid_a.n)`, c-product-normalized, taken from grid `a`.
    """
    E_a, V_a = eigen(kinetic(grid_a, mu) + np.diag(W_a))
    E_b, _ = eigen(kinetic(grid_b, mu) + np.diag(W_b))
    energies, residuals, idx = match_angle_stable(
        E_a, E_b, window, rel_tol=rel_tol, atol=atol
    )
    states = np.empty((idx.size, grid_a.n), dtype=np.complex128)
    for k, j in enumerate(idx):
        c = V_a[:, j].astype(np.complex128)
        states[k] = c / np.sqrt(c_product(c, c))
    return energies, residuals, states


def lcp_resonance_levels(
    nuclear_grid_a: FemDvrEcsGrid,
    nuclear_grid_b: FemDvrEcsGrid,
    mu: float,
    Vd_a: npt.NDArray[np.complex128],
    Vd_b: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    *,
    window: tuple[float, float, float, float] | None = None,
    n_levels: int | None = None,
    rel_tol: float = 1e-4,
    atol: float = 1e-8,
    golden_rule: bool = True,
) -> ResonanceLevels:
    """Quasi-bound levels of `H_N = T(mu) + V_d(R) - i Gamma(R)/2`.

    The Born-Oppenheimer approximation to the 2-D model's resonance energies:
    step 1 (the fixed-`R` electronic pole, `local_complex_potential`) supplies
    the complex curve; this is step 2, the nuclear eigenvalue problem in it. The
    thesis's `H_LCP` (Vana 2017 Eq. 1.65).

    `nuclear_grid_a`/`nuclear_grid_b` must share every real node and differ only
    in their ECS tail angle -- that is what makes the two spectra comparable.
    `Vd_a`/`Vd_b` are the curve laid onto each grid (identical on the real
    nodes, differing in the continued tail); `Gamma` is real and tail-zero, so
    the same array serves both.

    Physical levels are selected by two-angle stability (`match_angle_stable`);
    the rotated dissociative continuum fails that test and drops out. Levels
    with `Im E > atol` are unphysical and are dropped with a warning.

    See `docs/physics/lcp-resonance-levels.md`.
    """
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu}")
    for name, arr, grid in (
        ("Vd_a", Vd_a, nuclear_grid_a),
        ("Vd_b", Vd_b, nuclear_grid_b),
        ("Gamma", Gamma, nuclear_grid_a),
    ):
        if arr.shape != (grid.n,):
            raise ValueError(f"{name} has shape {arr.shape}, expected ({grid.n},)")
    ra, rb = nuclear_grid_a.points, nuclear_grid_b.points
    if ra.size != rb.size or not np.array_equal(ra[ra.imag == 0.0], rb[rb.imag == 0.0]):
        raise ValueError(
            "nuclear_grid_a and nuclear_grid_b must share their real nodes "
            "(same real segments and quadrature; only the ECS tail angle may "
            "differ) -- otherwise the two spectra are not comparable"
        )

    if window is None:
        window = _default_window(Vd_a, Gamma, nuclear_grid_a, atol)

    half_i_gamma = 0.5j * Gamma
    energies, residuals, states = _levels_from(
        nuclear_grid_a, nuclear_grid_b, mu,
        Vd_a - half_i_gamma, Vd_b - half_i_gamma,
        window, rel_tol, atol,
    )

    physical = energies.imag <= atol
    if not physical.all():
        warnings.warn(
            f"lcp_resonance_levels: dropped {int((~physical).sum())} level(s) with "
            "Im E > 0 (unphysical: a growing state). Usually an over-wide window "
            "or an under-resolved grid.",
            UserWarning,
            stacklevel=2,
        )
    energies, residuals, states = energies[physical], residuals[physical], states[physical]

    if n_levels is not None:
        energies, residuals, states = (
            energies[:n_levels], residuals[:n_levels], states[:n_levels]
        )

    widths = np.maximum(0.0, -2.0 * energies.imag)
    real_mask = nuclear_grid_a.points.imag == 0.0
    dens = np.abs(states) ** 2
    total = dens.sum(axis=1)
    real_weight = np.divide(
        dens[:, real_mask].sum(axis=1), total,
        out=np.zeros_like(total), where=total > 0.0,
    )

    gr = np.full(energies.size, np.nan + 1j * np.nan, dtype=np.complex128)
    if golden_rule and energies.size:
        E0, _resid0, chi0 = _levels_from(
            nuclear_grid_a, nuclear_grid_b, mu, Vd_a, Vd_b,
            (window[0], window[1], -atol, atol), rel_tol, atol,
        )
        if E0.size:
            g1 = np.array([c_product(c, Gamma * c).real for c in chi0])
            # Pair each complex level to the nearest comparator level in Re E.
            near = np.argmin(np.abs(energies.real[:, None] - E0.real[None, :]), axis=1)
            gr = E0[near].real - 0.5j * g1[near]

    return ResonanceLevels(
        energies=np.asarray(energies, dtype=np.complex128),
        widths=np.asarray(widths, dtype=np.float64),
        states=np.asarray(states, dtype=np.complex128),
        residuals=np.asarray(residuals, dtype=np.float64),
        real_weight=np.asarray(real_weight, dtype=np.float64),
        golden_rule=np.asarray(gr, dtype=np.complex128),
    )
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run --no-sync pytest libs/qscat/tests/test_lcp_resonance_levels.py -v`
Expected: PASS, all 6 tests.

If `test_zero_width_reproduces_the_analytic_morse_spectrum` fails on tolerance
rather than on shape, the grid is under-resolved — raise `n_real` from 30 to 50
in `grid_pair` and re-run. Do **not** loosen `rtol` past `1e-4`; if 1e-4 still
fails, something is wrong with the kinetic term or the mass, not the grid.

- [ ] **Step 5: Verify the core-purity guard still holds**

Run: `uv run --no-sync pytest libs/qscat/tests/test_core_no_model_import.py -q`
Expected: PASS. `lcp_resonance_levels` takes `mu` as a float and touches no
model — only Task 4's wrapper mentions `ResonanceModel`, under `TYPE_CHECKING`.

- [ ] **Step 6: Commit**

```bash
uv run ruff format libs/qscat/qscat/core/lcp.py libs/qscat/tests/test_lcp_resonance_levels.py
uv run ruff check libs/qscat/qscat/core/lcp.py libs/qscat/tests/test_lcp_resonance_levels.py
uv run mypy libs/qscat
git add libs/qscat/qscat/core/lcp.py libs/qscat/tests/test_lcp_resonance_levels.py
git commit -m "feat(lcp): lcp_resonance_levels -- complex quasi-bound nuclear levels"
```

---

### Task 4: `resonance_levels` — the model-level convenience wrapper

**Files:**
- Modify: `libs/qscat/qscat/core/lcp.py`
- Modify: `libs/qscat/qscat/model/ionic.py`
- Test: `libs/qscat/tests/test_lcp_resonance_levels.py` (extend)

**Interfaces:**
- Consumes: `_assemble_lcp` (Task 2), `lcp_resonance_levels` (Task 3),
  `resonance_pole_walk` + `anion_electronic_states` (existing).
- Produces: `resonance_levels(model, nuclear_grid_a, nuclear_grid_b, elec_grid_a,
  elec_grid_b, *, re_half_width=0.05, im_half_width=0.05, resid_tol=1e-3,
  **levels_kwargs) -> ResonanceLevels`. Task 7 calls this.
- Produces: `IonicResonanceModel.max_nuclear_ecs_angle_deg: float = 22.5`.

**Why the angle bound exists.** Hvizdoš et al. §II: for the H₂⁺ model the nuclear
ECS angle must satisfy `theta < pi/8` (22.5°) and the electronic `theta < pi/4`,
or `V(R,r)` **diverges** at large `R`/`r` — the quartic `a3 * R**4` term in the
coupling needs `4*theta < pi/2`. Neutral diatomics have no such bound (their
Morse and Gaussian forms are entire), so the attribute is optional and read with
`getattr`.

- [ ] **Step 1: Write the failing tests**

Append to `libs/qscat/tests/test_lcp_resonance_levels.py`:

```python
def _lcp_grids():
    """F2's LCP decks: fine nuclear grid at two tail angles + two electronic angles."""
    from qscat.core.grids import electronic_grid

    real = [(9, 1.8), (1, 2.0), (5, 2.5), (4, 2.596908), (4, 2.7), (4, 10.7)]
    cx = [(15, 30.0)]
    nuc_a = segmented_grid(real, cx, angle_deg=25.0, quadrature=12)
    nuc_b = segmented_grid(real, cx, angle_deg=15.0, quadrature=12)
    ea = electronic_grid(r_max=16.0, order=8, n_complex=6, angle_deg=35.0)
    eb = electronic_grid(r_max=16.0, order=8, n_complex=6, angle_deg=44.0)
    return nuc_a, nuc_b, ea, eb


@pytest.mark.slow
def test_f2_levels_are_bound_and_narrow():
    from qscat.core.lcp import resonance_levels
    from qscat.model import F2

    nuc_a, nuc_b, ea, eb = _lcp_grids()
    out = resonance_levels(F2, nuc_a, nuc_b, ea, eb, n_levels=6)

    assert out.energies.size > 0
    assert np.all(np.diff(out.energies.real) > 0)      # ascending, non-degenerate
    assert np.all(out.widths >= 0.0)                   # clamped
    assert np.all(out.residuals < 1e-3)                # genuinely angle-stable
    assert np.all(out.real_weight > 0.5)               # localized, not continuum
    # The comparator must agree with the complex result on the narrow levels:
    # for Gamma_v << level spacing the shift is first-order.
    narrow = out.widths < 1e-4
    assert narrow.any()
    np.testing.assert_allclose(
        out.energies.real[narrow], out.golden_rule.real[narrow], atol=1e-4
    )


@pytest.mark.slow
def test_gamma_support_condition_holds_for_f2():
    """Vana 2017 Sec. 1.5: Im V_res is nonzero ONLY where v0(R) < E_res(R)."""
    from qscat.core.lcp import local_complex_potential
    from qscat.model import F2

    nuc_a, _nuc_b, ea, eb = _lcp_grids()
    Vd, Gamma = local_complex_potential(F2, nuc_a, ea, eb)
    real = nuc_a.points.imag == 0.0
    R = nuc_a.points[real].real
    bound_region = Vd[real].real < F2.v0(R).real  # anion below neutral: no autodetachment
    assert np.all(Gamma[real][bound_region] < 1e-6)


def test_grid_angle_bound_is_enforced():
    from qscat.core.grids import electronic_grid
    from qscat.core.lcp import resonance_levels
    from qscat.model import H2P

    real, cx = [(5, 1.0), (90, 14.05)], [(3, 30.0)]
    too_steep = segmented_grid(real, cx, angle_deg=30.0, quadrature=8)  # > 22.5
    ok = segmented_grid(real, cx, angle_deg=20.0, quadrature=8)
    e = electronic_grid(r_max=16.0, order=8, n_complex=6, angle_deg=35.0)
    with pytest.raises(ValueError, match="max_nuclear_ecs_angle_deg"):
        resonance_levels(H2P, too_steep, ok, e, e)
```

Register the `slow` marker if it is not already in `pyproject.toml` (check with
`grep -n "markers" pyproject.toml`; the repo already uses `@pytest.mark.slow`
elsewhere, so it should be present).

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run --no-sync pytest libs/qscat/tests/test_lcp_resonance_levels.py -v -k "f2 or angle_bound or support"`
Expected: FAIL — `ImportError: cannot import name 'resonance_levels'`.

- [ ] **Step 3: Add the model attribute**

In `libs/qscat/qscat/model/ionic.py`, add to the `IonicResonanceModel` dataclass
body (after the existing fields):

```python
    # Hvizdos et al., Phys. Rev. A 97, 022704 (2018), Sec. II: the nuclear ECS
    # angle must stay below pi/8 or the quartic `a3 * R**4` term in `v_int`
    # diverges under the rotation (4*theta < pi/2). The electronic bound is
    # pi/4 (from exp(-r^2/3), 2*theta < pi/2). Neutral diatomics have no such
    # bound -- their Morse + Gaussian forms are entire -- so this attribute is
    # optional on the protocol and read with getattr.
    max_nuclear_ecs_angle_deg: float = 22.5
```

- [ ] **Step 4: Implement `resonance_levels`**

Append to `libs/qscat/qscat/core/lcp.py` (and add `"resonance_levels"` to
`__all__`):

```python
def _check_angle_bound(model: ResonanceModel, *grids: FemDvrEcsGrid) -> None:
    """Reject nuclear grids whose ECS tail angle exceeds the model's bound."""
    bound = getattr(model, "max_nuclear_ecs_angle_deg", None)
    if bound is None:
        return
    for g in grids:
        worst = max((el.angle_deg for el in g.spec.elements), default=0.0)
        if worst > bound:
            raise ValueError(
                f"nuclear grid ECS angle {worst} deg exceeds this model's "
                f"max_nuclear_ecs_angle_deg = {bound} deg; beyond it the "
                "interaction potential diverges under the rotation "
                "(Hvizdos et al., Phys. Rev. A 97, 022704 (2018), Sec. II)"
            )


def resonance_levels(
    model: ResonanceModel,
    nuclear_grid_a: FemDvrEcsGrid,
    nuclear_grid_b: FemDvrEcsGrid,
    elec_grid_a: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    *,
    re_half_width: float = 0.05,
    im_half_width: float = 0.05,
    resid_tol: float = 1e-3,
    window: tuple[float, float, float, float] | None = None,
    n_levels: int | None = None,
    rel_tol: float = 1e-4,
    atol: float = 1e-8,
    golden_rule: bool = True,
) -> ResonanceLevels:
    """Quasi-bound levels of `model`'s anion, straight from the model.

    Runs the electronic pole walk ONCE (`resonance_pole_walk`, seeded from the
    asymptotic anion bound state exactly as `local_complex_potential` does),
    lays the resulting curve onto BOTH nuclear grids with `_assemble_lcp`, and
    diagonalizes (`lcp_resonance_levels`). `E_res(R)` at real `R` does not
    depend on the nuclear tail angle, so the second grid costs one extra nuclear
    diagonalization and nothing else.

    `nuclear_grid_b` must share `nuclear_grid_a`'s real segments and quadrature
    and differ only in its ECS tail angle -- conventionally a SMALLER angle,
    which is always safe against the model's divergence bound.
    """
    _check_angle_bound(model, nuclear_grid_a, nuclear_grid_b)

    R_inf = nuclear_grid_a.R0
    eps_e, _ = anion_electronic_states(elec_grid_a, model, R_inf, 1)
    seed_window = (
        eps_e[0] - re_half_width, eps_e[0] + re_half_width, -im_half_width, im_half_width,
    )

    pts = nuclear_grid_a.points
    real_idx = np.flatnonzero(pts.imag == 0.0)
    walk = real_idx[np.argsort(pts[real_idx].real)[::-1]]  # descending R
    shift, gamma_w = resonance_pole_walk(
        model, pts[walk].real, elec_grid_a, elec_grid_b, seed_window,
        re_half_width=re_half_width, im_half_width=im_half_width, resid_tol=resid_tol,
    )

    Vd_a, Gamma = _assemble_lcp(model, nuclear_grid_a, shift, gamma_w)
    Vd_b, _ = _assemble_lcp(model, nuclear_grid_b, shift, gamma_w)

    real = pts.imag == 0.0
    bound_region = Vd_a[real].real < np.asarray(model.v0(pts[real].real)).real
    if np.any(Gamma[real][bound_region] > 1e-6):
        warnings.warn(
            "resonance_levels: Gamma(R) is nonzero where the anion curve lies "
            "BELOW the neutral (v0 > E_res), where autodetachment is closed "
            "(Vana 2017, Sec. 1.5). The widths downstream are suspect.",
            UserWarning,
            stacklevel=2,
        )

    return lcp_resonance_levels(
        nuclear_grid_a, nuclear_grid_b, model.mu, Vd_a, Vd_b, Gamma,
        window=window, n_levels=n_levels, rel_tol=rel_tol, atol=atol,
        golden_rule=golden_rule,
    )
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `uv run --no-sync pytest libs/qscat/tests/test_lcp_resonance_levels.py -v --runslow`
(if the repo has no `--runslow` flag, use `-m "slow or not slow"`; check
`grep -n "runslow\|addopts" pyproject.toml conftest.py 2>/dev/null`).
Expected: PASS. The F2 test takes ~1-3 minutes (the electronic walk dominates).

- [ ] **Step 6: Run the full library suite**

Run: `uv run --no-sync pytest libs/qscat/tests -q`
Expected: PASS, no regressions.

- [ ] **Step 7: Commit**

```bash
uv run ruff format libs/qscat/qscat libs/qscat/tests
uv run ruff check libs/qscat
uv run mypy libs/qscat
git add libs/qscat/qscat/core/lcp.py libs/qscat/qscat/model/ionic.py libs/qscat/tests/test_lcp_resonance_levels.py
git commit -m "feat(lcp): resonance_levels model wrapper + ECS angle bound"
```

---

### Task 5: Convergence and golden-rule validation

**Files:**
- Test: `libs/qscat/tests/test_lcp_resonance_levels.py` (extend)

**Interfaces:**
- Consumes: `lcp_resonance_levels`, `resonance_levels` (Tasks 3-4).
- Produces: no new API. This task exists because spec validations 3 and 4
  ("angle and grid independence", "golden-rule consistency") are the checks that
  decide whether a reported level is trustworthy, and they need their own gate.

- [ ] **Step 1: Write the failing tests**

Append to `libs/qscat/tests/test_lcp_resonance_levels.py`:

```python
def test_levels_are_independent_of_the_ecs_angle_pair():
    """A physical level does not move when the rotation angles change."""
    ga1, gb1 = grid_pair(angle_a=35.0, angle_b=25.0)
    ga2, gb2 = grid_pair(angle_a=45.0, angle_b=20.0)
    kw = dict(window=WINDOW, n_levels=5)
    one = lcp_resonance_levels(
        ga1, gb1, MU, morse(ga1.points), morse(gb1.points), np.zeros(ga1.n), **kw
    )
    two = lcp_resonance_levels(
        ga2, gb2, MU, morse(ga2.points), morse(gb2.points), np.zeros(ga2.n), **kw
    )
    np.testing.assert_allclose(one.energies, two.energies, rtol=1e-6, atol=1e-9)


def test_levels_converge_under_h_refinement():
    """Refining the real region must not move the levels -- and must shrink
    the two-angle residual, which is the per-level convergence metric."""
    coarse_a, coarse_b = grid_pair(n_real=20)
    fine_a, fine_b = grid_pair(n_real=40)
    kw = dict(window=WINDOW, n_levels=5)
    coarse = lcp_resonance_levels(
        coarse_a, coarse_b, MU, morse(coarse_a.points), morse(coarse_b.points),
        np.zeros(coarse_a.n), **kw
    )
    fine = lcp_resonance_levels(
        fine_a, fine_b, MU, morse(fine_a.points), morse(fine_b.points),
        np.zeros(fine_a.n), **kw
    )
    np.testing.assert_allclose(coarse.energies.real, fine.energies.real, rtol=1e-6)
    assert fine.residuals.max() <= coarse.residuals.max() * 10.0


def test_golden_rule_matches_the_complex_result_for_a_weak_constant_width():
    """First-order perturbation theory is EXACT for a constant Gamma, so the
    comparator must reproduce the complex levels to round-off."""
    ga, gb = grid_pair()
    g0 = 1e-5
    out = lcp_resonance_levels(
        ga, gb, MU, morse(ga.points), morse(gb.points), np.full(ga.n, g0),
        window=(-0.76, -0.69, -0.5 * g0 - 1e-6, -0.5 * g0 + 1e-6),
    )
    np.testing.assert_allclose(out.golden_rule, out.energies, atol=1e-9)


def test_golden_rule_can_be_switched_off():
    ga, gb = grid_pair()
    out = lcp_resonance_levels(
        ga, gb, MU, morse(ga.points), morse(gb.points), np.zeros(ga.n),
        window=WINDOW, golden_rule=False,
    )
    assert np.all(np.isnan(out.golden_rule.real))
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run --no-sync pytest libs/qscat/tests/test_lcp_resonance_levels.py -v -k "independent or refinement or golden"`
Expected: the two `golden_rule` tests may already pass (Task 3 implemented the
field); the angle/refinement tests are the new signal. Any FAIL here is a real
defect in Task 3, not a missing feature — fix `lcp_resonance_levels`, do not
loosen the tolerance.

- [ ] **Step 3: Fix any failures in `lcp_resonance_levels`**

There is no new implementation in this task. If a test fails, the likely causes,
in order: (a) `_default_window` being used instead of the explicit `WINDOW` —
check the kwargs actually reach the call; (b) the c-normalization dividing by a
near-zero `c_product` for a continuum state that slipped through the tolerance —
tighten `rel_tol`, do not special-case; (c) the golden-rule nearest-`Re` pairing
picking the wrong partner when levels are dense — verify `E0.size == energies.size`
in a scratch run before assuming a deeper bug.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run --no-sync pytest libs/qscat/tests/test_lcp_resonance_levels.py -v`
Expected: PASS (all non-slow tests).

- [ ] **Step 5: Commit**

```bash
uv run ruff format libs/qscat/tests/test_lcp_resonance_levels.py
uv run ruff check libs/qscat/tests/test_lcp_resonance_levels.py
git add libs/qscat/tests/test_lcp_resonance_levels.py
git commit -m "test(lcp): angle/grid independence + golden-rule consistency gates"
```

---

### Task 6: Config — the `resonance_levels` observable kind and artifact flag

**Files:**
- Modify: `apps/qscat-run/qscat_run/config.py`
- Modify: `apps/qscat-run/qscat_run/presets.py`
- Test: `apps/qscat-run/tests/test_config.py` (extend)

**Interfaces:**
- Consumes: nothing from the library tasks (pure config).
- Produces:
  - `ArtifactSpec.resonance_levels: bool = False`
  - `GridSpec`-side `nuclear_angle_b: float | None = None` on the config's grid
    block (read as `cfg.grid.nuclear_angle_b`).
  - `"resonance_levels"` accepted in `Observable.kind` and in
    `presets.VALIDITY` for the LCP-capable molecules (`NO`, `F2`).
  - `presets.nuclear_grid_at_angle(cfg, angle_deg) -> FemDvrEcsGrid`.
  Task 7 consumes all four.

- [ ] **Step 1: Write the failing tests**

Append to `apps/qscat-run/tests/test_config.py`:

```python
def test_resonance_levels_observable_needs_no_energies(tmp_path):
    from qscat_run import presets
    from qscat_run.config import load_config, validate_config

    cfg_path = tmp_path / "levels.yaml"
    cfg_path.write_text(
        "molecule: F2\n"
        "methods: [lcp]\n"
        "observables:\n"
        "  - kind: resonance_levels\n"
        "    channels: 6\n"
    )
    cfg = load_config(cfg_path)
    validate_config(cfg)
    resolved = presets.resolve_defaults(cfg)
    assert resolved.observables[0].kind == "resonance_levels"
    assert resolved.energies is None  # levels-only run: no sweep to resolve


def test_resonance_levels_is_rejected_for_n2(tmp_path):
    from qscat_run.config import ConfigError, load_config, validate_config

    cfg_path = tmp_path / "bad.yaml"
    cfg_path.write_text(
        "molecule: N2\nmethods: [lcp]\nobservables: [{kind: resonance_levels}]\n"
    )
    with pytest.raises(ConfigError, match="not valid for N2"):
        validate_config(load_config(cfg_path))


def test_lcp_without_da_is_still_rejected_when_no_levels_requested(tmp_path):
    from qscat_run.config import ConfigError, load_config, validate_config

    cfg_path = tmp_path / "bad2.yaml"
    cfg_path.write_text(
        "molecule: F2\nmethods: [lcp]\nobservables: [{kind: ve, channels: 2}]\n"
        "energies: {min: 0.01, max: 0.05, step: 0.02}\n"
    )
    with pytest.raises(ConfigError, match="no 'da' or 'resonance_levels' observable"):
        validate_config(load_config(cfg_path))


def test_artifacts_resonance_levels_flag_parses(tmp_path):
    from qscat_run.config import load_config

    cfg_path = tmp_path / "art.yaml"
    cfg_path.write_text(
        "molecule: F2\nmethods: [lcp]\n"
        "observables: [{kind: da, channels: 1}]\n"
        "energies: {min: 0.01, max: 0.05, step: 0.02}\n"
        "artifacts: {resonance_levels: true}\n"
    )
    assert load_config(cfg_path).artifacts.resonance_levels is True


def test_nuclear_angle_b_defaults_to_ten_degrees_below(tmp_path):
    from qscat_run import presets
    from qscat_run.config import load_config

    cfg_path = tmp_path / "ang.yaml"
    cfg_path.write_text(
        "molecule: F2\nmethods: [lcp]\nobservables: [{kind: resonance_levels}]\n"
    )
    cfg = presets.resolve_defaults(load_config(cfg_path))
    g_a, _ea, _eb = presets.resolve_lcp_grids(cfg)
    g_b = presets.nuclear_grid_at_angle(cfg, presets.nuclear_angle_b(cfg))
    ang_a = max(el.angle_deg for el in g_a.spec.elements)
    ang_b = max(el.angle_deg for el in g_b.spec.elements)
    assert ang_b == pytest.approx(ang_a - 10.0)
    # The real nodes must be shared -- that is what makes the two spectra
    # comparable in `lcp_resonance_levels`.
    ra, rb = g_a.points, g_b.points
    assert np.array_equal(ra[ra.imag == 0.0], rb[rb.imag == 0.0])
```

Add `import numpy as np` and `import pytest` at the top of that test file if not
already present.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run --no-sync pytest apps/qscat-run/tests/test_config.py -v -k "resonance_levels or nuclear_angle_b"`
Expected: FAIL — the kind is rejected as unknown / attributes missing.

- [ ] **Step 3: Implement the config changes**

In `apps/qscat-run/qscat_run/config.py`:

1. Add the artifact flag to `ArtifactSpec` (after `eigenstates`):

```python
    # Emit the quasi-bound vibrational levels of the anion in the LCP complex
    # potential (`qscat.core.lcp.resonance_levels`) -- the BO approximation to
    # the resonance energies. LCP path only.
    resonance_levels: bool = False
```

and in `_load_artifacts`'s `return ArtifactSpec(...)`:

```python
        resonance_levels=bool(raw.get("resonance_levels", False)),
```

2. Add `nuclear_angle_b` to `GridSpec` (`config.py:172`, the config-level grid
   block — not `qscat.dvr.GridSpec`), after the existing `nuclear` field:

```python
    # Tail ECS angle of the SECOND nuclear grid used for two-angle level
    # selection. Defaults to `angle_a - 10` degrees (eMoScat's decks pair
    # 44/35 and 40/30). Real segments and quadrature are always copied from
    # grid a -- the shared real region is what makes the comparison valid.
    nuclear_angle_b: float | None = None
```

and parse it in `_load_grid` (`config.py:182`):

```python
        nuclear_angle_b=(
            float(raw["nuclear_angle_b"]) if raw.get("nuclear_angle_b") is not None else None
        ),
```

3. In the observable validity loop, the existing `presets.VALIDITY` lookup does
   the work once `presets` lists the kind — no change needed there.

4. Relax the LCP `da` requirement so a levels-only run is legal. Replace the
   `if not any(obs.kind == "da" ...)` guard with:

```python
        kinds = {obs.kind for obs in cfg.observables}
        if "da" not in kinds and "resonance_levels" not in kinds:
            raise ConfigError(
                "methods includes 'lcp' but no 'da' or 'resonance_levels' observable "
                "is requested; LCP approximates the DA cross section -- add "
                "`{kind: da, channels: 1}` -- or ask for the quasi-bound levels with "
                "`{kind: resonance_levels}`"
            )
```

In `apps/qscat-run/qscat_run/presets.py`:

5. Make `energies` optional for levels-only runs. In `resolve_defaults`
   (`presets.py:534`), replace the energy-filling line
   (`presets.py:556`)

```python
    energies = cfg.energies if cfg.energies is not None else preset.default_energies
```

with

```python
    # A levels-only run has no energy sweep to fill: `resonance_levels` needs
    # a molecule and two nuclear grids, nothing else. Leave `energies` None so
    # the runner can tell "not requested" from "requested but unresolved".
    needs_energies = any(obs.kind != "resonance_levels" for obs in cfg.observables)
    energies = cfg.energies
    if energies is None and needs_energies:
        energies = preset.default_energies
```

   (`_run_lcp` already raises a clear `ConfigError` if an energy-consuming path
   finds `cfg.energies is None`, so nothing downstream silently misbehaves.)

6. Add `"resonance_levels"` to `VALIDITY["NO"]` and `VALIDITY["F2"]` (the two
   molecules with `lcp_grids`). Leave `N2` and `H2P` alone.

7. Parameterize the nuclear decks and add the two helpers:

```python
def _no_nuc_grid(angle_deg: float = _NO_NUC_ANGLE) -> FemDvrEcsGrid:
    return segmented_grid(
        _NO_NUC_REAL, _NO_NUC_COMPLEX, angle_deg=angle_deg, quadrature=_NO_NUC_QUAD
    )


def _f2_nuc_grid(angle_deg: float = _F2_NUC_ANGLE) -> FemDvrEcsGrid:
    return segmented_grid(
        _F2_NUC_REAL, _F2_NUC_COMPLEX, angle_deg=angle_deg, quadrature=_F2_NUC_QUAD
    )


_NUC_GRID_BUILDERS = {"NO": _no_nuc_grid, "F2": _f2_nuc_grid}

# How far below grid a's tail angle grid b sits, when not set explicitly.
# eMoScat's electronic LCP decks pair 44/35 and 40/30 -- about ten degrees.
_ANGLE_B_OFFSET = 10.0


def nuclear_angle_b(cfg: ExperimentConfig) -> float:
    """The second nuclear grid's tail angle: explicit, or `angle_a - 10` deg.

    Always moves DOWNWARD, which is unconditionally safe against the ionic
    model's `max_nuclear_ecs_angle_deg` divergence bound.
    """
    if cfg.grid.nuclear_angle_b is not None:
        return float(cfg.grid.nuclear_angle_b)
    g_a, _ea, _eb = resolve_lcp_grids(cfg)
    return max(el.angle_deg for el in g_a.spec.elements) - _ANGLE_B_OFFSET


def nuclear_grid_at_angle(cfg: ExperimentConfig, angle_deg: float) -> FemDvrEcsGrid:
    """This molecule's nuclear deck rebuilt at a different ECS tail angle.

    Same real segments and quadrature -- only the tail rotates -- so every real
    node is shared with `resolve_lcp_grids`'s nuclear grid, which is what
    `qscat.core.lcp.lcp_resonance_levels` requires of its two grids.
    """
    builder = _NUC_GRID_BUILDERS.get(cfg.molecule)
    if builder is None:
        raise ConfigError(
            f"no LCP nuclear deck for {cfg.molecule}; "
            f"available: {sorted(_NUC_GRID_BUILDERS)}"
        )
    return builder(angle_deg)
```

Export `nuclear_angle_b` and `nuclear_grid_at_angle` from `presets.__all__`.
Import `ConfigError` there if it is not already imported.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run --no-sync pytest apps/qscat-run/tests/test_config.py -v`
Expected: PASS, including every pre-existing config test.

- [ ] **Step 5: Commit**

```bash
uv run ruff format apps/qscat-run
uv run ruff check apps/qscat-run
git add apps/qscat-run/qscat_run/config.py apps/qscat-run/qscat_run/presets.py apps/qscat-run/tests/test_config.py
git commit -m "feat(qscat-run): resonance_levels observable kind + second-angle nuclear grid"
```

---

### Task 7: Runner — compute the levels on both entry points

**Files:**
- Modify: `apps/qscat-run/qscat_run/runner.py`
- Test: `apps/qscat-run/tests/test_runner_resonance_levels.py` (create)

**Interfaces:**
- Consumes: `qscat.core.lcp.resonance_levels` + `ResonanceLevels` (Task 4);
  `presets.nuclear_angle_b`/`nuclear_grid_at_angle`,
  `cfg.artifacts.resonance_levels` (Task 6).
- Produces:
  - `ResonanceLevelsRun` frozen dataclass: `label: str`, `levels: ResonanceLevels`,
    `R_axis: NDArray[float64]`, `Vd: NDArray[complex128]`, `Gamma: NDArray[float64]`.
  - `ExperimentResult.resonance_levels: list[ResonanceLevelsRun]`.
  Task 8 writes them out.

- [ ] **Step 1: Write the failing test**

Create `apps/qscat-run/tests/test_runner_resonance_levels.py`:

```python
"""The `resonance_levels` observable + artifact flag through the runner."""

from __future__ import annotations

import numpy as np
import pytest
from qscat_run import presets
from qscat_run.config import load_config, validate_config
from qscat_run.runner import run_experiment


def _resolved(path):
    cfg = load_config(path)
    validate_config(cfg)
    return presets.resolve_defaults(cfg)

LEVELS_ONLY = (
    "molecule: F2\n"
    "methods: [lcp]\n"
    "observables:\n"
    "  - kind: resonance_levels\n"
    "    channels: 4\n"
)


@pytest.mark.slow
def test_levels_only_run_produces_levels_and_no_cross_sections(tmp_path):
    cfg_path = tmp_path / "levels.yaml"
    cfg_path.write_text(LEVELS_ONLY)
    result = run_experiment(_resolved(cfg_path))

    assert len(result.resonance_levels) == 1
    run = result.resonance_levels[0]
    assert run.label == "lcp:resonance_levels"
    assert run.levels.energies.size == 4
    assert np.all(np.diff(run.levels.energies.real) > 0)
    assert run.levels.widths.size == 4
    assert run.R_axis.size == run.Vd.size == run.Gamma.size
    assert not result.cross_sections  # no sweep was requested or run


@pytest.mark.slow
def test_artifact_flag_adds_levels_to_a_da_run(tmp_path):
    cfg_path = tmp_path / "da.yaml"
    cfg_path.write_text(
        "molecule: F2\n"
        "methods: [lcp]\n"
        "observables: [{kind: da, channels: 1}]\n"
        "energies: {min: 0.02, max: 0.04, step: 0.02}\n"
        "artifacts: {resonance_levels: true}\n"
    )
    result = run_experiment(_resolved(cfg_path))
    assert "lcp:da:ch0" in result.cross_sections
    assert len(result.resonance_levels) == 1
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --no-sync pytest apps/qscat-run/tests/test_runner_resonance_levels.py -v`
Expected: FAIL — `AttributeError: 'ExperimentResult' object has no attribute 'resonance_levels'`.

- [ ] **Step 3: Implement the runner changes**

In `apps/qscat-run/qscat_run/runner.py`:

1. Extend the `qscat.core.lcp` import to include `resonance_levels as
   lcp_resonance_levels_for_model` and `ResonanceLevels`.

2. Add the dataclass next to `ResonanceState`:

```python
@dataclass(frozen=True)
class ResonanceLevelsRun:
    """The anion's quasi-bound vibrational levels in the LCP complex potential.

    The Born-Oppenheimer approximation to the 2-D model's resonance energies:
    the thesis's `omega_j`, promoted from real levels of `Re V_res` to genuine
    complex eigenvalues `E_v - i Gamma_v/2`. Wraps the library dataclass rather
    than restating its fields, and carries the curve the levels sit in so the
    artifact can draw them together.

    `R_axis` is the real nuclear coordinate (bohr); `Vd`/`Gamma` are the LCP
    curve on the FULL grid (real nodes plus ECS tail), as `local_complex_potential`
    returns them.
    """

    label: str
    levels: ResonanceLevels
    R_axis: npt.NDArray[np.float64]
    Vd: npt.NDArray[np.complex128]
    Gamma: npt.NDArray[np.float64]
```

3. Add the field to `ExperimentResult`:

```python
    resonance_levels: list[ResonanceLevelsRun] = field(default_factory=list)
```

4. In `_run_lcp`, change the signature's return type to include
   `list[ResonanceLevelsRun]`, and make the energies guard conditional:

```python
    kinds = {obs.kind for obs in cfg.observables}
    wants_sigma = "da" in kinds
    if wants_sigma and cfg.energies is None:
        raise ConfigError("no energies resolved for this config (missing 'energies' block?)")
    energies = cfg.energies.as_array() if cfg.energies is not None else np.empty(0)
```

Guard the `vibrational_states`, `local_complex_potential` and
`lcp_da_cross_section` work that only the cross section needs behind
`wants_sigma` where it is not also needed by the levels path — `local_complex_potential`
is needed by both, so leave it unconditional.

5. Add the levels computation near the end of `_run_lcp`, before its `return`:

```python
    levels_runs: list[ResonanceLevelsRun] = []
    if "resonance_levels" in kinds or cfg.artifacts.resonance_levels:
        n_req = next(
            (_n_channels(o) for o in cfg.observables if o.kind == "resonance_levels"),
            None,
        )
        t0 = time.time()
        g_R_b = presets.nuclear_grid_at_angle(cfg, presets.nuclear_angle_b(cfg))
        levels = lcp_resonance_levels_for_model(
            model, g_R, g_R_b, elec_a, elec_b, n_levels=n_req
        )
        timings["lcp:resonance_levels"] = time.time() - t0
        levels_runs.append(
            ResonanceLevelsRun(
                label="lcp:resonance_levels",
                levels=levels,
                R_axis=g_R.real_points,
                Vd=v_d,
                Gamma=gamma,
            )
        )
```

6. Update `_run_lcp`'s `return` to add `levels_runs`, and update
   `run_experiment`'s unpacking of the `_run_lcp` call plus the
   `ExperimentResult(...)` construction to pass `resonance_levels=`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run --no-sync pytest apps/qscat-run/tests/test_runner_resonance_levels.py -v -m "slow or not slow"`
Expected: PASS. These are slow (~2-5 min): they run the real F2 electronic walk.

- [ ] **Step 5: Run the app suite for regressions**

Run: `uv run --no-sync pytest apps/qscat-run/tests -q`
Expected: PASS. The existing LCP DA tests must be unaffected.

- [ ] **Step 6: Commit**

```bash
uv run ruff format apps/qscat-run
uv run ruff check apps/qscat-run
git add apps/qscat-run/qscat_run/runner.py apps/qscat-run/tests/test_runner_resonance_levels.py
git commit -m "feat(qscat-run): compute resonance levels on the LCP path"
```

---

### Task 8: Artifacts — CSV, NPZ and the levels-on-the-curve figure

**Files:**
- Modify: `apps/qscat-run/qscat_run/artifacts.py`
- Test: `apps/qscat-run/tests/test_artifacts_resonance_levels.py` (create)

**Interfaces:**
- Consumes: `ResonanceLevelsRun` (Task 7).
- Produces: `_write_resonance_levels(out_dir, run)` writing
  `resonance_levels_{label}.csv`, `.npz`, `.png`; called from `write_artifacts`.

- [ ] **Step 1: Write the failing test**

Create `apps/qscat-run/tests/test_artifacts_resonance_levels.py`:

```python
"""Artifact writers for the quasi-bound level table."""

from __future__ import annotations

import csv

import numpy as np
from qscat.core.lcp import ResonanceLevels
from qscat_run.artifacts import _write_resonance_levels
from qscat_run.runner import ResonanceLevelsRun


def _fake_run() -> ResonanceLevelsRun:
    n_grid = 12
    levels = ResonanceLevels(
        energies=np.array([-0.05 - 0.001j, -0.03 - 0.004j]),
        widths=np.array([0.002, 0.008]),
        states=np.ones((2, n_grid), dtype=np.complex128) / np.sqrt(n_grid),
        residuals=np.array([1e-9, 4e-9]),
        real_weight=np.array([0.999, 0.97]),
        golden_rule=np.array([-0.0501 - 0.0011j, -0.0299 - 0.0039j]),
    )
    return ResonanceLevelsRun(
        label="lcp:resonance_levels",
        levels=levels,
        R_axis=np.linspace(1.5, 6.0, n_grid),
        Vd=np.linspace(-0.1, 0.0, n_grid).astype(np.complex128),
        Gamma=np.linspace(0.01, 0.0, n_grid),
    )


def test_writes_csv_npz_and_png(tmp_path):
    _write_resonance_levels(tmp_path, _fake_run())
    stem = "resonance_levels_lcp_resonance_levels"
    assert (tmp_path / f"{stem}.csv").exists()
    assert (tmp_path / f"{stem}.npz").exists()
    assert (tmp_path / f"{stem}.png").exists()


def test_csv_columns_and_values(tmp_path):
    _write_resonance_levels(tmp_path, _fake_run())
    path = tmp_path / "resonance_levels_lcp_resonance_levels.csv"
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    assert list(rows[0]) == [
        "v", "Re_E", "Gamma_v", "residual", "real_weight", "Re_E0", "Gamma_v_1",
    ]
    assert len(rows) == 2
    assert float(rows[0]["v"]) == 0
    assert float(rows[0]["Re_E"]) == -0.05
    assert float(rows[1]["Gamma_v"]) == 0.008


def test_npz_round_trips_the_complex_energies_and_curve(tmp_path):
    run = _fake_run()
    _write_resonance_levels(tmp_path, run)
    data = np.load(tmp_path / "resonance_levels_lcp_resonance_levels.npz")
    np.testing.assert_allclose(data["energies"], run.levels.energies)
    np.testing.assert_allclose(data["states"], run.levels.states)
    np.testing.assert_allclose(data["R_axis"], run.R_axis)
    np.testing.assert_allclose(data["Vd"], run.Vd)
    np.testing.assert_allclose(data["Gamma"], run.Gamma)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --no-sync pytest apps/qscat-run/tests/test_artifacts_resonance_levels.py -v`
Expected: FAIL — `ImportError: cannot import name '_write_resonance_levels'`.

- [ ] **Step 3: Implement the writer**

In `apps/qscat-run/qscat_run/artifacts.py`, after `_write_resonance_state`:

```python
def _write_resonance_levels(out_dir: Path, run: ResonanceLevelsRun) -> None:
    """`resonance_levels_{label}.{csv,npz,png}` -- the quasi-bound level table.

    The png draws each level as a horizontal bar across the `V_d(R)` curve at
    its `Re E_v`, with bar thickness proportional to `Gamma_v` -- so a broad,
    short-lived level reads as a thick smear and a long-lived one as a hairline.
    """
    import csv

    stem = f"resonance_levels_{run.label.replace(':', '_')}"
    lv = run.levels

    with (out_dir / f"{stem}.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["v", "Re_E", "Gamma_v", "residual", "real_weight", "Re_E0", "Gamma_v_1"])
        for v in range(lv.energies.size):
            w.writerow([
                v,
                f"{lv.energies[v].real!r}",
                f"{lv.widths[v]!r}",
                f"{lv.residuals[v]!r}",
                f"{lv.real_weight[v]!r}",
                f"{lv.golden_rule[v].real!r}",
                f"{-2.0 * lv.golden_rule[v].imag!r}",
            ])

    np.savez(
        out_dir / f"{stem}.npz",
        energies=lv.energies,
        widths=lv.widths,
        states=lv.states,
        residuals=lv.residuals,
        real_weight=lv.real_weight,
        golden_rule=lv.golden_rule,
        R_axis=run.R_axis,
        Vd=run.Vd,
        Gamma=run.Gamma,
    )

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return  # plotting is the optional `plot` extra; csv/npz are the data

    n_real = run.R_axis.size
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(run.R_axis, run.Vd[:n_real].real, label=r"$V_d(R)$", color="tab:blue")
    ax.plot(run.R_axis, run.Gamma[:n_real], label=r"$\Gamma(R)$", color="tab:red", ls="--")
    for v in range(lv.energies.size):
        ax.axhline(
            lv.energies[v].real,
            color="k",
            lw=max(0.6, 400.0 * float(lv.widths[v])),
            alpha=0.55,
        )
        ax.annotate(
            rf"$\omega_{{{v}}}$",
            xy=(run.R_axis[-1], lv.energies[v].real),
            xytext=(-24, 2),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_xlabel(r"$R$ [bohr]")
    ax.set_ylabel("energy [hartree]")
    ax.set_title(f"quasi-bound levels -- {run.label}")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / f"{stem}.png", dpi=150)
    plt.close(fig)
```

Import `ResonanceLevelsRun` from `.runner` alongside the existing
`ResonanceState`/`EigenStates` imports, and call it from `write_artifacts`:

```python
    for run in result.resonance_levels:
        _write_resonance_levels(out_dir, run)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run --no-sync pytest apps/qscat-run/tests/test_artifacts_resonance_levels.py -v`
Expected: PASS, 3 tests.

- [ ] **Step 5: Commit**

```bash
uv run ruff format apps/qscat-run
uv run ruff check apps/qscat-run
git add apps/qscat-run/qscat_run/artifacts.py apps/qscat-run/tests/test_artifacts_resonance_levels.py
git commit -m "feat(qscat-run): resonance-levels csv/npz/png artifacts"
```

---

### Task 9: Example config, physics note, repo map

**Files:**
- Create: `apps/qscat-run/examples/f2-resonance-levels.yaml`
- Create: `docs/physics/lcp-resonance-levels.md`
- Modify: `CLAUDE.md`
- Modify: `apps/qscat-run/README.md`
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: everything above. Produces no code.

- [ ] **Step 1: Write the example config**

Create `apps/qscat-run/examples/f2-resonance-levels.yaml`:

```yaml
# Quasi-bound vibrational levels of the F2 anion in the LCP complex potential
# V_res(R) = E_res(R) - (i/2) Gamma(R) -- the Born-Oppenheimer approximation to
# the 2-D model's resonance energies (the thesis's omega_j, Vana 2017 Sec. 1.5).
#
# This is the cheapest run in the app: no energy sweep, no wavepacket, no
# cross section. It solves the fixed-R electronic pole walk once, lays the
# curve onto two nuclear grids differing only in ECS tail angle, and
# diagonalizes.
#
#   uv run qscat-run apps/qscat-run/examples/f2-resonance-levels.yaml -o out/
#
# Outputs: resonance_levels_lcp_resonance_levels.{csv,npz,png}

molecule: F2
methods: [lcp]

observables:
  # `channels` is the number of levels to report; omit it for every
  # angle-stable level inside the default window.
  - kind: resonance_levels
    channels: 6

grid:
  preset: default
  # Tail angle of the second nuclear grid. Omit to use `angle_a - 10` degrees.
  # Always keep it BELOW angle_a.
  # nuclear_angle_b: 15.0
```

- [ ] **Step 2: Verify the example runs**

Run: `uv run --no-sync pytest apps/qscat-run/tests/test_examples.py -q`
Expected: PASS — that test validates every example config parses and resolves.

Then a real run: `uv run qscat-run apps/qscat-run/examples/f2-resonance-levels.yaml -o /tmp/levels-check`
Expected: exits 0; `/tmp/levels-check/resonance_levels_lcp_resonance_levels.csv`
holds 6 rows with strictly increasing `Re_E` and non-negative `Gamma_v`.
**Record the actual six `(Re_E, Gamma_v)` pairs** — they go into the physics note
in Step 3 and are the first real numbers this capability has produced.

- [ ] **Step 3: Write the physics note**

Create `docs/physics/lcp-resonance-levels.md` covering, in this order:

1. **What this computes** — `H_N = T(mu) + V_d(R) - i Gamma(R)/2` on the nuclear
   FEM-DVR-ECS grid; eigenvalues `E_v - i Gamma_v/2`. Cite Váňa 2017 Eqs.
   1.63/1.65 and note qscat's `Vd` is the thesis's `E_res`.
2. **What eMoScat and the thesis actually did** — `TimeDependentModel2d.cpp:58-79`
   discards `Im(V_res)` and takes the 15 lowest `Re E` states blindly as a
   projection basis; those are the thesis's `omega_j`. Our `golden_rule` column
   reproduces them; the complex widths are the extension. The thesis reports
   widths only qualitatively and estimates lifetimes from TD peak-formation times.
3. **Terminology** — these are complex-scaled (ECS) resonance eigenstates, NOT
   Siegert pseudostates. Include the four-row comparison table from the spec and
   cite Hvizdoš et al., Phys. Rev. A **97**, 022704 (2018), App. A.
4. **Normalization** — bilinear c-product over the whole rotated grid; the
   pseudostates' surface term (their Eq. A5) is an artifact of truncating at a
   finite radius, which ECS does not do.
5. **Selection** — two nuclear ECS angles, `match_angle_stable`; why the real
   nodes must be shared; why the electronic walk runs only once.
6. **The ECS angle bound** — nuclear `< pi/8`, electronic `< pi/4` for H₂⁺, with
   the `a3 * R**4` / `exp(-r^2/3)` reason.
7. **Validation** — the two analytic oracles (Morse spectrum, constant-Γ rigid
   shift), angle/grid independence, golden-rule consistency; and the honest
   statement that **no prior computed width exists anywhere** to check `Gamma_v`
   against, so the imaginary parts are gated only by internal consistency until
   the M3 comparisons land.
8. **First results** — the F2 table from Step 2.
9. **Limits and what is next** — M3 (thesis `omega_j` reproduction, elastic-peak
   correspondence, NO lifetime bounds, both grid decks) and the deferred
   `Psi_d` projection / level populations.

- [ ] **Step 4: Update the repo map and changelog**

In `CLAUDE.md`, extend the `qscat.core` bullet's `lcp` sentence with:

```
Also `resonance_levels`/`lcp_resonance_levels` (+ `ResonanceLevels`): the
BORN-OPPENHEIMER approximation to the resonance energies -- the nuclear
eigenvalue problem IN the complex curve, `H_N = T(mu) + V_d - i*Gamma/2` on
the nuclear FEM-DVR-ECS grid, giving complex quasi-bound levels
`E_v - i*Gamma_v/2` (the thesis's `omega_j`, promoted from eMoScat's
real-part-only levels). Physical levels are picked by two NUCLEAR ECS angles
(`qscat.ecs.match_angle_stable`, the multi-state sibling of
`find_resonance_pole`); the electronic pole walk runs ONCE since `E_res(R)` at
real `R` is angle-independent. A golden-rule comparator (`Gamma=0` levels +
`<chi|Gamma|chi>`) rides along and reproduces what eMoScat/the thesis computed
-- its divergence from the complex result is the non-perturbative signal. NOT
Siegert pseudostates (see docs/physics/lcp-resonance-levels.md).
```

In `apps/qscat-run/README.md`, add `resonance_levels` to the observables matrix
with a note that it needs no `energies` block, and list the three artifacts.

In `CHANGELOG.md`, add an entry under the unreleased section describing the new
observable, the library API, and the H₂⁺ reduced-mass correction.

- [ ] **Step 5: Full verification**

```bash
uv run ruff check .
uv run ruff format --check .
uv run mypy libs/qscat
uv run --no-sync pytest -q
```
Expected: all clean, all pass. Report the actual test count.

- [ ] **Step 6: Commit**

```bash
git add apps/qscat-run/examples/f2-resonance-levels.yaml docs/physics/lcp-resonance-levels.md CLAUDE.md apps/qscat-run/README.md CHANGELOG.md
git commit -m "docs(lcp): resonance-levels example, physics note, repo map"
```

---

## After the plan

Run the `physics-reviewer` agent over the diff before merging — the repo's
promotion rule for anything landing in `qscat`. Ask it specifically about: the
complex-symmetric normalization, the golden-rule pairing, the default window, and
whether the `Gamma` support-condition warning has the right threshold.

Milestone 3 (reproducing the thesis `omega_j`, the elastic-peak correspondence,
the NO lifetime bounds, and running both grid decks as a convergence check) is a
**separate plan** — it needs data the author holds and the thesis-grid preset
variants, neither of which is a prerequisite for anything above.
