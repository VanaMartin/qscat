# N₂ Exact 2-D TI Cross-Section Implementation Plan (sub-project #6)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Solve the electron–N₂ problem exactly in 2-D via the driven Lippmann–Schwinger equation on a tensor-product FEM-DVR-ECS grid, verify it against Houfek's data, then use it as the oracle to measure where the 1-D LCP approximation fails.

**Architecture:** `projects/n2_2d_cross_section/` — a parametrized electronic grid, the energy-normalized Riccati–Bessel channel function, the 2-D Hamiltonian, one sparse LU per energy, and a T-matrix projection. Built entirely on sub-project #5's `qscat.dvr` tensor layer and `qscat.linalg`; reuses the vibrational states and model potentials from #2/#3. Nothing new is promoted into `qscat`.

**Tech Stack:** Python 3.12, numpy, `scipy.sparse` / `scipy.sparse.linalg.splu`, `scipy.special.spherical_jn`.

**Design spec:** `docs/superpowers/specs/2026-07-23-n2-2d-exact-cross-section-design.md`
**Source archaeology:** `.superpowers/sdd/n2-2d-exact-extraction.md`

## Global Constraints

- Python `>=3.12`. Everything through `uv` — never bare `python`/`pip`/conda. Conda `base` may be active; `uv run` correctly ignores it.
- **Atomic units throughout.** `μ = 12766.36`, electron mass `1`, `l = 2` (`PARAMS["impulsemomentum"]`). Cross sections in bohr².
- Package-absolute imports only. `projects/` may import `projects.*`; **`projects/` must NOT import `validation/`** (that dependency runs the other way). Use `projects.n2_resonance.potential` (`v0`, `lam`, `v_int`, `PARAMS`) — the projects-side copy.
- `uv run mypy libs/qscat` stays at **0 errors**. `uv run ruff check .` stays clean (line length 100; rules `E, F, I, UP, B, NPY`).
- The existing N₂ harness must not regress: **19 PASS / 0 PENDING / 2 NOTE / 0 FAIL**, exit 0. (Group E may *add* rows; the 2 NOTEs may become PASS if the exact model closes them — that would be the headline result, not a regression.)
- **`V_int` is the interaction ALONE:** `V_int(r,R) = -λ(R)·exp(-α_c r²)`. `v0(R)` and `l(l+1)/2r²` belong to the entrance channel, **not** the perturbation. Putting them in the driving term is a physics error that still produces plausible-looking numbers.
- **ECS masking:** `Ψ_i` and every channel projection must be zeroed outside the unscaled region via `TensorGrid.real_mask()`. A projection extending onto the scaled tail is meaningless.
- **DVR coefficients:** a *function* becomes basis coefficients as `c_j = f(x_j)·sqrt(w_j)` using `TensorGrid.sqrt_weights()` (bridge-summed, complex). With both sides in coefficient form, `c_product` *is* the quadrature integral — no extra weights.
- **`H` is complex symmetric, never Hermitian.** Use `qscat.linalg.c_product`, never `np.vdot`.
- Index order is numpy-native C-order (last axis fastest). Axis 0 = electronic `r`, axis 1 = nuclear `R`.
- **The model is a given testbed, not a description of real N₂.** Never tune a model parameter to improve agreement with anything. Houfek's data is the *gate* that certifies our solver; once it passes, the exact result is the **oracle** and the LCP is **under test**, so an LCP-vs-exact discrepancy is the *result*, not an error to minimize.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`

---

## File Structure

| File | Responsibility |
|---|---|
| `projects/n2_2d_cross_section/__init__.py` | **Create.** Empty package marker |
| `projects/n2_2d_cross_section/electronic_grid.py` | **Create.** Parametrized electronic FEM-DVR-ECS grid factory |
| `projects/n2_2d_cross_section/channels.py` | **Create.** Energy-normalized Riccati–Bessel `F_{E,l}`; masked channel vectors |
| `projects/n2_2d_cross_section/hamiltonian2d.py` | **Create.** The 2-D potential surface and `H_2D` assembly |
| `projects/n2_2d_cross_section/cross_section_2d.py` | **Create.** Driven solve, T-matrix, σ (the crux) |
| `projects/n2_2d_cross_section/convergence.py` | **Create.** Convergence study driver (runnable module) |
| `projects/n2_2d_cross_section/test_*.py` | **Create.** Per-task tests |
| `validation/n2/exact2d.py` | **Create.** Group E: anchors, LCP head-to-head, the two NOTEs |
| `validation/n2/experiment.py` | **Modify.** Add group E rows |
| `docs/physics/n2-2d-cross-section.md` | **Create.** Method, convergence table, LCP-vs-exact findings |
| `CLAUDE.md` | **Modify.** Add the sub-project to the repo map |

---

### Task 1: Electronic grid + energy-normalized channel function

**Files:**
- Create: `projects/n2_2d_cross_section/__init__.py` (empty), `projects/n2_2d_cross_section/electronic_grid.py`, `projects/n2_2d_cross_section/channels.py`
- Test: `projects/n2_2d_cross_section/test_channels.py`

**Interfaces:**
- Consumes: `qscat.dvr.{ElementSpec, GridSpec, FemDvrEcsGrid}`.
- Produces:
  - `n2_electronic_grid(*, r_max: float = 30.0, angle_deg: float = 35.0, order: int = 8, n_complex: int = 8, tail_alpha: float = 0.2, tail_skip: int = 2) -> FemDvrEcsGrid`
  - `riccati_bessel_en(r: NDArray[float64], k: float, l: int) -> NDArray[float64]`

**Background you need.** The electronic grid mirrors eMoScat's `N2-model.json` layout: real segments with element lengths `0.2` up to `r=1`, `1.0` to `r=5`, `2.0` to `r=7`, `3.0` to `r=10`, then `4.0` out to `r_max`; then `n_complex` ECS elements whose lengths grow exponentially (`skip` elements at the base length, then `base·exp(alpha·j)`), all at one angle. It is **parametrized** because Task 4 varies every one of these in a convergence study.

`F_{E,l}(r) = sqrt(2k/π)·r·j_l(kr)` is the **energy-normalized** regular free solution (`m = 1`), matching eMoScat's `sphBesselJEn` (`source/bessel.cpp:50`).

**Critical:** `scipy.special.spherical_jn` accepts **real arguments only**. That is fine and not a limitation — `F` is only ever evaluated on the **unscaled** region (it is masked to zero on the ECS tail, exactly as eMoScat does), where `grid.points == grid.real_points` and both are real. Evaluate on real points and zero the tail; never call `spherical_jn` with a complex argument.

- [ ] **Step 1: Write the failing tests — `projects/n2_2d_cross_section/test_channels.py`**

```python
"""Energy-normalized regular free radial function and the electronic grid.

`F_{E,l}(r) = sqrt(2k/pi) * r * j_l(k r)` is the energy-normalized regular
solution of the free radial equation at electron mass 1 (eMoScat
`sphBesselJEn`, source/bessel.cpp:50).
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.dvr import FemDvrEcsGrid

from projects.n2_2d_cross_section.channels import riccati_bessel_en
from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid


def test_l0_reduces_to_exact_normalized_sine() -> None:
    """j_0(x) = sin(x)/x, so F_{E,0}(r) == sqrt(2/(pi k)) sin(k r) EXACTLY.

    This pins the normalization CONSTANT, not just the shape -- the whole
    cross-section scale rides on it.
    """
    r = np.linspace(0.05, 40.0, 2000)
    k = 0.63
    got = riccati_bessel_en(r, k, 0)
    want = np.sqrt(2.0 / (np.pi * k)) * np.sin(k * r)
    assert np.abs(got - want).max() < 1e-12


def test_asymptotic_envelope_matches_energy_normalization() -> None:
    """F -> sqrt(2/(pi k)) sin(k r - l pi/2): the envelope fixes the constant."""
    k = 0.5
    r = np.linspace(200.0, 260.0, 20000)   # k r >> l(l+1), deep asymptotic
    f = riccati_bessel_en(r, k, 2)
    assert abs(np.abs(f).max() - np.sqrt(2.0 / (np.pi * k))) < 1e-4


def test_satisfies_free_radial_equation() -> None:
    """-F'' + l(l+1)/r^2 F = k^2 F, checked by finite differences."""
    k, ell = 0.7, 2
    h = 1e-4
    r = np.linspace(2.0, 12.0, 300)
    f = riccati_bessel_en(r, k, ell)
    fpp = (riccati_bessel_en(r + h, k, ell) - 2 * f + riccati_bessel_en(r - h, k, ell)) / h**2
    residual = -fpp + ell * (ell + 1) / r**2 * f - k**2 * f
    assert np.abs(residual).max() < 1e-5 * np.abs(f).max()


def test_regular_at_origin() -> None:
    """The REGULAR solution vanishes at r -> 0 like r^{l+1}."""
    assert abs(float(riccati_bessel_en(np.array([1e-6]), 0.5, 2)[0])) < 1e-15


def test_electronic_grid_shape_and_ecs_pivot() -> None:
    g = n2_electronic_grid(r_max=30.0, angle_deg=35.0, order=8, n_complex=8)
    assert isinstance(g, FemDvrEcsGrid)
    # R0 is x_min + sum(real element lengths), accumulated in floating point,
    # so compare approximately rather than exactly.
    assert g.R0 == pytest.approx(30.0)        # pivot at the end of the real region
    assert g.real_points.min() > 0.0          # Dirichlet endpoint at r=0 dropped
    # real region genuinely unscaled; tail genuinely scaled
    inside = g.real_points <= g.R0
    assert np.abs(g.points[inside].imag).max() < 1e-12
    assert np.abs(g.points[~inside].imag).max() > 1.0


def test_electronic_grid_is_parametrized() -> None:
    """Task 4's convergence study varies every one of these."""
    a = n2_electronic_grid(r_max=20.0, order=7, n_complex=6)
    b = n2_electronic_grid(r_max=45.0, order=9, n_complex=10)
    assert a.n != b.n
    assert a.R0 == 20.0 and b.R0 == 45.0
    assert n2_electronic_grid(angle_deg=25.0).points[-1] != n2_electronic_grid(
        angle_deg=40.0
    ).points[-1]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest projects/n2_2d_cross_section/test_channels.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'projects.n2_2d_cross_section'`

- [ ] **Step 3: Implement `projects/n2_2d_cross_section/electronic_grid.py`**

```python
"""Parametrized electronic FEM-DVR-ECS grid for the exact 2-D N2 model.

Layout follows eMoScat's `input/experimental/N2-model.json` `grids.electronic`:
a finely-resolved region near the origin where the interaction
`-lambda(R) exp(-alpha_c r^2)` lives, coarsening outward, then an ECS tail of
exponentially growing elements at a single angle.

EVERY parameter is exposed because sub-project #6's convergence study (Task 4)
varies all of them -- eMoScat asserted 35 degrees and a 98-bohr box without
ever documenting the study that justified them, so we redo it.
"""

from __future__ import annotations

import numpy as np
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec

__all__ = ["n2_electronic_grid"]

# (segment_end_bohr, element_length_bohr); the final segment runs to r_max.
_INNER_SEGMENTS: tuple[tuple[float, float], ...] = ((1.0, 0.2), (5.0, 1.0), (7.0, 2.0), (10.0, 3.0))
_OUTER_LENGTH = 4.0


def _ecs_tail(base: float, n: int, *, skip: int, alpha: float) -> list[float]:
    """eMoScat `uniform_increment`/`exp`: `skip` elements at `base`, then growing."""
    return [base if i < skip else base * float(np.exp(alpha * (i - skip + 1))) for i in range(n)]


def n2_electronic_grid(
    *,
    r_max: float = 30.0,
    angle_deg: float = 35.0,
    order: int = 8,
    n_complex: int = 8,
    tail_alpha: float = 0.2,
    tail_skip: int = 2,
) -> FemDvrEcsGrid:
    """Electronic radial grid: real region [0, r_max] + an ECS tail at `angle_deg`.

    The ECS pivot is `R0 == r_max` by construction.
    """
    if r_max <= _INNER_SEGMENTS[-1][0]:
        raise ValueError(f"r_max must exceed {_INNER_SEGMENTS[-1][0]} bohr, got {r_max}")

    elements: list[ElementSpec] = []
    start = 0.0
    for end, length in _INNER_SEGMENTS:
        k = round((end - start) / length)
        elements += [ElementSpec((end - start) / k) for _ in range(k)]
        start = end

    k_out = max(1, round((r_max - start) / _OUTER_LENGTH))
    elements += [ElementSpec((r_max - start) / k_out) for _ in range(k_out)]

    base = (r_max - start) / k_out
    elements += [
        ElementSpec(h, angle_deg)
        for h in _ecs_tail(base, n_complex, skip=tail_skip, alpha=tail_alpha)
    ]
    return FemDvrEcsGrid(GridSpec(quadrature=order, elements=elements, x_min=0.0))
```

- [ ] **Step 4: Implement `projects/n2_2d_cross_section/channels.py`**

```python
"""Asymptotic channel functions for the exact 2-D e-N2 scattering problem.

The entrance/exit channel of a VE transition is a free electron of momentum
`k` in partial wave `l`, times a neutral vibrational state. The electronic
factor is the ENERGY-NORMALIZED regular free radial solution

    F_{E,l}(r) = sqrt(2/(pi k)) (k r) j_l(k r) = sqrt(2 k / pi) r j_l(k r)

at electron mass 1 -- eMoScat's `sphBesselJEn` (`source/bessel.cpp:50`),
equivalently `sF_en` (`source/coulomb.cpp:75`) with charge 0 since N2 is
neutral. Energy normalization (`<F_E|F_E'> = delta(E-E')`) is what makes the
`sigma = 4 pi^3 |T|^2 / k^2` prefactor correct; getting the constant wrong
rescales every cross section.

`scipy.special.spherical_jn` takes REAL arguments only. That is not a
limitation here: `F` is only ever needed on the UNSCALED region, because a
channel projection that extends onto the exterior-complex-scaled tail is
meaningless and must be masked to zero there anyway (eMoScat zeroes it
explicitly -- `time_independent_model.cpp:149-151`). So `F` is evaluated on
real points and the tail is zeroed; no complex Bessel function is needed.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.special import spherical_jn

__all__ = ["riccati_bessel_en"]


def riccati_bessel_en(
    r: npt.NDArray[np.float64], k: float, l: int
) -> npt.NDArray[np.float64]:
    """`F_{E,l}(r) = sqrt(2k/pi) r j_l(k r)`, energy-normalized at mass 1.

    `r` must be REAL (see module docstring); `k = sqrt(2E) > 0`.
    """
    if k <= 0.0:
        raise ValueError(f"k must be positive, got {k}")
    rr = np.asarray(r, dtype=np.float64)
    out: npt.NDArray[np.float64] = np.sqrt(2.0 * k / np.pi) * rr * spherical_jn(l, k * rr)
    return out
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest projects/n2_2d_cross_section/test_channels.py -q`
Expected: PASS (6 tests)

- [ ] **Step 6: Lint**

Run: `uv run ruff check .`
Expected: clean.

- [ ] **Step 7: Commit**

```bash
git add projects/n2_2d_cross_section
git commit -m "$(cat <<'EOF'
feat(n2-2d): parametrized electronic grid + energy-normalized channel function

The electronic grid mirrors eMoScat's N2-model.json layout but exposes every
parameter (box size, ECS angle, DVR order, tail length/growth), because the
convergence study redoes the box/angle study eMoScat asserted but never
documented.

F_{E,l}(r) = sqrt(2k/pi) r j_l(kr) is the energy-normalized regular free
solution; the normalization constant is what makes the 4 pi^3 |T|^2 / k^2
prefactor correct, so it is pinned by an exact l=0 identity against
sqrt(2/(pi k)) sin(kr) rather than only by shape.

scipy's spherical_jn is real-argument-only, which costs nothing here: a
channel projection is meaningless on the complex-scaled tail and is masked to
zero there regardless, exactly as the reference does.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: The 2-D Hamiltonian

**Files:**
- Create: `projects/n2_2d_cross_section/hamiltonian2d.py`
- Test: `projects/n2_2d_cross_section/test_hamiltonian2d.py`

**Interfaces:**
- Consumes: `qscat.dvr.{TensorGrid, hamiltonian_nd}`, `projects.n2_resonance.potential.{v0, lam, v_int, PARAMS}`, Task 1's `n2_electronic_grid`.
- Produces:
  - `MU: float`, `ELL: int`
  - `potential_2d(r, R) -> ArrayLike` — the full surface
  - `interaction_2d(r, R) -> ArrayLike` — `V_int` ALONE
  - `build_h2d(tgrid: TensorGrid) -> sparse.csr_matrix`
  - `interaction_diag(tgrid: TensorGrid) -> NDArray[complex128]` — flat `V_int` on the tensor grid

**Background you need.** The full surface is
`V(r,R) = v0(R) + l(l+1)/(2r²) - λ(R)e^{-α_c r²}` (eMoScat `Potentials2d.cpp:18`) — exactly `v0` plus `v_eff_el` from the already-verified model module. **But the perturbation driving the scattering is `V_int = -λ(R)e^{-α_c r²}` alone**; `v0(R)` and the centrifugal term are entrance-channel potentials. Two different functions, both needed, easy to confuse — hence separate names.

`potential_2d` must **not** coerce to a real dtype: `r` and `R` are complex on the ECS tails, and discarding the imaginary part destroys the analytic continuation. (`validation/n2/model.py` documents this exact bug being found and fixed.)

- [ ] **Step 1: Write the failing tests — `projects/n2_2d_cross_section/test_hamiltonian2d.py`**

```python
"""The exact 2-D N2 potential surface and Hamiltonian."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import (
    ELL,
    MU,
    build_h2d,
    interaction_2d,
    interaction_diag,
    potential_2d,
)
from projects.n2_resonance.potential import lam, v0
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid


def _small_tgrid() -> TensorGrid:
    """Deliberately tiny: this task tests STRUCTURE, not converged physics."""
    return TensorGrid(
        [
            n2_electronic_grid(r_max=14.0, order=6, n_complex=4),
            n2_nuclear_grid(quadrature=8, r_max=20.0, n_complex=4),
        ]
    )


def test_potential_decomposes_into_channel_plus_interaction() -> None:
    """V = [v0(R) + l(l+1)/2r^2] + V_int -- the split the driven equation needs."""
    r, R = 1.7 + 0.0j, 2.1 + 0.0j
    channel = v0(R) + ELL * (ELL + 1) / (2 * r**2)
    assert complex(potential_2d(r, R)) == complex(channel + interaction_2d(r, R))


def test_interaction_excludes_v0_and_centrifugal() -> None:
    """The classic error: sweeping the channel potentials into the perturbation."""
    from projects.n2_resonance.potential import PARAMS

    alpha_c = PARAMS["potential"]["alpha_c"]
    r, R = 1.3 + 0.0j, 2.4 + 0.0j
    assert complex(interaction_2d(r, R)) == complex(-lam(R) * np.exp(-alpha_c * r**2))
    # V_int decays in r; v0(R) does not vanish where V_int does
    assert abs(complex(interaction_2d(30.0 + 0j, R))) < 1e-100
    assert abs(complex(v0(R))) > 1e-3


def test_potential_preserves_complex_points_on_the_ecs_tail() -> None:
    """Coercing to float here would silently destroy the analytic continuation."""
    tg = _small_tgrid()
    r, R = tg.points()
    vals = np.asarray(potential_2d(r, R))
    assert np.abs(np.broadcast_to(vals, tg.shape).imag).max() > 1e-6


def test_h2d_is_complex_symmetric_never_hermitian() -> None:
    tg = _small_tgrid()
    H = build_h2d(tg)
    assert isinstance(H, sp.csr_matrix)
    assert H.shape == (tg.size, tg.size)
    assert abs(H - H.T).max() < 1e-9 * abs(H).max()
    assert abs(H - H.conj().T).max() > 1e-3 * abs(H).max()


def test_interaction_diag_matches_pointwise_evaluation() -> None:
    tg = _small_tgrid()
    r, R = tg.points()
    want = np.broadcast_to(np.asarray(interaction_2d(r, R)), tg.shape).ravel()
    assert np.allclose(interaction_diag(tg), want, rtol=0, atol=1e-14)


def test_masses_are_on_the_right_axes() -> None:
    """Axis 0 is the ELECTRON (mass 1), axis 1 the nuclei (mass mu). Swapping
    them changes the spectrum by orders of magnitude."""
    tg = _small_tgrid()
    from qscat.dvr import kinetic_nd

    right = abs(kinetic_nd(tg, [1.0, MU])).max()
    wrong = abs(kinetic_nd(tg, [MU, 1.0])).max()
    assert right / wrong > 100.0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest projects/n2_2d_cross_section/test_hamiltonian2d.py -q`
Expected: FAIL — `ModuleNotFoundError: ... hamiltonian2d`

Note this task also needs `n2_nuclear_grid` to accept `quadrature`, `r_max` and `n_complex` keywords — it already does (`projects/n2_ti_cross_section/nuclear_grid.py`). Verify before implementing; if a keyword is missing, add it without changing the defaults.

- [ ] **Step 3: Implement `projects/n2_2d_cross_section/hamiltonian2d.py`**

```python
"""The exact 2-D electron-N2 potential surface and Hamiltonian.

    H_2D = -(1/2) d^2/dr^2 - (1/2 mu) d^2/dR^2
           + v0(R) + l(l+1)/(2 r^2) - lambda(R) exp(-alpha_c r^2)

verbatim eMoScat `Neutral2dPotential` (`source/Model2d/Potentials2d.cpp:18`),
built from the same model functions already ported and verified in
sub-project #2. No new potential physics.

TWO potentials live here and confusing them is a physics error, not a
convention choice:

- `potential_2d` is the FULL surface, which goes into `H_2D`.
- `interaction_2d` is `V_int = -lambda(R) exp(-alpha_c r^2)` ALONE -- the only
  perturbation relative to the entrance channel. `v0(R)` (the neutral
  molecule's own potential) and `l(l+1)/2r^2` (the centrifugal barrier) are
  CHANNEL potentials: they survive as `r -> infinity`, and the asymptotic
  channel function is an eigenfunction of the Hamiltonian containing them.
  Sweeping them into the driving term would produce a plausible-looking but
  wrong T-matrix.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from qscat.dvr import TensorGrid, hamiltonian_nd, potential_nd

from projects.n2_resonance.potential import PARAMS, v0, v_int

__all__ = [
    "MU",
    "ELL",
    "potential_2d",
    "interaction_2d",
    "build_h2d",
    "interaction_diag",
]

MU: float = 12766.36                       # N2 nuclear reduced mass (a.u.)
ELL: int = int(PARAMS["impulsemomentum"])  # fixed partial wave, l = 2


def interaction_2d(r: npt.ArrayLike, R: npt.ArrayLike) -> npt.ArrayLike:
    """`V_int(r,R) = -lambda(R) exp(-alpha_c r^2)` -- the perturbation ALONE."""
    return v_int(r, R)


def potential_2d(r: npt.ArrayLike, R: npt.ArrayLike) -> npt.ArrayLike:
    """The full surface `v0(R) + l(l+1)/(2 r^2) + V_int(r,R)`.

    Must not coerce to a real dtype: `r`/`R` are complex on the ECS tails.
    """
    rr = np.asarray(r)
    return v0(R) + ELL * (ELL + 1) / (2.0 * rr**2) + v_int(rr, R)


def interaction_diag(tgrid: TensorGrid) -> npt.NDArray[np.complex128]:
    """`V_int` evaluated on the tensor grid, flattened (C order)."""
    return potential_nd(tgrid, interaction_2d)


def build_h2d(tgrid: TensorGrid) -> sp.csr_matrix:
    """`H_2D` on `tgrid` (axis 0 = electronic r, axis 1 = nuclear R)."""
    return hamiltonian_nd(tgrid, [1.0, MU], potential_2d)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest projects/n2_2d_cross_section/test_hamiltonian2d.py -q`
Expected: PASS (6 tests)

- [ ] **Step 5: Commit**

```bash
git add projects/n2_2d_cross_section
git commit -m "$(cat <<'EOF'
feat(n2-2d): the exact 2-D potential surface and Hamiltonian

V(r,R) = v0(R) + l(l+1)/2r^2 - lambda(R) exp(-alpha_c r^2), verbatim
eMoScat Potentials2d.cpp:18, assembled with qscat.dvr.hamiltonian_nd from
the model functions already verified in sub-project #2.

Keeps potential_2d (the full surface, for H) and interaction_2d (V_int
alone, the perturbation) as separate named functions. v0(R) and the
centrifugal term are CHANNEL potentials that survive as r -> infinity;
sweeping them into the driving term would give a wrong T-matrix that still
looks plausible, so the split is enforced by name and by test.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3 (CRUX): Driven solve, T-matrix, cross section

**Files:**
- Create: `projects/n2_2d_cross_section/cross_section_2d.py`
- Test: `projects/n2_2d_cross_section/test_cross_section_2d.py`

**Interfaces:**
- Consumes: Tasks 1–2; `qscat.linalg.{SparseLU, c_product}`; `qscat.dvr.TensorGrid`; `projects.n2_ti_cross_section.vibrational.vibrational_states`.
- Produces:
  - `channel_vector(tgrid, k, chi_v, *, l=ELL) -> NDArray[complex128]` — masked DVR coefficients of `F_{E,l}(r)·χ_v(R)`
  - `ve_cross_section_2d(tgrid, eps, chi, v_init, vprimes, E, *, ordering="COLAMD", lam_scale=1.0, return_wavefunction=False) -> NDArray[float64]` (or `(sigma, psi_plus)`)

**Background you need — the method:**

```
E_tot = E + eps[v_init];   k = sqrt(2E)
Psi_i = mask * outer( F_{E,l}(r)*sqrt(w_r),  chi[v_init] )
rhs   = V_int * Psi_i
Psi_sc = SparseLU(E_tot*I - H2D).solve(rhs)
Psi_plus = Psi_i + Psi_sc
for each v':  k' = sqrt(2*(E_tot - eps[v'])),  skip if E_tot - eps[v'] <= 0
   Phi_f = mask * outer( F_{E',l}(r)*sqrt(w_r), chi[v'] )
   T = c_product(Phi_f, V_int * Psi_plus)
   sigma = 4 pi^3 |T|^2 / k^2
```

`chi` from `vibrational_states` is **already** a DVR coefficient vector (no extra `sqrt(w)`); `F` is a **function** and needs `·sqrt(w_r)` from `TensorGrid.sqrt_weights()[0]`. Mixing these up rescales σ. Renormalize `chi` by the c-product (`chi ← chi / sqrt(c_product(chi, chi))`) rather than relying on the Hermitian norm `eigen()` returns.

`lam_scale` multiplies `V_int` only — it exists for the two validation tests below and must NOT be used to tune physics.

**The two tests that make this task trustworthy without any reference data:**
- **V1, free-particle limit:** `lam_scale=0` ⇒ `V_int ≡ 0` ⇒ `rhs = 0` ⇒ `T = 0` ⇒ `σ = 0` exactly.
- **First Born limit (the strong one):** as `lam_scale → 0`, `Psi_sc → 0` so `T → c_product(Phi_f, V_int·Psi_i)`, the first Born amplitude, and `σ ∝ lam_scale²`. This tests the normalization, the masking, the coefficient convention and the T-matrix **all at once**, against a quantity computable directly. It is the single most valuable check in the sub-project.

- [ ] **Step 1: Write the failing tests — `projects/n2_2d_cross_section/test_cross_section_2d.py`**

```python
"""The exact 2-D driven-equation VE cross section (sub-project #6, crux).

Validated WITHOUT reference data: the free-particle limit and the first Born
limit together pin the normalization, the ECS masking, the DVR coefficient
convention and the T-matrix.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.dvr import TensorGrid
from qscat.linalg import c_product

from projects.n2_2d_cross_section.cross_section_2d import (
    channel_vector,
    ve_cross_section_2d,
)
from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import MU, interaction_diag
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states

# Small but physically sane: the interaction lives at r < ~3 bohr, so a modest
# box still supports a meaningful (if unconverged) T-matrix. Task 4 converges it.
TG = TensorGrid(
    [
        n2_electronic_grid(r_max=16.0, order=7, n_complex=5),
        n2_nuclear_grid(quadrature=10, r_max=22.0, n_complex=5),
    ]
)
EPS, CHI = vibrational_states(TG.grids[1], MU, 4)


def test_free_particle_limit_gives_exactly_zero() -> None:
    """lam_scale=0 removes the perturbation, so there is nothing to scatter off."""
    sigma = ve_cross_section_2d(TG, EPS, CHI, 0, [0, 1, 2], 0.2, lam_scale=0.0)
    assert np.all(sigma == 0.0)


def test_weak_coupling_matches_first_born() -> None:
    """As lam_scale -> 0, T -> <Phi_f|V_int|Psi_i>: the first Born amplitude.

    Pins normalization + masking + coefficient convention + T-matrix at once.
    """
    scale = 1e-4
    E, vp = 0.2, 1
    sigma = ve_cross_section_2d(TG, EPS, CHI, 0, [vp], E, lam_scale=scale)

    # First Born, computed directly here (independent of the solver's internals)
    k = np.sqrt(2.0 * E)
    e_tot = E + EPS[0]
    kp = np.sqrt(2.0 * (e_tot - EPS[vp]))
    psi_i = channel_vector(TG, k, CHI[0])
    phi_f = channel_vector(TG, kp, CHI[vp])
    v_int = scale * interaction_diag(TG)
    t_born = c_product(phi_f, v_int * psi_i)
    sigma_born = 4.0 * np.pi**3 * abs(t_born) ** 2 / (2.0 * E)

    assert sigma[0] == pytest.approx(sigma_born, rel=1e-3)


def test_sigma_scales_as_lambda_squared_in_the_born_regime() -> None:
    a = ve_cross_section_2d(TG, EPS, CHI, 0, [1], 0.2, lam_scale=1e-4)[0]
    b = ve_cross_section_2d(TG, EPS, CHI, 0, [1], 0.2, lam_scale=2e-4)[0]
    assert b / a == pytest.approx(4.0, rel=1e-3)


def test_sigma_is_real_and_non_negative() -> None:
    sigma = ve_cross_section_2d(TG, EPS, CHI, 0, [0, 1, 2, 3], 0.2)
    assert sigma.dtype == np.float64
    assert np.all(sigma >= 0.0)


def test_closed_channels_are_zero() -> None:
    """At E below a channel's threshold that channel cannot be populated."""
    e_small = 0.005
    sigma = ve_cross_section_2d(TG, EPS, CHI, 0, [0, 1, 2, 3], e_small)
    open_ = (e_small + EPS[0] - EPS) > 0.0
    assert np.all(sigma[~open_[:4]] == 0.0)


def test_channel_vector_is_masked_to_the_unscaled_region() -> None:
    """A channel projection on the complex-scaled tail is meaningless."""
    psi = channel_vector(TG, 0.6, CHI[0])
    assert np.all(psi[~TG.real_mask()] == 0.0)
    assert np.abs(psi[TG.real_mask()]).max() > 0.0


def test_array_of_energies_matches_scalar_calls() -> None:
    energies = [0.1, 0.2]
    both = ve_cross_section_2d(TG, EPS, CHI, 0, [1], energies)
    assert both.shape == (2, 1)
    for i, e in enumerate(energies):
        assert both[i, 0] == pytest.approx(
            ve_cross_section_2d(TG, EPS, CHI, 0, [1], e)[0], rel=1e-12
        )
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest projects/n2_2d_cross_section/test_cross_section_2d.py -q`
Expected: FAIL — `ModuleNotFoundError: ... cross_section_2d`

- [ ] **Step 3: Implement `projects/n2_2d_cross_section/cross_section_2d.py`**

```python
"""Exact 2-D VE cross section by the driven Lippmann-Schwinger equation.

    Psi_i    = F_{E,l}(r) chi_v(R)                     [masked to unscaled region]
    Psi_sc   = (E_tot - H_2D)^{-1} V_int Psi_i         [one sparse LU per energy]
    Psi^(+)  = Psi_i + Psi_sc
    T_{v->v'} = <chi_v' F_{E',l} | V_int | Psi^(+)>    [c-product, masked]
    sigma     = 4 pi^3 |T|^2 / k^2                     [bohr^2]

Conventions that must not be gotten wrong (each has bitten this repo or the
reference implementation):

- `chi` from `vibrational_states` is ALREADY a DVR coefficient vector; `F` is
  a FUNCTION and must be converted with `c_j = F(r_j) sqrt(w_j)` using the
  bridge-summed complex weight (`TensorGrid.sqrt_weights()`). Mixing the two
  rescales every cross section.
- Everything is paired with the C-PRODUCT (no conjugate): under ECS `H = H^T`,
  not `H^dagger`. With both sides in coefficient form the c-product IS the
  quadrature integral.
- `Psi_i` and every `Phi_f` are masked to the unscaled region. eMoScat uses a
  Hermitian dot here and is saved only by doing the same masking.

Elastic and inelastic share one formula: with `S = 1 - 2 pi i T`,
`|S - 1|^2 = 4 pi^2 |T|^2`, so Houfek's `pi |S-1|^2 / k^2` and our
`4 pi^3 |T|^2 / k^2` are the same expression. Unlike the 1-D LCP model, this
elastic T-matrix DOES contain the non-resonant background scattering.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from qscat.dvr import TensorGrid
from qscat.linalg import SparseLU, c_product

from projects.n2_2d_cross_section.channels import riccati_bessel_en
from projects.n2_2d_cross_section.hamiltonian2d import ELL, build_h2d, interaction_diag

__all__ = ["channel_vector", "ve_cross_section_2d"]


def channel_vector(
    tgrid: TensorGrid,
    k: float,
    chi_v: npt.NDArray[np.complex128],
    *,
    l: int = ELL,
) -> npt.NDArray[np.complex128]:
    """DVR coefficients of `F_{E,l}(r) chi_v(R)`, masked to the unscaled region.

    `chi_v` is already a coefficient vector; `F` is a function and picks up
    `sqrt(w_r)`.
    """
    g_r = tgrid.grids[0]
    f_vals = riccati_bessel_en(g_r.real_points, k, l)
    # sqrt_weights() is per-axis and broadcast-shaped ((n_r, 1) at D=2); ravel
    # it to pair elementwise with the 1-D electronic function values.
    sqrt_w_r = tgrid.sqrt_weights()[0].ravel()
    f_coeff = f_vals * sqrt_w_r

    chi = np.asarray(chi_v, dtype=np.complex128)
    chi = chi / np.sqrt(c_product(chi, chi))  # c-product normalization, not Hermitian

    psi = tgrid.outer([f_coeff, chi])
    psi[~tgrid.real_mask()] = 0.0
    return psi


def _sigma_at_one_energy(
    tgrid: TensorGrid,
    lu: SparseLU,
    v_diag: npt.NDArray[np.complex128],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float,
    *,
    want_psi: bool,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128] | None]:
    sigma = np.zeros(len(vprimes), dtype=np.float64)
    if E <= 0.0:
        return sigma, None

    e_tot = E + eps[v_init]
    k = float(np.sqrt(2.0 * E))

    psi_i = channel_vector(tgrid, k, chi[v_init])
    psi_plus = psi_i + lu.solve(v_diag * psi_i)
    v_psi = v_diag * psi_plus

    for j, vp in enumerate(vprimes):
        excess = e_tot - eps[vp]
        if excess <= 0.0:
            continue  # closed channel
        kp = float(np.sqrt(2.0 * excess))
        phi_f = channel_vector(tgrid, kp, chi[vp])
        t = c_product(phi_f, v_psi)
        sigma[j] = 4.0 * np.pi**3 * abs(t) ** 2 / (2.0 * E)

    return sigma, (psi_plus if want_psi else None)


def ve_cross_section_2d(
    tgrid: TensorGrid,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ordering: str = "COLAMD",
    lam_scale: float = 1.0,
    return_wavefunction: bool = False,
):
    """sigma_{v_init->v'}(E) in bohr^2, exact 2-D driven-equation solution.

    `E` may be scalar or an array; scalar returns shape `(len(vprimes),)`,
    array returns `(len(E), len(vprimes))`. One sparse LU per energy is
    reused across all `vprimes`.

    `lam_scale` scales `V_int` ONLY, for the free-particle and first-Born
    validation limits. It is a test lever, never a physics knob.
    """
    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    H = build_h2d(tgrid)
    v_diag = lam_scale * interaction_diag(tgrid)
    ident = sp.identity(tgrid.size, format="csc", dtype=np.complex128)

    out = []
    psis = []
    for e in e_arr:
        e_tot = float(e) + eps[v_init]
        lu = SparseLU((e_tot * ident - H).tocsc(), ordering=ordering)
        s, psi = _sigma_at_one_energy(
            tgrid, lu, v_diag, eps, chi, v_init, vprimes, float(e),
            want_psi=return_wavefunction,
        )
        out.append(s)
        psis.append(psi)

    sigma = np.stack(out)
    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    result = np.asarray(sigma[0], dtype=np.float64) if scalar else sigma
    if return_wavefunction:
        return result, (psis[0] if scalar else psis)
    return result
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest projects/n2_2d_cross_section/test_cross_section_2d.py -q -v`
Expected: PASS (7 tests).

If the first-Born test fails, debug **in this order** — do not loosen the tolerance:
1. **The `sqrt(w)` factor** — applied to `F` only, never to `chi`.
2. **The mask** — applied to both `Psi_i` and every `Phi_f`.
3. **`V_int` contents** — must exclude `v0(R)` and the centrifugal term.
4. **`k` vs `k'`** — `k` from the *incident* energy in the prefactor; `k'` in the exit channel function.
5. **c-product vs `vdot`** — no conjugation anywhere.
6. If `lam_scale=1e-4` still shows second-Born contamination, go smaller (1e-5) before suspecting a bug.

If it cannot be made to pass, report **BLOCKED** with the computed and Born values — do not proceed to Task 4 on an unvalidated solver.

- [ ] **Step 5: Commit**

```bash
git add projects/n2_2d_cross_section
git commit -m "$(cat <<'EOF'
feat(n2-2d): exact driven-equation VE cross section (crux)

Psi_sc = (E_tot - H_2D)^{-1} V_int Psi_i by sparse LU, one factorization per
energy reused across all final channels; T by c-product projection onto
masked channel functions; sigma = 4 pi^3 |T|^2 / k^2.

Validated with NO reference data: the free-particle limit (lam_scale=0 gives
exactly zero) and the first Born limit (as lam_scale -> 0, T approaches
<Phi_f|V_int|Psi_i> and sigma scales as lam^2). The Born check pins the
energy normalization, the ECS masking, the DVR coefficient convention and
the T-matrix simultaneously, against a quantity computed independently of
the solver.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Convergence study — choose the working grid

**Files:**
- Create: `projects/n2_2d_cross_section/convergence.py` (runnable: `python -m projects.n2_2d_cross_section.convergence`)
- Test: `projects/n2_2d_cross_section/test_convergence.py`

**Interfaces:**
- Consumes: Tasks 1–3.
- Produces:
  - `WORKING_GRID: dict[str, float | int]` — the chosen, documented grid parameters
  - `working_tgrid() -> TensorGrid`
  - `convergence_table(...) -> list[dict]` — the study results

**Background you need.** eMoScat asserts 35° and a 98-bohr electronic box without documenting the study that justified them. This task does that study and picks the **smallest** grid where σ is stable to ~1%, which is what makes the later comparison honest — and probably what makes the harness affordable.

Vary, one at a time, at a fixed anchor (**E = 0.2 Ha, v' = 1**):
- electronic `r_max` ∈ {16, 22, 30, 45}
- **ECS angle** ∈ {25, 30, 35, 40} — **the sharpest check**: a converged ECS result must not move when the contour rotates
- electronic `order` ∈ {7, 8, 9}; `n_complex` ∈ {5, 8, 11}
- nuclear `r_max` ∈ {20, 30, 40}, `quadrature` ∈ {10, 12, 14}

Record `N`, wall time, `lu.fill_factor`, and σ for each. **Also compare `ordering="COLAMD"` vs `"MMD_AT_PLUS_A"`** on the real Hamiltonian and use whichever is cheaper (a small random matrix suggested MMD roughly halves fill — confirm or refute on the real problem; it cannot affect correctness).

- [ ] **Step 1: Write `convergence.py`.** Structure it exactly like this; only `WORKING_GRID`'s values are yours to determine (Step 3).

```python
"""Convergence study for the exact 2-D N2 solver, and the working grid it picks.

eMoScat asserts a 35-degree ECS angle and a 98-bohr electronic box without
documenting the study that justified them. This module redoes that study and
records the numbers, so the grid used for the benchmark is EARNED rather than
inherited.

Run: `uv run python -m projects.n2_2d_cross_section.convergence`
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.cross_section_2d import ve_cross_section_2d
from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states

__all__ = ["BASELINE", "WORKING_GRID", "build_tgrid", "working_tgrid", "convergence_table"]

# The anchor the study is run at: the resonance region, a well-behaved VE channel.
STUDY_E, STUDY_VP, N_VIB = 0.2, 1, 4

BASELINE: dict[str, Any] = {
    "r_max": 30.0, "angle_deg": 35.0, "order": 8, "n_complex": 8,
    "nuc_r_max": 40.0, "nuc_quadrature": 14, "nuc_n_complex": 10,
}

# Filled in at Step 3 from the measured table below. Each value must carry a
# comment citing the number that justifies it.
WORKING_GRID: dict[str, Any] = dict(BASELINE)

SWEEPS: dict[str, list[Any]] = {
    "r_max": [16.0, 22.0, 30.0, 45.0],
    "angle_deg": [25.0, 30.0, 35.0, 40.0],   # the sharpest check
    "order": [7, 8, 9],
    "n_complex": [5, 8, 11],
    "nuc_r_max": [20.0, 30.0, 40.0],
    "nuc_quadrature": [10, 12, 14],
    "nuc_n_complex": [5, 10, 13],
}


def build_tgrid(params: dict[str, Any]) -> TensorGrid:
    """Both grids share `angle_deg` -- one ECS contour angle for the problem."""
    return TensorGrid(
        [
            n2_electronic_grid(
                r_max=params["r_max"], angle_deg=params["angle_deg"],
                order=params["order"], n_complex=params["n_complex"],
            ),
            n2_nuclear_grid(
                r_max=params["nuc_r_max"], angle_deg=params["angle_deg"],
                quadrature=params["nuc_quadrature"], n_complex=params["nuc_n_complex"],
            ),
        ]
    )


def working_tgrid() -> TensorGrid:
    return build_tgrid(WORKING_GRID)


def _one(params: dict[str, Any], ordering: str = "COLAMD") -> dict[str, Any]:
    t0 = time.perf_counter()
    tg = build_tgrid(params)
    eps, chi = vibrational_states(tg.grids[1], MU, N_VIB)
    sigma = float(
        ve_cross_section_2d(tg, eps, chi, 0, [STUDY_VP], STUDY_E, ordering=ordering)[0]
    )
    return {"N": tg.size, "sigma": sigma, "seconds": time.perf_counter() - t0}


def convergence_table() -> list[dict[str, Any]]:
    """Vary ONE axis at a time about BASELINE; report sigma and its drift."""
    rows: list[dict[str, Any]] = []
    for key, values in SWEEPS.items():
        prev: float | None = None
        for v in values:
            params = {**BASELINE, key: v}
            row = {"axis": key, "value": v, **_one(params)}
            row["pct_change"] = (
                None if prev is None else 100.0 * abs(row["sigma"] - prev) / abs(prev)
            )
            prev = row["sigma"]
            rows.append(row)
    return rows


def main() -> None:
    rows = convergence_table()
    print(f"| axis | value | N | sigma (bohr^2) | % change | s |")
    print(f"|---|---|---|---|---|---|")
    for r in rows:
        pc = "-" if r["pct_change"] is None else f"{r['pct_change']:.2f}%"
        print(
            f"| {r['axis']} | {r['value']} | {r['N']} | "
            f"{r['sigma']:.6e} | {pc} | {r['seconds']:.1f} |"
        )
    # Ordering comparison on the REAL Hamiltonian (a small random matrix
    # suggested MMD_AT_PLUS_A roughly halves fill; confirm or refute here).
    for ordering in ("COLAMD", "MMD_AT_PLUS_A"):
        r = _one(BASELINE, ordering=ordering)
        print(f"ordering={ordering:<14} N={r['N']} sigma={r['sigma']:.6e} {r['seconds']:.1f}s")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the study**

Run: `uv run python -m projects.n2_2d_cross_section.convergence 2>&1 | tee /tmp/conv.txt`
Record the full table — it is a deliverable, and it goes into the docs in Task 7.

- [ ] **Step 3: Choose `WORKING_GRID`** — the smallest parameter set where σ is stable to ~1% across *every* axis, and in particular is θ-independent. Write it into `convergence.py` as a module constant with a comment citing the measured numbers that justify each choice. Expose `working_tgrid()`.

- [ ] **Step 4: Write the regression test — `test_convergence.py`**

```python
"""The chosen working grid must actually be converged (V3)."""

from __future__ import annotations

import numpy as np
import pytest

from projects.n2_2d_cross_section.convergence import WORKING_GRID, working_tgrid
from projects.n2_2d_cross_section.cross_section_2d import ve_cross_section_2d
from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states
from qscat.dvr import TensorGrid


def _sigma(tgrid, E=0.2, vp=1) -> float:
    eps, chi = vibrational_states(tgrid.grids[1], MU, 4)
    return float(ve_cross_section_2d(tgrid, eps, chi, 0, [vp], E)[0])


@pytest.mark.slow
def test_working_grid_is_theta_independent() -> None:
    """THE decisive ECS check: bound-state/scattering results must not move
    when the complex-scaling contour rotates. Fill in the measured tolerance."""
    base = dict(WORKING_GRID)
    sigmas = []
    for theta in (30.0, 35.0, 40.0):
        params = {**base, "angle_deg": theta}
        tg = TensorGrid(
            [
                n2_electronic_grid(
                    r_max=params["r_max"], angle_deg=theta,
                    order=params["order"], n_complex=params["n_complex"],
                ),
                n2_nuclear_grid(
                    quadrature=params["nuc_quadrature"], r_max=params["nuc_r_max"],
                    n_complex=params["nuc_n_complex"], angle_deg=theta,
                ),
            ]
        )
        sigmas.append(_sigma(tg))
    spread = (max(sigmas) - min(sigmas)) / np.mean(sigmas)
    assert spread < 0.01, f"theta-dependence {spread:.2%}: grid is NOT converged"


@pytest.mark.slow
def test_working_grid_is_stable_under_refinement() -> None:
    """Refining past the working grid must not move sigma by more than ~1%."""
    coarse = _sigma(working_tgrid())
    fine_tg = TensorGrid(
        [
            n2_electronic_grid(
                r_max=WORKING_GRID["r_max"] * 1.5,
                angle_deg=WORKING_GRID["angle_deg"],
                order=WORKING_GRID["order"] + 1,
                n_complex=WORKING_GRID["n_complex"] + 3,
            ),
            n2_nuclear_grid(
                quadrature=WORKING_GRID["nuc_quadrature"] + 2,
                r_max=WORKING_GRID["nuc_r_max"] + 10.0,
                n_complex=WORKING_GRID["nuc_n_complex"] + 3,
            ),
        ]
    )
    assert abs(_sigma(fine_tg) - coarse) / coarse < 0.01
```

`n2_nuclear_grid` already exposes exactly the keywords used above — its signature is
`n2_nuclear_grid(*, angle_deg=35.0, r_max=40.0, n_complex=10, quadrature=14)`. Do not change
its defaults; sub-projects #3/#4 depend on them.

- [ ] **Step 5: Run, then tighten** both tolerances to just above the measured spread, recording the measured values in comments.

Run: `uv run pytest projects/n2_2d_cross_section/test_convergence.py -q -m slow`

- [ ] **Step 6: Commit** — include the convergence table in the commit message body.

---

### Task 5: The six anchors, LCP head-to-head, and the two NOTEs

**Files:**
- Create: `validation/n2/exact2d.py`
- Test: `projects/n2_2d_cross_section/test_anchors.py`

**Interfaces:**
- Consumes: Tasks 1–4; `validation.n2.reference.{ANCHOR_COORDS, anchors, ANCHOR_FACTOR, ANCHOR_MARGIN_HA}`; `validation.n2.cross_section.compute_anchor_results` (the LCP results).
- Produces: `Exact2dResult` dataclass and `compute_exact2d_results() -> list[Exact2dResult]`, `lru_cache`d (one factorization per distinct energy, reused across channels).

**Background you need.** `validation/n2/cross_section.py` is the model to follow: it computes σ at `reference.ANCHOR_COORDS`, classifies each anchor generally (never by hardcoding which coordinate is which) and caches the expensive setup. Mirror that structure. **`validation/` may import `projects/`** (the reverse is forbidden).

Report three ratios per anchor: `σ_exact/σ_houfek` (V4, the gate), `σ_LCP/σ_exact` (V5, **the primary deliverable**), and `σ_LCP/σ_houfek` (already known, for context). Group the anchors by energy so each distinct energy costs one factorization.

**On tolerances:** V4's bound must be **derived from the converged result and stated with its reasoning**, not picked in advance. If the exact solver agrees with Houfek *worse* than the LCP does, that is a finding to investigate and report — say so plainly rather than widening the bound.

- [ ] **Step 1: Implement `validation/n2/exact2d.py`** following `validation/n2/cross_section.py`'s structure (dataclass, `lru_cache`d setup, general classification, no hardcoded coordinates). The result type:

```python
@dataclass(frozen=True)
class Exact2dResult:
    """One anchor, compared three ways.

    The exact 2-D solver is the ORACLE here; `ratio_lcp_vs_exact` is the
    scientific deliverable, and `ratio_exact_vs_houfek` is the gate that
    certifies the oracle.
    """

    energy_ha: float
    channel: int                 # v' (0 = elastic)
    sigma_exact: float           # bohr^2, this sub-project
    sigma_lcp: float             # bohr^2, sub-project #3
    sigma_houfek: float          # bohr^2, CSVE.V00.J00
    ratio_exact_vs_houfek: float # V4 -- the GATE
    ratio_lcp_vs_exact: float    # V5 -- the DELIVERABLE
    ratio_lcp_vs_houfek: float   # context, already known
    gated: bool                  # same classification rule as C5
    mechanism: str               # empty if gated; else the known LCP limitation
```

Group the anchor coordinates by energy so each distinct energy costs exactly one
factorization, reused across that energy's channels — that is what makes 6 anchors cost 3
solves. Reuse `validation.n2.reference.anchors()` for the Houfek values and
`validation.n2.cross_section.compute_anchor_results()` for the LCP ones rather than
recomputing either.
- [ ] **Step 2: Run it and record every number.**
- [ ] **Step 3: Write `test_anchors.py`** asserting: σ real and ≥ 0; the four GATED anchors agree with Houfek within the tolerance you derived (documented in a comment with the measured ratios); and the exact result is at least as close to Houfek as the LCP at those anchors — **or, if it is not, an xfail with a written explanation of what was investigated.** Do not silently weaken the claim.
- [ ] **Step 4: Record V6** — what the exact model gives at the elastic (0.2, v'=0) and near-threshold (0.02, v'=1) anchors versus the LCP. This is a measurement; report it either way.
- [ ] **Step 5: Commit** with the full anchor table in the message body.

---

### Task 6: Nuclear dynamics — LCP vs exact (V7)

**Files:**
- Create: `projects/n2_2d_cross_section/nuclear_density.py`
- Test: `projects/n2_2d_cross_section/test_nuclear_density.py`

**Interfaces:**
- Consumes: Task 3's `return_wavefunction=True`; the LCP driven solution from `projects.n2_ti_cross_section.cross_section`.
- Produces: `nuclear_density(tgrid, psi) -> (R_real, density)`; `compare_to_lcp(E, v_init) -> dict` with centroid and width in `R` for both.

**Background you need.** An integrated cross section averages away much of what distinguishes an exact treatment from an approximate one. The nuclear-coordinate density of the driven solution is where a *local* width approximation is most likely to visibly break.

Project `|Ψ⁽⁺⁾(r,R)|²` over the electronic coordinate, restricted to the **unscaled** region in both coordinates (the ECS tail carries outgoing flux, not probability density), giving `ρ(R)`. Compare its shape against the LCP solver's 1-D driven solution `|ξ(R)|²` at the same energy — normalize both to unit area first, since they are not the same object dimensionally. Report centroid `⟨R⟩` and RMS width.

**This is exploratory: there is no pass/fail.** The tests should assert only that the machinery is sound (density real, non-negative, normalizable, supported where the molecule actually is), and the *comparison* is reported as data.

- [ ] **Step 1: Write the tests** — density real/non-negative; integrates to a finite value; peaks in the physically sensible range (`R ~ 1.5–3` bohr); restricted to the unscaled region.
- [ ] **Step 2: Implement `nuclear_density.py`.**
- [ ] **Step 3: Run `compare_to_lcp` at the resonance anchor (E = 0.2 Ha)** and record centroids, widths and the qualitative shape difference.
- [ ] **Step 4: Save a plot** to `docs/physics/figures/n2-2d-nuclear-density.png` (matplotlib is already a dev dependency; create the directory).
- [ ] **Step 5: Commit** with the measured comparison in the message body.

---

### Task 7: Harness group E and documentation

**Files:**
- Modify: `validation/n2/experiment.py`, `CLAUDE.md`
- Create: `docs/physics/n2-2d-cross-section.md`

- [ ] **Step 1: Wire group E** into `validation/n2/experiment.py`, following the existing C5/D1 pattern exactly: compute the exact-2D anchors, emit PASS/FAIL for gated ones and NOTE for documented-limited ones, and **guard the whole computation in try/except so a solver error becomes a labeled FAIL rather than crashing the table**. Preserve the existing 19 PASS / 0 PENDING / 2 NOTE / 0 FAIL; if the exact model closes either NOTE, that row may become PASS — report it prominently.
- [ ] **Step 2: Decide the harness grid.** If the working grid runs fast enough, group E uses it directly. If not, use a reduced grid with a **documented, looser** tolerance and say so in the row detail. Time it and decide from the measurement; state the decision in the docs.
- [ ] **Step 3: Write `docs/physics/n2-2d-cross-section.md`** covering: the driven-equation method and the exact formulas; why `V_int` excludes the channel potentials; the ECS masking and DVR coefficient conventions; the free-particle and first-Born validations; **the full convergence table** and the chosen working grid with its justification; the anchor comparison (exact vs Houfek, LCP vs exact); what happened to the two documented LCP limitations; the nuclear-density comparison with its figure; and the measured cost including the ordering choice. Be explicit that the model is a **given testbed**, that Houfek's data is the gate rather than the goal, and that the LCP-vs-exact discrepancy is the scientific result.
- [ ] **Step 4: Update `CLAUDE.md`** — add `projects/n2_2d_cross_section` to the repo map and mention group E.
- [ ] **Step 5: Full verification.**

```bash
uv run pytest projects/n2_2d_cross_section libs/qscat validation/n2 -q
uv run pytest -q -m "not slow"
uv run mypy libs/qscat
uv run ruff check .
uv run python -m validation.n2.experiment
docker/build.sh test
```
Expected: all pass; mypy 0; ruff clean; harness exit 0 with no regression; docker green.

- [ ] **Step 6: Commit.**

---

## Final verification

- [ ] Free-particle and first-Born limits pass — the solver is validated independently of any reference data.
- [ ] Convergence study complete, with a θ-independence check, and `WORKING_GRID` justified by measured numbers rather than inherited from eMoScat.
- [ ] The six anchors computed; V4 tolerance **derived and documented**, not assumed.
- [ ] LCP-vs-exact reported at every anchor — the primary deliverable.
- [ ] The two documented LCP limitations measured against the exact result, reported either way.
- [ ] Nuclear-density comparison recorded with a figure.
- [ ] Harness group E wired and guarded; existing 19 PASS / 0 PENDING / 2 NOTE / 0 FAIL not regressed; docker green.
- [ ] No model parameter was tuned to improve agreement with anything.
- [ ] Nothing new promoted into `qscat`; no `projects/` → `validation/` import.
