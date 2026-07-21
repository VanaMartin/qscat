"""Electronic ECS resonance pole for the N2 harness Group B check (B1).

Builds the two-angle FEM-DVR-ECS N2 electronic grids and computes E_res(R0)
via the promoted `qscat.ecs.find_resonance_pole` matcher, using this
package's own `model.v_eff_el` (the validated closed-form "single source"
potential -- see `model.py`).

Grid construction mirrors `projects/n2_resonance/grid_n2.n2_electronic_grid`
(the toy-model version, validated in `projects/n2_resonance/test_pole.py`
V1/V2), reimplemented in-line here rather than imported: `validation/` stays
self-contained and does not import from `projects/` (see CLAUDE.md's repo
map -- `projects/` is per-problem toy-model work, `validation/` benchmarks
against it independently).
"""

from __future__ import annotations

import model
import numpy as np
import numpy.typing as npt
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec, eigen, hamiltonian
from qscat.ecs import find_resonance_pole

R0: float = model.PARAMS["potential"]["R_0"]  # N2 equilibrium bond length (Bohr)

# Re in [0.04, 0.16] Ha, Im in [-0.05, 0] Ha -- brackets the R0 pole (~0.09 Ha
# region); see projects/n2_resonance/test_pole.py's identical WINDOW.
WINDOW: tuple[float, float, float, float] = (0.04, 0.16, -0.05, 0.0)


def n2_electronic_grid(
    angle_deg: float,
    *,
    r_pivot: float = 10.0,
    n_real: int = 8,
    r_max: float = 30.0,
    n_complex: int = 8,
    quadrature: int = 8,
) -> FemDvrEcsGrid:
    """Build the N2 electronic FEM-DVR-ECS grid: real region + ECS tail.

    `n_real` uniform real elements tile [0, r_pivot]; `n_complex` uniform
    complex elements at `angle_deg` tile [r_pivot, r_max]. Mirrors
    `projects/n2_resonance/grid_n2.n2_electronic_grid`.
    """
    real_length = r_pivot / n_real
    complex_length = (r_max - r_pivot) / n_complex

    elements = [ElementSpec(real_length) for _ in range(n_real)]
    elements += [ElementSpec(complex_length, angle_deg) for _ in range(n_complex)]

    spec = GridSpec(quadrature=quadrature, elements=elements, x_min=0.0)
    return FemDvrEcsGrid(spec)


def electronic_hamiltonian(R: float, grid: FemDvrEcsGrid) -> npt.NDArray[np.complex128]:
    """Assemble H_el(R) = kinetic(grid, 1.0) + diag(model.v_eff_el(grid.points, R))."""
    return hamiltonian(grid, lambda z: model.v_eff_el(z, R), 1.0)


def find_pole_at_R(
    R: float,
    angle_a: float = 35.0,
    angle_b: float = 44.0,
    window: tuple[float, float, float, float] = WINDOW,
) -> tuple[complex, float]:
    """Two-angle ECS resonance pole at fixed bond length `R` (Hartree).

    Diagonalizes `H_el(R)` on grids built at `angle_a`/`angle_b` and matches
    the angle-stable pole via `qscat.ecs.find_resonance_pole`. Returns
    `(E_pole, residual)`: `E_pole.real` is `E_res`, and
    `max(0, -2*E_pole.imag)` is `Gamma`. Mirrors the validated
    `projects/n2_resonance/pole.find_pole` (Task 2).
    """
    grid_a = n2_electronic_grid(angle_a)
    grid_b = n2_electronic_grid(angle_b)
    Ea, _ = eigen(electronic_hamiltonian(R, grid_a))
    Eb, _ = eigen(electronic_hamiltonian(R, grid_b))
    return find_resonance_pole(Ea, Eb, window)


def e_res_at_R0() -> tuple[complex, float]:
    """`(E_pole, residual)` at the N2 equilibrium bond length `R0`."""
    return find_pole_at_R(R0)
