"""Tests for the N2 nuclear vibrational states (Task 1, sub-project #3).
Builds the nuclear FEM-DVR-ECS grid + neutral vibrational states and
validates the FEM-DVR eigenvalues against the ANALYTIC Morse spectrum of
eMoScat's own potential (model-consistent check -- see module note below),
not against real N2 spectroscopy.

Maintainer decision (see the model-caveat analysis): eMoScat's
neutral N2 Morse (`D_0=0.75102` Ha =~ 20.4 eV) is accepted as-is -- it is a
MODEL potential for the resonance study, not a spectroscopic fit to real N2
(whose real dissociation energy is =~9.8 eV). Its vibrational spacing
(omega_e =~ 0.0125 Ha analytic, eps1-eps0 =~ 0.0124 Ha) is therefore =~16%
larger than real N2 (0.01074 Ha / 2358 cm^-1). That gap is a property of the
*model*, not a solver bug, so the internal correctness check below validates
the FEM-DVR solver against the closed-form Morse spectrum of THIS potential.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states

REAL_N2_SPACING_HA = 0.01074  # real N2 omega_e = 2358 cm^-1 = 0.2924 eV

_CONFIG = json.loads(
    (Path(__file__).resolve().parents[2] / "validation" / "n2" / "config.json").read_text()
)
MU = _CONFIG["reduced_mass"]  # N2 nuclear reduced mass (a.u.), 12766.36
D0 = _CONFIG["potential"]["D_0"]  # 0.75102 Ha
ALPHA0 = _CONFIG["potential"]["alpha_0"]  # 1.1535 bohr^-1


def _analytic_morse_levels(n: int) -> np.ndarray:
    """Closed-form bound-state spectrum of eMoScat's Morse potential.

    `eps_v = -D_0 + omega_e*(v+1/2) - (omega_e^2/(4*D_0))*(v+1/2)^2`, with
    `omega_e = alpha_0 * sqrt(2*D_0/mu)` (atomic units, hbar=1). This is the
    textbook Morse-oscillator eigenvalue formula applied to eMoScat's own
    `(D_0, alpha_0)`, independent of the FEM-DVR-ECS solver -- it is what the
    solver's `eps` must match if the grid/kinetic/eigen machinery is correct.
    """
    omega_e = ALPHA0 * math.sqrt(2 * D0 / MU)
    v = np.arange(n)
    return -D0 + omega_e * (v + 0.5) - (omega_e**2 / (4 * D0)) * (v + 0.5) ** 2


def test_eps_real_and_ascending():
    grid = n2_nuclear_grid()
    eps, chi = vibrational_states(grid, MU, 6)
    assert eps.shape == (6,)
    assert chi.shape == (6, grid.n)
    assert np.all(np.abs(eps.imag) < 1e-6)
    assert np.all(np.diff(eps.real) > 0.0)


def test_vibrational_eigenvalues_match_analytic_morse_spectrum():
    """Solver-correctness check: FEM-DVR eps vs. the closed-form Morse
    spectrum of eMoScat's own potential (model-consistent, not real N2).
    """
    grid = n2_nuclear_grid()
    eps, _chi = vibrational_states(grid, MU, 5)
    analytic = _analytic_morse_levels(5)
    np.testing.assert_allclose(eps.real, analytic, atol=1e-5)


def test_model_spacing_is_not_real_n2_spacing():
    """Documents the known model-vs-reality gap (not a solver defect).

    eMoScat's D_0 (~20.4 eV) is ~2x real N2's dissociation energy (~9.8 eV),
    so the model's omega_e (analytic: ~0.01251 Ha) and FEM-DVR eps1-eps0
    (~0.01241 Ha) sit ~16% above real N2's omega_e (0.01074 Ha). This is an
    intentional, accepted property of the model potential -- see the module
    docstring and `docs/physics/n2-resonance.md`'s "Model caveat" section.
    """
    omega_e = ALPHA0 * math.sqrt(2 * D0 / MU)
    assert math.isclose(omega_e, 0.01251, abs_tol=1e-5)

    grid = n2_nuclear_grid()
    eps, _chi = vibrational_states(grid, MU, 2)
    spacing = eps[1].real - eps[0].real
    assert math.isclose(spacing, 0.01241, abs_tol=1e-5)
    assert (spacing - REAL_N2_SPACING_HA) / REAL_N2_SPACING_HA > 0.10


def test_low_lying_levels_roughly_evenly_spaced():
    # Anharmonicity is small for low v: eps1-eps0 and eps2-eps1 shouldn't
    # differ by more than ~15%.
    grid = n2_nuclear_grid()
    eps, _chi = vibrational_states(grid, MU, 6)
    d01 = eps[1].real - eps[0].real
    d12 = eps[2].real - eps[1].real
    assert abs(d12 - d01) / d01 < 0.15
