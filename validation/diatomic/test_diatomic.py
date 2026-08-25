"""Tests for the NO/F2 exact-2D TI cross-section oracle (the model port).

There is no independent golden data for NO/F2 (only N2 has Houfek's data), so
these gate the exact solver's own well-posedness (real, finite, non-negative σ
of the right shape) plus the physical fact that each molecule has a low-lying
resonance producing a strong vibrational-excitation cross section -- computed
entirely through the promoted `qscat.core` + `qscat.model` library. The full
dense σ(E) curves + figures are generated from config through
`apps/qscat-run` (committed under `docs/physics/figures/`), not re-run here.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pytest
from qscat.core.driven import ve_cross_section
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid

from validation.diatomic.config import CONFIGS, MoleculeConfig

# The exact-2D VE channels [elastic, v'=1, v'=2] the dense curves (now produced
# by `apps/qscat-run`, e.g. `qscat-run run` on an F2/NO VE config) sweep.
VPRIMES = [0, 1, 2]


def _small_tgrid() -> TensorGrid:
    """Tiny, deliberately unconverged grid -- fast shape/finiteness checks only."""
    return TensorGrid(
        [
            electronic_grid(r_max=14.0, order=6, n_complex=4),
            nuclear_grid(r_max=20.0, quadrature=8, n_complex=4),
        ]
    )


def test_no_and_f2_are_configured() -> None:
    assert set(CONFIGS) == {"NO", "F2"}
    assert CONFIGS["NO"].model.ell == 1 and CONFIGS["F2"].model.ell == 1
    assert CONFIGS["NO"].model.mu == 13614.16
    assert CONFIGS["F2"].model.mu == 17315.99


@pytest.mark.parametrize("name", ["NO", "F2"])
def test_ti_curve_shape_and_physical(name: str) -> None:
    """On a small grid at a few energies: σ real, finite, ≥0, right shape."""
    cfg: MoleculeConfig = CONFIGS[name]
    tg = _small_tgrid()
    eps, chi = vibrational_states(tg.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)
    E = np.array([0.02, 0.04, 0.06])
    sigma = ve_cross_section(tg, cfg.model, eps, chi, 0, VPRIMES, E)
    assert sigma.shape == (len(E), len(VPRIMES))
    assert sigma.dtype == np.float64
    assert np.all(np.isfinite(sigma))
    assert np.all(sigma >= 0.0)
    # neutral vibrational spacing is the small, molecule-specific value
    assert 0.001 < (eps[1] - eps[0]) < 0.02


# The resonance window is SCANNED, not sampled at a few hand-picked energies:
# the claim under test is that a strong v'=1 excitation exists *somewhere in the
# window*, and three points cannot establish that. Windows, steps and floors are
# all calibrated against measurement on the fine `da_grid` (2026-08-17).
#
# The two molecules need DIFFERENT sampling, which is the reason this is a table
# rather than one shared rule:
#
#   F2 -- BROAD peak. A dense reference scan (97 energies, 0.004-0.100 Ha at
#         0.001) puts the v'=1 global maximum at 0.2340 bohr^2 at E=0.031.
#         Re-sampling that scan at every phase offset of coarser grids recovers,
#         worst case over phases: 0.999 of the peak at step 0.002, 0.996 at
#         0.004, 0.991 at 0.006. A 9-point scan therefore cannot miss it.
#
#   NO -- SHARP boomerang structure, the opposite regime. Measured at step 0.002
#         its v'=1 swings between adjacent points: 16.98 (0.027), 34.76 (0.029),
#         6.13 (0.031), 3.56 (0.033), 22.98 (0.035), 29.30 (0.037), 1.96 (0.039).
#         A 5.7x drop in one 0.002 Ha step. Coarsening to 0.004 makes the
#         recovered maximum swing 23.67-34.76 depending on phase, so NO keeps the
#         finer step. (This test finds a MAXIMUM; it does not claim to resolve
#         NO's curve, which 0.002 Ha does not.)
#
# The floors sit well below each measured peak (F2 0.15 vs 0.2340; NO 10.0 vs
# 34.76, and vs 23.67 even at the worst coarse phase), leaving margin for
# cross-architecture BLAS differences -- tiny here: MUMPS/x86_64 and
# SuperLU/arm64 agree on F2's value to ~11 significant figures.
#
# This test previously asserted a SINGLE shared bound, `max > 0.3`, from two
# probe energies. That bound was written on 2026-07-27 against the shared
# "working grid"; the 2026-08-15 curve-driver retirement repointed the test at
# the finer per-molecule `da_grid` without recalibrating it. On that deck F2
# never reaches 0.3 anywhere, so the test could not pass -- while for NO the
# same bound sat 100x below the measured value and tested nothing. It went
# unnoticed because the test is @slow and CI runs `-m "not slow"`; only a Docker
# `test` image build exercises it, and that build had been failing since.
_RESONANCE_WINDOW = {  # (E_min, E_max, step) in Hartree
    "F2": (0.010, 0.044, 0.004),
    "NO": (0.017, 0.045, 0.002),
}
_V1_FLOOR = {"F2": 0.15, "NO": 10.0}


def _window_energies(name: str) -> npt.NDArray[np.float64]:
    lo, hi, step = _RESONANCE_WINDOW[name]
    return np.round(np.arange(lo, hi + 0.5 * step, step), 6)


@pytest.mark.slow
@pytest.mark.parametrize("name", ["NO", "F2"])
def test_low_energy_resonance_drives_vibrational_excitation(name: str) -> None:
    """On the converged deck (the fine `da_grid`, the same grid `apps/qscat-run`
    uses for these molecules): each molecule's low-lying resonance produces a
    strong v'=1 excitation cross section somewhere in its resonance window --
    the physical signature the exact-2D solver must reproduce for the port to be
    meaningful.

    The window is SCANNED (see `_RESONANCE_WINDOW`), because "somewhere in the
    window" is the claim; a few hand-picked energies would assume the answer.
    Sampling and floors are per molecule and measured -- their cross sections
    differ by two orders of magnitude and their peaks by an order in width.
    Still not the full curve: the dense curves and figures come from
    `apps/qscat-run`.
    """
    cfg = CONFIGS[name]
    tg = cfg.da_grid()
    eps, chi = vibrational_states(tg.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)
    E = _window_energies(name)
    sigma = ve_cross_section(tg, cfg.model, eps, chi, 0, VPRIMES, E)
    assert np.all(sigma >= 0.0) and np.all(np.isfinite(sigma))
    peak_i = int(np.argmax(sigma[:, 1]))
    peak = float(sigma[peak_i, 1])
    assert peak > _V1_FLOOR[name], (
        f"{name}: strongest v'=1 excitation in the scanned window "
        f"[{E[0]:.3f}, {E[-1]:.3f}] Ha ({E.size} energies) was {peak:.4f} bohr^2 "
        f"at E={E[peak_i]:.3f}, below the resonant-excitation floor "
        f"{_V1_FLOOR[name]}"
    )
