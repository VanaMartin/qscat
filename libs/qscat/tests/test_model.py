"""Differential gate: `qscat.model.library.N2` against its independent
oracles -- `projects.n2_resonance.potential` (v0/v_int) and
`projects.n2_2d_cross_section.hamiltonian2d` (potential_2d/build_h2d/
interaction_diag), NOT against this module's own formula. `qscat.model` is
brand new here; these two N2 project modules are the pre-existing,
already-validated implementations it must reproduce to round-off before any
N2 project is rewired to consume it.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from qscat.dvr import TensorGrid
from qscat.model import ResonanceModel
from qscat.model.library import F2, N2, NO

from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import build_h2d, interaction_diag
from projects.n2_2d_cross_section.hamiltonian2d import potential_2d as oracle_surface
from projects.n2_resonance.potential import v0 as oracle_v0
from projects.n2_resonance.potential import v_int as oracle_v_int
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid

# Sample grids for the pointwise (non-tensor-grid) comparisons: real values
# spanning the well and the tail, deliberately NOT the DVR grid points.
_R_SAMPLE = np.array([1.5, 2.01943, 2.5, 3.5, 6.0], dtype=complex)
_R_SAMPLE_ROW = _R_SAMPLE.reshape(1, -1)
_r_SAMPLE = np.array([0.1, 0.5, 1.0, 1.7, 4.0], dtype=complex).reshape(-1, 1)


def _small_tgrid() -> TensorGrid:
    """Deliberately tiny: mirrors `n2_2d_cross_section/test_hamiltonian2d.py`'s
    `_small_tgrid` -- this tests STRUCTURE (does `N2.hamiltonian` reproduce
    `build_h2d`?), not converged physics."""
    return TensorGrid(
        [
            n2_electronic_grid(r_max=14.0, order=6, n_complex=4),
            n2_nuclear_grid(quadrature=8, r_max=20.0, n_complex=4),
        ]
    )


def test_n2_registry_constants_match_the_oracle() -> None:
    assert N2.mu == 12766.36
    assert N2.ell == 2


def test_n2_v0_matches_the_potential_oracle() -> None:
    got = np.asarray(N2.v0(_R_SAMPLE))
    want = np.asarray(oracle_v0(_R_SAMPLE))  # type: ignore[no-untyped-call]
    assert np.abs(got - want).max() < 1e-14


def test_n2_v_int_matches_the_potential_oracle() -> None:
    got = np.asarray(N2.v_int(_r_SAMPLE, _R_SAMPLE_ROW))
    want = np.asarray(oracle_v_int(_r_SAMPLE, _R_SAMPLE_ROW))  # type: ignore[no-untyped-call]
    assert np.abs(got - want).max() < 1e-14


def test_n2_surface_matches_the_hamiltonian2d_oracle() -> None:
    got = np.asarray(N2.surface(_r_SAMPLE, _R_SAMPLE_ROW))
    want = np.asarray(oracle_surface(_r_SAMPLE, _R_SAMPLE_ROW))
    assert np.abs(got - want).max() < 1e-14


def test_n2_hamiltonian_matches_build_h2d_to_round_off() -> None:
    tg = _small_tgrid()
    got = N2.hamiltonian(tg)
    want = build_h2d(tg)
    assert isinstance(got, sp.csr_matrix)
    assert got.shape == want.shape
    assert abs(got - want).max() < 1e-12 * abs(want).max()


def test_n2_interaction_diag_matches_the_oracle() -> None:
    tg = _small_tgrid()
    got = N2.interaction_diag(tg)
    want = interaction_diag(tg)
    assert np.abs(got - want).max() < 1e-12 * np.abs(want).max()


def test_n2_satisfies_the_resonance_model_protocol() -> None:
    assert isinstance(N2, ResonanceModel)


def test_no_and_f2_are_registered_with_the_deck_constants() -> None:
    """Not gated against an independent oracle (none exists for NO/F2 yet --
    see the design spec's Deliverable 2) -- just pins the transcribed deck
    values so a registry edit can't silently drift.
    """
    assert (NO.mu, NO.ell, NO.alpha_c) == (13614.16, 1, 1.0)
    assert (F2.mu, F2.ell, F2.alpha_c) == (17315.99, 1, 3.0)
    assert isinstance(NO, ResonanceModel)
    assert isinstance(F2, ResonanceModel)
