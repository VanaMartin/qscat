"""Differential gate: `qscat.model.library.N2` against a direct, independent
closed-form oracle -- the Morse/sigmoid/Gaussian expressions inlined here
from the eMoScat deck constants (`validation/n2/config.json`'s `potential`
block), NOT `N2`'s own methods and NOT `projects.n2_2d_cross_section
.hamiltonian2d` (which, after that module was rewired to delegate to
`qscat.model.N2`, would make such a comparison tautological -- N2 vs
itself). This test suite is self-contained: it imports nothing from
`projects`.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from qscat.dvr import TensorGrid, hamiltonian_nd, potential_nd
from qscat.model import ResonanceModel
from qscat.model.library import F2, N2, NO

# N2 deck constants (`validation/n2/config.json`), transcribed independently
# of `qscat.model.library.N2` for a direct-formula comparison.
_MU = 12766.36
_ELL = 2
_D0 = 0.75102
_ALPHA0 = 1.1535
_R0 = 2.01943
_LAMBDA_INF = 6.21066
_LAMBDA_1 = 1.05708
_R_LAMBDA = -27.9833
_LAMBDA_C = 5.38022
_R_C = 2.405
_ALPHA_C = 0.4

# Sample grids for the pointwise (non-tensor-grid) comparisons: real values
# spanning the well and the tail, deliberately NOT the DVR grid points.
_R_SAMPLE = np.array([1.5, 2.01943, 2.5, 3.5, 6.0], dtype=complex)
_R_SAMPLE_ROW = _R_SAMPLE.reshape(1, -1)
_r_SAMPLE = np.array([0.1, 0.5, 1.0, 1.7, 4.0], dtype=complex).reshape(-1, 1)


def _oracle_v0(R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
    """Direct Morse formula: `D0*(exp(-2*a*(R-R0)) - 2*exp(-a*(R-R0)))`."""
    Rc = np.asarray(R, dtype=np.complex128)
    a = _ALPHA0
    out = _D0 * (np.exp(-2 * a * (Rc - _R0)) - 2 * np.exp(-a * (Rc - _R0)))
    return np.asarray(out, dtype=np.complex128)


def _oracle_lam(R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
    """Direct sigmoid formula for `lambda(R)`, independent of `DiatomicResonanceModel.lam`."""
    Rc = np.asarray(R, dtype=np.complex128)
    lam0 = (_LAMBDA_C - _LAMBDA_INF) * (1 + np.exp(_LAMBDA_1 * (_R_C - _R_LAMBDA)))
    out = _LAMBDA_INF + lam0 / (1 + np.exp(_LAMBDA_1 * (Rc - _R_LAMBDA)))
    return np.asarray(out, dtype=np.complex128)


def _oracle_v_int(r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
    """Direct Gaussian-in-r interaction: `-lambda(R) * exp(-alpha_c * r**2)`."""
    rr = np.asarray(r, dtype=np.complex128)
    out = -_oracle_lam(R) * np.exp(-_ALPHA_C * rr**2)
    return np.asarray(out, dtype=np.complex128)


def _oracle_surface(r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
    """Direct full surface: `v0(R) + ell(ell+1)/(2 r^2) + v_int(r, R)`."""
    rr = np.asarray(r, dtype=np.complex128)
    out = _oracle_v0(R) + _ELL * (_ELL + 1) / (2.0 * rr**2) + _oracle_v_int(rr, R)
    return np.asarray(out, dtype=np.complex128)


def _small_tgrid() -> TensorGrid:
    """Deliberately tiny: tests STRUCTURE (does `N2.hamiltonian` assemble the
    surface correctly onto a `TensorGrid`?), not converged physics. Built
    from `qscat.core.grids.electronic_grid`/`nuclear_grid` (the promoted,
    project-independent grid builders), mirroring the shapes used by
    `n2_2d_cross_section/test_hamiltonian2d.py`'s `_small_tgrid` without
    importing the project shims that wrap them.
    """
    import qscat.core.grids as _grids

    electronic = _grids.electronic_grid(r_max=14.0, order=6, n_complex=4)
    nuclear = _grids.nuclear_grid(quadrature=8, r_max=20.0, n_complex=4)
    return TensorGrid([electronic, nuclear])


def test_n2_registry_constants_match_the_deck() -> None:
    assert N2.mu == _MU
    assert N2.ell == _ELL


def test_n2_v0_matches_the_direct_morse_formula() -> None:
    got = np.asarray(N2.v0(_R_SAMPLE))
    want = _oracle_v0(_R_SAMPLE)
    assert np.abs(got - want).max() < 1e-14


def test_n2_lam_matches_the_direct_sigmoid_formula() -> None:
    got = np.asarray(N2.lam(_R_SAMPLE))
    want = _oracle_lam(_R_SAMPLE)
    assert np.abs(got - want).max() < 1e-14
    # lam(R_c) == lambda_c, per the sigmoid's own construction.
    assert abs(complex(N2.lam(_R_C)) - _LAMBDA_C) < 1e-12


def test_n2_v_int_matches_the_direct_gaussian_formula() -> None:
    got = np.asarray(N2.v_int(_r_SAMPLE, _R_SAMPLE_ROW))
    want = _oracle_v_int(_r_SAMPLE, _R_SAMPLE_ROW)
    assert np.abs(got - want).max() < 1e-14


def test_n2_surface_matches_the_direct_closed_form() -> None:
    got = np.asarray(N2.surface(_r_SAMPLE, _R_SAMPLE_ROW))
    want = _oracle_surface(_r_SAMPLE, _R_SAMPLE_ROW)
    assert np.abs(got - want).max() < 1e-14


def test_n2_hamiltonian_matches_hamiltonian_nd_assembled_from_the_surface() -> None:
    """Self-consistent (built from the SAME `N2.surface`), gates the assembly
    wiring (`hamiltonian_nd(tgrid, [1, mu], surface)`), not the physics of
    the surface itself -- that's covered by the direct-formula tests above.
    The harness (`validation.n2.experiment`) gates the full end-to-end
    physics against Houfek's independent data.
    """
    tg = _small_tgrid()
    got = N2.hamiltonian(tg)
    want = hamiltonian_nd(tg, [1.0, N2.mu], N2.surface)
    assert isinstance(got, sp.csr_matrix)
    assert got.shape == want.shape
    assert abs(got - want).max() < 1e-12 * abs(want).max()


def test_n2_interaction_diag_matches_potential_nd_assembled_from_v_int() -> None:
    tg = _small_tgrid()
    got = N2.interaction_diag(tg)
    want = potential_nd(tg, N2.v_int)
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
