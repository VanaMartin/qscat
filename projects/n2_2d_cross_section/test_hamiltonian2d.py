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
    """Axis 0 is the ELECTRON (mass 1), axis 1 the nuclei (mass mu).

    Structural check, not an incidental-scale one: build the right and
    swapped kinetic operators explicitly from `kinetic_sparse` per-grid
    matrices and `kron_sum`, and assert `kinetic_nd(tg, [1.0, MU])` matches
    the explicit correct pairing to round-off while differing hugely from
    the explicit swapped pairing. This exercises `kinetic_nd`'s own
    axis-to-mass wiring rather than restating it, and is independent of
    which grid's finest FEM element happens to be numerically stiffer.

    Also ties the check to `build_h2d` itself (not just `kinetic_nd` in
    isolation): its kinetic part is recovered by subtracting the
    mass-independent potential diagonal from `H_2D`, and must equal the
    explicit correct pairing too -- so a mass argument swapped inside
    `build_h2d` is caught here, not just a swap inside `kinetic_nd`.
    """
    tg = _small_tgrid()
    from qscat.dvr import kinetic_nd, kinetic_sparse, potential_nd
    from qscat.linalg import kron_sum

    g_electronic, g_nuclear = tg.grids

    got = kinetic_nd(tg, [1.0, MU])

    right = kron_sum([kinetic_sparse(g_electronic, 1.0), kinetic_sparse(g_nuclear, MU)])
    wrong = kron_sum([kinetic_sparse(g_electronic, MU), kinetic_sparse(g_nuclear, 1.0)])

    assert abs(got - right).max() < 1e-9 * abs(right).max()

    # Wide, clearly-reported margin: `got` disagrees with the swapped pairing
    # by a large fraction of the swapped pairing's own scale -- not merely
    # more than round-off, but grossly different.
    diff_from_wrong = abs(got - wrong).max()
    assert diff_from_wrong > 0.5 * abs(wrong).max()

    H = build_h2d(tg)
    v_diag = potential_nd(tg, potential_2d)
    kinetic_from_h2d = H - sp.diags(v_diag, format="csr")
    assert abs(kinetic_from_h2d - right).max() < 1e-9 * abs(right).max()
