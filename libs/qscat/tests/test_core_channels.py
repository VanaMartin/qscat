"""Gates for `qscat.core.channels.channel_vector`'s charge-aware dispatch
(sub-project #D, Task 3): `charge=0` (the default) must stay the exact,
fast `riccati_bessel_en` path used by every existing VE/DA caller;
`charge != 0` routes the electronic factor through `coulomb_f_en` instead
(the H2+/DR Coulomb-incident-wave case).
"""

from __future__ import annotations

import numpy as np
from qscat.core.channels import channel_vector
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.dvr import TensorGrid


def test_channel_vector_charge_zero_unchanged() -> None:
    # charge=0 must be identical to the default (no Coulomb path)
    tg = TensorGrid([electronic_grid(r_max=14.0, order=6, n_complex=4),
                     nuclear_grid(r_max=20.0, n_complex=4, quadrature=8)])
    chi = np.zeros(tg.grids[1].n, dtype=np.complex128)
    chi[0] = 1.0
    a = channel_vector(tg, 0.6, chi, 1)
    b = channel_vector(tg, 0.6, chi, 1, charge=0)
    assert np.array_equal(a, b)


def test_channel_vector_coulomb_is_finite_and_differs() -> None:
    tg = TensorGrid([electronic_grid(r_max=14.0, order=6, n_complex=4),
                     nuclear_grid(r_max=20.0, n_complex=4, quadrature=8)])
    chi = np.zeros(tg.grids[1].n, dtype=np.complex128)
    chi[0] = 1.0
    free = channel_vector(tg, 0.6, chi, 1)
    coul = channel_vector(tg, 0.6, chi, 1, charge=-1)
    assert np.all(np.isfinite(coul)) and not np.allclose(coul, free)
