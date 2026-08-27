"""The O2 seed model: `lam(R)` is a `TailR` with the sign the physics needs."""

from __future__ import annotations

import numpy as np
from qscat.model import TailR

from validation.factory.targets.o2 import o2_seed


def test_o2_seed_lam_is_a_tail_form_with_the_right_asymptote_sign():
    m = o2_seed()
    assert isinstance(m.lam, TailR) and m.lam.q == 4
    # deeper at the equilibrium than at infinity (the anion binds MORE in the
    # molecule than O + O^- at -EA), and monotone beyond the well
    R = np.array([2.3, 4.0, 8.0, 30.0])
    lam = m.lam(R).real
    assert lam[0] > lam[-1] and np.all(np.diff(lam[1:]) < 0)
