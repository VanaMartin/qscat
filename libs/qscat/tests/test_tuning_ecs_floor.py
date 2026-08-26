"""`max_stable_angle` must not read round-off as growth: a potential that has
already decayed to 1e-17 on the tail (O2's EMO neutral at its pivot) returned
0 degrees at every angle before the absolute floor -- measured."""

from __future__ import annotations

import numpy as np
from qscat.model import O2
from qscat.tuning import max_stable_angle


def test_round_off_tail_is_not_growth():
    rng = np.random.default_rng(0)

    def noise(z):
        # flat at round-off, with a few-percent jitter -- "growth" of 1e-19
        return 1e-17 * (1.0 + 0.05 * rng.standard_normal(np.shape(z)))

    assert max_stable_angle(noise, 12.0, 20.0, angle_cap=35.0) == 35.0


def test_o2_neutral_has_a_stable_nuclear_angle():
    assert max_stable_angle(O2.v0, 12.0, 20.0, angle_cap=35.0) == 35.0


def test_genuine_growth_is_still_rejected():
    def grows(z):
        return 1e-3 * np.exp(0.1 * (np.asarray(z) - 12.0))  # grows along any ray

    assert max_stable_angle(grows, 12.0, 20.0, angle_cap=35.0) == 0.0
