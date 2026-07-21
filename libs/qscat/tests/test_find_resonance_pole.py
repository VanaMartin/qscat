"""Tests for the promoted general resonance-pole matcher (`qscat.ecs`).

Synthetic (no Hamiltonian involved): two eigenvalue sets that share exactly
one value ("pole", angle-stable) and otherwise differ ("continuum", rotates
with the angle) -- the matcher must pick out the shared value with a
near-zero residual, mirroring the physical N2 resonance search this was
promoted from (`projects/n2_resonance/pole.find_pole`,
`projects/n2_resonance/test_pole.py`).
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.ecs import find_resonance_pole


def test_shared_eigenvalue_is_matched_as_the_pole() -> None:
    # The shared pole value sits at index 1 in eigs_a but index 2 in eigs_b
    # (a different index in each array), so a naive same-index comparator
    # (rather than a true all-pairs |ea - eb| search) would fail this test.
    pole_val = 0.10 - 0.02j
    eigs_a = np.array([0.05 + 0.09j, pole_val, 0.14 - 0.30j])
    eigs_b = np.array([0.02 - 0.31j, 0.13 + 0.11j, pole_val])
    window = (0.0, 0.2, -0.05, 0.05)

    E_pole, residual = find_resonance_pole(eigs_a, eigs_b, window)

    assert abs(E_pole - pole_val) < 1e-12
    assert residual < 1e-12


def test_empty_window_raises_value_error() -> None:
    eigs_a = np.array([0.1 + 0.0j, 0.2 - 0.1j])
    eigs_b = np.array([0.1 + 0.0j, 0.2 - 0.1j])
    window = (5.0, 6.0, -1.0, 1.0)  # deliberately misses every eigenvalue

    with pytest.raises(ValueError):
        find_resonance_pole(eigs_a, eigs_b, window)
