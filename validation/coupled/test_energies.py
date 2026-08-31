"""The mesh has to resolve two scales an order of magnitude apart: cusps AT
each threshold, and interference spread across the range on the scale of a
resonance width."""

from __future__ import annotations

import numpy as np

from validation.coupled.energies import (
    BACKGROUND,
    CLUSTER_STEP,
    DEDUP_TOL,
    E_HI,
    E_LO,
    sweep_energies,
    vibrational_thresholds,
)


def test_the_mesh_spans_the_declared_window() -> None:
    E = sweep_energies()
    assert E[0] >= E_LO
    assert E[-1] <= E_HI
    assert np.all(np.diff(E) > 0.0), "the mesh must be strictly ascending"


def test_no_two_energies_are_wastefully_close() -> None:
    """Every solve costs ~15 s on the production deck. Two energies a
    ten-thousandth of a mHa apart are a wasted one -- dedup must use a
    tolerance, not exact equality."""
    E = sweep_energies()
    assert float(np.min(np.diff(E))) >= DEDUP_TOL


def test_every_threshold_in_range_is_bracketed_closely() -> None:
    """A cusp is non-analytic AT the threshold, so points must sit on both
    sides of it, close."""
    E = sweep_energies()
    inside = [t for t in vibrational_thresholds() if E_LO < t < E_HI]
    assert len(inside) >= 15, f"expected many thresholds in range, got {len(inside)}"
    for t in inside:
        below = E[E <= t]
        above = E[E >= t]
        assert below.size and above.size, f"threshold {t} not bracketed"
        assert t - below[-1] <= CLUSTER_STEP * 1.01
        assert above[0] - t <= CLUSTER_STEP * 1.01


def test_the_background_is_no_coarser_than_declared() -> None:
    E = sweep_energies()
    assert float(np.max(np.diff(E))) <= BACKGROUND * 1.01


def test_the_mesh_is_the_declared_size() -> None:
    """The campaign's wall clock is per-energy times this number; a mesh
    that silently doubled would silently double the run."""
    assert 900 <= sweep_energies().size <= 1150
