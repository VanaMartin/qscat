"""Load the declared finite boundary problem without opening a saved K/S grid."""

from pathlib import Path

import numpy as np

from projects.ukrmol_co.saved_matrices import array, read_interface, records


def load_boundary_model(source: Path) -> dict:
    """Retain all 298 channels; accept only the neutral 48-target/2048-root view."""
    return _load_boundary_model(source, minimum_excited_threshold_rydberg=0.3675)


def load_closed_boundary_model(source: Path) -> dict:
    """Load the same finite view for strictly negative incident Hartree energies.

    All six ground channels and positive excited thresholds remain present.
    This loader does not qualify a positive-energy grid or evaluate a bound
    state; the decaying solver must still enforce its negative-energy domain.
    """
    return _load_boundary_model(source, minimum_excited_threshold_rydberg=0.0)


def _load_boundary_model(source: Path, *, minimum_excited_threshold_rydberg: float) -> dict:
    interface = read_interface(source / "fort.10", source / "fort.21")
    if (
        interface["target_components"] != 48
        or interface["boundary_channels"] != 298
        or interface["boundary_poles"] != 2048
        or interface["multipole_orders"] != 2
        or interface["symmetry_charge"][3] != 0
        or abs(interface["radius_bohr"] - 18) > 1e-12
    ):
        raise ValueError("The declared neutral48/298/2048/18-bohr boundary view is required")
    labels = np.asarray(interface["channel_labels"])
    targets = np.asarray(interface["target_energies_hartree"])
    thresholds = np.asarray(interface["channel_thresholds_rydberg"])
    if not np.isfinite(targets).all() or not np.isfinite(thresholds).all():
        raise ValueError("Finite native target energies and channel thresholds required")
    np.testing.assert_allclose(
        thresholds, 2 * (targets[labels[:, 1] - 1] - targets[0]), atol=1e-10, rtol=0
    )
    if np.any(thresholds[:6] != 0) or np.any(thresholds[6:] <= minimum_excited_threshold_rydberg):
        domain = (
            "all-closed boundary"
            if minimum_excited_threshold_rydberg == 0
            else "six-open-channel grid"
        )
        raise ValueError(f"The declared {domain} domain differs")
    angular = labels[:, 2]
    raw = records(source / "fort.21")
    count = len(thresholds)
    cf = np.zeros((2, count, count))
    packed = array(raw[5], "<f8", count * (count + 1)).reshape(2, -1)
    indices = np.tril_indices(count)
    for matrix, values in zip(cf, packed, strict=True):
        matrix[indices] = values
        matrix.T[indices] = values
    return dict(
        interface=interface,
        thresholds=thresholds,
        angular=angular,
        ground=targets[0],
        poles=array(raw[6], "<f8", 2048),
        amplitudes=array(raw[7], "<f8", count * 2048).reshape(2048, count),
        c2=cf[0] + np.diag(angular * (angular + 1)),
        c3=cf[1],
    )
