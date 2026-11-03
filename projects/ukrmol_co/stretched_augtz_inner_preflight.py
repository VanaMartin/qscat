"""Saved-target threshold screen for the unchanged99-point/six-open-channel export."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

from qscat.units import HARTREE_TO_EV

from projects.ukrmol_co.phase_diagnostics import NATIVE_RYDBERG_PER_INPUT_EV
from projects.ukrmol_co.verify_stretched_properties import SECTORS, native_energies, require

DECK = "projects/ukrmol_co/compressed-cas12-inner-inputs-20261026/scattering.rsolve.doublet.B1.inp"


def threshold_screen(energies, deck):
    """No Hamiltonian/native action; do not silently narrow an incompatible legacy grid."""
    clean = "\n".join(s.split("!", 1)[0] for s in deck.splitlines())

    def values(name):
        found = re.findall(rf"\b{name}\s*=\s*([^\n/]+)", clean, re.I)
        require(len(found) == 1, "Unique native energy assignment " + name)
        return [float(x.strip()) for x in found[0].split(",") if x.strip()]

    require(
        values("nescat") == [99.0]
        and values("nerang") == [1.0]
        and values("ieunit") == [2.0]
        and values("einc") == [0.1, 0.05],
        "Only the unchanged legacy99-point grid is screened",
    )
    require(
        set(energies) == set(SECTORS) and all(len(v) == 8 for v in energies.values()),
        "Retained64 target inventory",
    )
    require(
        all(math.isfinite(e) for vs in energies.values() for e in vs),
        "Finite native target energies",
    )
    ground = energies["singlet.A1"][0]
    require(
        ground == min(e for vs in energies.values() for e in vs), "Own singlet-A1 ground threshold"
    )
    requested = [0.1 + 0.05 * i for i in range(99)]
    grid = [e * NATIVE_RYDBERG_PER_INPUT_EV / 2 for e in requested]
    targets = sorted(
        (e - ground, sector, i)
        for sector, vs in energies.items()
        for i, e in enumerate(vs, 1)
        if (sector, i) != ("singlet.A1", 1)
    )
    threshold = targets[0][0]
    crossed = [
        dict(sector=s, root=i, threshold_hartree=e, threshold_ev=e * HARTREE_TO_EV)
        for e, s, i in targets
        if e <= grid[-1]
    ]
    rows = [i for i, e in enumerate(grid) if e >= threshold]
    return dict(
        unchanged_six_open_grid_prerequisite_verdict=not crossed,
        first_excited_threshold_hartree=threshold,
        first_excited_threshold_ev=threshold * HARTREE_TO_EV,
        native_grid_maximum_hartree=grid[-1],
        native_grid_maximum_physical_ev=grid[-1] * HARTREE_TO_EV,
        crossed_retained_target_thresholds=crossed,
        affected_grid_indices_zero_based=rows,
        affected_requested_energies_ev=[round(requested[i], 12) for i in rows],
        new_native_calls=0,
        new_ci_calls=0,
        new_hamiltonian_actions=0,
        new_maps=0,
        route_holdout_s_matrices_analyzed=0,
        automatic_successor=False,
        expensive_inner_execution_released=False,
        full_sweep_released=False,
        scope="Saved N-electron thresholds versus existing diagnostic prescription;"
        "not a native variable-open-channel export or a rejection of the electronic model",
    )


def screen_saved(cidata, deck):
    """Decode saved native energies and identify the exact screened input bytes."""
    energies, norm_error = native_energies(cidata)
    return threshold_screen(energies, deck.decode()) | dict(
        cidata_sha256=hashlib.sha256(cidata).hexdigest(),
        deck_sha256=hashlib.sha256(deck).hexdigest(),
        maximum_native_norm_error=norm_error,
        native_rydberg_per_requested_ev=NATIVE_RYDBERG_PER_INPUT_EV,
        hartree_to_reported_ev=HARTREE_TO_EV,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cidata", type=Path)
    p.add_argument("--deck", type=Path, default=Path(DECK))
    a = p.parse_args()
    print(
        json.dumps(
            screen_saved(a.cidata.read_bytes(), a.deck.read_bytes()), indent=2, allow_nan=False
        )
    )
