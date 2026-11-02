"""Saved own-K phase reconnaissance for a pole design, never a pole measurement."""

from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path

import numpy as np
from qscat.units import hartree_to_ev

from projects.ukrmol_co.verify_stretched_inner import numbers, records, symmetric
from projects.ukrmol_co.verify_stretched_preparation import stream_files
from projects.ukrmol_co.verify_stretched_properties import require
from projects.ukrmol_co.verify_stretched_real_completion import verify


def phase_seed(energy_rydberg, matrices):
    energy = np.asarray(energy_rydberg, float)
    matrices = np.asarray(matrices, float)
    require(
        len(energy) >= 5 and np.all(np.diff(energy) > 0) and np.isfinite(matrices).all(),
        "Finite increasing saved phase grid required",
    )
    phase = np.unwrap(2 * np.arctan(np.linalg.eigvalsh(matrices)).sum(axis=1)) / 2
    gradient = np.gradient(phase, energy)
    index = int(np.argmax(gradient))
    require(0 < index < len(energy) - 1 and gradient[index] > 0, "No interior rising phase feature")
    width = 2 / gradient[index]
    return dict(
        index=index,
        energy_rydberg=float(energy[index]),
        phase_sum_rad=float(phase[index]),
        gradient_per_rydberg=float(gradient[index]),
        diagnostic_full_width_rydberg=float(width),
        diagnostic_seed_incident_hartree=[float(energy[index] / 2), float(-width / 4)],
        position_ev=float(hartree_to_ev(energy[index] / 2)),
        diagnostic_full_width_ev=float(hartree_to_ev(width / 2)),
        maximum_adjacent_phase_change_rad=float(abs(np.diff(phase)).max()),
        grid_energy_rydberg=energy.tolist(),
        phase_sum_rad_on_grid=phase.tolist(),
    )


def review(archive, collection):
    completion = verify(archive, collection)
    names = [n for n in collection["files"] if Path(n).name.endswith(("fort.921", "fort.10"))]
    files, _, _ = stream_files(archive, collection, capture_names=names)
    result = {}
    for case, date in [("dz", "20261031"), ("tz", "20261101")]:
        prefix = f"prepared/stretched-{case}-inner-source-{date}/outputs/inner/native/B1/"
        km = records(files[collection["reference_files"][prefix + "fort.921"]])
        energies = []
        matrices = []
        for row in km[5:]:
            values = numbers(row[24:], "<f8", 22)
            energies.append(values[0])
            matrices.append(symmetric(values[1:], 6))
        require(len(energies) == 99, "Own native99-point grid required")
        channel = records(files[collection["reference_files"][prefix + "fort.10"]])
        labels = [list(struct.unpack("<4q", row[:32])) for row in channel[68:74]]
        thresholds = [struct.unpack_from("<d", row, 32)[0] for row in channel[68:]]
        result[case] = phase_seed(energies, matrices) | dict(
            open_labels=labels,
            first_closed_threshold_rydberg=min(t for t in thresholds if t > max(energies)),
        )
    return dict(
        real_prerequisite_verdict=completion["combined_20degree_real_prerequisite_verdict"],
        cases=result,
        new_ode_maps=0,
        native_calls=0,
        pole_measurement=False,
        basis_pole_verdict=None,
        protected_route_holdout_s_matrices_analyzed=0,
        scope=(
            "Own saved native99-point phase gradients for seeds only;single-pole width heuristic "
            "with background/sampling uncertainty,not outgoing poles or a basis-sensitivity verdict"
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("collection", type=Path)
    args = parser.parse_args()
    print(
        json.dumps(
            review(args.archive, json.loads(args.collection.read_text())), indent=2, allow_nan=False
        )
    )
