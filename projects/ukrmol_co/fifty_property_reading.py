"""Independent saved-DENPROP reading for the explicit forty-state property export."""

from __future__ import annotations

import argparse
import json
import math
from itertools import pairwise
from pathlib import Path

IRREPS = ("A1", "B1", "B2", "A2")
SECTORS = tuple(f"{spin}.{irrep}" for irrep in IRREPS for spin in ("singlet", "triplet"))


def read_selected(text, reference, ground_dipole):
    """Check raw state/moment records without calling the producer inventory reader."""
    headers, states, moments = [], [], {}
    for line in text.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "6":
            headers.append(int(fields[4]))
        elif fields[0] == "5":
            states.append(
                (
                    int(fields[1]),
                    fields[-2].lower() + "." + fields[-1],
                    float(fields[8].replace("D", "E")),
                )
            )
        elif fields[0] == "1":
            key = tuple(int(v) for v in fields[1:8])
            value = float(fields[8].replace("D", "E"))
            if len(key) != 7 or key in moments or not math.isfinite(value):
                raise ValueError("Duplicate, incomplete or nonfinite moment")
            moments[key] = value
    if headers != [40] or [row[0] for row in states] != list(range(1, 41)):
        raise ValueError("Incomplete forty-state header/ID sequence")
    energies = [row[2] for row in states]
    if not all(math.isfinite(v) for v in energies) or any(
        b - a < -1e-10 for a, b in pairwise(energies)
    ):
        raise ValueError("Nonfinite or unordered property energies")
    if set(row[1] for row in states) != set(SECTORS) or set(reference) != set(SECTORS):
        raise ValueError("Wrong spin/irrep sectors")
    labels, mapping, errors = {}, [], []
    for sector in SECTORS:
        rows = [row for row in states if row[1] == sector]
        if len(rows) != 5 or len(reference[sector]) != 5:
            raise ValueError("Five target roots required in every sector")
        for root, (state_id, _, energy) in enumerate(rows, 1):
            want = float(reference[sector][root - 1])
            if not math.isfinite(want) or abs(energy - want) > 1e-7:
                raise ValueError("Property energy differs from independent CI")
            errors.append(abs(energy - want))
            labels[state_id] = f"{sector}.{root}"
            mapping.append(state_id)
    covered, ranks = set(), set()
    for key in moments:
        left, left_irrep, right, right_irrep, kind, rank, projection = key
        if (
            left not in labels
            or right not in labels
            or kind != 3
            or rank not in (1, 2)
            or abs(projection) > rank
        ):
            raise ValueError("Invalid state/operator/rank in moment")
        if left_irrep != IRREPS.index(labels[left].split(".")[1]) or right_irrep != IRREPS.index(
            labels[right].split(".")[1]
        ):
            raise ValueError("Moment symmetry indices differ from its state labels")
        covered.update((left, right))
        ranks.add(rank)
    if covered != set(range(1, 41)) or ranks != {1, 2}:
        raise ValueError("Incomplete state/rank moment coverage")
    ground = next(i for i, label in labels.items() if label == "singlet.A1.1")
    key = (ground, 0, ground, 0, 3, 1, 0)
    if key not in moments:
        raise ValueError("Missing ground-state dipole")
    dipole = -moments[key]
    if not math.isfinite(ground_dipole) or abs(dipole - ground_dipole) > 1e-5:
        raise ValueError("Ground-state dipole gate failed")
    return dict(
        selected_property_inventory_verified=True,
        target_count=40,
        target_map=mapping,
        state_labels_by_id={str(k): v for k, v in labels.items()},
        multipole_record_count=len(moments),
        maximum_energy_error_hartree=max(errors),
        ground_dipole_z_au=dipole,
        ground_dipole_difference_au=dipole - ground_dipole,
        physical_cross_frame_channel_correspondence_verdict=None,
        scattering_ensemble_verdict=None,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--properties", type=Path, required=True)
    parser.add_argument("--frame-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    frame = json.loads(args.frame_report.read_text())
    if frame["new_frame_import_verdict"] is not True:
        raise ValueError("New-frame qualification is missing")
    reference = {}
    for sector in SECTORS:
        spin, irrep = sector.split(".")
        count = 13 if irrep == "A1" else 8
        probe = next(
            p
            for p in frame["probes"]
            if (p["spin"], p["irrep"], p["roots"], p["space"])
            == (0 if spin == "singlet" else 2, irrep, count, 160)
        )
        if (
            probe["all_roots_pass"] is not True
            or probe["checkpoint_sha256"] != frame["checkpoint_sha256"]
        ):
            raise ValueError("Independent CI reference is unqualified or another frame")
        reference[sector] = [s["energy_hartree"] for s in probe["states"][:5]]
    result = read_selected(
        args.properties.read_text(),
        reference,
        frame["native_properties"]["ground_state_dipole_au"][2],
    )
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
