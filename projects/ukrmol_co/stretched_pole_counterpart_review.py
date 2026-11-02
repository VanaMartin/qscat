"""Saved DZ/TZ local-candidate comparison, not a cross-basis wavefunction matcher."""

from __future__ import annotations

import argparse
import io
import json
import struct
import tarfile
from pathlib import Path

import numpy as np
from qscat.units import hartree_to_ev

from projects.ukrmol_co.verify_stretched_dz_pole import difference
from projects.ukrmol_co.verify_stretched_dz_pole import verify as verify_dz
from projects.ukrmol_co.verify_stretched_inner import records
from projects.ukrmol_co.verify_stretched_properties import require, sha
from projects.ukrmol_co.verify_stretched_tz_pole import verify as verify_tz

READINGS = {
    "dz": "dc7819118262b9de02befe6f1766b28e4f74ec63747f711e3cbecd8230b5043a",
    "tz": "85390acc2aefa0ff1fe450e91e86447bca81d7dae1e5e180e490380a07d8f61f",
}


def compare(left, right, labels, residues):
    """Report physical shifts and phase-free diagnostics without an identity cutoff."""
    for reading in [left, right]:
        require(
            reading["original_exits"] == [0] * 6
            and reading["finite_reference_pole_verdict"] is True
            and reading["physical_pi_support_verdict"] is True,
            "Two qualified local poles",
        )
    require(labels[0] == labels[1], "Same open-channel labels necessary for comparison")
    poles = [complex(*r["roots"][-1]["z_hartree"]) for r in [left, right]]
    require(all(np.isfinite(z) and z.imag < 0 for z in poles), "Finite outgoing poles")
    widths = [float(hartree_to_ev(-2 * z.imag)) for z in poles]
    shifts = difference(poles[1], poles[0])
    budgets = dict(position_ev=0.05, full_width_ev=max(0.001, 0.05 * widths[0]))
    profiles = [np.asarray(r["character_witnesses"][1]["profile"]) for r in [left, right]]
    magnitudes = [abs(r) / np.linalg.norm(r) for r in residues]
    return dict(
        R_bohr=2.5,
        comparison="finite atomic basis plus separately optimized orbital frame",
        position_ev=[float(hartree_to_ev(z.real)) for z in poles],
        full_width_ev=widths,
        signed_tz_minus_dz_position_ev=float(hartree_to_ev((poles[1] - poles[0]).real)),
        signed_tz_minus_dz_full_width_ev=widths[1] - widths[0],
        difference=shifts,
        physical_budgets_dz_reference=budgets,
        full_width_relative_change_to_dz=shifts["full_width_ev"] / widths[0],
        full_width_relative_change_to_tz=shifts["full_width_ev"] / widths[1],
        alternate_tz_reference_width_budget_ev=max(0.001, 0.05 * widths[1]),
        candidate_budget_components={
            k: "reject" if v > budgets[k] else "pass" for k, v in shifts.items()
        },
        candidate_insensitivity_verdict=all(shifts[k] <= budgets[k] for k in shifts),
        same_model_numerical_budget_ev=0.0001,
        shifts_exceed_budgets_even_after_two_declared_numerical_budgets=all(
            shifts[k] - 2 * 0.0001 > budgets[k] for k in shifts
        ),
        open_labels=labels[0],
        profile_tv=float(abs(profiles[0] - profiles[1]).sum() / 2),
        profiles=[p.tolist() for p in profiles],
        normalized_residue_magnitude_distance=float(np.linalg.norm(magnitudes[0] - magnitudes[1])),
        residue_norm_ratio_tz_to_dz=float(
            np.linalg.norm(residues[1]) / np.linalg.norm(residues[0])
        ),
        counterpart_review="compatible local isolated ground-channel Pi-support candidates;"
        "no contradiction in these necessary observables, no unique state assignment",
        cross_basis_wavefunction_identity_verdict=None,
        global_uniqueness_or_geometry_identity_verdict=None,
        profile_or_residue_identity_threshold_applied=False,
        scientific_disposition=(
            "stop basis-insensitive/production inference;"
            "redirect to an explicit electronic-model decision"
        ),
        native_calls=0,
        new_ode_maps=0,
        route_holdout_s_matrices_analyzed=0,
        automatic_successor=False,
        full_sweep_released=False,
    )


def review(dz_bundle, tz_bundle):
    readings, labels, residues, input_hashes = [], [], [], {}
    for case, directory, verifier in [("dz", dz_bundle, verify_dz), ("tz", tz_bundle, verify_tz)]:
        raw = (directory / "reading.json").read_bytes()
        require(sha(raw) == READINGS[case], "Published independent reading " + case)
        collection = json.loads((directory / "collection.json").read_text())
        reading = verifier(directory / "evidence.tar.gz", collection)
        require(reading == json.loads(raw)["result"], "Independent reconstruction " + case)
        readings.append(reading)
        with tarfile.open(directory / "evidence.tar.gz") as archive:
            # The verifier authenticated every member above; recheck the few used here.
            def member(name, collection=collection):
                data = archive.extractfile(name).read()
                require(sha(data) == collection["files"][name]["sha256"], "Saved member " + name)
                return data

            model_file = next(n for n in collection["reference_files"] if n.endswith("/fort.10"))
            channel = records(member(collection["reference_files"][model_file]))
            labels.append([list(struct.unpack("<4q", row[:32])) for row in channel[68:74]])
            with np.load(
                io.BytesIO(member("outputs/check/moments-1.npz")), allow_pickle=False
            ) as z:
                residues.append(z["residue"])
        input_hashes[case] = dict(
            reading_sha256=sha(raw),
            archive_sha256=collection["archive"]["sha256"],
            source_commit=reading["source_commit"],
        )
    result = compare(*readings, labels, residues)
    result["input_evidence"] = input_hashes
    result["isolation_scope"] = "Unit winding and rank-one moments in each own 0.001Ha circle;"
    result["isolation_scope"] += "no pole count over the union of search boxes or global spectrum"
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dz_bundle", type=Path)
    parser.add_argument("tz_bundle", type=Path)
    args = parser.parse_args()
    print(json.dumps(review(args.dz_bundle, args.tz_bundle), indent=2, allow_nan=False))
