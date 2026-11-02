"""Authenticate the stopped guard and diagnose saved blocks/rates, with no new maps."""

from __future__ import annotations

import argparse
import io
import json
import struct
from pathlib import Path

import numpy as np

from projects.ukrmol_co.verify_stretched_inner import records
from projects.ukrmol_co.verify_stretched_preparation import stream_files
from projects.ukrmol_co.verify_stretched_properties import files_in, require, sha
from projects.ukrmol_co.verify_stretched_real_guard import verify


def rate_separation(energy, threshold, angle_degrees, length):
    """Asymptotic unwanted closed-dual decay relative to the desired open dual."""
    require(
        np.isfinite([energy, threshold, angle_degrees, length]).all()
        and 0 < energy < threshold
        and 0 < angle_degrees < 90
        and length > 0,
        "An open ground channel and closed coupled channel are required",
    )
    k = np.sqrt(energy)
    kappa = np.sqrt(threshold - energy)
    angle = np.deg2rad(angle_degrees)
    separation = kappa * np.cos(angle) - k * np.sin(angle)
    return dict(
        open_momentum=float(k),
        closed_decay_momentum=float(kappa),
        critical_angle_degrees=float(np.rad2deg(np.arctan(kappa / k))),
        relative_decay_rate_per_bohr=float(separation),
        finite_interval_damping_estimate=float(np.exp(-length * separation)),
    )


def saved_blocks(q, y, a):
    wall, _, rank, singular = np.linalg.lstsq(a.conj().T, a.T, rcond=None)
    require(rank == 6, "Saved open amplitudes are rank deficient")
    flux = (y - y.conj().T) / (2j) - a.conj().T @ a
    wronskian = a.conj().T @ wall - a.T
    maximum_flux = np.unravel_index(abs(flux).argmax(), flux.shape)
    maximum_wronskian = np.unravel_index(abs(wronskian).argmax(), wronskian.shape)
    return dict(
        amplitude_condition=float(singular[0] / singular[-1]),
        matching_condition=float(np.linalg.cond(np.eye(len(q)) - q @ y)),
        amplitude_maximum=float(abs(a).max()),
        open_closed_flux_defect=float(abs(flux[:6, 6:]).max()),
        open_open_flux_defect=float(abs(flux[:6, :6]).max()),
        closed_closed_flux_defect=float(abs(flux[6:, 6:]).max()),
        maximum_flux_indices_zero_based=[int(i) for i in maximum_flux],
        maximum_wronskian_indices_zero_based=[int(i) for i in maximum_wronskian],
        wronskian_residual=float(abs(wronskian).max()),
        wall_unitarity_defect=float(abs(wall.conj().T @ wall - np.eye(6)).max()),
    )


def diagnose(archive, collection, allocation):
    original = verify(archive, collection, allocation)
    require(
        original["original_exits"] == [2] * 6 and original["real_guard_verdict"] is None,
        "Diagnosis is only for the original stopped allocation",
    )
    prefix = "prepared/stretched-dz-inner-source-20261031/outputs/inner/native/B1/"
    channel_name = collection["reference_files"][prefix + "fort.10"]
    names = [channel_name, *[f"outputs/check/dz/real-{i}.npz" for i in [60, 80]]]
    files, hashes, _ = stream_files(archive, collection, capture_names=names)
    source = files_in(files["source.tar.gz"])
    channels = records(files[channel_name])
    labels = [struct.unpack("<4q", row[:32]) for row in channels[68:]]
    thresholds = [struct.unpack_from("<d", row, 32)[0] for row in channels[68:]]
    guard = json.loads(files["outputs/check/dz/guard.json"])
    rows = []
    for saved in guard["real_rows"]:
        index = saved["index"]
        name = f"outputs/check/dz/real-{index}.npz"
        with np.load(io.BytesIO(files[name]), allow_pickle=False) as z:
            q, y, a = [z[n] for n in ["q", "y", "a"]]
            energy = float(z["energy_rydberg"])
        blocks = saved_blocks(q, y, a)
        closed_index = blocks["maximum_wronskian_indices_zero_based"][0]
        require(closed_index >= 6, "Residual maximum is not on a closed-channel column")
        target = labels[closed_index][1]
        target_label = struct.unpack("<4q", channels[4 + target - 1][:32])
        diagnostic = saved["diagnostics"]
        rows.append(
            dict(
                index=index,
                energy_rydberg=energy,
                **blocks,
                closed_channel_label=list(labels[closed_index]),
                closed_target_label=list(target_label),
                closed_threshold_rydberg=thresholds[closed_index],
                rate=rate_separation(
                    energy,
                    thresholds[closed_index],
                    diagnostic["angle"],
                    diagnostic["endpoint"] - 18,
                ),
                saved_matrix_sha256=hashes[name],
            )
        )
    observed_ratio = rows[1]["wronskian_residual"] / rows[0]["wronskian_residual"]
    damping_ratio = (
        rows[1]["rate"]["finite_interval_damping_estimate"]
        / rows[0]["rate"]["finite_interval_damping_estimate"]
    )
    return dict(
        original_scientific_exits=original["original_exits"],
        rows=rows,
        observed_wronskian_residual_ratio=float(observed_ratio),
        asymptotic_relative_damping_ratio=float(damping_ratio),
        evidence_archive_sha256=collection["archive"]["sha256"],
        executed_outgoing_map_sha256=sha(source["projects/ukrmol_co/outgoing_map.py"]),
        scope=(
            "Saved matrices and asymptotic rate arithmetic only. Rate separation "
            "supports finite-tail closed-dual contamination; it is not a bound, "
            "a new integration, unique causal proof or a repaired guard."
        ),
        new_maps=0,
        native_calls=0,
        pole_searches=0,
        tz_guards_started=False,
        original_guard_requalified=False,
        full_sweep_released=False,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ["archive", "collection", "allocation"]:
        p.add_argument(name, type=Path)
    args = p.parse_args()
    print(
        json.dumps(
            diagnose(
                args.archive,
                json.loads(args.collection.read_text()),
                json.loads(args.allocation.read_text()),
            ),
            indent=2,
            allow_nan=False,
        )
    )
