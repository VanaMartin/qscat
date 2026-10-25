"""Inventory saved CO boundary/channel records and reconstruct real-axis S from K."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from pathlib import Path

import numpy as np
from scipy.io import FortranEOFError, FortranFile


def records(path: Path) -> list[bytes]:
    """Read complete little-endian native records; malformed trailing bytes raise."""
    result = []
    with FortranFile(path, "r", header_dtype="<u4") as stream:
        while True:
            try:
                result.append(stream.read_record("u1").tobytes())
            except FortranEOFError:
                break
    return result


def array(record: bytes, dtype: str, count: int) -> np.ndarray:
    if len(record) != count * 8:
        raise ValueError(f"Unexpected record length {len(record)}; expected {count * 8}")
    values = np.frombuffer(record, dtype)
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite native record")
    return values


def read_k(path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Read one native open-open K dataset; validate every energy record and EOF."""
    raw = records(path)
    key, dataset, _, _, _ = array(raw[0], "<i8", 5)
    if key != 11 or dataset != 1 or len(raw[1]) != 80:
        raise ValueError("Unsupported K dataset")
    symmetry = array(raw[2][:32], "<i8", 4)
    geometry = array(raw[2][32:], "<f8", 2)
    targets, vib, diss, channels, ranges = array(raw[3], "<i8", 5)
    if vib != 0 or diss != 0 or ranges != 1:
        raise ValueError("Only a single fixed-nuclei K energy range is supported")
    index, expected = array(raw[4][:16], "<i8", 2)
    start, step = array(raw[4][16:], "<f8", 2)
    if index != 1 or expected <= 0 or len(raw) != 5 + expected:
        raise ValueError("Incomplete K energy sequence or extra dataset")
    energies, matrices = [], []
    for record in raw[5:]:
        if len(record) < 32:
            raise ValueError("Truncated K record")
        nopen, ndopen, packed_size = np.frombuffer(record, "<i8", count=3)
        if ndopen != 0 or not 0 < nopen <= channels or packed_size != nopen * (nopen + 1) // 2:
            raise ValueError("Unsupported K dimensions")
        values = array(record[24:], "<f8", int(packed_size + 1))
        matrix = np.zeros((nopen, nopen))
        offset = 1
        for column in range(nopen):
            matrix[: column + 1, column] = values[offset : offset + column + 1]
            matrix[column, : column + 1] = values[offset : offset + column + 1]
            offset += column + 1
        energies.append(values[0])
        matrices.append(matrix)
    energies = np.asarray(energies)
    if np.any(np.diff(energies) <= 0):
        raise ValueError("K energies are not strictly increasing")
    np.testing.assert_allclose(energies, start + step * np.arange(expected), atol=1e-12, rtol=0)
    return (
        energies,
        np.stack(matrices),
        {
            "target_components": int(targets),
            "boundary_channels": int(channels),
            "symmetry_charge": symmetry.tolist(),
            "geometry_mass": geometry.tolist(),
            "actual_energy_records": len(energies),
            "open_channels": len(matrices[0]),
            "stored_triangle_only": True,
        },
    )


def scattering_matrix(k: np.ndarray) -> np.ndarray:
    """Cayley transform for a finite real symmetric flux-normalized open-channel K."""
    k = np.asarray(k, dtype=float)
    if k.ndim != 2 or k.shape[0] != k.shape[1] or not np.isfinite(k).all():
        raise ValueError("A finite square K matrix is required")
    np.testing.assert_allclose(k, k.T, atol=1e-12, rtol=0)
    identity = np.eye(len(k))
    return np.linalg.solve(identity - 1j * k, identity + 1j * k)


def read_interface(channels_path: Path, boundary_path: Path, *, selected_property_ids=None) -> dict:
    """Validate native metadata; original property IDs need an explicit selection.

    The default remains strict dense-row interpretation. Explicit selection
    preserves raw labels and exposes their dense lookup separately; it neither
    edits native bytes nor qualifies another consumer of the original pointers.
    """
    ch, bd = records(channels_path), records(boundary_path)
    if array(ch[0], "<i8", 5)[0] != 10 or array(bd[0], "<i8", 5)[0] != 11:
        raise ValueError("Wrong channel/boundary dataset key")
    if not np.array_equal(array(ch[2], "<i8", 4), array(bd[2], "<i8", 4)):
        raise ValueError("Channel/boundary dimensions differ")
    targets, vib, diss, count = array(ch[2], "<i8", 4)
    if vib != 0 or diss != 0 or len(ch) != 4 + targets + count:
        raise ValueError("Unsupported/incomplete fixed-nuclei channel file")
    if ch[3] != bd[3]:
        raise ValueError("Channel/boundary geometry or symmetry differs")
    symmetry = array(ch[3][:32], "<i8", 4)
    geometry = array(ch[3][32:], "<f8", 2)
    states = ch[4 : 4 + targets]
    target_labels = [array(v[:32], "<i8", 4).tolist() for v in states]
    target_energies = np.array([array(v[32:], "<f8", 1)[0] for v in states])
    selected = None
    if selected_property_ids is not None:
        selected = np.asarray(selected_property_ids)
        if (
            selected.shape != (targets,)
            or not np.issubdtype(selected.dtype, np.integer)
            or len(np.unique(selected)) != targets
            or np.any(selected < 1)
            or np.min(selected) != 1
        ):
            raise ValueError(
                "Unique positive integer property selection including ground1 required"
            )
        selected = np.sort(selected)
        if [row[0] for row in target_labels] != list(range(1, targets + 1)):
            raise ValueError("Selected target table must have consecutive written row IDs")
    lookup = (
        {int(original): i + 1 for i, original in enumerate(selected)}
        if selected is not None
        else {i: i for i in range(1, targets + 1)}
    )
    labels, thresholds = [], []
    target_rows = []
    for i, record in enumerate(ch[4 + targets :], 1):
        label = array(record[:32], "<i8", 4)
        threshold = array(record[32:], "<f8", 1)[0]
        if label[0] != i or label[1] not in lookup or label[2] < abs(label[3]):
            raise ValueError("Invalid channel label")
        labels.append(label.tolist())
        thresholds.append(float(threshold))
        target_rows.append(lookup[int(label[1])])
    if selected is not None:
        expected = 2 * (target_energies[np.asarray(target_rows) - 1] - target_energies[0])
        if np.max(abs(np.asarray(thresholds) - expected)) > 1e-10:
            raise ValueError("Mapped channel threshold differs from the selected target energy")
    multipoles, roots, vectors, buttle = array(bd[4][:32], "<i8", 4)
    if multipoles <= 0 or roots <= 0 or vectors != 0 or buttle != 0 or len(bd) != 8:
        raise ValueError("Unsupported boundary correction/vector layout")
    cf = array(bd[5], "<f8", int(multipoles * count * (count + 1) // 2)).reshape(multipoles, -1)
    poles = array(bd[6], "<f8", int(roots))
    array(bd[7], "<f8", int(roots * count))
    if np.any(np.diff(poles) < 0):
        raise ValueError("Boundary poles are not energy ordered")
    result = {
        "target_components": int(targets),
        "boundary_channels": int(count),
        "symmetry_charge": symmetry.tolist(),
        "geometry_mass": geometry.tolist(),
        "target_labels": target_labels,
        "target_energies_hartree": target_energies.tolist(),
        "channel_labels": labels,
        "channel_thresholds_rydberg": thresholds,
        "boundary_poles": int(roots),
        "radius_bohr": float(array(bd[4][32:], "<f8", 1)[0]),
        "boundary_pole_range_hartree": [float(poles.min()), float(poles.max())],
        "multipole_orders": int(multipoles),
        "multipole_packed_entries_per_order": cf.shape[1],
        "multipole_nonzero_entries": np.count_nonzero(cf, axis=1).tolist(),
        "multipole_maximum_absolute_coefficients_native": np.max(abs(cf), axis=1).tolist(),
        "buttle_flag": int(buttle),
        "saved_full_ci_vectors": int(vectors),
        "spectral_scope": "2048 exported poles; uncorrected truncated sum, not full spectrum",
    }
    if selected is not None:
        result.update(selected_property_ids=selected.tolist(), channel_target_rows=target_rows)
    return result


def unit(path: Path, name: str) -> str:
    text = "\n".join(line.split("!", 1)[0] for line in path.read_text().splitlines())
    matches = re.findall(rf"\b{name}\s*=\s*(\d+)", text, flags=re.I)
    if len(matches) != 1:
        raise ValueError(f"Missing/ambiguous unit {name}")
    return f"fort.{int(matches[0])}"


def inspect_case(directory: Path, model: str, contract: dict) -> dict:
    energy, matrices, header = read_k(directory / unit(directory / "rsolve.inp", "lukmt"))
    interface = read_interface(directory / "fort.10", directory / "fort.21")
    for key in ["target_components", "boundary_channels", "symmetry_charge", "geometry_mass"]:
        if header[key] != interface[key]:
            raise ValueError(f"K/interface {key} differs")
    assert len(energy) == contract["expected_energy_points"]
    assert header["open_channels"] == contract["expected_open_channels"]
    assert interface["target_components"] == contract["expected_target_components"][model]
    assert interface["boundary_channels"] == contract["expected_boundary_channels"][model]
    assert interface["boundary_poles"] == contract["expected_boundary_poles"]
    assert interface["multipole_orders"] == contract["expected_multipole_orders"]
    assert all(interface["multipole_nonzero_entries"])
    assert interface["symmetry_charge"][3] == 0
    np.testing.assert_allclose(
        interface["radius_bohr"], contract["expected_radius_bohr"], atol=1e-12, rtol=0
    )
    requested = contract["requested_energy_start_ev"] + contract[
        "requested_energy_step_ev"
    ] * np.arange(len(energy))
    np.testing.assert_allclose(
        energy / contract["native_rydberg_per_requested_ev"],
        requested,
        atol=contract["grid_requested_ev_absolute_budget"],
        rtol=0,
    )
    thresholds = np.asarray(interface["channel_thresholds_rydberg"])
    target_e = np.asarray(interface["target_energies_hartree"])
    pointers = np.asarray(interface["channel_labels"])[:, 1] - 1
    np.testing.assert_allclose(
        thresholds,
        2 * (target_e[pointers] - target_e[0]),
        atol=contract["metadata_absolute_budget"],
        rtol=0,
    )
    open_indices = np.flatnonzero(thresholds < energy.min())
    np.testing.assert_array_equal(open_indices, np.arange(matrices.shape[1]))
    assert np.count_nonzero(thresholds < energy.max()) == len(open_indices)
    phases = np.loadtxt(directory / unit(directory / "eigenp.inp", "luphso"), ndmin=2)
    assert phases.shape == (len(energy), 2) and np.isfinite(phases).all()
    np.testing.assert_allclose(
        phases[:, 0], requested, atol=contract["grid_requested_ev_absolute_budget"], rtol=0
    )
    metrics = {
        key: 0.0
        for key in [
            "unitarity_error",
            "reciprocity_error",
            "spectral_cayley_error",
            "determinant_phase_error_rad_mod_pi",
            "native_phase_error_rad_mod_pi",
            "cayley_condition_number",
        ]
    }
    phase_values = []
    for k, (_, native_phase) in zip(matrices, phases, strict=True):
        s = scattering_matrix(k)
        values, vectors = np.linalg.eigh(k)
        phase = np.sum(np.arctan(values))
        alternate = (vectors * np.exp(2j * np.arctan(values))) @ vectors.T
        measured = {
            "unitarity_error": np.max(abs(s.conj().T @ s - np.eye(len(s)))),
            "reciprocity_error": np.max(abs(s - s.T)),
            "spectral_cayley_error": np.max(abs(s - alternate)),
            "determinant_phase_error_rad_mod_pi": abs(
                np.angle(np.linalg.det(s) * np.exp(-2j * phase))
            )
            / 2,
            "native_phase_error_rad_mod_pi": abs(np.angle(np.exp(2j * (phase - native_phase)))) / 2,
            "cayley_condition_number": np.linalg.cond(np.eye(len(k)) - 1j * k),
        }
        for key, value in measured.items():
            metrics[key] = max(metrics[key], float(value))
        phase_values.append(float(phase))
    passed = metrics["native_phase_error_rad_mod_pi"] <= contract[
        "native_phase_reproduction_budget_rad_mod_pi"
    ] and all(
        value <= contract["matrix_algebra_absolute_budget"]
        for key, value in metrics.items()
        if key not in ["native_phase_error_rad_mod_pi", "cayley_condition_number"]
    )
    return {
        "header": header,
        "interface": interface,
        "maximum_metrics": metrics,
        "requested_energies_ev": requested.tolist(),
        "reconstructed_eigenphase_sums_rad": phase_values,
        "open_channel_labels": [interface["channel_labels"][i] for i in open_indices],
        "first_closed_threshold_ev": float(
            thresholds[len(open_indices) :].min() / contract["native_rydberg_per_requested_ev"]
        ),
        "saved_matrix_reproduction_verdict": bool(passed),
    }


def audit(directory: Path, contract: dict) -> dict:
    start = time.monotonic()
    manifest = json.loads((directory / "input-manifest.json").read_text())
    for name, digest in manifest.items():
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Input hash mismatch: {name}")
    cases = []
    for anchor in contract["anchors"]:
        for model in contract["models"]:
            for sector in contract["sectors"]:
                name = f"{anchor}/{model}/{sector}"
                row = {"case": name, "saved_matrix_reproduction_verdict": None}
                try:
                    case_path = directory / name
                    required = ["fort.10", "fort.21", "rsolve.inp", "eigenp.inp"]
                    required += [
                        unit(case_path / "rsolve.inp", "lukmt"),
                        unit(case_path / "eigenp.inp", "luphso"),
                    ]
                    if any(f"{name}/{file}" not in manifest for file in required):
                        raise ValueError("Case inputs are not completely checksum-bound")
                    row.update(inspect_case(directory / name, model, contract))
                except (ValueError, AssertionError, IndexError, OSError) as error:
                    row["failure"] = f"{type(error).__name__}: {error}"
                cases.append(row)
    return {
        "contract": contract,
        "input_sha256": manifest,
        "cases": cases,
        "saved_matrix_reproduction_verdict": len(cases) == contract["expected_cases"]
        and all(c["saved_matrix_reproduction_verdict"] is True for c in cases),
        "reader_wall_seconds": time.monotonic() - start,
        "independent_outer_propagation_verdict": None,
        "complex_pole_verdict": None,
        "native_solves": 0,
        "full_sweep_released": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.directory, json.loads(args.contract.read_text()))
    with args.output.open("x") as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if report["saved_matrix_reproduction_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
