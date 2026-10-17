"""Bounded synthetic timing and prospective saved-native exterior pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.linalg import LinAlgWarning

from projects.ukrmol_co.exterior_numerov import numerov_k
from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.saved_matrices import array, read_interface, read_k, records, unit


def cayley(k):
    identity = np.eye(len(k))
    return np.linalg.solve(identity - 1j * k, identity + 1j * k)


def matrix_error(left, right):
    return float(np.max(abs(left - right)))


def synthetic(contract):
    """Time dense 298-channel linear algebra without reading native data."""
    count = contract["channels"]
    rng = np.random.default_rng(contract["seed"])
    frame = rng.normal(size=(count, 8)) / np.sqrt(count)
    q = np.diag(np.linspace(0.1, 0.4, count)) + frame @ frame.T
    thresholds = np.r_[np.zeros(6), np.linspace(0.4, 4.0, count - 6)]
    angular = np.arange(count) % 6
    c2 = np.diag(angular * (angular + 1)) + 0.2 * (frame @ frame.T)
    c3 = 0.1 * (frame @ frame.T)
    start = time.monotonic()
    with warnings.catch_warnings():
        warnings.simplefilter("error", LinAlgWarning)
        k, diagnostics = numerov_k(
            q,
            thresholds,
            angular,
            c2,
            c3,
            contract["energy_rydberg"],
            contract["radius_bohr"],
            contract["endpoint_bohr"],
            contract["step_bohr"],
            contract["tail_order"],
        )
    return {
        "synthetic_timing_verdict": bool(np.isfinite(k).all()),
        "kernel_wall_seconds": time.monotonic() - start,
        "diagnostics": diagnostics,
        "native_data_read": False,
        "native_solves": 0,
        "native_full_k_reconstruction_verdict": None,
        "co_complex_pole_verdict": None,
        "full_sweep_released": False,
    }


def pilot(directory, contract, journal):
    manifest = json.loads((directory / "input-manifest.json").read_text())
    for name, want in contract["input_sha256"].items():
        if (
            manifest.get(name) != want
            or hashlib.sha256((directory / name).read_bytes()).hexdigest() != want
        ):
            raise ValueError(f"Input hash mismatch: {name}")
    case = directory / contract["case"]
    interface = read_interface(case / "fort.10", case / "fort.21")
    energies, native_k, header = read_k(case / unit(case / "rsolve.inp", "lukmt"))
    count = interface["boundary_channels"]
    if count != contract["expected_channels"] or interface["boundary_poles"] != 2048:
        raise ValueError("Unexpected native boundary dimensions")
    if interface["symmetry_charge"][3] != 0 or interface["radius_bohr"] != 18:
        raise ValueError("Only the saved neutral 18-bohr boundary is supported")
    if len(energies) != 99 or header["open_channels"] != 6:
        raise ValueError("Incomplete native energy/open-channel sequence")
    labels = np.asarray(interface["channel_labels"])
    angular, thresholds = labels[:, 2], np.asarray(interface["channel_thresholds_rydberg"])
    targets = np.asarray(interface["target_energies_hartree"])
    np.testing.assert_allclose(
        thresholds, 2 * (targets[labels[:, 1] - 1] - targets[0]), atol=1e-10, rtol=0
    )
    if np.any(thresholds[:6] != 0) or np.any(thresholds[6:] <= energies.max()):
        raise ValueError("Open-channel ordering/thresholds differ from the declared pilot")
    raw = records(case / "fort.21")
    packed = array(raw[5], "<f8", count * (count + 1)).reshape(2, -1)
    cf = np.zeros((2, count, count))
    indices = np.tril_indices(count)
    for coefficient, values in zip(cf, packed, strict=True):
        coefficient[indices] = values
        coefficient.T[indices] = values
    poles = array(raw[6], "<f8", 2048)
    amplitudes = array(raw[7], "<f8", count * 2048).reshape(2048, count)
    c2, c3 = cf[0] + np.diag(angular * (angular + 1)), cf[1]
    steps, endpoints, orders = (
        contract[key] for key in ["steps_bohr", "endpoints_bohr", "tail_orders"]
    )
    configurations = list(
        dict.fromkeys(
            [(h, endpoints[-1], orders[-1]) for h in steps]
            + [(steps[-1], r, orders[-1]) for r in endpoints]
            + [(steps[-1], endpoints[-1], p) for p in orders]
        )
    )
    rows = []
    start = time.monotonic()
    with journal.open("x") as stream:
        for index in contract["native_energy_indices_zero_based"]:
            energy, reference = energies[index], native_k[index]
            reference_s = cayley(reference)
            q = surface_matrix(poles, amplitudes, targets[0] + energy / 2)
            for h, endpoint, order in configurations:
                row = {
                    "energy_index": index,
                    "requested_energy_ev": float(energy / 0.0735),
                    "step_bohr": h,
                    "endpoint_bohr": endpoint,
                    "tail_order": order,
                    "completed": False,
                }
                begin = time.monotonic()
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("error", LinAlgWarning)
                        k, diagnostics = numerov_k(
                            q,
                            thresholds,
                            angular,
                            c2,
                            c3,
                            energy,
                            18,
                            endpoint,
                            h,
                            order,
                            contract["threshold_degeneracy_rydberg"],
                        )
                    s = cayley(k)
                    row.update(
                        {
                            "completed": True,
                            "k": k.tolist(),
                            "s_real": s.real.tolist(),
                            "s_imag": s.imag.tolist(),
                            "diagnostics": diagnostics,
                            "native_k_error": matrix_error(k, reference),
                            "native_s_error": matrix_error(s, reference_s),
                            "native_phase_error_rad_mod_pi": float(
                                abs(np.angle(np.linalg.det(s) / np.linalg.det(reference_s))) / 2
                            ),
                            "unitarity_error": matrix_error(s.conj().T @ s, np.eye(len(s))),
                            "reciprocity_error": matrix_error(s, s.T),
                        }
                    )
                except (ValueError, RuntimeError, np.linalg.LinAlgError, LinAlgWarning) as error:
                    row["failure"] = f"{type(error).__name__}: {error}"
                row["wall_seconds"] = time.monotonic() - begin
                rows.append(row)
                stream.write(json.dumps(row, allow_nan=False) + "\n")
                stream.flush()
    checks = {
        "complete_matrix": len(rows) == contract["expected_cases"]
        and all(r["completed"] for r in rows)
    }
    refinement = {}
    if checks["complete_matrix"]:
        finest = [
            r
            for r in rows
            if (r["step_bohr"], r["endpoint_bohr"], r["tail_order"])
            == (steps[-1], endpoints[-1], orders[-1])
        ]
        for metric, budget in contract["finest_budgets"].items():
            checks[metric] = max(r[metric] for r in finest) <= budget
        for field, family in [
            ("step_bohr", steps),
            ("endpoint_bohr", endpoints),
            ("tail_order", orders),
        ]:
            differences = [0.0, 0.0]
            for index in contract["native_energy_indices_zero_based"]:
                group = [
                    r
                    for r in rows
                    if r["energy_index"] == index
                    and all(
                        r[k] == v
                        for k, v in [
                            ("step_bohr", steps[-1]),
                            ("endpoint_bohr", endpoints[-1]),
                            ("tail_order", orders[-1]),
                        ]
                        if k != field
                    )
                ]
                group.sort(key=lambda r: family.index(r[field]))
                for i in range(2):
                    a = np.array(group[i]["s_real"]) + 1j * np.array(group[i]["s_imag"])
                    b = np.array(group[i + 1]["s_real"]) + 1j * np.array(group[i + 1]["s_imag"])
                    differences[i] = max(differences[i], matrix_error(a, b))
            refinement[field] = differences
            checks[field + "_refinement"] = (
                differences[-1] <= contract["refinement_s_absolute_budget"]
            )
        checks["spatial_fourth_order_trend"] = refinement["step_bohr"][1] <= max(
            refinement["step_bohr"][0] / 8, contract["roundoff_floor"]
        )
    return {
        "cases": rows,
        "checks": checks,
        "refinement_s": refinement,
        "pilot_verdict": bool(all(checks.values())),
        "wall_seconds": time.monotonic() - start,
        "native_data_read": True,
        "native_solves": 0,
        "native_full_k_reconstruction_verdict": None,
        "co_complex_pole_verdict": None,
        "full_sweep_released": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--directory", type=Path)
    args = parser.parse_args()
    raw = args.contract.read_bytes()
    contract = json.loads(raw)
    if contract["mode"] == "synthetic":
        report = synthetic(contract)
        verdict = report["synthetic_timing_verdict"]
    else:
        report = pilot(args.directory, contract, args.output.with_suffix(".jsonl"))
        verdict = report["pilot_verdict"]
    report.update({"contract": contract, "contract_sha256": hashlib.sha256(raw).hexdigest()})
    with args.output.open("x") as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if verdict is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
