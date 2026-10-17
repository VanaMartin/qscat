"""Read-only replay of saved six-case decisions and Cayley-coordinate diagnostics."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np

from projects.ukrmol_co.saved_matrices import read_k, unit


def scattering(k):
    """Right-solve the Cayley map, independently of the producer's left solve."""
    identity = np.eye(k.shape[-1])
    return np.linalg.solve((identity - 1j * k).T, (identity + 1j * k).T).T


def compare(left, right):
    ls, rs = scattering(left), scattering(right)
    # Eigenvalue arctangents provide a separate route to determinant phase.
    lp = np.arctan(np.linalg.eigvals(left)).sum().real
    rp = np.arctan(np.linalg.eigvals(right)).sum().real
    phase = abs((lp - rp + np.pi / 2) % np.pi - np.pi / 2)
    return {
        "k": float(np.max(abs(left - right))),
        "s": float(np.max(abs(ls - rs))),
        "phase_rad_mod_pi": float(phase),
    }


def coordinate_diagnostic(left, right):
    """Explain saved differences through an exact identity, without changing gates.

    S(A)-S(B) = 2i (I-iA)^-1 (A-B) (I-iB)^-1. Native B is
    symmetric; its orthonormal eigenframe locates the large-tangent channel.
    This does not locate a complex pole or infer a convergence bound.
    """
    identity = np.eye(len(left))
    delta = left - right
    predicted = 2j * np.linalg.solve(identity - 1j * left, delta)
    predicted = np.linalg.solve((identity - 1j * right).T, predicted.T).T
    actual = scattering(left) - scattering(right)
    values, vectors = np.linalg.eigh(right)
    index = int(np.argmax(abs(values)))
    vector = vectors[:, index]
    coefficient = float(vector @ delta @ vector)
    norm = float(np.linalg.norm(delta))
    return {
        "max_abs_native_k": float(np.max(abs(right))),
        "max_abs_independent_k": float(np.max(abs(left))),
        "largest_tangent_eigenvalue": float(values[index]),
        "distance_to_tangent_singularity_rad": float(np.pi / 2 - abs(np.arctan(values[index]))),
        "large_channel_local_cayley_derivative": float(2 / (1 + values[index] ** 2)),
        "relative_frobenius_k_change": norm / float(np.linalg.norm(right)),
        "large_channel_delta_frobenius_fraction": abs(coefficient) / norm if norm else 0.0,
        "cayley_difference_identity_max_residual": float(np.max(abs(predicted - actual))),
        "differences": compare(left, right),
    }


def expected_configurations(contract, probe):
    steps, endpoints, orders = (
        contract[k] for k in ["steps_bohr", "endpoints_bohr", "tail_orders"]
    )
    common = contract["comparison_endpoint_bohr"]
    result = [(h, common, orders[-1]) for h in steps[-(3 if probe else 2) :]]
    if probe:
        result += [(steps[-1], end, orders[-1]) for end in endpoints if end != common]
        result += [(steps[-1], common, order) for order in orders[:-1]]
    return result


def replay_energy(bundle, reference, contract):
    rows = bundle["rows"]
    keys = [(row["step_bohr"], row["endpoint_bohr"], row["tail_order"]) for row in rows]
    expected = expected_configurations(contract, bundle["probe"])
    if keys != expected or not all(row["completed"] for row in rows):
        raise ValueError("Incomplete or duplicate independent configurations")
    matrices = {key: np.asarray(row["k"]) for key, row in zip(keys, rows, strict=True)}
    steps, orders = contract["steps_bohr"], contract["tail_orders"]
    common = contract["comparison_endpoint_bohr"]
    families = {"step": [(h, common, orders[-1]) for h in steps[-(3 if bundle["probe"] else 2) :]]}
    if bundle["probe"]:
        families["endpoint"] = [(steps[-1], end, orders[-1]) for end in contract["endpoints_bohr"]]
        families["order"] = [(steps[-1], common, order) for order in orders]
    finest = matrices[steps[-1], common, orders[-1]]
    metrics = {"native_comparison": compare(finest, reference)}
    checks = {}
    budget = contract["matrix_and_phase_absolute_budget"]
    for name, sequence in families.items():
        changes = [
            compare(matrices[sequence[i]], matrices[sequence[i + 1]])
            for i in range(len(sequence) - 1)
        ]
        metrics[name + "_changes"] = changes
        checks.update({name + "_" + k: v <= budget for k, v in changes[-1].items()})
    checks.update({"native_" + k: v <= budget for k, v in metrics["native_comparison"].items()})
    if bundle["probe"]:
        changes = metrics["step_changes"]
        checks["spatial_trend"] = changes[-1]["s"] <= max(
            changes[0]["s"] / contract["spatial_improvement_factor"], contract["roundoff_floor"]
        )
    invariant_rows = []
    for row in rows:
        k = np.asarray(row["k"])
        if k.shape != (6, 6) or not np.isfinite(k).all():
            raise ValueError("Nonfinite or incompatible K matrix")
        s = scattering(k)
        np.testing.assert_allclose(s.real, row["s_real"], atol=1e-13, rtol=0)
        np.testing.assert_allclose(s.imag, row["s_imag"], atol=1e-13, rtol=0)
        invariant_rows.append(
            {
                "unitarity": float(np.max(abs(s @ s.conj().T - np.eye(6)))),
                "reciprocity": float(np.max(abs(s - s.T))),
            }
        )
    metrics["invariant_maxima"] = {
        k: max(row[k] for row in invariant_rows) for k in invariant_rows[0]
    }
    checks.update(
        {
            k: v <= contract["invariant_absolute_budget"]
            for k, v in metrics["invariant_maxima"].items()
        }
    )
    if checks != bundle["checks"] or all(checks.values()) is not bundle["verdict"]:
        raise ValueError("Independent decision replay differs")
    for name, values in metrics.items():
        old = bundle["metrics"][name]
        for a, b in zip(
            values if isinstance(values, list) else [values],
            old if isinstance(old, list) else [old],
            strict=True,
        ):
            for key in a:
                np.testing.assert_allclose(a[key], b[key], atol=1e-13, rtol=0)
    return metrics


def replay(root, raw):
    """Verify every saved native record, journal and finite-grid decision."""
    report = json.loads((root / "result/report.json").read_text())
    execution = json.loads((root / "execution.json").read_text())
    contract = report["contract"]
    for name, want in json.loads((root / "result-collection-manifest.json").read_text()).items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != want:
            raise ValueError(f"Collected output digest differs: {name}")
    for name, want in contract["input_sha256"].items():
        if hashlib.sha256((raw / name).read_bytes()).hexdigest() != want:
            raise ValueError(f"Native input digest differs: {name}")
    native, native_deltas = {}, []
    for case in contract["cases"]:
        output_unit = unit(raw / case / "rsolve.inp", "lukmt")
        grid, saved, header = read_k(raw / case / output_unit)
        matrices = []
        producer = next(row for row in report["native"] if row["case"] == case)
        for index, count in enumerate([8, 16, 32]):
            directory = root / "result/native" / case / str(count)
            deck = (directory / "rsolve.inp").read_bytes()
            if count != 8:
                inserted = f"  nrange = {count},\n".encode()
                if deck.count(inserted) != 1:
                    raise ValueError("Missing or duplicated one-line native mutation")
                deck = deck.replace(inserted, b"")
            if deck != (raw / case / "rsolve.inp").read_bytes():
                raise ValueError("Literal baseline deck differs")
            if not re.search(
                rf"PROPAGATED ACROSS\s+{count} SUBRANGES", (directory / "rsolve.out").read_text()
            ):
                raise ValueError("Native printed subrange count differs")
            energies, k, current_header = read_k(directory / output_unit)
            np.testing.assert_allclose(energies, grid, atol=1e-14, rtol=0)
            if current_header != header or k.shape != (99, 6, 6):
                raise ValueError("Native dimensions or complete grid differ")
            if (
                producer["stages"][index]["original_exit"] != 0
                or not producer["stages"][index]["completed"]
            ):
                raise ValueError("Incomplete native stage")
            matrices.append(k)
        if not np.array_equal(matrices[0], saved):
            raise ValueError("Native baseline is not exactly reproduced")
        deltas = [compare(a, b) for a, b in zip(matrices[1], matrices[2], strict=True)]
        maxima = {key: max(row[key] for row in deltas) for key in deltas[0]}
        checks = {
            "baseline": True,
            **{k: v <= contract["matrix_and_phase_absolute_budget"] for k, v in maxima.items()},
        }
        if checks != producer["checks"] or not producer["verdict"]:
            raise ValueError("Native gate replay differs")
        native[case] = matrices[-1]
        native_deltas.append(maxima)
    identities, count, diagnoses = set(), 0, []
    for bundle in report["probes"] + report["remainder"]:
        case, index = bundle["case"], bundle["energy_index"]
        identity = (case, index)
        if identity in identities or case not in contract["cases"] or not 0 <= index < 99:
            raise ValueError("Duplicate or unexpected energy identity")
        if bundle["probe"] != (index in contract["probe_energy_indices"]):
            raise ValueError("Probe identity differs")
        identities.add(identity)
        path = root / "result/independent" / case / f"{index:03d}"
        if json.loads(path.with_suffix(".json").read_text()) != bundle:
            raise ValueError("Per-energy report differs from aggregate")
        journal = [json.loads(line) for line in path.with_suffix(".jsonl").read_text().splitlines()]
        if journal != bundle["rows"]:
            raise ValueError("Persistent journal differs from report")
        for row in journal:
            np.testing.assert_allclose(
                row["requested_energy_ev"], 0.1 + 0.05 * index, atol=1e-11, rtol=0
            )
            np.testing.assert_allclose(
                row["energy_rydberg"], (0.1 + 0.05 * index) * 0.0735, atol=1e-14, rtol=0
            )
        replay_energy(bundle, native[case][index], contract)
        count += len(journal)
        if bundle["verdict"] is False:
            diagnoses.append(
                {
                    "case": case,
                    "energy_index": index,
                    "requested_energy_ev": journal[-1]["requested_energy_ev"],
                    "failed_checks": [k for k, v in bundle["checks"].items() if not v],
                    "metrics": bundle["metrics"],
                    "coordinate_diagnostic": coordinate_diagnostic(
                        np.asarray(journal[-1]["k"]), native[case][index]
                    ),
                    "adjacent_native_max_abs_k": {
                        str(j): float(np.max(abs(native[case][j])))
                        for j in [index - 1, index, index + 1]
                    },
                }
            )
    expected_probes = {
        (case, index) for case in contract["cases"] for index in contract["probe_energy_indices"]
    }
    if {(x["case"], x["energy_index"]) for x in report["probes"]} != expected_probes:
        raise ValueError("Incomplete probe coverage")
    if len(list((root / "result/independent").rglob("*.json"))) != len(identities):
        raise ValueError("Unreported independent output")
    if not all(x["verdict"] is True for x in report["probes"]) or not diagnoses:
        raise ValueError("Unexpected stopped-run verdict")
    if (
        report["stage_verdicts"] != {"native": True, "probes": True, "remainder": False}
        or report["reconstruction_verdict"] is not False
    ):
        raise ValueError("Stage verdicts differ")
    if execution["original_command_exit"] != 1 or execution["original_container_exit"] != 1:
        raise ValueError("Original rejected exits differ")
    return {
        "saved_decision_replay_passed": True,
        "reconstruction_verdict": False,
        "original_command_container_exits": [1, 1],
        "native_cases_verified": 6,
        "native_energy_records_verified": 1782,
        "native_baseline_matrices_exact": 594,
        "native16_to32_maxima": {
            key: max(row[key] for row in native_deltas) for key in native_deltas[0]
        },
        "independent_energy_identities_verified": len(identities),
        "independent_solves_journal_verified": count,
        "probes_passed": len(report["probes"]),
        "remainder_passed": sum(x["verdict"] is True for x in report["remainder"]),
        "remainder_rejected": len(diagnoses),
        "rejections": diagnoses,
        "original_pilot_verdict": False,
        "original_pilot_exits": [1, 1],
        "native_full_k_reconstruction_verdict": False,
        "co_complex_pole_verdict": None,
        "new_inner_solves": 0,
        "full_sweep_released": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "raw", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = replay(args.root, args.raw)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
