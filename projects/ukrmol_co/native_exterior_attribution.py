"""Fixed native-subrange and independent-step attribution; original pilot retained."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.linalg import LinAlgWarning

from projects.ukrmol_co.exterior_numerov import numerov_k
from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.native_exterior_study import cayley, matrix_error
from projects.ukrmol_co.saved_matrices import array, read_interface, read_k, records


def subrange_deck(original: bytes, subranges: int) -> bytes:
    """Keep the literal baseline, or insert one declared BPROP line only."""
    if subranges not in (0, 16, 32):
        raise ValueError("Only the frozen baseline/16/32 subrange controls are allowed")
    if re.search(rb"\bnrange\s*=", original, flags=re.I):
        raise ValueError("Original deck already specifies a subrange count")
    lines = original.splitlines(keepends=True)
    indices = [i for i, line in enumerate(lines) if line.lstrip().lower().startswith(b"&bprop")]
    if len(indices) != 1:
        raise ValueError("Exactly one BPROP namelist is required")
    if subranges:
        lines.insert(indices[0] + 1, f"  nrange = {subranges},\n".encode())
    return b"".join(lines)


def differences(left, right):
    """Maximum entrywise K/S and determinant-phase differences on matched grids."""
    lk, rk = np.asarray(left), np.asarray(right)
    ls, rs = (
        np.asarray([cayley(matrix) for matrix in k]) if k.ndim == 3 else cayley(k) for k in [lk, rk]
    )
    phase = abs(np.angle(np.linalg.det(ls) / np.linalg.det(rs))) / 2
    return {
        "k": matrix_error(lk, rk),
        "s": matrix_error(ls, rs),
        "phase_rad_mod_pi": float(np.max(phase)),
    }


def algebra(k):
    s = cayley(k)
    return {
        "unitarity": matrix_error(s.conj().T @ s, np.eye(len(k))),
        "reciprocity": matrix_error(s, s.T),
    }


def attribution_checks(native, independent, prior, contract):
    """Evaluate the fixed control matrix without changing comparison budgets."""
    index = contract["energy_index_zero_based"]
    budget = contract["matrix_and_phase_absolute_budget"]
    latest = independent[-1]
    step = differences(np.array(latest["k"]), np.array(prior["k"]))
    endpoint = differences(np.array(independent[0]["k"]), np.array(latest["k"]))
    native_change = differences(native[1][index], native[2][index])
    match = differences(native[2][index], np.array(latest["k"]))
    baseline_match = differences(native[0][index], np.array(latest["k"]))
    checks = {
        "native_subrange_refinement": all(v <= budget for v in native_change.values()),
        "independent_step_refinement": step["s"] <= budget and step["phase_rad_mod_pi"] <= budget,
        "independent_endpoint_refinement": endpoint["s"] <= budget
        and endpoint["phase_rad_mod_pi"] <= budget,
        "independent_spatial_trend": step["s"]
        <= max(
            contract["prior_step_s_change"] / contract["independent_spatial_improvement_factor"],
            contract["roundoff_floor"],
        ),
        "native_k_reproduction": match["k"] <= budget,
        "native_s_reproduction": match["s"] <= budget,
        "native_phase_reproduction": match["phase_rad_mod_pi"] <= budget,
    }
    invariant_values = [algebra(np.array(row["k"])) for row in independent]
    invariant_values.append(algebra(native[2][index]))
    checks["unitarity"] = max(v["unitarity"] for v in invariant_values) <= budget
    checks["reciprocity"] = max(v["reciprocity"] for v in invariant_values) <= budget
    reduction = 1 - match["s"] / baseline_match["s"] if baseline_match["s"] else None
    # The correction and residual are full matrices, not only their norms.
    old_s, new_s = cayley(native[0][index]), cayley(native[2][index])
    residual = new_s - cayley(np.array(latest["k"]))
    correction = new_s - old_s
    return checks, {
        "native_final_change": native_change,
        "independent_step_change": step,
        "independent_endpoint_change": endpoint,
        "refined_native_vs_independent": match,
        "baseline_native_vs_independent": baseline_match,
        "s_mismatch_fractional_reduction": reduction,
        "native_correction_s_real": correction.real.tolist(),
        "native_correction_s_imag": correction.imag.tolist(),
        "refined_residual_s_real": residual.real.tolist(),
        "refined_residual_s_imag": residual.imag.tolist(),
        "native_subrange_attribution": bool(
            all(checks.values())
            and reduction is not None
            and reduction >= contract["attribution_minimum_fractional_s_mismatch_reduction"]
        ),
    }


def run(directory, prior_path, replay_path, engine, output, contract, report):
    raw = prior_path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != contract["prior_pilot_sha256"]:
        raise ValueError("Prior rejected pilot digest differs")
    prior = json.loads(raw)
    replay = json.loads(replay_path.read_text())
    if any(replay.get(k) != v for k, v in contract["public_replay_required"].items()):
        raise ValueError("Required public reconstruction of the original rejection is missing")
    if prior["pilot_verdict"] is not False:
        raise ValueError("Original pilot must retain its false verdict")
    if hashlib.sha256(engine.read_bytes()).hexdigest() != contract["engine_rsolve_sha256"]:
        raise ValueError("Pinned rsolve executable differs")
    for name, want in contract["input_sha256"].items():
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != want:
            raise ValueError(f"Native input differs: {name}")
    source = directory / contract["case"]
    energies, saved, _ = read_k(source / "fort.921")
    if len(energies) != contract["native_energy_points"]:
        raise ValueError("Saved energy sequence is incomplete")
    native = []
    for subranges, printed in zip(
        contract["native_subranges"], contract["expected_printed_subranges"], strict=True
    ):
        case = output / f"native-{printed}"
        case.mkdir()
        for name in ["fort.10", "fort.21"]:
            shutil.copy2(source / name, case / name)
        (case / "rsolve.inp").write_bytes(
            subrange_deck((source / "rsolve.inp").read_bytes(), subranges)
        )
        row = {"subranges": printed, "completed": False}
        report["native_cases"].append(row)
        (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        begin = time.monotonic()
        with (case / "rsolve.inp").open("rb") as inp, (case / "rsolve.out").open("xb") as log:
            try:
                result = subprocess.run(
                    [str(engine)],
                    cwd=case,
                    stdin=inp,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=contract["resources"]["native_stage_timeout_seconds"],
                )
                row["original_exit"] = result.returncode
            except subprocess.TimeoutExpired:
                row.update(original_exit=124, timed_out=True)
        row["wall_seconds"] = time.monotonic() - begin
        row["input_sha256"] = {
            name: hashlib.sha256((case / name).read_bytes()).hexdigest()
            for name in ["fort.10", "fort.21", "rsolve.inp"]
        }
        (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        for name in ["fort.10", "fort.21"]:
            if row["input_sha256"][name] != contract["input_sha256"][contract["case"] + "/" + name]:
                raise ValueError(f"Native execution changed {name}")
        expected_deck = subrange_deck((source / "rsolve.inp").read_bytes(), subranges)
        if (case / "rsolve.inp").read_bytes() != expected_deck:
            raise ValueError("Native execution changed the declared deck")
        if row["original_exit"] != 0:
            raise ValueError(f"Native execution failed for {printed} subranges")
        text = (case / "rsolve.out").read_text()
        if not re.search(rf"PROPAGATED ACROSS\s+{printed} SUBRANGES", text):
            raise ValueError("Native printed subrange count differs")
        current_energies, k, _ = read_k(case / "fort.921")
        np.testing.assert_allclose(current_energies, energies, atol=1e-14, rtol=0)
        row.update(completed=True, saved_native_differences=differences(k, saved))
        (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        native.append(k)
        if subranges == 0:
            report["checks"]["baseline_reproduction"] = all(
                value <= contract["baseline_absolute_budget"]
                for value in row["saved_native_differences"].values()
            )
            if not report["checks"]["baseline_reproduction"]:
                raise ValueError("Literal same-engine baseline reproduction failed")
    interface = read_interface(source / "fort.10", source / "fort.21")
    count = interface["boundary_channels"]
    raw = records(source / "fort.21")
    cf = np.zeros((2, count, count))
    packed = array(raw[5], "<f8", count * (count + 1)).reshape(2, -1)
    indices = np.tril_indices(count)
    for matrix, values in zip(cf, packed, strict=True):
        matrix[indices] = values
        matrix.T[indices] = values
    labels = np.asarray(interface["channel_labels"])
    angular = labels[:, 2]
    thresholds = np.asarray(interface["channel_thresholds_rydberg"])
    index = contract["energy_index_zero_based"]
    energy = energies[index]
    poles = array(raw[6], "<f8", 2048)
    amplitudes = array(raw[7], "<f8", count * 2048).reshape(2048, count)
    total = interface["target_energies_hartree"][0] + energy / 2
    q = surface_matrix(poles, amplitudes, total)
    with (output / "independent.jsonl").open("x") as journal:
        for endpoint in contract["independent_endpoints_bohr"]:
            begin = time.monotonic()
            with warnings.catch_warnings():
                warnings.simplefilter("error", LinAlgWarning)
                k, diagnostics = numerov_k(
                    q,
                    thresholds,
                    angular,
                    cf[0] + np.diag(angular * (angular + 1)),
                    cf[1],
                    energy,
                    18,
                    endpoint,
                    contract["independent_step_bohr"],
                    contract["independent_tail_order"],
                )
            row = {
                "endpoint_bohr": endpoint,
                "step_bohr": contract["independent_step_bohr"],
                "tail_order": contract["independent_tail_order"],
                "k": k.tolist(),
                "diagnostics": diagnostics,
                "wall_seconds": time.monotonic() - begin,
            }
            report["independent_cases"].append(row)
            (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
            journal.write(json.dumps(row, allow_nan=False) + "\n")
            journal.flush()
    previous = next(
        row
        for row in prior["cases"]
        if (row["energy_index"], row["step_bohr"], row["endpoint_bohr"], row["tail_order"])
        == (index, 0.025, 400.0, 8)
    )
    checks, metrics = attribution_checks(native, report["independent_cases"], previous, contract)
    report["checks"].update(checks)
    report.update(metrics=metrics, attribution_control_verdict=bool(all(report["checks"].values())))
    for name, want in contract["input_sha256"].items():
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != want:
            report["attribution_control_verdict"] = None
            raise ValueError(f"Original input changed during attribution: {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["directory", "prior", "replay", "engine", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text())
    args.output.mkdir()
    report = {
        "contract": contract,
        "checks": {},
        "native_cases": [],
        "independent_cases": [],
        "attribution_control_verdict": None,
        "original_pilot_verdict": False,
        "native_full_k_reconstruction_verdict": None,
        "co_complex_pole_verdict": None,
        "new_inner_solves": 0,
        "full_sweep_released": False,
    }
    begin = time.monotonic()
    try:
        run(args.directory, args.prior, args.replay, args.engine, args.output, contract, report)
    except (
        ValueError,
        RuntimeError,
        OSError,
        AssertionError,
        np.linalg.LinAlgError,
        LinAlgWarning,
    ) as error:
        report["failure"] = f"{type(error).__name__}: {error}"
    report["wall_seconds"] = time.monotonic() - begin
    (args.output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if report["attribution_control_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
