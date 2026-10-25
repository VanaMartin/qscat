"""Gated six-case saved-native outer references and independent real-axis reconstruction."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import time
import warnings
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.linalg import LinAlgWarning

from projects.ukrmol_co.exterior_numerov import numerov_k
from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.native_exterior_attribution import algebra, differences, subrange_deck
from projects.ukrmol_co.native_exterior_study import cayley
from projects.ukrmol_co.native_reconstruction_checks import (
    complete_coverage,
    configurations,
    energy_checks,
)
from projects.ukrmol_co.saved_matrices import array, read_interface, read_k, records, unit


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, data):
    temporary = path.with_suffix(".partial")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def verify_inputs(directory, contract):
    for name, want in contract["input_sha256"].items():
        if digest(directory / name) != want:
            raise ValueError(f"Native input digest differs: {name}")


@lru_cache(maxsize=6)
def load_model(source):
    """Decode unchanged native boundary coefficients; retain every closed channel."""
    interface = read_interface(source / "fort.10", source / "fort.21")
    return _model_from_interface(source, interface, (48, 298))


def _model_from_interface(source, interface, expected):
    """Shared native assembly for the declared48/298 and explicit40/250 views."""
    if expected not in [(48, 298), (40, 250)]:
        raise ValueError("Only the two declared compressed target/channel layouts are supported")
    target_count, channel_count = expected
    output_unit = unit(source / "rsolve.inp", "lukmt")
    energies, saved, header = read_k(source / output_unit)
    for key in ["target_components", "boundary_channels", "symmetry_charge", "geometry_mass"]:
        if header[key] != interface[key]:
            raise ValueError(f"K/interface {key} differs")
    if (
        interface["boundary_channels"] != channel_count
        or interface["boundary_poles"] != 2048
        or interface["target_components"] != target_count
        or interface["multipole_orders"] != 2
        or interface["symmetry_charge"][3] != 0
        or abs(interface["radius_bohr"] - 18) > 1e-12
        or saved.shape != (99, 6, 6)
    ):
        raise ValueError("Only the declared six-open-channel/2048-pole model view is supported")
    np.testing.assert_allclose(energies / 0.0735, 0.1 + 0.05 * np.arange(99), atol=1e-11, rtol=0)
    labels = np.asarray(interface["channel_labels"])
    angular = labels[:, 2]
    thresholds = np.asarray(interface["channel_thresholds_rydberg"])
    targets = np.asarray(interface["target_energies_hartree"])
    target_rows = np.asarray(interface.get("channel_target_rows", labels[:, 1]))
    np.testing.assert_allclose(
        thresholds, 2 * (targets[target_rows - 1] - targets[0]), atol=1e-10, rtol=0
    )
    if np.any(thresholds[:6] != 0) or np.any(thresholds[6:] <= energies.max()):
        raise ValueError("Declared open ordering or closed-channel grid coverage differs")
    raw = records(source / "fort.21")
    count = len(thresholds)
    cf = np.zeros((2, count, count))
    packed = array(raw[5], "<f8", count * (count + 1)).reshape(2, -1)
    indices = np.tril_indices(count)
    for matrix, values in zip(cf, packed, strict=True):
        matrix[indices] = values
        matrix.T[indices] = values
    return {
        "energies": energies,
        "saved": saved,
        "header": header,
        "unit": output_unit,
        "thresholds": thresholds,
        "angular": angular,
        "ground": targets[0],
        "poles": array(raw[6], "<f8", 2048),
        "amplitudes": array(raw[7], "<f8", count * 2048).reshape(2048, count),
        "c2": cf[0] + np.diag(angular * (angular + 1)),
        "c3": cf[1],
    }


def native_case(task):
    directory, output, engine, case, contract = task
    source = directory / case
    destination = output / "native" / case
    destination.mkdir(parents=True)
    result = {"case": case, "verdict": None, "stages": [], "checks": {}}
    try:
        model = load_model(source)
        matrices = []
        for subranges, printed in zip(
            contract["native_subranges"], contract["expected_printed_subranges"], strict=True
        ):
            run = destination / str(printed)
            run.mkdir()
            for name in ["fort.10", "fort.21"]:
                shutil.copy2(source / name, run / name)
            deck = subrange_deck((source / "rsolve.inp").read_bytes(), subranges)
            (run / "rsolve.inp").write_bytes(deck)
            row = {
                "subranges": printed,
                "original_exit": None,
                "timed_out": False,
                "completed": False,
            }
            result["stages"].append(row)
            save(destination / "report.json", result)
            begin = time.monotonic()
            with (run / "rsolve.inp").open("rb") as inp, (run / "rsolve.out").open("xb") as log:
                with subprocess.Popen(
                    [str(engine)], cwd=run, stdin=inp, stdout=log, stderr=subprocess.STDOUT
                ) as process:
                    try:
                        row["original_exit"] = process.wait(
                            timeout=contract["resources"]["native_stage_timeout_seconds"]
                        )
                    except subprocess.TimeoutExpired:
                        process.kill()
                        row.update(original_exit=process.wait(), timed_out=True)
            row["wall_seconds"] = time.monotonic() - begin
            save(destination / "report.json", result)
            for name in ["fort.10", "fort.21"]:
                if digest(run / name) != contract["input_sha256"][case + "/" + name]:
                    raise ValueError("Native execution changed its input")
            if (run / "rsolve.inp").read_bytes() != deck:
                raise ValueError("Native execution changed its deck")
            if row["timed_out"] or row["original_exit"] != 0:
                raise ValueError("Native execution failed")
            if not re.search(
                rf"PROPAGATED ACROSS\s+{printed} SUBRANGES", (run / "rsolve.out").read_text()
            ):
                raise ValueError("Printed native subranges differ")
            energies, k, header = read_k(run / model["unit"])
            np.testing.assert_allclose(energies, model["energies"], atol=1e-14, rtol=0)
            if header != model["header"]:
                raise ValueError("Native output header differs")
            row.update(completed=True, saved_differences=differences(k, model["saved"]))
            row["output_sha256"] = {
                name: digest(run / name) for name in ["rsolve.inp", "rsolve.out", model["unit"]]
            }
            matrices.append(k)
            if not subranges:
                result["checks"]["baseline"] = all(
                    v <= contract["baseline_absolute_budget"]
                    for v in row["saved_differences"].values()
                )
                if not result["checks"]["baseline"]:
                    result["verdict"] = False
                    break
            save(destination / "report.json", result)
        if len(matrices) == 3:
            result["final_change"] = differences(matrices[1], matrices[2])
            result["per_energy_final_changes"] = [
                differences(a, b) for a, b in zip(matrices[1], matrices[2], strict=True)
            ]
            result["checks"].update(
                {
                    name: value <= contract["matrix_and_phase_absolute_budget"]
                    for name, value in result["final_change"].items()
                }
            )
            result["verdict"] = all(result["checks"].values())
    except (ValueError, RuntimeError, OSError, AssertionError, np.linalg.LinAlgError) as error:
        result["failure"] = f"{type(error).__name__}: {error}"
    save(destination / "report.json", result)
    return result


def independent_energy(task):
    directory, output, case, index, contract, probe = task
    destination = output / "independent" / case
    destination.mkdir(parents=True, exist_ok=True)
    result = {"case": case, "energy_index": index, "probe": probe, "verdict": None, "rows": []}
    try:
        model = load_model(directory / case)
        energy = model["energies"][index]
        _, native, _ = read_k(output / "native" / case / "32" / model["unit"])
        q = surface_matrix(model["poles"], model["amplitudes"], model["ground"] + energy / 2)
        with (destination / f"{index:03d}.jsonl").open("x") as stream:
            for step, endpoint, order in configurations(contract, probe):
                begin = time.monotonic()
                with warnings.catch_warnings():
                    warnings.simplefilter("error", LinAlgWarning)
                    k, diagnostics = numerov_k(
                        q,
                        model["thresholds"],
                        model["angular"],
                        model["c2"],
                        model["c3"],
                        energy,
                        18,
                        endpoint,
                        step,
                        order,
                        contract["threshold_degeneracy_rydberg"],
                    )
                s = cayley(k)
                row = {
                    "completed": True,
                    "requested_energy_ev": float(energy / 0.0735),
                    "energy_rydberg": float(energy),
                    "step_bohr": step,
                    "endpoint_bohr": endpoint,
                    "tail_order": order,
                    "k": k.tolist(),
                    "s_real": s.real.tolist(),
                    "s_imag": s.imag.tolist(),
                    "phase_rad_mod_pi": float(np.angle(np.linalg.det(s)) / 2),
                    "diagnostics": diagnostics,
                    "invariants": algebra(k),
                    "wall_seconds": time.monotonic() - begin,
                }
                result["rows"].append(row)
                stream.write(json.dumps(row, allow_nan=False) + "\n")
                stream.flush()
        checks, metrics = energy_checks(result["rows"], native[index], contract, probe)
        result.update(checks=checks, metrics=metrics, verdict=all(checks.values()))
    except (
        ValueError,
        RuntimeError,
        OSError,
        AssertionError,
        np.linalg.LinAlgError,
        LinAlgWarning,
    ) as error:
        result["failure"] = f"{type(error).__name__}: {error}"
    save(destination / f"{index:03d}.json", result)
    return result


def gated_pool(executor, worker, tasks, workers, receive):
    """Bound in-flight bundles; inspect a whole completed batch before replenishing."""
    tasks = iter(tasks)
    pending = set()
    verdict = True
    exhausted = False
    while pending or not exhausted:
        while verdict is True and not exhausted and len(pending) < workers:
            task = next(tasks, None)
            if task is None:
                exhausted = True
            else:
                pending.add(executor.submit(worker, task))
        if not pending:
            break
        done, pending = wait(pending, return_when=FIRST_COMPLETED)
        for future in done:
            row = future.result()
            receive(row)
            if row["verdict"] is not True:
                verdict = False if row["verdict"] is False or verdict is False else None
        if verdict is not True:
            exhausted = True
    return verdict


def run(args, contract, report):
    verify_inputs(args.directory, contract)
    if digest(args.engine) != contract["engine_rsolve_sha256"]:
        raise ValueError("Pinned native executable differs")
    for name, want in contract["prerequisite_sha256"].items():
        path = args.prerequisites / name
        if digest(path) != want:
            raise ValueError("Attribution prerequisite digest differs")
        previous = json.loads(path.read_text())
        if (
            previous["attribution_control_verdict"] is not True
            or previous["original_pilot_verdict"] is not False
        ):
            raise ValueError("Attribution prerequisite verdict differs")
        if "public" in name and previous["public_reconstruction_passed"] is not True:
            raise ValueError("Public attribution replay is missing")

    def receive(stage, row):
        report[stage].append(row)
        save(args.output / "report.json", report)

    with ProcessPoolExecutor(max_workers=contract["resources"]["workers"]) as executor:
        tasks = [
            (args.directory, args.output, args.engine, case, contract) for case in contract["cases"]
        ]
        verdict = gated_pool(executor, native_case, tasks, 6, lambda row: receive("native", row))
        report["stage_verdicts"]["native"] = verdict
        save(args.output / "report.json", report)
        if verdict is not True:
            report["reconstruction_verdict"] = verdict
            return
        for stage, probe in [("probes", True), ("remainder", False)]:
            tasks = [
                (args.directory, args.output, case, index, contract, probe)
                for index in (
                    contract["probe_energy_indices"]
                    if probe
                    else range(contract["expected_energy_points"])
                )
                if probe or index not in contract["probe_energy_indices"]
                for case in contract["cases"]
            ]
            verdict = gated_pool(
                executor,
                independent_energy,
                tasks,
                contract["resources"]["workers"],
                lambda row, stage=stage: receive(stage, row),
            )
            report["stage_verdicts"][stage] = verdict
            save(args.output / "report.json", report)
            if verdict is not True:
                report["reconstruction_verdict"] = verdict
                return
    report["reconstruction_verdict"] = complete_coverage(report, contract)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["directory", "output", "engine", "contract", "prerequisites"]:
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir()
    contract = json.loads(args.contract.read_text())
    report = {
        "contract": contract,
        "contract_sha256": digest(args.contract),
        "native": [],
        "probes": [],
        "remainder": [],
        "stage_verdicts": {},
        "reconstruction_verdict": None,
        "original_pilot_verdict": False,
        "original_pilot_exits": [1, 1],
        "original_saved_matrix_aggregate_verdict": False,
        "five_a2_interface_verdict": None,
        "co_complex_pole_verdict": None,
        "new_inner_solves": 0,
        "full_sweep_released": False,
    }
    begin = time.monotonic()
    save(args.output / "report.json", report)
    try:
        run(args, contract, report)
    except (ValueError, RuntimeError, OSError, AssertionError) as error:
        report.update(reconstruction_verdict=None, failure=f"{type(error).__name__}: {error}")
    finally:
        try:
            verify_inputs(args.directory, contract)
            if digest(args.engine) != contract["engine_rsolve_sha256"]:
                raise ValueError("Native executable changed")
            report["input_and_executable_digests_verified_after"] = True
        except (ValueError, OSError) as error:
            report.update(reconstruction_verdict=None, final_integrity_failure=str(error))
        report["wall_seconds"] = time.monotonic() - begin
        save(args.output / "report.json", report)
    if report["reconstruction_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
