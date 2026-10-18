"""Fixed two-identity native/spacing/endpoint/tail attribution of a preserved K rejection."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.linalg import LinAlgWarning

from projects.ukrmol_co.exterior_numerov import numerov_k
from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.native_exterior_attribution import algebra, differences, subrange_deck
from projects.ukrmol_co.native_exterior_study import cayley
from projects.ukrmol_co.native_full_reconstruction import digest, gated_pool, load_model, save
from projects.ukrmol_co.native_reconstruction_replay import coordinate_diagnostic
from projects.ukrmol_co.saved_matrices import read_k


def single_energy_deck(original, requested_literal, subranges):
    """Only the frozen count/first-energy/subrange mutations are permitted."""
    if subranges not in [32, 64] or requested_literal != "0.99999999999999978":
        raise ValueError("Only the frozen one-energy 32/64 controls are allowed")
    replacements = [
        (rb"\bnescat\s*=\s*99(?=\s*,)", b"nescat = 1"),
        (rb"\beinc\s*=\s*0\.1,0\.05,", b"einc = " + requested_literal.encode() + b",0.05,"),
    ]
    deck = original
    for pattern, replacement in replacements:
        if len(re.findall(pattern, deck, flags=re.I)) != 1:
            raise ValueError("Exactly one frozen energy count and range are required")
        deck = re.sub(pattern, replacement, deck, flags=re.I)
    deck = subrange_deck(deck, 32)
    return deck.replace(b"  nrange = 32,\n", b"  nrange = 64,\n") if subranges == 64 else deck


def native_case(task):
    raw, prior, output, engine, case, contract = task
    destination = output / "native" / case
    destination.mkdir(parents=True)
    result = {"case": case, "verdict": None, "stages": []}
    try:
        model = load_model(raw / case)
        index = contract["energy_index"]
        expected = model["energies"][index]
        matrices = []
        for n in [16, 32]:
            _, k, _ = read_k(prior / "result/native" / case / str(n) / model["unit"])
            matrices.append(k[index])
        for n in [32, 64]:
            folder = destination / str(n)
            folder.mkdir()
            for name in ["fort.10", "fort.21"]:
                shutil.copy2(raw / case / name, folder / name)
            deck = single_energy_deck(
                (raw / case / "rsolve.inp").read_bytes(), contract["requested_energy_literal"], n
            )
            (folder / "rsolve.inp").write_bytes(deck)
            stage = {"subranges": n, "completed": False, "original_exit": None, "timed_out": False}
            result["stages"].append(stage)
            save(destination / "report.json", result)
            begin = time.monotonic()
            with (
                (folder / "rsolve.inp").open("rb") as inp,
                (folder / "rsolve.out").open("xb") as log,
            ):
                with subprocess.Popen(
                    [str(engine)], cwd=folder, stdin=inp, stdout=log, stderr=subprocess.STDOUT
                ) as process:
                    try:
                        stage["original_exit"] = process.wait(
                            timeout=contract["resources"]["native_stage_timeout_seconds"]
                        )
                    except subprocess.TimeoutExpired:
                        process.kill()
                        stage.update(original_exit=process.wait(), timed_out=True)
            stage["wall_seconds"] = time.monotonic() - begin
            save(destination / "report.json", result)
            if stage["original_exit"] != 0 or stage["timed_out"]:
                raise ValueError("Native slice execution failed")
            if (folder / "rsolve.inp").read_bytes() != deck:
                raise ValueError("Native deck changed during execution")
            for name in ["fort.10", "fort.21"]:
                if digest(folder / name) != contract["input_sha256"][case + "/" + name]:
                    raise ValueError("Native input changed during execution")
            if not re.search(
                rf"PROPAGATED ACROSS\s+{n} SUBRANGES", (folder / "rsolve.out").read_text()
            ):
                raise ValueError("Printed native subrange count differs")
            energies, k, header = read_k(folder / model["unit"])
            if (
                k.shape != (1, 6, 6)
                or energies[0] != expected
                or any(
                    header[key] != model["header"][key]
                    for key in [
                        "target_components",
                        "boundary_channels",
                        "symmetry_charge",
                        "geometry_mass",
                    ]
                )
            ):
                raise ValueError("Single-energy identity or interface differs")
            stage.update(
                completed=True, differences=differences(k[0], matrices[-1]), k=k[0].tolist()
            )
            budget = (
                contract["slice_reproduction_budget"]
                if n == 32
                else contract["matrix_and_phase_absolute_budget"]
            )
            stage["verdict"] = all(value <= budget for value in stage["differences"].values())
            if not stage["verdict"]:
                result["verdict"] = False
                break
            if n == 64:
                matrices.append(k[0])
            save(destination / "report.json", result)
        if len(matrices) == 3:
            result.update(
                verdict=True,
                native_changes=[differences(matrices[i], matrices[i + 1]) for i in range(2)],
            )
    except (ValueError, RuntimeError, OSError, AssertionError, np.linalg.LinAlgError) as error:
        result["failure"] = f"{type(error).__name__}: {error}"
    save(destination / "report.json", result)
    return result


def independent_control(task):
    raw, output, case, configuration, contract = task
    step, endpoint, order = configuration
    result = {
        "case": case,
        "step_bohr": step,
        "endpoint_bohr": endpoint,
        "tail_order": order,
        "verdict": None,
    }
    folder = output / "independent" / case
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"h{step}-r{endpoint}-p{order}.json"
    save(path, result)
    try:
        model = load_model(raw / case)
        energy = model["energies"][contract["energy_index"]]
        q = surface_matrix(model["poles"], model["amplitudes"], model["ground"] + energy / 2)
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
                1e-10,
            )
        s = cayley(k)
        invariants = algebra(k)
        result.update(
            k=k.tolist(),
            s_real=s.real.tolist(),
            s_imag=s.imag.tolist(),
            diagnostics=diagnostics,
            energy_rydberg=float(energy),
            wall_seconds=time.monotonic() - begin,
            invariants=invariants,
            verdict=all(
                value <= contract["invariant_absolute_budget"] for value in invariants.values()
            ),
        )
    except (
        ValueError,
        RuntimeError,
        OSError,
        AssertionError,
        np.linalg.LinAlgError,
        LinAlgWarning,
    ) as error:
        result["failure"] = f"{type(error).__name__}: {error}"
    save(path, result)
    return result


def evaluate(old_rows, new_rows, reference, contract):
    rows = old_rows + new_rows
    keys = [(row["step_bohr"], row["endpoint_bohr"], row["tail_order"]) for row in rows]
    expected = [(0.025, 400.0, 8), (0.0125, 400.0, 8)] + [
        tuple(x) for x in contract["independent_configurations"]
    ]
    if len(set(keys)) != len(keys) or set(keys) != set(expected):
        raise ValueError("Missing, duplicate or undeclared independent control")
    matrices = {key: np.asarray(row["k"]) for key, row in zip(keys, rows, strict=True)}
    families = {
        "step": [(h, 400.0, 8) for h in [0.025, 0.0125, 0.00625]],
        "endpoint": [(0.0125, r, 8) for r in [200.0, 400.0, 800.0]],
        "order": [(0.0125, 400.0, p) for p in [4, 6, 8]],
    }
    metrics = {
        name + "_changes": [
            differences(matrices[sequence[i]], matrices[sequence[i + 1]]) for i in range(2)
        ]
        for name, sequence in families.items()
    }
    checks = {}
    for name in families:
        checks.update(
            {
                name + "_" + k: value <= contract["matrix_and_phase_absolute_budget"]
                for k, value in metrics[name + "_changes"][-1].items()
            }
        )
    metrics["native_comparison"] = differences(matrices[0.00625, 400.0, 8], reference)
    checks.update(
        {
            "native_" + k: value <= contract["matrix_and_phase_absolute_budget"]
            for k, value in metrics["native_comparison"].items()
        }
    )
    coarse, fine = metrics["step_changes"]
    checks["spatial_trend"] = fine["s"] <= max(coarse["s"] / 8, 1e-8)
    values, vectors = np.linalg.eigh(reference)
    vector = vectors[:, np.argmax(abs(values))]
    metrics["all_native_comparisons"] = [
        {
            "configuration": key,
            **differences(k, reference),
            "signed_large_channel_delta": float(vector @ (k - reference) @ vector),
            "coordinate_diagnostic": coordinate_diagnostic(k, reference),
        }
        for key, k in matrices.items()
    ]
    for name in ["unitarity", "reciprocity"]:
        checks[name] = all(
            algebra(k)[name] <= contract["invariant_absolute_budget"] for k in matrices.values()
        )
    return {"verdict": all(checks.values()), "checks": checks, "metrics": metrics}


def run(args, contract, report):
    for name, want in contract["input_sha256"].items():
        if digest(args.raw / name) != want:
            raise ValueError(f"Input digest differs: {name}")
    for name, want in contract["prior_sha256"].items():
        if digest(args.prior / name) != want:
            raise ValueError(f"Saved reference digest differs: {name}")
    if digest(args.engine) != contract["engine_rsolve_sha256"]:
        raise ValueError("Engine digest differs")
    if digest(args.prerequisite) != contract["public_prerequisite_sha256"]:
        raise ValueError("Public prerequisite digest differs")
    previous = json.loads(args.prerequisite.read_text())
    if (
        previous["public_reconstruction_passed"] is not True
        or previous["reconstruction_verdict"] is not False
    ):
        raise ValueError("Public prerequisite verdict differs")

    def receive(stage, row):
        report[stage].append(row)
        save(args.output / "report.json", report)

    with ProcessPoolExecutor(max_workers=4) as pool:
        tasks = [
            (args.raw, args.prior, args.output, args.engine, case, contract)
            for case in contract["cases"]
        ]
        verdict = gated_pool(pool, native_case, tasks, 2, lambda row: receive("native", row))
        report["stage_verdicts"]["native"] = verdict
        if verdict is not True:
            report["experiment_verdict"] = verdict
            return
        tasks = [
            (args.raw, args.output, case, config, contract)
            for config in contract["independent_configurations"]
            for case in contract["cases"]
        ]
        verdict = gated_pool(
            pool, independent_control, tasks, 4, lambda row: receive("independent", row)
        )
        report["stage_verdicts"]["independent_invariants"] = verdict
        if verdict is not True:
            report["experiment_verdict"] = verdict
            return
    for case in contract["cases"]:
        old = json.loads((args.prior / "result/independent" / case / "018.json").read_text())
        new = [row for row in report["independent"] if row["case"] == case]
        native = next(row for row in report["native"] if row["case"] == case)
        report["cases"].append(
            {
                "case": case,
                **evaluate(old["rows"], new, np.asarray(native["stages"][-1]["k"]), contract),
            }
        )
    report["experiment_verdict"] = all(row["verdict"] for row in report["cases"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["raw", "prior", "engine", "prerequisite", "contract", "output"]:
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text())
    args.output.mkdir()
    report = {
        "contract": contract,
        "native": [],
        "independent": [],
        "cases": [],
        "stage_verdicts": {},
        "experiment_verdict": None,
        "original_reconstruction_verdict": False,
        "original_reconstruction_exits": [1, 1],
        "native_full_k_reconstruction_verdict": False,
        "co_complex_pole_verdict": None,
        "new_inner_solves": 0,
        "full_sweep_released": False,
    }
    begin = time.monotonic()
    save(args.output / "report.json", report)
    try:
        run(args, contract, report)
    except (ValueError, RuntimeError, OSError, AssertionError, np.linalg.LinAlgError) as error:
        report.update(experiment_verdict=None, failure=f"{type(error).__name__}: {error}")
    finally:
        report["wall_seconds"] = time.monotonic() - begin
        save(args.output / "report.json", report)
    if report["experiment_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
