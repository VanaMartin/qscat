"""One saved CAS12 B1 outgoing-pole contrast; unchanged numerical pole machinery."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from copy import deepcopy
from pathlib import Path

import numpy as np

from projects.ukrmol_co.compressed_a2_contrast_run import real_guard, save_json
from projects.ukrmol_co.compressed_a2_pole_contrast import (
    load_pole_model,
    pole_insensitivity,
    require_real_guard,
)
from projects.ukrmol_co.pole_character import screen_witnesses
from projects.ukrmol_co.saved_pole_character import digest, verify_inputs


def resolved_contract(execution, base, proposal, reference, guard_path, guard_digest):
    """Change model/neutral reference, not numerical methods, budgets or identity rules."""
    z = complex(*reference["pole_incident_hartree"])
    half = proposal["b1_box_halfwidth_hartree"]
    result = deepcopy(base)
    result.update(
        case="compressed/cas12/B1",
        model_view="compressed-cas12-native40-v1",
        raw_directory=execution["raw_directory"],
        real_guard_result=guard_path,
        real_guard_budgets={
            k: proposal[k]
            for k in [
                "real_s_absolute_budget",
                "real_invariant_absolute_budget",
                "real_phase_sum_budget_rad",
                "real_flux_absolute_budget",
            ]
        },
        search_box_incident_hartree=[z.real - half, z.real + half, z.imag - half, z.imag + half],
        seeds_incident_hartree=[
            [z.real, z.imag],
            [z.real + proposal["second_seed_delta_hartree"], z.imag],
        ],
        root_seed_refinement_delta_hartree=proposal["second_seed_delta_hartree"],
        root_evaluation_caps=proposal["root_caps"]["B1"],
        settings=proposal["settings"],
        claim_boundary="Finite CAS12/five B1 contrast; no Pi pair, electronic accuracy or curve",
    )
    for key in [
        "reference_energy_rydberg",
        "same_root_position_full_width_budget_ev",
        "contour_radius_hartree",
        "contour_nodes",
        "projection_columns",
        "projection_seed",
        "derivative_steps_hartree",
        "derivative_residue_relative_budget",
    ]:
        result[key] = proposal[key]
    result["resources"] = {
        k: proposal["future_resources"][k]
        for k in [
            "cpuset_cpus",
            "threads_per_worker",
            "contour_workers",
            "memory_gib",
            "per_map_seconds",
            "rhs_cap",
        ]
    }
    result["resources"]["total_map_cap"] = proposal["future_resources"]["b1_map_cap"] - 2
    for key in ["cpuset_cpus", "threads_per_worker", "contour_workers", "memory_gib"]:
        result["resources"][key] = execution["resources"][key]
    result["input_sha256"].update(execution["input_sha256"])
    result["input_sha256"][guard_path] = guard_digest
    return result


def run(root, source, output, execution):
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    report = dict(
        model_insensitivity_verdict=None,
        disposition="inconclusive",
        scientific_exit=2,
        native_calls=0,
        new_inner_solves=0,
        real_maps_reserved=0,
        pi_pair_verdict=None,
        physical_co_accuracy_verdict=None,
        width_curve_verdict=None,
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
        timed_out=False,
    )

    def save():
        report["wall_seconds"] = time.monotonic() - started
        save_json(output / "result.json", report)

    try:
        report["inputs_verified_before"] = verify_inputs(root, execution)
        for name, want in execution["source_sha256"].items():
            if digest(source / name) != want:
                raise ValueError(f"Changed frozen source prerequisite: {name}")

        def read(name):
            return json.loads((source / execution[name]).read_text())

        base, proposal, policy, verified = [
            read(k)
            for k in [
                "base_pole_contract",
                "numerical_proposal",
                "character_policy",
                "inner_verification",
            ]
        ]
        native = json.loads((root / execution["inner_result"]).read_text())
        reference = json.loads((root / execution["cas11_five_reference"]).read_text())
        if (
            verified["inner_export_verdict"] is not True
            or native["inner_export_verdict"] is not True
            or native["scientific_exit"] != 0
            or reference["disposition"] != "pass"
            or not reference["completed"]
        ):
            raise ValueError(
                "Independent native export and completed local CAS11/five reference required"
            )
        contract = resolved_contract(
            execution, base, proposal, reference, "not-yet-measured", "unset"
        )
        model = load_pole_model(root / execution["raw_directory"] / contract["case"], contract)
        for name, want in verified["array_sha256"].items():
            if hashlib.sha256(np.asarray(model[name], dtype="<f8").tobytes()).hexdigest() != want:
                raise ValueError(f"Own-model array differs from independent native replay: {name}")
        if model["interface"]["channel_labels"][:6] != policy["open_labels"]["B1"]:
            raise ValueError("Own open-channel identity differs")
        save()
        guard = real_guard(model, output, contract["case"], proposal, report)
        guard_name = execution["output_relative"] + "/check/guard.json"
        contract = resolved_contract(
            execution, base, proposal, reference, guard_name, digest(guard)
        )
        require_real_guard(root, contract)
        frozen = output / "contract.resolved.json"
        save_json(frozen, contract)
        remaining = execution["per_case_hard_cap_seconds"] - (time.monotonic() - started)
        if remaining <= 0:
            raise TimeoutError("Hard cap exhausted before complex maps")
        command = [
            sys.executable,
            "-m",
            "projects.ukrmol_co.compressed_pole_reference",
            "--root",
            str(root),
            "--output",
            str(output / "pole"),
            "--contract",
            str(frozen),
        ]
        with (output / "pole.log").open("x") as log:
            process = subprocess.Popen(
                command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True
            )
            try:
                process.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                report["timed_out"] = True
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        report["original_pole_command_exit"] = process.returncode
        pole = json.loads((output / "pole/result.json").read_text())
        report["original_pole_scientific_exit"] = pole["scientific_exit"]
        if (
            report["timed_out"]
            or process.returncode != 0
            or pole["finite_reference_pole_verdict"] is not True
            or pole["scientific_exit"] != 0
        ):
            raise ValueError("Numerical pole prerequisite failed; model sensitivity inconclusive")
        z = complex(*pole["pole_incident_hartree"])
        residues = []
        for name in ["moments-0", "moments-1", "derivative-0", "derivative-1"]:
            with np.load(output / f"pole/{name}.npz", allow_pickle=False) as data:
                residues.append(data["residue"])
        screen = screen_witnesses(
            residues, policy["open_labels"]["B1"], -2 * z.imag, policy["character"]
        )
        decision = pole_insensitivity(
            z, complex(*reference["pole_incident_hartree"]), screen, reference["screen"], proposal
        )
        disposition = decision["disposition"]
        report.update(
            pole_incident_hartree=pole["pole_incident_hartree"],
            screen=screen,
            contrast=decision,
            disposition=disposition,
            model_insensitivity_verdict=True
            if disposition == "pass"
            else False
            if disposition == "reject"
            else None,
            scientific_exit=0 if disposition == "pass" else 1 if disposition == "reject" else 2,
        )
    except Exception as error:
        report["error"] = repr(error)
    finally:
        pole = output / "pole/result.json"
        report["pole_maps_reserved"] = (
            len(json.loads(pole.read_text())["evaluations"]) if pole.is_file() else 0
        )
        report["new_maps_reserved"] = report["real_maps_reserved"] + report["pole_maps_reserved"]
        try:
            report["inputs_verified_after"] = verify_inputs(root, execution)
            if (
                report["new_maps_reserved"] > execution["total_map_cap"]
                or time.monotonic() - started > execution["per_case_hard_cap_seconds"]
            ):
                raise ValueError("Frozen map/time cap exceeded")
        except Exception as error:
            report.update(
                integrity_or_cap_error=repr(error),
                scientific_exit=2,
                disposition="inconclusive",
                model_insensitivity_verdict=None,
            )
        for name in ["memory.peak", "cpu.stat"]:
            path = Path("/sys/fs/cgroup") / name
            if path.is_file():
                report["cgroup_" + name] = path.read_text()
        save()
    return report["scientific_exit"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "source", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(
        run(args.root, args.source, args.output, json.loads(args.contract.read_text()))
    )
