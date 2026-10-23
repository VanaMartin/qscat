"""Finite B2 adapter around the frozen compressed pole machinery, not a new solver."""

from __future__ import annotations

import argparse
import json
import resource
import time
from copy import deepcopy
from pathlib import Path

import numpy as np

from projects.ukrmol_co.compressed_pole_reference import run as run_pole
from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.native_full_reconstruction import load_model
from projects.ukrmol_co.outgoing_map import outgoing_map, physical_s
from projects.ukrmol_co.pole_character import compare_pair, screen_witnesses
from projects.ukrmol_co.saved_matrices import scattering_matrix
from projects.ukrmol_co.saved_pole_character import digest, verified_labels, verify_inputs


def b2_contract(policy, base, execution):
    """Only translate the already frozen pair allocation to the existing pole API."""
    result = deepcopy(base)
    z = complex(*policy["b1_pole_incident_hartree"])
    half = policy["b2_box_halfwidth_hartree"]
    result.update(
        case="compressed/nine/B2",
        search_box_incident_hartree=[z.real - half, z.real + half, z.imag - half, z.imag + half],
        seeds_incident_hartree=[
            [z.real, z.imag],
            [z.real + policy["b2_second_seed_delta_hartree"], z.imag],
        ],
        settings=policy["b2_settings"],
        root_evaluation_caps=policy["b2_root_evaluation_caps"],
        root_seed_refinement_delta_hartree=policy["b2_second_seed_delta_hartree"],
        resources=deepcopy(policy["b2_resources"]),
        claim_boundary=(
            "One witnessed B2 counterpart of the compressed finite reference, "
            "not absolute CO accuracy, exact Lambda or full-S equivalence"
        ),
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
        result[key] = policy[key]
    result["input_sha256"].update(policy["input_sha256"])
    result["input_sha256"].update(execution["input_sha256"])
    result["resources"]["total_map_cap"] -= len(policy["b2_development_indices"])
    if result["case"] not in policy["cases"] or result["resources"]["total_map_cap"] != 65:
        raise ValueError("B2 translation leaves frozen67-map allocation")
    return result


def real_checks(actual, reference, invariants, policy):
    error = float(np.max(abs(actual - reference)))
    phase = float(abs(np.angle(np.linalg.slogdet(actual)[0] / np.linalg.slogdet(reference)[0])) / 2)
    budget = policy["real_invariant_absolute_budget"]
    return dict(
        s_absolute_error=error,
        phase_sum_error_mod_pi_rad=phase,
        checks=dict(
            s=error <= policy["real_s_absolute_budget"],
            phase=phase <= policy["real_phase_sum_budget_rad"],
            flux=invariants["flux_defect"] <= policy["real_flux_absolute_budget"],
            wronskian=invariants["wronskian_residual"] <= budget,
            unitarity=invariants["unitarity_defect"] <= budget,
            reciprocity=invariants["reciprocity_defect"] <= budget,
        ),
    )


def run(root, source, output, execution):
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    report = dict(
        paired_verdict=None,
        paired_disposition="inconclusive",
        scientific_exit=2,
        real_rows=[],
        new_propagations=0,
        native_calls=0,
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
        physical_co_accuracy_verdict=None,
        width_curve_verdict=None,
    )

    def save():
        usage = resource.getrusage(resource.RUSAGE_SELF)
        report.update(
            wall_seconds=time.monotonic() - started,
            parent_cpu_seconds=usage.ru_utime + usage.ru_stime,
            peak_parent_rss_platform_units=usage.ru_maxrss,
        )
        for name in ["memory.peak", "cpu.stat"]:
            path = Path("/sys/fs/cgroup") / name
            if path.is_file():
                report["cgroup_" + name] = path.read_text()
        temporary = output / "result.partial"
        temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        temporary.replace(output / "result.json")

    policy = None
    try:
        for key, sha in [
            ("policy_contract", "policy_sha256"),
            ("pole_contract", "pole_contract_sha256"),
        ]:
            if digest(source / execution[key]) != execution[sha]:
                raise ValueError("Changed frozen policy/base contract")
        policy = json.loads((source / execution["policy_contract"]).read_text())
        base = json.loads((source / execution["pole_contract"]).read_text())
        for name, want in execution["input_sha256"].items():
            if digest(root / name) != want:
                raise ValueError(f"Changed saved-screen prerequisite: {name}")
        report["inputs_verified_before"] = verify_inputs(root, policy)
        left = json.loads((root / execution["b1_screen"]).read_text())
        owner = json.loads((root / execution["b1_screen_execution"]).read_text())
        if (
            left["scientific_exit"] != 0
            or left["b1_character_disposition"] != "pass"
            or left["paired_verdict"] is not None
            or not left["b2_calculation_eligible"]
            or [
                owner[k]
                for k in ["original_command_exit", "original_container_exit", "original_owner_exit"]
            ]
            != [0, 0, 0]
            or not owner["source_integrity_after"]
            or owner["timed_out"]
            or owner["container_oom_killed"]
        ):
            raise ValueError("Saved B1 screen is not an authenticated eligible pass")
        labels = verified_labels(root, policy)
        contract = b2_contract(policy, base, execution)
        report.update(b1_screen=left["b1_screen"], labels=labels)
        (output / "b2-contract.resolved.json").write_text(json.dumps(contract, indent=2) + "\n")
        model = load_model(root / policy["raw_directory"] / contract["case"])
        if policy["b2_development_indices"] != [60, 80]:
            raise ValueError("Only two declared B2 development identities allowed")
        endpoint, angle, order, rtol = policy["b2_settings"][-1]
        for index in policy["b2_development_indices"]:
            energy = model["energies"][index]
            row = dict(index=index, completed=False)
            report["real_rows"].append(row)
            report["new_propagations"] += 1
            save()
            q = surface_matrix(model["poles"], model["amplitudes"], model["ground"] + energy / 2)
            y, a, diagnostics = outgoing_map(
                model["thresholds"],
                model["angular"],
                model["c2"],
                model["c3"],
                energy,
                18,
                endpoint,
                angle,
                order,
                rtol=rtol,
                amplitudes=True,
                timeout=policy["b2_resources"]["per_map_seconds"],
                max_nfev=policy["b2_resources"]["rhs_cap"],
            )
            actual, invariants = physical_s(q, y, a)
            reference = scattering_matrix(model["saved"][index])
            np.savez(
                output / f"real-{index}.npz",
                q=q,
                y=y,
                a=a,
                s=actual,
                reference=reference,
                energy_rydberg=energy,
            )
            row.update(
                completed=True,
                diagnostics=diagnostics,
                invariants=invariants,
                **real_checks(actual, reference, invariants, policy),
            )
            save()
            if not all(row["checks"].values()):
                report.update(
                    scientific_exit=1, stopped_at="B2 real-axis gate; no complex propagation"
                )
                return 1
        run_pole(root, output / "b2", contract)
        right = json.loads((output / "b2/result.json").read_text())
        report.update(
            b2_original_scientific_exit=right["scientific_exit"],
            new_propagations=2 + len(right["evaluations"]),
        )
        if right["scientific_exit"] != 0 or right["finite_reference_pole_verdict"] is not True:
            report.update(
                scientific_exit=right["scientific_exit"], stopped_at="B2 pole prerequisite"
            )
            return report["scientific_exit"]
        z = complex(*right["pole_incident_hartree"])
        residues = []
        for name in ["moments-0", "moments-1", "derivative-0", "derivative-1"]:
            with np.load(output / f"b2/{name}.npz", allow_pickle=False) as capture:
                residues.append(capture["residue"])
        screen = screen_witnesses(residues, labels["B2"], -2 * z.imag, policy["character"])
        paired = compare_pair(
            left["b1_screen"],
            screen,
            complex(*policy["b1_pole_incident_hartree"]),
            z,
            policy["character"],
            policy["same_root_position_full_width_budget_ev"],
        )
        disposition = paired["disposition"]
        report.update(
            b2_screen=screen,
            counterpart=paired,
            paired_disposition=disposition,
            paired_verdict={"pass": True, "reject": False, "inconclusive": None}[disposition],
            scientific_exit={"pass": 0, "reject": 1, "inconclusive": 2}[disposition],
            b2_pole_incident_hartree=right["pole_incident_hartree"],
        )
        return report["scientific_exit"]
    except Exception as error:
        report.update(error=repr(error), scientific_exit=2)
        return 2
    finally:
        try:
            if policy is not None:
                report["inputs_verified_after"] = verify_inputs(root, policy)
            for name, want in execution["input_sha256"].items():
                if digest(root / name) != want:
                    raise ValueError(f"Changed prerequisite after use: {name}")
        except Exception as error:
            report.update(
                integrity_error=repr(error),
                scientific_exit=2,
                paired_verdict=None,
                paired_disposition="inconclusive",
            )
        save()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "source", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    run(args.root, args.source, args.output, json.loads(args.contract.read_text()))
    raise SystemExit(json.loads((args.output / "result.json").read_text())["scientific_exit"])
