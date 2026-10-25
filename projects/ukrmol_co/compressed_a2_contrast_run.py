"""Finite five/nine retained-A2 experiment owner; numerical methods are unchanged."""

from __future__ import annotations

import argparse
import json
import os
import resource
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.compressed_a2_pole_contrast import (
    contrast_contract,
    load_pole_model,
    pole_insensitivity,
    require_real_guard,
)
from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.outgoing_map import outgoing_map, physical_s
from projects.ukrmol_co.paired_pole_character import real_checks
from projects.ukrmol_co.pole_character import compare_pair, screen_witnesses
from projects.ukrmol_co.saved_matrices import scattering_matrix
from projects.ukrmol_co.saved_pole_character import digest, verified_labels, verify_inputs


def save_json(path, value):
    """Persist progress atomically without hiding partial scientific failures."""
    partial = path.with_suffix(".partial")
    partial.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    partial.replace(path)


def original_pass(owner):
    return (
        [
            owner[k]
            for k in ["original_command_exit", "original_container_exit", "original_owner_exit"]
        ]
        == [0, 0, 0]
        and owner["source_integrity_after"]
        and not owner["timed_out"]
        and not owner["container_oom_killed"]
    )


def load_prerequisites(root, source, execution):
    """Only authenticated saved passes can release this frozen contrast."""
    for name, want in execution["source_sha256"].items():
        if digest(source / name) != want:
            raise ValueError(f"Changed source contract/pointer: {name}")

    def read_source(name):
        return json.loads((source / execution[name]).read_text())

    proposal, base, policy = [
        read_source(key) for key in ["proposal", "pole_contract", "character_policy"]
    ]
    gate = json.loads((root / execution["view_result"]).read_text())
    view_owner = json.loads((root / execution["view_execution"]).read_text())
    pair = json.loads((root / execution["nine_pair_result"]).read_text())
    pair_owner = json.loads((root / execution["nine_pair_execution"]).read_text())
    if (
        gate["mapped_view_verdict"] is not True
        or gate["scientific_exit"] != 0
        or len(gate["cases"]) != 2
        or not all(row["completed"] and row["own_w_d_byte_exact"] for row in gate["cases"])
        or not original_pass(view_owner)
        or pair["paired_verdict"] is not True
        or pair["scientific_exit"] != 0
        or not original_pass(pair_owner)
    ):
        raise ValueError("Saved view/nine pair is not an authenticated pass")
    labels = verified_labels(root, policy)
    references = {}
    for sector in ["B1", "B2"]:
        pole = json.loads((root / execution["nine_pole_results"][sector]).read_text())
        screen = pair["b1_screen" if sector == "B1" else "b2_screen"]
        if (
            pole["finite_reference_pole_verdict"] is not True
            or pole["scientific_exit"] != 0
            or pole["case"] != "compressed/nine/" + sector
            or screen["disposition"] != "pass"
            or not all(pole["contour_witness"]["checks"].values())
            or not all(all(row["checks"].values()) for row in pole["derivative_witnesses"])
        ):
            raise ValueError("Saved nine pole/character prerequisite does not pass")
        references[sector] = dict(pole=pole["pole_incident_hartree"], screen=screen)
    return proposal, base, policy, labels, references


def real_guard(model, directory, case, proposal, progress):
    """Exactly two own-model real S comparisons, never route-holdout indices."""
    report = dict(case=case, scientific_exit=2, real_rows=[])
    path = directory / "guard.json"
    save_json(path, report)
    if proposal["real_development_indices"] != [60, 80]:
        raise ValueError("Only the frozen two development identities are allowed")
    endpoint, angle, order, rtol = proposal["settings"][-1]
    resources = proposal["future_resources"]
    for index in [60, 80]:
        row = dict(index=index, completed=False)
        report["real_rows"].append(row)
        save_json(path, report)
        progress["real_maps_reserved"] += 1
        energy = model["energies"][index]
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
            timeout=resources["per_map_seconds"],
            max_nfev=resources["rhs_cap"],
        )
        actual, invariants = physical_s(q, y, a)
        reference = scattering_matrix(model["saved"][index])
        np.savez(
            directory / f"real-{index}.npz",
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
            **real_checks(actual, reference, invariants, proposal),
        )
        save_json(path, report)
        if not all(row["checks"].values()):
            raise ValueError("Own-model real-axis prerequisite failed; no complex map")
    report["scientific_exit"] = 0
    save_json(path, report)
    return path


def run_sector(
    root,
    directory,
    execution,
    proposal,
    base,
    policy,
    labels,
    sector,
    reference,
    record,
    *,
    seed=None,
):
    """One sector's guards and unchanged pole CLI within its original hard cap."""
    directory.mkdir(exist_ok=False)
    started = time.monotonic()
    try:
        # Placeholder digest is never used for propagation; replace after measured guards.
        contract = contrast_contract(proposal, base, sector, "not-yet-measured", "unset", seed=seed)
        model = load_pole_model(root / proposal["raw_directory"] / contract["case"], contract)
        if model["interface"]["channel_labels"][:6] != labels:
            raise ValueError("Own-model open-channel identity changed")
        guard = real_guard(model, directory, contract["case"], proposal, record)
        guard_name = execution["output_relative"] + f"/check/{sector}/guard.json"
        contract = contrast_contract(proposal, base, sector, guard_name, digest(guard), seed=seed)
        contract["input_sha256"].update(execution["input_sha256"])
        require_real_guard(root, contract)
        resolved = directory / "contract.resolved.json"
        save_json(resolved, contract)
        resources = proposal["future_resources"]
        remaining = resources["per_sector_hard_cap_seconds"] - (time.monotonic() - started)
        if remaining <= 0:
            raise TimeoutError("Sector hard cap exhausted before pole work")
        command = [
            sys.executable,
            "-m",
            "projects.ukrmol_co.compressed_pole_reference",
            "--root",
            str(root),
            "--output",
            str(directory / "pole"),
            "--contract",
            str(resolved),
        ]
        with (directory / "pole.log").open("x") as log:
            process = subprocess.Popen(
                command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True
            )
            try:
                process.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                record["timed_out"] = True
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            record["original_pole_command_exit"] = process.returncode
        pole = json.loads((directory / "pole/result.json").read_text())
        record["original_pole_scientific_exit"] = pole["scientific_exit"]
        if (
            record["timed_out"]
            or process.returncode != 0
            or pole["scientific_exit"] != 0
            or pole["finite_reference_pole_verdict"] is not True
        ):
            raise ValueError("Pole prerequisite failed; model sensitivity is inconclusive")
        z = complex(*pole["pole_incident_hartree"])
        residues = []
        for name in ["moments-0", "moments-1", "derivative-0", "derivative-1"]:
            with np.load(directory / f"pole/{name}.npz", allow_pickle=False) as capture:
                residues.append(capture["residue"])
        screen = screen_witnesses(residues, labels, -2 * z.imag, policy["character"])
        decision = pole_insensitivity(
            z, complex(*reference["pole"]), screen, reference["screen"], proposal
        )
        record.update(
            pole_incident_hartree=pole["pole_incident_hartree"],
            screen=screen,
            contrast=decision,
            disposition=decision["disposition"],
            completed=True,
        )
    except Exception as error:
        record.update(error=repr(error), disposition="inconclusive")
    finally:
        path = directory / "pole/result.json"
        record["pole_maps_reserved"] = (
            len(json.loads(path.read_text())["evaluations"]) if path.is_file() else 0
        )
        record["new_maps_reserved"] = record["real_maps_reserved"] + record["pole_maps_reserved"]
        record["wall_seconds"] = time.monotonic() - started
        cap = proposal["future_resources"]["b1_map_cap" if sector == "B1" else "b2_map_cap"]
        if record["new_maps_reserved"] > cap:
            record.update(disposition="inconclusive", error="Frozen map cap exceeded")
        if record["wall_seconds"] > proposal["future_resources"]["per_sector_hard_cap_seconds"]:
            record.update(
                disposition="inconclusive", timed_out=True, error="Sector hard cap exceeded"
            )
        save_json(directory / "result.json", record)
    return record


def run(root, source, output, execution):
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    report = dict(
        model_insensitivity_verdict=None,
        disposition="inconclusive",
        scientific_exit=2,
        sectors=[],
        new_maps_reserved=0,
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
            new_maps_reserved=sum(row["new_maps_reserved"] for row in report["sectors"]),
        )
        for name in ["memory.peak", "cpu.stat"]:
            path = Path("/sys/fs/cgroup") / name
            if path.is_file():
                report["cgroup_" + name] = path.read_text()
        save_json(output / "result.json", report)

    try:
        report["inputs_verified_before"] = verify_inputs(root, execution)
        proposal, base, policy, labels, references = load_prerequisites(root, source, execution)
        for sector in ["B1", "B2"]:
            record = dict(
                sector=sector,
                completed=False,
                disposition="inconclusive",
                real_maps_reserved=0,
                pole_maps_reserved=0,
                new_maps_reserved=0,
                timed_out=False,
            )
            report["sectors"].append(record)
            save()
            seed = (
                complex(*report["sectors"][0]["pole_incident_hartree"]) if sector == "B2" else None
            )
            run_sector(
                root,
                output / sector,
                execution,
                proposal,
                base,
                policy,
                labels[sector],
                sector,
                references[sector],
                record,
                seed=seed,
            )
            save()
            if record["disposition"] == "inconclusive":
                report["stopped_at"] = sector + " prerequisite/identity/numerical gate"
                break
            if sector == "B2":
                left = report["sectors"][0]
                pair = compare_pair(
                    left["screen"],
                    record["screen"],
                    complex(*left["pole_incident_hartree"]),
                    complex(*record["pole_incident_hartree"]),
                    policy["character"],
                    proposal["same_root_position_full_width_budget_ev"],
                )
                report["five_pair"] = pair
                if pair["disposition"] != "pass":
                    report["stopped_at"] = "Five-model counterpart prerequisite"
                    break
            if record["disposition"] == "reject":
                report.update(
                    disposition="reject",
                    model_insensitivity_verdict=False,
                    scientific_exit=1,
                    stopped_at=sector + " resolved pole sensitivity",
                )
                break
        else:
            report.update(disposition="pass", model_insensitivity_verdict=True, scientific_exit=0)
    except Exception as error:
        report["error"] = repr(error)
    finally:
        try:
            report["inputs_verified_after"] = verify_inputs(root, execution)
        except Exception as error:
            report.update(
                integrity_error=repr(error),
                disposition="inconclusive",
                scientific_exit=2,
                model_insensitivity_verdict=None,
            )
        save()
    return report["scientific_exit"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "source", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    run(args.root, args.source, args.output, json.loads(args.contract.read_text()))
    raise SystemExit(json.loads((args.output / "result.json").read_text())["scientific_exit"])
