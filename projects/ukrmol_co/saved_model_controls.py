"""Read retained CO model differences and extraction sensitivity without native solves."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import re
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.a2_scattering import phases
from projects.ukrmol_co.phase_fit import fit_phase, predict_phase
from projects.ukrmol_co.sparse_comparison import compare_phases, read_run
from projects.ukrmol_co.sparse_pilot import audit_sector, target_map


def extraction_control(table: np.ndarray, seed: dict, windows: list, contract: dict) -> dict:
    """Keep every declared background/window/start and predict five unused point folds."""
    settings = contract["extraction"]
    conversion = contract["hartree_to_ev"]
    rows = []
    for low, high in windows:
        mask = (table[:, 0] >= low - 5e-5) & (table[:, 0] <= high + 5e-5)
        energy = table[mask, 0] * contract["native_rydberg_per_requested_ev"] / 2
        eta = table[mask, 1]
        for order in settings["background_orders"]:
            row = {
                "window_requested_ev": [low, high],
                "nback": order,
                "starts": [],
                "folds": [],
                "valid": False,
            }
            rows.append(row)
            try:
                for displacement, factor in settings["starts"]:
                    fit = fit_phase(
                        energy,
                        eta,
                        order,
                        seed["energy_hartree"] + displacement * np.ptp(energy),
                        seed["width_hartree"] * factor,
                    )
                    row["starts"].append(fit)
                if any(
                    not f["optimizer_success"] or f["at_parameter_bound"] for f in row["starts"]
                ):
                    continue
                best = min(row["starts"], key=lambda f: f["residual_sum_squares_rad2"])
                row["best_fit"] = best
                for fold in range(settings["held_point_folds"]):
                    held = np.arange(len(energy)) % settings["held_point_folds"] == fold
                    fitted = fit_phase(
                        energy[~held],
                        eta[~held],
                        order,
                        best["position_hartree"],
                        best["full_width_hartree"],
                    )
                    difference = (
                        predict_phase(energy[held], fitted) - eta[held] + np.pi / 2
                    ) % np.pi - np.pi / 2
                    row["folds"].append(
                        {
                            "fold": fold,
                            "held_indices": np.flatnonzero(held).tolist(),
                            "fit": fitted,
                            "maximum_phase_error_rad_mod_pi": float(max(abs(difference))),
                        }
                    )
                row["valid"] = all(
                    f["fit"]["optimizer_success"] and not f["fit"]["at_parameter_bound"]
                    for f in row["folds"]
                )
            except (ValueError, np.linalg.LinAlgError) as error:
                row["error"] = str(error)
    fits = [
        f
        for row in rows
        for f in row["starts"]
        if f["optimizer_success"] and not f["at_parameter_bound"]
    ]
    result = {
        "rows": rows,
        "extraction_stability_verdict": None,
        "scope": "Phase representation sensitivity; no pole or statistical uncertainty claim",
    }
    if not fits:
        return result
    centers = np.asarray([f["position_hartree"] for f in fits]) * conversion
    widths = np.asarray([f["full_width_hartree"] for f in fits]) * conversion
    gates = settings["gates"]
    width_budget = max(
        gates["absolute_width_floor_ev"],
        gates["maximum_full_width_span_fraction_of_median"] * np.median(widths),
    )
    held_max = max(
        (f["maximum_phase_error_rad_mod_pi"] for r in rows for f in r["folds"]), default=None
    )
    result.update(
        position_range_ev=[float(min(centers)), float(max(centers))],
        full_width_range_ev=[float(min(widths)), float(max(widths))],
        position_span_ev=float(np.ptp(centers)),
        full_width_span_ev=float(np.ptp(widths)),
        full_width_span_fraction_of_median=float(np.ptp(widths) / np.median(widths)),
        full_width_budget_ev=float(width_budget),
        maximum_held_phase_error_rad_mod_pi=held_max,
    )
    per_order = {}
    for order in settings["background_orders"]:
        selected = [
            f
            for row in rows
            if row["nback"] == order
            for f in row["starts"]
            if f["optimizer_success"] and not f["at_parameter_bound"]
        ]
        if selected:
            per_order[str(order)] = {
                "position_span_ev": float(
                    np.ptp([f["position_hartree"] for f in selected]) * conversion
                ),
                "full_width_span_ev": float(
                    np.ptp([f["full_width_hartree"] for f in selected]) * conversion
                ),
                "maximum_held_phase_error_rad_mod_pi": max(
                    (
                        f["maximum_phase_error_rad_mod_pi"]
                        for r in rows
                        if r["nback"] == order
                        for f in r["folds"]
                    ),
                    default=None,
                ),
            }
    result["per_order_window_diagnosis"] = per_order
    if all(row["valid"] for row in rows):
        result["extraction_stability_verdict"] = bool(
            np.ptp(centers) <= gates["maximum_position_span_ev"]
            and np.ptp(widths) <= width_budget
            and held_max <= gates["maximum_held_phase_error_rad_mod_pi"]
        )
    return result


def read_model(root: Path, model: dict, contract: dict) -> tuple[dict, dict]:
    """Read the original CAS11 record directly, preserving its historical output schema."""
    path = root / model["directory"]
    settings = contract | model
    if model["record"] == "pilot.json":
        return read_run(path, model["roots"], settings)
    saved = json.loads((path / model["record"]).read_text())
    resources = json.loads((path / "resources.json").read_text())
    if resources["exit_code"] != 0 or not saved["all_stages_complete"]:
        raise ValueError("Original CAS11 run failed or incomplete")
    result, tables = {"resources": resources, "sectors": {}}, {}
    for sector in contract["sectors"]:
        stages = saved["sectors"][sector]["stages"]
        if [(s["name"], s["original_exit"]) for s in stages] != [
            (name, 0) for name in ("scatci", "rsolve", "eigenp", "reson")
        ] or any(s["timed_out"] for s in stages):
            raise ValueError("Original CAS11 stage failed")
        directory = path / sector
        raw = audit_sector(
            directory,
            model["roots"],
            model["dimension"],
            model["continuum_coefficients"],
            model["boundary_channels"],
        )
        telemetry = json.loads((directory / "eps-telemetry.jsonl").read_text())
        if telemetry["which"] != contract["raw_checks"]["native_selection_which"]:
            raise ValueError("CAS11 native selection changed")
        result["sectors"][sector] = raw
        tables[sector] = phases(directory)
    delta = float(
        np.max(abs((tables["B1"][:, 1] - tables["B2"][:, 1] + np.pi / 2) % np.pi - np.pi / 2))
    )
    if delta > contract["raw_checks"]["pi_phase_rad_mod_pi"]:
        raise ValueError("CAS11 Pi symmetry failed")
    result["pi_phase_difference_rad_mod_pi"] = delta
    return result, tables


def equilibrium_inputs(root: Path, contract: dict) -> dict:
    """Resolve own-model target IDs and actual preparation rather than target-only defaults."""
    targets = contract["equilibrium_targets"]
    paths = {key: root / value for key, value in targets.items()}
    configs = {key: json.loads((path / "config.json").read_text()) for key, path in paths.items()}
    common = (
        "bond_length",
        "basis",
        "model",
        "orbitals",
        "frozen_orbitals",
        "sa_singlet_roots",
        "sa_triplet_roots",
        "target_states_used",
        "target_singlet_roots",
        "target_triplet_roots",
        "radius",
        "maxl",
        "deletion_threshold",
        "propagation_radius",
        "virtual_orbitals",
    )
    for key in common:
        if configs["eq11"][key] != configs["eq12_preparation"][key]:
            raise ValueError(f"Equilibrium preparation mismatch: {key}")
    maps = {}
    for name, active in (("eq11", [5, 3, 3, 0]), ("eq12", [6, 3, 3, 0])):
        target = json.loads((paths[name] / "target.json").read_text())
        resources = json.loads((paths[name] / "resources.json").read_text())
        imported = json.loads((paths[name] / "result.json").read_text())
        if not target["converged"] or resources["exit_code"] != 0:
            raise ValueError("Equilibrium target prerequisite failed")
        if (
            imported["target_state_max_energy_difference_hartree"] > 1e-7
            or abs(imported["target_dipole_difference_au"]) > 1e-5
        ):
            raise ValueError("Original target energy/dipole import failed")
        if configs[name]["active_orbitals"] != active:
            raise ValueError("Equilibrium active-space model changed")
        ground = min(s["energy_hartree"] for s in target["states"])
        first_excited = min(
            s["energy_hartree"] - ground
            for s in target["states"]
            if s["energy_hartree"] > ground + 1e-7
        )
        if first_excited <= 5 * contract["native_rydberg_per_requested_ev"] / 2:
            raise ValueError("Equilibrium grid reaches an excited target threshold")
        directory = root / contract["models"][name]["directory"]
        mapping = target_map(
            (paths[name] / "output/CO/geom1/prop.out").read_text(), target["states"]
        )
        for sector in contract["sectors"]:
            own = target_map((directory / sector / "fort.24").read_text(), target["states"])
            text = (directory / sector / "scatci.inp").read_text()
            match = re.search(r"\bidtarg\s*=\s*([\d,\s]+)", text)
            actual = list(map(int, re.findall(r"\d+", match[1]))) if match else []
            if own != mapping or actual != mapping:
                raise ValueError("Own-model native target map changed")
        maps[name] = mapping
    for name in ("eq12", "eq12refined"):
        for sector in contract["sectors"]:
            for stage in ("rsolve", "eigenp", "reson"):
                left = root / contract["models"]["eq11"]["directory"] / sector / f"{stage}.inp"
                right = root / contract["models"][name]["directory"] / sector / f"{stage}.inp"
                if left.read_bytes() != right.read_bytes():
                    raise ValueError("Equilibrium actual outer deck differs")
    return {
        "common_preparation_settings": {key: configs["eq11"][key] for key in common},
        "own_target_maps": maps,
        "target_only_vs_preparation_differences": {
            key: [configs["eq12"].get(key), configs["eq12_preparation"].get(key)]
            for key in configs["eq12"].keys() | configs["eq12_preparation"].keys()
            if configs["eq12"].get(key) != configs["eq12_preparation"].get(key)
        },
        "execution_scope": (
            "CAS11 selected run used4 ranks/release engine; CAS12 used8 ranks/audited engine. "
            "Existing exact CAS12 engine reproduction supports consistency; no matched CAS11 "
            "engine rerun. Target diagonalizers differ; original import/action gates apply."
        ),
    }


def study(root: Path, contract: dict) -> dict:
    started = time.monotonic()
    if json.loads((root / "contract.json").read_text()) != contract:
        raise ValueError("Fetched and declared study contracts differ")
    manifest = json.loads((root / "manifest.json").read_text())
    for name, row in manifest.items():
        path = root / name
        if (
            path.stat().st_size != row["bytes"]
            or hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]
        ):
            raise ValueError(f"Fetched checksum failed: {name}")
    inputs = equilibrium_inputs(root, contract)
    raw, tables = {}, {}
    for name, model in contract["models"].items():
        raw[name], tables[name] = read_model(root, model, contract)
    comparisons = {}
    for candidate in ("eq12", "eq12refined"):
        comparisons[candidate] = {}
        for sector in contract["sectors"]:
            row = compare_phases(tables["eq11"][sector], tables[candidate][sector], contract)
            row["observable_insensitivity_verdict"] = row.pop("root_pair_stability_verdict")
            if any(
                not raw[m]["sectors"][sector]["native_resonances"]
                or raw[m]["sectors"][sector]["native_maxfit_diagnostic"]
                for m in ("eq11", candidate)
            ):
                row["observable_insensitivity_verdict"] = None
            row["maximum_phase_requested_ev"] = float(
                tables[candidate][sector][
                    np.argmax(
                        abs(
                            (
                                tables[candidate][sector][:, 1]
                                - tables["eq11"][sector][:, 1]
                                + np.pi / 2
                            )
                            % np.pi
                            - np.pi / 2
                        )
                    ),
                    0,
                ]
            )
            comparisons[candidate][sector] = row
    extraction = {}
    for name in contract["extraction"]["models"]:
        windows = contract["extraction"][
            "equilibrium_windows_requested_ev"
            if name.startswith("eq")
            else "stretched_windows_requested_ev"
        ]
        extraction[name] = {}
        for sector in contract["sectors"]:
            native = raw[name]["sectors"][sector]
            if len(native["native_resonances"]) != 1 or native["native_maxfit_diagnostic"]:
                extraction[name][sector] = {
                    "extraction_stability_verdict": None,
                    "reason": "Empty/multiple/MAXFIT native candidates",
                }
                continue
            extraction[name][sector] = extraction_control(
                tables[name][sector], native["native_resonances"][0], windows, contract
            )
    dense = {}
    for sector in contract["sectors"]:
        path = (
            root
            / contract["equilibrium_targets"]["eq11"]
            / f"output/CO/geom1/eigenph.doublet.{sector}"
        )
        dense[sector] = compare_phases(np.loadtxt(path), tables["eq11"][sector], contract)
    refinement = {
        s: compare_phases(tables["eq12"][s], tables["eq12refined"][s], contract)
        for s in contract["sectors"]
    }
    numerical_controls = all(
        row["root_pair_stability_verdict"] is True
        for group in (dense, refinement)
        for row in group.values()
    ) and all(
        raw[name]["sectors"][s]["native_resonances"]
        and not raw[name]["sectors"][s]["native_maxfit_diagnostic"]
        for name in ("eq11", "eq12", "eq12refined")
        for s in contract["sectors"]
    )
    if not numerical_controls:
        for group in comparisons.values():
            for row in group.values():
                row["observable_insensitivity_verdict"] = None
                row["unset_reason"] = (
                    "Saved numerical baseline/refinement prerequisite failed or unset"
                )
    return {
        "verified_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "input_manifest": manifest,
        "files_verified": len(manifest),
        "equilibrium_input_reading": inputs,
        "raw_reading": raw,
        "equilibrium_active_space_comparisons": comparisons,
        "cas11_selected_vs_dense": dense,
        "cas12_refinement": refinement,
        "numerical_sensitivity_controls_passed": numerical_controls,
        "extraction_controls": extraction,
        "reader_wall_seconds": time.monotonic() - started,
        "electronic_model_verdict": None,
        "spectral_truncation_verdict": None,
        "global_continuity_verdict": None,
        "full_sweep_released": False,
        "scope": (
            "Practical saved-output model/extraction sensitivity after accepted root controls; "
            "no absolute omission bound, full model convergence or resonance-pole identity"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("root", "contract", "output"):
        parser.add_argument(f"--{field}", type=Path, required=True)
    args = parser.parse_args()
    result = study(args.root, json.loads(args.contract.read_text()))
    result["source_sha256"] = {
        str(p): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (
            Path(__file__).relative_to(Path.cwd()),
            Path(__file__).with_name("phase_fit.py").relative_to(Path.cwd()),
            args.contract,
        )
    }
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
