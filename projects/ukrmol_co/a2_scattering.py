"""Matched retained-A2 channel sensitivity, with a dense-observable baseline gate."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np

from projects.ukrmol_co.a2_channel_preflight import (
    SECTORS,
    expected_states,
    property_inventory,
    replace_assignment,
    stage,
)
from projects.ukrmol_co.analyze import native_resonances
from projects.ukrmol_co.run import run_profiled
from projects.ukrmol_co.sparse_comparison import compare_phases, read_run
from projects.ukrmol_co.sparse_comparison_owner import sha
from projects.ukrmol_co.sparse_pilot import audit_sector
from projects.ukrmol_co.sparse_scattering import boundary_spectrum, mpi_boundary_deck, selected_deck


def save(path: Path, value: dict) -> None:
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def channel_map(inventory: dict, a2_roots: int) -> list[int]:
    """Select explicit sector roots from the complete48-state property representation."""
    if a2_roots not in (5, 9):
        raise ValueError("Declared A2 channel counts are five and nine")
    labels = {label: int(i) for i, label in inventory["state_labels_by_id"].items()}
    return [
        labels[f"{spin}.{irrep}.{root}"]
        for spin, irrep in SECTORS
        for root in range(1, (a2_roots if irrep == "A2" else 5) + 1)
    ]


def channel_deck(text: str, inventory: dict, a2_roots: int, roots: int) -> str:
    """Change contraction counts/sequences and the corresponding native interface map."""
    counts = [5] * 6 + [a2_roots] * 2
    text = replace_assignment(text, "numtgt", ",".join(map(str, counts)))
    text = replace_assignment(
        text, "ntgtf", ",".join(str(i) for i, count in enumerate(counts, 1) for _ in range(count))
    )
    text = replace_assignment(
        text, "ntgts", ",".join(str(i) for count in counts for i in range(1, count + 1))
    )
    mapping = channel_map(inventory, a2_roots)
    return mpi_boundary_deck(
        selected_deck(text, roots, 1e-12, 1000),
        f"ntarg={len(mapping)},idtarg={','.join(map(str, mapping))},rmatr=18,",
    )


def input_gate(root: Path, contract: dict) -> dict:
    """Require the completed exact-frame target/property import before any scattering."""
    expected = expected_states(root, contract)
    preflight = root / contract["target_preflight_directory"]
    owner = json.loads((preflight.parent.parent / "execution.json").read_text())
    saved = json.loads((preflight / "pilot.json").read_text())
    if (
        owner["original_controller_exit"] != 0
        or saved["channel_import_verdict"] is not True
        or saved["unchanged_orbital_hashes_passed"] is not True
    ):
        raise ValueError("Completed A2 import prerequisite failed")
    inventory = property_inventory((preflight / "fort.24").read_text(), expected)
    if inventory["target_map"] != saved["property_inventory"]["target_map"]:
        raise ValueError("Fresh48-state property map disagrees with saved import")
    target = root / contract["target_reference_directory"]
    target_saved = json.loads((target / "target.json").read_text())
    if abs(inventory["ground_dipole_z_au"] - target_saved["ground_dipole_au"][2]) > 1e-5:
        raise ValueError("Native target ground dipole failed")
    if json.loads((target / "resources.json").read_text())["exit_code"] != 0:
        raise ValueError("Original dense reference failed")
    first_excitation = (
        min(
            value
            for values in expected.values()
            for value in values
            if value > expected["singlet.A1"][0] + 1e-7
        )
        - expected["singlet.A1"][0]
    )
    if first_excitation <= 5 * contract["native_rydberg_per_requested_ev"] / 2:
        raise ValueError("Requested grid reaches an excited target threshold")
    for sector in contract["sectors"]:
        retained = target / "output/CO/geom1" / f"doublet.{sector}"
        for name in ("fort.16", "fort.26", "fort.24"):
            original = target / "output/CO/geom1" / name
            if sha(retained / name) != sha(original):
                raise ValueError(f"Dense retained input differs from original frame: {name}")
        if sha(retained / "fort.16") != sha(preflight / "fort.16"):
            raise ValueError("Scattering and qualified target integrals differ")
    return inventory


def outer_decks(geometry: Path, directory: Path, sector: str) -> None:
    """Keep the original outer method; explicitly request the common99-point grid."""
    for name in ("rsolve", "eigenp", "reson"):
        text = (geometry / "inputs" / f"scattering.{name}.doublet.{sector}.inp").read_text()
        if name == "rsolve":
            text = replace_assignment(text, "nescat", "99")
            text = replace_assignment(text, "einc", "0.1,0.05")
        (directory / f"{name}.inp").write_text(text)
    (directory / "reson_message").write_text("")


def phases(directory: Path) -> np.ndarray:
    unit = int(re.search(r"luphso\s*=\s*(\d+)", (directory / "eigenp.inp").read_text())[1])
    result = np.loadtxt(directory / f"fort.{unit}", ndmin=2)
    if result.shape != (99, 2) or not np.isfinite(result).all():
        raise ValueError("Incomplete99-point phase output")
    np.testing.assert_allclose(result[:, 0], 0.1 + 0.05 * np.arange(99), atol=5e-5, rtol=0)
    return result


def verify_model_inputs(
    root: Path, work: Path, contract: dict, inventory: dict, model: str
) -> None:
    """Resolve actual inputs/decks against the frozen common frame and declared count."""
    geometry = root / contract["target_reference_directory"] / "output/CO/geom1"
    imported = root / contract["target_preflight_directory"]
    for sector in contract["sectors"]:
        directory = work / model / sector
        text = (geometry / "inputs" / f"scattering.scatci.doublet.{sector}.inp").read_text()
        expected_deck = channel_deck(
            text, inventory, 5 if model == "five" else 9, contract["scattering_roots"]
        )
        if (directory / "scatci.inp").read_text() != expected_deck:
            raise ValueError("Actual scattering deck differs from declared channel model")
        unit = int(re.search(r"megul\s*=\s*(\d+)", text)[1])
        for name, path in {
            "fort.16": geometry / f"doublet.{sector}" / "fort.16",
            "fort.24": imported / "fort.24",
            "fort.26": imported / "fort.26",
            f"fort.{unit}": geometry / f"doublet.{sector}" / f"fort.{unit}",
        }.items():
            if sha(directory / name) != sha(path):
                raise ValueError(f"Actual scattering input changed: {model}/{sector}/{name}")


def compute(root: Path, work: Path, contract: dict, model: str) -> None:
    inventory = input_gate(root, contract)
    geometry = root / contract["target_reference_directory"] / "output/CO/geom1"
    imported = root / contract["target_preflight_directory"]
    a2_roots = 9 if model == "nine" else 5
    report = {"model": model, "target_map": channel_map(inventory, a2_roots), "sectors": {}}
    for sector in contract["sectors"]:
        directory = work / sector
        directory.mkdir()
        outer_decks(geometry, directory, sector)
        native = {"stages": []}
        if model == "dense":
            collected = geometry.parent / "collected_scattering_data"
            inputs = {
                "fort.10": collected / "channels" / f"channels.geom1.doublet.{sector}",
                "fort.21": collected / "rmat_amplitudes" / f"ramps.geom1.doublet.{sector}",
            }
            poles, amplitudes = boundary_spectrum(inputs["fort.21"])
            if (
                len(poles) != contract["models"]["five"]["dimension"]
                or amplitudes.shape != (len(poles), contract["models"]["five"]["boundary_channels"])
                or np.any(np.diff(poles) < 0)
            ):
                raise ValueError("Saved dense boundary does not contain the full spectrum")
        else:
            text = (geometry / "inputs" / f"scattering.scatci.doublet.{sector}.inp").read_text()
            unit = int(re.search(r"megul\s*=\s*(\d+)", text)[1])
            inputs = {
                "fort.16": geometry / f"doublet.{sector}" / "fort.16",
                "fort.24": imported / "fort.24",
                "fort.26": imported / "fort.26",
                f"fort.{unit}": geometry / f"doublet.{sector}" / f"fort.{unit}",
            }
            (directory / "scatci.inp").write_text(
                channel_deck(text, inventory, a2_roots, contract["scattering_roots"])
            )
        hashes = {}
        for name, source in inputs.items():
            shutil.copy2(source, directory / name)
            hashes[name] = sha(directory / name)
        stages = (
            ("rsolve", "eigenp", "reson")
            if model == "dense"
            else ("scatci", "rsolve", "eigenp", "reson")
        )
        for name in stages:
            env = os.environ.copy()
            command = [f"/opt/ukrmolp/bin.double/{name}"]
            if name == "scatci":
                command = [
                    "mpiexec",
                    "-n",
                    str(contract["resources"]["mpi_ranks"]),
                    "/opt/ukrmolp/bin.double/mpi-scatci",
                    "scatci.inp",
                ]
                env.update(
                    LD_PRELOAD=str(
                        root / contract["telemetry_directory"] / "scattering_telemetry.so"
                    ),
                    UKRMOL_EPS_AUDIT=str(directory / "eps-telemetry.jsonl"),
                )
            stage(
                directory,
                native,
                name,
                command,
                directory / f"{name}.inp",
                contract["resources"]["native_stage_timeout_seconds"],
                env,
            )
        if model != "dense":
            settings = contract["models"][model]
            native.update(
                audit_sector(
                    directory,
                    contract["scattering_roots"],
                    settings["dimension"],
                    settings["continuum_coefficients"],
                    settings["boundary_channels"],
                )
            )
        native["input_sha256"] = hashes
        report["sectors"][sector] = native
        save(work / "pilot.json", report)
    pi = float(
        np.max(
            abs(
                (phases(work / "B1")[:, 1] - phases(work / "B2")[:, 1] + np.pi / 2) % np.pi
                - np.pi / 2
            )
        )
    )
    report.update(all_stages_complete=True, pi_phase_difference_rad_mod_pi=pi)
    save(work / "pilot.json", report)
    if pi > contract["raw_checks"]["pi_phase_rad_mod_pi"]:
        raise ValueError("Pi symmetry failed")


def compare(root: Path, work: Path, contract: dict, baseline_gate: bool) -> dict:
    """Compare observables across channel models; never require equal Hamiltonian poles."""
    inventory = input_gate(root, contract)
    candidate_name = "five" if baseline_gate else "nine"
    settings = contract | contract["models"][candidate_name]
    verify_model_inputs(root, work, contract, inventory, candidate_name)
    right_path = work / candidate_name
    saved = json.loads((right_path / "pilot.json").read_text())
    resources = json.loads((right_path / "resources.json").read_text())
    if not saved.get("all_stages_complete") or resources["exit_code"] != 0:
        raise ValueError("Candidate native run failed/incomplete")
    raw, right_phases = {}, {}
    expected_map = channel_map(inventory, 5 if baseline_gate else 9)
    if saved["target_map"] != expected_map:
        raise ValueError("Wrong explicit retained channel map")
    for sector in contract["sectors"]:
        directory = right_path / sector
        stages = saved["sectors"][sector]["stages"]
        if [(s["name"], s["original_exit"]) for s in stages] != [
            (n, 0) for n in ("scatci", "rsolve", "eigenp", "reson")
        ] or any(s["timed_out"] for s in stages):
            raise ValueError("Candidate original stage failed")
        raw[sector] = audit_sector(
            directory,
            contract["scattering_roots"],
            settings["dimension"],
            settings["continuum_coefficients"],
            settings["boundary_channels"],
        )
        telemetry = json.loads((directory / "eps-telemetry.jsonl").read_text())
        if telemetry["which"] != contract["raw_checks"]["native_selection_which"]:
            raise ValueError("Wrong native spectral selection")
        right_phases[sector] = phases(directory)
    if baseline_gate:
        dense = work / "dense"
        dense_saved = json.loads((dense / "pilot.json").read_text())
        if (
            not dense_saved["all_stages_complete"]
            or json.loads((dense / "resources.json").read_text())["exit_code"] != 0
        ):
            raise ValueError("Dense outer replay failed")
        left_phases = {s: phases(dense / s) for s in contract["sectors"]}
        for sector in contract["sectors"]:
            stages = dense_saved["sectors"][sector]["stages"]
            if [(s["name"], s["original_exit"]) for s in stages] != [
                (n, 0) for n in ("rsolve", "eigenp", "reson")
            ] or any(s["timed_out"] for s in stages):
                raise ValueError("Dense outer original stage failed")
            collected = (
                root
                / contract["target_reference_directory"]
                / "output/CO/collected_scattering_data"
            )
            for name, original in {
                "fort.10": collected / "channels" / f"channels.geom1.doublet.{sector}",
                "fort.21": collected / "rmat_amplitudes" / f"ramps.geom1.doublet.{sector}",
            }.items():
                if sha(dense / sector / name) != sha(original):
                    raise ValueError("Dense outer input changed")
            for name in ("rsolve", "eigenp", "reson"):
                if (dense / sector / f"{name}.inp").read_bytes() != (
                    right_path / sector / f"{name}.inp"
                ).read_bytes():
                    raise ValueError("Dense and selected baseline outer decks differ")
    else:
        prerequisite = json.loads((work / "baseline-qualification.json").read_text())
        if prerequisite["dense_baseline_verdict"] is not True:
            raise ValueError("Five-channel dense prerequisite failed or unset")
        verify_model_inputs(root, work, contract, inventory, "five")
        _, left_phases = read_run(
            work / "five", contract["scattering_roots"], contract | contract["models"]["five"]
        )
        for sector in contract["sectors"]:
            left, right = work / "five" / sector, right_path / sector
            for name in saved["sectors"][sector]["input_sha256"]:
                if sha(left / name) != sha(right / name):
                    raise ValueError("Matched target/integral/CONGEN input changed")
            for name in ("rsolve", "eigenp", "reson"):
                if (left / f"{name}.inp").read_bytes() != (right / f"{name}.inp").read_bytes():
                    raise ValueError("Matched outer deck changed")
    comparisons = {}
    for sector in contract["sectors"]:
        comparison = compare_phases(left_phases[sector], right_phases[sector], contract)
        verdict = comparison.pop("root_pair_stability_verdict")
        if raw[sector]["native_maxfit_diagnostic"] or not raw[sector]["native_resonances"]:
            verdict = None
        if baseline_gate:
            output = work / "dense" / sector / "reson.out"
            if "MAXFIT" in output.read_text().upper() or not native_resonances(output):
                verdict = None
        comparison["observable_budget_verdict"] = verdict
        comparisons[sector] = comparison
    verdicts = [r["observable_budget_verdict"] for r in comparisons.values()]
    for values in (left_phases, right_phases):
        pi = np.max(abs((values["B1"][:, 1] - values["B2"][:, 1] + np.pi / 2) % np.pi - np.pi / 2))
        if pi > contract["raw_checks"]["pi_phase_rad_mod_pi"]:
            raise ValueError("Independent Pi check failed")
    return {
        "dense_baseline_verdict" if baseline_gate else "a2_channel_sensitivity_verdict": None
        if None in verdicts
        else all(verdicts),
        "sectors": comparisons,
        "candidate_raw_reading": raw,
        "independent_raw_reading_passed": True,
        "root_pair_stability_verdict": None,
        "spectral_truncation_verdict": None,
        "global_continuity_verdict": None,
        "electronic_model_verdict": None,
        "scope": "Dense-observable qualification of stretched five-A2 baseline"
        if baseline_gate
        else "Five-versus-nine A2 sensitivity; channel convergence and pole identity remain open",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("check-inputs", "model", "baseline-gate", "compare"))
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--model", choices=("dense", "five", "nine"))
    parser.add_argument("--profiled", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    root, work = args.root.resolve(), args.work.resolve()
    contract = json.loads(args.contract.read_text())
    if args.action == "check-inputs":
        save(work / "channel-input-gate.json", input_gate(root, contract))
    elif args.action == "model":
        if args.model is None:
            parser.error("model action requires --model")
        if (
            args.model == "nine"
            and json.loads((work / "baseline-qualification.json").read_text())[
                "dense_baseline_verdict"
            ]
            is not True
        ):
            raise ValueError("Expanded channel run requires passing dense baseline")
        directory = work / args.model
        if args.profiled:
            compute(root, directory, contract, args.model)
        else:
            directory.mkdir()
            (directory / "scratch").mkdir()
            env = os.environ | {
                "TMPDIR": str(directory / "scratch"),
                "PSI_SCRATCH": str(directory / "scratch"),
            }
            raise SystemExit(
                run_profiled(
                    directory,
                    env,
                    [
                        "python3",
                        "-m",
                        "projects.ukrmol_co.a2_scattering",
                        *sys.argv[1:],
                        "--profiled",
                    ],
                )
            )
    else:
        baseline = args.action == "baseline-gate"
        result = compare(root, work, contract, baseline)
        save(
            work / ("baseline-qualification.json" if baseline else "independent-comparison.json"),
            result,
        )
        if (
            result["dense_baseline_verdict" if baseline else "a2_channel_sensitivity_verdict"]
            is not True
        ):
            raise SystemExit(1)


if __name__ == "__main__":
    main()
