"""Import explicit A2 channels at fixed orbitals before declaring scattering."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import signal
import subprocess
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.run import run_profiled
from projects.ukrmol_co.sparse_comparison_owner import sha

SECTORS = [(spin, irrep) for irrep in ("A1", "B1", "B2", "A2") for spin in ("singlet", "triplet")]


def replace_assignment(text: str, name: str, value: str) -> str:
    """Replace exactly one scalar/list assignment, retaining its original comment."""
    pattern = rf"(\b{re.escape(name)}\s*=\s*)[^!\n/]+"
    result, count = re.subn(pattern, lambda m: m[1] + value + ", ", text)
    if count != 1:
        raise ValueError(f"Expected one native assignment: {name}")
    return result


def property_inventory(text: str, expected: dict[str, list[float]]) -> dict:
    """Resolve nonuniform target sectors and every multipole's own energy-ordered ID."""
    states, moments, header = [], {}, []
    for line in text.splitlines():
        tokens = line.split()
        if not tokens:
            continue
        if tokens[0] == "6":
            header.append(int(tokens[4]))
        elif tokens[0] == "5":
            states.append(
                (
                    int(tokens[1]),
                    tokens[-2].lower() + "." + tokens[-1],
                    float(tokens[8].replace("D", "E")),
                )
            )
        elif tokens[0] == "1":
            key = tuple(map(int, tokens[1:8]))
            value = float(tokens[8].replace("D", "E"))
            if key in moments or not np.isfinite(value):
                raise ValueError("Duplicate or nonfinite multipole")
            moments[key] = value
    total = sum(map(len, expected.values()))
    if header != [total] or [s[0] for s in states] != list(range(1, total + 1)):
        raise ValueError("Incomplete DENPROP header/state IDs")
    if not np.isfinite([s[2] for s in states]).all() or np.any(
        np.diff([s[2] for s in states]) < -1e-10
    ):
        raise ValueError("Nonfinite or unordered DENPROP energies")
    if set(s[1] for s in states) != set(expected):
        raise ValueError("Unexpected DENPROP spin/irrep")
    labels, mapping, maximum_error = {}, [], 0.0
    for spin, irrep in SECTORS:
        sector = f"{spin}.{irrep}"
        group = [s for s in states if s[1] == sector]
        if len(group) != len(expected[sector]):
            raise ValueError(f"Incomplete retained target sector: {sector}")
        error = float(np.max(abs(np.asarray([s[2] for s in group]) - expected[sector])))
        if error > 1e-7:
            raise ValueError(f"Native target energy gate failed: {sector}")
        maximum_error = max(maximum_error, error)
        for root, state in enumerate(group, 1):
            labels[state[0]] = f"{sector}.{root}"
            mapping.append(state[0])
    referenced = set()
    angular_ranks = set()
    for key in moments:
        left, _, right, _, kind, rank, projection = key
        if left not in labels or right not in labels or kind != 3 or rank not in (1, 2):
            raise ValueError("Multipole refers to an absent state or unexpected rank")
        if abs(projection) > rank:
            raise ValueError("Invalid multipole projection")
        referenced.update((left, right))
        angular_ranks.add(rank)
    if referenced != set(labels) or angular_ranks != {1, 2}:
        raise ValueError("Incomplete dipole/quadrupole coverage")
    ground_id = next(i for i, label in labels.items() if label == "singlet.A1.1")
    dipole = -moments[(ground_id, 0, ground_id, 0, 3, 1, 0)]
    return {
        "target_map": mapping,
        "state_labels_by_id": labels,
        "target_count": total,
        "maximum_energy_error_hartree": maximum_error,
        "multipole_record_count": len(moments),
        "ground_dipole_z_au": dipole,
    }


def expected_states(root: Path, contract: dict) -> dict[str, list[float]]:
    """Require already validated twelve-root actions on this exact frozen checkpoint."""
    run = root / contract["target_reference_directory"]
    config = json.loads((run / "config.json").read_text())
    for key, expected_config in contract["fixed_target_settings"].items():
        if config.get(key) != expected_config:
            raise ValueError(f"Frozen target setting changed: {key}")
    checkpoint = run / "output/CO/geom1/co.casscf.chk"
    target = json.loads((run / "target.json").read_text())
    geometry = run / "output/CO/geom1"
    if sha(geometry / "co.molden") != target["molden_sha256"] or not target["converged"]:
        raise ValueError("Frozen orbital export gate failed")
    reference = root / contract["independent_reference_directory"]
    expected = {}
    for spin, irrep in SECTORS:
        sector = f"{spin}.{irrep}"
        saved = json.loads((reference / f"{sector}.json").read_text())
        count = 9 if irrep == "A2" else 5
        if (
            saved["checkpoint_sha256"] != sha(checkpoint)
            or not saved["gate_passed"]
            or not all(saved["ci_converged"])
        ):
            raise ValueError(f"Independent fixed-orbital prerequisite failed: {sector}")
        for field in ("physical_residual_norms_hartree", "penalized_residual_norms_hartree"):
            values = np.asarray(saved[field][:count])
            if len(values) != count or not np.isfinite(values).all() or max(values) > 1e-9:
                raise ValueError(f"Independent action gate failed: {sector}")
        spins = np.asarray(saved["spin_square"][:count])
        if (
            len(spins) != count
            or not np.isfinite(spins).all()
            or max(abs(spins - (0 if spin == "singlet" else 2))) > 1e-6
        ):
            raise ValueError(f"Independent spin gate failed: {sector}")
        expected[sector] = saved["energy_hartree"][:count]
        original = [
            s["energy_hartree"]
            for s in target["states"]
            if s["spin"] == spin and s["irrep"] == irrep
        ]
        np.testing.assert_allclose(expected[sector][:5], original, atol=1e-7, rtol=0)
    return expected


def stage(
    work: Path, report: dict, name: str, command: list[str], deck: Path, timeout: float, env: dict
) -> None:
    started = time.time()
    with deck.open() as inp, (work / f"{name}.out").open("x") as out:
        p = subprocess.Popen(
            command,
            cwd=work,
            env=env,
            stdin=inp,
            stdout=out,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        timed_out = False
        try:
            p.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(p.pid, signal.SIGTERM)
            try:
                p.wait(timeout=20)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
    report["stages"].append(
        {
            "name": name,
            "original_exit": p.returncode,
            "timed_out": timed_out,
            "wall_seconds": time.time() - started,
        }
    )
    (work / "pilot.json").write_text(json.dumps(report, indent=2) + "\n")
    text = (work / f"{name}.out").read_text(errors="replace")
    if (
        p.returncode
        or timed_out
        or any(
            s in text
            for s in ("MPI_ABORT", "***DUE TO AN ERROR", "CO_PTR_FAIL", "Fortran runtime error")
        )
    ):
        raise ValueError(f"Native {name} failed")


def compute(root: Path, work: Path, contract: dict) -> None:
    expected = expected_states(root, contract)
    reference = root / contract["target_reference_directory"]
    geometry = reference / "output/CO/geom1"
    target = json.loads((reference / "target.json").read_text())
    report = {
        "scope": contract["scientific_question"],
        "stages": [],
        "target_sectors": {},
        "channel_import_verdict": None,
        "scattering_channel_verdict": None,
        "scattering_roots": None,
    }
    names = ["fort.16", "co.molden", "co.casscf.chk"] + [
        f"fort.{unit}" for unit in (710, 730, 711, 731, 712, 732, 713, 733)
    ]
    report["reused_input_sha256"] = {}
    for name in names:
        shutil.copy2(geometry / name, work / name)
        report["reused_input_sha256"][name] = sha(work / name)
    for index, (spin, irrep) in enumerate(SECTORS, 1):
        count = len(expected[f"{spin}.{irrep}"])
        name = f"target.scatci.{spin}.{irrep}"
        text = (geometry / "inputs" / f"{name}.inp").read_text()
        text = replace_assignment(text, "nstat", str(count))
        text = replace_assignment(text, "maxiter", "1000")
        text = replace_assignment(text, "nciset", str(index))
        deck = work / f"{name}.inp"
        deck.write_text(text)
        telemetry = work / f"{name}.telemetry.jsonl"
        env = os.environ | {
            "LD_PRELOAD": str(root / contract["telemetry_directory"] / "scattering_telemetry.so"),
            "UKRMOL_EPS_AUDIT": str(telemetry),
        }
        stage(
            work,
            report,
            name,
            [
                "mpiexec",
                "-n",
                str(contract["resources"]["mpi_ranks"]),
                "/opt/ukrmolp/bin.double/mpi-scatci",
                str(deck),
            ],
            deck,
            contract["resources"]["target_sector_timeout_seconds"],
            env,
        )
        rows = [json.loads(line) for line in telemetry.read_text().splitlines()]
        if len(rows) != 1:
            raise ValueError("Missing unique native target residual telemetry")
        row = rows[0]
        residual = np.asarray([v["absolute_residual"] for v in row["eigenpairs"]])
        if (
            row["requested"] != count
            or row["converged"] < count
            or row["reason"] <= 0
            or row["solver_type"] != "krylovschur"
            or row["which"] != 1
            or len(residual) < count
            or not np.isfinite(residual).all()
            or max(residual) > 1e-9
        ):
            raise ValueError("Native target solver/action gate failed")
        report["target_sectors"][f"{spin}.{irrep}"] = row
        (work / "pilot.json").write_text(json.dumps(report, indent=2) + "\n")
    denprop = replace_assignment(
        (geometry / "inputs/target.denprop.inp").read_text(), "ntgtl", "5,5,5,5,5,5,9,9"
    )
    deck = work / "target.denprop.inp"
    deck.write_text(denprop)
    stage(
        work,
        report,
        "target.denprop",
        ["/opt/ukrmolp/bin.double/denprop"],
        deck,
        contract["resources"]["denprop_timeout_seconds"],
        os.environ.copy(),
    )
    if "has completed execution" not in (work / "target.denprop.out").read_text():
        raise ValueError("DENPROP completion marker missing")
    report["property_inventory"] = property_inventory((work / "fort.24").read_text(), expected)
    difference = report["property_inventory"]["ground_dipole_z_au"] - target["ground_dipole_au"][2]
    if abs(difference) > 1e-5:
        raise ValueError("Native ground dipole gate failed")
    report.update(
        ground_dipole_difference_au=difference,
        channel_import_verdict=True,
        unchanged_orbital_hashes_passed=all(
            sha(work / n) == h for n, h in report["reused_input_sha256"].items()
        ),
    )
    if not report["unchanged_orbital_hashes_passed"]:
        raise ValueError("Frozen inputs changed during native preflight")
    (work / "pilot.json").write_text(json.dumps(report, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--profiled", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    root, work = args.root.resolve(), args.work.resolve()
    contract = json.loads(args.contract.read_text())
    if args.profiled:
        compute(root, work, contract)
    else:
        import sys

        work.mkdir()
        (work / "scratch").mkdir()
        env = os.environ | {"TMPDIR": str(work / "scratch"), "PSI_SCRATCH": str(work / "scratch")}
        raise SystemExit(
            run_profiled(
                work,
                env,
                [
                    "python3",
                    "-m",
                    "projects.ukrmol_co.a2_channel_preflight",
                    *sys.argv[1:],
                    "--profiled",
                ],
            )
        )


if __name__ == "__main__":
    main()
