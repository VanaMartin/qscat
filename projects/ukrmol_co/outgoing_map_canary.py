"""One-anchor, bounded real-S / complex-DtN feasibility application."""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.native_full_reconstruction import load_model
from projects.ukrmol_co.outgoing_map import outgoing_map, physical_s
from projects.ukrmol_co.saved_matrices import scattering_matrix


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def run(root, output, contract):
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    report = dict(
        case=contract["case"],
        rows=[],
        canary_verdict=None,
        scientific_exit=2,
        native_calls=0,
        pole_searches=0,
        qualified_pole_verdict=None,
        full_sweep_released=False,
        route_holdout_s_matrices_analyzed=0,
    )
    raw = root / contract["raw_directory"]

    def save():
        report["wall_seconds"] = time.monotonic() - started
        usage = resource.getrusage(resource.RUSAGE_SELF)
        report["peak_process_rss_platform_units"] = usage.ru_maxrss
        report["process_cpu_seconds"] = usage.ru_utime + usage.ru_stime
        for name in ["memory.peak", "cpu.stat"]:
            path = Path("/sys/fs/cgroup") / name
            if path.is_file():
                report["cgroup_" + name] = path.read_text()
        temporary = output / "result.partial"
        temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        temporary.replace(output / "result.json")

    def verify():
        for name, want in contract["input_sha256"].items():
            if digest(root / name) != want:
                raise ValueError(f"Changed prerequisite: {name}")
        manifest = json.loads((raw / "input-manifest.json").read_text())
        if len(manifest) != 123:
            raise ValueError("Changed original input inventory")
        for name, want in manifest.items():
            if digest(raw / name) != want:
                raise ValueError(f"Changed original input: {name}")
        return len(manifest)

    try:
        report["inputs_verified_before"] = verify()
        with (output / "controls.log").open("x") as log:
            result = subprocess.run(
                [sys.executable, "-m", "projects.ukrmol_co.test_outgoing_map"],
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=120,
            )
        report["controls_original_exit"] = result.returncode
        if result.returncode:
            raise ValueError("Toy controls failed; no CO application")
        previous = json.loads((root / contract["qualified_report"]).read_text())
        if previous["saved_matrix_reproduction_verdict"] is not False:
            raise ValueError("Original aggregate failure must remain false")
        prior = next(row for row in previous["cases"] if row["case"] == contract["case"])
        if prior["saved_matrix_reproduction_verdict"] is not True:
            raise ValueError("Unqualified prerequisite case")
        if contract["case"] != "compressed/nine/B1" or contract["development_indices"] != [
            60,
            70,
            80,
        ]:
            raise ValueError("Scope differs from the one-anchor declaration")
        model = load_model(raw / contract["case"])
        common = (model["thresholds"], model["angular"], model["c2"], model["c3"])
        report["original_aggregate_verdict"] = False
        report["boundary_channels"] = len(model["thresholds"])
        save()
        for index in contract["development_indices"]:
            if index % 5 == 4:
                raise ValueError("Reserved S index is forbidden")
            energy = model["energies"][index]
            q = surface_matrix(model["poles"], model["amplitudes"], model["ground"] + energy / 2)
            reference = scattering_matrix(model["saved"][index])
            for endpoint in contract["real_endpoints_bohr"]:
                row = dict(
                    kind="real",
                    index=index,
                    energy_rydberg=float(energy),
                    endpoint=endpoint,
                    completed=False,
                )
                report["rows"].append(row)
                save()
                y, a, diagnostics = outgoing_map(
                    *common,
                    energy,
                    18,
                    endpoint,
                    contract["real_angle_deg"],
                    amplitudes=True,
                    timeout=contract["resources"]["propagation_cap_seconds"],
                    max_nfev=contract["resources"]["rhs_cap"],
                )
                actual, invariants = physical_s(q, y, a)
                np.savez(
                    output / f"real-{index}-{endpoint}.npz", y=y, a=a, s=actual, reference=reference
                )
                error = float(np.max(abs(actual - reference)))
                checks = dict(
                    s=error <= contract["s_absolute_budget"],
                    flux=invariants["flux_defect"] <= contract["flux_absolute_budget"],
                    wronskian=invariants["wronskian_residual"]
                    <= contract["invariant_absolute_budget"],
                    unitarity=invariants["unitarity_defect"]
                    <= contract["invariant_absolute_budget"],
                    reciprocity=invariants["reciprocity_defect"]
                    <= contract["invariant_absolute_budget"],
                )
                row.update(
                    completed=True,
                    diagnostics=diagnostics,
                    invariants=invariants,
                    s_absolute_error=error,
                    checks=checks,
                )
                save()
                if not all(checks.values()):
                    report.update(
                        canary_verdict=False, scientific_exit=1, stopped_at="real accuracy gate"
                    )
                    return 1
        maps = []
        for ordinal, (endpoint, angle, order, rtol) in enumerate(contract["complex_settings"]):
            row = dict(
                kind="complex",
                ordinal=ordinal,
                endpoint=endpoint,
                angle=angle,
                order=order,
                rtol=rtol,
                completed=False,
            )
            report["rows"].append(row)
            save()
            energy = complex(*contract["complex_incident_rydberg"])
            y, _, diagnostics = outgoing_map(
                *common,
                energy,
                18,
                endpoint,
                angle,
                order,
                reference=energy.real,
                rtol=rtol,
                timeout=contract["resources"]["propagation_cap_seconds"],
                max_nfev=contract["resources"]["rhs_cap"],
            )
            np.savez(output / f"complex-{ordinal}.npz", y=y)
            row.update(completed=True, diagnostics=diagnostics)
            maps.append(y)
            save()
        for row, y in zip(report["rows"][-len(maps) :], maps, strict=True):
            error = float(18 * np.max(abs(y - maps[-1])))
            symmetry = float(18 * np.max(abs(y - y.T)))
            row.update(
                radius_scaled_map_difference=error,
                radius_scaled_transpose_defect=symmetry,
                checks=dict(
                    map=error <= contract["radius_scaled_y_budget"],
                    transpose=symmetry <= contract["radius_scaled_y_budget"],
                ),
            )
        verdict = all(all(row["checks"].values()) for row in report["rows"])
        report.update(canary_verdict=verdict, scientific_exit=0 if verdict else 1)
        return report["scientific_exit"]
    except BaseException as error:
        report.update(error=repr(error), scientific_exit=2)
        if report["rows"] and not report["rows"][-1]["completed"]:
            report["rows"][-1]["error"] = repr(error)
        return 2
    finally:
        try:
            report["inputs_verified_after"] = verify()
        except BaseException as error:
            report.update(integrity_error=repr(error), scientific_exit=2)
        save()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = run(args.root, args.output, json.loads(args.contract.read_text()))
    # The persisted integrity result is authoritative even if finally changed it.
    raise SystemExit(json.loads((args.output / "result.json").read_text())["scientific_exit"])
