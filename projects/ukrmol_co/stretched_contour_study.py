"""Authorized fixed three-endpoint DZ80 diagnostic; never retries the original guard."""

from __future__ import annotations

import argparse
import json
import os
import signal
import sys
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.compressed_a2_contrast_run import save_json
from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.outgoing_map import outgoing_map, physical_s
from projects.ukrmol_co.paired_pole_character import real_checks
from projects.ukrmol_co.saved_matrices import scattering_matrix
from projects.ukrmol_co.stretched_basis_properties import digest
from projects.ukrmol_co.stretched_real_guard import own_model


def authenticate(root, source, contract):
    if not contract["execution_authorized"] or (
        contract["development_index"],
        contract["boundary_bohr"],
        contract["angle_degrees"],
        contract["endpoints_bohr"],
        contract["tail_order"],
        contract["rtol"],
    ) != (80, 18, 20, [60, 90, 120], 10, 1e-12):
        raise ValueError("Outside the explicitly authorized fixed three-map study")
    for base, key in [(root, "input_sha256"), (source, "source_sha256")]:
        for name, want in contract[key].items():
            path = base / name
            if not path.resolve().is_relative_to(base.resolve()) or digest(path) != want:
                raise ValueError("Changed/out-of-root prerequisite: " + name)
    receipt = json.loads((root / contract["model"]["receipt"]).read_text())
    if not (
        receipt["case"] == "dz"
        and receipt["inner_export_verdict"] is True
        and receipt["public_reader_exit"] == 0
        and receipt["original_exits"] == [0] * 6
        and receipt["contract_sha256"] == contract["model"]["inner_contract_sha256"]
        and (receipt["roots"], receipt["target_count"], receipt["boundary_channels"])
        == (2048, 64, 400)
        and receipt["spectral_truncation_verdict"] is None
        and not receipt["full_sweep_released"]
    ):
        raise ValueError("Own-DZ public export receipt required")
    coupled = json.loads((root / contract["coupled_reading"]).read_text())
    if not (
        coupled["reader_exit"] == 0
        and coupled["result"]["coupled_control_verdict"] is True
        and coupled["result"]["integrations_verified"] == 54
        and coupled["result"]["original_exits"] == [0] * 6
        and coupled["result"]["co_guard_repaired"] is False
    ):
        raise ValueError("Independent public coupled-control pass required")
    original = json.loads((root / contract["original_guard_result"]).read_text())
    if original["scientific_exit"] != 2 or original["real_guard_verdict"] is not None:
        raise ValueError("Original failed guard must remain unchanged")
    return receipt


def compute(root, source, output, contract):
    started = time.monotonic()
    result = dict(
        scientific_exit=2,
        contour_study_verdict=None,
        real_rows=[],
        real_maps_reserved=0,
        native_calls=0,
        new_inner_solves=0,
        new_complex_maps=0,
        tz_started=False,
        original_guard_requalified=False,
        route_holdout_s_matrices_analyzed=0,
        automatic_successor=False,
        full_sweep_released=False,
    )

    def save():
        result["science_wall_seconds"] = time.monotonic() - started
        for name in ["memory.peak", "cpu.stat"]:
            path = Path("/sys/fs/cgroup") / name
            if path.exists():
                result["cgroup_" + name] = path.read_text()
        save_json(output / "result.json", result)

    try:
        receipt = authenticate(root, source, contract)
        model = own_model(root / contract["model"]["native_directory"], receipt)
        index = contract["development_index"]
        energy = model["energies"][index]
        q = surface_matrix(model["poles"], model["amplitudes"], model["ground"] + energy / 2)
        reference = scattering_matrix(model["saved"][index])
        captures = []
        for endpoint in contract["endpoints_bohr"]:
            if time.monotonic() - started > contract["science_timeout_seconds"]:
                raise TimeoutError("Authorized combined science cap")
            row = dict(index=index, endpoint=endpoint, completed=False)
            result["real_rows"].append(row)
            result["real_maps_reserved"] += 1
            save()
            y, a, diagnostics = outgoing_map(
                model["thresholds"],
                model["angular"],
                model["c2"],
                model["c3"],
                energy,
                contract["boundary_bohr"],
                endpoint,
                contract["angle_degrees"],
                contract["tail_order"],
                rtol=contract["rtol"],
                amplitudes=True,
                timeout=contract["per_map_seconds"],
                max_nfev=contract["rhs_cap"],
            )
            actual, invariants = physical_s(q, y, a)
            np.savez(
                output / f"endpoint-{endpoint}.npz",
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
                **real_checks(actual, reference, invariants, contract),
            )
            captures.append((y, a))
            save()
        stability = dict(
            y_absolute=float(abs(captures[-1][0] - captures[-2][0]).max()),
            a_absolute=float(abs(captures[-1][1] - captures[-2][1]).max()),
        )
        result["last_two_stability"] = stability
        if not all(result["real_rows"][-1]["checks"].values()) or any(
            value > contract["last_two_y_and_a_absolute_stability_budget"]
            for value in stability.values()
        ):
            raise ValueError("Fixed study longest real gates or last-two map stability failed")
        authenticate(root, source, contract)
        result.update(scientific_exit=0, contour_study_verdict=True, inputs_unchanged=True)
    except Exception as error:
        result["error"] = repr(error)
    finally:
        save()
    return result["scientific_exit"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "source", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--profiled", action="store_true")
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text())
    if args.profiled:

        def timeout(*_):
            raise TimeoutError("Authorized study science wall cap")

        signal.signal(signal.SIGALRM, timeout)
        signal.alarm(contract["science_timeout_seconds"])
        return compute(args.root, args.source, args.output, contract)
    from projects.ukrmol_co.run import run_profiled

    args.output.mkdir(exist_ok=False)
    (args.output / "scratch").mkdir()
    return run_profiled(
        args.output,
        os.environ | dict(TMPDIR=str(args.output / "scratch")),
        [
            sys.executable,
            "-m",
            "projects.ukrmol_co.stretched_contour_study",
            *sys.argv[1:],
            "--profiled",
        ],
    )


if __name__ == "__main__":
    raise SystemExit(main())
