"""Four real-axis development checks on qualified own64/400 exports; no pole search."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import signal
import sys
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.compressed_a2_contrast_run import real_guard, save_json
from projects.ukrmol_co.native_full_reconstruction import _model_from_interface
from projects.ukrmol_co.saved_matrices import read_interface
from projects.ukrmol_co.stretched_basis_properties import digest


def authenticate(root, source, contract):
    for base, key in [(root, "input_sha256"), (source, "source_sha256")]:
        for name, want in contract[key].items():
            path = base / name
            if not path.resolve().is_relative_to(base.resolve()) or digest(path) != want:
                raise ValueError("Changed/out-of-root real guard prerequisite: " + name)
    if list(contract["models"]) != ["dz", "tz"]:
        raise ValueError("Only the frozen DZ then TZ own64/400 exports are supported")
    receipts = {}
    for case, spec in contract["models"].items():
        receipt = json.loads((root / spec["receipt"]).read_text())
        if not (
            receipt["case"] == case
            and receipt["inner_export_verdict"] is True
            and receipt["public_reader_exit"] == 0
            and receipt["original_exits"] == [0] * 6
            and receipt["contract_sha256"] == spec["inner_contract_sha256"]
            and (receipt["roots"], receipt["target_count"], receipt["boundary_channels"])
            == (2048, 64, 400)
            and receipt["spectral_truncation_verdict"] is None
            and not receipt["full_sweep_released"]
        ):
            raise ValueError("Independent public own export required: " + case)
        receipts[case] = receipt
    return receipts


def own_model(directory, receipt):
    interface = read_interface(directory / "fort.10", directory / "fort.21")
    model = _model_from_interface(directory, interface, (64, 400))
    for name, want in receipt["array_sha256"].items():
        value = np.asarray(model[name], dtype="<f8")
        if not np.isfinite(value).all() or hashlib.sha256(value.tobytes()).hexdigest() != want:
            raise ValueError("Own array differs from independent public reading: " + name)
    return model


def compute(root, source, output, contract):
    started = time.monotonic()
    report = dict(
        scientific_exit=2,
        real_guard_verdict=None,
        real_maps_reserved=0,
        models={},
        native_calls=0,
        new_inner_solves=0,
        new_complex_maps=0,
        pole_verdict=None,
        model_insensitivity_verdict=None,
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
        automatic_successor=False,
    )

    def save():
        report["wall_seconds"] = time.monotonic() - started
        save_json(output / "result.json", report)

    try:
        receipts = authenticate(root, source, contract)
        save()
        for case, spec in contract["models"].items():
            model = own_model(root / spec["native_directory"], receipts[case])
            directory = output / case
            directory.mkdir()
            path = real_guard(model, directory, "stretched/cas12/" + case + "/B1", contract, report)
            report["models"][case] = dict(guard_sha256=digest(path), completed=True)
            save()
        authenticate(root, source, contract)
        report.update(scientific_exit=0, real_guard_verdict=True, inputs_unchanged=True)
    except Exception as error:
        report["error"] = repr(error)
    finally:
        save()
    return report["scientific_exit"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "source", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--profiled", action="store_true")
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text())
    if args.profiled:

        def timeout(*_):
            raise TimeoutError("Frozen four-map guard wall cap exhausted")

        signal.signal(signal.SIGALRM, timeout)
        signal.alarm(contract["guard_timeout_seconds"])
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
            "projects.ukrmol_co.stretched_real_guard",
            *sys.argv[1:],
            "--profiled",
        ],
    )


if __name__ == "__main__":
    raise SystemExit(main())
