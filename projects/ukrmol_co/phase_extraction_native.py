"""Apply the calibrated extraction rule to qualified saved nine-A2 anchor outputs."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import time
from pathlib import Path

from projects.ukrmol_co.a2_scattering import phases
from projects.ukrmol_co.a2_scattering_reading import reconstruct
from projects.ukrmol_co.phase_extraction import estimate_phase
from projects.ukrmol_co.saved_model_controls import study


def read_anchors(root: Path, calibration: dict, contract: dict) -> dict:
    """Gate on unchanged analytic calibration and original raw evidence before fitting."""
    if calibration.get("calibration_verdict") is not True:
        raise ValueError("Passing known-pole calibration is required before native application")
    if calibration["contract"] != contract:
        raise ValueError("Calibration and native extraction contracts differ")
    for name, digest in calibration["source_sha256"].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Calibrated source changed: {name}")
    start = time.monotonic()
    qualified = {name: reconstruct(root / name) for name in ("compressed", "equilibrium")}
    for reading in qualified.values():
        assert (
            reading["dense_baseline_verdict"] is reading["a2_channel_sensitivity_verdict"] is True
        )
    saved = root / "saved"
    old_contract = json.loads((saved / "contract.json").read_text())
    old = study(saved, old_contract)
    assert old["numerical_sensitivity_controls_passed"] is True
    for sector in contract["native_application"]["sectors"]:
        raw = old["raw_reading"]["stretchednine"]["sectors"][sector]
        assert raw["native_resonances"] and not raw["native_maxfit_diagnostic"]
    settings = contract["estimator"] | {"hartree_to_ev": contract["hartree_to_ev"]}
    native, verdicts = {}, []
    directories = {
        "compressed-nine": root / "compressed/nine",
        "equilibrium-nine": root / "equilibrium/nine",
        "stretched-nine": saved / old_contract["models"]["stretchednine"]["directory"],
    }
    for model in contract["native_application"]["models"]:
        windows = [
            [v * contract["native_rydberg_per_requested_ev"] / 2 for v in window]
            for window in contract["native_application"]["windows_requested_ev"][model]
        ]
        native[model] = {}
        for sector in contract["native_application"]["sectors"]:
            table = phases(directories[model] / sector)
            energy = table[:, 0] * contract["native_rydberg_per_requested_ev"] / 2
            row = estimate_phase(energy, table[:, 1], windows, settings)
            native[model][sector] = row
            verdicts.append(row["estimator_verdict"])
    return {
        "verified_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "source_sha256": {
            str(Path(__file__).relative_to(Path(__file__).resolve().parents[2])): hashlib.sha256(
                Path(__file__).read_bytes()
            ).hexdigest()
        }
        | calibration["source_sha256"],
        "contract": contract,
        "calibration_verdict": calibration["calibration_verdict"],
        "original_raw_file_manifests": {
            "compressed": qualified["compressed"]["raw_file_manifest"],
            "equilibrium": qualified["equilibrium"]["raw_file_manifest"],
            "saved": old["input_manifest"],
        },
        "original_controller_exits": {
            name: reading["original_controller_exit"] for name, reading in qualified.items()
        },
        "original_stretched_channel_sensitivity_verdict": False,
        "original_extraction_verdicts": {
            name: {s: row["extraction_stability_verdict"] for s, row in rows.items()}
            for name, rows in old["extraction_controls"].items()
        },
        "native_extraction": native,
        "native_extraction_verdict": None if None in verdicts else all(verdicts),
        "reader_wall_seconds": time.monotonic() - start,
        "native_solves": 0,
        "full_sweep_released": False,
        "scope": "New prospective estimator on saved nine-A2 phases after finite analytic "
        "calibration. Original all-order/window and channel rejections remain preserved. "
        "No native S-matrix pole, global identifiability or electronic-model certification.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "calibration", "contract", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = read_anchors(
        args.root, json.loads(args.calibration.read_text()), json.loads(args.contract.read_text())
    )
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if result["native_extraction_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
