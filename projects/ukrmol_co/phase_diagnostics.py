"""Reconstruct native residuals and check independent, held-point phase fits.

The selected pinned engine converts requested eV to Rydberg using 0.0735.
The independent fitter consumes those actual Hartree energies. Modern physical
eV conversions use qscat.units only when reporting the result.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from qscat.units import HARTREE_TO_EV

from projects.ukrmol_co.phase_fit import fit_phase, predict_phase

# Pinned rsolve.f/eigenp.f RYD literal, verified against supplied source below.
NATIVE_RYDBERG_PER_INPUT_EV = 0.0735


def numbers(text: str) -> list[float]:
    """Read the native decimal/exponent fields, preserving their printed precision."""
    return [float(x.replace("D", "E")) for x in re.findall(r"[-+]?\d+\.\d+(?:D[-+]?\d+)?", text)]


def diagnose(run: Path, replays: Path, engine_source: Path) -> dict:
    """Compare every completed window replay with raw points and held-point fits."""
    for name in ("rsolve.f", "eigenp.f"):
        if not re.search(r"RYD/0\.073500D0/", (engine_source / name).read_text()):
            raise ValueError("Pinned native energy conversion differs; resolve the current source")
    native = json.loads((replays / "execution.json").read_text())
    target = json.loads((run / "result.json").read_text())["target_excitations_hartree"]
    first_excited_threshold = min(value for value in target.values() if value > 0)
    rows = []
    for entry in native["replays"]:
        directory = replays / (
            f"{entry['irrep']}-nback{entry['nback']}-window{entry['window_label']}"
        )
        raw = (directory / "reson.out").read_text()
        grid = np.array(
            [
                [float(e), float(eta.replace("D", "E"))]
                for e, eta in re.findall(
                    r"Grid point\s+\d+\s+E=\s*([-\d.]+)\s+Ryd, Eta=\s*([-\d.D+]+)", raw
                )
            ]
        )
        table = np.loadtxt(run / f"output/CO/geom1/eigenph.doublet.{entry['irrep']}")
        energy_rydberg = table[:, 0] * NATIVE_RYDBERG_PER_INPUT_EV
        indices = np.argmin(np.abs(grid[:, 0, None] - energy_rydberg), axis=1)
        np.testing.assert_allclose(energy_rydberg[indices], grid[:, 0], atol=5.01e-6, rtol=0)
        np.testing.assert_allclose(table[indices, 1], grid[:, 1], atol=5.1e-8, rtol=0)
        if len(set(indices)) != len(indices):
            raise ValueError("Native printed points do not map one-to-one to the saved grid")
        energy = energy_rydberg[indices] / 2
        if energy[-1] >= first_excited_threshold:
            raise ValueError("Phase fits must remain below the first excited threshold")
        phase = grid[:, 1]
        fitted = raw.split("Fitted resonance parameters")[-1]
        positions = numbers(re.search(r"Positions / Ryd([^\n]+)", fitted)[1])
        widths = numbers(re.search(r"Widths\s+/ Ryd([^\n]+)", fitted)[1])
        if len(positions) != 1 or len(widths) != 1:
            raise ValueError("The phase diagnostic requires exactly one native candidate")
        position, width = positions[0] / 2, widths[0] / 2
        coefficients = numbers(re.search(r"Background([^\n]+)", fitted)[1])
        residuals = np.array(numbers(fitted.split("Residues")[-1].split("Goodness factor")[0]))
        goodness = numbers(re.search(r"Goodness factor:([^\n]+)", fitted)[1])[0]
        if residuals.size != len(energy) or len(coefficients) != entry["nback"]:
            raise ValueError("Incomplete native residue/background record")
        rounding_bound = 2 * np.sum(np.abs(residuals)) * 5e-8 + len(residuals) * (5e-8) ** 2
        np.testing.assert_allclose(
            residuals @ residuals, goodness, atol=rounding_bound + 1e-10, rtol=0
        )
        prediction = np.arctan2(width / 2, position - energy) + (
            np.polynomial.polynomial.polyval(2 * energy, coefficients)
        )
        np.testing.assert_allclose(phase - prediction, residuals, atol=1e-6, rtol=0)
        span = energy[-1] - energy[0]
        starts = []
        for displacement, width_factor in ((0, 1), (-0.05, 0.5), (0.05, 2)):
            trial = fit_phase(
                energy, phase, entry["nback"], position + displacement * span, width * width_factor
            )
            if not trial["optimizer_success"] or trial["at_parameter_bound"]:
                raise ValueError("Independent phase fit failed or reached its bounds")
            starts.append(trial)
        best = min(starts, key=lambda x: x["residual_sum_squares_rad2"])
        folds = []
        for fold in range(5):
            held = np.arange(len(energy)) % 5 == fold
            trial = fit_phase(
                energy[~held],
                phase[~held],
                entry["nback"],
                best["position_hartree"],
                best["full_width_hartree"],
            )
            if not trial["optimizer_success"] or trial["at_parameter_bound"]:
                raise ValueError("Held-point phase fit failed or reached its bounds")
            difference = predict_phase(energy[held], trial) - phase[held]
            folds.append(
                {
                    "fold": fold,
                    "held_indices": np.flatnonzero(held).tolist(),
                    "fit": trial,
                    "held_residual_rad": difference.tolist(),
                }
            )
        held_residuals = np.concatenate([f["held_residual_rad"] for f in folds])
        rows.append(
            {
                "irrep": entry["irrep"],
                "nback": entry["nback"],
                "window": entry["window_label"],
                "native_position_hartree": position,
                "native_full_width_hartree": width,
                "native_goodness_sum_squares_rad2": goodness,
                "native_printed_residual_rms_rad": float(np.sqrt(np.mean(residuals**2))),
                "energy_hartree": energy.tolist(),
                "phase_rad": phase.tolist(),
                "independent_fit": best,
                "start_controls": starts,
                "held_point_folds": folds,
                "maximum_start_position_difference_ev": max(
                    abs(s["position_hartree"] - best["position_hartree"]) * HARTREE_TO_EV
                    for s in starts
                ),
                "maximum_start_width_fractional_difference": max(
                    abs(s["full_width_hartree"] - best["full_width_hartree"])
                    / best["full_width_hartree"]
                    for s in starts
                ),
                "native_position_difference_ev": (best["position_hartree"] - position)
                * HARTREE_TO_EV,
                "native_width_difference_ev": (best["full_width_hartree"] - width) * HARTREE_TO_EV,
                "held_point_rms_rad": float(np.sqrt(np.mean(held_residuals**2))),
                "held_point_max_abs_rad": float(np.max(np.abs(held_residuals))),
                "native_output_sha256": hashlib.sha256(
                    (directory / "reson.out").read_bytes()
                ).hexdigest(),
            }
        )
    return {
        "run": run.name,
        "native_replay_diagnostic": replays.name,
        "rows": rows,
        "scope": "Native residue audit and independent phase fits; pole assignment remains open",
        "native_rydberg_per_input_ev": NATIVE_RYDBERG_PER_INPUT_EV,
        "physical_ev_per_requested_ev": NATIVE_RYDBERG_PER_INPUT_EV * HARTREE_TO_EV / 2,
        "energy_grid_precision": "Native conversion; printed grid agreement within 5.01e-6 Ryd",
        "input_sha256": {
            str(Path(directory.name) / name): hashlib.sha256(
                (directory / name).read_bytes()
            ).hexdigest()
            for directory, name in [
                (run, "config.json"),
                (run, "result.json"),
                (replays, "execution.json"),
                *[(run, f"output/CO/geom1/eigenph.doublet.{irrep}") for irrep in ("B1", "B2")],
            ]
        },
        "source_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                Path(__file__),
                Path(__file__).with_name("phase_fit.py"),
                *[engine_source / name for name in ("rsolve.f", "eigenp.f")],
            ]
        },
    }


def main() -> None:
    """Preserve the independent diagnostic in a fresh output directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("replays", type=Path)
    parser.add_argument("--engine-source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    result = diagnose(args.run, args.replays, args.engine_source)
    (args.output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "fits": len(result["rows"]),
                "maximum_native_position_difference_ev": max(
                    abs(r["native_position_difference_ev"]) for r in result["rows"]
                ),
                "maximum_native_width_difference_ev": max(
                    abs(r["native_width_difference_ev"]) for r in result["rows"]
                ),
                "physical_ev_per_requested_ev": result["physical_ev_per_requested_ev"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
