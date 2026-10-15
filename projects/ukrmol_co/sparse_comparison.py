"""Independent raw reading and fixed-window sensitivity of two selected-root runs.

A stable root-count pair is evidence of sensitivity, not an absolute bound on the
uncorrected omitted pole sum or certification of a resonance pole.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

from projects.ukrmol_co.phase_fit import fit_phase
from projects.ukrmol_co.sparse_pilot import audit_sector
from projects.ukrmol_co.sparse_scattering import ci_spectrum, interface_digest


def compare_phases(reference: np.ndarray, candidate: np.ndarray, contract: dict) -> dict:
    """Compare complete matched grids and identical-window fits at each background order."""
    checks = contract["raw_checks"]
    for phase in (reference, candidate):
        if phase.shape != (99, 2) or not np.isfinite(phase).all():
            raise ValueError("Incomplete or nonfinite comparison grid")
        np.testing.assert_allclose(
            phase[:, 0],
            0.1 + 0.05 * np.arange(99),
            atol=checks["phase_grid_atol_requested_ev"],
            rtol=0,
        )
    np.testing.assert_allclose(reference[:, 0], candidate[:, 0], atol=5e-5, rtol=0)
    delta = float(np.max(abs((candidate[:, 1] - reference[:, 1] + np.pi / 2) % np.pi - np.pi / 2)))
    budgets = contract["observable_budgets"]
    low, high = budgets["fixed_requested_energy_window_ev"]
    window = (reference[:, 0] >= low - 5e-5) & (reference[:, 0] <= high + 5e-5)
    energy = reference[window, 0] * contract["native_rydberg_per_requested_ev"] / 2
    rows = []
    for nback in budgets["fit_background_orders"]:
        row = {"nback": nback, "reference_fit": None, "candidate_fit": None, "passed": None}
        try:
            fits = [
                fit_phase(
                    energy,
                    phase[window, 1],
                    nback,
                    contract["fit_initial_position_hartree"],
                    contract["fit_initial_full_width_hartree"],
                )
                for phase in (reference, candidate)
            ]
            row.update(reference_fit=fits[0], candidate_fit=fits[1])
            if all(f["optimizer_success"] and not f["at_parameter_bound"] for f in fits):
                conversion = contract["hartree_to_ev"]
                position = (
                    abs(fits[1]["position_hartree"] - fits[0]["position_hartree"]) * conversion
                )
                width = (
                    abs(fits[1]["full_width_hartree"] - fits[0]["full_width_hartree"]) * conversion
                )
                width_budget = max(
                    budgets["full_width_floor_ev"],
                    budgets["full_width_relative"] * fits[0]["full_width_hartree"] * conversion,
                )
                row.update(
                    position_difference_ev=position,
                    full_width_difference_ev=width,
                    full_width_budget_ev=width_budget,
                    passed=bool(position <= budgets["position_ev"] and width <= width_budget),
                )
        except (ValueError, np.linalg.LinAlgError) as error:
            row["error"] = str(error)
        rows.append(row)
    valid = all(row["passed"] is not None for row in rows) and bool(rows)
    return {
        "maximum_phase_difference_rad_mod_pi": delta,
        "phase_budget_passed": delta <= budgets["phase_rad_mod_pi"],
        "fixed_window_fits": rows,
        "root_pair_stability_verdict": (
            bool(delta <= budgets["phase_rad_mod_pi"] and all(row["passed"] for row in rows))
            if valid
            else None
        ),
    }


def read_run(path: Path, roots: int, contract: dict) -> tuple[dict, dict]:
    """Verify original exits, solver coverage, raw exports and symmetry before fitting."""
    saved = json.loads((path / "pilot.json").read_text())
    resources = json.loads((path / "resources.json").read_text())
    if not saved.get("all_stages_complete") or resources["exit_code"] != 0:
        raise ValueError("Native run is incomplete or failed")
    phases, result = {}, {"resources": resources, "sectors": {}}
    for sector in contract["sectors"]:
        directory = path / sector
        stages = saved["sectors"][sector]["stages"]
        if [(s["name"], s["original_exit"]) for s in stages] != [
            (stage, 0) for stage in ("scatci", "rsolve", "eigenp", "reson")
        ] or any(s["timed_out"] for s in stages):
            raise ValueError("Native stage exit or timeout failed")
        raw = audit_sector(
            directory,
            roots,
            contract["dimension"],
            contract["continuum_coefficients"],
            contract.get("boundary_channels", 250),
        )
        telemetry = json.loads((directory / "eps-telemetry.jsonl").read_text())
        if telemetry["which"] != contract["raw_checks"]["native_selection_which"]:
            raise ValueError("Native selection changed")
        raw["solver_selection"] = {k: v for k, v in telemetry.items() if k != "eigenpairs"}
        raw["original_stages"] = stages
        unit = int(re.search(r"luphso\s*=\s*(\d+)", (directory / "eigenp.inp").read_text())[1])
        phases[sector] = np.loadtxt(directory / f"fort.{unit}", ndmin=2)
        result["sectors"][sector] = raw
    pi = float(
        np.max(abs((phases["B1"][:, 1] - phases["B2"][:, 1] + np.pi / 2) % np.pi - np.pi / 2))
    )
    if pi > contract["raw_checks"]["pi_phase_rad_mod_pi"]:
        raise ValueError("Independent Pi check failed")
    result["pi_phase_difference_rad_mod_pi"] = pi
    result["target_map"] = saved["target_map"]
    return result, phases


def compare_runs(reference: Path, candidate: Path, contract: dict) -> dict:
    """Reconstruct the two-sector comparison without trusting saved producer verdicts."""
    left, left_phases = read_run(reference, contract["baseline_roots"], contract)
    right, right_phases = read_run(candidate, contract["candidate_roots"], contract)
    if left["target_map"] != right["target_map"]:
        raise ValueError("Native target mapping changed")
    result = {
        "reference": str(reference),
        "candidate": str(candidate),
        "reference_raw_reading": left,
        "candidate_raw_reading": right,
        "sectors": {},
        "spectral_truncation_verdict": None,
        "electronic_model_verdict": None,
        "scope": "Root-pair sensitivity only; omitted-spectrum and pole uncertainty unqualified",
    }
    for sector in contract["sectors"]:
        a, b = reference / sector, candidate / sector
        # Only the requested root count may differ in the native Hamiltonian deck.
        decks = [
            re.sub(r"\bnstat\s*=\s*\d+", "nstat = ROOTS", (p / "scatci.inp").read_text())
            for p in (a, b)
        ]
        if decks[0] != decks[1]:
            raise ValueError("Native Hamiltonian input changed beyond root count")
        for stage in ("rsolve", "eigenp", "reson"):
            if (a / f"{stage}.inp").read_bytes() != (b / f"{stage}.inp").read_bytes():
                raise ValueError("Outer input changed")
        energies_a, _ = ci_spectrum(a / "fort.25")
        energies_b, _ = ci_spectrum(b / "fort.25")
        np.testing.assert_allclose(
            energies_a,
            energies_b[: len(energies_a)],
            atol=contract["raw_checks"]["shared_prefix_energy_atol_hartree"],
            rtol=0,
        )
        digests = [interface_digest(p / "fort.21", boundary=True) for p in (a, b)]
        if digests[0] != digests[1]:
            raise ValueError("Native channel/threshold/multipole interface changed")
        comparison = compare_phases(left_phases[sector], right_phases[sector], contract)
        comparison.update(
            shared_prefix_energy_error_hartree=float(
                np.max(abs(energies_a - energies_b[: len(energies_a)]))
            ),
            interface_sha256=digests[0],
        )
        raw = [reading["sectors"][sector] for reading in (left, right)]
        if any(r["native_maxfit_diagnostic"] or not r["native_resonances"] for r in raw):
            comparison["root_pair_stability_verdict"] = None
            comparison["unset_reason"] = "Empty or MAXFIT native fit"
        result["sectors"][sector] = comparison
    verdicts = [r["root_pair_stability_verdict"] for r in result["sectors"].values()]
    result["root_pair_stability_verdict"] = None if None in verdicts else all(verdicts)
    result["independent_raw_reading_passed"] = True
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("reference", "candidate", "contract", "output"):
        parser.add_argument(f"--{field}", type=Path, required=True)
    args = parser.parse_args()
    result = compare_runs(args.reference, args.candidate, json.loads(args.contract.read_text()))
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
