"""Reconstruct a completed matched A2 control from checksum-fetched output files."""

from __future__ import annotations

import argparse
import datetime
import json
from pathlib import Path

import numpy as np

from projects.ukrmol_co.a2_scattering import channel_deck, channel_map, phases
from projects.ukrmol_co.analyze import native_resonances
from projects.ukrmol_co.sparse_comparison import compare_phases, read_run
from projects.ukrmol_co.sparse_comparison_owner import sha


def reconstruct(directory: Path) -> dict:
    """Recompute raw checks and both observable gates without a live campaign mount."""
    manifest = json.loads((directory / "manifest.json").read_text())
    for name, row in manifest.items():
        path = directory / name
        if path.stat().st_size != row["bytes"] or sha(path) != row["sha256"]:
            raise ValueError(f"Fetched output checksum failed: {name}")
    contract = json.loads((directory / "contract.json").read_text())
    owner = json.loads((directory / "execution.json").read_text())
    observer = json.loads((directory / "observer-execution.json").read_text())
    assert owner["original_controller_exit"] == observer["original_exit"] == 0
    assert observer["finish_reason"] == "Finite owner ended"
    assert [(s["label"], s["original_exit"]) for s in owner["steps"]] == [
        (n, 0)
        for n in (
            "runtime-gate",
            "check-inputs",
            "model-dense",
            "model-five",
            "baseline-gate",
            "model-nine",
            "compare",
        )
    ]
    inventory = json.loads((directory / "channel-input-gate.json").read_text())
    target = json.loads((directory / "independent-target-reading.json").read_text())
    assert target["channel_import_verdict"] is True
    assert inventory["target_map"] == target["target_map"]
    shared = json.loads((directory / "shared-input-audit.json").read_text())["actual_input_files"]
    raw, tables = {}, {}
    for model, count in (("five", 5), ("nine", 9)):
        settings = contract | contract["models"][model]
        raw[model], tables[model] = read_run(
            directory / model, contract["scattering_roots"], settings
        )
        assert raw[model]["target_map"] == channel_map(inventory, count)
        for sector in contract["sectors"]:
            path = directory / model / sector
            original = (directory / "original" / sector / "scatci.inp").read_text()
            assert (path / "scatci.inp").read_text() == channel_deck(
                original, inventory, count, contract["scattering_roots"]
            )
            for name in (
                "fort.16",
                "fort.24",
                "fort.26",
                "fort." + ("721" if sector == "B1" else "722"),
            ):
                assert shared[f"five/{sector}/{name}"] == shared[f"nine/{sector}/{name}"]
            assert sha(path / "fort.24") == target["input_manifest"]["fort.24"]["sha256"]
            assert (
                shared[f"{model}/{sector}/fort.26"]["sha256"]
                == target["input_manifest"]["fort.26"]["sha256"]
            )
            assert shared[f"{model}/{sector}/fort.24"]["sha256"] == sha(path / "fort.24")
            for name in ("rsolve", "eigenp", "reson"):
                assert (path / (name + ".inp")).read_bytes() == (
                    directory / "five" / sector / (name + ".inp")
                ).read_bytes()
    dense = directory / "dense"
    pilot = json.loads((dense / "pilot.json").read_text())
    resources = json.loads((dense / "resources.json").read_text())
    assert pilot["all_stages_complete"] and resources["exit_code"] == 0
    tables["dense"] = {}
    for sector in contract["sectors"]:
        stages = pilot["sectors"][sector]["stages"]
        assert [(s["name"], s["original_exit"]) for s in stages] == [
            (n, 0) for n in ("rsolve", "eigenp", "reson")
        ]
        assert not any(s["timed_out"] for s in stages)
        tables["dense"][sector] = phases(dense / sector)
        for name in ("rsolve", "eigenp", "reson"):
            assert (dense / sector / (name + ".inp")).read_bytes() == (
                directory / "five" / sector / (name + ".inp")
            ).read_bytes()
    pi = {}
    for model, values in tables.items():
        pi[model] = float(
            np.max(abs((values["B1"][:, 1] - values["B2"][:, 1] + np.pi / 2) % np.pi - np.pi / 2))
        )
        assert pi[model] <= contract["raw_checks"]["pi_phase_rad_mod_pi"]
    comparisons, fit_parity = {}, []
    for name, left, right, filename in (
        ("baseline", "dense", "five", "baseline-qualification.json"),
        ("channel", "five", "nine", "independent-comparison.json"),
    ):
        saved = json.loads((directory / filename).read_text())
        rows = {}
        for sector in contract["sectors"]:
            a, b = tables[left][sector], tables[right][sector]
            row = compare_phases(a, b, contract)
            verdict = row.pop("root_pair_stability_verdict")
            for model in (left, right):
                output = directory / model / sector / "reson.out"
                if "MAXFIT" in output.read_text().upper() or not native_resonances(output):
                    verdict = None
            row["observable_budget_verdict"] = verdict
            assert verdict is saved["sectors"][sector]["observable_budget_verdict"]
            delta = abs((b[:, 1] - a[:, 1] + np.pi / 2) % np.pi - np.pi / 2)
            row["maximum_phase_difference_requested_energy_ev"] = float(a[np.argmax(delta), 0])
            row["points_exceeding_phase_budget"] = int(
                np.count_nonzero(delta > contract["observable_budgets"]["phase_rad_mod_pi"])
            )
            np.testing.assert_allclose(
                row["maximum_phase_difference_rad_mod_pi"],
                saved["sectors"][sector]["maximum_phase_difference_rad_mod_pi"],
                atol=1e-14,
                rtol=0,
            )
            for fit, original in zip(
                row["fixed_window_fits"], saved["sectors"][sector]["fixed_window_fits"], strict=True
            ):
                assert fit["passed"] is original["passed"]
                if fit["passed"] is None:
                    continue
                for side in ("reference_fit", "candidate_fit"):
                    for key in ("position_hartree", "full_width_hartree"):
                        error = (
                            abs(fit[side][key] - original[side][key]) * contract["hartree_to_ev"]
                        )
                        assert error <= 1e-8
                        fit_parity.append(error)
                for field in ("position", "full_width"):
                    fit["signed_" + field + "_change_ev"] = (
                        fit["candidate_fit"][field + "_hartree"]
                        - fit["reference_fit"][field + "_hartree"]
                    ) * contract["hartree_to_ev"]
            rows[sector] = row
        verdicts = [r["observable_budget_verdict"] for r in rows.values()]
        comparisons[name] = {
            "sectors": rows,
            "verdict": None if None in verdicts else all(verdicts),
        }
    return {
        "verified_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "reader_source_sha256": sha(Path(__file__)),
        "checksum_verified_files": len(manifest),
        "raw_file_manifest": manifest,
        "original_controller_exit": owner["original_controller_exit"],
        "owner_finished_utc": owner["finished_utc"],
        "observer_original_exit": observer["original_exit"],
        "observer_capture_count": len(observer["captures"]),
        "dense_resources": resources,
        "selected_raw_reading": raw,
        "pi_phase_difference_rad_mod_pi": pi,
        "comparisons": comparisons,
        "maximum_fit_platform_parity_difference_ev": max(fit_parity, default=None),
        "dense_baseline_verdict": comparisons["baseline"]["verdict"],
        "a2_channel_sensitivity_verdict": comparisons["channel"]["verdict"],
        "independent_local_raw_reconstruction_passed": True,
        "spectral_truncation_verdict": None,
        "electronic_model_verdict": None,
        "global_continuity_verdict": None,
        "full_sweep_released": False,
        "scope": "Checksum-fetched raw export/grid/Pi/deck/fit reconstruction using shared public "
        "readers and fitter,not an independent scattering engine. Actual large shared inputs are "
        "checksum-audited remotely. Fixed-orbital five-versus-nine A2 sensitivity at2048only;"
        "full channel convergence,absolute omission,pole identity and extraction stability "
        "remain unqualified.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = reconstruct(args.directory)
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
