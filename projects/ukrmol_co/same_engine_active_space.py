"""Recheck saved active-space sensitivity with the audited CAS11 baseline."""

from __future__ import annotations

import argparse
import datetime
import json
from pathlib import Path

import numpy as np

from projects.ukrmol_co.a2_scattering_reading import reconstruct
from projects.ukrmol_co.saved_model_controls import read_model, study
from projects.ukrmol_co.sparse_comparison import compare_phases
from projects.ukrmol_co.sparse_comparison_owner import sha


def compare(cas11: Path, saved: Path, provenance: Path, contract: dict) -> dict:
    """Verify original outputs and compare every sector/order on the frozen window."""
    manifest = json.loads((provenance / "manifest.json").read_text())
    for name, row in manifest.items():
        path = provenance / name
        if sha(path) != row["sha256"] or path.stat().st_size != row["bytes"]:
            raise ValueError(f"Native-owner provenance checksum failed: {name}")
    owners = {}
    for name in ("cas11", "cas12", "cas12refined"):
        declaration = json.loads((provenance / name / "contract.json").read_text())
        owner = json.loads((provenance / name / "execution.json").read_text())
        assert owner["original_controller_exit"] == 0
        assert declaration["resources"]["mpi_ranks"] == 8
        assert owner["image_id"] == declaration["image_id"]
        assert owner["engine_verified"] == owner["engine_verified_after_native"]
        assert owner["engine_verified"]["installed_files"] == 242
        owners[name] = {
            "image_id": owner["image_id"],
            "engine": owner["engine_verified"],
            "engine_evidence_sha256": declaration["engine_evidence_sha256"],
            "mpi_ranks": declaration["resources"]["mpi_ranks"],
        }
    assert owners["cas11"] == owners["cas12"] == owners["cas12refined"]
    assert (provenance / "cas11/execution.json").read_bytes() == (
        cas11 / "execution.json"
    ).read_bytes()
    fresh = reconstruct(cas11)
    assert fresh["dense_baseline_verdict"] is True
    previous_contract = json.loads((saved / "contract.json").read_text())
    previous = study(saved, previous_contract)
    assert previous["numerical_sensitivity_controls_passed"] is True
    for key in ("observable_budgets", "hartree_to_ev", "native_rydberg_per_requested_ev"):
        assert contract[key] == previous_contract[key]
    fresh_contract = json.loads((cas11 / "contract.json").read_text())
    assert fresh_contract == json.loads((provenance / "cas11/contract.json").read_text())
    for key in ("observable_budgets", "hartree_to_ev", "native_rydberg_per_requested_ev"):
        assert contract[key] == fresh_contract[key]
    old11, old_phases = read_model(saved, previous_contract["models"]["eq11"], contract)
    groups, reproduction, deck_audit = {}, {}, {}
    for candidate in contract["cas12_models"]:
        model = previous_contract["models"][candidate]
        raw, tables = read_model(saved, model, contract)
        rows = {}
        for sector in contract["sectors"]:
            left = cas11 / "five" / sector
            right = saved / model["directory"] / sector
            for stage in ("rsolve", "eigenp", "reson"):
                a, b = left / (stage + ".inp"), right / (stage + ".inp")
                text = a.read_text()
                if stage == "rsolve":
                    for fresh_line, old_line in contract["rsolve_formatting_only_line_pairs"]:
                        assert "".join(fresh_line.split("!", 1)[0].split()) == "".join(
                            old_line.split("!", 1)[0].split()
                        )
                        assert text.splitlines().count(fresh_line) == 1
                        assert b.read_text().splitlines().count(old_line) == 1
                        text = text.replace(fresh_line + "\n", old_line + "\n")
                assert text.encode() == b.read_bytes()
                deck_audit[f"{candidate}/{sector}/{stage}"] = {
                    "fresh_sha256": sha(a),
                    "candidate_sha256": sha(b),
                    "strict_byte_equal": a.read_bytes() == b.read_bytes(),
                    "supplementary_exact_line_audit_passed": True,
                }
            table = np.loadtxt(left / ("fort.102" + sector[-1]), ndmin=2)
            row = compare_phases(table, tables[sector], contract)
            row["observable_insensitivity_verdict"] = row.pop("root_pair_stability_verdict")
            row["maximum_phase_requested_ev"] = float(
                table[
                    np.argmax(
                        abs((tables[sector][:, 1] - table[:, 1] + np.pi / 2) % np.pi - np.pi / 2)
                    ),
                    0,
                ]
            )
            for original in (fresh["selected_raw_reading"]["five"], raw):
                s = original["sectors"][sector]
                if not s["native_resonances"] or s["native_maxfit_diagnostic"]:
                    row["observable_insensitivity_verdict"] = None
            rows[sector] = row
            reproduction[sector] = compare_phases(old_phases[sector], table, contract)
        groups[candidate] = rows
    return {
        "verified_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "reader_source_sha256": sha(Path(__file__)),
        "contract": contract,
        "native_owner_manifest": manifest,
        "matched_engine_provenance": owners,
        "fresh_cas11_reading": fresh,
        "previous_raw_file_manifest": previous["input_manifest"],
        "previous_input_reading": previous["equilibrium_input_reading"],
        "previous_numerical_controls_passed": previous["numerical_sensitivity_controls_passed"],
        "old_cas11_pi_phase_difference_rad_mod_pi": old11["pi_phase_difference_rad_mod_pi"],
        "old_to_fresh_cas11": reproduction,
        "outer_deck_audit": deck_audit,
        "original_strict_byte_deck_prerequisite_verdict": False,
        "original_strict_observable_insensitivity_verdict": None,
        "same_engine_active_space_comparisons": groups,
        "new_native_solves": 0,
        "electronic_model_verdict": None,
        "spectral_truncation_verdict": None,
        "global_continuity_verdict": None,
        "full_sweep_released": False,
        "scope": "Saved-output active-space sensitivity with matched audited scattering engine "
        "and eight MPI ranks. Original target diagonalizers/imports, optimized orbital frames "
        "and native target-file representations differ; this is not a pure orbital-count "
        "effect, full model convergence, absolute omission bound or pole identity.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("cas11", "saved", "provenance", "contract", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.cas11, args.saved, args.provenance, json.loads(args.contract.read_text()))
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
