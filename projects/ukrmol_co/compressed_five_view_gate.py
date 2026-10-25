"""One saved-only qualification of the compressed mapped view; no S/pole evaluation."""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.compressed_five_view import load_compressed_five
from projects.ukrmol_co.native_full_reconstruction import _model_from_interface, load_model
from projects.ukrmol_co.saved_matrices import read_interface, records
from projects.ukrmol_co.saved_pole_character import digest, verify_inputs


def run(root, output, contract):
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    report = dict(
        mapped_view_verdict=None,
        scientific_exit=2,
        cases=[],
        new_propagations=0,
        new_pole_evaluations=0,
        native_calls=0,
        route_holdout_s_matrices_analyzed=0,
        pole_calculation_released=False,
        model_insensitivity_verdict=None,
        full_sweep_released=False,
    )

    def save():
        usage = resource.getrusage(resource.RUSAGE_SELF)
        report.update(
            wall_seconds=time.monotonic() - started,
            process_cpu_seconds=usage.ru_utime + usage.ru_stime,
            peak_process_rss_platform_units=usage.ru_maxrss,
        )
        for name in ["memory.peak", "cpu.stat"]:
            path = Path("/sys/fs/cgroup") / name
            if path.is_file():
                report["cgroup_" + name] = path.read_text()
        temporary = output / "result.partial"
        temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        temporary.replace(output / "result.json")

    try:
        proposal = Path(__file__).with_name(contract["source_proposal"])
        if digest(proposal) != contract["source_proposal_sha256"]:
            raise ValueError("Changed source proposal")
        report["inputs_verified_before"] = verify_inputs(root, contract)
        previous = json.loads((root / contract["original_strict_report"]).read_text())
        if previous["saved_matrix_reproduction_verdict"] is not False:
            raise ValueError("Original strict aggregate failure must be preserved")
        raw = root / contract["raw_directory"]
        for sector in ["B1", "B2"]:
            case = "compressed/five/" + sector
            prior = next(row for row in previous["cases"] if row["case"] == case)
            if prior["saved_matrix_reproduction_verdict"] is not None:
                raise ValueError("Original five-case verdict must remain unset")
            source = raw / case
            row = dict(case=case, original_strict_case_verdict=None, completed=False)
            report["cases"].append(row)
            save()
            try:
                read_interface(source / "fort.10", source / "fort.21")
            except ValueError as error:
                if str(error) != "Invalid channel label":
                    raise
                row["strict_default_failure"] = repr(error)
            else:
                raise ValueError("The default reader unexpectedly accepted original property IDs")
            model = load_compressed_five(source)
            boundary = records(source / "fort.21")
            np.testing.assert_array_equal(model["poles"], np.frombuffer(boundary[6], "<f8"))
            np.testing.assert_array_equal(
                model["amplitudes"], np.frombuffer(boundary[7], "<f8").reshape(2048, 250)
            )
            reference = raw / "compressed/nine" / sector
            strict = load_model(reference)
            interface = read_interface(
                reference / "fort.10",
                reference / "fort.21",
                selected_property_ids=list(range(1, 49)),
            )
            identity = _model_from_interface(reference, interface, (48, 298))
            for key in strict:
                if isinstance(strict[key], np.ndarray):
                    np.testing.assert_array_equal(strict[key], identity[key])
                elif strict[key] != identity[key]:
                    raise ValueError(f"Nine-model identity view changed {key}")
            row.update(
                completed=True,
                view_evidence=model["view_evidence"],
                original_property_ids=model["interface"]["selected_property_ids"],
                channel_target_rows=model["interface"]["channel_target_rows"],
                raw_channel_labels=model["interface"]["channel_labels"],
                boundary_channels=len(model["thresholds"]),
                closed_channels=len(model["thresholds"]) - 6,
                roots=len(model["poles"]),
                own_w_d_byte_exact=True,
                nine_identity_all_model_fields_exact=True,
                model_array_sha256={
                    key: hashlib.sha256(model[key].tobytes()).hexdigest()
                    for key in ["poles", "amplitudes", "thresholds", "angular", "c2", "c3"]
                },
            )
            save()
        report.update(
            mapped_view_verdict=True, scientific_exit=0, original_strict_aggregate_verdict=False
        )
    except Exception as error:
        report.update(error=repr(error), scientific_exit=2)
    finally:
        try:
            report["inputs_verified_after"] = verify_inputs(root, contract)
        except Exception as error:
            report.update(integrity_error=repr(error), scientific_exit=2, mapped_view_verdict=None)
        save()
    return report["scientific_exit"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(run(args.root, args.output, json.loads(args.contract.read_text())))
