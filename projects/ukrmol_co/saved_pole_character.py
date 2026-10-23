"""Saved-only first stop of the compressed pair contract; never propagates/searches."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import resource
import tarfile
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.pole_character import screen_witnesses
from projects.ukrmol_co.saved_matrices import read_interface


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def verify_inputs(root, contract):
    for name, want in contract["input_sha256"].items():
        if digest(root / name) != want:
            raise ValueError(f"Changed prerequisite: {name}")
    raw = root / contract["raw_directory"]
    manifest = json.loads((raw / "input-manifest.json").read_text())
    if len(manifest) != 123:
        raise ValueError("Changed original input inventory")
    for name, want in manifest.items():
        if digest(raw / name) != want:
            raise ValueError(f"Changed original input: {name}")
    return len(manifest)


def load_b1(root, contract):
    """Authenticate archived evidence; do not extract paths or execute producer code."""
    collection = json.loads((root / contract["b1_collection"]).read_text())
    if digest(root / contract["b1_archive"]) != collection["archive_sha256"]:
        raise ValueError("B1 archive/collection disagree")
    names = ["outputs/check/result.json", "execution.json"]
    names += [
        f"outputs/check/{name}.npz"
        for name in ["moments-0", "moments-1", "derivative-0", "derivative-1"]
    ]
    selected = {}
    seen = set()
    with tarfile.open(root / contract["b1_archive"]) as archive:
        for member in archive:
            name = member.name
            if not member.isfile() or name in seen or name not in collection["files"]:
                raise ValueError("Unexpected/duplicate B1 archive member")
            seen.add(name)
            data = archive.extractfile(member).read()
            want = collection["files"][name]
            if len(data) != want["bytes"] or hashlib.sha256(data).hexdigest() != want["sha256"]:
                raise ValueError(f"Changed B1 member: {name}")
            if name in names:
                selected[name] = data
    if seen != set(collection["files"]) or set(selected) != set(names):
        raise ValueError("Incomplete B1 archive")
    result = json.loads(selected[names[0]])
    owner = json.loads(selected[names[1]])
    if (
        result["finite_reference_pole_verdict"] is not True
        or result["case"] != "compressed/nine/B1"
        or result["scientific_exit"] != 0
        or not all(result["contour_witness"]["checks"].values())
        or not all(all(row["checks"].values()) for row in result["derivative_witnesses"])
        or [
            owner[n]
            for n in ["original_command_exit", "original_container_exit", "original_owner_exit"]
        ]
        != [0, 0, 0]
        or not owner["source_integrity_after"]
        or owner["timed_out"]
        or owner["container_oom_killed"]
    ):
        raise ValueError("B1 finite pole/residue prerequisite does not pass")
    pole = complex(*result["pole_incident_hartree"])
    np.testing.assert_allclose(
        pole, complex(*contract["b1_pole_incident_hartree"]), atol=1e-15, rtol=0
    )
    residues = []
    for name in names[2:]:
        with np.load(io.BytesIO(selected[name]), allow_pickle=False) as capture:
            residues.append(capture["residue"])
    return pole, residues, len(seen)


def verified_labels(root, contract):
    labels = {}
    for sector in ["B1", "B2"]:
        source = root / contract["raw_directory"] / "compressed/nine" / sector
        interface = read_interface(source / "fort.10", source / "fort.21")
        thresholds = np.array(interface["channel_thresholds_rydberg"])
        opened = np.flatnonzero(thresholds == 0)
        if (
            interface["target_labels"][0] != contract["ground_target_label"]
            or opened.tolist() != list(range(6))
            or np.any(thresholds[6:] <= 0)
            or interface["channel_labels"][:6] != contract["open_labels"][sector]
            or interface["boundary_channels"] != 298
            or interface["boundary_poles"] != 2048
        ):
            raise ValueError(f"Unqualified ground-target/open-label scope: {sector}")
        labels[sector] = interface["channel_labels"][:6]
    return labels


def run(root, output, contract):
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    report = dict(
        stage="saved B1 necessary character screen only",
        paired_verdict=None,
        b1_character_disposition=None,
        scientific_exit=2,
        new_propagations=0,
        native_calls=0,
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
    )
    try:
        report["inputs_verified_before"] = verify_inputs(root, contract)
        labels = verified_labels(root, contract)
        pole, residues, count = load_b1(root, contract)
        np.savez(
            output / "saved-residues.npz",
            **{f"residue_{i}": value for i, value in enumerate(residues)},
        )
        result = screen_witnesses(residues, labels["B1"], -2 * pole.imag, contract["character"])
        report.update(
            archived_files_verified=count,
            labels=labels,
            b1_screen=result,
            b1_character_disposition=result["disposition"],
            scientific_exit={"pass": 0, "reject": 1, "inconclusive": 2}[result["disposition"]],
            b2_calculation_eligible=result["disposition"] == "pass",
            b2_execution_released=False,
        )
    except Exception as error:
        report.update(error=repr(error), scientific_exit=2)
    finally:
        try:
            report["inputs_verified_after"] = verify_inputs(root, contract)
        except Exception as error:
            report.update(
                integrity_error=repr(error), scientific_exit=2, b2_calculation_eligible=False
            )
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
        (output / "result.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report["scientific_exit"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(run(args.root, args.output, json.loads(args.contract.read_text())))
