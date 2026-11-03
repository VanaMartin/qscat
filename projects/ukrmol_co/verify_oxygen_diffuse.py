"""Reconstruct the original failed diffuse diagnostic; never infer an affinity."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

import numpy as np

from projects.ukrmol_co.verify_stretched_properties import files_in, require, sha

SOURCE = "2b21b43a33d94a0fb4157dc7d8a4c4a1ace4da24"
ARCHIVE = "7ba5b429be6a64a84b32c63498e126c7167d69b0cc45c9360000a7e4b4e50645"
SOURCE_ARCHIVE = "3925cccad754db528128d1bc141b6de396ddf7a71748147205e9465d38858fe8"
CONTRACT = "projects/ukrmol_co/oxygen-affinity-diffuse-release-20261103.json"


def verify(path, collection):
    raw = path.read_bytes()
    require(sha(raw) == ARCHIVE == collection["archive"]["sha256"], "Original failed archive")
    kept = files_in(raw)
    require(
        set(kept) == set(collection["files"]) and len(kept) == 26, "Full partial-evidence inventory"
    )
    for name, pin in collection["files"].items():
        require(sha(kept[name]) == pin["sha256"] and len(kept[name]) == pin["bytes"], name)
    require(sha(kept["source.tar.gz"]) == SOURCE_ARCHIVE, "Frozen source archive")
    source = files_in(kept["source.tar.gz"])
    pins = json.loads(kept["source-manifest.json"])
    require(set(source) == set(pins) and len(source) == 15, "Executed source inventory")
    for name, pin in pins.items():
        require(
            sha(source[name]) == pin
            and source[name] == subprocess.check_output(["git", "show", SOURCE + ":" + name]),
            "Executed Git source " + name,
        )
    c = json.loads(source[CONTRACT])
    for name, pin in c["source_sha256"].items():
        require(sha(source[name]) == pin, "Frozen calculation prerequisite " + name)
    owner = json.loads(kept["execution.json"])
    container = json.loads(kept["container-inspect.json"])
    science = json.loads(kept["outputs/result.json"])
    exits = [
        json.loads(kept["actual-image-controls.json"])["original_exit"],
        json.loads(kept["outputs/controls.json"])["original_exit"],
        science["original_science_exit"],
        owner["original_command_exit"],
        owner["original_container_exit"],
        owner["original_owner_exit"],
        json.loads(kept["launcher-exit.json"])["original_exit"],
    ]
    require(
        exits == [0, 0, 2, 2, 2, 2, 2]
        and not owner["timed_out"]
        and not owner["container_oom_killed"]
        and owner["source_integrity_after"]
        and container["State"]["ExitCode"] == 2
        and not container["State"]["OOMKilled"],
        "Original failed exits/time/OOM/integrity",
    )
    require(
        owner["source_commit"] == SOURCE
        and owner["source_archive_sha256"] == SOURCE_ARCHIVE
        and container["Image"] == c["image"],
        "Executed image/source",
    )
    require(
        "Ran 4 tests" in kept["outputs/controls.log"].decode()
        and [r["original_profile_exit"] for r in science["states"]] == [0, 1],
        "Controls/workers",
    )
    neutral = json.loads(kept["outputs/O/atom.json"])
    failed = json.loads(kept["outputs/O-/diagnostics.json"])
    require(
        neutral["nao"] == failed["nao"] == 105 and neutral["basis"] == failed["basis"],
        "Same extended basis",
    )
    require(
        neutral["rohf_converged"]
        and neutral["rohf_internally_stable"]
        and neutral["uccsd_converged"]
        and neutral["finite_method_energy_complete"],
        "Completed neutral only",
    )
    log = kept["outputs/O/run.log"].decode()
    scf = float(re.search(r"converged SCF energy =\s*([-+.\deE]+)", log)[1])
    cc = float(re.search(r"E\(UCCSD\) =\s*([-+.\deE]+)", log)[1])
    triples = float(re.search(r"UCCSD\(T\) correction =\s*([-+.\deE]+)", log)[1])
    np.testing.assert_allclose(
        [scf, cc, cc + triples],
        [
            neutral["rohf_energy_hartree"],
            neutral["uccsd_energy_hartree"],
            neutral["total_energy_hartree"],
        ],
        atol=1e-11,
        rtol=0,
    )
    failed_log = kept["outputs/O-/run.log"].decode()
    require(
        not failed["rohf_converged"]
        and len(failed["stages"]) == 1
        and "cycle= 100" in failed_log
        and "SCF not converged." in failed_log
        and "Atomic ROHF did not converge" in failed_log
        and "UCCSD" not in failed_log
        and "outputs/O-/atom.json" not in kept
        and not any("affinity" in name for name in science),
        "No anion correlation or affinity",
    )
    resources = [json.loads(kept["outputs/" + s + "/resources.json"]) for s in ["O", "O-"]]
    require(
        [r["exit_code"] for r in resources] == [0, 1]
        and max(r["kernel_peak_container_memory_bytes"] for r in resources) < 16 * 2**30,
        "Original resource records",
    )
    observations = [json.loads(row) for row in kept["monitor.jsonl"].decode().splitlines()]
    require(
        len(observations) == 1
        and observations[0]["original_command_exit"] == 2
        and owner["wall_seconds"] < 600,
        "Immediate terminal failure observation",
    )
    return dict(
        source_commit=SOURCE,
        archive_sha256=ARCHIVE,
        raw_files_verified=26,
        executed_git_sources_verified=15,
        original_exits=exits,
        original_worker_exits=[0, 1],
        neutral_total_energy_hartree=neutral["total_energy_hartree"],
        anion_rohf_converged=False,
        anion_cc_calls=0,
        electron_affinity_ev=None,
        diffuse_shift_ev=None,
        failure_requalified=False,
        owner_wall_seconds=owner["wall_seconds"],
        worker_cpu_seconds=sum(r["container_cpu_seconds"] for r in resources),
        kernel_peak_container_memory_bytes=max(
            r["kernel_peak_container_memory_bytes"] for r in resources
        ),
        no_quantum_calls_in_reconstruction=True,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive", type=Path)
    p.add_argument("collection", type=Path)
    a = p.parse_args()
    print(json.dumps(verify(a.archive, json.loads(a.collection.read_text())), indent=2))
