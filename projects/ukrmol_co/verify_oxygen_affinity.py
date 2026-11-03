"""Reconstruct saved oxygen orbital-control evidence without any quantum calls."""

from __future__ import annotations

import argparse
import io
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from qscat.units import HARTREE_TO_EV

from projects.ukrmol_co.verify_stretched_properties import files_in, require, sha

SOURCE = "fb063d9065df7f329a801cd9b4cf4150b6303d93"
ARCHIVE = "862cbe2aecc1843aa497b6a03641f761ebea071954e4158f525d34cc65459262"
SOURCE_ARCHIVE = "39b48c8653ae06880f227d22272c548511725b9bd302782525f022a4a8e7081a"
CONTRACT = "projects/ukrmol_co/oxygen-affinity-orbital-release-20261103.json"


def verify(path, collection, original_path):
    raw = path.read_bytes()
    require(sha(raw) == ARCHIVE == collection["archive"]["sha256"], "Original archive digest")
    kept = files_in(raw)
    require(set(kept) == set(collection["files"]) and len(kept) == 24, "Full raw inventory")
    for name, pin in collection["files"].items():
        require(sha(kept[name]) == pin["sha256"] and len(kept[name]) == pin["bytes"], name)
    require(sha(kept["source.tar.gz"]) == SOURCE_ARCHIVE, "Frozen source archive")
    source = files_in(kept["source.tar.gz"])
    pins = json.loads(kept["source-manifest.json"])
    require(set(source) == set(pins) and len(source) == 14, "Full executed source inventory")
    for name, pin in pins.items():
        require(
            sha(source[name]) == pin
            and source[name] == subprocess.check_output(["git", "show", SOURCE + ":" + name]),
            "Executed Git source " + name,
        )
    c = json.loads(source[CONTRACT])
    for name, pin in c["source_sha256"].items():
        require(sha(source[name]) == pin, "Calculation prerequisite " + name)
    original_raw = original_path.read_bytes()
    require(
        sha(original_raw) == "3a2521b4b8d12b6c4a816cfb51e15917f5dafbfce1777ddbeb5a1626cfe5498d",
        "Original qualified 5Z archive",
    )
    old = files_in(original_raw)
    for name, pin in c["input_sha256"].items():
        require(sha(old[name]) == pin, "Original immutable input " + name)
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
        exits == [0] * 7
        and not owner["timed_out"]
        and not owner["container_oom_killed"]
        and owner["source_integrity_after"]
        and container["State"]["ExitCode"] == 0
        and not container["State"]["OOMKilled"],
        "Original exits/time/OOM/integrity",
    )
    require(
        owner["source_commit"] == SOURCE
        and owner["source_archive_sha256"] == SOURCE_ARCHIVE
        and container["Image"] == c["image"],
        "Executed source/image pins",
    )
    require("Ran 4 tests" in kept["outputs/controls.log"].decode(), "Runtime controls")
    atoms, resources = {}, []
    for species in ["O", "O-"]:
        prefix = "outputs/" + species + "/"
        a = json.loads(kept[prefix + "result.json"])
        original = json.loads(old[prefix + "atom.json"])
        log = kept[prefix + "run.log"].decode()
        triples = [float(x) for x in re.findall(r"UCCSD\(T\) correction =\s*([-+.\deE]+)", log)]
        np.testing.assert_allclose(
            triples,
            [a["original_triples_hartree"], a["semicanonical_triples_hartree"]],
            atol=1e-14,
            rtol=0,
        )
        np.testing.assert_allclose(
            a["original_total_hartree"], original["total_energy_hartree"], atol=1e-10, rtol=0
        )
        np.testing.assert_allclose(
            a["original_ccsd_correlation_hartree"],
            a["rotated_ccsd_correlation_hartree"],
            atol=1e-8,
            rtol=0,
        )
        for kind, corr in [("original", "original"), ("semicanonical", "rotated")]:
            total = (
                a["reference_energy_hartree"]
                + a[corr + "_ccsd_correlation_hartree"]
                + a[kind + "_triples_hartree"]
            )
            np.testing.assert_allclose(total, a[kind + "_total_hartree"], atol=1e-12, rtol=0)
        with np.load(io.BytesIO(kept[prefix + "orbital-witnesses.npz"]), allow_pickle=False) as w:
            np.testing.assert_allclose(w["density_new"], w["density_old"], atol=1e-10, rtol=0)
            for spin, n in zip("ab", w["nocc"], strict=True):
                u, f = w["rotation" + spin], w["fock" + spin]
                np.testing.assert_allclose(u.T @ u, np.eye(len(u)), atol=1e-10, rtol=0)
                np.testing.assert_allclose(u[:n, n:], 0, atol=1e-14, rtol=0)
                rotated = u.T @ f @ u
                for block in [rotated[:n, :n], rotated[n:, n:]]:
                    np.testing.assert_allclose(block, np.diag(np.diag(block)), atol=1e-9, rtol=0)
                np.testing.assert_allclose(
                    w["t1" + spin + "_new"],
                    u[:n, :n].T @ w["t1" + spin + "_old"] @ u[n:, n:],
                    atol=1e-12,
                    rtol=0,
                )
        require(a["new_scf_calls"] == a["new_ccsd_iterations"] == 0, "Saved-state-only scope")
        r = json.loads(kept[prefix + "resources.json"])
        require(
            r["exit_code"] == 0 and r["kernel_peak_container_memory_bytes"] < 16 * 2**30,
            "Original worker exit/resources",
        )
        resources.append(r)
        atoms[species] = a
    affinity = {
        kind: (atoms["O"][kind + "_total_hartree"] - atoms["O-"][kind + "_total_hartree"])
        * HARTREE_TO_EV
        for kind in ["original", "semicanonical"]
    }
    for kind, value in affinity.items():
        np.testing.assert_allclose(
            value, science[kind + "_electron_affinity_ev"], atol=1e-11, rtol=0
        )
    shift = affinity["semicanonical"] - affinity["original"]
    np.testing.assert_allclose(shift, science["orbital_convention_ea_shift_ev"], atol=1e-11, rtol=0)
    observations = [json.loads(row) for row in kept["monitor.jsonl"].decode().splitlines()]
    require(
        len(observations) == 1
        and observations[0]["original_command_exit"] == 0
        and owner["wall_seconds"] < 600,
        "Immediate terminal observation before first periodic check",
    )
    return dict(
        source_commit=SOURCE,
        archive_sha256=ARCHIVE,
        raw_files_verified=24,
        executed_git_sources_verified=14,
        original_exits=exits,
        original_electron_affinity_ev=affinity["original"],
        semicanonical_electron_affinity_ev=affinity["semicanonical"],
        orbital_convention_ea_shift_ev=shift,
        meaningful_at_0_1_mev=abs(shift) > 0.0001,
        owner_wall_seconds=owner["wall_seconds"],
        worker_cpu_seconds=sum(r["container_cpu_seconds"] for r in resources),
        kernel_peak_container_memory_bytes=max(
            r["kernel_peak_container_memory_bytes"] for r in resources
        ),
        no_quantum_calls_in_reconstruction=True,
        absolute_accuracy_qualified=False,
        double_amplitude_rotation_checked_by_executed_ccsd_energy_only=True,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive", type=Path)
    p.add_argument("collection", type=Path)
    p.add_argument("original_5z_archive", type=Path)
    a = p.parse_args()
    print(
        json.dumps(
            verify(a.archive, json.loads(a.collection.read_text()), a.original_5z_archive), indent=2
        )
    )
