"""Saved-only atomic log/ledger/resource reconstruction; never import the producer."""

from __future__ import annotations

import argparse
import csv
import datetime
import io
import itertools
import json
import math
import re
import subprocess
from pathlib import Path

from qscat.units import HARTREE_TO_EV

from projects.ukrmol_co.verify_stretched_properties import files_in, require, sha

SOURCE = "493856bd83751e36d66596157e8cbfc9a74efb5a"
SOURCE_ARCHIVE = "eccbc700815716d20dae586553b064df130c009bb9cca1940ce80ebbdd11570b"
ARCHIVE = "2908950e1ded0c448fc83daf3162ce5b2ae0df52887cb77708fbfbd60bb0a842"
NUMBER = r"([-+]?\d+(?:\.\d*)?(?:[eEdD][-+]?\d+)?)"


def log_energies(text):
    values = []
    for prefix in [r"^converged SCF energy = ", r"^E\(UCCSD\) = ", r"^UCCSD\(T\) correction = "]:
        matches = re.findall(prefix + NUMBER, text, re.M)
        require(len(matches) == 1, "Unique converged native atomic energy " + prefix)
        values.append(float(matches[0].replace("D", "E")))
    require(all(math.isfinite(v) for v in values), "Finite raw atomic energies")
    require(re.search(r"^UCCSD converged$", text, re.M) is not None, "Raw CC convergence")
    require(
        re.search(
            r"^<class .+ wavefunction is stable in the internal stability analysis$", text, re.M
        )
        is not None,
        "Raw internal ROHF stability",
    )
    return dict(
        rohf_energy_hartree=values[0],
        uccsd_energy_hartree=values[1],
        triples_correction_hartree=values[2],
        total_energy_hartree=values[1] + values[2],
    )


def verify(
    path,
    collection,
    *,
    source_commit=SOURCE,
    source_archive_sha256=SOURCE_ARCHIVE,
    archive_sha256=ARCHIVE,
    contract_path="projects/ukrmol_co/atomic-asymptote-release-20261103.json",
    basis="aug-cc-pVQZ",
    source_entries=14,
    raw_entries=47,
    runtime_controls=6,
):
    raw = path.read_bytes()
    require(
        sha(raw) == archive_sha256 == collection["archive"]["sha256"],
        "Exact original atomic archive",
    )
    kept = files_in(raw)
    require(
        set(kept) == set(collection["files"]) and len(kept) == raw_entries,
        "Full original file inventory",
    )
    for n, pin in collection["files"].items():
        require(
            len(kept[n]) == pin["bytes"] and sha(kept[n]) == pin["sha256"],
            "Full original bytes " + n,
        )

    def load(name):
        return json.loads(kept[name])

    source = files_in(kept["source.tar.gz"])
    pins = load("source-manifest.json")
    require(
        set(source) == set(pins)
        and len(source) == source_entries
        and sha(kept["source.tar.gz"]) == source_archive_sha256,
        "Exact executed source archive",
    )
    for n, h in pins.items():
        require(
            sha(source[n]) == h
            and source[n] == subprocess.check_output(["git", "show", source_commit + ":" + n]),
            "Original executed Git source " + n,
        )
    c = json.loads(source[contract_path])
    inputs = json.loads(source[c["inputs"]])
    owner = load("execution.json")
    container = load("container-inspect.json")
    result = load("outputs/result.json")
    controls = load("outputs/controls.json")
    exits = [
        controls["original_exit"],
        result["scientific_exit"],
        owner["original_command_exit"],
        owner["original_container_exit"],
        owner["original_owner_exit"],
        container["State"]["ExitCode"],
        load("service-return.json")["original_service_exit"],
        load("launcher-exit.json")["original_exit"],
    ]
    require(
        exits == [0] * 8 and not owner["timed_out"] and not container["State"]["OOMKilled"],
        "Original exits/time/OOM",
    )
    require(
        container["Image"] == c["image"]
        and owner["source_commit"] == source_commit
        and owner["source_integrity_after"]
        and f"Ran {runtime_controls} tests".encode() in kept["outputs/controls.log"],
        "Image/source/runtime controls",
    )
    host = container["HostConfig"]
    mounts = {m["Destination"]: m for m in container["Mounts"]}
    require(
        host["CpusetCpus"] == "0-3"
        and host["Memory"] == host["MemorySwap"] == 16 * 2**30
        and host["NetworkMode"] == "none"
        and not mounts["/opt/qmodeling"]["RW"],
        "Original allocation",
    )
    for n, h in c["source_sha256"].items():
        require(sha(source[n]) == h, "Original source prerequisite " + n)
    totals = {}
    resources = {}
    samples = 0
    for state, row in zip(inputs["states"], result["atoms"], strict=True):
        species = state["id"]
        prefix = "outputs/" + species + "/"
        require(
            row == dict(id=species, state="finished", original_profile_exit=0),
            "Original atom exit " + species,
        )
        reported = load(prefix + "atom.json")
        native = log_energies(kept[prefix + "run.log"].decode())
        for k, value in native.items():
            require(abs(reported[k] - value) <= 1e-10, "Raw log arithmetic " + species + "/" + k)
        require(
            reported["basis"] == basis
            and reported["pyscf_version"] == "2.11.0"
            and reported["external_rohf_stability"] is None
            and reported["rohf_internally_stable"]
            and reported["frozen_spatial_orbitals"] == 1
            and [s["name"] for s in reported["stages"]]
            == ["rohf", "rohf_internal_stability", "uccsd", "triples"]
            and all(s["original_exit"] == 0 for s in reported["stages"]),
            "Recorded state/method gates " + species,
        )
        spin = state["spin_2s"] / 2
        require(
            abs(reported["reference_s2"] - spin * (spin + 1)) <= 1e-8,
            "Recorded reference spin " + species,
        )
        require(
            all(r["target_weight"] > 0.5 for r in reported["occupied_angular_character"])
            and reported["minimum_overlap_eigenvalue"] > 1e-9,
            "Recorded angular/rank gates " + species,
        )
        r = load(prefix + "resources.json")
        rows = list(csv.DictReader(io.StringIO(kept[prefix + "resources.csv"].decode())))
        require(
            r["exit_code"] == 0 and r["wall_seconds"] < 1800 and len(rows) > 2,
            "Resource/time bound " + species,
        )
        peak = max(int(x["container_memory_bytes"]) for x in rows)
        require(
            peak == r["sampled_peak_container_memory_bytes"]
            and r["kernel_peak_container_memory_bytes"] >= peak,
            "Raw resource samples " + species,
        )
        totals[species] = native["total_energy_hartree"]
        resources[species] = r
        samples += len(rows)
    sums = {
        "C+O": totals["C"] + totals["O"],
        "C+O-": totals["C"] + totals["O-"],
        "C-+O": totals["C-"] + totals["O"],
    }
    eas = {"O": totals["O"] - totals["O-"], "C": totals["C"] - totals["C-"]}
    for k, v in sums.items():
        require(abs(v - result["fragment_sum_hartree"][k]) <= 2e-10, "Fragment identity " + k)
    for k, v in eas.items():
        require(abs(v - result["electron_affinity_hartree"][k]) <= 2e-10, "EA identity " + k)
    observations = [json.loads(s) for s in kept["monitor.jsonl"].decode().splitlines()]
    times = [datetime.datetime.fromisoformat(owner["started_utc"])] + [
        datetime.datetime.fromisoformat(r["observed_utc"]) for r in observations
    ]
    gaps = [(b - a).total_seconds() for a, b in itertools.pairwise(times)]
    require(
        gaps and min(gaps) > 0 and observations[-1]["container_command_exit"] == 0,
        "Terminal monitor observation",
    )
    return dict(
        finite_method_endpoint_reconstruction_verdict=True,
        source_commit=source_commit,
        archived_files_verified=len(kept),
        executed_git_sources_verified=len(source),
        source_prerequisites_verified=len(c["source_sha256"]),
        resource_samples_verified=samples,
        original_exits=exits,
        original_atom_profile_exits=[0] * 4,
        total_energy_hartree=totals,
        fragment_sum_hartree=sums,
        electron_affinity_hartree=eas,
        electron_affinity_ev={k: v * HARTREE_TO_EV for k, v in eas.items()},
        resources=resources,
        science_wall_seconds=result["wall_seconds"],
        owner_wall_seconds=owner["wall_seconds"],
        kernel_peak_container_memory_bytes=max(
            r["kernel_peak_container_memory_bytes"] for r in resources.values()
        ),
        maximum_observation_gap_seconds=max(gaps),
        within_requested_900second_observation_bound=max(gaps) <= 900,
        finished_before_first_600second_periodic_check=owner["wall_seconds"] < 600,
        atomic_absolute_accuracy_qualified=False,
        molecular_dissociation_qualified=False,
        full_correlated_L2_purity_proved=False,
        new_quantum_calls_during_verification=0,
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
        independence="Full payload hashes and executed Git bytes;"
        "raw SCF/CC/triples parsing and independent endpoint arithmetic;resource CSV replay. "
        "State/spin/occupation gates are source-authenticated producer diagnostics,"
        "not fresh wavefunction or integral calculations.",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive", type=Path)
    p.add_argument("collection", type=Path)
    a = p.parse_args()
    print(
        json.dumps(
            verify(a.archive, json.loads(a.collection.read_text())), indent=2, allow_nan=False
        )
    )
