"""Authenticate the active property launch without turning placeholders into exits."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

from projects.ukrmol_co.verify_stretched_properties import files_in, require, sha

SOURCE = "58a61ba052b29a583cf8360d9b2d2620ae60e989"
ARCHIVE = "e7ea1f68fe6a4ec9134dc00c35c1155c14e972e3efc2f316c5f6274b066a211a"
CONTRACT = "projects/ukrmol_co/stretched-augtz-properties-contract-20261102.json"


def launch_disposition(owner, science, container):
    require(
        owner["state"] == "running"
        and container["State"]["Running"]
        and not container["State"]["OOMKilled"]
        and science["property_inventory_verdict"] is None
        and science["native_denprop_calls"] == 1
        and science["stages"] == [],
        "Active one-call launch",
    )
    require(
        all(
            science[n] == 0
            for n in [
                "qc_calls",
                "target_diagonalizations",
                "scattering_calls",
                "new_maps",
                "route_holdout_s_matrices_analyzed",
            ]
        )
        and not science["full_sweep_released"],
        "No automatic successor/protected work",
    )
    return dict(
        property_inventory_verdict=None,
        original_terminal_exits=None,
        scope=(
            "Active launch only;producer/owner initial2 fields are not "
            "terminal failures or successes"
        ),
    )


def verify(archive, collection):
    require(
        len(archive) == collection["archive"]["bytes"]
        and sha(archive) == collection["archive"]["sha256"],
        "Launch archive bytes/hash",
    )
    files = files_in(archive)
    require(
        len(files) == 14 and set(files) == set(collection["files"]),
        "Exact launch evidence inventory",
    )
    for name, pin in collection["files"].items():
        require(
            sha(files[name]) == pin["sha256"] and len(files[name]) == pin["bytes"],
            "Launch bytes " + name,
        )

    def load(name):
        return json.loads(files[name])

    source = files_in(files["source.tar.gz"])
    pins = load("source-manifest.json")
    require(
        sha(files["source.tar.gz"]) == ARCHIVE and len(source) == 20 and set(source) == set(pins),
        "Exact executing20-source capsule",
    )
    for name, pin in pins.items():
        require(
            sha(source[name]) == pin
            and source[name] == subprocess.check_output(["git", "show", SOURCE + ":" + name]),
            "Executing Git source " + name,
        )
    c = json.loads(source[CONTRACT])
    model = c["models"]["augtz"]
    require(
        c["execution_authorized"]
        and len(model["input_sha256"]) == 18
        and len(c["source_sha256"]) == 6,
        "One reviewed18-input/six-prerequisite release",
    )
    for name, pin in c["source_sha256"].items():
        require(sha(source[name]) == pin, "Source prerequisite " + name)
    preflight = load("input-preflight.json")
    require(
        preflight["original_exit"] == 0
        and json.loads(preflight["stdout"])
        == dict(source_prerequisites_verified=6, raw_input_pins_verified=18, new_native_calls=0),
        "Original host input hash-gate attestation",
    )
    owner, science, container = [
        load(n)
        for n in ["execution.json", "outputs/properties/result.json", "container-inspect.json"]
    ]
    disposition = launch_disposition(owner, science, container)
    require(
        owner["source_commit"] == SOURCE
        and owner["source_archive_sha256"] == ARCHIVE
        and owner["owner_sha256"]
        == sha(source["projects/ukrmol_co/stretched_basis_properties_owner.py"]),
        "Original executing owner provenance",
    )
    require(
        load("outputs/controls.json")["original_exit"] == 0
        and b"Ran 10 tests" in files["outputs/controls.log"],
        "All10 original runtime controls",
    )
    require(len(science["copied_input_sha256"]) == 12, "Original12 copied input pins")
    for name, pin in science["copied_input_sha256"].items():
        require(
            pin == model["input_sha256"]["output/CO/geom1/" + name], "Own unchanged input " + name
        )
    deck = files["outputs/properties/target.denprop.inp"]
    require(
        sha(deck) == science["deck_sha256"] == c["matched_deck_sha256"],
        "Exact executed64-state deck",
    )
    active = "\n".join(line.split("!", 1)[0] for line in deck.decode().splitlines())
    for name, want in [
        ("ntgtl", [8] * 8),
        ("ntgtf", list(range(1, 9))),
        ("nftsor", [710, 730, 711, 731, 712, 732, 713, 733]),
    ]:
        fields = re.findall(rf"\b{name}\s*=\s*([\d, ]+)", active, re.I)
        require(
            len(fields) == 1 and list(map(int, re.findall(r"\d+", fields[0]))) == want,
            "Native property selection " + name,
        )
    host = container["HostConfig"]
    mounts = {m["Destination"]: m for m in container["Mounts"]}
    require(
        container["Image"] == c["image"]
        and host["CpusetCpus"] == "4-11"
        and host["Memory"] == host["MemorySwap"] == 8 * 2**30
        and host["NetworkMode"] == "none"
        and set(mounts) == {"/work", "/opt/qmodeling", "/out"}
        and not mounts["/work"]["RW"]
        and not mounts["/opt/qmodeling"]["RW"],
        "Immutable8-core/8GiB launch",
    )
    return dict(
        observed_utc=collection["observed_utc"],
        started_utc=owner["started_utc"],
        source_commit=SOURCE,
        source_archive_sha256=ARCHIVE,
        launch_raw_files_verified=14,
        executed_git_sources_verified=20,
        source_prerequisites_verified=6,
        raw_input_hash_gate_attested=18,
        copied_input_pins_reconciled=12,
        runtime_controls_passed=10,
        native_denprop_calls=1,
        qc_calls=0,
        target_diagonalizations=0,
        scattering_calls=0,
        new_maps=0,
        route_holdout_s_matrices_analyzed=0,
        automatic_successor=False,
        full_sweep_released=False,
        **disposition,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive", type=Path)
    p.add_argument("collection", type=Path)
    a = p.parse_args()
    print(
        json.dumps(verify(a.archive.read_bytes(), json.loads(a.collection.read_text())), indent=2)
    )
