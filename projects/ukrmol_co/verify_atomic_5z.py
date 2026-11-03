"""Reconstruct the matched5Z result or original zero-atom stop without quantum calls."""

from __future__ import annotations

import argparse
import ast
import json
import subprocess
from pathlib import Path

from projects.ukrmol_co.atomic_basis_comparison import compare
from projects.ukrmol_co.verify_atomic_references import verify as atomic_verify
from projects.ukrmol_co.verify_stretched_properties import files_in, require, sha

SOURCE = "e5f8e74e3323cba2986ac407f17d6691857f8a7e"
SOURCE_ARCHIVE = "80e7dec680e35faf822c8a548d2cacecf0d84042f681da25ce374e369043ec77"
ARCHIVE = "3a2521b4b8d12b6c4a816cfb51e15917f5dafbfce1777ddbeb5a1626cfe5498d"
CONTRACT = "projects/ukrmol_co/atomic-asymptote-5z-source-compat-release-20261103.json"
STOP_SOURCE = "875e3785496df4f81eda7cd61bd7ccf70cf85569"
STOP_ARCHIVE = "69a94cd4509aa530dac03a89e50399cdc9497ba0e2dc6bd3533f2bc391763085"


def calculation_source(raw):
    text = raw.decode()
    node = next(
        n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef) and n.name == "calculate"
    )
    return ("\n".join(text.splitlines()[node.lineno - 1 : node.end_lineno]) + "\n").encode()


def verify(path, collection, qz_archive, qz_collection):
    checked = atomic_verify(
        path,
        collection,
        source_commit=SOURCE,
        source_archive_sha256=SOURCE_ARCHIVE,
        archive_sha256=ARCHIVE,
        contract_path=CONTRACT,
        basis="aug-cc-pV5Z",
        source_entries=21,
        raw_entries=48,
        runtime_controls=8,
    )
    kept = files_in(path.read_bytes())
    source = files_in(kept["source.tar.gz"])
    c = json.loads(source[CONTRACT])
    preflight = json.loads(kept["actual-image-controls.json"])
    require(
        preflight["original_exit"] == 0 and "Ran 8 tests" in preflight["stderr"],
        "Actual-image no-quantum controls before science",
    )
    old = subprocess.check_output(
        [
            "git",
            "show",
            "493856bd83751e36d66596157e8cbfc9a74efb5a:projects/ukrmol_co/atomic_references.py",
        ]
    )
    require(
        calculation_source(old)
        == calculation_source(source["projects/ukrmol_co/atomic_references.py"])
        and sha(calculation_source(old)) == c["qz_calculation_source_sha256"],
        "Exact calculation-source bytes preserved,independent of AST metadata version",
    )
    qz = atomic_verify(qz_archive, qz_collection)
    predecessor = json.loads(source[c["predecessor_result"]])
    for name in ["total_energy_hartree", "fragment_sum_hartree", "electron_affinity_ev"]:
        require(qz[name] == predecessor[name], "Replayed qualified QZ predecessor " + name)
    comparison = compare(
        qz,
        checked,
        json.loads(source["projects/ukrmol_co/neutral-results.json"]),
        c["diagnostic_budgets"],
    )
    return checked | dict(
        basis="aug-cc-pV5Z",
        qz_predecessor_independently_replayed=True,
        calculation_source_byte_identical_to_qz=True,
        actual_image_preflight_exit=0,
        finite_basis_comparison=comparison,
        original_startup_stop_requalified=False,
        preserved_startup_stop_archive_sha256=STOP_ARCHIVE,
    )


def verify_startup_stop(path, collection):
    raw = path.read_bytes()
    require(sha(raw) == STOP_ARCHIVE == collection["archive"]["sha256"], "Original stopped archive")
    kept = files_in(raw)
    require(set(kept) == set(collection["files"]) and len(kept) == 14, "Stopped raw inventory")
    for name, pin in collection["files"].items():
        require(
            len(kept[name]) == pin["bytes"] and sha(kept[name]) == pin["sha256"],
            "Stopped original bytes " + name,
        )
    source = files_in(kept["source.tar.gz"])
    pins = json.loads(kept["source-manifest.json"])
    require(set(source) == set(pins) and len(source) == 19, "Stopped source inventory")
    for name, h in pins.items():
        require(
            sha(source[name]) == h
            and source[name] == subprocess.check_output(["git", "show", STOP_SOURCE + ":" + name]),
            "Stopped executed Git source " + name,
        )
    owner = json.loads(kept["execution.json"])
    container = json.loads(kept["container-inspect.json"])
    exits = [
        json.loads(kept["outputs/controls.json"])["original_exit"],
        owner["original_command_exit"],
        owner["original_container_exit"],
        owner["original_owner_exit"],
        json.loads(kept["service-return.json"])["original_service_exit"],
        json.loads(kept["launcher-exit.json"])["original_exit"],
    ]
    require(
        exits == [1, 2, 2, 2, 2, 2]
        and container["State"]["ExitCode"] == 2
        and not container["State"]["OOMKilled"]
        and not owner["timed_out"]
        and owner["source_integrity_after"],
        "Original startup exits/time/OOM/integrity",
    )
    log = kept["outputs/controls.log"].decode()
    require(
        "Ran 8 tests" in log
        and "FAILED (failures=1)" in log
        and "test_calculation_algorithm_is_identical_to_the_executed_qz_implementation" in log
        and not any(
            n.startswith("outputs/") and n not in ["outputs/controls.log", "outputs/controls.json"]
            for n in kept
        ),
        "Control failure before quantum outputs",
    )
    return dict(
        original_startup_stop_reconstructed=True,
        source_commit=STOP_SOURCE,
        archived_files_verified=14,
        executed_git_sources_verified=19,
        original_exits=exits,
        original_scientific_exit=None,
        new_integral_calls=0,
        new_scf_calls=0,
        new_cc_calls=0,
        original_stop_requalified=False,
        finite_method_endpoint_reconstruction_verdict=None,
        owner_wall_seconds=owner["wall_seconds"],
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive", type=Path)
    p.add_argument("collection", type=Path)
    p.add_argument("qz_archive", type=Path, nargs="?")
    p.add_argument("qz_collection", type=Path, nargs="?")
    p.add_argument("--startup-stop", action="store_true")
    a = p.parse_args()
    collection = json.loads(a.collection.read_text())
    if a.startup_stop:
        result = verify_startup_stop(a.archive, collection)
    elif a.qz_archive is None or a.qz_collection is None:
        p.error("Qualified original QZ archive/collection required")
    else:
        result = verify(
            a.archive, collection, a.qz_archive, json.loads(a.qz_collection.read_text())
        )
    print(json.dumps(result, indent=2, allow_nan=False))
