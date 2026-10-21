"""Finish only the seven unrun CSF sectors after a saved CONGEN trailer-check defect."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path

from projects.ukrmol_co.a2_channel_preflight import stage
from projects.ukrmol_co.fifty_import_verification import digest
from projects.ukrmol_co.immutable_continuum import (
    IRREPS,
    congen_completion,
    congen_deck,
    csf_records,
    native_inventory,
)


def complete(root, output, contract):
    parent = json.loads((root / contract["parent_contract"]).read_text())

    def verify():
        for name, want in (parent["input_sha256"] | contract["input_sha256"]).items():
            if digest(root / name) != want:
                raise ValueError(f"Changed completion input: {name}")

    verify()
    original = root / contract["original_preparation"]
    owner = json.loads((root / contract["original_execution"]).read_text())
    failed = json.loads((original / "failure.json").read_text())
    assert owner["original_owner_exit"] == 1
    assert (
        owner["steps"][0]["original_command_exit"]
        == owner["steps"][0]["original_container_exit"]
        == 1
    )
    assert [(s["name"], s["original_exit"], s["timed_out"]) for s in failed["stages"]] == [
        ("integrals", 0, False),
        ("congen.singlet.A1", 0, False),
    ]
    assert failed["error"] == "AssertionError()"
    assert 'assert "has completed execution"' in (root / contract["original_traceback"]).read_text()
    inventory = native_inventory((original / "integrals.native.out").read_text())
    dimensions = contract["target_csf_dimensions"]
    congen_completion(
        (original / "congen.singlet.A1.out").read_text(), 710, dimensions["singlet.A1"]
    )
    saved = csf_records(
        original / "fort.710", inventory["total"], "singlet.A1", dimensions["singlet.A1"]
    )
    output.mkdir(exist_ok=False)
    for name in (
        "moints",
        "co.molden",
        "fort.710",
        "integrals.native.out",
        "congen.singlet.A1.out",
        "target.congen.singlet.A1.inp",
        "integrals.inp",
    ):
        shutil.copy2(original / name, output / name)
    (output / "fort.16").symlink_to("moints")
    report = dict(
        stages=[],
        completed_layout_generation_verdict=None,
        original_preparation_exits=[1, 1, 1],
        original_native_stage_exits=[0, 0],
        retained_subspace_dominance_verdict=False,
        native_target_vector_reuse_verdict=None,
        scattering_ensemble_verdict=None,
        new_integral_calls=0,
        new_orbital_optimizations=0,
        target_diagonalizations=0,
        property_exports=0,
        scattering_inner_solves=0,
        native_inventory=inventory,
        saved_first_sector_reconstructed=True,
        target_csf_records={"singlet.A1": saved},
    )
    try:
        for irrep in IRREPS:
            for spin in ("singlet", "triplet"):
                label = f"{spin}.{irrep}"
                if label == "singlet.A1":
                    continue
                before = root / contract["baseline_inputs"] / f"target.congen.{label}.inp"
                card = output / f"target.congen.{label}.inp"
                card.write_text(congen_deck(before.read_text(), inventory["total"]))
                stage(
                    output,
                    report,
                    f"congen.{label}",
                    ["/opt/ukrmolp/bin.double/congen"],
                    card,
                    15,
                    os.environ.copy(),
                )
                unit = (710 if spin == "singlet" else 730) + IRREPS.index(irrep)
                congen_completion(
                    (output / f"congen.{label}.out").read_text(), unit, dimensions[label]
                )
                report["target_csf_records"][label] = csf_records(
                    output / f"fort.{unit}", inventory["total"], label, dimensions[label]
                )
        for name in ("moints", "co.molden", "fort.710", "integrals.inp"):
            assert digest(output / name) == digest(original / name)
        report["completed_layout_generation_verdict"] = True
        report["scope"] = (
            "Generation-only completion: seven previously unrun target-CSF calls. "
            "Original failure preserved; no integral rerun, target-vector conversion "
            "or Hamiltonian/property/scattering claim."
        )
        verify()
        (output / "result.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    except BaseException as error:
        report["error"] = repr(error)
        (output / "failure.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        raise
    finally:
        verify()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    args = parser.parse_args()
    complete(args.root, args.output, json.loads(args.contract.read_text()))


if __name__ == "__main__":
    main()
