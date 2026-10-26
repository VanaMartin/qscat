"""Preparation routing/card/dimension controls; never calls a native program."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from projects.ukrmol_co import compressed_cas12_preparation as prep
from projects.ukrmol_co.compressed_cas12_preparation_owner import container_command
from projects.ukrmol_co.immutable_continuum import native_inventory

TEMPLATES = Path(__file__).parent / "compressed-cas12-preparation-inputs-20261026"


def inventory():
    return {
        "target": prep.TARGET,
        "PCOs": [0] * 8,
        "TGT+PCOs": prep.TARGET,
        "continuum": [51, 34, 34, 22, 0, 0, 0, 0],
        "total": [59, 37, 37, 22, 0, 0, 0, 0],
    }


def log():
    return (
        "\n".join(
            f"WFN GROUP {i + 1}\n{s} {irrep} x continuum^1"
            for i, (irrep, s) in enumerate(
                (irrep, s) for irrep in ["A1", "B1", "B2", "A2"] for s in ["singlet", "triplet"]
            )
        )
        + """
WFN GROUP 9
L^2: ORMAS distribution 1/1
CSF NUMBERS 906049 TO 990990 ( 84942 CSFS)
TOTAL NUMBER OF CSF'S GENERATED IS 990990
"""
    )


def cards(sector="B1", inv=None):
    return prep.resolved_cards(
        (TEMPLATES / f"scattering.congen.doublet.{sector}.inp").read_text(),
        (TEMPLATES / f"scattering.scatci.doublet.{sector}.inp").read_text(),
        inv or inventory(),
        sector,
    )


def test_cas12_inventory_is_explicit_and_cas11_default_is_not_relaxed():
    text = "\n".join(
        f"Final number of orbitals for {kind}\n" + " ".join(map(str, row))
        for kind, row in inventory().items()
    )
    text += "\nTwo_electron_integrals took [s]: 1\n" + "NUMBER OF ERROR MESSAGES: 0\n" * 4
    assert native_inventory(text, target_counts=prep.TARGET) == inventory()
    with unittest.TestCase().assertRaises(ValueError):
        native_inventory(text)
    with unittest.TestCase().assertRaises(ValueError):
        native_inventory(text.replace("59 37 37 22", "58 37 37 22"), target_counts=prep.TARGET)


def test_resolved_cards_preserve_geometry_free_algebra_and_own_target_order():
    for sector in ["B1", "B2"]:
        congen, scatci = cards(sector)
        assert prep.values(congen, "nob0") == [8, 3, 3, 0]
        assert prep.values(scatci, "ntgtf") == [i for i in range(1, 9) for _ in range(5)]
        assert prep.values(scatci, "ntgts") == list(range(1, 6)) * 8
        # Only native address list whitespace may change at the planning anchor.
        original = (TEMPLATES / f"scattering.congen.doublet.{sector}.inp").read_text()
        assert congen.split("&wfngrp") == original.split("&wfngrp")
        inv = inventory()
        inv["continuum"] = [52, 35, 35, 23, 0, 0, 0, 0]
        inv["total"] = [60, 38, 38, 23, 0, 0, 0, 0]
        updated, changed = cards(sector, inv)
        assert prep.values(updated, "nob") == [60, 38, 38, 23]
        assert sum(prep.values(changed, "notgt")) == 290


def test_contracted_dimension_comes_from_actual_continuum_and_final_l2_ranges():
    _, scatci = cards()
    result = prep.dimensions(log(), scatci)
    assert result["raw_congen_csfs"] == 990990
    assert result["uncontracted_l2_csfs"] == 84942
    assert result["contracted_continuum_dimension"] == 1410
    assert result["contracted_dimension"] == 86352
    for text in [
        log().replace("84942 CSFS", "84941 CSFS"),
        log().replace("WFN GROUP 9", "no ninth group"),
        log().replace("triplet A2", "triplet B1"),
        log() + "TOTAL NUMBER OF CSF'S GENERATED IS 990990\n",
    ]:
        with unittest.TestCase().assertRaises(ValueError):
            prep.dimensions(text, scatci)


def test_failed_property_prerequisite_starts_no_native_stage():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        with (
            patch.object(prep, "authenticate", side_effect=ValueError("unverified properties")),
            patch.object(prep, "stage") as stage,
        ):
            exit_code = prep.compute(root, root, root, {})
        stage.assert_not_called()
        report = json.loads((root / "result.json").read_text())
        assert exit_code == 2 and report["preparation_verdict"] is None
        assert report["integral_calls"] == report["congen_calls"] == 0


def test_failed_integral_call_never_starts_congen_or_retries():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        work = root / "new"
        work.mkdir()
        target = root / "target"
        geometry = target / "output/CO/geom1"
        geometry.mkdir(parents=True)
        for name in ["co.molden", "co.casscf.chk", "fort.26"] + [
            f"fort.{n}" for n in [710, 730, 711, 731, 712, 732, 713, 733]
        ]:
            (geometry / name).write_bytes(b"immutable saved target input")
        (geometry / "fort.24").write_text("saved matched property output")
        (target / "target.json").write_text("{}")
        contract = dict(
            templates=str(TEMPLATES.resolve()),
            resources=dict(integral_timeout_seconds=10),
            input_sha256={},
            target_directory="target",
            property_directory="target/output/CO/geom1",
        )
        with (
            patch.object(prep, "authenticate", return_value=(target, geometry, {})),
            patch.object(prep, "stage", side_effect=ValueError("native integral exit1")) as stage,
        ):
            exit_code = prep.compute(root, root, work, contract)
        stage.assert_called_once()
        report = json.loads((work / "result.json").read_text())
        assert exit_code == 2 and report["integral_calls"] == 1 and report["congen_calls"] == 0
        assert report["preparation_verdict"] is None and not report["sectors"]


def test_exact_owner_command_has_only_one_writable_mount_and_four_cores():
    command = container_command(
        Path("/persistent"), Path("/frozen"), Path("/new"), dict(image="sha256:pinned")
    )
    assert "--cpuset-cpus=0-3" in command and "--memory=16g" in command
    mounts = [command[i + 1] for i, arg in enumerate(command) if arg == "-v"]
    assert mounts == ["/persistent:/work:ro", "/frozen:/opt/qmodeling:ro", "/new:/out"]
    compile(command[-1], "<container-controls-and-preparation>", "exec")
    assert "test_compressed_cas12_preparation" in command[-1]
    assert (
        "sparse_pilot" not in command[-1]
        and "compressed_cas12_properties','--root" not in command[-1]
    )


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
