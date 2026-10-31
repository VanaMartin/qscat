"""Exact64-state preparation controls with no native/quantum calculation."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co import stretched_basis_preparation as prep
from projects.ukrmol_co.a2_channel_preflight import replace_assignment
from projects.ukrmol_co.compressed_cas12_preparation import dimensions, values
from projects.ukrmol_co.stretched_basis_preparation_owner import container_command, own
from projects.ukrmol_co.test_compressed_cas12_preparation import inventory, log

SOURCE = Path(__file__).resolve().parents[2]


def test_cards_retain_all64_native_roots_without_changing_congen_groups():
    congen, scatci = prep.cards(inventory(), SOURCE)
    original = (SOURCE / prep.TEMPLATES / "scattering.congen.doublet.B1.inp").read_text()
    assert congen.split("&wfngrp") == original.split("&wfngrp")
    assert values(scatci, "numtgt") == [8] * 8
    assert values(scatci, "ntgtf") == [i for i in range(1, 9) for _ in range(8)]
    assert values(scatci, "ntgts") == list(range(1, 9)) * 8
    assert values(congen, "nob0") == [8, 3, 3, 0]


def test_actual_contracted_dimension_includes_eight_roots_and_unchanged_l2():
    _, scatci = prep.cards(inventory(), SOURCE)
    d = dimensions(log(), scatci, roots_per_sector=8)
    assert d["uncontracted_l2_csfs"] == 84942 and d["contracted_continuum_dimension"] == 2256
    assert d["contracted_dimension"] == 87198 and d["retained_roots_per_sector"] == [8] * 8
    with unittest.TestCase().assertRaises(ValueError):
        dimensions(log(), scatci)
    with unittest.TestCase().assertRaises(ValueError):
        dimensions(
            log(), replace_assignment(scatci, "numtgt", "7,8,8,8,8,8,8,8"), roots_per_sector=8
        )


def test_bad_property_prerequisite_makes_no_native_call():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        with (
            patch.object(prep, "authenticate", side_effect=ValueError("no public64-state input")),
            patch("projects.ukrmol_co.a2_channel_preflight.stage") as stage,
        ):
            assert prep.compute(root, root, root, "dz", {}) == 2
        stage.assert_not_called()
        assert json.loads((root / "result.json").read_text())["integral_calls"] == 0


def test_tz_preparation_requires_separate_independent_dz_preparation_receipt():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        contract = root / prep.CONTRACT
        contract.parent.mkdir(parents=True)
        contract.write_text("{}")
        receipt = root / "receipt.json"
        value = dict(
            case="dz",
            preparation_verdict=True,
            contract_sha256=prep.digest(contract),
            original_scientific_command_container_owner_exits=[0] * 4,
            public_reader_exit=0,
            roots_per_sector=[8] * 8,
            moints_sha256="a" * 64,
            csf_sha256="b" * 64,
        )
        for changes in [{}, dict(public_reader_exit=1), dict(preparation_verdict=False)]:
            receipt.write_text(json.dumps(value | changes))
            want = prep.digest(receipt)
            if changes:
                with unittest.TestCase().assertRaises(ValueError):
                    prep.prior_preparation(root, root, "tz", "receipt.json", want)
            else:
                assert (
                    prep.prior_preparation(root, root, "tz", "receipt.json", want)["sha256"] == want
                )
        with unittest.TestCase().assertRaises(ValueError):
            prep.prior_preparation(root, root, "tz")
        assert prep.prior_preparation(root, root, "dz") is None


def test_failed_integral_never_starts_congen_or_retries():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        target = root / "target"
        geometry = target / "output/CO/geom1"
        geometry.mkdir(parents=True)
        work = root / "work"
        work.mkdir()
        for n in prep.COPY:
            (geometry / n).write_bytes(b"saved input")
        (geometry / "fort.24").write_text("saved properties")
        (target / "target.json").write_text("{}")
        native = {
            s: dict(energies=np.arange(8) - 100.0)
            for s in [
                "singlet.A1",
                "triplet.A1",
                "singlet.B1",
                "triplet.B1",
                "singlet.B2",
                "triplet.B2",
                "singlet.A2",
                "triplet.A2",
            ]
        }
        contract = dict(engine_sha256={}, resources=dict(integral_timeout_seconds=10))
        with (
            patch.object(prep, "authenticate", return_value=(target, geometry, None)),
            patch(
                "projects.ukrmol_co.stretched_basis_feasibility.read_cidata", return_value=native
            ),
            patch("projects.ukrmol_co.a2_channel_preflight.property_inventory", return_value={}),
            patch(
                "projects.ukrmol_co.a2_channel_preflight.stage",
                side_effect=ValueError("native exit1"),
            ) as stage,
        ):
            assert prep.compute(root, SOURCE, work, "dz", contract) == 2
        stage.assert_called_once()
        report = json.loads((work / "result.json").read_text())
        assert (
            report["integral_calls"] == 1
            and report["congen_calls"] == 0
            and report["preparation_verdict"] is None
        )


def test_owner_is_one_case_readonly_four_core_no_network_job():
    command = container_command(
        Path("/root"), Path("/source"), Path("/out"), "dz", dict(image="sha256:pinned")
    )
    assert (
        "--cpuset-cpus=0-3" in command and "--memory=32g" in command and "--network=none" in command
    )
    assert [command[i + 1] for i, v in enumerate(command) if v == "-v"] == [
        "/root:/work:ro",
        "/source:/opt/qmodeling:ro",
        "/out:/out",
    ]
    assert "tz" not in command[-1] and "inner" not in command[-1]
    compile(command[-1], "<preparation-owner>", "exec")


def test_owner_failed_authentication_does_not_start_scientific_container():
    with tempfile.TemporaryDirectory() as tmp:
        capsule = Path(tmp) / "prepared/capsule"
        source = capsule / "source"
        path = source / prep.CONTRACT
        path.parent.mkdir(parents=True)
        path.write_text("{}")
        (capsule / "source.tar.gz").write_bytes(b"frozen")
        (capsule / "source-manifest.json").write_text("{}")
        with (
            patch(
                "projects.ukrmol_co.stretched_basis_preparation_owner.authenticate",
                side_effect=ValueError("stop"),
            ),
            patch("subprocess.Popen") as process,
        ):
            assert own(capsule, prep.digest(capsule / "source.tar.gz"), "commit", "dz") == 2
            process.assert_not_called()


def test_reuse_requires_frozen_tz_release_and_same_public_csf():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        path = root / prep.CONTRACT
        path.parent.mkdir(parents=True)
        path.write_text("{}")
        prepared = root / "prepared/dz"
        prepared.mkdir(parents=True)
        (prepared / "fort.721").write_bytes(b"algebraic CSFs")
        receipt = root / "receipt.json"
        receipt.write_text(json.dumps(dict(csf_sha256=prep.digest(prepared / "fort.721"))))
        previous = dict(path="receipt.json", sha256=prep.digest(receipt))
        release = dict(
            predecessor=previous,
            contract_sha256=prep.digest(path),
            preparation_directory="prepared/dz",
            input_sha256={"fort.721": prep.digest(prepared / "fort.721")},
        )
        (root / prep.REUSE).write_text(json.dumps(release))
        assert prep.reuse_prerequisites(root, root, "tz", previous)[0] == prepared
        with unittest.TestCase().assertRaises(ValueError):
            prep.reuse_prerequisites(root, root, "dz", previous)
        (prepared / "fort.721").write_bytes(b"changed")
        with unittest.TestCase().assertRaises(ValueError):
            prep.reuse_prerequisites(root, root, "tz", previous)


def test_reuse_copies_only_algebra_and_rejects_deck_or_target_changes():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        prepared = root / "dz"
        prepared.mkdir()
        for name, data in [("fort.721", b"CSFs"), ("congen.B1.out", b"completed log")]:
            (prepared / name).write_bytes(data)
        (prepared / "scattering.congen.doublet.B1.inp").write_text("exact deck")
        work = root / "tz"
        work.mkdir()
        (work / "fort.710").write_bytes(b"target CSFs")
        release = dict(
            preparation_directory="dz",
            native_inventory=inventory(),
            input_sha256={n: prep.digest(prepared / n) for n in ["fort.721", "congen.B1.out"]},
            target_csf_sha256={"fort.710": prep.digest(work / "fort.710")},
        )
        for deck in ["changed deck", "exact deck"]:
            if deck == "changed deck":
                with unittest.TestCase().assertRaises(ValueError):
                    prep.reuse_csf(work, prepared, release, inventory(), deck)
            else:
                assert (
                    prep.reuse_csf(work, prepared, release, inventory(), deck)[
                        "original_congen_exit"
                    ]
                    == 0
                )
        assert not (work / "moints").exists() and not (work / "fort.26").exists()
        (work / "fort.710").write_bytes(b"changed target")
        with unittest.TestCase().assertRaises(ValueError):
            prep.reuse_csf(work, prepared, release, inventory(), "exact deck")


def test_reuse_owner_passes_explicit_opt_in_without_automatic_successor():
    command = container_command(
        Path("/root"),
        Path("/source"),
        Path("/out"),
        "tz",
        dict(image="sha256:pinned"),
        "receipt",
        "a" * 64,
        True,
    )
    assert "--reuse-dz-csf" in command[-1] and "--previous-receipt" in command[-1]
    assert "inner" not in command[-1]
    compile(command[-1], "<reuse-owner>", "exec")


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
