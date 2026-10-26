"""Freeze controls:64-state inventory, strict predecessor and no automatic peer job."""

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co import a2_channel_preflight as native
from projects.ukrmol_co import stretched_basis_feasibility as reader
from projects.ukrmol_co import stretched_basis_properties as properties
from projects.ukrmol_co.stretched_basis_properties_owner import (
    container_command,
    own,
    service_command,
)

COUNTS = {"dz": [6, 5, 5, 6, 5, 6, 3, 4], "tz": [6, 4, 5, 6, 5, 7, 3, 4]}


def deck(case):
    return (
        "&input\n ntgt=8,\n nftsor=710,730,711,731,712,732,713,733,\n"
        " ntgtf=1,2,3,4,5,6,7,8,\n ntgtl="
        + ",".join(map(str, COUNTS[case]))
        + ", ! unchanged comment\n nftint=16, nftg=26, nftmt=24,\n/\n"
    )


def table():
    expected = {
        sector: [-100 + 0.1 * (index * 8 + j) for j in range(8)]
        for index, sector in enumerate(reader.SECTORS)
    }
    rows = ["6 0 0 0 64"]
    for index in range(64):
        sector = reader.SECTORS[index // 8]
        spin, irrep = sector.split(".")
        rows.append(f"5 {index + 1} 0 0 0 0 0 0 {expected[sector][index % 8]:.17g} {spin} {irrep}")
        for rank in [1, 2]:
            rows.append(f"1 {index + 1} 0 {index + 1} 0 3 {rank} 0 0.0")
    return "\n".join(rows), expected


def test_both_decks_change_only_the_exact_original_counts():
    for case in ["dz", "tz"]:
        original = deck(case)
        changed = properties.matched_deck(original, COUNTS[case])
        assert (
            changed.replace("8,8,8,8,8,8,8,8, ", ",".join(map(str, COUNTS[case])) + ", ")
            == original
        )
        for corrupt in [
            original.replace("nftg=26", "nftg=27"),
            original.replace("ntgt=8", "ntgt=7"),
            original.replace("ntgtf=1,2", "ntgtf=2,1"),
        ]:
            with unittest.TestCase().assertRaises(ValueError):
                properties.matched_deck(corrupt, COUNTS[case])


def test_64_inventory_rejects_wrong_count_bad_energy_and_missing_moments():
    text, expected = table()
    assert native.property_inventory(text, expected)["target_count"] == 64
    for broken in [
        text.replace("6 0 0 0 64", "6 0 0 0 40"),
        text.replace("-100 ", "-99 ", 1),
        "\n".join(line for line in text.splitlines() if not line.startswith("1 64 ")),
    ]:
        with unittest.TestCase().assertRaises(ValueError):
            native.property_inventory(broken, expected)


def receipt_value():
    return dict(
        case="dz",
        contract_sha256="frozen",
        independent_property_inventory_verdict=True,
        target_count=64,
        roots_per_sector=[8] * 8,
        original_scientific_command_container_owner_exits=[0] * 4,
        public_reader_exit=0,
        property_sha256="a" * 64,
    )


def test_tz_needs_independent_hash_pinned_dz_success_not_a_producer_claim():
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        path = root / "receipt.json"
        for changes in [
            {},
            dict(public_reader_exit=1),
            dict(independent_property_inventory_verdict=False),
            dict(original_scientific_command_container_owner_exits=[2] * 4),
            dict(target_count=40),
        ]:
            path.write_text(json.dumps(receipt_value() | changes))
            want = properties.digest(path)
            if changes:
                with unittest.TestCase().assertRaises(ValueError):
                    properties.predecessor(root, "tz", "frozen", "receipt.json", want)
            else:
                assert (
                    properties.predecessor(root, "tz", "frozen", "receipt.json", want)["case"]
                    == "dz"
                )
        for args in [(None, None), ("receipt.json", "bad")]:
            with unittest.TestCase().assertRaises(ValueError):
                properties.predecessor(root, "tz", "frozen", *args)
        assert properties.predecessor(root, "dz", "frozen") is None


def simulated(fail):
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        work = root / "work"
        work.mkdir()
        geometry = root / "output/CO/geom1"
        (geometry / "inputs").mkdir(parents=True)
        for name in properties.COPY:
            (geometry / name).write_bytes(b"saved-input")
        (geometry / "inputs/target.denprop.inp").write_text(deck("dz"))
        (root / "target.json").write_text(json.dumps(dict(ground_dipole_au=[0, 0, 0])))
        text, expected = table()
        actual_digest = properties.digest

        def measured(path):
            return (
                "pinned-executable"
                if str(path).startswith("/opt/ukrmolp/bin")
                else actual_digest(path)
            )

        def call(*args):
            if fail == "native":
                raise ValueError("original native exit failure")
            (work / "target.denprop.out").write_text("has completed execution")
            (work / "fort.24").write_text(text)

        auth = ValueError("prerequisite rejected") if fail == "prerequisite" else None
        contract = dict(
            models=dict(dz=dict(original_property_counts=COUNTS["dz"])),
            denprop_sha256="pinned-executable",
            resources=dict(native_timeout_seconds=64800),
            matched_deck_sha256=hashlib.sha256(
                properties.matched_deck(deck("dz"), COUNTS["dz"]).encode()
            ).hexdigest(),
        )
        if fail == "deck":
            contract["matched_deck_sha256"] = "wrong-frozen-deck"
        with (
            patch.object(properties, "authenticate", return_value=(root, None), side_effect=auth),
            patch.object(properties, "digest", side_effect=measured),
            patch.object(
                reader,
                "read_cidata",
                return_value={k: dict(energies=np.array(v)) for k, v in expected.items()},
            ),
            patch.object(native, "stage", side_effect=call) as stage,
            patch("subprocess.Popen") as process,
        ):
            code = properties.compute(root, root, work, "dz", contract)
            process.assert_not_called()
        return code, json.loads((work / "result.json").read_text()), stage.call_count


def test_success_is_one_property_call_and_no_peer_release():
    code, report, calls = simulated(None)
    assert code == 0 and calls == report["native_denprop_calls"] == 1
    assert report["property_inventory"]["target_count"] == 64
    assert report["qc_calls"] == report["target_diagonalizations"] == report["new_maps"] == 0


def test_failed_prerequisite_or_native_call_stops_without_retry():
    for failure, expected in [("prerequisite", 0), ("deck", 0), ("native", 1)]:
        code, report, calls = simulated(failure)
        assert code == 2 and calls == expected and report["property_inventory_verdict"] is None


def test_owner_is_one_named_readonly_no_network_eight_core_job():
    command = container_command(
        Path("/old"), Path("/source"), Path("/new"), "dz", dict(image="sha256:pinned")
    )
    assert (
        "--cpuset-cpus=4-11" in command
        and "--memory=64g" in command
        and "--network=none" in command
    )
    assert [command[i + 1] for i, v in enumerate(command) if v == "-v"] == [
        "/old:/work:ro",
        "/source:/opt/qmodeling:ro",
        "/new:/out",
    ]
    assert "tz" not in command[-1]
    compile(command[-1], "<container-owner>", "exec")


def test_host_predecessor_failure_never_creates_a_scientific_container():
    with tempfile.TemporaryDirectory() as temporary:
        capsule = Path(temporary) / "prepared/capsule"
        source = capsule / "source"
        path = source / properties.CONTRACT
        path.parent.mkdir(parents=True)
        path.write_text(json.dumps(dict(resources={})))
        (capsule / "source.tar.gz").write_bytes(b"frozen")
        (capsule / "source-manifest.json").write_text("{}")
        with (
            patch(
                "projects.ukrmol_co.stretched_basis_properties_owner.authenticate",
                side_effect=ValueError("No independent DZ receipt"),
            ),
            patch("subprocess.Popen") as process,
        ):
            assert own(capsule, properties.digest(capsule / "source.tar.gz"), "commit", "tz") == 2
            process.assert_not_called()
        assert (
            json.loads((capsule / "execution.json").read_text())["peer_job_automatically_released"]
            is False
        )


def test_one_durable_service_has_the_frozen_total_owner_cap():
    command = service_command(Path("/capsule"), dict(owner_timeout_seconds=68400), ["--case", "dz"])
    assert "--property=RuntimeMaxSec=68400" in command and "--wait" in command
    assert "--under-service" in command and not any("tz" in word for word in command)


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
