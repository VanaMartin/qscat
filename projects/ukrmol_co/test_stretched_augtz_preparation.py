"""No-native release controls for one own-integral/conditional-CSF reuse step."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from projects.ukrmol_co import stretched_basis_preparation as prep
from projects.ukrmol_co import test_stretched_basis_preparation as old
from projects.ukrmol_co.stretched_basis_preparation_owner import container_command

CONTRACT = "projects/ukrmol_co/stretched-augtz-preparation-contract-20261103.json"


def test_one_explicit_aug_case_new_name_and_readonly_four_core_allocation():
    c = json.loads(Path(CONTRACT).read_text())
    command = container_command(
        Path("/root"),
        Path("/source"),
        Path("/out"),
        "augtz",
        c,
        "receipt",
        "a" * 64,
        True,
        CONTRACT,
    )
    assert command[3] == "co-stretched-augtz-preparation-20261103"
    assert (
        "--cpuset-cpus=0-3" in command and "--memory=32g" in command and "--network=none" in command
    )
    assert CONTRACT in command[-1] and "--reuse-dz-csf" in command[-1]
    assert c["runtime_controls_module"] in command[-1] and "inner" not in command[-1]
    compile(command[-1], "<aug-preparation-owner>", "exec")


def test_aug_cannot_fall_back_to_a_fresh_congen_even_with_no_other_inputs():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        with patch("projects.ukrmol_co.a2_channel_preflight.stage") as stage:
            assert prep.compute(root, root, root, "augtz", {}) == 2
        stage.assert_not_called()
        r = json.loads((root / "result.json").read_text())
        assert r["integral_calls"] == r["congen_calls"] == 0
        assert "no fresh/fallback CONGEN" in r["error"]


def test_aug_authentication_uses_its_own_contract_and_public_property_receipt():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        props = root / "properties"
        props.mkdir()
        (props / "result.json").write_text(
            json.dumps(dict(scientific_exit=0, property_inventory_verdict=True))
        )
        c = dict(
            execution_authorized=True,
            source_sha256={},
            property_input_contract="aug-contract.json",
            reuse_release="reuse.json",
            reuse_source_contract_sha256="old-dz-contract",
            models=dict(
                augtz=dict(
                    property_receipt="receipt.json",
                    property_directory="properties",
                    input_sha256={},
                )
            ),
        )
        (root / "aug-contract.json").write_text(
            json.dumps(dict(models=dict(augtz=dict(target_directory="target", input_sha256={}))))
        )
        receipt = dict(
            case="augtz",
            contract_sha256="pin",
            independent_property_inventory_verdict=True,
            public_reader_exit=0,
            target_count=64,
            roots_per_sector=[8] * 8,
            original_scientific_command_container_owner_exits=[0] * 4,
            property_sha256="pin",
        )
        for changes in [{}, dict(public_reader_exit=1), dict(case="tz"), dict(target_count=40)]:
            (root / "receipt.json").write_text(json.dumps(receipt | changes))
            with (
                patch.object(prep, "digest", return_value="pin"),
                patch.object(prep, "prior_preparation", return_value="qualified-dz"),
            ):
                if changes:
                    with unittest.TestCase().assertRaises(ValueError):
                        prep.authenticate(root, root, "augtz", c)
                else:
                    assert prep.authenticate(root, root, "augtz", c) == (
                        root / "target",
                        props,
                        "qualified-dz",
                    )


def test_aug_contract_source_and_reuse_overlay_are_exact():
    c = json.loads(Path(CONTRACT).read_text())
    release = json.loads(Path(c["reuse_release"]).read_text())
    assert len(c["source_sha256"]) >= 8
    for n, pin in c["source_sha256"].items():
        assert prep.digest(Path(n)) == pin, n
    assert release["case"] == "augtz" and release["preparation_contract"] == CONTRACT
    assert release["contract_sha256"] == prep.digest(Path(CONTRACT))
    assert c["native_integral_calls_max"] == 1 and c["native_congen_calls_max"] == 0
    assert c["automatic_successor"] is c["full_sweep_released"] is False


if __name__ == "__main__":
    functions = [
        fn for name, fn in list(globals().items()) if name.startswith("test_") and callable(fn)
    ]
    functions += [
        old.test_cards_retain_all64_native_roots_without_changing_congen_groups,
        old.test_actual_contracted_dimension_includes_eight_roots_and_unchanged_l2,
        old.test_bad_property_prerequisite_makes_no_native_call,
        old.test_tz_preparation_requires_separate_independent_dz_preparation_receipt,
        old.test_failed_integral_never_starts_congen_or_retries,
        old.test_owner_is_one_case_readonly_four_core_no_network_job,
        old.test_owner_failed_authentication_does_not_start_scientific_container,
        old.test_reuse_requires_frozen_tz_release_and_same_public_csf,
        old.test_reuse_copies_only_algebra_and_rejects_deck_or_target_changes,
        old.test_reuse_owner_passes_explicit_opt_in_without_automatic_successor,
    ]
    suite = unittest.TestSuite(unittest.FunctionTestCase(fn) for fn in functions)
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
