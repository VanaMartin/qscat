"""No-native controls for the separate aug-TZ64 property release and old routing."""

import json
import unittest
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from projects.ukrmol_co import stretched_basis_properties as properties
from projects.ukrmol_co import test_stretched_basis_properties as original
from projects.ukrmol_co.stretched_basis_properties_owner import container_command, service_command

CONTRACT = "projects/ukrmol_co/stretched-augtz-properties-contract-20261102.json"


def test_aug_owner_has_only_the_exact_reviewed_case_and_new_names():
    c = json.loads(Path(CONTRACT).read_text())
    cmd = container_command(
        Path("/work"), Path("/source"), Path("/out"), "augtz", c, contract_path=CONTRACT
    )
    assert "--cpuset-cpus=4-11" in cmd and "--memory=8g" in cmd and "--memory-swap=8g" in cmd
    assert "--network=none" in cmd and cmd[3] == c["container_name"]
    assert CONTRACT in cmd[-1] and c["runtime_controls_module"] in cmd[-1]
    assert "'--case', 'augtz'" in cmd[-1] and "dz-receipt" not in cmd[-1]
    compile(cmd[-1], "<augtz-property-owner>", "exec")
    service = service_command(
        Path("/capsule"), c["resources"], ["--case", "augtz", "--contract", CONTRACT]
    )
    assert "--unit=co-stretched-augtz-properties-20261102" in service
    assert "--property=RuntimeMaxSec=46800" in service


def fixture(root):
    target = root / "native"
    target.mkdir()
    (target / "config.json").write_text(
        json.dumps(
            dict(
                basis="aug-cc-pVTZ",
                bond_length=2.5,
                active_orbitals=[6, 3, 3, 0],
                sa_singlet_roots=[5] * 4,
                sa_triplet_roots=[5] * 4,
            )
        )
    )
    proof = dict(
        reuse_feasibility_verdict=True,
        corrected_reader_exit=0,
        fresh_public_payloads_verified=13,
        native_root_count=64,
        original_import_analysis_verifier_exits=[0] * 3,
        imported_checkpoint_sha256="checkpoint",
        original_property_counts=[5, 3, 6, 7, 6, 7, 3, 3],
    )
    c = dict(
        execution_authorized=True,
        source_sha256={},
        verification="proof.json",
        models=dict(
            augtz=dict(
                basis="aug-cc-pVTZ",
                target_directory="native",
                original_property_counts=proof["original_property_counts"],
                input_sha256={"output/CO/geom1/co.casscf.chk": "checkpoint"},
            )
        ),
    )
    return proof, c


def test_aug_requires_its_own_public_reuse_proof_and_explicit_release():
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        proof, c = fixture(root)
        for changes in [
            {},
            dict(reuse_feasibility_verdict=False),
            dict(corrected_reader_exit=1),
            dict(fresh_public_payloads_verified=0),
            dict(native_root_count=40),
            dict(imported_checkpoint_sha256="wrong"),
            dict(original_import_analysis_verifier_exits=[0, 2, 0]),
            dict(original_property_counts=[8] * 8),
        ]:
            (root / "proof.json").write_text(json.dumps(proof | changes))
            with patch.object(properties, "digest", return_value="checkpoint"):
                if changes:
                    with unittest.TestCase().assertRaises(ValueError):
                        properties.authenticate(root, root, "augtz", c)
                else:
                    assert properties.authenticate(root, root, "augtz", c) == (
                        root / "native",
                        None,
                    )
        (root / "proof.json").write_text(json.dumps(proof))
        bad = deepcopy(c)
        bad["execution_authorized"] = False
        with unittest.TestCase().assertRaises(ValueError):
            properties.authenticate(root, root, "augtz", bad)
        with unittest.TestCase().assertRaises(ValueError):
            properties.authenticate(root, root, "augtz", c, "receipt", "sha")


def test_aug_deck_changes_only_counts_and_preserves_the_native_set_order():
    counts = [5, 3, 6, 7, 6, 7, 3, 3]
    old = original.deck("dz").replace("6,5,5,6,5,6,3,4", ",".join(map(str, counts)))
    new = properties.matched_deck(old, counts)
    assert new.replace("8,8,8,8,8,8,8,8, ", ",".join(map(str, counts)) + ", ") == old
    with unittest.TestCase().assertRaises(ValueError):
        properties.matched_deck(old, [8] * 8)


def test_exact_frozen_source_pins_and_no_automatic_scientific_successor():
    c = json.loads(Path(CONTRACT).read_text())
    assert len(c["models"]["augtz"]["input_sha256"]) == 18
    assert c["native_denprop_calls_max"] == 1 and c["automatic_successor"] is False
    assert (
        c["qc_calls"] == c["target_diagonalizations"] == c["scattering_calls"] == c["new_maps"] == 0
    )
    assert c["route_holdout_s_matrices_analyzed"] == 0 and c["full_sweep_released"] is False
    assert len(c["source_sha256"]) >= 5
    for name, pin in c["source_sha256"].items():
        assert properties.digest(Path(name)) == pin, name


if __name__ == "__main__":
    functions = [
        fn for name, fn in list(globals().items()) if name.startswith("test_") and callable(fn)
    ]
    functions += [
        original.test_both_decks_change_only_the_exact_original_counts,
        original.test_64_inventory_rejects_wrong_count_bad_energy_and_missing_moments,
        original.test_tz_needs_independent_hash_pinned_dz_success_not_a_producer_claim,
        original.test_success_is_one_property_call_and_no_peer_release,
        original.test_failed_prerequisite_or_native_call_stops_without_retry,
        original.test_owner_is_one_named_readonly_no_network_eight_core_job,
    ]
    suite = unittest.TestSuite(unittest.FunctionTestCase(fn) for fn in functions)
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
