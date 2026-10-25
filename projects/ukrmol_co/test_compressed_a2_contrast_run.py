"""Owner routing controls: failed guards/identities never spend on another sector."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co import compressed_a2_contrast_run as owner
from projects.ukrmol_co.pole_character import screen_witnesses
from projects.ukrmol_co.test_compressed_a2_pole_contrast import BASE, CHARACTER, PROPOSAL


def screen(sector):
    residue = np.zeros((6, 6), complex)
    residue[0, 0] = 0.1j
    return screen_witnesses(
        [residue] * 4, CHARACTER["open_labels"][sector], 0.08, CHARACTER["character"]
    )


def staged(outcomes, *, pair_mismatch=False, corrupt_after=False):
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        seen = []

        def sector(*args, **kwargs):
            name, record = args[7], args[9]
            seen.append(name)
            if name == "B2":
                assert kwargs["seed"] == 0.14 - 0.04j
            record.update(
                disposition=outcomes[name],
                completed=True,
                new_maps_reserved=42,
                pole_incident_hartree=[
                    0.14 + (0.01 if pair_mismatch and name == "B2" else 0),
                    -0.04,
                ],
                screen=screen(name),
            )
            return record

        with (
            patch.object(
                owner,
                "load_prerequisites",
                return_value=(
                    PROPOSAL,
                    BASE,
                    CHARACTER,
                    CHARACTER["open_labels"],
                    {"B1": {}, "B2": {}},
                ),
            ),
            patch.object(
                owner,
                "verify_inputs",
                side_effect=[123, ValueError("changed")] if corrupt_after else [123, 123],
            ),
            patch.object(owner, "run_sector", side_effect=sector),
        ):
            owner.run(root, root, root / "output", {})
        return seen, json.loads((root / "output/result.json").read_text())


def test_resolved_b1_rejection_never_opens_b2():
    seen, report = staged({"B1": "reject"})
    assert seen == ["B1"] and report["scientific_exit"] == 1
    assert report["model_insensitivity_verdict"] is False and report["new_maps_reserved"] == 42
    assert "five_pair" not in report


def test_inconclusive_b1_never_opens_b2_or_becomes_physical_rejection():
    seen, report = staged({"B1": "inconclusive"})
    assert seen == ["B1"] and report["scientific_exit"] == 2
    assert report["model_insensitivity_verdict"] is None


def test_qualified_pass_runs_b2_with_new_b1_seed_and_requires_pair():
    seen, report = staged({"B1": "pass", "B2": "pass"})
    assert seen == ["B1", "B2"] and report["scientific_exit"] == 0
    assert report["model_insensitivity_verdict"] is True
    assert report["five_pair"]["disposition"] == "pass"
    assert report["new_maps_reserved"] == 84


def test_failed_counterpart_and_after_integrity_override_a_model_rejection():
    _, report = staged({"B1": "pass", "B2": "reject"}, pair_mismatch=True)
    assert report["scientific_exit"] == 2 and report["model_insensitivity_verdict"] is None
    _, report = staged({"B1": "reject"}, corrupt_after=True)
    assert report["scientific_exit"] == 2 and report["model_insensitivity_verdict"] is None
    assert report["sectors"][0]["disposition"] == "reject"  # Original result is not erased.


def test_failed_real_guard_never_starts_complex_pole_process():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        record = dict(real_maps_reserved=0, timed_out=False)
        model = dict(interface=dict(channel_labels=CHARACTER["open_labels"]["B1"]))
        with (
            patch.object(owner, "load_pole_model", return_value=model),
            patch.object(owner, "real_guard", side_effect=ValueError("real guard failed")),
            patch.object(owner.subprocess, "Popen") as process,
        ):
            owner.run_sector(
                root,
                root / "B1",
                {},
                PROPOSAL,
                BASE,
                CHARACTER,
                CHARACTER["open_labels"]["B1"],
                "B1",
                {},
                record,
            )
        process.assert_not_called()
        assert record["disposition"] == "inconclusive"
        assert record["new_maps_reserved"] == 0


def test_measured_phase_failure_keeps_first_capture_and_skips_second_guard():
    with TemporaryDirectory() as temporary:
        directory = Path(temporary)
        model = dict(
            energies=np.ones(99) * 0.2,
            saved=np.zeros((99, 6, 6)),
            ground=-100,
            poles=np.ones(2),
            amplitudes=np.ones((2, 6)),
            thresholds=np.zeros(6),
            angular=np.ones(6),
            c2=np.eye(6),
            c3=np.eye(6),
        )
        invariants = dict(
            flux_defect=0, wronskian_residual=0, unitarity_defect=0, reciprocity_defect=0
        )
        progress = dict(real_maps_reserved=0)
        actual = np.eye(6) * np.exp(4e-7j)  # Each entry passes; phase sum exceeds1e-6.
        with (
            patch.object(owner, "surface_matrix", return_value=np.eye(6)),
            patch.object(owner, "outgoing_map", return_value=(np.eye(6), np.eye(6), {})) as maps,
            patch.object(owner, "physical_s", return_value=(actual, invariants)),
        ):
            with unittest.TestCase().assertRaisesRegex(ValueError, "real-axis"):
                owner.real_guard(model, directory, "compressed/five/B1", PROPOSAL, progress)
        maps.assert_called_once()
        data = json.loads((directory / "guard.json").read_text())
        assert data["scientific_exit"] == 2 and len(data["real_rows"]) == 1
        assert data["real_rows"][0]["checks"]["s"] and not data["real_rows"][0]["checks"]["phase"]
        assert (directory / "real-60.npz").is_file() and progress["real_maps_reserved"] == 1


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
