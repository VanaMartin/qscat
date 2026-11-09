"""Synthetic boundary-only assembly and local-view guard controls; no CO data."""

import hashlib
import json
import unittest
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co import compressed_a2_pole_contrast as dispatch
from projects.ukrmol_co import local_pole_inputs


def fixture():
    labels = [[i + 1, 1 if i < 6 else 2, 1, 1] for i in range(298)]
    info = dict(
        target_components=48,
        boundary_channels=298,
        boundary_poles=2048,
        multipole_orders=2,
        symmetry_charge=[1, 2, 0, 0],
        radius_bohr=18,
        channel_labels=labels,
        target_energies_hartree=[-112, -111.75] + [-111] * 46,
        channel_thresholds_rydberg=[0] * 6 + [0.5] * 292,
    )
    packed = np.zeros((2, 298 * 299 // 2), dtype="<f8")
    packed[0, :3] = [2, 3, 5]
    packed[1, :3] = [7, 11, 13]
    raw = [b""] * 5 + [
        packed.tobytes(),
        np.arange(2048, dtype="<f8").tobytes(),
        np.zeros((2048, 298), dtype="<f8").tobytes(),
    ]
    return info, raw


class LocalInputs(unittest.TestCase):
    def test_full_boundary_assembly_needs_no_k_grid(self):
        info, raw = fixture()
        with (
            patch.object(local_pole_inputs, "read_interface", return_value=info),
            patch.object(local_pole_inputs, "records", return_value=raw) as read,
        ):
            model = local_pole_inputs.load_boundary_model(Path("no-k-file"))
        self.assertEqual(model["amplitudes"].shape, (2048, 298))
        self.assertNotIn("saved", model)
        self.assertNotIn("energies", model)
        read.assert_called_once_with(Path("no-k-file/fort.21"))
        np.testing.assert_allclose(model["c2"][:2, :2], [[4, 3], [3, 7]], atol=0, rtol=0)
        np.testing.assert_allclose(model["c3"][:2, :2], [[7, 11], [11, 13]], atol=0, rtol=0)
        np.testing.assert_allclose(model["thresholds"][6:], 0.5, atol=0, rtol=0)

    def test_wrong_layout_radius_or_threshold_rejects(self):
        info, raw = fixture()
        for key, value in [
            ("target_components", 64),
            ("boundary_poles", 4096),
            ("radius_bohr", 17),
            ("channel_thresholds_rydberg", [0] * 298),
        ]:
            changed = deepcopy(info)
            changed[key] = value
            with (
                patch.object(local_pole_inputs, "read_interface", return_value=changed),
                patch.object(local_pole_inputs, "records", return_value=raw),
                self.assertRaises((ValueError, AssertionError)),
            ):
                local_pole_inputs.load_boundary_model(Path("unused"))

    def test_dispatch_requires_explicit_view_and_never_equilibrium(self):
        c = dict(case="stretched/nine/B1", model_view="local-cas11-nine-2048-v1")
        with patch.object(dispatch, "load_boundary_model", return_value="boundary") as load:
            self.assertEqual(dispatch.load_pole_model(Path("raw/stretched/nine/B1"), c), "boundary")
            load.assert_called_once()
        for changed in [dict(c, model_view=None), dict(c, case="equilibrium/nine/B1")]:
            with self.assertRaises(ValueError):
                dispatch.load_pole_model(Path("raw/stretched/nine/B1"), changed)

    def test_local_guard_requires_own_two_development_rows(self):
        c = dict(
            case="stretched/nine/B1",
            model_view="local-cas11-nine-2048-v1",
            real_guard_result="guard.json",
            input_sha256={},
            real_guard_budgets=dict(
                real_s_absolute_budget=1e-6,
                real_phase_sum_budget_rad=1e-6,
                real_invariant_absolute_budget=1e-6,
                real_flux_absolute_budget=1e-7,
            ),
        )
        with TemporaryDirectory() as directory:
            root = Path(directory)
            for indices in [[16, 20], [14, 20], [60, 80], [16]]:
                data = dict(
                    case=c["case"],
                    scientific_exit=0,
                    real_rows=[
                        dict(
                            index=i,
                            completed=True,
                            s_absolute_error=0,
                            phase_sum_error_mod_pi_rad=0,
                            invariants=dict(
                                flux_defect=0,
                                wronskian_residual=0,
                                unitarity_defect=0,
                                reciprocity_defect=0,
                            ),
                        )
                        for i in indices
                    ],
                )
                payload = json.dumps(data).encode()
                (root / "guard.json").write_bytes(payload)
                c["input_sha256"]["guard.json"] = hashlib.sha256(payload).hexdigest()
                if indices == [16, 20]:
                    dispatch.require_real_guard(root, c)
                else:
                    with self.assertRaises(ValueError):
                        dispatch.require_real_guard(root, c)


if __name__ == "__main__":
    unittest.main()
