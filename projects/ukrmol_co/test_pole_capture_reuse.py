"""Exact saved-map reuse rejects changed energy, model identity and method bytes."""

import hashlib
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from projects.ukrmol_co.compressed_pole_reference import digest, reuse_evaluation


def test_exact_pinned_map_reuse_and_no_nearest_point_substitution():
    with TemporaryDirectory() as directory:
        root = Path(directory)
        source = root / "saved.npz"
        output = root / "out"
        output.mkdir()
        z = 0.2 - 0.04j
        np.savez(
            source,
            q=np.eye(298),
            y=1j * np.eye(298),
            a=np.empty((0, 298)),
            z_hartree=z,
            total_hartree=-112 + z,
        )
        row = dict(
            kind="root",
            amplitudes=False,
            z_hartree=[z.real, z.imag],
            setting=[80, 30, 8, 1e-10],
            file="capture.npz",
        )
        model_pin = hashlib.sha256(b"synthetic native inventory").hexdigest()
        methods = {
            n: digest(Path(__file__).with_name(n))
            for n in ["outgoing_map.py", "native_boundary.py"]
        }
        sign, logabs = np.linalg.slogdet((1 - 1j) * np.eye(298))
        prior = dict(
            case="interior/nine/B1",
            reference_energy_rydberg=0.4,
            image_id="synthetic-runtime",
            path=str(source),
            z_hartree=row["z_hartree"],
            setting=row["setting"],
            amplitudes=False,
            input_manifest_sha256=model_pin,
            method_sha256=methods,
            diagnostics=dict(wall_seconds=3),
            det_phase=float(np.angle(sign)),
            log_abs_det_m=float(logabs),
        )
        c = dict(
            case="interior/nine/B1",
            reference_energy_rydberg=0.4,
            image="synthetic-runtime",
            raw_directory="raw",
            input_sha256={"raw/input-manifest.json": model_pin, str(source): digest(source)},
            reused_evaluations=[prior],
        )
        result = reuse_evaluation(row, c, output, -112)
        assert result["completed"] and result["new_ode_call"] is False
        assert (output / row["file"]).read_bytes() == source.read_bytes()
        assert (
            reuse_evaluation(dict(row, z_hartree=[z.real + 1e-12, z.imag]), c, output, -112) is None
        )
        assert reuse_evaluation(dict(row, kind="contour"), c, output, -112) is None
        assert reuse_evaluation(dict(row, amplitudes=True), c, output, -112) is None
        assert reuse_evaluation(row, dict(c, case="interior/nine/B2"), output, -112) is None
        for change in ["map", "model", "method", "diagnostics", "reference", "image"]:
            altered = dict(prior)
            pins = dict(c["input_sha256"])
            if change == "map":
                pins[str(source)] = "0" * 64
            elif change == "model":
                altered["input_manifest_sha256"] = "0" * 64
            elif change == "method":
                altered["method_sha256"] = dict(methods, **{"outgoing_map.py": "0" * 64})
            elif change == "diagnostics":
                altered["det_phase"] = 1
            elif change == "reference":
                altered["reference_energy_rydberg"] = 0.5
            else:
                altered["image_id"] = "different-runtime"
            with unittest.TestCase().assertRaises(ValueError):
                reuse_evaluation(
                    row, dict(c, input_sha256=pins, reused_evaluations=[altered]), output, -112
                )


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
