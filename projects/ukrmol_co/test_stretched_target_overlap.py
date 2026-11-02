"""CAS12 native-address/sign controls and finite-subspace decision counterexamples."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from scipy.io import FortranFile

from projects.ukrmol_co.ci_overlap import frozen_core_overlaps
from projects.ukrmol_co.ensemble_frame_overlaps import full_minor_overlaps
from projects.ukrmol_co.native_target_determinants import determinant_index, expand, parity
from projects.ukrmol_co.stretched_dz_pole_owner import command
from projects.ukrmol_co.stretched_target_overlap import CONTRACT, canonical_order, diagnose


def fixture(path, bad_phase=False):
    header = bytearray(304)
    words = np.zeros(23, dtype="<i8")
    for i, v in {0: 0, 5: 14, 6: 28, 7: 8, 8: 14, 9: 8, 10: 1, 11: 3, 12: 2, 21: 8, 22: -1}.items():
        words[i] = v
    reference = list(range(1, 15))
    packed = [0]
    for particle in range(16, 29, 2):
        packed += [1, 14, particle]
    words[13] = len(packed)
    header[120:] = words.tobytes()
    phases = np.ones(8, dtype="<i8")
    if bad_phase:
        phases[0] = -1
    # The tail has 1+5*nsym words: zero, target, inventory/zeros, target/zeros.
    metadata = np.r_[
        [8, 3, 3],
        reference,
        np.ones(8, dtype=int),
        0,
        [8, 3, 3],
        [8, 3, 3],
        [0, 0, 0],
        [8, 3, 3],
        [0, 0, 0],
    ]
    offsets = np.r_[np.arange(1, 10), [1, 2, 5, 8, 11, 14, 17, 20, 23]]
    with FortranFile(path, "w") as f:
        for value in [
            np.frombuffer(header, "u1"),
            phases,
            metadata.astype("<i8"),
            offsets.astype("<i8"),
            np.array(packed, dtype="<i8"),
            np.ones(8, dtype="<f8"),
        ]:
            f.write_record(value)
    return dict(phases=phases, vectors=np.eye(8))


def test_fermion_sorting_signs_and_alpha_beta_reordering():
    assert parity([2, 0, 1]) == 1 and parity([1, 0, 2]) == -1
    a, b, sign = determinant_index(list(range(1, 15)), 0)
    assert a == b == 31 and sign == (-1) ** 21
    with unittest.TestCase().assertRaises(ValueError):
        determinant_index([1] * 14, 0)
    with unittest.TestCase().assertRaises(ValueError):
        determinant_index(list(range(3, 17)), 0)


def test_positive_exact_native_layout_reexpands_without_normalization_or_solves():
    with TemporaryDirectory() as tmp:
        p = Path(tmp) / "csf"
        native = fixture(p)
        ci, row = expand(p, native, 0, "A1")
        assert row == dict(csfs=8, determinant_terms=8, determinant_shape=[792, 792])
        np.testing.assert_allclose(
            ci.reshape(8, -1) @ ci.reshape(8, -1).T, np.eye(8), atol=1e-14, rtol=0
        )
        assert np.count_nonzero(ci) == 8


def test_native_phase_and_wrong_sector_cannot_be_silently_repaired():
    with TemporaryDirectory() as tmp:
        p = Path(tmp) / "csf"
        native = fixture(p, True)
        with unittest.TestCase().assertRaisesRegex(ValueError, "canonical phase"):
            expand(p, native, 0, "A1")
        with unittest.TestCase().assertRaises(ValueError):
            expand(p, native, 0, "B1")


def test_dominance_screen_handles_permutations_and_lost_manifold_directions():
    p = np.eye(8)[:, [1, 0, 2, 3, 4, 5, 6, 7]]
    assert diagnose(p)["retained_eight_dominance_passed"] is True
    assert diagnose(0.7 * np.eye(8))["retained_eight_dominance_passed"] is False
    with unittest.TestCase().assertRaises(ValueError):
        diagnose(1.1 * np.eye(8))


def test_explicit_orbital_inventory_drives_canonical_metric_order():
    labels = [
        ("A1", 1),
        ("A1", 2),
        ("A1", 3),
        ("B1", 1),
        ("B2", 1),
        ("A1", 4),
        ("B1", 2),
        ("B2", 2),
        ("A1", 5),
        ("B1", 3),
        ("B2", 3),
        ("A1", 6),
        ("A1", 7),
        ("A1", 8),
    ]
    target = dict(orbitals=[dict(irrep=s, index_in_irrep=i) for s, i in labels])
    assert canonical_order(target) == [0, 1, 2, 5, 8, 11, 12, 13, 3, 6, 9, 4, 7, 10]
    target["orbitals"][0]["irrep"] = "A2"
    with unittest.TestCase().assertRaises(ValueError):
        canonical_order(target)


def test_four_core_owner_has_no_native_or_network_successor():
    c = json.loads(Path(CONTRACT).read_text())
    cmd = command(Path("/root"), Path("/source"), Path("/out"), c, CONTRACT)
    assert "--cpuset-cpus=0-3" in cmd and "--memory=8g" in cmd and "--network=none" in cmd
    assert "projects.ukrmol_co.stretched_target_overlap" in cmd[-1]
    compile(cmd[-1], "<target-overlap-owner>", "exec")


def test_nonorthogonal_frozen_core_formula_against_full_occupied_minors():
    rng = np.random.default_rng(78)
    frames = [
        np.linalg.qr(rng.normal(size=(8, 5)) + 1j * rng.normal(size=(8, 5)))[0] for _ in range(2)
    ]
    metric = frames[0].conj().T @ frames[1]
    left = rng.normal(size=(2, 6, 4)) + 1j * rng.normal(size=(2, 6, 4))
    right = rng.normal(size=(3, 6, 4)) + 1j * rng.normal(size=(3, 6, 4))
    actual = frozen_core_overlaps(left, right, metric, 1, (2, 1))
    expected = full_minor_overlaps(left, right, metric, 1, (2, 1))
    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(
        actual,
        frozen_core_overlaps(right, left, metric.conj().T, 1, (2, 1)).conj().T,
        atol=1e-12,
        rtol=1e-12,
    )


def test_core_factor_and_singular_core_stop():
    metric = np.diag([0.8, 0.9, 1.0])
    ci = np.ones((1, 1, 1))
    np.testing.assert_allclose(
        frozen_core_overlaps(ci, ci, metric, 2, (0, 0)), [[(0.8 * 0.9) ** 2]], atol=1e-14, rtol=0
    )
    with unittest.TestCase().assertRaises(ValueError):
        frozen_core_overlaps(ci, ci, np.diag([0.0, 1.0, 1.0]), 2, (0, 0))


def test_root_phase_covariance_does_not_normalize_lost_projection():
    roots = np.eye(4).reshape(4, 2, 2)
    phase = np.exp(1j * np.arange(4))
    actual = frozen_core_overlaps(
        roots * phase[:, None, None], roots, np.diag([0.8, 1.0, 1.0]), 1, (1, 1)
    )
    np.testing.assert_allclose(actual, 0.8**2 * np.diag(phase.conj()), atol=1e-13, rtol=0)


def test_frozen_source_prerequisites_exist_and_match():
    from projects.ukrmol_co.stretched_basis_properties import digest

    c = json.loads(Path(CONTRACT).read_text())
    for name, pin in c["source_sha256"].items():
        assert digest(Path(name)) == pin, name


if __name__ == "__main__":
    functions = [
        fn for name, fn in list(globals().items()) if name.startswith("test_") and callable(fn)
    ]
    suite = unittest.TestSuite(unittest.FunctionTestCase(fn) for fn in functions)
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
