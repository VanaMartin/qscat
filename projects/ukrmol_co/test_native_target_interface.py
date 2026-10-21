"""Analytic operator-bound controls for saved target integral differences."""

import copy

import numpy as np
import pytest

from projects.ukrmol_co.native_target_interface import compare_integrals, compare_target


def blocks():
    geometry = np.array([[0, 0, 0, 6], [0, 0, 2.1323, 8]], dtype=float)
    source = dict(
        total=13,
        counts=np.array([7, 3, 3, 0, 0, 0, 0, 0]),
        mapping=np.arange(1, 14),
        geometry=geometry,
        overlap=np.eye(13),
        h=np.zeros((13, 13)),
        eri=np.zeros((13,) * 4),
        nuclear=48 / 2.1323,
    )
    full = copy.deepcopy(source)
    full.update(
        total=154,
        counts=np.array([58, 37, 37, 22, 0, 0, 0, 0]),
        mapping=np.array([*range(1, 8), 59, 60, 61, 96, 97, 98]),
        geometry=np.r_[geometry, [[0, 0, 0, 0]]],
    )
    return source, full


def test_spin_counted_triangle_bound_not_just_largest_integral():
    source, full = blocks()
    full["h"][0, 1] = full["h"][1, 0] = 1e-12
    full["eri"][0, 0, 0, 0] = 2e-12
    result = compare_integrals(source, full)
    assert result["operator_bound_hartree"] == pytest.approx(8e-12, rel=1e-15, abs=0)
    assert result["target_integral_equivalence_verdict"] is True


@pytest.mark.parametrize(
    "defect", ["aggregate", "element", "overlap", "geometry", "nan", "symmetry", "mapping"]
)
def test_invalid_saved_target_block_does_not_pass(defect):
    source, full = blocks()
    if defect == "aggregate":
        full["eri"][:] = 1e-11
    elif defect == "element":
        full["h"][0, 0] = 1e-6
    elif defect == "overlap":
        full["overlap"][0, 0] += 1e-8
    elif defect == "geometry":
        full["geometry"][0, 0] += 1e-5
    elif defect == "nan":
        full["eri"][0, 0, 0, 0] = np.nan
    elif defect == "symmetry":
        full["eri"][0, 1, 0, 0] = 1e-5
    else:
        full["mapping"][7] = 58
    with pytest.raises((ValueError, AssertionError)):
        compare_integrals(source, full)


def target_fixture(tmp_path, defect=None):
    import hashlib

    from scipy.io import FortranFile

    from projects.ukrmol_co.test_native_csf_verification import cidata_fixture

    cidata = tmp_path / "fort.26"
    pairs = cidata_fixture(cidata)
    native = tmp_path / "target.bin"
    expected = []
    signs = np.ones(16, dtype="<i8")
    signs[0] = -1
    with FortranFile(native, "w") as stream:
        for index, (sector, _) in enumerate(pairs.items()):
            unit = (710, 730, 711, 731, 712, 732, 713, 733)[index]
            with FortranFile(tmp_path / f"fort.{unit}", "w") as csf:
                csf.write_record(np.zeros(304, dtype="u1"))
                csf.write_record(signs)
            roots = 13 if index < 2 else 8
            expected.append(
                dict(
                    sector=sector,
                    dimension=16,
                    roots=roots,
                    sign_int64_sha256=hashlib.sha256(signs.tobytes()).hexdigest(),
                )
            )
            stream.write_record(np.array([index + 1, 16, roots], dtype="<i8"))
            energies = np.arange(roots, dtype="<f8")
            if defect == "energies":
                energies[0] += 1e-12
            stream.write_record(energies)
            vectors = np.eye(roots, 16) * signs
            if defect == "vectors":
                vectors[0, 0] = 1
            stream.write_record(vectors.reshape(-1))
        if defect == "trailing":
            stream.write_record(np.zeros(1))
    return native, cidata, tmp_path, expected


def test_native_target_output_compared_byte_exactly_for_all_roots(tmp_path):
    result = compare_target(*target_fixture(tmp_path))
    assert sum(s["roots"] for s in result) == 74
    assert all(s["exact_native_coefficients"] and s["exact_native_energies"] for s in result)


@pytest.mark.parametrize("defect", ["energies", "vectors", "trailing"])
def test_malformed_native_import_output_rejects(tmp_path, defect):
    with pytest.raises(ValueError):
        compare_target(*target_fixture(tmp_path, defect))
