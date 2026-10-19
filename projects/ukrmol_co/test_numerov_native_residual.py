"""Primitive precision, array-layout and cancellation checks against mpmath."""

import mpmath as mp
import numpy as np
import pytest

from projects.ukrmol_co.numerov_native_residual import native_residual, numpy_residual


def require_native():
    native = pytest.importorskip("qscat_extended_kernels")
    if not native.capability()[3] or np.finfo(np.longdouble).nmant + 1 != 64:
        pytest.skip("Host has no qualified64-bit-significand native longdouble")


def deterministic_case():
    dtype = np.longdouble
    a = ((np.arange(77).reshape(7, 11) * 7 % 23) - 11).astype(dtype) / 8
    x = ((np.arange(55).reshape(11, 5) * 11 % 19) - 9).astype(dtype) / 16
    a += dtype(2) ** -60
    b = np.full((7, 5), dtype(3) / 7)
    return a, b, x


def mp_reference(a, b, x):
    with mp.workdps(100):
        aa, bb, xx = [
            mp.matrix(
                [
                    [
                        mp.mpf(int(v.as_integer_ratio()[0])) / int(v.as_integer_ratio()[1])
                        for v in row
                    ]
                    for row in value
                ]
            )
            for value in [a, b, x]
        ]
        return bb - aa * xx


def scaled_error(result, reference, a, b, x):
    with mp.workdps(100):
        error = max(
            abs(mp.mpf(str(result[i, j])) - reference[i, j])
            for i in range(len(result))
            for j in range(result.shape[1])
        )

        def norm(value):
            return mp.mpf(str(np.max(np.sum(abs(value), axis=1))))

        denominator = norm(a) * norm(x) + norm(b)
        return float(error / denominator) if denominator else float(error)


def test_numpy_reference_cancellation_signal():
    dtype = np.longdouble
    a = np.array([[1]], dtype=dtype) + dtype(2) ** -60
    b = x = np.ones((1, 1), dtype=dtype)
    expected = -(dtype(2) ** -60) if np.finfo(dtype).nmant + 1 >= 64 else dtype(0)
    np.testing.assert_allclose(numpy_residual(a, b, x), [[expected]], atol=0, rtol=0)


@pytest.mark.parametrize("layout", ["contiguous", "fortran", "transpose", "reversed"])
def test_native_dense_layout_and_inputs(layout):
    require_native()
    a, b, x = deterministic_case()
    if layout == "fortran":
        a, b, x = [np.asfortranarray(v) for v in [a, b, x]]
    if layout == "transpose":
        a, b, x = [v.T.copy().T for v in [a, b, x]]
    if layout == "reversed":
        a, b, x = a[:, ::-1], b[::-1, ::-1], x[::-1, ::-1]
    before = [v.copy() for v in [a, b, x]]
    actual = native_residual(a, b, x)
    reference = mp_reference(a, b, x)
    assert scaled_error(actual, reference, a, b, x) <= 32 * np.finfo(np.longdouble).eps
    assert actual.dtype == np.dtype(np.longdouble) and actual.flags.owndata
    for actual_input, old in zip([a, b, x], before, strict=True):
        np.testing.assert_array_equal(actual_input, old)


def test_native_low_bits_and_cancellation():
    require_native()
    dtype = np.longdouble
    a = np.array([[1]], dtype=dtype) + dtype(2) ** -60
    b = x = np.ones((1, 1), dtype=dtype)
    np.testing.assert_allclose(native_residual(a, b, x), [[-(dtype(2) ** -60)]], atol=0, rtol=0)
    a = np.array([[dtype(2) ** 60, 1, -(dtype(2) ** 60)]], dtype=dtype)
    np.testing.assert_allclose(
        native_residual(a, np.zeros((1, 1), dtype=dtype), np.ones((3, 1), dtype=dtype)),
        [[-1]],
        atol=0,
        rtol=0,
    )


@pytest.mark.parametrize("bad", ["shape", "nonfinite", "dtype"])
def test_invalid_operands_rejected_before_native_call(bad):
    a, b, x = deterministic_case()
    if bad == "shape":
        b = b[:, :-1]
    if bad == "nonfinite":
        a[0, 0] = np.nan
    if bad == "dtype":
        x = x.astype(np.float32)
    with pytest.raises(ValueError, match="Finite compatible"):
        native_residual(a, b, x)
