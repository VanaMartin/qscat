"""Property-based tests (hypothesis) for the dimension-general primitives.

The D-general library layer (`qscat.linalg`, `qscat.dvr`, `qscat.ecs`) has
invariants that hold for *any* valid input, which is exactly what property-based
testing checks. These complement the fixed-input differential/analytic tests and
justify hypothesis as a dependency.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from hypothesis import given, settings
from hypothesis import strategies as st
from qscat.dvr import gll_nodes_weights
from qscat.ecs import ecs_map
from qscat.linalg import kron_sum


@given(
    x=st.lists(st.floats(min_value=0.0, max_value=200.0), min_size=1, max_size=20),
    r0=st.floats(min_value=1.0, max_value=100.0),
    theta=st.floats(min_value=0.0, max_value=45.0),
)
@settings(max_examples=100, deadline=None)
def test_ecs_map_identity_on_interior_and_at_zero_angle(
    x: list[float], r0: float, theta: float
) -> None:
    xa = np.asarray(x, dtype=np.float64)
    z = ecs_map(xa, r0, theta)
    # Interior (x <= R0): the map is the identity — real, equal to x.
    inside = xa <= r0
    assert np.allclose(z[inside].imag, 0.0)
    assert np.allclose(z[inside].real, xa[inside])
    # theta == 0: identity everywhere.
    z0 = ecs_map(xa, r0, 0.0)
    assert np.allclose(z0.imag, 0.0)
    assert np.allclose(z0.real, xa)


@given(n=st.integers(min_value=2, max_value=40))
@settings(max_examples=50, deadline=None)
def test_gll_weights_sum_to_interval_and_nodes_bracketed(n: int) -> None:
    nodes, weights = gll_nodes_weights(n)
    assert nodes.shape == weights.shape == (n,)
    # GLL quadrature on (-1, 1): weights sum to the interval length, 2.
    assert np.isclose(weights.sum(), 2.0)
    # Nodes are ordered and lie within the closed interval, endpoints included.
    assert np.all(np.diff(nodes) > 0)
    assert np.isclose(nodes[0], -1.0) and np.isclose(nodes[-1], 1.0)


@given(
    a=st.lists(st.floats(-5, 5), min_size=1, max_size=5),
    b=st.lists(st.floats(-5, 5), min_size=1, max_size=5),
)
@settings(max_examples=100, deadline=None)
def test_kron_sum_trace_identity(a: list[float], b: list[float]) -> None:
    # For diagonal blocks, kron_sum(A,B) = A (x) I + I (x) B, so
    # trace = n_B * trace(A) + n_A * trace(B), for A,B of any size.
    A = sp.diags(np.asarray(a, dtype=np.float64))
    B = sp.diags(np.asarray(b, dtype=np.float64))
    ks = kron_sum([A, B])
    assert ks.shape == (len(a) * len(b), len(a) * len(b))
    expected = len(b) * sum(a) + len(a) * sum(b)
    assert np.isclose(ks.diagonal().sum(), expected)
